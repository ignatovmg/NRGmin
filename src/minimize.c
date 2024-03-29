#include <stdio.h>
#include <stdlib.h>
#include <jansson.h>
#include <signal.h>
#include <setjmp.h>
#include <fenv.h>

#include "mol2/gbsa.h"
#include "mol2/genborn.h"
#include "mol2/minimize.h"
#include "mol2/pdb.h"

#include "utils.h"
#include "setup.h"
#include "parse_options.h"
#include "energy.h"

#ifdef OPENMP
#include <omp.h>
#endif


jmp_buf FE_RETURN_LOCATION;

/**
 * Handler for SIGFPE (floating point exception)
 */
void overflow_handler(__attribute__((unused)) int signal_number) {
    if (feclearexcept(FE_ALL_EXCEPT) != 0){
        WRN_MSG("Couldn't clear exceptions");
    }

    longjmp(FE_RETURN_LOCATION, 1);
}


int main(int argc, char **argv) {
    VERBOSITY = DEBUG;

    #ifndef OPENMP
    // For serial version, SIGFPE will be caught and "exception" field will be added
    // to the corresponding model entry in output json file. In OMP version SIGFPE
    // will be switched off
    mol_enable_floating_point_exceptions();
    signal(SIGFPE, overflow_handler);
    #endif

    bool opts_error;
    struct options opts = options_populate_from_argv(argc, argv, &opts_error);
    if (opts_error) {
        ERR_MSG("Couldn't parse options");
        options_free(opts);
        exit(EXIT_FAILURE);
    }
    if (opts.help) {
        options_free(opts);
        usage_message(argv);
        exit(EXIT_SUCCESS);
    }

    struct mol_atom_group_list* ag_list = mol_atom_group_list_from_options(&opts);
    if (!ag_list) {
        ERR_MSG("Couldn't create atom group for minimization");
        options_free(opts);
        exit(EXIT_FAILURE);
    }

    // Create separate json array for each model
    json_t** json_log_total = calloc(ag_list->size, sizeof(json_t *));
    for (size_t i = 0; i < ag_list->size; i++) {
        json_log_total[i] = json_array();
    }

    #ifdef OPENMP
    if (opts.num_threads > 0) {
        omp_set_num_threads(opts.num_threads);
    }

    #pragma omp parallel shared(ag_list, json_log_total, opts)
    {
        int num_threads = omp_get_num_threads();
        int thread_id = omp_get_thread_num();

        if (thread_id == 0) {
            DEBUG_MSG("Using %i threads", num_threads);
        }
    #else
        int thread_id = 0;
    #endif

        // Create local copy of energy prms for each thread
        // TODO: move this outside of omp parallel and implement
        //       energy_prms_copy() function in setup.c
        struct energy_prms* min_prms;
        size_t nstages;
        if (!energy_prms_populate_from_options(&min_prms, &nstages, opts)) {
            if (thread_id == 0) {
                ERR_MSG("Couldn't populate energy parameters");
                options_free(opts);
                mol_atom_group_list_free(ag_list);
            }
            exit(EXIT_FAILURE);
        }

        #ifdef OPENMP
        // Wait till all threads successfully create min_prms
        #pragma omp barrier
        #pragma omp for
        #endif
        for (size_t modeli = 0; modeli < ag_list->size; modeli++) {
            INFO_MSG("Started model %zu", modeli);

            json_t *json_log_model = json_log_total[modeli];

            struct mol_atom_group *ag = &ag_list->members[modeli];

            // Setup fixed atoms and nblists
            mol_fixed_init(ag);

            struct agsetup ag_setup;
            init_nblst(ag, &ag_setup);

            // Run stages of minimization
            for (size_t stage_id = 0; stage_id < nstages; stage_id++) {
                json_t *json_log_stage = json_object();

                struct energy_prms *stage_prms = calloc(1, sizeof(struct energy_prms));
                *stage_prms = min_prms[stage_id];
                stage_prms->ag = ag;
                ag_setup.nblst->nbcut = stage_prms->nbcut;

                if (stage_prms->fixed) {
                    mol_fixed_update(ag, stage_prms->fixed->setups[modeli]->natoms, stage_prms->fixed->setups[modeli]->atoms);
                } else {
                    mol_fixed_update(ag, 0, NULL);
                }
                mol_nblst_update_fixed(ag, &ag_setup);

                // Set up ACE
                struct acesetup ace_setup;
                if (stage_prms->ace) {
                    ace_setup.efac = stage_prms->ace_efac;
                    ace_ini(ag, &ace_setup);
                    ace_fixedupdate(ag, &ag_setup, &ace_setup);
                    ace_updatenblst(&ag_setup, &ace_setup);
                    stage_prms->ace_setup = &ace_setup;
                } else {
                    stage_prms->ace_setup = NULL;
                }

                // Set up GBSA
                if (stage_prms->gbsa) {
                    if (!mol_init_gb(ag)) {
                        ERR_MSG("Failed to initialize GBSA for model %zu", modeli);
                        exit(EXIT_FAILURE);
                    }
                }

                stage_prms->ag_setup = &ag_setup;
                bool sigfpe_caught = false;

                // Minimize energy
                if (!stage_prms->score_only) {
                    stage_prms->json_log = json_array();

                    if (setjmp(FE_RETURN_LOCATION) == 1) {
                        // If minimization failed with SIGFPE the execution jumps to this location
                        WRN_MSG("Floating point exception detected for model %zu during minimization", modeli);
                        sigfpe_caught = true;
                        if (stage_prms->json_log && (json_array_size(stage_prms->json_log) > 0)) {
                            // Add "exception" to the last recorded energy set
                            json_t *last_energy_dict = json_array_get(stage_prms->json_log, json_array_size(stage_prms->json_log) - 1);
                            json_object_set_new(last_energy_dict, "exception", json_string("Floating point exception during minimization"));
                        }
                    } else {
                        mol_minimize_ag(MOL_LBFGS, stage_prms->nsteps, stage_prms->tol, ag, (void *) stage_prms, energy_func);
                    }

                    // Record energy every time it was evaluated
                    if (stage_prms->json_log_setup.print_step) {
                        json_object_set(json_log_stage, "steps", stage_prms->json_log);
                    }
                    json_decref(stage_prms->json_log);
                    stage_prms->json_log = NULL;
                }

                struct mol_atom_group* ag_copy;
                if (stage_prms->pull_ligand_away) {
                    ag_copy = mol_atom_group_copy(stage_prms->ag);
                    pull_ligand_away(stage_prms->ag, opts.rec_natoms, 200);
                }

                // Record final energy
                if (stage_prms->json_log_setup.print_stage) {
                    stage_prms->json_log = json_array();

                    if (sigfpe_caught) {
                        // if SIGFPE was caught previously during minimization, do not compute the final energy and set "exception" instead
                        json_t* energy_dict = json_object();
                        json_object_set_new(
                                energy_dict,
                                "exception",
                                json_string("Floating point exception during minimization"));
                        json_array_append_new(stage_prms->json_log, energy_dict);
                    } else {
                        if (setjmp(FE_RETURN_LOCATION) == 1) {
                            // If energy evalution failed with SIGFPE, the execution jumps to this location
                            WRN_MSG("Floating point exception detected for model %zu during final energy evaluation", modeli);
                            sigfpe_caught = true;

                            if (stage_prms->json_log) {
                                json_t *last_energy_dict = json_array_get(stage_prms->json_log, 0);
                                json_object_set_new(
                                        last_energy_dict,
                                        "exception",
                                        json_string("Floating point exception during final energy evaluation"));
                            }
                        } else {
                            energy_func((void *) stage_prms, NULL, NULL, 0, 0);
                        }
                    }

                    json_t *final_energy = json_deep_copy(json_array_get(stage_prms->json_log, 0));
                    json_object_set_new(json_log_stage, "final", final_energy);
                    json_decref(stage_prms->json_log);
                    stage_prms->json_log = NULL;
                }

                if (stage_prms->pull_ligand_away) {
                    memcpy(stage_prms->ag->coords, ag_copy->coords, ag_copy->natoms * sizeof(struct mol_vector3));
                    mol_atom_group_free(ag_copy);
                }

                json_array_append_new(json_log_model, json_log_stage);

                if (stage_prms->ace_setup) {
                    destroy_acesetup(stage_prms->ace_setup);
                }
                if (stage_prms->gbsa) {
                    if (!mol_gb_delete_metadata(ag)) {
                        ERR_MSG("Failed to delete GenBorn metadata for model %zu", modeli);
                        exit(EXIT_FAILURE);
                    }
                }
                free(stage_prms);
            }

            destroy_agsetup(&ag_setup);
            INFO_MSG("Finished model %zu", modeli);
        }

        energy_prms_free(&min_prms, nstages);

    #ifdef OPENMP
    }
    #endif

    // Write minimized models
    if (!opts.score_only) {
        FILE *out_pdb;
        FOPEN_ELSE(out_pdb, opts.out_pdb, "w") {
            options_free(opts);
            mol_atom_group_list_free(ag_list);
            exit(EXIT_FAILURE);
        }

        DEBUG_MSG("Writing minimized models to %s", opts.out_pdb);
        for (size_t modeli = 0; modeli < ag_list->size; modeli++) {
            if (ag_list->size > 1) {
                fprintf(out_pdb, "MODEL %zu\n", (modeli + 1));
            }
            mol_fwrite_pdb(out_pdb, &ag_list->members[modeli]);
            if (ag_list->size > 1) {
                fprintf(out_pdb, "ENDMDL\n");
            }
        }
    }

    DEBUG_MSG("Merging energy logs for all models");
    json_t* json_log_merged = json_array();
    for (size_t i = 0; i < ag_list->size; i++) {
        json_array_append_new(json_log_merged, json_log_total[i]);
    }
    DEBUG_MSG("Writing energies to %s", opts.out_json);
    if (json_dump_file(json_log_merged, opts.out_json, JSON_INDENT(4)) != 0) {
        ERR_MSG("Couldn't write to file %s", opts.out_json);
        return EXIT_FAILURE;
    }
    json_decref(json_log_merged);
    free(json_log_total);

    // Free everything
    mol_atom_group_list_free(ag_list);
    options_free(opts);

    INFO_MSG("Completed");
    return EXIT_SUCCESS;
}
