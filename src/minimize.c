#include <stdio.h>
#include <stdlib.h>
#include <jansson.h>

#include "mol2/gbsa.h"
#include "mol2/minimize.h"
#include "mol2/pdb.h"

#include "utils.h"
#include "setup.h"
#include "parse_options.h"
#include "energy.h"

#define __TOL__ 5E-4


int main(int argc, char **argv) {
    VERBOSITY = DEBUG;

    mol_enable_floating_point_exceptions();

    bool opts_error;
    struct options opts = options_populate_from_argv(argc, argv, &opts_error);
    if (opts_error) {
        ERR_MSG("Couldn't parse options");
        options_free(opts);
        usage_message(argv);
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

    struct energy_prms* min_prms;
    size_t nstages;
    if (!energy_prms_populate_from_options(&min_prms, &nstages, opts)) {
        ERR_MSG("Couldn't populate energy parameters");
        options_free(opts);
        mol_atom_group_list_free(ag_list);
        exit(EXIT_FAILURE);
    }

    FILE* out_pdb;
    FOPEN_ELSE(out_pdb, opts.out_pdb, "w") {
        options_free(opts);
        mol_atom_group_list_free(ag_list);
        energy_prms_free(&min_prms, nstages);
        exit(EXIT_FAILURE);
    }

    json_t* json_log_total = json_array();

    for (size_t modeli = 0; modeli < ag_list->size; modeli++) {
        INFO_MSG("Started model %zu", modeli);

        if (ag_list->size > 1) {
            fprintf(out_pdb, "MODEL %zu\n", (modeli + 1));
        }

        json_t* json_log_model = json_array();

        struct mol_atom_group *ag = &ag_list->members[modeli];

        // Setup fixed atoms and nblists
        struct agsetup ag_setup;
        mol_fixed_init(ag);
        init_nblst(ag, &ag_setup);

        // Run stages of minimization
        for (size_t stage_id = 0; stage_id < nstages; stage_id++) {
            INFO_MSG("Stage %zu", stage_id);
            json_t* json_log_stage = json_object();

            struct energy_prms *stage_prms = &min_prms[stage_id];
            stage_prms->ag = ag;

            if (stage_prms->fixed) {
                mol_fixed_update(ag, stage_prms->fixed->natoms, stage_prms->fixed->atoms);
            } else {
                mol_fixed_update(ag, 0, NULL);
            }

            update_nblst(ag, &ag_setup);

            // Set up GBSA
            if (stage_prms->gbsa) {
                struct acesetup ace_setup;
                ace_setup.efac = 0.5;
                ace_ini(ag, &ace_setup);
                ace_fixedupdate(ag, &ag_setup, &ace_setup);
                ace_updatenblst(&ag_setup, &ace_setup);
                stage_prms->ace_setup = &ace_setup;
            } else {
                stage_prms->ace_setup = NULL;
            }

            stage_prms->ag_setup = &ag_setup;

            // Minimize energy
            if (!stage_prms->score_only) {
                stage_prms->json_log = json_array();
                mol_minimize_ag(MOL_LBFGS, stage_prms->nsteps, __TOL__, ag, (void *) stage_prms, energy_func);

                // Record energy every time it was evaluated
                if (stage_prms->json_log_setup.print_step) {
                    json_object_set_new(json_log_stage, "steps", stage_prms->json_log);
                } else {
                    json_decref(stage_prms->json_log);
                }
                stage_prms->json_log = NULL;
            }

            // Record final energy
            if (stage_prms->json_log_setup.print_stage) {
                stage_prms->json_log = json_array();
                energy_func((void *) stage_prms, NULL, NULL, 0, 0);
                json_t *final_energy = json_deep_copy(json_array_get(stage_prms->json_log, 0));
                json_object_set_new(json_log_stage, "final", final_energy);
                json_decref(stage_prms->json_log);
                stage_prms->json_log = NULL;
            }

            json_array_append_new(json_log_model, json_log_stage);

            if (stage_prms->ace_setup) {
                destroy_acesetup(stage_prms->ace_setup);
            }
        }

        mol_fwrite_pdb(out_pdb, ag);
        json_array_append_new(json_log_total, json_log_model);

        if (ag_list->size > 1) {
            fprintf(out_pdb, "ENDMDL\n");
        }

        destroy_agsetup(&ag_setup);

        INFO_MSG("Finished model %zu", modeli);
    }

    if (json_dump_file(json_log_total, opts.out_json, JSON_INDENT(4)) != 0) {
        ERR_MSG("Couldn't write to file %s", opts.out_json);
        return EXIT_FAILURE;
    }
    json_decref(json_log_total);

    // Free everything
    mol_atom_group_list_free(ag_list);
    energy_prms_free(&min_prms, nstages);
    options_free(opts);

    INFO_MSG("Completed\n");
    return EXIT_SUCCESS;
}
