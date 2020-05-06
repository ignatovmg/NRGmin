#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <jansson.h>

#include "mol2/json.h"
#include "mol2/benergy.h"
#include "mol2/gbsa.h"
#include "mol2/icharmm.h"
#include "mol2/minimize.h"
#include "mol2/nbenergy.h"
#include "mol2/pdb.h"
#include "mol2/fitting.h"
#include "mol2/noe.h"

#include "utils.h"
#include "potentials.h"
#include "parse_options.h"

#define __TOL__ 5E-4


static lbfgsfloatval_t energy_func(
        void *restrict prm,
        const double *restrict array,
        double *restrict gradient,
        const int array_size,
        const lbfgsfloatval_t step);


int main(int argc, char **argv) {
    mol_enable_floating_point_exceptions();

    bool parse_error;
    struct options parsed = parse_args(argc, argv, &parse_error);
    if (parse_error) {
        ERR_MSG("Parse error");
        exit(EXIT_FAILURE);
    }
    if (parsed.help) {
        exit(EXIT_SUCCESS);
    }

    struct mol_atom_group_list* ag_list = mol_atom_group_list_from_options(&parsed);
    if (!ag_list) {
        ERR_MSG("Couldn't read atom groups");
        exit(EXIT_FAILURE);
    }

    struct energy_prm* min_prms;
    size_t nstages;
    if (!energy_prm_read(&min_prms, &nstages, parsed)) {
        ERR_MSG("Couldn't fill params");
        exit(EXIT_FAILURE);
    }

    INFO_MSG("Started the main loop\n");

    FILE* out_pdb;
    FOPEN_ELSE(out_pdb, parsed.out_pdb, "w") {
        ERR_MSG("Can't open %s", parsed.out_pdb);
        mol_atom_group_list_free(ag_list);
        exit(EXIT_FAILURE);
    }

    json_t* json_log_total = json_array();

    for (size_t modeli = 0; modeli < ag_list->size; modeli++) {
        INFO_MSG("Started model %zu\n", modeli);

        if (ag_list->size > 1) {
            fprintf(out_pdb, "MODEL %zu\n", (modeli + 1));
        }

        json_t* json_log_model = json_array();

        struct mol_atom_group *ag = &ag_list->members[modeli];
        ag->gradients = calloc(ag->natoms, sizeof(struct mol_vector3));

        for (size_t stage_id = 0; stage_id < nstages; stage_id++) {
            INFO_MSG("Stage %zu\n", stage_id);
            json_t* json_log_stage = json_object();

            struct energy_prm *stage_prms = &min_prms[stage_id];
            stage_prms->ag = ag;

            struct agsetup ag_setup;
            struct acesetup ace_setup;

            // Setup fixed atoms and lists
            mol_fixed_init(ag);

            if (stage_prms->fixed) {
                mol_fixed_update(ag, stage_prms->fixed->natoms, stage_prms->fixed->atoms);
            } else {
                mol_fixed_update(ag, 0, NULL);
            }

            init_nblst(ag, &ag_setup);
            update_nblst(ag, &ag_setup);

            if (stage_prms->gbsa) {
                ace_setup.efac = 0.5;
                ace_ini(ag, &ace_setup);
                ace_fixedupdate(ag, &ag_setup, &ace_setup);
                ace_updatenblst(&ag_setup, &ace_setup);
                stage_prms->ace_setup = &ace_setup;
            } else {
                stage_prms->ace_setup = NULL;
            }

            stage_prms->ag_setup = &ag_setup;

            mol_minimize_ag(MOL_LBFGS, stage_prms->nsteps, __TOL__, ag, (void *) stage_prms, energy_func);

            //free_agsetup(stage_prms->ag_setup);
            //if (stage_prms->ace_setup) {
            //    free_acesetup(stage_prms->ace_setup);
            //}

            // Record energy every time it was evaluated
            json_object_set_new(json_log_stage, "steps", stage_prms->json_log);
            stage_prms->json_log = NULL;

            // Record final energy
            energy_func((void *) stage_prms, NULL, NULL, 0, 0);
            json_t* final_energy = json_deep_copy(json_array_get(stage_prms->json_log, 0));
            json_object_set_new(json_log_stage, "final", final_energy);
            json_decref(stage_prms->json_log);
            stage_prms->json_log = NULL;

            json_array_append_new(json_log_model, json_log_stage);
        }

        mol_fwrite_pdb(out_pdb, ag);
        json_array_append_new(json_log_total, json_log_model);

        if (ag_list->size > 1) {
            fprintf(out_pdb, "ENDMDL\n");
        }
    }

    json_dump_file(json_log_total, parsed.out_json, JSON_INDENT(4));
    json_decref(json_log_total);

    mol_atom_group_list_free(ag_list);

    INFO_MSG("Completed\n");
    return EXIT_SUCCESS;
}


static lbfgsfloatval_t energy_func(
        void *restrict prm,
        const double *restrict array,
        double *restrict gradient,
        const int array_size,
        __attribute__((unused)) const lbfgsfloatval_t step) {

    lbfgsfloatval_t total_energy = 0.0;
    lbfgsfloatval_t term_energy;

    struct energy_prm *energy_prm = (struct energy_prm *) prm;
    if (energy_prm->json_log == NULL) {
        energy_prm->json_log = json_array();
    }

    json_t* energy_dict = json_object();

    if (array != NULL) {
        assert((size_t)array_size == energy_prm->ag->active_atoms->size * 3);
        mol_atom_group_set_actives(energy_prm->ag, array);
    }

    bool updated = check_clusterupdate(energy_prm->ag, energy_prm->ag_setup);
    if (updated) {
        if (energy_prm->ace_setup != NULL) {
            ace_updatenblst(energy_prm->ag_setup, energy_prm->ace_setup);
        }
    }

    mol_zero_gradients(energy_prm->ag);

    if (energy_prm->ace_setup != NULL) {
        term_energy = 0.0;
        aceeng(energy_prm->ag, &term_energy, energy_prm->ace_setup, energy_prm->ag_setup);
        json_object_set_new(energy_dict, "gbsa", json_real(term_energy));
        total_energy += term_energy;
    }

    if (energy_prm->vdw) {
        term_energy = 0.0;
        vdweng(energy_prm->ag, &term_energy, energy_prm->ag_setup->nblst);
        json_object_set_new(energy_dict, "vdw", json_real(term_energy));
        total_energy += term_energy;
    }

    if (energy_prm->vdw03) {
        term_energy = 0.0;
        vdwengs03(
                1.0,
                energy_prm->ag_setup->nblst->nbcof,
                energy_prm->ag, &term_energy,
                energy_prm->ag_setup->nf03,
                energy_prm->ag_setup->listf03);
        json_object_set_new(energy_dict, "vdw03", json_real(term_energy));
        total_energy += term_energy;
    }

    if (energy_prm->bonds) {
        term_energy = 0.0;
        beng(energy_prm->ag, &term_energy);
        json_object_set_new(energy_dict, "bonds", json_real(term_energy));
        total_energy += term_energy;
    }

    if (energy_prm->angles) {
        term_energy = 0.0;
        aeng(energy_prm->ag, &term_energy);
        json_object_set_new(energy_dict, "angles", json_real(term_energy));
        total_energy += term_energy;
    }

    if (energy_prm->dihedrals) {
        term_energy = 0.0;
        teng(energy_prm->ag, &term_energy);
        json_object_set_new(energy_dict, "dihedrals", json_real(term_energy));
        total_energy += term_energy;
    }

    if (energy_prm->impropers) {
        term_energy = 0.0;
        ieng(energy_prm->ag, &term_energy);
        json_object_set_new(energy_dict, "impropers", json_real(term_energy));
        total_energy += term_energy;
    }

    if (energy_prm->sprst_pairs != NULL) {
        term_energy = 0.0;
        pairspring_energy(energy_prm->sprst_pairs, energy_prm->ag, &term_energy);
        json_object_set_new(energy_dict, "pairsprings", json_real(term_energy));
        total_energy += term_energy;
    }

    if (energy_prm->sprst_points != NULL) {
        term_energy = 0.0;
        pointspring_energy(energy_prm->sprst_points, energy_prm->ag, &term_energy);
        json_object_set_new(energy_dict, "pointsprings", json_real(term_energy));
        total_energy += term_energy;
    }

    if (energy_prm->nmr != NULL) {
        mol_noe_calc_peaks(energy_prm->nmr->spec, energy_prm->ag, true);
        mol_noe_calc_energy(
                energy_prm->nmr->spec,
                energy_prm->ag->gradients,
                energy_prm->nmr->weight,
                energy_prm->nmr->power);

        term_energy = energy_prm->nmr->spec->energy;
        total_energy += term_energy;

        json_object_set_new(energy_dict, "noe", json_real(term_energy));
        json_object_set_new(energy_dict, "noe_details", mol_noe_to_json_object(energy_prm->nmr->spec));
    }

    if (energy_prm->fit_prms != NULL) {
        term_energy = mol_fitting_score_aglist(energy_prm->ag,
                                           energy_prm->fit_prms->ag_list,
                                           energy_prm->fit_prms->ag_count,
                                           &energy_prm->fit_prms->prms,
                                           energy_prm->fit_prms->weight);
        json_object_set_new(energy_dict, "density", json_real(term_energy));
        total_energy += term_energy;
    }

    if (gradient != NULL) {
        for (int i = 0; i < array_size / 3; i++) {
            int atom_i = energy_prm->ag->active_atoms->members[i];
            gradient[3 * i] = -energy_prm->ag->gradients[atom_i].X;
            gradient[3 * i + 1] = -energy_prm->ag->gradients[atom_i].Y;
            gradient[3 * i + 2] = -energy_prm->ag->gradients[atom_i].Z;
        }
    }

    json_object_set_new(energy_dict, "total", json_real(total_energy));
    json_array_append_new(energy_prm->json_log, energy_dict);

    return total_energy;
}

