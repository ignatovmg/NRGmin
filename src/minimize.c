#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <jansson.h>

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

    struct options parsed = parse_args(argc, argv);

    struct mol_atom_group_list* ag_list = mol_atom_group_list_from_options(parsed);
    if (!ag_list) {
        ERR_MSG("Couldn't read atom groups");
        exit(EXIT_FAILURE);
    }

    struct energy_prm* min_prms;
    size_t nstages;
    if (!energy_prm_read(&min_prms, &nstages, parsed, ag_list)) {
        ERR_MSG("Couldn't fill params");
        mol_atom_group_list_free(ag_list);
        exit(EXIT_FAILURE);
    }

    INFO_MSG("Started the main loop\n");
    FILE* out_pdb = fopen(parsed.out_pdb, "w");
    if (!out_pdb) {
        ERR_MSG("Can't open %s", parsed.out_pdb);
        mol_atom_group_list_free(ag_list);
        exit(EXIT_FAILURE);
    }

    for (int modeli = 0; modeli < ag_list->size; modeli++) {
        INFO_MSG("Started model %i\n", modeli);

        if (ag_list->size > 1) {
            if (out_pdb != NULL) {
                fprintf(out_pdb, "MODEL %i\n", (modeli + 1));
            }
        }

        struct mol_atom_group *ag = &ag_list->members[modeli];
        ag->gradients = calloc(ag->natoms, sizeof(struct mol_vector3));

        for (size_t stage_id; stage_id < nstages; stage_id++) {
            struct energy_prm *stage_prms = &min_prms[stage_id];

            stage_prms->ag = ag;

            // Setup fixed atoms and lists
            mol_fixed_init(ag);

            if (stage_prms->fixed) {
                mol_fixed_update(ag, stage_prms->fixed->natoms, stage_prms->fixed->atoms);
            } else {
                mol_fixed_update(ag, 0, NULL);
            }

            init_nblst(ag, stage_prms->ag_setup);
            update_nblst(ag, stage_prms->ag_setup);

            if (stage_prms->gbsa) {
                stage_prms->ace_setup->efac = 0.5;
                ace_ini(ag, stage_prms->ace_setup);
                ace_fixedupdate(ag, stage_prms->ag_setup, stage_prms->ace_setup);
                ace_updatenblst(stage_prms->ag_setup, stage_prms->ace_setup);
            }

            mol_minimize_ag(MOL_LBFGS, stage_prms->nsteps, __TOL__, ag, (void *) (&stage_prms), energy_func);
        }
    }

    FILE* out_json = fopen(parsed.out_json, "w");
    if (!out_json) {
        ERR_MSG("Can't open %s", parsed.out_json);
        mol_atom_group_list_free(ag_list);
        exit(EXIT_FAILURE);
    }

    INFO_MSG("Completed\n");
    return EXIT_SUCCESS;
}


static lbfgsfloatval_t energy_func(
        void *restrict prm,
        const double *restrict array,
        double *restrict gradient,
        const int array_size,
        const lbfgsfloatval_t step) {
    lbfgsfloatval_t energy = 0.0;
    struct energy_prm *energy_prm = (struct energy_prm *) prm;

    if (array != NULL) {
        assert(array_size == energy_prm->ag->active_atoms->size * 3);
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
        aceeng(energy_prm->ag, &energy, energy_prm->ace_setup, energy_prm->ag_setup);
    }

    if (energy_prm->vdw) {
        vdweng(energy_prm->ag, &energy, energy_prm->ag_setup->nblst);
    }

    if (energy_prm->vdw03) {
        vdwengs03(
                1.0,
                energy_prm->ag_setup->nblst->nbcof,
                energy_prm->ag, &energy,
                energy_prm->ag_setup->nf03,
                energy_prm->ag_setup->listf03);
    }

    if (energy_prm->bonds) {
        beng(energy_prm->ag, &energy);
    }

    if (energy_prm->angles) {
        aeng(energy_prm->ag, &energy);
    }

    if (energy_prm->dihedrals) {
        teng(energy_prm->ag, &energy);
    }

    if (energy_prm->impropers) {
        ieng(energy_prm->ag, &energy);
    }

    if (energy_prm->sprst_pairs != NULL) {
        pairspring_energy(energy_prm->sprst_pairs, energy_prm->ag, &energy);
    }

    if (energy_prm->sprst_points != NULL) {
        pointspring_energy(energy_prm->sprst_points, energy_prm->ag, &energy);
    }

    if (energy_prm->nmr != NULL) {
        mol_noe_calc_peaks(energy_prm->nmr->spec, energy_prm->ag, true);
        mol_noe_calc_energy(energy_prm->nmr->spec, energy_prm->ag->gradients, energy_prm->nmr->weight, 1. / 6.);
        energy += energy_prm->nmr->spec->energy;
    }

    if (energy_prm->fit_prms != NULL) {
        energy += mol_fitting_score_aglist(energy_prm->ag,
                                           energy_prm->fit_prms->ag_list,
                                           energy_prm->fit_prms->ag_count,
                                           &energy_prm->fit_prms->prms,
                                           energy_prm->fit_prms->weight);
    }

    if (gradient != NULL) {
        for (int i = 0; i < array_size / 3; i++) {
            int atom_i = energy_prm->ag->active_atoms->members[i];
            gradient[3 * i] = -energy_prm->ag->gradients[atom_i].X;
            gradient[3 * i + 1] = -energy_prm->ag->gradients[atom_i].Y;
            gradient[3 * i + 2] = -energy_prm->ag->gradients[atom_i].Z;
        }
    }

    return energy;
}

