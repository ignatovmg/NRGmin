#include <assert.h>

#include "energy.h"

#include "mol2/benergy.h"
#include "mol2/gbsa.h"
#include "mol2/nbenergy.h"
#include "mol2/fitting.h"

#ifdef NOE
#include "mol2/noe.h"
#endif

#include "utils.h"

lbfgsfloatval_t energy_func(
        void *restrict prm,
        const double *restrict array,
        double *restrict gradient,
        const int array_size,
        __attribute__((unused)) const lbfgsfloatval_t step) {

    static const double scale_vdw_s03 = 1.0;
    static const double scale_coul_s03 = 0.8333333;
    static const double eeps = 80.0;

    lbfgsfloatval_t total_energy = 0.0;
    lbfgsfloatval_t term_energy;

    struct energy_prms *energy_prm = (struct energy_prms *) prm;
    json_t* energy_dict = json_object();
    if (energy_prm->json_log) {
        json_array_append_new(energy_prm->json_log, energy_dict);
    }

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

    if (energy_prm->eleng) {
        term_energy = 0.0;
        eleng(energy_prm->ag, eeps, &term_energy, energy_prm->ag_setup->nblst);
        json_object_set_new(energy_dict, "eleng", json_real(term_energy));
        total_energy += term_energy;
    }

    if (energy_prm->elengs03) {
        term_energy = 0.0;
        elengs03(scale_coul_s03,
                energy_prm->ag_setup->nblst->nbcof,
                energy_prm->ag,
                eeps,
                &term_energy,
                energy_prm->ag_setup->nf03,
                energy_prm->ag_setup->listf03);
        json_object_set_new(energy_dict, "elengs03", json_real(term_energy));
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
                scale_vdw_s03,
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

#ifdef NOE
    if (energy_prm->nmr != NULL) {
        mol_noe_calc_peaks(energy_prm->nmr->spec, energy_prm->ag, true, true);
        mol_noe_calc_energy(
                energy_prm->nmr->spec,
                energy_prm->ag->gradients,
                energy_prm->nmr->weight,
                energy_prm->nmr->power);

        term_energy = energy_prm->nmr->spec->energy;
        total_energy += term_energy;

        json_object_set_new(energy_dict, "noe", json_real(term_energy));

        if (energy_prm->json_log_setup.print_noe_matrix) {
            json_object_set_new(energy_dict, "noe_details", mol_noe_to_json_object(energy_prm->nmr->spec));
        }
    }
#endif

    if (energy_prm->density != NULL) {
        term_energy = mol_fitting_score(energy_prm->ag,
                energy_prm->density->ag,
                &energy_prm->density->prms,
                energy_prm->density->weight);
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

    if (!energy_prm->json_log) {
        json_decref(energy_dict);
    }

    return total_energy;
}


void pointspring_energy(const struct pointsprings_setup *sprst, struct mol_atom_group *ag, double *een) {
    for (size_t i = 0; i < sprst->nsprings; i++) {
        size_t nat = sprst->springs[i].natoms;

        if (nat > 0) {
            struct mol_vector3 tot_vec = {0, 0, 0};
            struct mol_vector3 g;

            for (size_t i1 = 0; i1 < nat; i1++) {
                size_t i2 = sprst->springs[i].atoms[i1];
                MOL_VEC_ADD(tot_vec, tot_vec, ag->coords[i2]);
            }

            MOL_VEC_DIV_SCALAR(tot_vec, tot_vec, nat);
            tot_vec.X -= sprst->springs[i].X0;
            tot_vec.Y -= sprst->springs[i].Y0;
            tot_vec.Z -= sprst->springs[i].Z0;

            double fk = sprst->springs[i].weight;
            (*een) += fk * MOL_VEC_SQ_NORM(tot_vec);

            fk = 2 * fk / nat;
            MOL_VEC_MULT_SCALAR(g, tot_vec, fk);

            for (size_t i1 = 0; i1 < nat; i1++) {
                size_t i2 = sprst->springs[i].atoms[i1];
                MOL_VEC_SUB(ag->gradients[i2], ag->gradients[i2], g);
            }
        }
    }
}


void pairspring_energy(const struct pairsprings_setup *sprst, struct mol_atom_group *ag, double *een) {
    size_t i, j, k, idx1, idx2, lni1, lni2, nm;
    double *xtot_a, *ytot_a, *ztot_a, *d_a;
    double xtot, ytot, ztot, d, d2, aved, sumd;
    double delta, coef, ln, ler, rer, fk, hk, gradx, grady, gradz;
    int potential, averaging;
    struct mol_vector3 g;

    for (i = 0; i < sprst -> nsprings; i++) {
        ln = sprst->springs[i].distance;
        ler = sprst->springs[i].lerror;
        rer = sprst->springs[i].rerror;
        fk = sprst->springs[i].weight;

        // start to calculate the average distance over two groups
        lni1 = sprst->springs[i].group_size1;
        lni2 = sprst->springs[i].group_size2;
        xtot_a = calloc(lni1 * lni2, sizeof(double));
        ytot_a = calloc(lni1 * lni2, sizeof(double));
        ztot_a = calloc(lni1 * lni2, sizeof(double));
        d_a = calloc(lni1 * lni2, sizeof(double));

        aved = 0;
        for (j = 0; j < lni1; j++) {
            idx1 = sprst->springs[i].group1[j];
            for (k = 0; k < lni2; k++) {
                idx2 = sprst->springs[i].group2[k];
                xtot = ag->coords[idx2].X - ag->coords[idx1].X;
                ytot = ag->coords[idx2].Y - ag->coords[idx1].Y;
                ztot = ag->coords[idx2].Z - ag->coords[idx1].Z;
                d2 = xtot * xtot + ytot * ytot + ztot * ztot;
                d = sqrt(d2);
                aved += pow(d, -6.0);
                xtot_a[lni2 * j + k] = xtot;
                ytot_a[lni2 * j + k] = ytot;
                ztot_a[lni2 * j + k] = ztot;
                d_a[lni2 * j + k] = d;
            }
        }
        sumd = aved;
        averaging = sprst->springs[i].average;
        if (averaging == 0) {
            aved = pow(aved, -1.0 / 6.0);
            nm = 1;
        } else if (averaging == 1) {  // R6 average
            aved = pow(aved / lni1 / lni2, -1.0 / 6);
            nm = lni1 * lni2;
        }
        delta = aved - ln;

        potential = sprst->springs[i].potential;
        if (potential == 0) {  // square-well
            if (delta < 0) {
                delta = (delta < -ler) ? (delta + ler) : 0.0;
            } else {
                delta = (delta > rer) ? (delta - rer) : 0.0;
            }
            (*een) += fk * delta * delta;
            coef = fk * 2.0 * delta * pow(nm, 1.0 / 6.0) * pow(sumd, -7.0 / 6.0);
        } else if (potential == 1) {  // biharmonic, temperature=300
            if (delta < 0) {
                hk = fk * 300 * 0.0019872041 / (2 * ler * ler);
            } else {
                hk = fk * 300 * 0.0019872041 / (2 * rer * rer);
            }
            hk = fmin(1000, hk);
            (*een) += hk * delta * delta;
            coef = hk * 2.0 * delta * pow(nm, 1.0 / 6.0) * pow(sumd, -7.0 / 6.0);
        } else if (potential == 2) {  // soft-square
            if (delta < 0) {
                delta = (delta < -ler) ? (delta + ler) : 0.0;
            } else {
                delta = (delta > rer) ? (delta - rer) : 0.0;
            }

            if (delta <= 3.0) {  // here the switch bound is 3.0
                (*een) += fk * delta * delta;
                coef = fk * 2.0 * delta * pow(nm, 1.0 / 6.0) * pow(sumd, -7.0 / 6.0);
            } else {  //delta > 3.0
                (*een) += fk * (delta - 45 / delta + 21);
                coef = fk * (1 + 45 / (delta * delta)) * pow(nm, 1.0 / 6) * pow(sumd, -7.0 / 6);
            }
        }
        // calculate and update gradients first groups
        for (j = 0; j < lni1; j++) {
            idx1 = sprst->springs[i].group1[j];
            gradx = 0;
            grady = 0;
            gradz = 0;
            for (k = 0; k < lni2; k++) {
                gradx += coef * xtot_a[j * lni2 + k] * pow(d_a[j * lni2 + k], -8.0);
                grady += coef * ytot_a[j * lni2 + k] * pow(d_a[j * lni2 + k], -8.0);
                gradz += coef * ztot_a[j * lni2 + k] * pow(d_a[j * lni2 + k], -8.0);
            }
            g.X = -gradx;
            g.Y = -grady;
            g.Z = -gradz;
            MOL_VEC_SUB(ag->gradients[idx1], ag->gradients[idx1], g);
        }
        // calculate and update gradients second groups
        for (k = 0; k < lni2; k++) {
            idx2 = sprst->springs[i].group2[k];
            gradx = 0;
            grady = 0;
            gradz = 0;
            for (j = 0; j < lni1; j++) {
                gradx += coef * xtot_a[j * lni2 + k] * pow(d_a[j * lni2 + k], -8.0);
                grady += coef * ytot_a[j * lni2 + k] * pow(d_a[j * lni2 + k], -8.0);
                gradz += coef * ztot_a[j * lni2 + k] * pow(d_a[j * lni2 + k], -8.0);
            }
            g.X = -gradx;
            g.Y = -grady;
            g.Z = -gradz;
            MOL_VEC_ADD(ag->gradients[idx2], ag->gradients[idx2], g);
        }
        free(xtot_a);
        free(ytot_a);
        free(ztot_a);
        free(d_a);
    }
}