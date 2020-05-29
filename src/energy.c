#include "energy.h"

#include "mol2/benergy.h"
#include "mol2/gbsa.h"
#include "mol2/nbenergy.h"
#include "mol2/fitting.h"
#include "mol2/noe.h"


lbfgsfloatval_t energy_func(
        void *restrict prm,
        const double *restrict array,
        double *restrict gradient,
        const int array_size,
        __attribute__((unused)) const lbfgsfloatval_t step) {

    lbfgsfloatval_t total_energy = 0.0;
    lbfgsfloatval_t term_energy;

    struct energy_prms *energy_prm = (struct energy_prms *) prm;
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

    double ele_eps = 1.0;
    if (energy_prm->eleng) {
        term_energy = 0.0;
        eleng(energy_prm->ag, ele_eps, &term_energy, energy_prm->ag_setup->nblst);
        json_object_set_new(energy_dict, "eleng", json_real(term_energy));
        total_energy += term_energy;
    }

    if (energy_prm->elengs03) {
        term_energy = 0.0;
        elengs03(1.0,
                energy_prm->ag_setup->nblst->nbcof,
                energy_prm->ag,
                ele_eps,
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

        if (energy_prm->json_log_setup.print_noe_matrix) {
            json_object_set_new(energy_dict, "noe_details", mol_noe_to_json_object(energy_prm->nmr->spec));
        }
    }

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

    if (energy_prm->json_log) {
        json_array_append_new(energy_prm->json_log, energy_dict);
    } else {
        json_decref(energy_dict);
    }

    return total_energy;
}


void pointspring_energy(const struct pointsprings_setup *sprst, struct mol_atom_group *ag, double *een) {
    size_t i, i1, i2, nat;
    double xtot, ytot, ztot, fk;
    struct mol_vector3 g;

    for (i = 0; i < sprst->nsprings; i++) {
        nat = sprst->springs[i].natoms;

        if (nat > 0) {
            xtot = 0.0;
            ytot = 0.0;
            ztot = 0.0;

            for (i1 = 0; i1 < nat; i1++) {
                i2 = sprst->springs[i].atoms[i1];
                xtot += ag->coords[i2].X;
                ytot += ag->coords[i2].Y;
                ztot += ag->coords[i2].Z;
            }

            xtot = xtot / nat - sprst->springs[i].X0;
            ytot = ytot / nat - sprst->springs[i].Y0;
            ztot = ztot / nat - sprst->springs[i].Z0;

            fk = sprst->springs[i].weight;
            (*een) += fk * (xtot * xtot + ytot * ytot + ztot * ztot);

            fk = 2 * fk / nat;
            g.X = xtot * fk;
            g.Y = ytot * fk;
            g.Z = ztot * fk;

            for (i1 = 0; i1 < nat; i1++) {
                i2 = sprst->springs[i].atoms[i1];
                MOL_VEC_SUB(ag->gradients[i2], ag->gradients[i2], g);
            }
        }
    }
}


void pairspring_energy(const struct pairsprings_setup *sprst, struct mol_atom_group *ag, double *een) {
    size_t i, i1, i2;
    double xtot, ytot, ztot, fk, d, d2, ln, er, coef, delta;
    struct mol_vector3 g;

    for (i = 0; i < sprst->nsprings; i++) {
        ln = sprst->springs[i].length;
        er = sprst->springs[i].error;
        fk = sprst->springs[i].weight / 2.0;

        i1 = sprst->springs[i].atoms[0];
        i2 = sprst->springs[i].atoms[1];

        xtot = ag->coords[i2].X - ag->coords[i1].X;
        ytot = ag->coords[i2].Y - ag->coords[i1].Y;
        ztot = ag->coords[i2].Z - ag->coords[i1].Z;

        d2 = xtot * xtot + ytot * ytot + ztot * ztot;
        d = sqrt(d2);

        delta = fabs(d - ln);
        delta = (delta > er) ? ((delta - er) * delta / (d - ln)) : 0.0;

        //(*een) += fk * (d - ln) * (d - ln);
        (*een) += fk * delta * delta;
        //coef = fk * 2 * (1.0 - ln / d);
        coef = fk * 2.0 * delta / d;

        g.X = -coef * xtot;
        g.Y = -coef * ytot;
        g.Z = -coef * ztot;

        MOL_VEC_SUB(ag->gradients[i1], ag->gradients[i1], g);
        MOL_VEC_ADD(ag->gradients[i2], ag->gradients[i2], g);
    }
}