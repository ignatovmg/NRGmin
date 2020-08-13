/**
 * \file Functions for energy calculation
 */

#ifndef ENERGYMIN_ENERGY_H
#define ENERGYMIN_ENERGY_H

#include "mol2/minimize.h"

#include "setup.h"


/**
 * Master function which computes the enery terms and the gradients
 *
 * @param prm Energy setup "struct energy_prms" from setup.h
 * @param array Internal use
 * @param gradient Atom group gradient to fill
 * @param array_size Internal use
 * @param step Minimization step size
 * @return Energy value
 */
lbfgsfloatval_t energy_func(
        void *restrict prm,
        const double *restrict array,
        double *restrict gradient,
        const int array_size,
        const lbfgsfloatval_t step);


/**
 * Compute pointwise distance restraints
 *
 * @param sprst Setup
 * @param ag Atom group to fill in the gradients
 * @param een Result is added to this variable
 */
void pointspring_energy(const struct pointsprings_setup *sprst, struct mol_atom_group *ag, double *een);


/**
 * Compute pairwise distance restraints
 *
 * @param sprst Setup
 * @param ag Atom group to fill in the gradients
 * @param een Result is added to this variable
 */
void pairspring_energy(const struct pairsprings_setup *sprst, struct mol_atom_group *ag, double *een);


#endif //ENERGYMIN_ENERGY_H
