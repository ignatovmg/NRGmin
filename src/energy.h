#ifndef ENERGYMIN_ENERGY_H
#define ENERGYMIN_ENERGY_H

#include "mol2/minimize.h"

#include "setup.h"


lbfgsfloatval_t energy_func(
        void *restrict prm,
        const double *restrict array,
        double *restrict gradient,
        const int array_size,
        const lbfgsfloatval_t step);

void pointspring_energy(const struct pointsprings_setup *sprst, struct mol_atom_group *ag, double *een);

void pairspring_energy(const struct pairsprings_setup *sprst, struct mol_atom_group *ag, double *een);


#endif //ENERGYMIN_ENERGY_H
