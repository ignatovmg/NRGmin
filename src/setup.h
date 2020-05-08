#ifndef ENERGYMIN_SETUP_H
#define ENERGYMIN_SETUP_H

#include <stdio.h>
#include <stdlib.h>
#include <jansson.h>

#include "mol2/benergy.h"
#include "mol2/icharmm.h"
#include "mol2/pdb.h"
#include "mol2/fitting.h"
#include "mol2/noe.h"

#include "parse_options.h"


struct json_log_setup {
    bool print_step;
    bool print_stage;
    bool print_noe_matrix;
};

struct energy_prms {
    struct mol_atom_group *ag;
    struct agsetup *ag_setup;
    struct acesetup *ace_setup;
    json_t* json_log;

    struct pairsprings_setup *sprst_pairs;
    struct pointsprings_setup *sprst_points;
    struct noe_setup *nmr;
    struct density_setup *fit_prms;
    struct fixed_setup *fixed;

    int nsteps;
    bool bonds;
    bool angles;
    bool dihedrals;
    bool impropers;
    bool vdw;
    bool vdw03;
    bool gbsa;

    bool score_only;
    bool verbose;

    struct json_log_setup json_log_setup;
};

struct noe_setup {
    double weight;
    double power;
    struct mol_noe *spec;
};

struct density_setup {
    double weight;
    size_t ag_count;
    struct mol_atom_group **ag_list;
    struct mol_fitting_params prms;
};

struct pairspring {
    size_t laspr[2];      /**< list of atoms */
    double lnspr;
    double erspr;
    double fkspr;      /**< force constant */
};

struct pointspring {
    size_t naspr;      /**< number of affected atoms */
    size_t *laspr;      /**< list of atoms */
    double fkspr;      /**< force constant */
    double X0, Y0, Z0; /**< anchor point */
};

struct pairsprings_setup {
    size_t nsprings;       /**< number of springs */
    struct pairspring *springs;  /**< array of springs */
};

struct pointsprings_setup {
    size_t nsprings;       /**< number of springs */
    struct pointspring *springs;  /**< array of springs */
};

struct fixed_setup {
    size_t natoms;
    size_t *atoms;
};


bool energy_prms_populate_from_options(
        struct energy_prms **result_energy_prm,
        size_t *result_nstages,
        const struct options opts);

void energy_prms_free(struct energy_prms **prms, size_t nstages);

struct mol_atom_group_list* mol_atom_group_list_from_options(struct options *opts);


#endif //ENERGYMIN_SETUP_H