/**
 * \file Functions for parsing variable minimization setup formats
 */

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


/**
 * Regulates the output in json log file
 */
struct json_log_setup {
    bool print_step; ///< Output energy for every step
    bool print_stage; ///< Output final energy for every minimization stage
    bool print_noe_matrix; ///< Matrix can be bulky, so this switches "noe_details" field off in the json log
};


/**
 * The structure sets up different
 * options and terms computed by energy function
 */
struct energy_prms {
    struct mol_atom_group *ag; ///< Minimized atom group
    struct agsetup *ag_setup; ///< Atom group setup
    struct acesetup *ace_setup; ///< GBSA setup
    json_t* json_log; ///< Energy terms are written here every time energy function is envoked
    struct json_log_setup json_log_setup; ///< json log file parameters

    struct pairsprings_setup *sprst_pairs; ///< Pairwise distance restraints
    struct pointsprings_setup *sprst_points; ///< Pointwise distance restraints
    struct noe_setup *nmr; ///< NOE matrix restraints
    struct density_setup *density; ///< Density fitting
    struct fixed_setup *fixed; ///< Fixed atoms

    int nsteps; ///< Max number of steps during minimization
    bool bonds; ///< Flags
    bool angles;
    bool dihedrals;
    bool impropers;
    bool vdw;
    bool vdw03;
    bool gbsa;

    bool score_only; ///< Don't perform minimizition and only output energy terms
};

struct noe_setup {
    double weight;
    double power;
    struct mol_noe *spec;
};

struct density_setup {
    double weight;
    struct mol_atom_group *ag;
    struct mol_fitting_params prms;
};

struct pairspring {
    size_t *group1;      /**< list of first atoms */
    size_t *group2;      /**< list of second atoms */
    size_t leng1;        /**< length of 1st group */
    size_t leng2;        /**< length of 2nd group */
    double weight;      /**< force constant */
    double distance;    /**< distance constraint */
    double lerror;      /**< left error */
    double rerror;      /**< right error */
    char potential[12]; /**< potential for penalty function */
    char average[4];    /**< way to average for group atoms */
};
// write the free function to free the spaces

struct pointspring {
    size_t natoms; ///< Number of atom in a group being pulled
    size_t *atoms; ///< Atom IDs
    double weight; ///< Weight
    double X0, Y0, Z0; ///< Attachment point
};

struct pairsprings_setup {
    size_t nsprings;       /**< number of springs */
    struct pairspring *springs;  /**< array of springs */
};

struct pointsprings_setup {
    size_t nsprings;
    struct pointspring *springs;
};

struct fixed_setup {
    size_t natoms;
    size_t *atoms;
};


/**
 * Fill energy minimization setups from parsed command line arguments
 *
 * @param result_energy_prm Structure being filled in
 * @param result_nstages Number of minimization stages
 * @param opts Parsed command line options
 * @return True on success, false on failure
 */
bool energy_prms_populate_from_options(
        struct energy_prms **result_energy_prm,
        size_t *result_nstages,
        const struct options opts);

void energy_prms_free(struct energy_prms **prms, size_t nstages);


/**
 * Create atom group list from command line arguments. Gradient is already allocated
 *
 * @param opts Parsed command line options
 * @return List of atom groups to minimize
 */
struct mol_atom_group_list* mol_atom_group_list_from_options(struct options *opts);


#endif //ENERGYMIN_SETUP_H
