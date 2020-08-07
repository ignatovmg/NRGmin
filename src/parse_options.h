/**
 * \file Parsing command line arguments
 */

#ifndef ENERGYMIN_PARSE_OPTIONS_H
#define ENERGYMIN_PARSE_OPTIONS_H

#include <stdio.h>
#include <stdbool.h>
#include <jansson.h>


/**
 * Parsed command line options
 *
 * Structure can be passed using two modes: single file mode or rec/lig mode
 * Rec/lig mode is used when structure is passed using rec-* lig-* arguments, and
 * single file mode is enabled for the opposite.
 *
 * All flags in this structure can be filled by passing a file to --setup-json,
 * so this one argument is enough to specify all minimization parameters.
 */
struct options {
    char* out_pdb; ///< Where to write the minimized molecule(s)
    char* out_json; ///< Log file with energy terms

    int print_step; ///< Log energies for every step
    int print_stage; ///< Log energies for every stage
    int print_noe_matrix; ///< Log NOE matrix is NOE calculation is on

    bool separate; ///< Receptor and ligand provided separately
    size_t rec_natoms; ///< Number of atoms in receptor (set to zero in single file mode)
    size_t lig_natoms; ///< Number of atoms in ligand (set to zero in single file mode)

    ///< Parameters for single file mode
    char* psf;
    char* prm;
    char* rtf;
    char* pdb;
    char* json;

    ///< Parameters for rec/lig mode
    char* rec_psf;
    char* rec_prm;
    char* rec_rtf;
    char* rec_pdb;
    char* rec_json;

    char* lig_psf;
    char* lig_prm;
    char* lig_rtf;
    char* lig_pdb;
    char* lig_json;

    char* fixed_pdb; ///< PDB file with fixed atoms

    char* pair_springs_txt; ///< Text file setup for pairwise restraints
    char* point_springs_txt; ///< Text file setup for pointwise restraints

    char* noe_txt; ///< Text file setup for NOE
    char* noe_json; ///< Json file setup for NOE

    char* density_json; ///< Json file setup for density fitting

    char* setup_json; ///< This file duplicates all the other options present in this structure. It also allows for
///< multistage minimization, where each stage can have different minimization terms
    json_t* setup_json_root; ///< Contents of setup_json. This has to persist, if setup_json was used, othewise
///< all char* options in this structure will be erased

    int nsteps; ///< Number of minimization steps

    int bonds; ///< 1 - on, 0 - off
    int angles;
    int dihedrals;
    int impropers;
    int vdw;
    int vdw03;
    int eleng;
    int elengs03;
    int ace;
    int gbsa;

    int fix_receptor; ///< Works only rec/lig mode
    int fix_ligand; ///< Works only rec/lig mode
    int score_only; ///< Score without doing minimization

    int verbosity; ///< Levels 0 - QUIET, 1 - ERROR, 2 - WARNING, 3- INFO, 4 - DEBUG

    int num_threads;
    int help; ///< Print help and exit
};


struct options options_get_default(void);

/**
 * Fill options structure from command line args
 * @param argc Main function argc
 * @param argv Main function argv
 * @param error True if error occured
 * @return Parsed options
 */
struct options options_populate_from_argv(const int argc, char *const *argv, bool *error);

void options_free(struct options opts);

void usage_message(char *const *argv);


#endif //ENERGYMIN_PARSE_OPTIONS_H
