//
// Created by Mikhail Ignatov on 2020-05-02.
//

#ifndef ENERGYMIN_PARSE_OPTIONS_H
#define ENERGYMIN_PARSE_OPTIONS_H

#include <stdio.h>
#include <stdbool.h>
#include <jansson.h>


struct options {
    char* out_pdb;
    char* out_json;

    int print_step;
    int print_stage;
    int print_noe_matrix;

    bool separate;
    size_t rec_natoms;
    size_t lig_natoms;

    char* psf;
    char* prm;
    char* rtf;
    char* pdb;
    char* json;

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

    char* fixed_pdb;

    char* pair_springs_txt;
    char* point_springs_txt;

    char* noe_txt;
    char* noe_json;

    char* density_json;

    char* setup_json;
    json_t* setup_json_root;

    int nsteps;

    int bonds;
    int angles;
    int dihedrals;
    int impropers;
    int vdw;
    int vdw03;
    int gbsa;

    int fix_receptor;
    int fix_ligand;
    int score_only;

    int verbosity;

    int help;
};


struct options get_defaut_options(void);

struct options parse_args(const int argc, char *const *argv, bool *error);

void free_options(struct options prms);

void usage_message(char *const *argv);


#endif //ENERGYMIN_PARSE_OPTIONS_H
