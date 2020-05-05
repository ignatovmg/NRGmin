//
// Created by Mikhail Ignatov on 2020-05-02.
//

#ifndef ENERGYMIN_PARSE_OPTIONS_H
#define ENERGYMIN_PARSE_OPTIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <stdbool.h>


struct options {
    char* out_pdb;
    char* out_json;

    bool separate;

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

    int nsteps;

    int bonds;
    int angles;
    int dihedrals;
    int impropers;
    int vdw;
    int vdw03;
    int gbsa;

    int verbose;
    int fix_receptor;
    int score_only;
    int help;
};


struct options parse_args(const int argc, const char** argv, bool *error);

void usage_message(const char **argv);


#endif //ENERGYMIN_PARSE_OPTIONS_H
