#include <string.h>
#include <stdlib.h>
#include <getopt.h>

#include "parse_options.h"
#include "utils.h"


/**
 * Stringize "arg" find the dot, replace "-" with "_" and compare to "value"
 */
#define _FILL_PARAM(arg, value) do {                                   \
    char str_arg[] = #arg;                                             \
    char* dot_pos;                                                     \
    if ((dot_pos = strchr(str_arg, '.')) != NULL) {                    \
        dot_pos = dot_pos + 1;                                         \
    }                                                                  \
    char* dash_pos = dot_pos;                                          \
    while ((dash_pos = strchr(dash_pos, '_')) != NULL) {               \
        *dash_pos = '-';                                               \
    };                                                                 \
    if (strcmp(dot_pos, long_options[option_index].name) == 0) {       \
            (arg) = (value);                                           \
            break;                                                     \
    }                                                                  \
} while (0)


struct options options_get_default()
{
    struct options prms = {
            .out_pdb = "out.pdb",
            .out_json = "out.json",

            .print_step = false,
            .print_stage = true,
            .print_noe_matrix = true,

            .separate = false,
            .rec_natoms = 0,
            .lig_natoms = 0,

            .psf = NULL,
            .prm = NULL,
            .rtf = NULL,
            .pdb = NULL,
            .json = NULL,

            .rec_psf = NULL,
            .rec_prm = NULL,
            .rec_rtf = NULL,
            .rec_pdb = NULL,
            .rec_json = NULL,

            .lig_psf = NULL,
            .lig_prm = NULL,
            .lig_rtf = NULL,
            .lig_pdb = NULL,
            .lig_json = NULL,

            .fixed_pdb = NULL,

            .pair_springs_txt = NULL,
            .point_springs_txt = NULL,

            .noe_txt = NULL,
            .noe_json = NULL,

            .density_json = NULL,

            .setup_json = NULL,
            .setup_json_root = NULL,

            .nsteps = 1000,

            .bonds = 1,
            .angles = 1,
            .dihedrals = 1,
            .impropers = 1,
            .vdw = 1,
            .vdw03 = 1,
            .eleng = 1,
            .elengs03 = 1,
            .ace = 0,
            .gbsa = 0,

            .verbosity = DEBUG,
            .fix_receptor = 0,
            .fix_ligand = 0,
            .score_only = 0,
            .help = 0,

            .num_threads = 0
    };

    return prms;
}


void options_free(struct options opts)
{
    if (opts.setup_json_root) {
        json_decref(opts.setup_json_root);
    }
}


static bool _fill_prms_from_json(struct options* opts, const json_t* root)
{
    if (!json_is_object(root)) {
        return false;
    }

    json_t* dict = json_object_get(root, "options");
    if (!dict) {
        return true;
    }

    if (!json_is_object(dict)) {
        ERR_MSG("Key 'options' must point to a dictionary");
        return false;
    }

    json_error_t error;

    int code = json_unpack_ex(
            dict, &error, JSON_STRICT,

            "{s?i, s?i, s?i "
            " s?b, s?b, s?b, s?b, s?b, s?b, s?b, s?b, s?b, s?b, s?b, s?b, s?b, s?b, s?b, s?b, "
            " s?s, s?s, s?s, s?s, s?s, s?s, s?s, s?s, s?s, s?s, s?s, s?s, s?s, s?s, "
            " s?s, s?s, s?s, s?s, s?s, s?s, s?s, s?s, s?s}",

            // integer options
            "nsteps", &opts->nsteps,
            "verbosity", &opts->verbosity,
            "num_threads", &opts->num_threads,

            // binary options
            "bonds", &opts->bonds,
            "angles", &opts->angles,
            "dihedrals", &opts->dihedrals,
            "impropers", &opts->impropers,
            "vdw", &opts->vdw,
            "vdw03", &opts->vdw03,
            "eleng", &opts->eleng,
            "elengs03", &opts->elengs03,
            "ace", &opts->ace,
            "gbsa", &opts->gbsa,
            "fix_receptor", &opts->fix_receptor,
            "fix_ligand", &opts->fix_ligand,
            "score_only", &opts->score_only,
            "print_step", &opts->print_step,
            "print_stage", &opts->print_stage,
            "print_noe_matrix", &opts->print_noe_matrix,

            // char options
            "out_pdb", &opts->out_pdb,
            "out_json", &opts->out_json,
            "psf", &opts->psf,
            "prm", &opts->prm,
            "rtf", &opts->rtf,
            "pdb", &opts->pdb,
            "json", &opts->json,
            "rec_psf", &opts->rec_psf,
            "rec_prm", &opts->rec_prm,
            "rec_rtf", &opts->rec_rtf,
            "rec_pdb", &opts->rec_pdb,
            "rec_json", &opts->rec_json,
            "lig_psf", &opts->lig_psf,
            "lig_prm", &opts->lig_prm,
            "lig_rtf", &opts->lig_rtf,
            "lig_pdb", &opts->lig_pdb,
            "lig_json", &opts->lig_json,
            "fixed_pdb", &opts->fixed_pdb,
            "pair_springs_txt", &opts->pair_springs_txt,
            "point_springs_txt", &opts->point_springs_txt,
            "noe_txt", &opts->noe_txt,
            "noe_json", &opts->noe_json,
            "density_json", &opts->density_json
    );

    if (code != 0) {
        JSON_ERR_MSG(error, "Couldn't parse setup from json");
        return false;
    }

    return true;
}


static int _check_prms(struct options *opts)
{
    bool full = false;
    bool rec_sep = false;
    bool lig_sep = false;

    VERBOSITY = opts->verbosity >= 0 ? opts->verbosity : 0;

    // Unpack setup json first and fill the global options if there are any
    if (opts->setup_json) {
        json_t* setup = read_json_file(opts->setup_json);
        if (!setup) {
            return 1;
        }

        if (!_fill_prms_from_json(opts, setup)) {
            json_decref(setup);
            return 1;
        }
        opts->setup_json_root = setup;

        VERBOSITY = opts->verbosity >= 0 ? opts->verbosity : 0;
    }

    if (opts->num_threads < 0) {
        ERR_MSG("Number of threads can't be negative");
        return 1;
    }

    if (opts->gbsa && opts->ace) {
        WRN_MSG("You're using GBSA and ACE at the same time");
    }

    if (opts->json || opts->pdb || opts->psf || opts->prm || opts->rtf) {
        full = true;
        opts->separate = false;

        if (!(opts->json || (opts->pdb && opts->psf && opts->prm && opts->rtf))) {
            ERR_MSG("Json and/or pdb, psf, prm and rtf must be provided");
            return 1;
        }

        if (opts->fix_receptor) {
            ERR_MSG("Can't fix the receptor when using single file mode");
            return 1;
        }

        if (opts->fix_ligand) {
            ERR_MSG("Can't fix the ligand when using single file mode");
            return 1;
        }
    }

    if (opts->rec_pdb || opts->rec_psf || opts->rec_prm || opts->rec_rtf || opts->rec_json) {
        if (full) {
            ERR_MSG("You can't provide options for rec/lig mode and single file mode at the same time");
            return 1;
        }

        if (!(opts->rec_json || (opts->rec_pdb && opts->rec_psf && opts->rec_prm && opts->rec_rtf))) {
            ERR_MSG("Json and/or pdb, psf, prm and rtf must be provided for receptor");
            return 1;
        }

        rec_sep = true;

        if (!(opts->lig_json || (opts->lig_pdb && opts->lig_psf && opts->lig_prm && opts->lig_rtf))) {
            ERR_MSG("Json and/or pdb, psf, prm and rtf must be provided for ligand");
            return 1;
        }

        lig_sep = true;
        opts->separate = true;
    }

    if (!(full || (rec_sep && lig_sep))) {
        ERR_MSG("The molecule wasn't provided");
        return 1;
    }

    if (opts->nsteps < 0) {
        ERR_MSG("Number of steps can't be negative");
        return 1;
    }

#ifndef NOE
    if (opts->noe_txt || opts->noe_json) {
        ERR_MSG("Rebuild NRGmin with -DNOE=ON to use NOE spectrum fitting");
        return 1;
    }
#endif

    return 0;
}


struct options options_populate_from_argv(const int argc, char *const *argv, bool *error)
{
    struct options prms = options_get_default();

    struct option long_options[] =
            {
                    {"out-pdb",              required_argument, 0,                 0},
                    {"out-json",             required_argument, 0,                 0},
                    {"verbosity",            required_argument, 0,                 0},
                    {"num-threads",          required_argument, 0,                 0},

                    {"print-step",           no_argument,       &prms.print_step,  1},
                    {"print-stage",          no_argument,       &prms.print_stage, 1},
                    {"print-noe-matrix",     no_argument,  &prms.print_noe_matrix, 1},

                    {"psf",                  required_argument, 0,                 0},
                    {"prm",                  required_argument, 0,                 0},
                    {"rtf",                  required_argument, 0,                 0},
                    {"pdb",                  required_argument, 0,                 0},
                    {"json",                 required_argument, 0,                 0},

                    {"rec-psf",              required_argument, 0,                 0},
                    {"rec-prm",              required_argument, 0,                 0},
                    {"rec-rtf",              required_argument, 0,                 0},
                    {"rec-pdb",              required_argument, 0,                 0},
                    {"rec-json",             required_argument, 0,                 0},

                    {"lig-psf",              required_argument, 0,                 0},
                    {"lig-prm",              required_argument, 0,                 0},
                    {"lig-rtf",              required_argument, 0,                 0},
                    {"lig-pdb",              required_argument, 0,                 0},
                    {"lig-json",             required_argument, 0,                 0},

                    {"fixed-pdb",            required_argument, 0,                 0},

                    {"pair-springs-txt",     required_argument, 0,                 0},
                    {"point-springs-txt",    required_argument, 0,                 0},

                    {"noe-txt",              required_argument, 0,                 0},
                    {"noe-json",             required_argument, 0,                 0},

                    {"density-json",         required_argument, 0,                 0},

                    {"setup-json",           required_argument, 0,                 0},

                    {"nsteps",               required_argument, 0,                 0},

                    {"bonds-on",             no_argument,  &prms.bonds,            1},
                    {"bonds-off",            no_argument,  &prms.bonds,            0},
                    {"angles-on",            no_argument,  &prms.angles,           1},
                    {"angles-off",           no_argument,  &prms.angles,           0},
                    {"dihedrals-on",         no_argument,  &prms.dihedrals,        1},
                    {"dihedrals-off",        no_argument,  &prms.dihedrals,        0},
                    {"impropers-on",         no_argument,  &prms.impropers,        1},
                    {"impropers-off",        no_argument,  &prms.impropers,        0},
                    {"vdw-on",               no_argument,  &prms.vdw,              1},
                    {"vdw-off",              no_argument,  &prms.vdw,              0},
                    {"vdw03-on",             no_argument,  &prms.vdw03,            1},
                    {"vdw03-off",            no_argument,  &prms.vdw03,            0},
                    {"eleng-on",             no_argument,  &prms.eleng,            1},
                    {"eleng-off",            no_argument,  &prms.eleng,            0},
                    {"elengs03-on",          no_argument,  &prms.elengs03,         1},
                    {"elengs03-off",         no_argument,  &prms.elengs03,         0},
                    {"gbsa-on",              no_argument,  &prms.gbsa,             1},
                    {"gbsa-off",             no_argument,  &prms.gbsa,             0},
                    {"ace-on",               no_argument,  &prms.ace,              1},
                    {"ace-off",              no_argument,  &prms.ace,              0},

                    {"fix-receptor",         no_argument,  &prms.fix_receptor,     1},
                    {"fix-ligand",           no_argument,  &prms.fix_ligand,       1},
                    {"score-only",           no_argument,  &prms.score_only,       1},
                    {"help",                 no_argument,  0,                    'h'},

                    {0,                      0,            0,                      0}
            };

    int option_index = 0;
    int opt;
    while (1) {
        opt = getopt_long_only(argc, argv, "h", long_options, &option_index);

        if (opt == -1) {
            break;
        }

        switch (opt) {
            case 0:
                if (long_options[option_index].flag != 0)
                    break;
                //DEBUG_MSG("Option %s = %s", long_options[option_index].name, optarg);
                break;

            case 'h':
                *error = false;
                prms.help = true;
                return prms;

            case '?':
                *error = true;
                return prms;

            default:
                *error = true;
                return prms;
        }

        _FILL_PARAM(prms.out_pdb, optarg);
        _FILL_PARAM(prms.out_json, optarg);
        _FILL_PARAM(prms.psf, optarg);
        _FILL_PARAM(prms.prm, optarg);
        _FILL_PARAM(prms.rtf, optarg);
        _FILL_PARAM(prms.pdb, optarg);
        _FILL_PARAM(prms.json, optarg);
        _FILL_PARAM(prms.rec_psf, optarg);
        _FILL_PARAM(prms.rec_prm, optarg);
        _FILL_PARAM(prms.rec_rtf, optarg);
        _FILL_PARAM(prms.rec_pdb, optarg);
        _FILL_PARAM(prms.rec_json, optarg);
        _FILL_PARAM(prms.lig_psf, optarg);
        _FILL_PARAM(prms.lig_prm, optarg);
        _FILL_PARAM(prms.lig_rtf, optarg);
        _FILL_PARAM(prms.lig_pdb, optarg);
        _FILL_PARAM(prms.lig_json, optarg);
        _FILL_PARAM(prms.fixed_pdb, optarg);
        _FILL_PARAM(prms.pair_springs_txt, optarg);
        _FILL_PARAM(prms.point_springs_txt, optarg);
        _FILL_PARAM(prms.noe_txt, optarg);
        _FILL_PARAM(prms.noe_json, optarg);
        _FILL_PARAM(prms.density_json, optarg);
        _FILL_PARAM(prms.setup_json, optarg);

        if (strcmp("nsteps", long_options[option_index].name) == 0) {
            prms.nsteps = atoi(optarg);
        }
        if (strcmp("verbosity", long_options[option_index].name) == 0) {
            prms.verbosity = atoi(optarg);
        }
        if (strcmp("num-threads", long_options[option_index].name) == 0) {
            prms.num_threads = atoi(optarg);
        }
    }

    if (_check_prms(&prms) != 0) {
        *error = true;
        return prms;
    }

    *error = false;
    return prms;
}


void usage_message(char *const *argv) {
    fprintf(stderr, "\nUsage %s [ options ]\n\n", argv[0]);
    fprintf(stderr, "Energy minimization using LBFGS.\n"
           "\n"
           "A molecule can be passed using two modes: single file mode or rec/lig mode.\n"
           "Rec/lig mode is used, when structure is passed using rec-* lig-* arguments,\n"
           "and single file mode is enabled for the opposite. Structure can be passed\n"
           "using a set of 4 files pdb, psf, rtf, prm or as a single json file which\n"
           "contains coordinates as well as force field parameters. Also and pair of\n"
           "json and pdb can be passed, if one want to use geometry from json and\n"
           "coordinates from pdb.\n"
           "\n"
           "All flags can be as well filled by passing a global setup file to --setup-json,\n"
           "so this one argument is enough to run minimization.\n"
           "\n"
           "Examples of different minimization setups can be found in examples/energy_setup.\n"
           "\n"
           "    Output control:\n"
           "\n"
           "    --out-pdb Where to write the minimized molecule(s)\n"
           "    --out-json Log file with energy terms\n"
           "    --print-step Log energies for every step\n"
           "    --print-stage Log energies for every stage\n"
           "    --print-noe-matrix Log NOE matrix is NOE calculation is on\n"
           "    --verbosity 0 - QUIET, 1 - ERROR, 2 - WARNING, 3 - INFO, >=4 - DEBUG\n"
           "\n"
           "\n"
           "    Single file mode:\n"
           "\n"
           "    --psf Geometry\n"
           "    --prm Forcefield parameters\n"
           "    --rtf Topology file with residues description\n"
           "    --pdb Atom coordinates (can contain multiple models)\n"
           "    --json Json file with coordinates, geometry and force field parameters\n"
           "\n"
           "    Parameters for rec/lig mode\n"
           "\n"
           "    --rec-psf Receptor geometry\n"
           "    --rec-prm Receptor forcefield parameters\n"
           "    --rec-rtf Receptor topology file\n"
           "    --rec-pdb Receptor atom coordinates (can contain multiple models)\n"
           "    --rec-json Receptor json file with coordinates, geometry and force field parameters\n"
           "    --lig-psf Ligand geometry\n"
           "    --lig-prm Ligand forcefield parameters\n"
           "    --lig-rtf Ligand topology file with residues description\n"
           "    --lig-pdb Ligand atom coordinates (can contain multiple models)\n"
           "    --lig-json Ligand json file with coordinates, geometry and force field parameters\n"
           "\n"
           "\n"
           "    Minimization setup\n"
           "\n"
           "    --nsteps Number of minimization steps\n"
           "\n"
           "    Energy terms switches. Everything is on by default except for GBSA and ACE\n"
           "\n"
           "    --bonds-on/--bonds-off Bonds energy term\n"
           "    --angles-on/--angles-off Angles\n"
           "    --dihedrals-on/--dihedrals-off Dihedrals\n"
           "    --impropers-on/--impropers-on Impropers\n"
           "    --vdw-on/--vdw-off VDW\n"
           "    --vdw03-on/--vdw03-off 1-4 VDW\n"
           "    --eleng-on/--eleng-off Coulomb electrostatics\n"
           "    --elengs03-on/--elengs03-off 1-4 Coulomb electrostatics\n"
           "    --gbsa-on/gbsa-off GBSA\n"
           "    --ace-on/ace-off ACE\n"
           "\n"
           "    Fixed atoms (in rec/lig mode ligand atom IDs must be increased by the number of receptor atoms)\n"
           "\n"
           "    --fixed-pdb PDB file with fixed atoms. Only atom ID is read (lines starting with ATOM)\n"
           "    --fix-receptor Fix receptor atoms (works only in rec/lig mode)\n"
           "    --fix-ligand Fix ligand atoms (works only in rec/lig mode)\n"
           "\n"
           "    Distance restraints\n"
           "\n"
           "    --pair-springs-txt Text file setup for pairwise restraints\n"
           "    --point-springs-txt Text file setup for pointwise restraints\n"
           "\n"
#ifdef NOE
           "    NOE setup\n"
           "\n"
           "    --noe-txt Text file setup for NOE\n"
           "    --noe-json Json file setup for NOE\n"
           "\n"
#endif
           "    Density setup\n"
           "\n"
           "    --density-json Json file setup for density fitting\n"
           "\n"
           "    Global setup in json format\n"
           "\n"
           "    --setup-json This file duplicates all the other options. It also allows for\n"
           "                 multistage minimization protocol, where each stage can have\n"
           "                 different minimization terms\n"
           "\n"
           "    Miscellaneous\n"
           "\n"
           "    --num-threads Number of OpenMP threads to use. 0 (default) means using maximum number of threads\n"
           "    --score-only Score only without doing minimization\n"
           "    --help Show this message and exit\n\n");
}
