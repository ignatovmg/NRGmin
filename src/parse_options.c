#include <string.h>
#include <stdlib.h>
#include <getopt.h>

#include "parse_options.h"
#include "utils.h"


#define _STRINGIZE_AND_REPLACE_UNDERSCORE(arg, sep, dst) do { \
    char str_arg[] = #arg;  \
    char* dot_pos; \
    if ((dot_pos = strchr(str_arg, (sep))) != NULL) { \
        dot_pos = dot_pos + 1; \
    }  \
    char* dash_pos = dot_pos; \
    while ((dash_pos = strchr(dash_pos, '_')) != NULL) { \
        *dash_pos = '-'; \
    };  \
    dst = dot_pos; \
} while (0)


#define _FILL_PARAM(arg, value) do {      \
    char str_arg[] = #arg;  \
    char* dot_pos; \
    if ((dot_pos = strchr(str_arg, '.')) != NULL) { \
        dot_pos = dot_pos + 1; \
    }  \
    char* dash_pos = dot_pos; \
    while ((dash_pos = strchr(dash_pos, '_')) != NULL) { \
        *dash_pos = '-'; \
    };  \
    if (strcmp(dot_pos, long_options[option_index].name) == 0) { \
            (arg) = (value);  \
            break; \
    }      \
} while (0)


struct options get_defaut_options()
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
            .gbsa = 0,

            .verbosity = DEBUG,
            .fix_receptor = 0,
            .fix_ligand = 0,
            .score_only = 0,
            .help = 0
    };

    return prms;
}


void free_options(struct options opts)
{
    if (opts.setup_json_root) {
        json_decref(opts.setup_json_root);
    }
}


static bool _fill_prms_from_json(struct options* prms, json_t* root)
{
    if (!json_is_object(root)) {
        return false;
    }

    json_t* opts = json_object_get(root, "options");
    if (!opts) {
        return true;
    }

    if (!json_is_object(opts)) {
        ERR_MSG("Key 'options' must point to a dictionary");
        json_decref(opts);
        return false;
    }

    json_error_t error;

    int code = json_unpack_ex(
            opts, &error, JSON_STRICT,

            "{s?i, s?i "
            " s?b, s?b, s?b, s?b, s?b, s?b, s?b, s?b, s?b, s?b, s?b, s?b, s?b, s?b "
            " s?s, s?s, s?s, s?s, s?s, s?s, s?s, s?s, s?s, s?s, s?s, s?s, s?s, s?s "
            " s?s, s?s, s?s, s?s, s?s, s?s, s?s, s?s, s?s, s?s}",

            // integer options
            "nsteps", &prms->nsteps,
            "verbosity", &prms->verbosity,

            // binary options
            "bonds", &prms->bonds,
            "angles", &prms->angles,
            "dihedrals", &prms->dihedrals,
            "impropers", &prms->impropers,
            "vdw", &prms->vdw,
            "vdw03", &prms->vdw03,
            "gbsa", &prms->gbsa,
            "fix-receptor", &prms->fix_receptor,
            "fix-ligand", &prms->fix_ligand,
            "score-only", &prms->score_only,
            "print-step", &prms->print_step,
            "print-stage", &prms->print_stage,
            "print-noe-matrix", &prms->print_noe_matrix,

            // char options
            "out-pdb", &prms->out_pdb,
            "out-json", &prms->out_json,
            "psf", &prms->psf,
            "prm", &prms->prm,
            "rtf", &prms->rtf,
            "pdb", &prms->pdb,
            "json", &prms->json,
            "rec-psf", &prms->rec_psf,
            "rec-prm", &prms->rec_prm,
            "rec-rtf", &prms->rec_rtf,
            "rec-pdb", &prms->rec_pdb,
            "rec-json", &prms->rec_json,
            "lig-psf", &prms->lig_psf,
            "lig-prm", &prms->lig_prm,
            "lig-rtf", &prms->lig_rtf,
            "lig-pdb", &prms->lig_pdb,
            "lig-json", &prms->lig_json,
            "fixed-pdb", &prms->fixed_pdb,
            "pair-springs-txt", &prms->pair_springs_txt,
            "point-springs-txt", &prms->point_springs_txt,
            "noe-txt", &prms->noe_txt,
            "noe-json", &prms->noe_json,
            "density-json", &prms->density_json
    );

    if (code != 0) {
        JSON_ERR_MSG(error, "Couldn't parse setup from json");
        json_decref(opts);
        return false;
    }

    return true;
}


static int _check_prms(struct options *prms)
{
    bool full = false;
    bool rec_sep = false;
    bool lig_sep = false;

    VERBOSITY = prms->verbosity;

    // Unpack setup json first and fill the global options if there are any
    if (prms->setup_json) {
        json_t* setup = read_json_file(prms->setup_json);
        if (!setup) {
            return 1;
        }

        if (!_fill_prms_from_json(prms, setup)) {
            json_decref(setup);
            ERR_MSG("Couldn't parse options");
            return 1;
        }
        //json_decref(setup);
        prms->setup_json_root = setup;
    }

    VERBOSITY = prms->verbosity;

    if (prms->json || prms->pdb || prms->psf || prms->prm || prms->rtf) {
        full = true;
        prms->separate = false;

        if (!(prms->json || (prms->pdb && prms->psf && prms->prm && prms->rtf))) {
            ERR_MSG("Wrong sfds");
            return 1;
        }

        if (prms->fix_receptor) {
            ERR_MSG("Can't fix receptor when provided as a single file");
            return 1;
        }

        if (prms->fix_ligand) {
            ERR_MSG("Can't fix ligand when provided as a single file");
            return 1;
        }

        INFO_MSG("Using prm %s", prms->prm);
    }

    if (prms->rec_pdb || prms->rec_psf || prms->rec_prm || prms->rec_rtf || prms->rec_json) {
        if (full) {
            ERR_MSG("Cant sep and full at the same time");
            return 1;
        }

        if (!(prms->rec_json || (prms->rec_pdb && prms->rec_psf && prms->rec_prm && prms->rec_rtf))) {
            ERR_MSG("Rec error");
            return 1;
        }

        rec_sep = true;

        if (!(prms->lig_json || (prms->lig_pdb && prms->lig_psf && prms->lig_prm && prms->lig_rtf))) {
            ERR_MSG("Lig error");
            return 1;
        }

        lig_sep = true;
        prms->separate = true;
    }

    if (!(full || (rec_sep && lig_sep))) {
        ERR_MSG("Nothing provided");
        return 1;
    }

    if (prms->nsteps < 0) {
        ERR_MSG("nsteps < 0");
        return 1;
    }

    return 0;
}


struct options parse_args(const int argc, char *const *argv, bool *error)
{
    struct options prms = get_defaut_options();

    struct option long_options[] =
            {
                    {"out-pdb",              required_argument, 0,                 0},
                    {"out-json",              required_argument, 0,                 0},

                    {"print-step", no_argument, &prms.print_step, 1},
                    {"print-stage", no_argument, &prms.print_stage, 1},
                    {"print-noe-matrix", no_argument, &prms.print_noe_matrix, 1},

                    {"psf",              required_argument, 0,                 0},
                    {"prm",              required_argument, 0,                 0},
                    {"rtf",              required_argument, 0,                 0},
                    {"pdb",              required_argument, 0,                 0},
                    {"json",             required_argument, 0,                 0},

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

                    {"fixed-pdb",             required_argument, 0,                 0},

                    {"pair-springs-txt",             required_argument, 0,                 0},
                    {"point-springs-txt",             required_argument, 0,                 0},

                    {"noe-txt",             required_argument, 0,                 0},
                    {"noe-json",             required_argument, 0,                 0},

                    {"density-json",             required_argument, 0,                 0},
                    {"setup-json",             required_argument, 0,                 0},

                    {"nsteps",           required_argument, 0,                 0},

                    {"bonds-on",             no_argument,  &prms.bonds,         1},
                    {"bonds-off",             no_argument,  &prms.bonds,         0},
                    {"angles-on",     no_argument,       &prms.angles,     1},
                    {"angles-off",     no_argument,       &prms.angles,     0},
                    {"dihedrals-on",          no_argument,       &prms.dihedrals,          1},
                    {"dihedrals-off",          no_argument,       &prms.dihedrals,          0},
                    {"impropers-on", no_argument,       &prms.impropers, 1},
                    {"impropers-off", no_argument,       &prms.impropers, 0},
                    {"vdw-on", no_argument,       &prms.vdw, 1},
                    {"vdw-off", no_argument,       &prms.vdw, 0},
                    {"vdw03-on", no_argument,       &prms.vdw03, 1},
                    {"vdw03-off", no_argument,       &prms.vdw03, 0},
                    {"gbsa-on", no_argument,       &prms.gbsa, 1},
                    {"gbsa-off", no_argument,       &prms.gbsa, 0},

                    {"fix-receptor",           no_argument,        &prms.fix_receptor,                 1},
                    {"fix-ligand",           no_argument,        &prms.fix_ligand,                 1},
                    {"score-only",           no_argument,        &prms.score_only,                 1},
                    {"verbosity",             no_argument,       &prms.verbosity,                 1},
                    {"help",             no_argument,       0,                 'h'},

                    {0, 0,                                  0,                 0}
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
                INFO_MSG("Option %s = %s", long_options[option_index].name, optarg);
                break;

            case 'h':
                usage_message(argv);
                *error = false;
                prms.help = true;
                return prms;

            case '?':
                usage_message(argv);
                *error = true;
                return prms;

            default:
                usage_message(argv);
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
    }

    if (_check_prms(&prms) != 0) {
        usage_message(argv);
        *error = true;
        return prms;
    }

    *error = false;
    return prms;
}


void usage_message(char *const *argv) {
    fprintf(stderr, "Usage msg\n");
    return;

    printf("\nUsage %s [ options ]\n\n", argv[0]);
    printf("\t--prm ------------------ Parameter file (required)\n"
           "\t--rtf ------------------ Topology file (required)\n"
           "\t--pdb ------------------ Full molecule pdb file (can contain multiple models)\n"
           "\t--psf ------------------ Full molecule psf file\n"
           "\t--json ----------------- Full molecule json file\n"
           "\t--rec-pdb -------------- Receptor pdb file\n"
           "\t--rec-psf -------------- Receptor psf file\n"
           "\t--rec-json ------------- Receptor json file\n"
           "\t--lig-pdb -------------- Ligand pdb file\n"
           "\t--lig-psf -------------- Ligand psf file\n"
           "\t--lig-json ------------- Ligand json file\n"
           "\t--out ------------------ Minimized structure (default: min.pdb)\n"
           "\t--protocol ------------- Protocol file\n"
           "\t--torsions-off --------- Turns torsional forces off\n"
           "\t--fitting-pdb ---------- PDB file for atomic density fitting\n"
           "\t--fitting-weight ------- Weight of atomic density fit in energy function\n"
           "\t--noe-params ----------- File with 2D NOE spectrum params\n"
           "\t--pairsprings ---------- File with pairwise distance restraints\n"
           "\t--pointsprings --------- File with point distance restraints\n"
           "\t--score_only ----------- Don't run minimization, just score the models\n"
           "\t--nsteps --------------- Number of minimization steps (default: 1000)\n"
           "\t--gbsa ----------------- Turn GBSA on (default: off)\n"
           "\t--json-log ------------- Log starting and final energy terms to specified json file\n"
           "\t--help ----------------- This message\n\n"

           "Protocol file format:\n\n"
           "\t> number_of_steps1 fixed_atoms_file1 pair_springs_file1 point_springs_file1\n"
           "\t> number_of_steps2 fixed_atoms_file2 pair_springs_file2 point_springs_file2\n"
           "\t> ...\n\n"

           "Pairwise springs file format:\n\n"
           "\t> number_of_springs\n"
           "\t> length absolute_error force_constant atom1_id(PDB) atom1_name atom2_id(PDB) atom2_name\n"
           "\t> ...\n\n"

           "Point springs (attached to a single point) file format:\n\n"
           "\t> number_of_springs\n"
           "\t> number_of_atoms_attached force_constant X0 Y0 Z0\n"
           "\t> atom1_id(from PDB) atom1_name\n"
           "\t> atom2_id(from PDB) atom2_name\n"
           "\t> ...\n\n"

           "NOE params file format (lines starting with # are skipped):\n"
           "\tAtom group file path\n"
           "\tExperimental matrix file path\n"
           "\tFrequency (in GHz): ex. 0.00006\n"
           "\tCorrelation time (in nanoseconds): ex. 0.1\n"
           "\tMixing time (in seconds): ex. 0.3\n"
           "\tDistance cutoff for atom groups (Angstroms): ex. 10\n"
           "\tMask everything except for the peaks on experimental matrix, when computing NMR energy (on/off)\n"
           "\tNMR energy weight in energy function\n\n"

           "\tNOE params file example:\n"
           "\t> # groups file\n"
           "\t> 104.fake_groups\n"
           "\t> # experimental matrix\n"
           "\t> 104.fake_mat\n"
           "\t> # frequency\n"
           "\t> 0.00006\n"
           "\t> # cor time\n"
           "\t> 0.1\n"
           "\t> # mix time\n"
           "\t> 0.3\n"
           "\t> # distance cutoff\n"
           "\t> 1000\n"
           "\t> # mask\n"
           "\t> off\n"
           "\t> # nmr energy weight\n"
           "\t> 100.0\n\n"

           "Fixed atoms file format: simply a slice of pdb containing fixed atoms\n\n"

           "If you wish to skip any of the files in the protocol file,\n"
           "just replace it with a dot \".\". For instance if no fixed atoms\n"
           "are needed then the line should be \n"
           "\t1000 . spring_pairs spring_points\n\n"

           "WARNING: If ligand and receptor are specified separately,\n"
           "\tthen ligand atom ids need to be increased by the number of receptor atoms.\n"
           "\tConsistency is verified by comparing atom names.\n\n");
}
