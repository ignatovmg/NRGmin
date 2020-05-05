#include <string.h>


#include "parse_options.h"
#include "utils.h"


#define _FILL_PARAM(arg, value) do {                       \
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
            arg = value;                               \
    }                                                         \
} while (0)


static struct options _get_defaut_options()
{
    struct options prms = {
            .out_pdb = "out.pdb",
            .out_json = "out.json",

            .separate = false,

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

            .nsteps = 1000,

            .bonds = 1,
            .angles = 1,
            .dihedrals = 1,
            .impropers = 1,
            .vdw = 1,
            .vdw03 = 1,
            .gbsa = 0,

            .verbose = 0,
            .fix_receptor = 0,
            .score_only = 0,
            .help = 0
    };

    return prms;
}


static int _check_prms(struct options prms)
{
    if (prms.json) {
        prms.separate = false;
        INFO_MSG("Using json %s", prms.json);

    } else if (prms.pdb && prms.psf && prms.prm && prms.rtf) {
        prms.separate = false;
        INFO_MSG("Using prm %s", prms.prm);

    } else if (prms.rec_pdb || prms.rec_psf || prms.rec_prm || prms.rec_rtf || prms.rec_json) {
        prms.separate = true;

        if (prms.rec_json) {
            INFO_MSG("Using json %s", prms.rec_json);

        } else if (prms.rec_pdb && prms.rec_psf && prms.rec_prm && prms.rec_rtf) {
            INFO_MSG("Using prm %s", prms.rec_prm);

        } else {
            ERR_MSG("kjhlkjh");
            return 1;
        }

        if (prms.lig_json) {
            INFO_MSG("Using json %s", prms.lig_json);

        } else if (prms.lig_pdb && prms.lig_psf && prms.lig_prm && prms.lig_rtf) {
            INFO_MSG("Using prm %s", prms.lig_prm);

        } else {
            ERR_MSG("sgsg");
            return 1;
        }
    } else {
        ERR_MSG("kljh");
        return 1;
    }

    if (prms.nsteps < 0) {
        ERR_MSG("sdf");
        return 1;
    }

    return 0;
}


struct options parse_args(const int argc, const char** argv, bool *error)
{
    struct options prms = _get_defaut_options();

    struct option long_options[] =
            {
                    {"out-pdb",              required_argument, 0,                 0},
                    {"out-json",              required_argument, 0,                 0},

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
                    {"score-only",           no_argument,        &prms.score_only,                 1},
                    {"verbose",             no_argument,       &prms.verbose,                 1},
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
                printf("Option %s", long_options[option_index].name);
                if (optarg)
                    printf(" = %s", optarg);
                printf("\n");
                break;

            case 'h':
                usage_message(argv);
                printf("ASDF");
                exit(EXIT_SUCCESS);
                break;

            case '?':
                usage_message(argv);
                printf("ASDF");
                exit(EXIT_FAILURE);
                break;

            default:
                usage_message(argv);
                printf("ASDF");
                exit(EXIT_FAILURE);
                break;
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

    if (_check_prms(prms) != 0) {
        *error = true;
    }

    return prms;
}


void usage_message(const char **argv) {
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
