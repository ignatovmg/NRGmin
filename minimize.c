#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <jansson.h>
#include <getopt.h>

#include "mol2/json.h"
#include "mol2/benergy.h"
#include "mol2/gbsa.h"
#include "mol2/icharmm.h"
#include "mol2/minimize.h"
#include "mol2/nbenergy.h"
#include "mol2/pdb.h"
#include "mol2/fitting.h"

#include "nmrgrad/noe.h"

#define __TOL__ 5E-4

#define ERR_MSG(fmt, ...) do {                                              \
    fprintf(stderr,"[Error] (file %s, %s, line %i):\n" fmt "\n"             \
            "Exiting ...\n", __FILE__, __func__, __LINE__, ##__VA_ARGS__);  \
    exit(EXIT_FAILURE);                                                     \
} while(0)

#define WRN_MSG(fmt, ...) do {                                         \
    fprintf(stderr,"[Warning] (file %s, %s, line %i):\n" fmt "\n",     \
            __FILE__, __func__, __LINE__, ##__VA_ARGS__);              \
} while(0)

#define INFO_MSG(fmt, ...) do {                               \
    fprintf(stderr,"[Info]: " fmt);                           \
} while(0)

#define MYFOPEN(fp, path, spec)  fopen(path, spec); {         \
    if (fp == NULL) {                                         \
        ERR_MSG("ERROR opening file %s\n", path);              \
    }                                                         \
}

#define MYCALLOC(ptr, n, size)  calloc(n, size); {            \
    if (ptr == NULL) {                                        \
        ERR_MSG("ERROR allocating %i bytes for variable %s\n",\
            (int)((n))*(int)((size)), #ptr);                   \
    }                                                         \
}

#define READ_WORD(f, word, line) do { \
    word = NULL; \
    while (fgets(line, 512, f) != NULL) { \
        if (line[0] != '#') {\
            word = strtok(line, " \t\n"); \
            if (word == NULL || word[0] == '#') { \
                ERR_MSG("Line cannot be empty\n"); \
            } \
            break; \
        } \
    } \
} while(0)

#define MAGIC_ARGS(arg) {                                     \
    char name[] = #arg;                                       \
    char* pos = strchr(name, '_');                            \
    if (pos != NULL) *pos = '-';                              \
    if (strcmp(name, long_options[option_index].name) == 0) { \
            arg = optarg;                                     \
    }                                                         \
}

struct energy_prm {
    struct mol_atom_group *ag;
    struct agsetup *ag_setup;
    struct acesetup *ace_setup;
    struct springset_pairs *sprst_pairs;
    struct springset_points *sprst_points;
    struct noesy_spectrum *nmr;
    struct fitting_params *fit_prms;
    bool no_geom;
    bool bonds;
    bool angles;
    bool torsions;
    bool vdw;
};

struct pair_spring {
    struct mol_atom_group *ag;  /**< affected atomgroup */
    int laspr[2];      /**< list of atoms */
    double lnspr;
    double erspr;
    double fkspr;      /**< force constant */
};

struct point_spring {
    struct mol_atom_group *ag;  /**< affected atomgroup */
    int naspr;      /**< number of affected atoms */
    int *laspr;      /**< list of atoms */
    double fkspr;      /**< force constant */
    double X0, Y0, Z0; /**< anchor point */
};

struct springset_pairs {
    int nsprings;       /**< number of springs */
    struct pair_spring *springs;  /**< array of springs */
};

struct springset_points {
    int nsprings;       /**< number of springs */
    struct point_spring *springs;  /**< array of springs */
};

struct noesy_spectrum {
    struct nmr_noe *spec;
    bool print_noe_matrix;
    double weight;
};

struct fitting_params {
    struct mol_atom_group **ag_list;
    int ag_count;
    double weight;
    struct mol_fitting_params prms;
};

void read_fix(char *ffile, int *nfix, size_t **fix);

struct springset_pairs *read_springset_pairs(struct mol_atom_group *ag, char *sfile);

struct springset_points *read_springset_points(struct mol_atom_group *ag, char *sfile);

struct noesy_spectrum *read_noesy_spectrum(struct mol_atom_group *ag, char *sfile);

struct fitting_params *make_fitting_prms(char *fitting_pdblist, double weight, double radius);

void free_fitting_prms(struct fitting_params *prms);

void free_noesy_spectrum(struct noesy_spectrum *nmr);

void free_springset_pairs(struct springset_pairs **sprst);

void free_springset_points(struct springset_points **sprst);

void springeng_pair(struct springset_pairs *sprst, double *een);

void springeng_point(struct springset_points *sprst, double *een);

static lbfgsfloatval_t energy_func(
        void *restrict prm,
        const double *restrict array,
        double *restrict gradient,
        const int array_size,
        const lbfgsfloatval_t step);

static void fprint_energy_terms(FILE *stream, void *restrict prm, char *prefix);

void usage_message(char **argv);

void help_message(void);

int main(int argc, char **argv) {
    mol_enable_floating_point_exceptions();

    //static int verbose_flag = 0;
    static int ace_flag = 0;
    static int torsions_off = 0;
    static int print_noe_matrix = 0;
    static int score_only = 0;
    static int fix_rec = 0;
    static int nonpar_terms_only = 0; // if geometry is not provided
    double fitting_weight = 1.;
    char *protocol = NULL;
    char *pdb = NULL;
    char *psf = NULL;
    char *prm = NULL;
    char *json = NULL;
    char *rtf = NULL;
    char *out = NULL;
    char *rec_pdb = NULL;
    char *rec_psf = NULL;
    char *rec_json = NULL;
    char *lig_pdb = NULL;
    char *lig_psf = NULL;
    char *lig_json = NULL;
    char *noe_params = NULL;
    char *pairsprings = NULL;
    char *pointsprings = NULL;
    char *fixed_pdb = NULL;
    char *fitting_pdblist = NULL;
    int nsteps = 1000;

    static struct option long_options[] =
            {
                    //{"verbose", no_argument, &verbose_flag, 1},
                    {"rtf",              required_argument, 0,                 0},
                    {"prm",              required_argument, 0,                 0},
                    {"out",              required_argument, 0,                 0},
                    {"protocol",         required_argument, 0,                 0},
                    {"pdb",              required_argument, 0,                 0},
                    {"psf",              required_argument, 0,                 0},
                    {"json",             required_argument, 0,                 0},
                    {"rec-pdb",          required_argument, 0,                 0},
                    {"rec-psf",          required_argument, 0,                 0},
                    {"rec-json",         required_argument, 0,                 0},
                    {"lig-pdb",          required_argument, 0,                 0},
                    {"lig-psf",          required_argument, 0,                 0},
                    {"lig-json",         required_argument, 0,                 0},
                    {"nsteps",           required_argument, 0,                 0},
                    {"noe-params",       required_argument, 0,                 0},
                    {"pairsprings",      required_argument, 0,                 0},
                    {"pointsprings",     required_argument, 0,                 0},
                    {"fitting-pdblist",  required_argument, 0,                 0},
                    {"fitting-weight",   required_argument, 0,                 0},
                    {"fixed-pdb",   required_argument, 0,                 0},
                    {"gbsa",             no_argument,       &ace_flag,         1},
                    {"torsions-off",     no_argument,       &torsions_off,     1},
                    {"fix-rec",          no_argument,       &fix_rec,          1},
                    {"print_noe_matrix", no_argument,       &print_noe_matrix, 1},
                    {"score_only",       no_argument,       &score_only,       1},
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
                exit(EXIT_SUCCESS);
                break;

            case '?':
                usage_message(argv);
                exit(EXIT_FAILURE);
                break;

            default:
                usage_message(argv);
                exit(EXIT_FAILURE);
                break;
        }

        MAGIC_ARGS(rtf);
        MAGIC_ARGS(prm);
        MAGIC_ARGS(out);
        MAGIC_ARGS(pdb);
        MAGIC_ARGS(json);
        MAGIC_ARGS(psf);
        MAGIC_ARGS(protocol);
        MAGIC_ARGS(rec_pdb);
        MAGIC_ARGS(rec_psf);
        MAGIC_ARGS(rec_json);
        MAGIC_ARGS(lig_pdb);
        MAGIC_ARGS(lig_psf);
        MAGIC_ARGS(lig_json);
        MAGIC_ARGS(noe_params);
        MAGIC_ARGS(pairsprings);
        MAGIC_ARGS(pointsprings);
        MAGIC_ARGS(fitting_pdblist);
        MAGIC_ARGS(fixed_pdb);

        if (strcmp("nsteps", long_options[option_index].name) == 0) {
            nsteps = atoi(optarg);
        }
        if (strcmp("fitting-weight", long_options[option_index].name) == 0) {
            fitting_weight = atof(optarg);
        }
    }

    struct mol_atom_group *ag = NULL;
    struct mol_atom_group *ag_json = NULL;
    struct mol_atom_group_list *aglist = NULL;
    struct mol_atom_group *rec_ag = NULL;
    struct mol_atom_group *lig_ag = NULL;

    if (protocol != NULL && score_only != 0) {
        ERR_MSG("--protocol is not effective, when --score_only flag is provided\n");
    }

    if ((rtf == NULL || prm == NULL) && json == NULL && (lig_json == NULL || rec_json)) {
        if (score_only != 0) {
            nonpar_terms_only = 1;
            WRN_MSG("Without parameter files or json files score_only can compute only NOE, Pairsprings, Pointsprings and Density score\n");
        } else {
            ERR_MSG("A pair of RTF and PRM files, a JSON file or --scory_only flag is required\n");
        }
    }

    if (nsteps < 0) {
        ERR_MSG("Number of steps must be non-negative (nsteps = %i)\n", nsteps);
    }

    if (out == NULL) {
        out = "min.pdb";
    }

    if (json != NULL) {
        INFO_MSG("Reading json file\n");
        ag_json = mol_read_json(json);

        if  (psf != NULL) {
            WRN_MSG("Parameter --psf is not effective, since --json was provided\n");
        }
    }

    if (rec_json != NULL) {
        INFO_MSG("Reading rec json file\n");
        rec_ag = mol_read_json(rec_json);

        if  (rec_psf != NULL) {
            WRN_MSG("Parameter --rec-psf is not effective, since --rec-json was provided\n");
        }
    }

    if (lig_json != NULL) {
        INFO_MSG("Reading lig json file\n");
        lig_ag = mol_read_json(lig_json);

        if  (lig_psf != NULL) {
            WRN_MSG("Parameter --lig-psf is not effective, since --lig-json was provided\n");
        }
    }

    if (rec_pdb != NULL) {
        if (rec_json != NULL) {
            INFO_MSG("Reading receptor coordinates from pdb file, using parameters from --rec-json\n");
            struct mol_atom_group* rec_ag_crd = mol_read_pdb(rec_pdb);
            if (rec_ag->natoms != rec_ag_crd->natoms) {
                ERR_MSG("--rec-json and --rec-pdb have different atom numbers (%i, %i)", (int)rec_ag->natoms, (int)rec_ag_crd->natoms);
            }
            memcpy(rec_ag->coords, rec_ag_crd->coords, sizeof(struct mol_vector3) * rec_ag->natoms);
            mol_atom_group_free(rec_ag_crd);

        } else {
            INFO_MSG("Reading rec pdb file\n");
            rec_ag = mol_read_pdb(rec_pdb);

            if  (rec_psf != NULL) {
                INFO_MSG("Reading rec psf file\n");
                mol_atom_group_read_geometry(rec_ag, rec_psf, prm, rtf);

            } else if (score_only != 0 && nonpar_terms_only != 0) {
                WRN_MSG("Geometry for the receptor is not provided, computing only non-parametric terms\n");
            } else {
                ERR_MSG("--rec-psf or --rec-json must be provided with --rec-pdb\n");
            }
        }
    }

    if (lig_pdb != NULL) {
        if (lig_json != NULL) {
            INFO_MSG("Reading ligand coordinates from pdb file, using parameters from json\n");
            struct mol_atom_group* lig_ag_crd = mol_read_pdb(lig_pdb);
            if (lig_ag->natoms != lig_ag_crd->natoms) {
                ERR_MSG("--lig-json and --lig-pdb have different atom numbers (%i, %i)", (int)lig_ag->natoms, (int)lig_ag_crd->natoms);
            }
            memcpy(lig_ag->coords, lig_ag_crd->coords, sizeof(struct mol_vector3) * lig_ag->natoms);
            mol_atom_group_free(lig_ag_crd);

        } else {
            INFO_MSG("Reading lig pdb file\n");
            lig_ag = mol_read_pdb(lig_pdb);

            if  (lig_psf != NULL) {
                INFO_MSG("Reading lig psf file\n");
                mol_atom_group_read_geometry(lig_ag, lig_psf, prm, rtf);
            } else if (score_only != 0 && nonpar_terms_only != 0) {
                WRN_MSG("Geometry for the ligand is not provided, computing only non-parametric terms\n");
            } else {
                ERR_MSG("--lig-psf or --lig-json must be provided with --lig-pdb\n");
            }
        }
    }

    if (pdb != NULL) {
        INFO_MSG("Trying to read the pdb file as a multimodel one (with MODEL records)\n");
        aglist = mol_read_pdb_models(pdb);
        if (aglist == NULL) {
            INFO_MSG("The pdb file doesn't have MODEL records. Reading as a regular pdb file\n");
            ag = mol_read_pdb(pdb);
            if (ag == NULL) {
                ERR_MSG("Failed reading pdb file");
            }
            aglist = mol_atom_group_list_create(1);
            aglist->members[0] = *ag;
        }

        if (ag_json != NULL) {
            INFO_MSG("Reading coordinates from --pdb and geometry from --json\n");
            // Copy geometry from json to aglist
            struct mol_atom_group_list *fin_aglist = mol_atom_group_list_create(aglist->size);
            for (size_t i = 0; i < fin_aglist->size; i++) {
                fin_aglist->members[i] = *mol_atom_group_copy(ag_json);
                if (ag_json->natoms != aglist->members[i].natoms) {
                    ERR_MSG("Model %i in --pdb and --json have different atom numbers (%i, %i)", (int)i, (int)aglist->members[i].natoms, (int)ag_json->natoms);
                }
                memcpy(fin_aglist->members[i].coords,
                       aglist->members[i].coords,
                       sizeof(struct mol_vector3) * aglist->members[i].natoms);
            }
            mol_atom_group_list_free(aglist);
            mol_atom_group_free(ag_json);
            aglist = fin_aglist;
            ag_json = NULL;

        } else if (psf != NULL) {
            INFO_MSG("Reading geometry from --psf\n");
            for (size_t i = 0; i < aglist->size; i++) {
                mol_atom_group_read_geometry(&aglist->members[i], psf, prm, rtf);
            }
        } else if (score_only != 0 && nonpar_terms_only != 0) {
            WRN_MSG("Geometry for the molecule is not provided, computing only non-parametric terms\n");
        } else {
            ERR_MSG("--json or --psf must be provided with --pdb");
        }
    }

    // If rec and lig were supplied separately
    if (aglist == NULL) {
        aglist = mol_atom_group_list_create(1);

        if (ag_json != NULL) {
            INFO_MSG("Using single model from file provided with --json\n");
            aglist->members[0] = *ag_json;

        } else if (rec_ag != NULL && lig_ag != NULL) {
            INFO_MSG("Using single model assembled from receptor and ligand provided separately\n");
            ag = mol_atom_group_join(rec_ag, lig_ag);
            aglist->members[0] = *ag;
            mol_atom_group_free(lig_ag);
            mol_atom_group_free(rec_ag);

        } else {
            ERR_MSG("Couldn't create an atom group\n");
        }
    } else {
        INFO_MSG("Using model(s) from file provided with --pdb or --json\n");
    }

    // Fill minimization parameters
    struct energy_prm engpar;
    engpar.no_geom = false;
    engpar.bonds = true;
    engpar.angles = true;
    engpar.torsions = true;
    engpar.vdw = true;

    // Torsions
    if (torsions_off == 1) {
        engpar.torsions = false;
    }

    // If geometry is not provided
    if (nonpar_terms_only != 0) {
        engpar.no_geom = true;
        engpar.bonds = false;
        engpar.angles = false;
        engpar.torsions = false;
        engpar.vdw = false;
    }

    // NMR 2D spectrum
    if (noe_params != NULL) {
        engpar.nmr = read_noesy_spectrum(&aglist->members[0], noe_params);
        engpar.nmr->print_noe_matrix = print_noe_matrix;
    } else {
        engpar.nmr = NULL;
    }

    // Pairsprings
    if (pairsprings != NULL) {
        engpar.sprst_pairs = read_springset_pairs(ag, pairsprings);
    } else {
        engpar.sprst_pairs = NULL;
    }

    // Point springs
    if (pointsprings != NULL) {
        engpar.sprst_points = read_springset_points(ag, pointsprings);
    } else {
        engpar.sprst_points = NULL;
    };

    // Density fitting
    if (fitting_pdblist != NULL) {
        engpar.fit_prms = make_fitting_prms(fitting_pdblist, fitting_weight, 2.0);
    } else {
        engpar.fit_prms = NULL;
    }

    // Fixed part
    int nfix_glob = 0;
    size_t *fix_glob = NULL;
    if (fixed_pdb != NULL) {
        read_fix(fixed_pdb, &nfix_glob, &fix_glob);
    }

    FILE *outfile = MYFOPEN(outfile, out, "w");
    for (int modeli = 0; modeli < aglist->size; modeli++) {
        if (aglist->size > 1) {
            fprintf(outfile, "MODEL %i\n", (modeli + 1));
        }

        ag = &aglist->members[modeli];
        engpar.ag = ag;

        struct agsetup ags;
        struct acesetup ace_setup;

        if (!engpar.no_geom) {
            ag->gradients = MYCALLOC(ag->gradients, ag->natoms, sizeof(struct mol_vector3))
            mol_fixed_init(ag);
            mol_fixed_update(ag, nfix_glob, fix_glob);

            init_nblst(ag, &ags);
            update_nblst(ag, &ags);

            if (ace_flag == 1) {
                ace_setup.efac = 0.5;
                ace_ini(ag, &ace_setup);
                ace_fixedupdate(ag, &ags, &ace_setup);
                ace_updatenblst(&ags, &ace_setup);
            }
        }

        if (score_only == 0) {
            if (protocol != NULL) {
                FILE *prot_file = MYFOPEN(prot_file, protocol, "r");
                int cur_nsteps;
                char fix_path[1024];
                char spr_pair_path[1024];
                char spr_point_path[1024];

                int c;
                while ((c = fscanf(prot_file, "%i %s %s %s", &cur_nsteps, fix_path,
                                   spr_pair_path, spr_point_path)) != EOF) {
                    if (c != 4) {
                        ERR_MSG("Wrong protocol file format (%i words read)\n", c);
                    }

                    if (cur_nsteps < 0) {
                        ERR_MSG("Wrong protocol file format (number of steps must be non-negative)\n");
                    }

                    int nfix = 0;
                    size_t *fix = NULL;

                    if (strcmp(fix_path, ".") != 0) {
                        read_fix(fix_path, &nfix, &fix);
                        mol_fixed_update(ag, nfix, fix);
                        update_nblst(ag, &ags);

                        if (ace_flag == 1) {
                            ace_fixedupdate(ag, &ags, &ace_setup);
                            ace_updatenblst(&ags, &ace_setup);
                        }

                        free(fix);
                        fix = NULL;
                    }

                    struct springset_pairs *sprst_pairs = NULL;
                    struct springset_points *sprst_points = NULL;

                    if (strcmp(spr_pair_path, ".") != 0) {
                        sprst_pairs = read_springset_pairs(ag, spr_pair_path);
                    }

                    if (strcmp(spr_point_path, ".") != 0) {
                        sprst_points = read_springset_points(ag, spr_point_path);
                    }

                    engpar.ag_setup = &ags;
                    engpar.sprst_pairs = sprst_pairs;
                    engpar.sprst_points = sprst_points;

                    if (ace_flag == 1) {
                        engpar.ace_setup = &ace_setup;
                    } else {
                        engpar.ace_setup = NULL;
                    }

                    fprint_energy_terms(outfile, &engpar, "REMARK START ");

                    if (cur_nsteps > 0) {
                        mol_minimize_ag(MOL_LBFGS, cur_nsteps, __TOL__, ag, (void *) (&engpar), energy_func);
                    }

                    free_springset_pairs(&sprst_pairs);
                    free_springset_points(&sprst_points);
                }
            } else {
                engpar.ag_setup = &ags;

                if (ace_flag == 1) {
                    engpar.ace_setup = &ace_setup;
                } else {
                    engpar.ace_setup = NULL;
                }

                fprint_energy_terms(outfile, &engpar, "REMARK START ");

                if (nsteps > 0) {
                    mol_minimize_ag(MOL_LBFGS, nsteps, __TOL__, ag, (void *) (&engpar), energy_func);
                }
            }
        }

        fprintf(outfile, "REMARK\n");
        fprint_energy_terms(outfile, &engpar, "REMARK FINAL ");

        mol_fwrite_pdb(outfile, ag);

        if (aglist->size > 1) {
            fprintf(outfile, "ENDMDL\n");
        }
        fflush(outfile);
    }

    free_noesy_spectrum(engpar.nmr);
    free_fitting_prms(engpar.fit_prms);

    if (engpar.sprst_points != NULL) {
        free_springset_points(&engpar.sprst_points);
    }
    if (engpar.sprst_pairs != NULL) {
        free_springset_pairs(&engpar.sprst_pairs);
    }

    if (fix_glob != NULL) {
        free(fix_glob);
    }

    fclose(outfile);
    mol_atom_group_list_free(aglist);

    return EXIT_SUCCESS;
}


static lbfgsfloatval_t energy_func(
        void *restrict prm,
        const double *restrict array,
        double *restrict gradient,
        const int array_size,
        const lbfgsfloatval_t step) {
    lbfgsfloatval_t energy = 0.0;
    struct energy_prm *energy_prm = (struct energy_prm *) prm;

    if (array != NULL) {
        assert(array_size == energy_prm->ag->active_atoms->size * 3);
        mol_atom_group_set_actives(energy_prm->ag, array);
    }
    bool updated = check_clusterupdate(energy_prm->ag, energy_prm->ag_setup);
    if (updated) {
        if (energy_prm->ace_setup != NULL) {
            ace_updatenblst(energy_prm->ag_setup, energy_prm->ace_setup);
        }
    }

    mol_zero_gradients(energy_prm->ag);

    if (energy_prm->ace_setup != NULL) {
        aceeng(energy_prm->ag, &energy, energy_prm->ace_setup, energy_prm->ag_setup);
    }

    if (energy_prm->vdw) {
        vdweng(energy_prm->ag, &energy, energy_prm->ag_setup->nblst);
        vdwengs03(1.0, energy_prm->ag_setup->nblst->nbcof, energy_prm->ag, &energy,
                  energy_prm->ag_setup->nf03, energy_prm->ag_setup->listf03);
    }

    if (energy_prm->bonds) {
        beng(energy_prm->ag, &energy);
    }

    if (energy_prm->angles) {
        aeng(energy_prm->ag, &energy);
    }

    if (energy_prm->torsions) {
        teng(energy_prm->ag, &energy);
        ieng(energy_prm->ag, &energy);
    }

    if (energy_prm->sprst_pairs != NULL) {
        springeng_pair(energy_prm->sprst_pairs, &energy);
    }

    if (energy_prm->sprst_points != NULL) {
        springeng_point(energy_prm->sprst_points, &energy);
    }

    if (energy_prm->nmr != NULL) {
        struct nmr_noe *spec = energy_prm->nmr->spec;
        nmr_r2_mat(spec->r2, energy_prm->ag, spec->grps);
        nmr_compute_peaks(spec, energy_prm->ag);
        nmr_energy(spec, energy_prm->ag->gradients, energy_prm->nmr->weight);
        energy += spec->energy;
    }

    if (energy_prm->fit_prms != NULL) {
        energy += mol_fitting_score_aglist(energy_prm->ag,
                                           energy_prm->fit_prms->ag_list,
                                           energy_prm->fit_prms->ag_count,
                                           &energy_prm->fit_prms->prms,
                                           energy_prm->fit_prms->weight);
    }

    if (gradient != NULL) {
        for (int i = 0; i < array_size / 3; i++) {
            int atom_i = energy_prm->ag->active_atoms->members[i];
            gradient[3 * i] = -energy_prm->ag->gradients[atom_i].X;
            gradient[3 * i + 1] = -energy_prm->ag->gradients[atom_i].Y;
            gradient[3 * i + 2] = -energy_prm->ag->gradients[atom_i].Z;
        }
    }

    return energy;
}

static void fprint_energy_terms(FILE *stream, void *restrict prm, char *prefix) {
    lbfgsfloatval_t energy = 0.0;
    lbfgsfloatval_t total = 0.0;
    struct energy_prm *energy_prm = (struct energy_prm *) prm;

    char fmt[1024];
    strcpy(fmt, prefix);

    if (!energy_prm->no_geom) {
        bool updated = check_clusterupdate(energy_prm->ag, energy_prm->ag_setup);
        if (updated) {
            if (energy_prm->ace_setup != NULL) {
                ace_updatenblst(energy_prm->ag_setup, energy_prm->ace_setup);
            }
        }

        if (energy_prm->ace_setup != NULL) {
            aceeng(energy_prm->ag, &energy, energy_prm->ace_setup, energy_prm->ag_setup);
            strcpy(fmt, prefix);
            fprintf(stream, strcat(fmt, "ACE: % .3f\n"), energy);
            total += energy;
            energy = 0.0;
        }

        if (energy_prm->vdw) {
            vdweng(energy_prm->ag, &energy, energy_prm->ag_setup->nblst);
            strcpy(fmt, prefix);
            fprintf(stream, strcat(fmt, "VWD: % .3f\n"), energy);
            total += energy;
            energy = 0.0;

            vdwengs03(1.0, energy_prm->ag_setup->nblst->nbcof, energy_prm->ag, &energy,
                      energy_prm->ag_setup->nf03, energy_prm->ag_setup->listf03);
            strcpy(fmt, prefix);
            fprintf(stream, strcat(fmt, "VWD03: % .3f\n"), energy);
            total += energy;
            energy = 0.0;
        }

        if (energy_prm->bonds) {
            beng(energy_prm->ag, &energy);
            strcpy(fmt, prefix);
            fprintf(stream, strcat(fmt, "Bonded: % .3f\n"), energy);
            total += energy;
            energy = 0.0;
        }

        if (energy_prm->angles) {
            aeng(energy_prm->ag, &energy);
            strcpy(fmt, prefix);
            fprintf(stream, strcat(fmt, "Angles: % .3f\n"), energy);
            total += energy;
            energy = 0.0;
        }

        if (energy_prm->torsions) {
            teng(energy_prm->ag, &energy);
            strcpy(fmt, prefix);
            fprintf(stream, strcat(fmt, "Torsions: % .3f\n"), energy);
            total += energy;
            energy = 0.0;

            ieng(energy_prm->ag, &energy);
            strcpy(fmt, prefix);
            fprintf(stream, strcat(fmt, "Impropers: % .3f\n"), energy);
            total += energy;
            energy = 0.0;
        }
    }

    if (energy_prm->sprst_pairs != NULL) {
        springeng_pair(energy_prm->sprst_pairs, &energy);
        strcpy(fmt, prefix);
        fprintf(stream, strcat(fmt, "Pairsprings: % .3f\n"), energy);
        total += energy;
        energy = 0.0;
    }

    if (energy_prm->sprst_points != NULL) {
        springeng_point(energy_prm->sprst_points, &energy);
        strcpy(fmt, prefix);
        fprintf(stream, strcat(fmt, "Pointsprings: % .3f\n"), energy);
        total += energy;
        energy = 0.0;
    }

    if (energy_prm->nmr != NULL) {
        nmr_r2_mat(energy_prm->nmr->spec->r2, energy_prm->ag, energy_prm->nmr->spec->grps);
        nmr_compute_peaks_no_grad(energy_prm->nmr->spec, energy_prm->ag);
        energy = nmr_energy(energy_prm->nmr->spec, NULL, energy_prm->nmr->weight);

        if (energy_prm->nmr->print_noe_matrix) {
            strcpy(fmt, prefix);
            fprintf(stdout, strcat(fmt, "NOESY_MATRIX_SCALE: %e\n"), energy_prm->nmr->spec->scale);
            nmr_matrix_fwrite(stdout, energy_prm->nmr->spec->in, energy_prm->nmr->spec->size);
        }

        strcpy(fmt, prefix);
        fprintf(stream, strcat(fmt, "NOE: % .6f\n"), energy);
        total += energy;
        energy = 0.0;
    }

    if (energy_prm->fit_prms != NULL) {
        energy = mol_fitting_score_aglist(energy_prm->ag,
                                          energy_prm->fit_prms->ag_list,
                                          energy_prm->fit_prms->ag_count,
                                          &energy_prm->fit_prms->prms,
                                          energy_prm->fit_prms->weight);
        strcpy(fmt, prefix);
        fprintf(stream, strcat(fmt, "Density: % .4f\n"), energy);
        total += energy;
        energy = 0.0;
    }

    strcpy(fmt, prefix);
    fprintf(stream, strcat(fmt, "Total: % .3f\n"), total);
}


struct fitting_params *make_fitting_prms(char *fitting_pdblist, double weight, double radius) {
    struct fitting_params *fit_prms = calloc(1, sizeof(struct fitting_params));

    /*struct mol_atom_group_list* ag_list_fit = mol_read_pdb_models(fitting_pdb);
    if (ag_list_fit == NULL) {
        ag_list_fit = mol_atom_group_list_create(1);
        mol_atom_group_free(ag_list_fit->members);
        ag_list_fit->members = mol_read_pdb(fitting_pdb);
    }*/
    FILE *f = fopen(fitting_pdblist, "r");
    char pdb_file[512];
    struct mol_atom_group **aglist = calloc(1, sizeof(struct mol_atom_group *));
    int ag_count = 0;
    int cur_size = 1;
    while (fgets(pdb_file, 512, f) != NULL) {
        struct mol_atom_group *ag = mol_read_pdb(pdb_file);
        ag_count++;

        if (ag_count > cur_size) {
            cur_size *= 2;
            aglist = realloc(aglist, cur_size * sizeof(struct mol_atom_group *));
        }

        aglist[ag_count - 1] = ag;
    }

    fit_prms->ag_list = aglist;
    fit_prms->ag_count = ag_count;
    fit_prms->prms.radius = radius;
    fit_prms->weight = weight;
    return fit_prms;
}


void free_fitting_prms(struct fitting_params *prms) {
    if (prms != NULL) {
        for (int i = 0; i < prms->ag_count; i++) {
            mol_atom_group_free(prms->ag_list[i]);
        }
        free(prms);
    }
}


struct noesy_spectrum *read_noesy_spectrum(struct mol_atom_group *ag, char *sfile) {
    struct noesy_spectrum *nmr = calloc(1, sizeof(struct noesy_spectrum));

    char line[512];
    FILE *f = MYFOPEN(f, sfile, "r");

    char *word;
    READ_WORD(f, word, line);
    char groups_path[512];
    strcpy(groups_path, word);

    READ_WORD(f, word, line);
    char matrix_path[512];
    strcpy(matrix_path, word);

    READ_WORD(f, word, line);
    double freq = atof(word);

    READ_WORD(f, word, line);
    double corr_time = atof(word);

    READ_WORD(f, word, line);
    double mix_time = atof(word);

    READ_WORD(f, word, line);
    double dist_cutoff = atof(word);

    READ_WORD(f, word, line);
    bool mask_on = false;
    if (strcmp(word, "on\0") == 0) {
        mask_on = true;
    } else if (strcmp(word, "off\0") != 0) {
        ERR_MSG("Wrong value (on/off)");
    }

    READ_WORD(f, word, line);
    double weight = atof(word);

    fclose(f);

    struct nmr_group_list *groups = nmr_group_list_read(groups_path);
    nmr->spec = nmr_noe_create(groups, freq, corr_time, mix_time, dist_cutoff);
    nmr->spec->in_grad = nmr_grad_create(ag->natoms, groups->ngroups);
    nmr->spec->rx_grad = nmr_grad_create(ag->natoms, groups->ngroups);

    nmr->spec->omega = freq;
    nmr->spec->t_mix = mix_time;
    nmr->spec->t_cor = corr_time;
    nmr->spec->cutoff = dist_cutoff;

    int *mask = NULL;
    if (mask_on) {
        mask = calloc(groups->ngroups * groups->ngroups, sizeof(int));
    }
    nmr->spec->exp = nmr_matrix_read(matrix_path, groups->ngroups, mask);

    nmr->weight = weight;

    nmr->print_noe_matrix = false;

    return nmr;
}

void free_noesy_spectrum(struct noesy_spectrum *nmr) {
    if (nmr != NULL) {
        nmr_noe_free(nmr->spec);
        free(nmr);
    }
}

void read_fix(char *ffile, int *nfix, size_t **fix) {
    int linesz = 91;
    char *buffer = (char *) MYCALLOC(buffer, linesz, sizeof(char));

    *nfix = 0;
    FILE *fp = MYFOPEN(fp, ffile, "r");

    while (fgets(buffer, linesz - 1, fp) != NULL) {
        if (!strncmp(buffer, "ATOM", 4))(*nfix)++;
    }

    rewind(fp);
    *fix = MYCALLOC(*fix, *nfix, sizeof(size_t));
    //fp = MYFOPEN(fp, ffile, "r");
    int na = 0;

    while (fgets(buffer, linesz - 1, fp) != NULL) {
        if (!strncmp(buffer, "ATOM", 4)) {
            (*fix)[na] = atoi(buffer + 4) - 1;
            na++;
        }
    }

    free(buffer);
    fclose(fp);
}

struct springset_pairs *read_springset_pairs(struct mol_atom_group *ag, char *sfile) {
    FILE *fp = MYFOPEN(fp, sfile, "r");

    struct springset_pairs *sprst;
    sprst = MYCALLOC(sprst, 1, sizeof(struct springset_pairs));

    int c;
    if (fscanf(fp, "%i", &sprst->nsprings) != 1) {
        ERR_MSG("Wrong spring file format\n");
    }

    sprst->springs = MYCALLOC(sprst->springs, sprst->nsprings, sizeof(struct pair_spring));
    struct pair_spring *sprs = sprst->springs;

    int id = 0;
    char name1[8], name2[8];
    int aid1, aid2;
    while (id < sprst->nsprings) {
        //             id1 n1 id2 n2 len err fk
        c = fscanf(fp, "%lf %lf %lf %i %s %i %s", &sprs[id].lnspr, &sprs[id].erspr,
                   &sprs[id].fkspr, &aid1, name1, &aid2, name2);

        //printf("%lf %lf %lf %i %s %i %s\n", sprs[id].lnspr, sprs[id].erspr, sprs[id].fkspr, aid1, name1, aid2, name2);
        if (c != 7) {
            ERR_MSG("Wrong spring file format\n");
        }

        sprs[id].ag = ag;
        //printf("%s %s\n", ag->atom_name[aid1-1], ag->atom_name[aid2-1]);
        if (strstr(ag->atom_name[aid1 - 1], name1) == NULL || strstr(ag->atom_name[aid2 - 1], name2) == NULL) {
            ERR_MSG("Inconsistent numbering in file %s\n", sfile);
        }

        sprs[id].laspr[0] = aid1 - 1;
        sprs[id].laspr[1] = aid2 - 1;

        id++;
    }

    fclose(fp);
    return sprst;
}

struct springset_points *read_springset_points(struct mol_atom_group *ag, char *sfile) {
    FILE *fp = MYFOPEN(fp, sfile, "r");

    struct springset_points *sprst;
    sprst = MYCALLOC(sprst, 1, sizeof(struct springset_points))

    int c;
    if (fscanf(fp, "%i", &sprst->nsprings) != 1) {
        ERR_MSG("Wrong spring file format\n");
    }

    sprst->springs = MYCALLOC(sprst->springs, sprst->nsprings, sizeof(struct point_spring));
    struct point_spring *sprs = sprst->springs;

    int id = 0;
    char name[8];
    int aid, naspr;
    double fkspr, X0, Y0, Z0;

    while (id < sprst->nsprings) {
        c = fscanf(fp, "%i %lf %lf %lf %lf", &naspr, &fkspr, &X0, &Y0, &Z0);
        if (c != 5) {
            ERR_MSG("Wrong spring file format\n");
        }

        sprs[id].naspr = naspr;
        sprs[id].fkspr = fkspr;
        sprs[id].X0 = X0;
        sprs[id].Y0 = Y0;
        sprs[id].Z0 = Z0;
        sprs[id].laspr = MYCALLOC(sprs[id].laspr, naspr, sizeof(int));
        sprs[id].ag = ag;

        for (int i = 0; i < naspr; i++) {
            c = fscanf(fp, "%i %s", &aid, name);
            if (c != 2) {
                ERR_MSG("Wrong spring file format\n");
            }

            if (strstr(ag->atom_name[aid - 1], name) == NULL) {
                ERR_MSG("Inconsistent numbering in file %s\n"
                        "Provided atom id %i has name %s instead of %s\n",
                        sfile, aid, ag->atom_name[aid - 1], name);
            }

            sprs[id].laspr[i] = aid - 1;
        }

        id++;
    }

    fclose(fp);
    return sprst;
}

void free_springset_pairs(struct springset_pairs **sprst) {
    if (*sprst != NULL) {
        free((*sprst)->springs);
        free(*sprst);
        *sprst = NULL;
    }
}

void free_springset_points(struct springset_points **sprst) {
    if (*sprst != NULL) {
        for (int i = 0; i < (*sprst)->nsprings; i++) {
            free((*sprst)->springs[i].laspr);
        }
        free((*sprst)->springs);
        free(*sprst);
        *sprst = NULL;
    }
}

void springeng_point(struct springset_points *sprst, double *een) {
    int i, i1, i2, nat;
    double xtot, ytot, ztot, fk;
    struct mol_vector3 g;
    struct mol_atom_group *ag;

    for (i = 0; i < sprst->nsprings; i++) {
        nat = sprst->springs[i].naspr;

        if (nat > 0) {
            xtot = 0.0;
            ytot = 0.0;
            ztot = 0.0;
            ag = sprst->springs[i].ag;

            for (i1 = 0; i1 < nat; i1++) {
                i2 = sprst->springs[i].laspr[i1];
                xtot += ag->coords[i2].X;
                ytot += ag->coords[i2].Y;
                ztot += ag->coords[i2].Z;
            }

            xtot = xtot / nat - sprst->springs[i].X0;
            ytot = ytot / nat - sprst->springs[i].Y0;
            ztot = ztot / nat - sprst->springs[i].Z0;

            fk = sprst->springs[i].fkspr;
            (*een) += fk * (xtot * xtot + ytot * ytot + ztot * ztot);

            fk = 2 * fk / nat;
            g.X = xtot * fk;
            g.Y = ytot * fk;
            g.Z = ztot * fk;

            for (i1 = 0; i1 < nat; i1++) {
                i2 = sprst->springs[i].laspr[i1];
                MOL_VEC_SUB(ag->gradients[i2], ag->gradients[i2], g);
            }
        }
    }
}

void springeng_pair(struct springset_pairs *sprst, double *een) {
    int i, i1, i2;
    double xtot, ytot, ztot, fk, d, d2, ln, er, coef, delta;
    struct mol_vector3 g;
    struct mol_atom_group *ag;

    for (i = 0; i < sprst->nsprings; i++) {
        ag = sprst->springs[i].ag;

        ln = sprst->springs[i].lnspr;
        er = sprst->springs[i].erspr;
        fk = sprst->springs[i].fkspr / 2.0;

        i1 = sprst->springs[i].laspr[0];
        i2 = sprst->springs[i].laspr[1];

        xtot = ag->coords[i2].X - ag->coords[i1].X;
        ytot = ag->coords[i2].Y - ag->coords[i1].Y;
        ztot = ag->coords[i2].Z - ag->coords[i1].Z;

        d2 = xtot * xtot + ytot * ytot + ztot * ztot;
        d = sqrt(d2);

        delta = fabs(d - ln);
        delta = (delta > er) ? ((delta - er) * delta / (d - ln)) : 0.0;

        //(*een) += fk * (d - ln) * (d - ln);
        (*een) += fk * delta * delta;
        //coef = fk * 2 * (1.0 - ln / d);
        coef = fk * 2.0 * delta / d;

        g.X = -coef * xtot;
        g.Y = -coef * ytot;
        g.Z = -coef * ztot;

        MOL_VEC_SUB(ag->gradients[i1], ag->gradients[i1], g);
        MOL_VEC_ADD(ag->gradients[i2], ag->gradients[i2], g);
    }
}

void usage_message(char **argv) {
    printf("\nUsage %s [ prm ] [ rtf ] [ -pdb ] [ -psf ]\n"
           "         [ -rec-pdb ] [ -rec-psf  ] [ -lig-pdb  ] [ -lig-psf ]\n"
           "         [ -lig-json ] [ -protocol ] [ -gbsa     ] [ -verbose ]\n"
           "         [ -help    ] [ -nsteps   ] [ -out      ]\n\n", argv[0]);

    help_message();
}

void help_message(void) {
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
           "\t--help ----------------- This message\n\n"

           "Prm and rft are required. One of the following must be provided:\n\n"
           "\t * pdb and psf or json of the full molecule (--pdb --psf)\n"
           "\t * OR ((rec-pdb + rec-psf) or rec-json) and ((lig-pdb + lig-psf) or lig-json)\n\n"

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
           "\tCorellation time (in nanoseconds): ex. 0.1\n"
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

