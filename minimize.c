#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <jansson.h>

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
    fprintf(stderr,"[Info]: " fmt, ##__VA_ARGS__);            \
} while(0)

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
    struct pairsprings_setup *sprst_pairs;
    struct pairsprings_js_setup *sprst_js_pairs;
    struct pointsprings_setup *sprst_points;
    struct noe_setup *nmr;
    struct density_setup *fit_prms;
    bool no_geom;
    bool bonds;
    bool angles;
    bool torsions;
    bool vdw;
};

struct pairspring {
    int laspr[2];      /**< list of atoms */
    double lnspr;
    double erspr;
    double fkspr;      /**< force constant */
};

struct pointspring {
    int naspr;      /**< number of affected atoms */
    int *laspr;      /**< list of atoms */
    double fkspr;      /**< force constant */
    double X0, Y0, Z0; /**< anchor point */
};

struct pairsprings_setup {
    int nsprings;       /**< number of springs */
    struct pairspring *springs;  /**< array of springs */
};

struct pairsprings_js_setup {
    int nsprings;       /**< number of springs */
    json_t *springs;  /**< json array of springs */
};

struct pointsprings_setup {
    int nsprings;       /**< number of springs */
    struct pointspring *springs;  /**< array of springs */
};

struct noe_setup {
    struct nmr_noe *spec;
    double weight;
    bool print_noe_matrix;
};

struct density_setup {
    struct mol_atom_group **ag_list;
    int ag_count;
    double weight;
    struct mol_fitting_params prms;
};

void fixed_atoms_read(char *ffile, size_t *nfix, size_t **fix);

struct pairsprings_setup *pairsprings_setup_read(struct mol_atom_group *ag, char *sfile);

struct pairsprings_js_setup *pairsprings_js_setup_read(struct mol_atom_group *ag, char *sfile);

void _pairspring_js_energy(struct pairsprings_js_setup *sprst, struct mol_atom_group *ag, double *een);

void pairsprings_js_setup_free(struct pairsprings_js_setup **sprst);

struct pointsprings_setup *pointsprings_setup_read(struct mol_atom_group *ag, char *sfile);

struct noe_setup *noe_setup_read(struct mol_atom_group *ag, char *sfile);

struct density_setup *density_setup_create(char *fitting_pdblist, double weight, double radius);

void noe_setup_free(struct noe_setup *nmr);

void density_setup_free(struct density_setup *prms);

void pairsprings_setup_free(struct pairsprings_setup **sprst);

void pointsprings_setup_free(struct pointsprings_setup **sprst);

void _pairspring_energy(struct pairsprings_setup *sprst, struct mol_atom_group *ag, double *een);

void _pointspring_energy(struct pointsprings_setup *sprst, struct mol_atom_group *ag, double *een);

static lbfgsfloatval_t energy_func(
        void *restrict prm,
        const double *restrict array,
        double *restrict gradient,
        const int array_size,
        const lbfgsfloatval_t step);

static void fprint_energy_terms(FILE *stream, void *restrict prm, char *prefix, json_t *json_model_dict);

void usage_message(char **argv);

void help_message(void);

struct mol_atom_group_list *read_ag_list(char *prm, char *rtf, char *pdb, char *psf, char *json, int score_only);

struct mol_atom_group_list* merge_ag_lists(struct mol_atom_group_list* ag1, struct mol_atom_group_list* ag2);

FILE* _fopen_err(char* file, char* mode) {
    FILE* f = fopen(file, mode);
    if (f == NULL) {
        ERR_MSG("Cannot open file %s\n", file);
    }
    return f;
}

int main(int argc, char **argv) {
    mol_enable_floating_point_exceptions();

    //static int verbose_flag = 0;
    static int ace_flag = 0;
    static int torsions_off = 0;
    static int print_noe_matrix = 0;
    static int score_only = 0;
    static int fix_rec = 0;
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
    char *pairsprings_json = NULL;
    char *pointsprings = NULL;
    char *fixed_pdb = NULL;
    char *fitting_pdblist = NULL;
    char *json_log = NULL;
    int nsteps = 1000;
    static int nonpar_terms_only = 0; // if geometry is not provided

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
                    {"pairsprings-json", required_argument, 0,                 0},
                    {"pointsprings",     required_argument, 0,                 0},
                    {"fitting-pdblist",  required_argument, 0,                 0},
                    {"fitting-weight",   required_argument, 0,                 0},
                    {"fixed-pdb",        required_argument, 0,                 0},
                    {"json-log",        required_argument, 0,                 0},
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
        MAGIC_ARGS(pairsprings_json);
        MAGIC_ARGS(pointsprings);
        MAGIC_ARGS(fitting_pdblist);
        MAGIC_ARGS(fixed_pdb);
        MAGIC_ARGS(json_log);

        if (strcmp("nsteps", long_options[option_index].name) == 0) {
            nsteps = atoi(optarg);
        }
        if (strcmp("fitting-weight", long_options[option_index].name) == 0) {
            fitting_weight = atof(optarg);
        }
    }

    struct mol_atom_group_list *aglist = NULL;

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

    if (pdb != NULL || json != NULL) {
        aglist = read_ag_list(prm, rtf, pdb, psf, json, score_only);
    } else {
        struct mol_atom_group_list *rec_aglist, *lig_aglist;
        rec_aglist = read_ag_list(prm, rtf, rec_pdb, rec_psf, rec_json, score_only);
        lig_aglist = read_ag_list(prm, rtf, lig_pdb, lig_psf, lig_json, score_only);
        aglist = merge_ag_lists(rec_aglist, lig_aglist);
        mol_atom_group_list_free(rec_aglist);
        mol_atom_group_list_free(lig_aglist);
    }
    INFO_MSG("Finished creating atom group\n");

    // Fill minimization parameters
    INFO_MSG("Filling energy parameters\n");
    struct energy_prm engpar;
    engpar.no_geom = false;
    engpar.bonds = true;
    engpar.angles = true;
    engpar.torsions = true;
    engpar.vdw = true;

    // Torsions
    if (torsions_off == 1) {
        INFO_MSG("Torsions OFF\n");
        engpar.torsions = false;
    }

    // If geometry is not provided
    if (nonpar_terms_only != 0) {
        INFO_MSG("Geometry OFF\n");
        engpar.no_geom = true;
        engpar.bonds = false;
        engpar.angles = false;
        engpar.torsions = false;
        engpar.vdw = false;
    }

    // NMR 2D spectrum
    if (noe_params != NULL) {
        INFO_MSG("NMR ON\n");
        engpar.nmr = noe_setup_read(&aglist->members[0], noe_params);
        engpar.nmr->print_noe_matrix = print_noe_matrix;
    } else {
        engpar.nmr = NULL;
    }

    // Pairsprings
    if (pairsprings != NULL) {
        INFO_MSG("Pairsprings ON\n");
        engpar.sprst_pairs = pairsprings_setup_read(&aglist->members[0], pairsprings);
    } else {
        engpar.sprst_pairs = NULL;
    }
    
    // pairsprings_json
    if (pairsprings_json != NULL) {
        INFO_MSG("Pairsprings-Json ON\n");
        engpar.sprst_js_pairs = pairsprings_js_setup_read(&aglist->members[0], pairsprings_json);
    } else {
        engpar.sprst_js_pairs = NULL;
    }

    // Point springs
    if (pointsprings != NULL) {
        INFO_MSG("Pointsprings ON\n");
        engpar.sprst_points = pointsprings_setup_read(&aglist->members[0], pointsprings);
    } else {
        engpar.sprst_points = NULL;
    };

    // Density fitting
    if (fitting_pdblist != NULL) {
        INFO_MSG("Density ON\n");
        engpar.fit_prms = density_setup_create(fitting_pdblist, fitting_weight, 2.0);
    } else {
        engpar.fit_prms = NULL;
    }

    // Fixed part
    size_t nfix_glob = 0;
    size_t *fix_glob = NULL;
    if (fixed_pdb != NULL) {
        INFO_MSG("Fixed atoms ON\n");
        fixed_atoms_read(fixed_pdb, &nfix_glob, &fix_glob);
    }

    // Init json array to record energy terms fro each model
    json_t *json_log_root = NULL;
    if (json_log != NULL) {
        INFO_MSG("Energy terms to JSON file ON\n");
        json_log_root = json_array();
    }

    // Start main loop
    FILE *outfile = NULL;
    if (score_only == 0) {
        outfile = _fopen_err(out, "w");
    }

    INFO_MSG("Started the main loop\n");
    for (int modeli = 0; modeli < aglist->size; modeli++) {
        INFO_MSG("Started model %i\n", modeli);

        if (aglist->size > 1) {
            if (outfile != NULL) {
                fprintf(outfile, "MODEL %i\n", (modeli + 1));
            }
        }

        struct mol_atom_group *ag = &aglist->members[modeli];
        engpar.ag = ag;

        // Init json log for current model
        json_t *json_start_dict = NULL, *json_final_dict = NULL, *json_model_dict = NULL;
        if (json_log_root != NULL) {
            json_model_dict = json_object();
            json_array_append(json_log_root, json_model_dict);
        }

        struct agsetup ags;
        struct acesetup ace_setup;

        // If no parameters was provided skip atom group setup
        if (!engpar.no_geom) {
            ag->gradients = calloc(ag->natoms, sizeof(struct mol_vector3));

            // Setup fixed atoms and lists
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

            // Setup GBSA
            engpar.ag_setup = &ags;
            if (ace_flag == 1) {
                engpar.ace_setup = &ace_setup;
            } else {
                engpar.ace_setup = NULL;
            }
        }

        // If not score_only, enter minimization function
        if (score_only == 0) {
            if (protocol != NULL) {
                FILE *prot_file = _fopen_err(protocol, "r");
                int cur_nsteps;
                char fix_path[1024];
                char spr_pair_path[1024];
                char spr_point_path[1024];
                char spr_js_pair_path[1024];

                int c;
                while ((c = fscanf(prot_file, "%i %s %s %s %s", &cur_nsteps, fix_path,
                                   spr_pair_path, spr_point_path, spr_js_pair_path)) != EOF) {
                    if (c != 5) {
                        ERR_MSG("Wrong protocol file format (%i words read)\n", c);
                    }

                    if (cur_nsteps < 0) {
                        ERR_MSG("Wrong protocol file format (number of steps must be non-negative)\n");
                    }

                    size_t nfix = 0;
                    size_t *fix = NULL;

                    // Setup fixed atoms
                    if (strcmp(fix_path, ".") != 0) {
                        fixed_atoms_read(fix_path, &nfix, &fix);
                        mol_fixed_update(ag, nfix, fix);
                        update_nblst(ag, &ags);

                        if (ace_flag == 1) {
                            ace_fixedupdate(ag, &ags, &ace_setup);
                            ace_updatenblst(&ags, &ace_setup);
                        }

                        free(fix);
                        fix = NULL;
                    }

                    // Setup springs
                    struct pairsprings_setup *sprst_pairs = NULL;
                    struct pointsprings_setup *sprst_points = NULL;

                    if (strcmp(spr_pair_path, ".") != 0) {
                        sprst_pairs = pairsprings_setup_read(ag, spr_pair_path);
                    }

                    if (strcmp(spr_point_path, ".") != 0) {
                        sprst_points = pointsprings_setup_read(ag, spr_point_path);
                    }
                    
                    struct pairsprings_js_setup *sprst_js_pairs=NULL;
                    if (strcmp(spr_js_pair_path, ".") != 0) {
                        sprst_js_pairs = pairsprings_js_setup_read(ag, spr_js_pair_path);
                    }

                    engpar.ag_setup = &ags;
                    engpar.sprst_pairs = sprst_pairs;
                    engpar.sprst_points = sprst_points;
                    engpar.sprst_js_pairs = sprst_js_pairs;

                    // GBSA
                    if (ace_flag == 1) {
                        engpar.ace_setup = &ace_setup;
                    } else {
                        engpar.ace_setup = NULL;
                    }

                    // Record starting energy terms
                    if (json_model_dict != NULL) {
                        json_start_dict = json_object();
                        json_object_set(json_model_dict, "START", json_start_dict);
                    }
                    fprint_energy_terms(outfile, &engpar, "REMARK START ", json_start_dict);

                    // Minimize energy
                    if (cur_nsteps > 0) {
                        mol_minimize_ag(MOL_LBFGS, cur_nsteps, __TOL__, ag, (void *) (&engpar), energy_func);
                    }

                    pairsprings_setup_free(&sprst_pairs);
                    pointsprings_setup_free(&sprst_points);
                    pairsprings_js_setup_free(&sprst_js_pairs);
                }
            } else {
                // Record starting energy terms
                if (json_model_dict != NULL) {
                    json_start_dict = json_object();
                    json_object_set(json_model_dict, "START", json_start_dict);
                }
                fprint_energy_terms(outfile, &engpar, "REMARK START ", json_start_dict);

                // Minimize energy
                if (nsteps > 0) {
                    mol_minimize_ag(MOL_LBFGS, nsteps, __TOL__, ag, (void *) (&engpar), energy_func);
                }
            }
        }

        // Record final energy terms
        if (json_model_dict != NULL) {
            json_final_dict = json_object();
            json_object_set(json_model_dict, "FINAL", json_final_dict);
        }
        if (outfile != NULL) {
            fprintf(outfile, "REMARK\n");
        }
        fprint_energy_terms(outfile, &engpar, "REMARK FINAL ", json_final_dict);

        // Write final model
        if (outfile != NULL) {
            mol_fwrite_pdb(outfile, ag);
            if (aglist->size > 1) {
                fprintf(outfile, "ENDMDL\n");
            }
            fflush(outfile);
        }

        if (ags.nblst != NULL) {
            destroy_agsetup(&ags);
        }
    }
    INFO_MSG("Finished minimization loop\n");

    // Free everything
    noe_setup_free(engpar.nmr);
    density_setup_free(engpar.fit_prms);

    if (engpar.sprst_points != NULL) {
        pointsprings_setup_free(&engpar.sprst_points);
    }
    if (engpar.sprst_js_pairs != NULL) {
        pairsprings_js_setup_free(&engpar.sprst_js_pairs);
    }
    if (engpar.sprst_pairs != NULL) {
        pairsprings_setup_free(&engpar.sprst_pairs);
    }
    if (fix_glob != NULL) {
        free(fix_glob);
    }
    if (json_log_root != NULL) {
        json_dump_file(json_log_root, json_log, JSON_INDENT(4));
        json_decref(json_log_root);
    }

    if (outfile != NULL) {
        fclose(outfile);
    }
    mol_atom_group_list_free(aglist);

    INFO_MSG("Completed\n");
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
        _pairspring_energy(energy_prm->sprst_pairs, energy_prm->ag, &energy);
    }
    
    if (energy_prm->sprst_js_pairs != NULL) {
        _pairspring_js_energy(energy_prm->sprst_js_pairs, energy_prm->ag, &energy);
    }

    if (energy_prm->sprst_points != NULL) {
        _pointspring_energy(energy_prm->sprst_points, energy_prm->ag, &energy);
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

static void fprint_energy_terms(FILE *stream, void *restrict prm, char *prefix, json_t *json_model_dict) {
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
            if (stream != NULL) {
                fprintf(stream, strcat(fmt, "ACE: % .3f\n"), energy);
            }
            if (json_model_dict != NULL) {
                json_object_set_new(json_model_dict, "ACE", json_real(energy));
            }
            total += energy;
            energy = 0.0;
        }

        if (energy_prm->vdw) {
            vdweng(energy_prm->ag, &energy, energy_prm->ag_setup->nblst);
            strcpy(fmt, prefix);
            if (stream != NULL) {
                fprintf(stream, strcat(fmt, "VDW: % .3f\n"), energy);
            }
            if (json_model_dict != NULL) {
                json_object_set_new(json_model_dict, "VDW", json_real(energy));
            }
            total += energy;
            energy = 0.0;

            vdwengs03(1.0, energy_prm->ag_setup->nblst->nbcof, energy_prm->ag, &energy,
                      energy_prm->ag_setup->nf03, energy_prm->ag_setup->listf03);
            strcpy(fmt, prefix);
            if (stream != NULL) {
                fprintf(stream, strcat(fmt, "VDW03: % .3f\n"), energy);
            }
            if (json_model_dict != NULL) {
                json_object_set_new(json_model_dict, "VDW03", json_real(energy));
            }
            total += energy;
            energy = 0.0;
        }

        if (energy_prm->bonds) {
            beng(energy_prm->ag, &energy);
            strcpy(fmt, prefix);
            if (stream != NULL) {
                fprintf(stream, strcat(fmt, "Bonds % .3f\n"), energy);
            }
            if (json_model_dict != NULL) {
                json_object_set_new(json_model_dict, "Bonds", json_real(energy));
            }
            total += energy;
            energy = 0.0;
        }

        if (energy_prm->angles) {
            aeng(energy_prm->ag, &energy);
            strcpy(fmt, prefix);
            if (stream != NULL) {
                fprintf(stream, strcat(fmt, "Angles: % .3f\n"), energy);
            }
            if (json_model_dict != NULL) {
                json_object_set_new(json_model_dict, "Angles", json_real(energy));
            }
            total += energy;
            energy = 0.0;
        }

        if (energy_prm->torsions) {
            teng(energy_prm->ag, &energy);
            strcpy(fmt, prefix);
            if (stream != NULL) {
                fprintf(stream, strcat(fmt, "Torsions: % .3f\n"), energy);
            }
            if (json_model_dict != NULL) {
                json_object_set_new(json_model_dict, "Torsions", json_real(energy));
            }
            total += energy;
            energy = 0.0;

            ieng(energy_prm->ag, &energy);
            strcpy(fmt, prefix);
            if (stream != NULL) {
                fprintf(stream, strcat(fmt, "Impropers: % .3f\n"), energy);
            }
            if (json_model_dict != NULL) {
                json_object_set_new(json_model_dict, "Impropers", json_real(energy));
            }
            total += energy;
            energy = 0.0;
        }
    }

    if (energy_prm->sprst_pairs != NULL) {
        _pairspring_energy(energy_prm->sprst_pairs, energy_prm->ag, &energy);
        strcpy(fmt, prefix);
        if (stream != NULL) {
            fprintf(stream, strcat(fmt, "Pairsprings: % .3f\n"), energy);
        }
        if (json_model_dict != NULL) {
            json_object_set_new(json_model_dict, "Pairsprings", json_real(energy));
        }
        total += energy;
        energy = 0.0;
    }
    
    if (energy_prm->sprst_js_pairs != NULL) {
        _pairspring_js_energy(energy_prm->sprst_js_pairs, energy_prm->ag, &energy);
        strcpy(fmt, prefix);
        if (stream != NULL) {
            fprintf(stream, strcat(fmt, "Pairsprings-Json: % .3f\n"), energy);
        }
        if (json_model_dict != NULL) {
            json_object_set_new(json_model_dict, "Pairsprings-Json", json_real(energy));
        }
        total += energy;
        energy = 0.0;
    }
    
    if (energy_prm->sprst_points != NULL) {
        _pointspring_energy(energy_prm->sprst_points, energy_prm->ag, &energy);
        strcpy(fmt, prefix);
        if (stream != NULL) {
            fprintf(stream, strcat(fmt, "Pointsprings: % .3f\n"), energy);
        }
        if (json_model_dict != NULL) {
            json_object_set_new(json_model_dict, "Pointsprings", json_real(energy));
        }
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
        if (stream != NULL) {
            fprintf(stream, strcat(fmt, "NOE: % .6f\n"), energy);
        }
        if (json_model_dict != NULL) {
            json_object_set_new(json_model_dict, "NOE", json_real(energy));
        }
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
        if (stream != NULL) {
            fprintf(stream, strcat(fmt, "Density: % .4f\n"), energy);
        }
        if (json_model_dict != NULL) {
            json_object_set_new(json_model_dict, "Density", json_real(energy));
        }
        total += energy;
        energy = 0.0;
    }

    strcpy(fmt, prefix);
    if (stream != NULL) {
        fprintf(stream, strcat(fmt, "Total: % .3f\n"), total);
    }
    if (json_model_dict != NULL) {
        json_object_set_new(json_model_dict, "Total", json_real(total));
    }
}




/*
 * Create mol_atom_group_list with geometry (if score_only == 1 can be without)
 */
struct mol_atom_group_list *read_ag_list(char *prm, char *rtf, char *pdb, char *psf, char *json, int score_only) {
    INFO_MSG("Started reading a new atom group\n");
    struct mol_atom_group_list* ag_list = NULL;
    struct mol_atom_group* ag_json = NULL;

    if (json != NULL) {
        INFO_MSG("Reading json file %s\n", json);
        ag_json = mol_read_json(json);

        if (psf != NULL) {
            WRN_MSG("Parameter --psf is not effective, since --json was provided\n");
        }
    }

    if (pdb != NULL) {
        INFO_MSG("Reading %s\n", pdb);
        INFO_MSG("Trying to read %s as a multimodel one (with MODEL records)\n", pdb);
        ag_list = mol_read_pdb_models(pdb);

        if (ag_list == NULL) {
            INFO_MSG("File %s doesn't have MODEL records. Reading as a regular pdb file\n", pdb);
            struct mol_atom_group* ag = mol_read_pdb(pdb);
            if (ag == NULL) {
                ERR_MSG("Failed reading %s", pdb);
            }

            ag_list = mol_atom_group_list_create(1);
            ag_list->members[0] = *ag;
            free(ag);
        }

        if (ag_json != NULL) {
            INFO_MSG("Reading coordinates from %s and geometry from %s\n", pdb, json);
            // Copy geometry from json to aglist
            struct mol_atom_group_list *fin_aglist = mol_atom_group_list_create(ag_list->size);

            for (size_t i = 0; i < fin_aglist->size; i++) {
                struct mol_atom_group* tmp_copy = mol_atom_group_copy(ag_json);
                fin_aglist->members[i] = *tmp_copy;
                free(tmp_copy); // Do I need this?

                if (fin_aglist->members[i].natoms != ag_list->members[i].natoms) {
                    ERR_MSG("Model %i in %s and %s have different atom numbers (%i, %i)",
                            (int) i, pdb, json,
                            (int) ag_list->members[i].natoms,
                            (int) fin_aglist->members[i].natoms);
                }
                memcpy(fin_aglist->members[i].coords,
                       ag_list->members[i].coords,
                       sizeof(struct mol_vector3) * ag_list->members[i].natoms);
            }
            mol_atom_group_list_free(ag_list);
            mol_atom_group_free(ag_json);
            ag_list = fin_aglist;

        } else if (psf != NULL) {
            INFO_MSG("Reading geometry from %s\n", psf);
            for (size_t i = 0; i < ag_list->size; i++) {
                mol_atom_group_read_geometry(&ag_list->members[i], psf, prm, rtf);
            }
        } else if (score_only != 0) {
            WRN_MSG("Geometry for the molecule is not provided, computing only non-parametric terms\n");
        } else {
            ERR_MSG("--json or --psf must be provided with --pdb");
        }
    } else if (ag_json != NULL) {
        INFO_MSG("Using geometry and coordinates from %s\n", json);
        ag_list = mol_atom_group_list_create(1);
        ag_list->members[0] = *ag_json;
        free(ag_json);

    } else {
        ERR_MSG("Json or pdb file must be provided");
    }

    return ag_list;
}


struct mol_atom_group_list* merge_ag_lists(struct mol_atom_group_list* ag1, struct mol_atom_group_list* ag2) {
    if (ag1->size != ag2->size) {
        ERR_MSG("Receptor and ligand have different numbers of models (%i, %i)", (int)ag1->size, (int)ag2->size);
    }

    INFO_MSG("Using models assembled from receptor and ligand models provided separately\n");
    struct mol_atom_group_list *ag_list = mol_atom_group_list_create(ag1->size);
    for (size_t i = 0; i < ag_list->size; i++) {
        struct mol_atom_group* _join = mol_atom_group_join(&ag1->members[i], &ag2->members[i]);
        ag_list->members[i] = *_join;
        free(_join);
        //ag_list->members[i] = *mol_atom_group_join(&ag1->members[i], &ag2->members[i]);
    }

    return ag_list;
}


struct density_setup *density_setup_create(char *fitting_pdblist, double weight, double radius) {
    struct density_setup *fit_prms = calloc(1, sizeof(struct density_setup));

    FILE *f = _fopen_err(fitting_pdblist, "r");
    char pdb_file[512];
    struct mol_atom_group **aglist = calloc(1, sizeof(struct mol_atom_group *));
    int ag_count = 0;
    int cur_size = 1;
    while (fgets(pdb_file, 512, f) != NULL) {
        // Remove new line character
        for (int i = 0; i < 512; i++) {
            if (pdb_file[i] == '\n') {
                pdb_file[i] = '\0';
                break;
            }
        }

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


void density_setup_free(struct density_setup *prms) {
    if (prms != NULL) {
        for (int i = 0; i < prms->ag_count; i++) {
            mol_atom_group_free(prms->ag_list[i]);
        }
        free(prms->ag_list);
        free(prms);
    }
}


struct noe_setup *noe_setup_read(struct mol_atom_group *ag, char *sfile) {
    struct noe_setup *nmr = calloc(1, sizeof(struct noe_setup));

    char line[512];
    FILE *f = _fopen_err(sfile, "r");

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

void noe_setup_free(struct noe_setup *nmr) {
    if (nmr != NULL) {
        nmr_noe_free(nmr->spec);
        free(nmr);
    }
}

void fixed_atoms_read(char *ffile, size_t *nfix, size_t **fix) {
    int linesz = 91;
    char *buffer = calloc(linesz, sizeof(char));

    *nfix = 0;
    FILE *fp = _fopen_err(ffile, "r");

    while (fgets(buffer, linesz - 1, fp) != NULL) {
        if (!strncmp(buffer, "ATOM", 4))(*nfix)++;
    }

    rewind(fp);
    *fix = calloc(*nfix, sizeof(size_t));
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

struct pairsprings_setup *pairsprings_setup_read(struct mol_atom_group *ag, char *sfile) {
    FILE *fp = _fopen_err(sfile, "r");

    struct pairsprings_setup *sprst;
    sprst = calloc(1, sizeof(struct pairsprings_setup));

    int c;
    if (fscanf(fp, "%i", &sprst->nsprings) != 1) {
        ERR_MSG("Wrong spring file format\n");
    }

    sprst->springs = calloc(sprst->nsprings, sizeof(struct pairspring));
    struct pairspring *sprs = sprst->springs;

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

        if (strstr(ag->atom_name[aid1-1], name1) == NULL) {
            ERR_MSG("Inconsistent numbering in file %s\n"
                    "Provided atom %i has name (%i, %s, %s) instead of %s\n",
                    sfile, aid1, ag->residue_id[aid1-1].residue_seq, ag->residue_name[aid1-1], ag->atom_name[aid1-1], name1);
        }
        if (strstr(ag->atom_name[aid2 - 1], name2) == NULL) {
            ERR_MSG("Inconsistent numbering in file %s\n"
                    "Provided atom %i has name (%i, %s, %s) instead of %s\n",
                    sfile, aid2, ag->residue_id[aid2-1].residue_seq, ag->residue_name[aid2-1], ag->atom_name[aid2-1], name2);
        }

        sprs[id].laspr[0] = aid1 - 1;
        sprs[id].laspr[1] = aid2 - 1;

        id++;
    }

    fclose(fp);
    return sprst;
}

struct pointsprings_setup *pointsprings_setup_read(struct mol_atom_group *ag, char *sfile) {
    FILE *fp = _fopen_err(sfile, "r");

    struct pointsprings_setup *sprst;
    sprst = calloc(1, sizeof(struct pointsprings_setup));

    int c;
    if (fscanf(fp, "%i", &sprst->nsprings) != 1) {
        ERR_MSG("Wrong spring file format\n");
    }

    sprst->springs = calloc(sprst->nsprings, sizeof(struct pointspring));
    struct pointspring *sprs = sprst->springs;

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
        sprs[id].laspr = calloc(naspr, sizeof(int));

        for (int i = 0; i < naspr; i++) {
            c = fscanf(fp, "%i %s", &aid, name);
            if (c != 2) {
                ERR_MSG("Wrong spring file format\n");
            }

            if (strstr(ag->atom_name[aid - 1], name) == NULL) {
                ERR_MSG("Inconsistent numbering in file %s\n"
                        "Provided atom %i has name (%i, %s, %s) instead of %s\n",
                        sfile, aid, ag->residue_id[aid-1].residue_seq, ag->residue_name[aid-1], ag->atom_name[aid-1], name);
            }

            sprs[id].laspr[i] = aid - 1;
        }

        id++;
    }

    fclose(fp);
    return sprst;
}

void pairsprings_setup_free(struct pairsprings_setup **sprst) {
    if (*sprst != NULL) {
        free((*sprst)->springs);
        free(*sprst);
        *sprst = NULL;
    }
}

void pointsprings_setup_free(struct pointsprings_setup **sprst) {
    if (*sprst != NULL) {
        for (int i = 0; i < (*sprst)->nsprings; i++) {
            free((*sprst)->springs[i].laspr);
        }
        free((*sprst)->springs);
        free(*sprst);
        *sprst = NULL;
    }
}

void _pointspring_energy(struct pointsprings_setup *sprst, struct mol_atom_group* ag, double *een) {
    int i, i1, i2, nat;
    double xtot, ytot, ztot, fk;
    struct mol_vector3 g;

    for (i = 0; i < sprst->nsprings; i++) {
        nat = sprst->springs[i].naspr;

        if (nat > 0) {
            xtot = 0.0;
            ytot = 0.0;
            ztot = 0.0;

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

void _pairspring_energy(struct pairsprings_setup *sprst, struct mol_atom_group* ag, double *een) {
    int i, i1, i2;
    double xtot, ytot, ztot, fk, d, d2, ln, er, coef, delta;
    struct mol_vector3 g;

    for (i = 0; i < sprst->nsprings; i++) {
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

void _pairspring_js_energy(struct pairsprings_js_setup *sprst, struct mol_atom_group* ag, double *een) {
    int i, j, k, idx2, idx1, m, n;
    size_t lni1, lni2;
    double xtot, ytot, ztot, coef, delta, d2, d, ln, ler, rer, fk, hk;
    json_t *i1, *i2;
    const char *potential, *averaging;
    struct mol_vector3 g;
    
    for (i = 0; i< sprst -> nsprings; i++) {
        json_t *line_1 = json_array_get(sprst->springs, i);
        ln = json_number_value(json_object_get(line_1, "distance"));
        ler = json_number_value(json_object_get(line_1, "lerror"));
        rer = json_number_value(json_object_get(line_1, "rerror"));
        fk = json_number_value(json_object_get(line_1, "weight"));

        // start to calculate the average distance over two groups
        i1 = json_object_get(line_1, "group1");
        i2 = json_object_get(line_1, "group2");
        
        lni1 = json_array_size(i1);
        lni2 = json_array_size(i2);
        double aved = 0;
        for(j = 0; j < lni1; j++) {
            idx1 = json_integer_value(json_array_get(i1,j));
            for (k = 0; k < lni2; k++) {
                idx2 = json_integer_value(json_array_get(i2,k));
                xtot = ag->coords[idx2].X - ag->coords[idx1].X;
                ytot = ag->coords[idx2].Y - ag->coords[idx1].Y;
                ztot = ag->coords[idx2].Z - ag->coords[idx1].Z;
                d2 = xtot*xtot + ytot*ytot + ztot*ztot;
                d = sqrt(d2);
                aved += 1/(pow(d, 6));
            }
        }
        averaging = json_string_value(json_object_get( line_1, "average"));
        if (!(strcmp(averaging, "SUM"))) {
            aved = 1/pow(aved,1.0/6);
        } else {  // R6 average
            aved = 1/pow(aved/lni1/lni2,1.0/6);
        }
        delta = aved - ln;
        potential = json_string_value(json_object_get( line_1, "potential"));

        if (!(strncmp(potential, "SQUARE-WELL", 4))) {  // square-well
            if (delta < 0) {
                delta = (delta < -ler) ? (delta+ler) : 0.0;
            } else {
                delta = (delta > rer) ? (delta-rer) : 0.0;
            }
            (*een) += fk * delta * delta;
            coef = fk * 2.0 * delta / d;
        }
        else if (!(strncmp(potential, "BIHARMONIC", 4))) {  // biharmonic, temperature=300
            if (delta < 0) {
                hk = fk * 300*0.0019872041/(2*ler*ler);
            } else {
                hk = fk * 300*0.0019872041/(2*rer*rer);
            }
            
            hk = fmin(1000, hk);
            (*een) += hk * delta * delta;
            coef = hk*2.0*delta/d;
        }
        else if (!(strncmp(potential, "SOFT-SQUARE", 4))) {  // soft-square
            if (delta < 0) {
                delta = (delta < -ler) ? (delta+ler) : 0.0;
            } else {
                delta = (delta > rer) ? (delta-rer) : 0.0;
            }
            
            if (delta <= 0.5) {  // here the switch bound is 0.5
                (*een) += fk* delta*delta;
                coef = fk*2.0*delta/d;
            }
            else {  //delta > 0.5
                (*een) += fk*(delta-0.25);
                coef = fk;
            }
        }
        else {
            ERR_MSG("The potential should be one of SQUARE-WEll, BIHARMONIC, SOFT-SQUARE.\n");
        }
        g.X = -coef * xtot;
        g.Y = -coef * ytot;
        g.Z = -coef * ztot;
        
        // how to add gradients to all the atoms
        for(m = 0; m < lni1; m++) {
            idx1 = json_integer_value(json_array_get(i1,m));
            MOL_VEC_SUB(ag->gradients[idx1], ag->gradients[idx1], g);
        }
        for(n = 0; n < lni2; n++) {
            idx2 = json_integer_value(json_array_get(i2,n));
            MOL_VEC_ADD(ag->gradients[idx2], ag->gradients[idx2], g);
        }
    }
}

struct pairsprings_js_setup *pairsprings_js_setup_read(struct mol_atom_group *ag, char *sfile) {
    FILE *fp = _fopen_err(sfile, "r");
    
    struct pairsprings_js_setup *sprst;
    sprst = calloc(1, sizeof(struct pairsprings_js_setup));
    sprst->springs = calloc(300, sizeof(struct pairspring));
    
    char c[200];
    json_error_t error;
    int id = 0;
    
    json_t *jsprings = json_array();
    json_t *spr_line;
    while (fgets(c, 200, fp)!=NULL) {
        spr_line = json_loads(c, 0, &error);
        if (spr_line==NULL) {
            ERR_MSG("Wrong spring file format, should be in json format.\n");
        } else if (json_is_array(spr_line)) {
            jsprings = spr_line;
            continue;
        }
        else {
            json_array_append(jsprings, spr_line);
        }
        
        // message to say the file format is wrong.
        if (json_array_get( json_object_get(spr_line, "group1"), 0)==NULL) {
            ERR_MSG("Wrong spring file format\n");
        }
        if (json_array_get( json_object_get(spr_line, "group2"), 0)==NULL) {
            ERR_MSG("Wrong spring file format\n");
        }
        if ( json_object_get(spr_line, "potential") ==NULL) {
            ERR_MSG("Wrong spring file format\n");
        }
        if ( json_object_get(spr_line, "weight") ==NULL) {
            ERR_MSG("Wrong spring file format\n");
        }
        if ( json_object_get(spr_line, "lerror") ==NULL) {
            ERR_MSG("Wrong spring file format\n");
        }
        if ( json_object_get(spr_line, "rerror") ==NULL) {
            ERR_MSG("Wrong spring file format\n");
        }
        if ( json_object_get(spr_line, "distance") ==NULL) {
            ERR_MSG("Wrong spring file format\n");
        }
        if ( json_object_get(spr_line, "average") ==NULL) {
            ERR_MSG("Wrong spring file format\n");
        }

        id++;
    }
    json_decref(spr_line);
    sprst->nsprings = id;
    sprst->springs = jsprings;
    fclose(fp);
    return sprst;
}

// free function may not be used
void pairsprings_js_setup_free(struct pairsprings_js_setup **sprst) {
    if (*sprst != NULL) {
        free((*sprst)->springs);
        free(*sprst);
        *sprst = NULL;
    }
}

void usage_message(char **argv) {
    printf("\nUsage %s [ options ]\n\n", argv[0]);
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
           "\t--pairsprings-json ----- Json File with pairwise distance restraints\n"
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

