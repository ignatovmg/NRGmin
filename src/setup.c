#include "utils.h"
#include "setup.h"
#include "parse_options.h"

#include "mol2/json.h"
#include "mol2/icharmm.h"


static struct mol_atom_group_list *_read_ag_list(
        const char *prm,
        const char *rtf,
        const char *pdb,
        const char *psf,
        const char *json,
        const int score_only) {

    DEBUG_MSG("Started reading a new atom group\n");
    struct mol_atom_group_list* ag_list = NULL;
    struct mol_atom_group* ag_json = NULL;

    // Read molecule from json
    if (json != NULL) {
        DEBUG_MSG("Reading json file %s", json);
        ag_json = mol_read_json(json);

        if (ag_json == NULL) {
            ERR_MSG("Failed reading %s", json);
            return NULL;
        }
    }

    // Read molecule from pdb
    if (pdb != NULL) {
        DEBUG_MSG("Reading %s", pdb);
        DEBUG_MSG("Trying to read %s as a multimodel one (with MODEL records)", pdb);
        ag_list = mol_read_pdb_models(pdb);

        // If not multimodel, try to read as a single model
        if (ag_list == NULL) {
            DEBUG_MSG("File %s doesn't have MODEL records. Reading as a regular pdb file", pdb);
            struct mol_atom_group* ag_single = mol_read_pdb(pdb);

            // If both single and multi-model failed - give up
            if (ag_single == NULL) {
                ERR_MSG("Failed reading %s", pdb);
                mol_atom_group_free(ag_json);
                return NULL;
            }

            ag_list = mol_atom_group_list_create(1);
            ag_list->members[0] = *ag_single;
            free(ag_single);
        }

        // If json was provied too then merge geometry from json and coords from pdb
        if (ag_json != NULL) {
            DEBUG_MSG("Reading coordinates from %s and geometry from %s", pdb, json);

            struct mol_atom_group_list *fin_aglist = mol_atom_group_list_create(ag_list->size);

            for (size_t i = 0; i < fin_aglist->size; i++) {
                struct mol_atom_group* tmp_copy = mol_atom_group_copy(ag_json);
                fin_aglist->members[i] = *tmp_copy;
                free(tmp_copy); // Do I need this?

                if (fin_aglist->members[i].natoms != ag_list->members[i].natoms) {
                    DEBUG_MSG("Model %i in %s and %s have different atom numbers (%i, %i)",
                            (int) i, pdb, json,
                            (int) ag_list->members[i].natoms,
                            (int) fin_aglist->members[i].natoms);
                    mol_atom_group_list_free(fin_aglist);
                    mol_atom_group_list_free(ag_list);
                    mol_atom_group_free(ag_json);
                    return NULL;
                }
                memcpy(fin_aglist->members[i].coords,
                        ag_list->members[i].coords,
                        sizeof(struct mol_vector3) * ag_list->members[i].natoms);
            }

            mol_atom_group_list_free(ag_list);
            mol_atom_group_free(ag_json);
            ag_list = fin_aglist;

        } else if (psf && prm && rtf) {
            // If json wasn't provided read geometry from prm rtf psf
            DEBUG_MSG("Reading geometry from %s", psf);
            for (size_t i = 0; i < ag_list->size; i++) {
                if (!mol_atom_group_read_geometry(&ag_list->members[i], psf, prm, rtf)) {
                    ERR_MSG("Couldn't fill geometry from psf rtf and prm");
                    mol_atom_group_list_free(ag_list);
                }
            }

        } else if (score_only != 0) {
            // if they weren't provided, check the score_only flag
            WRN_MSG("Geometry for the molecule is not provided, computing only non-parametric terms");

        } else {
            // Give up if nothing worked
            ERR_MSG("Force field and geometry are not provided");
            mol_atom_group_list_free(ag_list);
            return NULL;
        }

    } else if (ag_json != NULL) {
        DEBUG_MSG("Using geometry and coordinates from %s", json);
        ag_list = mol_atom_group_list_create(1);
        ag_list->members[0] = *ag_json;
        free(ag_json);

    } else {
        ERR_MSG("Json or pdb file must be provided");
        return NULL;
    }

    return ag_list;
}


static struct mol_atom_group_list* _merge_ag_lists(const struct mol_atom_group_list* ag1, const struct mol_atom_group_list* ag2) {
    if (ag1->size != ag2->size) {
        ERR_MSG("Receptor and ligand have different numbers of models (%zu, %zu)", ag1->size, ag2->size);
        return NULL;
    }

    DEBUG_MSG("Using models assembled from receptor and ligand models provided separately\n");
    struct mol_atom_group_list *ag_list = mol_atom_group_list_create(ag1->size);

    for (size_t i = 0; i < ag_list->size; i++) {
        struct mol_atom_group* _join = mol_atom_group_join(&ag1->members[i], &ag2->members[i]);
        ag_list->members[i] = *_join;
        free(_join);
    }

    return ag_list;
}


struct mol_atom_group_list* mol_atom_group_list_from_options(struct options *opts)
{
    struct mol_atom_group_list* ag_list;

    if (opts->separate) {
        struct mol_atom_group_list* rec_list = _read_ag_list(
                opts->rec_prm,
                opts->rec_rtf,
                opts->rec_pdb,
                opts->rec_psf,
                opts->rec_json,
                opts->score_only);
        struct mol_atom_group_list* lig_list = _read_ag_list(
                opts->lig_prm,
                opts->lig_rtf,
                opts->lig_pdb,
                opts->lig_psf,
                opts->lig_json,
                opts->score_only);

        if (!rec_list || !lig_list) {
            if (rec_list) {
                mol_atom_group_list_free(rec_list);
            }
            if (lig_list) {
                mol_atom_group_list_free(lig_list);
            }
            return NULL;
        }

        opts->rec_natoms = rec_list->members[0].natoms;
        opts->lig_natoms = lig_list->members[0].natoms;

        ag_list = _merge_ag_lists(rec_list, lig_list);

    } else {
        ag_list = _read_ag_list(
                opts->prm,
                opts->rtf,
                opts->pdb,
                opts->psf,
                opts->json,
                opts->score_only);

        opts->rec_natoms = 0;
        opts->lig_natoms = 0;
    }
    return ag_list;
}


static void _fixed_setup_free(struct fixed_setup **fixed)
{
    if (*fixed != NULL) {
        if ((*fixed)->atoms) {
            free((*fixed)->atoms);
        }
        free(*fixed);
        *fixed = NULL;
    }
}


static struct fixed_setup* _fixed_setup_read_txt(const char *file) {
    FILE* fp;
    FOPEN_ELSE(fp, file, "r") {
        return NULL;
    }

    size_t nfix = 0;

    int linesz = 91;
    char *buffer = calloc(linesz, sizeof(char));

    while (fgets(buffer, linesz - 1, fp) != NULL) {
        if (!strncmp(buffer, "ATOM", 4)) nfix++;
    }

    rewind(fp);
    size_t* fix = calloc(nfix, sizeof(size_t));
    size_t na = 0;

    while (fgets(buffer, linesz - 1, fp) != NULL) {
        if (!strncmp(buffer, "ATOM", 4)) {
            fix[na] = atoi(buffer + 4) - 1;
            na++;
        }
    }

    free(buffer);
    fclose(fp);

    struct fixed_setup* result = calloc(1, sizeof(struct fixed_setup));
    result->atoms = calloc(nfix, sizeof(size_t));
    result->natoms = nfix;

    return result;
}


static struct fixed_setup* _fixed_setup_read_json(const json_t *root)
{
    if (!json_is_array(root)) {
        ERR_MSG("Fixed atoms json setup must be an array");
        return NULL;
    }

    struct fixed_setup* result = calloc(1, sizeof(struct fixed_setup));
    result->natoms = json_array_size(root);
    result->atoms = calloc(result->natoms, sizeof(size_t));

    size_t counter;
    json_t* atom_id;
    json_array_foreach(root, counter, atom_id) {
        if (!json_is_integer(atom_id)) {
            _fixed_setup_free(&result);
            json_decref(atom_id);
            return NULL;
        }
        result->atoms[counter] = json_integer_value(atom_id);
    }

    return result;
}


static void _pairsprings_setup_free(struct pairsprings_setup **sprst) {
    if (*sprst != NULL) {
        free((*sprst)->springs);
        free(*sprst);
        *sprst = NULL;
    }
}


static struct fixed_setup* _fixed_setup_atom_range(const size_t start_atom, const size_t end_atom)
{
    struct fixed_setup* fixed = calloc(1, sizeof(struct fixed_setup));
    fixed->natoms = end_atom - start_atom;
    fixed->atoms = calloc(fixed->natoms, sizeof(size_t));
    for (size_t i = 0; i < fixed->natoms; i++) {
        fixed->atoms[i] = start_atom + i;
    }
    return fixed;
}


static struct pairsprings_setup *_pairsprings_setup_read_txt(const char *path) {
    FILE *fp;
    FOPEN_ELSE(fp, path, "r") {
        return NULL;
    }

    struct pairsprings_setup *sprst;
    sprst = calloc(1, sizeof(struct pairsprings_setup));

    if (fscanf(fp, "%zu", &sprst->nsprings) != 1) {
        ERR_MSG("First line is pairsprings setup must be the number of springs");
        free(sprst);
        fclose(fp);
        return NULL;
    }

    sprst->springs = calloc(sprst->nsprings, sizeof(struct pairspring));
    struct pairspring *sprs = sprst->springs;

    size_t id = 0;
    char name1[8], name2[8];
    size_t aid1, aid2;
    while (id < sprst->nsprings) {
        size_t c = fscanf(fp,
                "%lf %lf %lf %zu %s %zu %s",
                &sprs[id].lnspr,
                &sprs[id].erspr,
                &sprs[id].fkspr,
                &aid1,
                name1,
                &aid2,
                name2);

        if (c != 7) {
            ERR_MSG("Each line in pairsprings setup must have 7 tokens");
            free(sprst->springs);
            free(sprst);
            fclose(fp);
            return NULL;
        }

        sprs[id].laspr[0] = aid1 - 1;
        sprs[id].laspr[1] = aid2 - 1;

        id++;
    }

    fclose(fp);
    return sprst;
}


static struct pairsprings_setup *_pairsprings_setup_read_json(const json_t *root) {
    if (!json_is_array(root)) {
        ERR_MSG("Pairsprings json setup must be an array");
        return NULL;
    }

    size_t nsprings = json_array_size(root);
    struct pairspring* spring_set = calloc(nsprings, sizeof(struct pairspring));

    bool error = false;
    size_t counter;
    json_t* spring;

    json_array_foreach(root, counter, spring) {
        json_error_t j_error;

        int result = json_unpack_ex(
                spring, &j_error, 0,
                "{s:F, s:F, s:F, s:i, s:i}",
                "length", &spring_set[counter].lnspr,
                "error", &spring_set[counter].erspr,
                "weight", &spring_set[counter].fkspr,
                "atom1", &spring_set[counter].laspr[0],
                "atom2", &spring_set[counter].laspr[1]);

        if (result != 0) {
            JSON_ERR_MSG(j_error, "Wrong pairspring setup json format");
            error = true;
            break;
        }

        /*double dv;
        size_t iv;

        json_t* value = json_object_get(spring, "length");
        if (!value || !(dv = json_number_value(value))) {
            ERR_MSG("sf");
            error = true;
            break;
        }
        spring_set[counter].lnspr = dv;

        value = json_object_get(spring, "error");
        if (!value || !(dv = json_number_value(value))) {
            ERR_MSG("sf");
            error = true;
            break;
        }
        spring_set[counter].erspr= dv;

        value = json_object_get(spring, "weight");
        if (!value || !(dv = json_number_value(value))) {
            ERR_MSG("sf");
            error = true;
            break;
        }
        spring_set[counter].fkspr = dv;

        value = json_object_get(spring, "atom1");
        if (!value || !(iv = json_integer_value(value))) {
            ERR_MSG("sf");
            error = true;
            break;
        }
        spring_set[counter].laspr[0] = iv - 1;

        value = json_object_get(spring, "atom2");
        if (!value || !(iv = json_integer_value(value))) {
            ERR_MSG("sf");
            error = true;
            break;
        }
        spring_set[counter].laspr[1] = iv - 1;*/
    }

    if (error) {
        free(spring_set);
        return NULL;
    }

    struct pairsprings_setup *sprst;
    sprst = calloc(1, sizeof(struct pairsprings_setup));
    sprst->springs = spring_set;
    sprst->nsprings = nsprings;

    return sprst;
}


static void _pointsprings_setup_free(struct pointsprings_setup **sprst) {
    if (*sprst != NULL) {
        for (size_t i = 0; i < (*sprst)->nsprings; i++) {
            if ((*sprst)->springs[i].laspr != NULL) {
                free((*sprst)->springs[i].laspr);
            }
        }
        free((*sprst)->springs);
        free(*sprst);
        *sprst = NULL;
    }
}


static struct pointsprings_setup *_pointsprings_setup_read_txt(const char *path) {
    FILE *fp;
    FOPEN_ELSE(fp, path, "r") {
        return NULL;
    }

    struct pointsprings_setup *sprst;
    sprst = calloc(1, sizeof(struct pointsprings_setup));

    if (fscanf(fp, "%zu", &sprst->nsprings) != 1) {
        ERR_MSG("Wrong pointsprings setup json format");
        free(sprst);
        fclose(fp);
        return NULL;
    }

    sprst->springs = calloc(sprst->nsprings, sizeof(struct pointspring));
    struct pointspring *sprs = sprst->springs;

    char name[8];
    size_t aid, naspr;
    double fkspr, X0, Y0, Z0;

    for (size_t id = 0; id < sprst->nsprings; id++) {
        size_t c = fscanf(fp, "%zu %lf %lf %lf %lf", &naspr, &fkspr, &X0, &Y0, &Z0);
        if (c != 5) {
            ERR_MSG("First line for each pointspring must have 5 tokens");
            fclose(fp);
            _pointsprings_setup_free(&sprst);
            return NULL;
        }

        sprs[id].naspr = naspr;
        sprs[id].fkspr = fkspr;
        sprs[id].X0 = X0;
        sprs[id].Y0 = Y0;
        sprs[id].Z0 = Z0;
        sprs[id].laspr = calloc(naspr, sizeof(size_t));

        for (size_t i = 0; i < naspr; i++) {
            c = fscanf(fp, "%zu %s", &aid, name);

            if (c != 2) {
                ERR_MSG("Wrong pointsprings setup format");
                fclose(fp);
                _pointsprings_setup_free(&sprst);
                return NULL;
            }

            sprs[id].laspr[i] = aid - 1;
        }
    }

    fclose(fp);
    return sprst;
}


static struct pointsprings_setup *_pointsprings_setup_read_json(const json_t *root)
{
    if (!json_is_array(root)) {
        ERR_MSG("Pointsprings json setup must be an array");
        return NULL;
    }

    size_t nsprings = json_array_size(root);
    struct pointspring *sprs = calloc(nsprings, sizeof(struct pointspring));

    size_t spring_counter;
    json_t* spring;
    json_array_foreach(root, spring_counter, spring) {
        struct pointspring* cur_spring = &sprs[spring_counter];

        json_error_t j_error;
        int code = json_unpack_ex(
                spring, &j_error, 0,
                "{s:F, s:[F,F,F]}",
                "weight", &cur_spring->fkspr,
                "coords", &cur_spring->X0, &cur_spring->Y0,  &cur_spring->Z0);
        if (code != 0) {
            JSON_ERR_MSG(j_error, "Couldn't parse json pointspring");
            free(sprs);
            return NULL;
        }

        json_t* atoms = json_object_get(spring, "atoms");
        if (!atoms || !json_is_array(atoms)) {
            ERR_MSG("Can't read pointspring atoms from json");
            free(sprs);
            return NULL;
        }

        cur_spring->naspr = json_array_size(atoms);
        cur_spring->laspr = calloc(cur_spring->naspr, sizeof(size_t));
        size_t atom_counter;
        json_t* atom_id;

        json_array_foreach(atoms, atom_counter, atom_id) {
            if (!json_is_integer(atom_id)) {
                for (size_t i = 0; i < nsprings; i++) {
                    free(sprs[i].laspr);
                }
                free(sprs);
                return NULL;
            }
            cur_spring->laspr[atom_counter] = json_integer_value(atom_id);
        }
    }

    struct pointsprings_setup* output = calloc(1, sizeof(struct pointsprings_setup));
    output->nsprings = nsprings;
    output->springs = sprs;
    return output;
}


/*struct density_setup *density_setup_create(char *fitting_pdblist, double weight, double radius) {
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


void density_setup_free(struct density_setup **prms) {
    if (prms != NULL) {
        for (int i = 0; i < prms->ag_count; i++) {
            mol_atom_group_free(prms->ag_list[i]);
        }
        free(prms->ag_list);
        free(prms);
    }
}*/


static void _noe_setup_free(struct noe_setup **nmr) {
    if (*nmr != NULL) {
        mol_noe_free((*nmr)->spec);
        free(*nmr);
        *nmr = NULL;
    }
}


#define READ_WORD(f, word, line) do { \
    (word) = NULL; \
    while (fgets(line, 512, f) != NULL) { \
        if ((line)[0] != '#') { \
            (word) = strtok((line), " \t\n"); \
            if ((word) == NULL || (word)[0] == '#') { \
                ERR_MSG("Line cannot be empty"); \
            } \
            break; \
        } \
    } \
} while(0)


static struct noe_setup *_noe_setup_read_txt(const char *sfile) {
    char line[512];
    FILE *f;
    if (!(f = fopen(sfile, "r"))) {
        ERR_MSG("jhkjh");
        return NULL;
    }

    struct noe_setup *nmr = calloc(1, sizeof(struct noe_setup));

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

    struct mol_noe_group_list *groups = mol_noe_group_list_read_txt(groups_path);
    nmr->spec = mol_noe_create(groups, freq, corr_time, mix_time, dist_cutoff);
    mol_noe_alloc_grad(nmr->spec);

    int *mask = NULL;
    if (mask_on) {
        mask = calloc(groups->ngroups * groups->ngroups, sizeof(int));
    }
    nmr->spec->exp = mol_noe_matrix_read_txt_stacked(matrix_path, groups->ngroups, mask);
    nmr->weight = weight;

    return nmr;
}


static struct noe_setup *_noe_setup_read_json(json_t *root) {
    struct mol_noe* noe = mol_noe_from_json_object(root);
    if (!noe) {
        ERR_MSG("Couldn't parse NOE setup from json");
        return NULL;
    }

    struct noe_setup* result = calloc(1, sizeof(struct noe_setup));
    result->spec = noe;

    mol_noe_alloc_grad(noe);

    json_error_t j_error;
    int code = json_unpack_ex(root, &j_error, 0, "{s:F, s:F}", "weight", &result->weight, "power", &result->power);
    if (code != 0) {
        JSON_ERR_MSG(j_error, "Couldn't unpack 'weight' and 'power' in json NOE setup");
        _noe_setup_free(&result);
        return NULL;
    }

    return result;
}


static struct noe_setup *_noe_setup_read_json_file(const char *file) {
    json_t* root = read_json_file(file);
    if (!root) {
        return NULL;
    }

    struct noe_setup* result = _noe_setup_read_json(root);
    if (!result) {
        json_decref(root);
        return NULL;
    }

    return result;
}


void energy_prms_free(struct energy_prms **prms, size_t nstages)
{
    if (*prms != NULL) {
        for (size_t i = 0; i < nstages; i++) {
            _pairsprings_setup_free(&((*prms)->sprst_pairs));
            _pointsprings_setup_free(&((*prms)->sprst_points));
            //density_setup_free(&((*prms)->fit_prms));
            _noe_setup_free(&((*prms)->nmr));
            _fixed_setup_free(&((*prms)->fixed));
        }
        free(*prms);
        *prms = NULL;
    }
}


bool energy_prms_populate_from_options(
        struct energy_prms **result_energy_prm,
        size_t *result_nstages,
        const struct options opts)
{
    struct energy_prms* all_stage_prms = calloc(1, sizeof(struct energy_prms));
    size_t nstages = 1;

    *result_energy_prm = NULL;
    *result_nstages = 0;

    all_stage_prms->ag = NULL;
    all_stage_prms->ag_setup = NULL;
    all_stage_prms->ace_setup = NULL;
    all_stage_prms->json_log = NULL;

    all_stage_prms->sprst_pairs = NULL;
    all_stage_prms->sprst_points = NULL;
    all_stage_prms->nmr = NULL;
    all_stage_prms->fit_prms = NULL;
    all_stage_prms->fixed = NULL;

    all_stage_prms->nsteps = opts.nsteps;

    all_stage_prms->bonds = opts.bonds;
    all_stage_prms->angles = opts.angles;
    all_stage_prms->dihedrals = opts.dihedrals;
    all_stage_prms->impropers = opts.impropers;
    all_stage_prms->vdw = opts.vdw;
    all_stage_prms->vdw03 = opts.vdw03;
    all_stage_prms->gbsa = opts.gbsa;

    all_stage_prms->score_only = opts.score_only;
    all_stage_prms->verbose = opts.verbosity;

    all_stage_prms->json_log_setup.print_step = opts.print_step;
    all_stage_prms->json_log_setup.print_stage = opts.print_stage;
    all_stage_prms->json_log_setup.print_noe_matrix = opts.print_noe_matrix;

    // Setup global fixed atoms
    if ((int)(opts.fixed_pdb != NULL) + (int)(opts.fix_receptor) + (int)(opts.fix_ligand) > 1) {
        ERR_MSG("Different fixed atoms flags can't be combined");
        energy_prms_free(&all_stage_prms, nstages);
        return false;
    }

    if (opts.rec_natoms == 0 && opts.lig_natoms == 0 && opts.separate) {
        ERR_MSG("rec_natoms and lig_natoms can't be zero in rec/lig mode, populate them first");
        energy_prms_free(&all_stage_prms, nstages);
        return false;
    }

    if (opts.fixed_pdb) {
        all_stage_prms->fixed = _fixed_setup_read_txt(opts.fixed_pdb);
        if (!all_stage_prms->fixed) {
            ERR_MSG("Couldn't parse fixed atoms from %s", opts.fixed_pdb);
            energy_prms_free(&all_stage_prms, nstages);
            return false;
        }
    } else if (opts.fix_receptor) {
        if (all_stage_prms->fixed) {
            ERR_MSG("Cannot use fix-receptor flag with another one");
            energy_prms_free(&all_stage_prms, nstages);
            return false;
        }
        all_stage_prms->fixed = _fixed_setup_atom_range(0, opts.rec_natoms);
    } else if (opts.fix_ligand) {
        if (all_stage_prms->fixed) {
            ERR_MSG("Cannot use fix-ligand flag with another one");
            energy_prms_free(&all_stage_prms, nstages);
            return false;
        }
        all_stage_prms->fixed = _fixed_setup_atom_range(opts.rec_natoms, opts.rec_natoms + opts.lig_natoms);
    }

    if (opts.pair_springs_txt) {
        all_stage_prms->sprst_pairs = _pairsprings_setup_read_txt(opts.pair_springs_txt);

        if (!all_stage_prms->sprst_pairs) {
            ERR_MSG("Couldn't parse pairsprings from %s", opts.pair_springs_txt);
            energy_prms_free(&all_stage_prms, nstages);
            return false;
        }
    }

    if (opts.point_springs_txt) {
        all_stage_prms->sprst_points = _pointsprings_setup_read_txt(opts.point_springs_txt);

        if (!all_stage_prms->sprst_points) {
            ERR_MSG("Couldn't parse pointsprings from %s", opts.point_springs_txt);
            energy_prms_free(&all_stage_prms, nstages);
            return false;
        }
    }

    if (opts.noe_json) {
        all_stage_prms->nmr = _noe_setup_read_json_file(opts.noe_json);

        if (!all_stage_prms->nmr) {
            ERR_MSG("Couldn't parse NOE from %s", opts.noe_json);
            energy_prms_free(&all_stage_prms, nstages);
            return false;
        }
    }

    if (opts.noe_txt) {
        if (all_stage_prms->nmr) {
            ERR_MSG("Cannot use NOE txt and json formats at the same time");
            energy_prms_free(&all_stage_prms, nstages);
            return false;
        }

        all_stage_prms->nmr = _noe_setup_read_txt(opts.noe_txt);

        if (!all_stage_prms->nmr) {
            ERR_MSG("Couldn't parse NOE from %s", opts.noe_txt);
            energy_prms_free(&all_stage_prms, nstages);
            return false;
        }
    }

    /*
    if (opts.density_json) {
        all_stage_prms->fit_prms = density_setup_read_json_file(opts.density_json);

        if (!all_stage_prms->fit_prms) {
            ERR_MSG("SDF");
            energy_prms_free(&all_stage_prms, nstages);
            return false;
        }
    }*/

    // read json
    json_t* setup = NULL;

    if (opts.setup_json) {
        json_t *setup_root = read_json_file(opts.setup_json);
        if (!setup_root) {
            energy_prms_free(&all_stage_prms, nstages);
            return false;
        }

        setup = json_object_get(setup_root, "stages");
        if (setup && !json_is_array(setup)) {
            ERR_MSG("Key 'stages' must point to a dictionary");
            json_decref(setup);
            energy_prms_free(&all_stage_prms, nstages);
            return false;
        }
    }

    if (setup) {
        nstages = json_array_size(setup);

        // copy default params to each stage
        struct energy_prms* buffer = calloc(nstages, sizeof(struct energy_prms));
        for (size_t i = 0; i < nstages; i++) {
            // TODO: need to recursive copy of setups, because later on the
            //       pointers can be overwritten and never be freed
            memcpy(buffer + i, all_stage_prms, sizeof(struct energy_prms));
        }
        free(all_stage_prms);
        all_stage_prms = buffer;

        json_t* stage_desc;
        size_t stage_id;
        bool error = false;

        json_array_foreach(setup, stage_id, stage_desc) {
            struct energy_prms* stage_prms = &all_stage_prms[stage_id];

            bool stage_fix_rec = false;
            bool stage_fix_lig = false;

            json_error_t j_error;
            int code = json_unpack_ex(
                    stage_desc, &j_error, 0,
                    "{s?b, s?b, s?b, s?b, s?b, s?b, s?b, s?b, s?b, s?b, s?b, s:i}",
                    "bonds", &stage_prms->bonds,
                    "angles", &stage_prms->angles,
                    "dihedrals", &stage_prms->dihedrals,
                    "impropers", &stage_prms->impropers,
                    "vdw", &stage_prms->vdw,
                    "vdw03", &stage_prms->vdw03,
                    "gbsa", &stage_prms->gbsa,
                    "fix_receptor", &stage_fix_rec,
                    "fix_ligand", &stage_fix_lig,
                    "score_only", &stage_prms->score_only,
                    "verbosity", &stage_prms->verbose,
                    "nsteps", &stage_prms->nsteps);

            if (code != 0) {
                JSON_ERR_MSG(j_error, "Couldn't parse stage flags in json");
                error = true;
                break;
            }

            if ((stage_fix_rec || stage_fix_lig) && !opts.separate) {
                ERR_MSG("You can't provide fix-? flags in a single file mode");
                error = true;
                break;
            }

            json_t* stage_fixed = json_object_get(stage_desc, "fixed");
            if ((int)(stage_fixed != NULL) + (int)(stage_fix_rec) + (int)(stage_fix_lig) > 1) {
                ERR_MSG("Different fixed atoms flags can't be combined");
                error = true;
                break;
            }

            // Read fixed atoms
            if (stage_fixed) {
                _fixed_setup_free(&stage_prms->fixed);
                stage_prms->fixed = _fixed_setup_read_json(stage_fixed);
                if (!stage_prms->fixed) {
                    ERR_MSG("Couldn't parse fixed atoms from %s", opts.setup_json);
                    error = true;
                    break;
                }
            } else if (stage_fix_rec) {
                _fixed_setup_free(&stage_prms->fixed);
                stage_prms->fixed = _fixed_setup_atom_range(0, opts.rec_natoms);
            } else if (stage_fix_lig) {
                _fixed_setup_free(&stage_prms->fixed);
                stage_prms->fixed = _fixed_setup_atom_range(opts.rec_natoms, opts.rec_natoms + opts.lig_natoms);
            }

            // Pairsprings
            json_t* stage_pairsprings = json_object_get(stage_desc, "pairsprings");
            if (stage_pairsprings) {
                stage_prms->sprst_pairs = _pairsprings_setup_read_json(stage_pairsprings);
                if (!stage_prms->sprst_pairs) {
                    ERR_MSG("Couldn't parse pairsprings from %s", opts.setup_json);
                    error = true;
                    break;
                }
            }

            // Pointsprings
            json_t* stage_pointsprings = json_object_get(stage_desc, "pointsprings");
            if (stage_pointsprings) {
                stage_prms->sprst_points = _pointsprings_setup_read_json(stage_pointsprings);
                if (!stage_prms->sprst_points) {
                    ERR_MSG("Couldn't parse poinsprings from %s", opts.setup_json);
                    error = true;
                    break;
                }
            }

            /*json_t* stage_density = json_object_get(stage_desc, "density");
            if (stage_density) {
                stage_prms->fit_prms = density_setup_read_json(stage_density);
                json_decref(stage_density);
                if (!stage_prms->fit_prms) {
                    error = true;
                    break;
                }
            }*/

            // NOE
            json_t* stage_noe = json_object_get(stage_desc, "noe");
            if (stage_noe) {
                stage_prms->nmr = _noe_setup_read_json(stage_noe);
                if (!stage_prms->nmr) {
                    ERR_MSG("Couldn't parse NOE from %s", opts.setup_json);
                    error = true;
                    break;
                }
            }
        }

        json_decref(setup);

        if (error) {
            energy_prms_free(&all_stage_prms, nstages);
            return false;
        }
    }

    *result_energy_prm = all_stage_prms;
    *result_nstages = nstages;

    return true;
}
