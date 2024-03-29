#include "utils.h"
#include "setup.h"
#include "parse_options.h"

#include "mol2/json.h"
#include "mol2/icharmm.h"


/**
 * defines function _{POTENTIAL_NAME}_setup_read_json_file(const char *path)
 */
#define _POTENTIAL_SETUP_READ_JSON_FILE(NAME)                                                   \
    static struct NAME ## _setup* _ ## NAME ## _setup_read_json_file(const char *path) {        \
        json_t* root = read_json_file(path);                                                    \
        if (!root) {                                                                            \
            return NULL;                                                                        \
        }                                                                                       \
        struct NAME ## _setup* result = _ ## NAME ## _setup_read_json(root);                    \
        if (!result) {                                                                          \
            json_decref(root);                                                                  \
            return NULL;                                                                        \
        }                                                                                       \
        json_decref(root);                                                                      \
        return result;                                                                          \
    }


/* **************************************
 *          Creating atom groups        *
 ****************************************/


static struct mol_atom_group_list *_mol_atom_group_list_from_pdb(const char* pdb)
{
    DEBUG_MSG("Reading %s", pdb);
    DEBUG_MSG("Trying to read %s as a multimodel one (with MODEL records)", pdb);
    struct mol_atom_group_list *ag_list = mol_read_pdb_models(pdb);

    // If not multimodel, try to read as a single model
    if (ag_list == NULL) {
        DEBUG_MSG("File %s doesn't have MODEL records, reading as single model", pdb);
        struct mol_atom_group* ag_single = mol_read_pdb(pdb);

        // If both single and multi-model failed - give up
        if (ag_single == NULL) {
            ERR_MSG("Failed parsing %s", pdb);
            return NULL;
        }

        ag_list = mol_atom_group_list_create(1);
        ag_list->members[0] = *ag_single;
        free(ag_single);
    }

    DEBUG_MSG("Succesfully parsed %s (%zu models)", pdb, ag_list->size);
    return ag_list;
}


static struct mol_atom_group_list *_read_ag_list(
        const char *prm,
        const char *rtf,
        const char *pdb,
        const char *psf,
        const char *json,
        const int score_only) {

    DEBUG_MSG("Started reading a new atom group");
    struct mol_atom_group_list* ag_list = NULL;
    struct mol_atom_group* ag_json = NULL;

    const char* mol_name;

    // Read molecule from json
    if (json != NULL) {
        DEBUG_MSG("Reading json file %s", json);
        mol_name = json;
        ag_json = mol_read_json(json);

        if (ag_json == NULL) {
            ERR_MSG("Failed reading %s", json);
            return NULL;
        }
    }

    // Read molecule from pdb
    if (pdb != NULL) {
        mol_name = pdb;
        ag_list = _mol_atom_group_list_from_pdb(pdb);
        if (!ag_list) {
            mol_atom_group_free(ag_json);
            return NULL;
        }

        // If json was provied too then merge geometry from json and coords from pdb
        if (ag_json != NULL) {
            DEBUG_MSG("Reading coordinates from %s and geometry from %s", pdb, json);

            struct mol_atom_group_list *final_aglist = mol_atom_group_list_create(ag_list->size);

            for (size_t i = 0; i < final_aglist->size; i++) {
                struct mol_atom_group* tmp_copy = mol_atom_group_copy(ag_json);

                // TODO: These are initialized from the beginning for some reason. Create issue in libmol
                //free_if_not_null(tmp_copy->active_atoms);
                //free_if_not_null(tmp_copy->active_bonds);
                //free_if_not_null(tmp_copy->active_angles);
                //free_if_not_null(tmp_copy->active_dihedrals);
                //free_if_not_null(tmp_copy->active_impropers);

                final_aglist->members[i] = *tmp_copy;
                free(tmp_copy);

                if (final_aglist->members[i].natoms != ag_list->members[i].natoms) {
                    ERR_MSG("Model %zu in %s and %s have different atom numbers (%zu, %zu)",
                            i, pdb, json,
                            ag_list->members[i].natoms,
                            final_aglist->members[i].natoms);
                    mol_atom_group_list_free(final_aglist);
                    mol_atom_group_list_free(ag_list);
                    mol_atom_group_free(ag_json);
                    return NULL;
                }
                memcpy(final_aglist->members[i].coords,
                        ag_list->members[i].coords,
                        sizeof(struct mol_vector3) * ag_list->members[i].natoms);
            }

            mol_atom_group_list_free(ag_list);
            mol_atom_group_free(ag_json);
            ag_list = final_aglist;

        } else if (psf && prm && rtf) {
            // If json wasn't provided read geometry from prm rtf psf
            DEBUG_MSG("Reading geometry from %s", psf);
            struct mol_atom_group_list* ag_list_copy =  mol_atom_group_list_create(ag_list->size);

            // Fill the geometry for the first model
            if (!mol_atom_group_read_geometry(&ag_list->members[0], psf, prm, rtf)) {
                ERR_MSG("Couldn't fill geometry from psf rtf and prm");
                mol_atom_group_list_free(ag_list);
                mol_atom_group_list_free(ag_list_copy);
                return NULL;
            }

            // Copy the first model into each member of ag_list_copy and fill the coordinates from ag_list
            // this is done for speed up, because reading geomerty for each model separately can be very slow
            for (size_t i = 0; i < ag_list->size; i++) {
                struct mol_atom_group* buf = mol_atom_group_copy(&ag_list->members[0]);
                ag_list_copy->members[i] = *buf;
                memcpy(ag_list_copy->members[i].coords,
                        ag_list->members[i].coords,
                        ag_list->members[i].natoms * sizeof(struct mol_vector3));
                free(buf);
            }
            mol_atom_group_list_free(ag_list);
            ag_list = ag_list_copy;
            DEBUG_MSG("Done reading geometry");

        } else if (score_only != 0) {
            // if they weren't provided, check the score_only flag
            WRN_MSG("Geometry for %s is not provided, computing only non-parametric terms", mol_name);

        } else {
            // Give up if nothing worked
            ERR_MSG("Force field and geometry are not provided for %s", mol_name);
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


/*
 * Copy single model in ag_list N times
 */
static bool _expand_mol_atom_group_list(struct mol_atom_group_list* ag_list, const size_t num_models)
{
    if (ag_list->size != 1) {
        ERR_MSG("Cannot expand atom group list of size != 1");
        return false;
    }
    struct mol_atom_group* model = mol_atom_group_copy(ag_list->members);
    mol_atom_group_free(ag_list->members);

    ag_list->size = num_models;
    ag_list->members = calloc(num_models, sizeof(struct mol_atom_group));
    for (size_t i = 0; i < num_models; i++) {
        struct mol_atom_group* buf = mol_atom_group_copy(model);
        ag_list->members[i] = *buf;
        free(buf);
    }
    mol_atom_group_free(model);
    return true;
}


static struct mol_atom_group_list* _merge_ag_lists(
        struct mol_atom_group_list* ag1,
        struct mol_atom_group_list* ag2) {
    if (ag1->size != ag2->size) {
        size_t num_models = MAX(ag1->size, ag2->size);
        if (ag1->size == 1) {
            DEBUG_MSG("Using the same receptor model with %zu ligand models", num_models);
            _expand_mol_atom_group_list(ag1, num_models);
        } else if (ag2->size == 1) {
            DEBUG_MSG("Using the same ligand model with %zu receptor models", num_models);
            _expand_mol_atom_group_list(ag2, num_models);
        } else {
            ERR_MSG("Receptor and ligand have different numbers of models (%zu, %zu)", ag1->size, ag2->size);
            return NULL;
        }
    }

    DEBUG_MSG("Using models assembled from receptor and ligand models provided separately");
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

        mol_atom_group_list_free(rec_list);
        mol_atom_group_list_free(lig_list);
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

    opts->ag_list = ag_list;
    opts->num_models = 0;
    if (ag_list) {
        opts->num_models = ag_list->size;
    }
    return ag_list;
}


/* **************************************
 *               Fixed atoms            *
 ****************************************/


static void _fixed_setup_free(struct fixed_setup **fixed)
{
    if (*fixed != NULL) {
        (*fixed)->ref_count -= 1;
        if ((*fixed)->ref_count == 0) {
            if ((*fixed)->atoms) {
                free((*fixed)->atoms);
                (*fixed)->atoms = NULL;
            }
            free(*fixed);
        }
        *fixed = NULL;
    }
}


static struct fixed_setup_multi* _fixed_setup_multi_create(size_t size)
{
    struct fixed_setup_multi* out = calloc(1, sizeof(struct fixed_setup_multi));
    out->size = size;
    out->setups = calloc(size, sizeof(struct fixed_setup*));
    out->ref_count = 1;
    return out;
}


static void _fixed_setup_multi_free(struct fixed_setup_multi **fixed)
{
    if (*fixed != NULL) {
        (*fixed)->ref_count -= 1;
        if ((*fixed)->ref_count == 0) {
            if ((*fixed)->setups != NULL) {
                for (size_t i = 0; i < (*fixed)->size; i++) {
                    _fixed_setup_free((*fixed)->setups + i);
                }
                free((*fixed)->setups);
                (*fixed)->setups = NULL;
            }
            free(*fixed);
            *fixed = NULL;
        }
    }
}


static struct fixed_setup* _fixed_atoms_from_ag(struct mol_atom_group* ag, struct mol_atom_group* fix_ag, double cutoff) {
    bool *fix_mask = calloc(ag->natoms, sizeof(bool));
    for (size_t i = 0; i < ag->natoms; i++) {
        for (size_t j = 0; j < fix_ag->natoms; j++) {
            if (MOL_VEC_EUCLIDEAN_DIST(ag->coords[i], fix_ag->coords[j]) < cutoff) {
                fix_mask[i] = true;
            }
        }
    }

    struct fixed_setup* out = calloc(1, sizeof(struct fixed_setup));
    out->ref_count = 1;
    out->natoms = 0;
    for (size_t i = 0; i < ag->natoms; i++) {
        if (fix_mask[i]) {
            out->natoms += 1;
        }
    }

    out->atoms = calloc(out->natoms, sizeof(out->natoms));
    size_t counter = 0;
    for (size_t i = 0; i < ag->natoms; i++) {
        if (fix_mask[i]) {
            out->atoms[counter++] = i;
        }
    }

    free(fix_mask);
    return out;
}


static struct fixed_setup_multi* _fixed_setup_multi_read_from_pdb(const char* fix_pdb, const struct mol_atom_group_list* ag) {
    if (!ag) {
        ERR_MSG("Minimized atom group is empty");
        return NULL;
    }

    struct mol_atom_group_list* fix_ag = _mol_atom_group_list_from_pdb(fix_pdb);
    if (!fix_ag) {
        ERR_MSG("Could not parse %s", fix_pdb);
        return NULL;
    }

    if ((fix_ag->size != 1) && (fix_ag->size != ag->size)) {
        mol_atom_group_list_free(fix_ag);
        ERR_MSG("Number of models in fixed pdb must equal 1 or match the input pdb");
        return NULL;
    }

    if ((fix_ag->size == 1) && (ag->size > 1)) {
        DEBUG_MSG("Found one model in fixed pdb and multiple models in minimized pdb. "
            "Applying fixed pdb to all models");
    }
    struct fixed_setup_multi* fix = _fixed_setup_multi_create(ag->size);
    for (size_t i = 0; i < ag->size; i++) {
        struct mol_atom_group* ref_model = ag->members + i;
        struct mol_atom_group* fix_model = fix_ag->members;
        if (fix_ag->size > 1) {
            fix_model += i;
        }
        fix->setups[i] = _fixed_atoms_from_ag(ref_model, fix_model, 0.01);
        if (fix->setups[i]->natoms == 0) {
            WRN_MSG("No atoms were fixed for model %zu", i);
        }
    }
    mol_atom_group_list_free(fix_ag);
    return fix;
}


static struct fixed_setup_multi* _fixed_setup_multi_read_json(const json_t *root, const size_t num_models) {
    if (!json_is_array(root)) {
        ERR_MSG("Fixed atoms json setup must be an array");
        return NULL;
    }

    struct fixed_setup *result = calloc(1, sizeof(struct fixed_setup));
    result->ref_count = 1;
    result->natoms = json_array_size(root);
    result->atoms = calloc(result->natoms, sizeof(size_t));

    size_t counter;
    json_t *atom_id;
    json_array_foreach(root, counter, atom_id) {
        if (!json_is_integer(atom_id)) {
            _fixed_setup_free(&result);
            json_decref(atom_id);
            return NULL;
        }
        result->atoms[counter] = json_integer_value(atom_id);
    }

    struct fixed_setup_multi* out = _fixed_setup_multi_create(num_models);
    for (size_t i = 0; i < num_models; i++) {
        out->setups[i] = result;
        result->ref_count += 1;
    }
    result->ref_count -= 1;
    return out;
}


static struct fixed_setup* _fixed_setup_atom_range(const size_t start_atom, const size_t end_atom)
{
    struct fixed_setup* fixed = calloc(1, sizeof(struct fixed_setup));
    fixed->ref_count = 1;
    fixed->natoms = end_atom - start_atom;
    fixed->atoms = calloc(fixed->natoms, sizeof(size_t));
    for (size_t i = 0; i < fixed->natoms; i++) {
        fixed->atoms[i] = start_atom + i;
    }
    return fixed;
}


static struct fixed_setup_multi* _fixed_setup_multi_atom_range(size_t size, const size_t start_atom, const size_t end_atom)
{
    struct fixed_setup_multi* out = _fixed_setup_multi_create(size);
    for (size_t i = 0; i < size; i++) {
        struct fixed_setup* buf = _fixed_setup_atom_range(start_atom, end_atom);
        out->setups[i] = buf;
    }
    return out;
}


/* **************************************
 *               Pairsprings            *
 ****************************************/


static void _pairsprings_setup_free(struct pairsprings_setup **sprst) {
    if (*sprst != NULL) {
        size_t i;
        for (i = 0;  i < (*sprst)->nsprings; i++) {
            if ((*sprst)->springs[i].group1 != NULL) {
                free((*sprst)->springs[i].group1);
                free((*sprst)->springs[i].group2);
                (*sprst)->springs[i].group1 = NULL;
                (*sprst)->springs[i].group2 = NULL;
            }
        }
        free((*sprst)->springs);
        (*sprst)->springs = NULL;
        free(*sprst);
        *sprst = NULL;
    }
}


static struct pairsprings_setup *_pairsprings_setup_read_txt(const char *path) {
    WRN_MSG("You are reading pairsprings in txt format. This is a legacy "
            "format and will be removed in the future versions.");

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
                          &sprs[id].distance,
                          &sprs[id].lerror,
                          &sprs[id].weight,
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
        sprs[id].rerror = sprs[id].lerror;
        sprs[id].group_size1 = 1;
        sprs[id].group_size2 = 1;
        sprs[id].group1 = calloc(1, sizeof(size_t));
        sprs[id].group2 = calloc(1, sizeof(size_t));
        sprs[id].group1[0] = aid1 - 1;
        sprs[id].group2[0] = aid2 - 1;
        sprs[id].average = 0;
        sprs[id].potential = 2;

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
    json_t* g1;
    json_t* g2;

    size_t len1, len2, i;
    json_array_foreach(root, counter, spring) {
        json_error_t j_error;
        int result = json_unpack_ex(
                spring, &j_error, 0,
                "{s:F, s:F, s:F, s:F, s:i, s:i}",
                "distance", &spring_set[counter].distance,
                "lerror", &spring_set[counter].lerror,
                "rerror", &spring_set[counter].rerror,
                "weight", &spring_set[counter].weight,
                "potential", &spring_set[counter].potential,
                "average", &spring_set[counter].average);

        if (result != 0) {
            JSON_ERR_MSG(j_error, "Wrong pairspring setup json format");
            error = true;
            break;
        }
        if (spring_set[counter].average >1 || spring_set[counter].average <0) {
            ERR_MSG("The average should be 0 (SUM), or 1 (R-6).");
            error = true;
            break;
        }

        if (spring_set[counter].potential>2 || spring_set[counter].potential <0) {
            ERR_MSG("Potential in pairsprings should be one of 0 (SQUARE-WEll), 1 (BIHARMONIC), 2 (SOFT-SQUARE).");
            error = true;
            break;
        }
        g1 = json_object_get(spring, "group1");
        g2 = json_object_get(spring, "group2");
        len1 = json_array_size(g1);
        len2 = json_array_size(g2);
        if (len1==0 || len2==0) {
            ERR_MSG("Wrong pairspring setup group size.");
            error = true;
            break;
        }
        spring_set[counter].group_size1 = len1;
        spring_set[counter].group_size2 = len2;
        spring_set[counter].group1 = calloc(len1, sizeof(size_t));
        spring_set[counter].group2 = calloc(len2, sizeof(size_t));

        for (i = 0; i < len1; ++i) {
            if (json_is_integer(json_array_get(g1, i))) {
                spring_set[counter].group1[i] = json_integer_value(json_array_get(g1, i));
            } else {
                ERR_MSG("Index in group1 must be interger");
            }
        }
        for (i = 0; i < len2; ++i) {
            if (json_is_integer(json_array_get(g2, i))) {
                spring_set[counter].group2[i] = json_integer_value(json_array_get(g2, i));
            } else {
                ERR_MSG("Index in group2 must be interger");
            }
        }
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


//_POTENTIAL_SETUP_READ_JSON_FILE(pairsprings)


/* **************************************
 *            Pointsprings              *
 ****************************************/


static void _pointsprings_setup_free(struct pointsprings_setup **sprst) {
    if (*sprst != NULL) {
        for (size_t i = 0; i < (*sprst)->nsprings; i++) {
            if ((*sprst)->springs[i].atoms != NULL) {
                free((*sprst)->springs[i].atoms);
                (*sprst)->springs[i].atoms = NULL;
            }
        }
        free((*sprst)->springs);
        (*sprst)->springs = NULL;
        free(*sprst);
        *sprst = NULL;
    }
}


static struct pointsprings_setup *_pointsprings_setup_read_txt(const char *path) {
    WRN_MSG("You are reading pointsprings in txt format. This is a legacy "
            "format and will be removed in the future versions.");

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

        sprs[id].natoms = naspr;
        sprs[id].weight = fkspr;
        sprs[id].X0 = X0;
        sprs[id].Y0 = Y0;
        sprs[id].Z0 = Z0;
        sprs[id].atoms = calloc(naspr, sizeof(size_t));

        for (size_t i = 0; i < naspr; i++) {
            c = fscanf(fp, "%zu %s", &aid, name);

            if (c != 2) {
                ERR_MSG("Wrong pointsprings setup format");
                fclose(fp);
                _pointsprings_setup_free(&sprst);
                return NULL;
            }

            sprs[id].atoms[i] = aid - 1;
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
                "weight", &cur_spring->weight,
                "coords", &cur_spring->X0, &cur_spring->Y0,  &cur_spring->Z0);
        if (code != 0) {
            JSON_ERR_MSG(j_error, "Couldn't parse json pointspring");
            free(sprs);
            return NULL;
        }

        json_t* atoms = json_object_get(spring, "atoms");
        if (!json_is_array(atoms)) {
            ERR_MSG("Can't read pointspring atoms from json");
            free(sprs);
            return NULL;
        }

        cur_spring->natoms = json_array_size(atoms);
        cur_spring->atoms = calloc(cur_spring->natoms, sizeof(size_t));
        size_t atom_counter;
        json_t* atom_id;

        json_array_foreach(atoms, atom_counter, atom_id) {
            if (!json_is_integer(atom_id)) {
                for (size_t i = 0; i < nsprings; i++) {
                    free(sprs[i].atoms);
                }
                free(sprs);
                return NULL;
            }
            cur_spring->atoms[atom_counter] = json_integer_value(atom_id);
        }
    }

    struct pointsprings_setup* output = calloc(1, sizeof(struct pointsprings_setup));
    output->nsprings = nsprings;
    output->springs = sprs;
    return output;
}


//_POTENTIAL_SETUP_READ_JSON_FILE(pointsprings)


/* **************************************
 *               Density                *
 ****************************************/


static void _density_setup_free(struct density_setup** density)
{
    if (*density != NULL) {
        if ((*density)->ag != NULL) {
            mol_atom_group_free((*density)->ag);
        }
        free(*density);
        *density = NULL;
    }
}


static struct density_setup* _density_setup_read_json(json_t* root)
{
    json_error_t j_error;
    int code;

    struct density_setup *density = calloc(1, sizeof(struct density_setup));
    density->ag = NULL;
    char* pdb;

    code = json_unpack_ex(
            root, &j_error, JSON_STRICT,
            "{s:s, s:F, s:F}",
            "pdb", &pdb,
            "weight", &density->weight,
            "atom_radius", &density->prms.radius);

    if (code != 0) {
        JSON_ERR_MSG(j_error, "Couldn't parse density setup");
        _density_setup_free(&density);
        return NULL;
    }

    struct mol_atom_group* ag = mol_read_pdb(pdb);
    if (!ag) {
        ERR_MSG("Couldn't read atom group in density setup");
        _density_setup_free(&density);
        return NULL;
    }
    density->ag = ag;

    return density;
}


_POTENTIAL_SETUP_READ_JSON_FILE(density)


/* **************************************
 *                 NOE                  *
 ****************************************/


#ifdef NOE

static void _noe_setup_free(struct noe_setup **noe) {
    if (*noe != NULL) {
        mol_noe_free((*noe)->spec);
        free(*noe);
        *noe = NULL;
    }
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


_POTENTIAL_SETUP_READ_JSON_FILE(noe)

#endif


/* **************************************
 *      Populate energy parameters      *
 ****************************************/


void energy_prms_free(struct energy_prms **prms, size_t nstages)
{
    if (*prms != NULL) {
        for (size_t i = 0; i < nstages; i++) {
            _pairsprings_setup_free(&((*prms + i)->sprst_pairs));
            _pointsprings_setup_free(&((*prms + i)->sprst_points));
            _density_setup_free(&((*prms + i)->density));
#ifdef NOE
            _noe_setup_free(&((*prms + i)->nmr));
#endif
            _fixed_setup_multi_free(&((*prms + i)->fixed));
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
    all_stage_prms->density = NULL;
    all_stage_prms->fixed = NULL;
#ifdef NOE
    all_stage_prms->nmr = NULL;
#endif

    all_stage_prms->nsteps = opts.nsteps;

    all_stage_prms->bonds = opts.bonds;
    all_stage_prms->angles = opts.angles;
    all_stage_prms->dihedrals = opts.dihedrals;
    all_stage_prms->impropers = opts.impropers;
    all_stage_prms->vdw = opts.vdw;
    all_stage_prms->vdw03 = opts.vdw03;
    all_stage_prms->eleng = opts.eleng;
    all_stage_prms->elengs03 = opts.elengs03;
    all_stage_prms->ace = opts.ace;
    all_stage_prms->gbsa = opts.gbsa;

    all_stage_prms->tol = opts.tol;
    all_stage_prms->nbcut = opts.nbcut;
    all_stage_prms->ace_efac = opts.ace_efac;
    all_stage_prms->scale_vdw_s03 = opts.scale_vdw_s03;
    all_stage_prms->scale_coul_s03 = opts.scale_coul_s03;
    all_stage_prms->eeps = opts.eeps;
    all_stage_prms->gbcut = opts.gbcut;

    all_stage_prms->score_only = opts.score_only;
    all_stage_prms->pull_ligand_away = opts.pull_ligand_away;

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
        DEBUG_MSG("Reading %s", opts.fixed_pdb);
        all_stage_prms->fixed = _fixed_setup_multi_read_from_pdb(opts.fixed_pdb, opts.ag_list);
        if (!all_stage_prms->fixed) {
            ERR_MSG("Couldn't get fixed atoms from %s", opts.fixed_pdb);
            energy_prms_free(&all_stage_prms, nstages);
            return false;
        }
    } else if (opts.fix_receptor) {
        DEBUG_MSG("Fixing receptor atoms");
        if (all_stage_prms->fixed) {
            ERR_MSG("Cannot use fix-receptor flag with another one");
            energy_prms_free(&all_stage_prms, nstages);
            return false;
        }
        all_stage_prms->fixed = _fixed_setup_multi_atom_range(opts.num_models, 0, opts.rec_natoms);
    } else if (opts.fix_ligand) {
        DEBUG_MSG("Fixing ligand atoms");
        if (all_stage_prms->fixed) {
            ERR_MSG("Cannot use fix-ligand flag with another one");
            energy_prms_free(&all_stage_prms, nstages);
            return false;
        }
        all_stage_prms->fixed = _fixed_setup_multi_atom_range(opts.num_models, opts.rec_natoms, opts.rec_natoms + opts.lig_natoms);
    }

    if (opts.pair_springs_txt) {
        DEBUG_MSG("Parsing %s", opts.pair_springs_txt);
        all_stage_prms->sprst_pairs = _pairsprings_setup_read_txt(opts.pair_springs_txt);

        if (!all_stage_prms->sprst_pairs) {
            ERR_MSG("Couldn't parse pairsprings from %s", opts.pair_springs_txt);
            energy_prms_free(&all_stage_prms, nstages);
            return false;
        }
    }

    if (opts.point_springs_txt) {
        DEBUG_MSG("Parsing %s", opts.point_springs_txt);
        all_stage_prms->sprst_points = _pointsprings_setup_read_txt(opts.point_springs_txt);

        if (!all_stage_prms->sprst_points) {
            ERR_MSG("Couldn't parse pointsprings from %s", opts.point_springs_txt);
            energy_prms_free(&all_stage_prms, nstages);
            return false;
        }
    }

    if (opts.noe_json) {
        DEBUG_MSG("Parsing %s", opts.noe_json);
#ifdef NOE
        all_stage_prms->nmr = _noe_setup_read_json_file(opts.noe_json);

        if (!all_stage_prms->nmr) {
            ERR_MSG("Couldn't parse NOE from %s", opts.noe_json);
            energy_prms_free(&all_stage_prms, nstages);
            return false;
        }
#else
        ERR_MSG("Rebuild NRGmin with -DNOE=ON to use NOE spectrum fitting");
        energy_prms_free(&all_stage_prms, nstages);
        return false;
#endif
    }

    if (opts.density_json) {
        DEBUG_MSG("Parsing %s", opts.density_json);
        all_stage_prms->density = _density_setup_read_json_file(opts.density_json);

        if (!all_stage_prms->density) {
            ERR_MSG("Couldn't parse density from %s", opts.density_json);
            energy_prms_free(&all_stage_prms, nstages);
            return false;
        }
    }

    // read json
    json_t* setup = NULL;
    json_t *setup_root = NULL;

    if (opts.setup_json) {
        DEBUG_MSG("Parsing %s", opts.setup_json);
        setup_root = read_json_file(opts.setup_json);
        if (!setup_root) {
            energy_prms_free(&all_stage_prms, nstages);
            return false;
        }

        setup = json_object_get(setup_root, "stages");
        if (!json_is_array(setup)) {
            ERR_MSG("Key 'stages' must point to a dictionary");
            json_decref(setup_root);
            energy_prms_free(&all_stage_prms, nstages);
            return false;
        }
    }

    if (setup) {
        nstages = json_array_size(setup);
        DEBUG_MSG("Found %zu stages", nstages);

        // copy global params to each stage
        struct energy_prms* buffer = calloc(nstages, sizeof(struct energy_prms));
        for (size_t i = 0; i < nstages; i++) {
            // TODO: need to recursive copy of setups, because later on the
            //       pointers can be overwritten and never be freed
            memcpy(buffer + i, all_stage_prms, sizeof(struct energy_prms));
            if (all_stage_prms->fixed != NULL) {
                all_stage_prms->fixed->ref_count += 1;
            }
        }
        if (all_stage_prms->fixed != NULL) {
            all_stage_prms->fixed->ref_count -= 1;
        }
        free(all_stage_prms);
        all_stage_prms = buffer;

        json_t* stage_desc;
        size_t stage_id;
        bool error = false;

        json_array_foreach(setup, stage_id, stage_desc) {
            DEBUG_MSG("Parsing stage %zu", stage_id);
            struct energy_prms* stage_prms = &all_stage_prms[stage_id];

            bool stage_fix_rec = false;
            bool stage_fix_lig = false;

            json_error_t j_error;
            DEBUG_MSG("Unpacking options");
            int code = json_unpack_ex(
                    stage_desc, &j_error, 0,
                    "{s?b, s?b, s?b, s?b, s?b, s?b, s?b, s?b, s?b, s?b, s?b, s?b, s?b, s?b, s:i, "
                    " s?f, s?f, s?f, s?f, s?f, s?f, s?f}",
                    "bonds", &stage_prms->bonds,
                    "angles", &stage_prms->angles,
                    "dihedrals", &stage_prms->dihedrals,
                    "impropers", &stage_prms->impropers,
                    "vdw", &stage_prms->vdw,
                    "vdw03", &stage_prms->vdw03,
                    "eleng", &stage_prms->eleng,
                    "elengs03", &stage_prms->elengs03,
                    "ace", &stage_prms->ace,
                    "gbsa", &stage_prms->gbsa,
                    "fix_receptor", &stage_fix_rec,
                    "fix_ligand", &stage_fix_lig,
                    "score_only", &stage_prms->score_only,
                    "pull_ligand_away", &stage_prms->pull_ligand_away,
                    "nsteps", &stage_prms->nsteps,

                    // float options
                    "tol", &stage_prms->tol,
                    "nbcut", &stage_prms->nbcut,
                    "ace_efac", &stage_prms->ace_efac,
                    "scale_vdw_s03", &stage_prms->scale_vdw_s03,
                    "scale_coul_s03", &stage_prms->scale_coul_s03,
                    "eeps", &stage_prms->eeps,
                    "gbcut", &stage_prms->gbcut
                    );

            if (code != 0) {
                JSON_ERR_MSG(j_error, "Couldn't parse stage flags in json");
                error = true;
                break;
            }
            if (stage_prms->nsteps < 1) {
                ERR_MSG("Field nsteps must be positive (%i <= 0)", stage_prms->nsteps);
                error = true;
                break;
            }

            if ((stage_fix_rec || stage_fix_lig) && !opts.separate) {
                ERR_MSG("You can't provide fix-* flags in a single file mode");
                error = true;
                break;
            }

            if ((!opts.separate) && stage_prms->pull_ligand_away) {
                ERR_MSG("Option 'pull_ligand_away' can be used only in rec/lig mode");
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
            // TODO: Dangling pointers can be left here, needs fixing
            if (stage_fixed) {
                DEBUG_MSG("Creating fixed atoms for stage");
                _fixed_setup_multi_free(&stage_prms->fixed);
                stage_prms->fixed = _fixed_setup_multi_read_json(stage_fixed, opts.num_models);
                if (!stage_prms->fixed) {
                    ERR_MSG("Couldn't parse fixed atoms from %s", opts.setup_json);
                    error = true;
                    break;
                }
            } else if (stage_fix_rec) {
                _fixed_setup_multi_free(&stage_prms->fixed);
                stage_prms->fixed = _fixed_setup_multi_atom_range(opts.num_models, 0, opts.rec_natoms);
            } else if (stage_fix_lig) {
                _fixed_setup_multi_free(&stage_prms->fixed);
                stage_prms->fixed = _fixed_setup_multi_atom_range(opts.num_models, opts.rec_natoms, opts.rec_natoms + opts.lig_natoms);
            }

            // Pairsprings
            json_t* stage_pairsprings = json_object_get(stage_desc, "pairsprings");
            if (stage_pairsprings) {
                DEBUG_MSG("Creating pairsprings for stage");
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
                DEBUG_MSG("Creating pointsprings for stage");
                stage_prms->sprst_points = _pointsprings_setup_read_json(stage_pointsprings);
                if (!stage_prms->sprst_points) {
                    ERR_MSG("Couldn't parse poinsprings from %s", opts.setup_json);
                    error = true;
                    break;
                }
            }

            // Density
            json_t* stage_density = json_object_get(stage_desc, "density");
            if (stage_density) {
                DEBUG_MSG("Creating density for stage");
                stage_prms->density = _density_setup_read_json(stage_density);
                if (!stage_prms->density) {
                    ERR_MSG("Couldn't parse density from %s", opts.setup_json);
                    error = true;
                    break;
                }
            }

            // NOE
            json_t* stage_noe = json_object_get(stage_desc, "noe");
            if (stage_noe) {
#ifndef NOE
                ERR_MSG("Rebuild NRGmin with -DNOE=ON to use NOE spectrum fitting");
                error = true;
                break;
#else

                DEBUG_MSG("Creating NOE for stage");
                stage_prms->nmr = _noe_setup_read_json(stage_noe);
                if (!stage_prms->nmr) {
                    ERR_MSG("Couldn't parse NOE from %s", opts.setup_json);
                    error = true;
                    break;
                }
#endif
            }
        }

        if (error) {
            energy_prms_free(&all_stage_prms, nstages);
            json_decref(setup_root);
            return false;
        }
    }

    if (setup_root) {
        json_decref(setup_root);
    }

    *result_energy_prm = all_stage_prms;
    *result_nstages = nstages;

    return true;
}


void pull_ligand_away(struct mol_atom_group* ag, size_t rec_size, double dist)
{
    struct mol_vector3 rec_com = {0., 0., 0.};
    for (size_t i = 0; i < rec_size; i++) {
        MOL_VEC_ADD(rec_com, rec_com, ag->coords[i]);
    }
    MOL_VEC_DIV_SCALAR(rec_com, rec_com, rec_size);

    struct mol_vector3 lig_com = {0., 0., 0.};
    for (size_t i = rec_size; i < ag->natoms; i++) {
        MOL_VEC_ADD(lig_com, lig_com, ag->coords[i]);
    }
    MOL_VEC_DIV_SCALAR(lig_com, lig_com, ag->natoms - rec_size);

    struct mol_vector3 tv;
    MOL_VEC_SUB(tv, lig_com, rec_com);
    double mult = dist / (sqrt(MOL_VEC_SQ_NORM(tv)) + 1e-6);
    MOL_VEC_MULT_SCALAR(tv, tv, mult);

    for (size_t i = rec_size; i < ag->natoms; i++) {
        MOL_VEC_ADD(ag->coords[i], ag->coords[i], tv);
    }
}