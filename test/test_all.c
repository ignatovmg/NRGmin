#include <stdio.h>
#include <check.h>
#include <jansson.h>

#include "mol2/benergy.h"
#include "mol2/pdb.h"
#include "mol2/vector.h"
#include "mol2/atom_group.h"
#include "utils.h"
#include "setup.h"
#include "parse_options.h"
#include "energy.h"

// These functions are defined only in newer versions of check
#ifndef ck_assert_double_eq_tol
#define ck_assert_double_eq_tol(x, y, tol) do { \
	sprintf(msg, "Calculated: %lf, Reference: %lf, Tolerance %lf (line %i)\n", x, y, tol, __LINE__); \
	ck_assert_msg(fabs((x) - (y)) <= (tol), msg); \
	} while(0);
#endif

#ifndef ck_assert_ptr_nonnull
#define ck_assert_ptr_nonnull(x) do { \
	sprintf(msg, "Pointer is NULL (line %i)\n", __LINE__); \
	ck_assert_msg((x) != NULL, msg); \
	} while(0);
#endif

#ifndef ck_assert_ptr_null
#define ck_assert_ptr_null(x) do { \
	sprintf(msg, "Pointer is not NULL (line %i)\n", __LINE__); \
	ck_assert_msg((x) == NULL, msg); \
	} while(0);
#endif

struct mol_atom_group *test_ag1, *test_ag2;
struct pairsprings_setup *sprst, *sj;
json_t *spj;
json_error_t *err_j;
/*static void _compare_arrays_int(const int* x, const int* y, const size_t len)
{
    for (size_t i = 0; i < len; i++) {
	    ck_assert_int_eq(x[i], y[i]);
	}
}*/

static void _compare_arrays_size_t(const size_t* x, const size_t* y, const size_t len)
{
    for (size_t i = 0; i < len; i++) {
        ck_assert_int_eq((int)x[i], (int)y[i]);
    }
}

/*static void _compare_arrays_double(const double* x, const double* y, const size_t len, const double tol)
{
    for (size_t i = 0; i < len; i++) {
        ck_assert_double_eq_tol(x[i], y[i], tol);
    }
}*/

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
                spring_set[counter].group1[i] = json_integer_value(json_array_get(g1, i)) - 1;
            } else {
                ERR_MSG("Index in group1 must be interger");
            }
        }
        for (i = 0; i < len2; ++i) {
            if (json_is_integer(json_array_get(g2, i))) {
                spring_set[counter].group2[i] = json_integer_value(json_array_get(g2, i)) - 1;
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

static void test_pairspring_gradients(
        struct mol_atom_group *mob_ag,
        const struct pairsprings_setup *sps,
        const double delta,
        const double tol)
{
    struct mol_vector3 *num_grad = calloc(mob_ag->natoms, sizeof(struct mol_vector3));

    mol_zero_gradients(mob_ag);
    lbfgsfloatval_t een_orig;
    pairspring_energy(sps, mob_ag, &een_orig);
    struct mol_vector3 *old_grad_ptr = mob_ag->gradients;
    double *vect = calloc(3*mob_ag->natoms, sizeof(double));
    for (int i = 0; i < mob_ag->natoms; i++)
    {
        vect[3*i] = old_grad_ptr[i].X;
        vect[3*i+1] = old_grad_ptr[i].Y;
        vect[3*i+2] = old_grad_ptr[i].Z;
    }

    mol_zero_gradients(mob_ag);

    double oldcrd;
    lbfgsfloatval_t een_x, een_y, een_z;
    for (size_t i = 0; i < mob_ag->natoms; i++) {
        een_x = 0;
        een_y = 0;
        een_z = 0;
        oldcrd = mob_ag->coords[i].X;
        mob_ag->coords[i].X = oldcrd + delta;
        pairspring_energy(sps, mob_ag, &een_x);
        num_grad[i].X = (een_orig - een_x) / delta;
        mob_ag->coords[i].X = oldcrd;

        oldcrd = mob_ag->coords[i].Y;
        mob_ag->coords[i].Y = oldcrd + delta;
        pairspring_energy(sps, mob_ag, &een_y);
        num_grad[i].Y = (een_orig - een_y) / delta;
        mob_ag->coords[i].Y = oldcrd;

        oldcrd = mob_ag->coords[i].Z;
        mob_ag->coords[i].Z = oldcrd + delta;
        pairspring_energy(sps, mob_ag, &een_z);
        num_grad[i].Z = (een_orig - een_z) / delta;
        mob_ag->coords[i].Z = oldcrd;
    }

    char msg[256];
    for (size_t i = 0; i < mob_ag->natoms; i++) {
        sprintf(msg,
                "\n(atom: %ld) calc: (%lf, %lf, %lf): numerical: (%lf, %lf, %lf)\n",
                i,
                vect[3*i], vect[3*i+1], vect[3*i+2],
                num_grad[i].X, num_grad[i].Y, num_grad[i].Z);

        ck_assert_msg(fabs(vect[3*i] - num_grad[i].X) < tol, msg);
        ck_assert_msg(fabs(vect[3*i+1] - num_grad[i].Y) < tol, msg);
        ck_assert_msg(fabs(vect[3*i+2] - num_grad[i].Z) < tol, msg);
    }
    free(num_grad);
    free(vect);
}

START_TEST(test_check_getopt_success)
{
	bool error;
	struct options opts;

	switch (_i) {
	    case 0:
            opts = options_populate_from_argv(3, (char *[]) {"sham", "sham", "-h"}, &error);
            ck_assert(!error);
            ck_assert(opts.help);
            break;
        case 1:
            opts = options_populate_from_argv(10, (char *[]) {"sham", "sham",
                                                              "--pdb", "test.pdb",
                                                              "--psf", "test.psf",
                                                              "--rtf", "test.rtf",
                                                              "--prm", "test.prm"}, &error);
            ck_assert(!error);
            ck_assert_str_eq(opts.pdb, "test.pdb");
            break;
        case 2:
            opts = options_populate_from_argv(4, (char *[]) {"sham", "sham",
                                                             "--json", "test.json"}, &error);
            ck_assert(!error);
            ck_assert_str_eq(opts.json, "test.json");
            break;
        case 3:
            opts = options_populate_from_argv(18, (char *[]) {"sham", "sham",
                                                              "--rec-pdb", "test.pdb",
                                                              "--rec-psf", "test.psf",
                                                              "--rec-rtf", "test.rtf",
                                                              "--rec-prm", "test.prm",
                                                              "--lig-pdb", "test.pdb",
                                                              "--lig-psf", "test.psf",
                                                              "--lig-rtf", "test.rtf",
                                                              "--lig-prm", "test.prm"}, &error);
            ck_assert(!error);
            ck_assert_str_eq(opts.rec_pdb, "test.pdb");
            ck_assert_str_eq(opts.lig_pdb, "test.pdb");
            ck_assert(opts.separate);
            break;
        case 4:
            opts = options_populate_from_argv(12, (char *[]) {"sham", "sham",
                                                              "--rec-pdb", "test.pdb",
                                                              "--rec-psf", "test.psf",
                                                              "--rec-rtf", "test.rtf",
                                                              "--rec-prm", "test.prm",
                                                              "--lig-json", "test.json"}, &error);
            ck_assert(!error);
            ck_assert_str_eq(opts.rec_pdb, "test.pdb");
            ck_assert_str_eq(opts.lig_json, "test.json");
            ck_assert(opts.separate);
            break;
        case 5:
            opts = options_populate_from_argv(12, (char *[]) {"sham", "sham",
                                                              "--rec-json", "test.json",
                                                              "--lig-pdb", "test.pdb",
                                                              "--lig-psf", "test.psf",
                                                              "--lig-rtf", "test.rtf",
                                                              "--lig-prm", "test.prm"}, &error);
            ck_assert(!error);
            ck_assert_str_eq(opts.rec_json, "test.json");
            ck_assert_str_eq(opts.lig_pdb, "test.pdb");
            ck_assert(opts.separate);
            break;
        case 6:
            opts = options_populate_from_argv(11, (char *[]) {"sham", "sham",
                                                              "--rec-json", "test.json",
                                                              "--lig-json", "test.json",
                                                              "--nsteps", "500",
                                                              "--pair-springs-txt", "test",
                                                              "--fix-receptor"}, &error);
            ck_assert(!error);
            ck_assert_str_eq(opts.rec_json, "test.json");
            ck_assert_str_eq(opts.lig_json, "test.json");
            ck_assert_int_eq(opts.nsteps, 500);
            ck_assert_str_eq(opts.pair_springs_txt, "test");
            ck_assert(opts.fix_receptor);
            ck_assert(!opts.fix_ligand);
            break;
        case 7:
            opts = options_populate_from_argv(4, (char *[]) {"sham", "sham",
                                                             "--setup-json", "setup_global.json"}, &error);
            ck_assert(!error);
            ck_assert_str_eq(opts.json, "test.json");
            ck_assert_str_eq(opts.pdb, "test.pdb");
            ck_assert_int_eq(opts.nsteps, 100);
            ck_assert(!opts.dihedrals);
            ck_assert(!opts.bonds);
            break;
    }
}
END_TEST


START_TEST(test_check_getopt_failure)
{
    bool error;

    switch (_i) {
        case 0:
            options_populate_from_argv(3, (char *[]) {"sham", "sham", "-a"}, &error);
            ck_assert(error);
            break;
        case 1:
            options_populate_from_argv(8, (char *[]) {"sham", "sham",
                                                      "--pdb", "test.pdb",
                                                      "--rtf", "test.rtf",
                                                      "--prm", "test.prm"}, &error);
            ck_assert(error);
            break;
        case 2:
            options_populate_from_argv(6, (char *[]) {"sham", "sham",
                                                      "--json", "test.json",
                                                      "--rec-pdb", "test"}, &error);
            ck_assert(error);
            break;
        case 3:
            options_populate_from_argv(6, (char *[]) {"sham", "sham",
                                                      "--json", "test.json",
                                                      "--nsteps", "-1"}, &error);
            ck_assert(error);
            break;
        case 4:
            options_populate_from_argv(10, (char *[]) {"sham", "sham",
                                                       "--rec-pdb", "test.pdb",
                                                       "--rec-rtf", "test.rtf",
                                                       "--rec-prm", "test.prm",
                                                       "--lig-json", "test.json"}, &error);
            ck_assert(error);
            break;
        case 5:
            options_populate_from_argv(10, (char *[]) {"sham", "sham",
                                                       "--rec-json", "test.json",
                                                       "--lig-pdb", "test.pdb",
                                                       "--lig-psf", "test.psf",
                                                       "--lig-prm", "test.prm"}, &error);
            ck_assert(error);
            break;
        case 6:
            options_populate_from_argv(5, (char *[]) {"sham", "sham",
                                                      "--json", "test.json",
                                                      "--fix-receptor"}, &error);
            ck_assert(error);
            break;
        case 7:
            options_populate_from_argv(4, (char *[]) {"sham", "sham",
                                                      "--setup-json", "setup_global_unknown_args.json"}, &error);
            ck_assert(error);
            break;
    }
}
END_TEST


START_TEST(test_mol_atom_group_list_from_options)
{
    struct options opts;
    struct mol_atom_group_list *ag_list;

	switch (_i) {
	    case 0:
	        // Pass single pdb
            opts = options_get_default();
            opts.pdb = "BACE_4_rec.pdb";
            opts.psf = "BACE_4_rec.psf";
            opts.prm = "BACE_4_rec_parm.prm";
            opts.rtf = "BACE_4_rec_pdbamino.rtf";

            ag_list = mol_atom_group_list_from_options(&opts);
            ck_assert_ptr_nonnull(ag_list);
            ck_assert_int_eq(opts.rec_natoms, 0);
            ck_assert_int_eq(opts.lig_natoms, 0);
            ck_assert_int_eq(ag_list->members[0].natoms, 3598);
            mol_atom_group_list_free(ag_list);
            break;

	    case 1:
	        // Pass merge rec and lig
            opts = options_get_default();
            opts.rec_pdb = "BACE_4_rec.pdb";
            opts.rec_psf = "BACE_4_rec.psf";
            opts.rec_prm = "BACE_4_rec_parm.prm";
            opts.rec_rtf = "BACE_4_rec_pdbamino.rtf";
            opts.lig_json = "BACE_4_lig.json";
            opts.separate = true;

            ag_list = mol_atom_group_list_from_options(&opts);
            ck_assert_ptr_nonnull(ag_list);
            ck_assert_int_eq(opts.rec_natoms, 3598);
            ck_assert_int_eq(opts.lig_natoms, 42);
            ck_assert_int_eq(ag_list->members[0].natoms, 3598 + 42);
            mol_atom_group_list_free(ag_list);
            break;

	    case 2:
	        // Pass single json
            opts = options_get_default();
            opts.json = "BACE_4_lig.json";

            ag_list = mol_atom_group_list_from_options(&opts);
            ck_assert_ptr_nonnull(ag_list);
            ck_assert_int_eq(opts.rec_natoms, 0);
            ck_assert_int_eq(opts.lig_natoms, 0);
            ck_assert_int_eq(ag_list->members[0].natoms, 42);
            mol_atom_group_list_free(ag_list);
            break;

        case 3:
            // Fail trying to merge geometry from json and pdb with different atom numbers
            opts = options_get_default();
            opts.pdb = "BACE_4_rec.pdb";
            opts.json = "BACE_4_lig.json";
            opts.separate = false;

            ag_list = mol_atom_group_list_from_options(&opts);
            ck_assert_ptr_null(ag_list);
            break;

        case 4:
            // Fail with missing prm file
            opts = options_get_default();
            opts.pdb = "BACE_4_rec.pdb";
            opts.psf = "BACE_4_rec.psf";
            opts.rtf = "BACE_4_rec_pdbamino.rtf";

            ag_list = mol_atom_group_list_from_options(&opts);
            ck_assert_ptr_null(ag_list);
            break;

        case 5:
            // Pass with missing prms and score_only flag set
            opts = options_get_default();
            opts.pdb = "BACE_4_rec.pdb";
            opts.score_only = true;

            ag_list = mol_atom_group_list_from_options(&opts);
            ck_assert_ptr_nonnull(ag_list);
            break;

        case 6:
            // Pass with two pdb models and json geometry
            opts = options_get_default();
            opts.pdb = "BACE_4_lig_2models_close.pdb";
            opts.json = "BACE_4_lig.json";

            ag_list = mol_atom_group_list_from_options(&opts);
            ck_assert_ptr_nonnull(ag_list);
            ck_assert_int_eq(ag_list->size, 2);
            break;

        /*case 6:
            // Fail with wrong prm file (this produces segfault, because of buggy mol_atom_group_read_geometry)
            opts = options_get_default();
            opts.pdb = "BACE_4_rec.pdb";
            opts.psf = "BACE_4_rec.pdb";
            opts.prm = "BACE_4_rec_parm.prm";
            opts.rtf = "BACE_4_rec_pdbamino.rtf";
            opts.separate = false;

            ag_list = mol_atom_group_list_from_options(&opts);
            ck_assert_ptr_null(ag_list);
            */
    }
}


START_TEST(test_energy_prm_from_flags)
{
    struct options opts;
    struct energy_prms* prms;
    size_t nstages;

    switch (_i) {
        case 0:
            // Pass fix rec
            opts = options_get_default();
            opts.separate = true;
            opts.rec_natoms = 1000;
            opts.lig_natoms = 100;
            opts.fix_receptor = true;

            ck_assert(energy_prms_populate_from_options(&prms, &nstages, opts));
            ck_assert_ptr_nonnull(prms);
            ck_assert_int_eq(prms[0].fixed->natoms, 1000);
            _compare_arrays_size_t(prms[0].fixed->atoms, (size_t[]){0, 1, 2, 3, 4}, 5);

            energy_prms_free(&prms, nstages);
            break;

        case 1:
            // Pass fix lig
            opts = options_get_default();
            opts.separate = true;
            opts.rec_natoms = 1000;
            opts.lig_natoms = 100;
            opts.fix_ligand = true;

            ck_assert(energy_prms_populate_from_options(&prms, &nstages, opts));
            ck_assert_ptr_nonnull(prms);
            ck_assert_int_eq(prms[0].fixed->natoms, 100);
            _compare_arrays_size_t(prms[0].fixed->atoms, (size_t[]){1000, 1001, 1002, 1003, 1004}, 5);

            energy_prms_free(&prms, nstages);
            break;

        case 2:
            // Fail fix lig and rec
            opts = options_get_default();
            opts.separate = true;
            opts.rec_natoms = 1000;
            opts.lig_natoms = 100;
            opts.fix_ligand = true;
            opts.fix_receptor = true;

            ck_assert(!energy_prms_populate_from_options(&prms, &nstages, opts));
            ck_assert_ptr_null(prms);
            break;

        case 3:
            // Pass with fix pdb
            opts = options_get_default();
            opts.fixed_pdb = "BACE_4_rec.pdb";

            ck_assert(energy_prms_populate_from_options(&prms, &nstages, opts));
            ck_assert_ptr_nonnull(prms);
            ck_assert_int_eq(prms->fixed->natoms, 3598);
            energy_prms_free(&prms, nstages);
            break;

        case 4:
            // Pass with pointsprings
            opts = options_get_default();
            opts.point_springs_txt = "pointsprings_xl.txt";

            ck_assert(energy_prms_populate_from_options(&prms, &nstages, opts));
            ck_assert_ptr_nonnull(prms);
            ck_assert_int_eq(prms->sprst_points->nsprings, 1);
            ck_assert_int_eq(prms->sprst_points->springs[0].natoms, 24);
            ck_assert_double_eq_tol(prms->sprst_points->springs[0].weight, 8, 10e-3);
            ck_assert_double_eq_tol(prms->sprst_points->springs[0].X0, 30.381, 10e-3);
            ck_assert_double_eq_tol(prms->sprst_points->springs[0].Y0, 5.877, 10e-3);
            ck_assert_double_eq_tol(prms->sprst_points->springs[0].Z0, 14.905, 10e-3);

            // Compare first 6 IDs
            _compare_arrays_size_t(
                    prms->sprst_points->springs->atoms,
                    (size_t[]){3598, 3599, 3600, 3601, 3602, 3603},
                    6);

            energy_prms_free(&prms, nstages);
            break;

        case 5:
            // Fail with wrong pointsprings
            opts = options_get_default();
            opts.point_springs_txt = "pointsprings_xl_invalid.txt";

            ck_assert(!energy_prms_populate_from_options(&prms, &nstages, opts));
            ck_assert_ptr_null(prms);
            break;

        case 6:
            // Pass with pairsprings
            opts = options_get_default();
            opts.pair_springs_txt = "pairsprings.txt";

            ck_assert(energy_prms_populate_from_options(&prms, &nstages, opts));
            ck_assert_ptr_nonnull(prms);
            ck_assert_int_eq(prms->sprst_pairs->nsprings, 1);
            ck_assert_int_eq(prms->sprst_pairs->springs->group_size1, 1);
            ck_assert_int_eq(prms->sprst_pairs->springs->group_size2, 1);
            ck_assert_double_eq_tol(prms->sprst_pairs->springs->weight, 5.0, 10e-3);
            ck_assert_double_eq_tol(prms->sprst_pairs->springs->lerror, 2.0, 10e-3);
            ck_assert_double_eq_tol(prms->sprst_pairs->springs->rerror, 2.0, 10e-3);
            ck_assert_double_eq_tol(prms->sprst_pairs->springs->distance, 5.0, 10e-3);
            ck_assert_int_eq(prms->sprst_pairs->springs->average, 0);
            ck_assert_int_eq(prms->sprst_pairs->springs->potential, 2);
            ck_assert_int_eq(prms->sprst_pairs->springs->group1[0], 101);
            ck_assert_int_eq(prms->sprst_pairs->springs->group2[0], 3598);

            energy_prms_free(&prms, nstages);
            break;

        case 7:
            // Fail with pairsprings
            opts = options_get_default();
            opts.pair_springs_txt = "pairsprings_invalid.txt";

            ck_assert(!energy_prms_populate_from_options(&prms, &nstages, opts));
            ck_assert_ptr_null(prms);
            break;

        case 8:
            // Pass with NOE
            opts = options_get_default();
            opts.noe_json = "noe.json";

            ck_assert(energy_prms_populate_from_options(&prms, &nstages, opts));
            ck_assert_ptr_nonnull(prms);
            ck_assert_double_eq_tol(prms->nmr->weight, 1000, 10e-3);
            ck_assert_double_eq_tol(prms->nmr->power, 1. / 6., 10e-3);
            ck_assert_int_eq(prms->nmr->spec->size, 2);
            ck_assert_double_eq_tol(prms->nmr->spec->exp[1], -0.074589, 10e-3);
            energy_prms_free(&prms, nstages);
            break;

        case 9:
            // Pass density
            opts = options_get_default();
            opts.density_json = "density.json";

            ck_assert(energy_prms_populate_from_options(&prms, &nstages, opts));
            ck_assert_ptr_nonnull(prms);
            ck_assert_double_eq_tol(prms[0].density->weight, 1000, 0.001);
            ck_assert_double_eq_tol(prms[0].density->prms.radius, 2, 0.001);
            ck_assert_int_eq(prms[0].density->ag->natoms, 42);
            energy_prms_free(&prms, nstages);
            break;
    }
}


START_TEST(test_energy_prm_from_json)
{
    struct options opts;
    struct energy_prms* prms;
    size_t nstages;

    switch (_i) {
        case 0:
            // Pass fix rec
            opts = options_get_default();
            opts.separate = true;
            opts.rec_natoms = 1000;
            opts.setup_json = "setup_fixed_flag.json";

            ck_assert(energy_prms_populate_from_options(&prms, &nstages, opts));
            ck_assert_ptr_nonnull(prms);
            ck_assert_int_eq(prms->fixed->natoms, opts.rec_natoms);
            energy_prms_free(&prms, nstages);
            break;

        case 1:
            // Fail fix rec
            opts = options_get_default();
            opts.separate = false;
            opts.setup_json = "setup_fixed_flag.json";

            ck_assert(!energy_prms_populate_from_options(&prms, &nstages, opts));
            ck_assert_ptr_null(prms);
            break;

        case 2:
            // Pass fix rec
            opts = options_get_default();
            opts.setup_json = "setup_fixed_json.json";

            ck_assert(energy_prms_populate_from_options(&prms, &nstages, opts));
            ck_assert_ptr_nonnull(prms);
            ck_assert_int_eq(prms->fixed->natoms, 7);
            _compare_arrays_size_t((size_t[]){0,1,2,3,4,8,2000}, prms->fixed->atoms, 7);
            energy_prms_free(&prms, nstages);
            break;

        case 3:
            // Pass pairsprings
            opts = options_get_default();
            opts.setup_json = "setup_pairsprings.json";

            ck_assert(energy_prms_populate_from_options(&prms, &nstages, opts));
            ck_assert_ptr_nonnull(prms);

            energy_prms_free(&prms, nstages);
            break;

        case 4:
            // Fail pairsprings
            opts = options_get_default();
            opts.setup_json = "setup_pairsprings_invalid.json";

            ck_assert(!energy_prms_populate_from_options(&prms, &nstages, opts));
            ck_assert_ptr_null(prms);
            energy_prms_free(&prms, nstages);
            break;

        case 5:
            // Pass pointsprings
            opts = options_get_default();
            opts.setup_json = "setup_pointsprings.json";

            ck_assert(energy_prms_populate_from_options(&prms, &nstages, opts));
            ck_assert_ptr_nonnull(prms);
            energy_prms_free(&prms, nstages);
            break;

        case 6:
            // Fail pointsprings
            opts = options_get_default();
            opts.setup_json = "setup_pointsprings_invalid.json";

            ck_assert(!energy_prms_populate_from_options(&prms, &nstages, opts));
            ck_assert_ptr_null(prms);
            energy_prms_free(&prms, nstages);
            break;

        case 7:
            // Pass noe
            opts = options_get_default();
            opts.setup_json = "setup_noe.json";

            ck_assert(energy_prms_populate_from_options(&prms, &nstages, opts));
            ck_assert_ptr_nonnull(prms);
            ck_assert_int_eq(prms[0].nmr->spec->size, 2);
            energy_prms_free(&prms, nstages);
            break;

        case 8:
            // Fail noe
            opts = options_get_default();
            opts.setup_json = "setup_noe_invalid.json";

            ck_assert(!energy_prms_populate_from_options(&prms, &nstages, opts));
            ck_assert_ptr_null(prms);
            energy_prms_free(&prms, nstages);
            break;

        case 9:
            // Pass density
            opts = options_get_default();
            opts.setup_json = "setup_density.json";

            ck_assert(energy_prms_populate_from_options(&prms, &nstages, opts));
            ck_assert_ptr_nonnull(prms);
            ck_assert_double_eq_tol(prms[0].density->weight, 1000, 0.001);
            ck_assert_double_eq_tol(prms[0].density->prms.radius, 2, 0.001);
            ck_assert_int_eq(prms[0].density->ag->natoms, 42);
            energy_prms_free(&prms, nstages);
            break;

        case 10:
            // Fail density
            opts = options_get_default();
            opts.setup_json = "setup_density_invalid.json";

            ck_assert(!energy_prms_populate_from_options(&prms, &nstages, opts));
            ck_assert_ptr_null(prms);
            energy_prms_free(&prms, nstages);
            break;

        case 11:
            // Pass generic setup
            opts = options_get_default();
            opts.separate = true;
            opts.rec_natoms = 1000;
            opts.lig_natoms = 100;
            opts.setup_json = "setup_all.json";

            ck_assert(energy_prms_populate_from_options(&prms, &nstages, opts));
            ck_assert_ptr_nonnull(prms);
            ck_assert_int_eq(nstages, 3);

            // Check some pointers for different stages
            ck_assert(!prms[0].vdw);
            ck_assert(prms[1].vdw);
            ck_assert(prms[2].vdw);

            ck_assert_ptr_nonnull(prms[0].fixed);
            ck_assert_ptr_nonnull(prms[1].fixed);
            ck_assert_ptr_null(prms[2].fixed);

            ck_assert_ptr_null(prms[0].nmr);
            ck_assert_ptr_nonnull(prms[1].nmr);
            ck_assert_ptr_null(prms[2].nmr);

            ck_assert_ptr_null(prms[0].sprst_pairs);
            ck_assert_ptr_nonnull(prms[1].sprst_pairs);
            ck_assert_ptr_nonnull(prms[2].sprst_pairs);

            // Check some random values from the setup to make sure it was parsed properly
            ck_assert_double_eq_tol(prms[1].nmr->power, 1./3., 10e-3);
            ck_assert_double_eq_tol(prms[1].nmr->weight, 1000, 10e-3);
            ck_assert_ptr_nonnull(prms[1].nmr->spec);
            ck_assert_int_eq(prms[1].nmr->spec->size, 2);
            ck_assert_int_eq(prms[2].sprst_points->springs[0].atoms[3], 4);
            ck_assert_int_eq(prms[1].sprst_pairs->springs[0].group2[0], 4);

            energy_prms_free(&prms, nstages);
            break;

        default:
            ck_assert(false);
    }
}


START_TEST(test_pairspring_penalty)
{
    mol_zero_gradients(test_ag1);
    test_pairspring_gradients(test_ag1, sprst, 10e-5, 10e-3);
}
END_TEST

START_TEST(test_pairspring_penalty2)
{   // test wit json input
    mol_zero_gradients(test_ag2);
    test_pairspring_gradients(test_ag2, sj, 10e-5, 10e-3);
}
END_TEST

void setup_real(void)
{
    sprst = _pairsprings_setup_read_txt("355w1.txt");
    spj = json_load_file("355w1.json", 0, err_j);
    sj = _pairsprings_setup_read_json(spj);
    test_ag1 = mol_read_pdb("30355_best.pdb");
    test_ag2 = mol_read_pdb("30355_best.pdb");
    test_ag1->gradients = calloc(test_ag1->natoms, sizeof(struct mol_vector3));
    test_ag2->gradients = calloc(test_ag2->natoms, sizeof(struct mol_vector3));
}

void teardown_real(void)
{
    mol_atom_group_free(test_ag1);
    mol_atom_group_free(test_ag2);
}

Suite *lists_suite(void)
{
    Suite *suite = suite_create("functionality");
    TCase *tcase_real = tcase_create("real");
    tcase_add_checked_fixture(tcase_real, setup_real, teardown_real);
    tcase_add_loop_test(tcase_real, test_check_getopt_success, 0, 8);
    tcase_add_loop_test(tcase_real, test_check_getopt_failure, 0, 8);
    tcase_add_loop_test(tcase_real, test_mol_atom_group_list_from_options, 0, 7);
    tcase_add_loop_test(tcase_real, test_energy_prm_from_flags, 0, 10);
    tcase_add_loop_test(tcase_real, test_energy_prm_from_json, 0, 12);
    tcase_add_test(tcase_real, test_pairspring_penalty);
    tcase_add_test(tcase_real, test_pairspring_penalty2);
    suite_add_tcase(suite, tcase_real);

    return suite;
}

int main(void)
{
    Suite *suite = lists_suite();
    SRunner *runner = srunner_create(suite);
    srunner_run_all(runner, CK_ENV);

    int number_failed = srunner_ntests_failed(runner);
    srunner_free(runner);
    return number_failed;
}
