#include <stdio.h>
#include <check.h>

#include "parse_options.h"
#include "setup.h"


char msg[512];

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
            ck_assert_double_eq_tol(prms->sprst_pairs->springs->weight, 5.0, 10e-3);
            ck_assert_double_eq_tol(prms->sprst_pairs->springs->error, 2.0, 10e-3);
            ck_assert_double_eq_tol(prms->sprst_pairs->springs->length, 5.0, 10e-3);

            _compare_arrays_size_t(
                    prms->sprst_pairs->springs->atoms,
                    (size_t[]){101, 3598},
                    2);

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

#ifdef NOE
            ck_assert(energy_prms_populate_from_options(&prms, &nstages, opts));
            ck_assert_ptr_nonnull(prms);
            ck_assert_double_eq_tol(prms->nmr->weight, 1000, 10e-3);
            ck_assert_double_eq_tol(prms->nmr->power, 1. / 6., 10e-3);
            ck_assert_int_eq(prms->nmr->spec->size, 2);
            ck_assert_double_eq_tol(prms->nmr->spec->exp[1], -0.074589, 10e-3);
            energy_prms_free(&prms, nstages);
#else
            ck_assert(!energy_prms_populate_from_options(&prms, &nstages, opts));
            ck_assert_ptr_null(prms);
#endif
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

        case 10:
            // Pass fixed pdb
            opts = options_get_default();
            opts.fixed_pdb = "BACE_4_lig_far.pdb";

            ck_assert(energy_prms_populate_from_options(&prms, &nstages, opts));
            ck_assert_ptr_nonnull(prms);
            ck_assert_int_eq(prms[0].fixed->natoms, 42);
            _compare_arrays_size_t(prms[0].fixed->atoms, (size_t[]){3598, 3599, 3600, 3601, 3602, 3603}, 6);

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

#ifdef NOE
            ck_assert(energy_prms_populate_from_options(&prms, &nstages, opts));
            ck_assert_ptr_nonnull(prms);
            ck_assert_int_eq(prms[0].nmr->spec->size, 2);
            energy_prms_free(&prms, nstages);
#else
            ck_assert(!energy_prms_populate_from_options(&prms, &nstages, opts));
            ck_assert_ptr_null(prms);
#endif
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
#ifdef NOE
            opts.setup_json = "setup_all.json";
#else
            opts.setup_json = "setup_all_no_noe.json";
#endif
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

#ifdef NOE
            ck_assert_ptr_null(prms[0].nmr);
            ck_assert_ptr_nonnull(prms[1].nmr);
            ck_assert_ptr_null(prms[2].nmr);
#endif

            ck_assert_ptr_null(prms[0].sprst_pairs);
            ck_assert_ptr_nonnull(prms[1].sprst_pairs);
            ck_assert_ptr_nonnull(prms[2].sprst_pairs);

            // Check some random values from the setup to make sure it was parsed properly
#ifdef NOE
            ck_assert_double_eq_tol(prms[1].nmr->power, 1./3., 10e-3);
            ck_assert_double_eq_tol(prms[1].nmr->weight, 1000, 10e-3);
            ck_assert_ptr_nonnull(prms[1].nmr->spec);
            ck_assert_int_eq(prms[1].nmr->spec->size, 2);
#endif
            ck_assert_int_eq(prms[2].sprst_points->springs[0].atoms[3], 4);
            ck_assert_int_eq(prms[1].sprst_pairs->springs[0].atoms[1], 5);

            energy_prms_free(&prms, nstages);
            break;

        default:
            ck_assert(false);
    }
}


Suite *lists_suite(void)
{
    Suite *suite = suite_create("functionality");
    TCase *tcase_real = tcase_create("real");
    //tcase_add_checked_fixture(tcase_real, setup_real, teardown_real);
    tcase_add_loop_test(tcase_real, test_check_getopt_success, 0, 8);
    tcase_add_loop_test(tcase_real, test_check_getopt_failure, 0, 8);
    tcase_add_loop_test(tcase_real, test_mol_atom_group_list_from_options, 0, 7);
    tcase_add_loop_test(tcase_real, test_energy_prm_from_flags, 0, 11);
    tcase_add_loop_test(tcase_real, test_energy_prm_from_json, 0, 12);
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