#include <stdio.h>
#include <check.h>
#include <errno.h>
#include <math.h>

#include "parse_options.h"
#include "potentials.h"


const double tolerance_strict = 1e-9;

#ifndef ck_assert_double_eq_tol
#define ck_assert_double_eq_tol(val, ref, tol) \
	ck_assert(fabs(((val))-((ref))) < (tol))
#endif

#ifndef ck_assert_double_eq_def
#define ck_assert_double_eq_def(val, ref) \
    ck_assert_double_eq_tol(val, ref, tolerance_strict)
#endif


START_TEST(test_check_getopt_success)
{
	bool error;
	struct options opts;

	switch (_i) {
	    case 0:
            opts = parse_args(3, (char *[]) {"sham", "sham", "-h"}, &error);
            ck_assert(!error);
            ck_assert(opts.help);
            break;
        case 1:
            opts = parse_args(10, (char *[]) {"sham", "sham",
                                          "--pdb", "test.pdb",
                                          "--psf", "test.psf",
                                          "--rtf", "test.rtf",
                                          "--prm", "test.prm"}, &error);
            ck_assert(!error);
            ck_assert_str_eq(opts.pdb, "test.pdb");
            break;
        case 2:
            opts = parse_args(4, (char *[]) {"sham", "sham",
                                              "--json", "test.json"}, &error);
            ck_assert(!error);
            ck_assert_str_eq(opts.json, "test.json");
            break;
        case 3:
            opts = parse_args(18, (char *[]) {"sham", "sham",
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
            opts = parse_args(12, (char *[]) {"sham", "sham",
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
            opts = parse_args(12, (char *[]) {"sham", "sham",
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
            opts = parse_args(11, (char *[]) {"sham", "sham",
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
    }
}
END_TEST


START_TEST(test_check_getopt_failure)
{
    bool error;
    struct options opts;

    switch (_i) {
        case 0:
            opts = parse_args(3, (char *[]) {"sham", "sham", "-a"}, &error);
            ck_assert(error);
            break;
        case 1:
            opts = parse_args(8, (char *[]) {"sham", "sham",
                                              "--pdb", "test.pdb",
                                              "--rtf", "test.rtf",
                                              "--prm", "test.prm"}, &error);
            ck_assert(error);
            break;
        case 2:
            opts = parse_args(6, (char *[]) {"sham", "sham",
                                             "--json", "test.json",
                                             "--rec-pdb", "test"}, &error);
            ck_assert(error);
            break;
        case 3:
            opts = parse_args(6, (char *[]) {"sham", "sham",
                                              "--json", "test.json",
                                              "--nsteps", "-1"}, &error);
            ck_assert(error);
            break;
        case 4:
            opts = parse_args(10, (char *[]) {"sham", "sham",
                                              "--rec-pdb", "test.pdb",
                                              "--rec-rtf", "test.rtf",
                                              "--rec-prm", "test.prm",
                                              "--lig-json", "test.json"}, &error);
            ck_assert(error);
            break;
        case 5:
            opts = parse_args(10, (char *[]) {"sham", "sham",
                                              "--rec-json", "test.json",
                                              "--lig-pdb", "test.pdb",
                                              "--lig-psf", "test.psf",
                                              "--lig-prm", "test.prm"}, &error);
            ck_assert(error);
            break;
        case 6:
            opts = parse_args(5, (char *[]) {"sham", "sham",
                                              "--json", "test.json",
                                              "--fix-receptor"}, &error);
            ck_assert(error);
            break;
    }
}
END_TEST


void setup_real(void)
{

}

void teardown_real(void)
{

}

Suite *lists_suite(void)
{
    Suite *suite = suite_create("getopt");
    TCase *tcase_real = tcase_create("real");
    tcase_add_checked_fixture(tcase_real, setup_real, teardown_real);
    tcase_add_loop_test(tcase_real, test_check_getopt_success, 0, 7);
    tcase_add_loop_test(tcase_real, test_check_getopt_failure, 0, 7);
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