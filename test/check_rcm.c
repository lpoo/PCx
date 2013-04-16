/* Test for Reverse Cuthill-McKee
 *
 * Authors: Raniere Silva <ra092767@ime.unicamp.br>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <check.h>
#include "../SRC/rcm.h"

/* Test for
 *
 * comp = [1 2 3 4 5];
 * */
START_TEST (test_CC_create00) {
  CC *comps;
  int comp[] = {1, 2, 3, 4, 5};
  int begin[] = {1, 2, 3, 4, 5, 6};
  int node[] = {1, 2, 3, 4, 5};
  comps = create_CC(5, 5, comp);
  fail_unless(comps->N == 5,
      "`N` not set correctly on creation");
  fail_unless(memcmp(comps->begin, begin, sizeof(begin)) == 0,
      "`begin` not set correctly on creation");
  fail_unless(memcmp(comps->node, node, sizeof(node)) == 0,
      "`node` not set correctly on creation");
  fail_unless(memcmp(comps->comp, comp, sizeof(comp)) == 0,
      "`comp` not set correctly on creation");
  free_CC(comps);
}
END_TEST

/* Test for
 *
 * comp = [1 1 1 1 1];
 * */
START_TEST (test_CC_create01) {
  CC *comps;
  int comp[] = {1, 1, 1, 1, 1};
  int begin[] = {1, 6};
  int node[] = {1, 2, 3, 4, 5};
  comps = create_CC(5, 1, comp);
  fail_unless(comps->N == 1,
      "`N` not set correctly on creation");
  fail_unless(memcmp(comps->begin, begin, sizeof(begin)) == 0,
      "`begin` not set correctly on creation");
  fail_unless(memcmp(comps->node, node, sizeof(node)) == 0,
      "`node` not set correctly on creation");
  fail_unless(memcmp(comps->comp, comp, sizeof(comp)) == 0,
      "`comp` not set correctly on creation");
  free_CC(comps);
}
END_TEST

/* Test for
 *
 * comp = [1 1 1 2 2];
 * */
START_TEST (test_CC_create02) {
  CC *comps;
  int comp[] = {1, 1, 1, 2, 2};
  int begin[] = {1, 4, 6};
  int node[] = {1, 2, 3, 4, 5};
  comps = create_CC(5, 2, comp);
  fail_unless(comps->N == 2,
      "`N` not set correctly on creation");
  fail_unless(memcmp(comps->begin, begin, sizeof(begin)) == 0,
      "`begin` not set correctly on creation");
  fail_unless(memcmp(comps->node, node, sizeof(node)) == 0,
      "`node` not set correctly on creation");
  fail_unless(memcmp(comps->comp, comp, sizeof(comp)) == 0,
      "`comp` not set correctly on creation");
  free_CC(comps);
}
END_TEST

/* Test for
 *
 * comp = [1 2 1 2 2];
 * */
START_TEST (test_CC_create03) {
  CC *comps;
  int comp[] = {1, 2, 1, 2, 2};
  int begin[] = {1, 3, 6};
  int node[] = {1, 3, 2, 4, 5};
  comps = create_CC(5, 2, comp);
  fail_unless(comps->N == 2,
      "`N` not set correctly on creation");
  fail_unless(memcmp(comps->begin, begin, sizeof(begin)) == 0,
      "`begin` not set correctly on creation");
  fail_unless(memcmp(comps->node, node, sizeof(node)) == 0,
      "`node` not set correctly on creation");
  fail_unless(memcmp(comps->comp, comp, sizeof(comp)) == 0,
      "`comp` not set correctly on creation");
  free_CC(comps);
}
END_TEST

Suite *CC_suite(void) {
  Suite *s = suite_create("CC");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_CC_create00);
  tcase_add_test(tc_core, test_CC_create01);
  tcase_add_test(tc_core, test_CC_create02);
  tcase_add_test(tc_core, test_CC_create03);
  suite_add_tcase(s, tc_core);

  return s;
}

/* Test for
 *
 * A = [1 0 0 0 0;
 *      0 1 0 0 0;
 *      0 0 1 0 0;
 *      0 0 0 1 0;
 *      0 0 0 0 1]
 */
START_TEST (test_fcc00) {
  int *xadj = NULL;
  int *adncy = NULL;
  CC *comps = fcc(5, xadj, adncy);
  int comp[] = {1, 2, 3, 4, 5};
  int begin[] = {1, 2, 3, 4, 5, 6};
  int node[] = {1, 2, 3, 4, 5};
  fail_unless(comps->N == 5,
      "`N` not set correctly on creation");
  fail_unless(memcmp(comps->begin, begin, sizeof(begin)) == 0,
      "`begin` not set correctly on creation");
  fail_unless(memcmp(comps->node, node, sizeof(node)) == 0,
      "`node` not set correctly on creation");
  fail_unless(memcmp(comps->comp, comp, sizeof(comp)) == 0,
      "`comp` not set correctly on creation");
  free_CC(comps);
}
END_TEST

/* Test for
 *
 * A = [1 1 1 1 1;
 *      1 1 0 0 0;
 *      1 0 1 0 0;
 *      1 0 0 1 0;
 *      1 0 0 0 1]
 */
START_TEST (test_fcc01) {
  int xadj[] = {1, 5, 6, 7, 8, 9};
  int adncy[] = {2, 3, 4, 5, 1, 1, 1, 1};
  CC *comps = fcc(5, xadj, adncy);
  int comp[] = {1, 1, 1, 1, 1};
  int begin[] = {1, 6};
  int node[] = {1, 2, 3, 4, 5};
  fail_unless(comps->N == 1,
      "`N` not set correctly on creation");
  fail_unless(memcmp(comps->begin, begin, sizeof(begin)) == 0,
      "`begin` not set correctly on creation");
  fail_unless(memcmp(comps->node, node, sizeof(node)) == 0,
      "`node` not set correctly on creation");
  fail_unless(memcmp(comps->comp, comp, sizeof(comp)) == 0,
      "`comp` not set correctly on creation");
  free_CC(comps);
}
END_TEST

/* Test for
 *
 * A = [1 1 1 0 0;
 *      1 1 0 0 0;
 *      1 0 1 0 0;
 *      0 0 0 1 1;
 *      0 0 0 1 1]
 */
START_TEST (test_fcc02) {
  int xadj[] = {1, 3, 4, 5, 6, 7};
  int adncy[] = {2, 3, 1, 1, 5, 4};
  CC *comps = fcc(5, xadj, adncy);
  int comp[] = {1, 1, 1, 2, 2};
  int begin[] = {1, 4, 6};
  int node[] = {1, 2, 3, 4, 5};
  fail_unless(comps->N == 2,
      "`N` not set correctly on creation");
  fail_unless(memcmp(comps->begin, begin, sizeof(begin)) == 0,
      "`begin` not set correctly on creation");
  fail_unless(memcmp(comps->node, node, sizeof(node)) == 0,
      "`node` not set correctly on creation");
  fail_unless(memcmp(comps->comp, comp, sizeof(comp)) == 0,
      "`comp` not set correctly on creation");
  free_CC(comps);
}
END_TEST

/* Test for
 *
 * A = [1 0 1 0 0;
 *      0 1 0 1 1;
 *      1 0 1 0 0;
 *      0 1 0 1 0;
 *      0 1 0 0 1]
 */
START_TEST (test_fcc03) {
  int xadj[] = {1, 2, 4, 5, 6, 7};
  int adncy[] = {3, 4, 5, 1, 2, 2};
  CC *comps = fcc(5, xadj, adncy);
  int comp[] = {1, 2, 1, 2, 2};
  int begin[] = {1, 3, 6};
  int node[] = {1, 3, 2, 4, 5};
  fail_unless(comps->N == 2,
      "`N` not set correctly on creation");
  fail_unless(memcmp(comps->begin, begin, sizeof(begin)) == 0,
      "`begin` not set correctly on creation");
  fail_unless(memcmp(comps->node, node, sizeof(node)) == 0,
      "`node` not set correctly on creation");
  fail_unless(memcmp(comps->comp, comp, sizeof(comp)) == 0,
      "`comp` not set correctly on creation");
  free_CC(comps);
}
END_TEST

Suite *fcc_suite(void) {
  Suite *s = suite_create("fcc");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_fcc00);
  tcase_add_test(tc_core, test_fcc01);
  tcase_add_test(tc_core, test_fcc02);
  tcase_add_test(tc_core, test_fcc03);
  suite_add_tcase(s, tc_core);

  return s;
}

/* Test for
 *
 * A = [1 0 0 0 0;
 *      0 1 0 0 0;
 *      0 0 1 0 0;
 *      0 0 0 1 0;
 *      0 0 0 0 1]
 */
START_TEST (test_calc_degree00) {
  int *xadj = NULL;
  int *d =calc_degree(5, xadj);
  int degree[] = {0, 0, 0, 0, 0};
  fail_unless(memcmp(degree, d, sizeof(degree)) == 0,
      "`degree` not computed correctly");
  free(d);
  d = NULL;
}
END_TEST

/* Test for
 *
 * A = [1 1 1 1 1;
 *      1 1 0 0 0;
 *      1 0 1 0 0;
 *      1 0 0 1 0;
 *      1 0 0 0 1]
 */
START_TEST (test_calc_degree01) {
  int xadj[] = {1, 5, 6, 7, 8, 9};
  int *d =calc_degree(5, xadj);
  int degree[] = {4, 1, 1, 1, 1};
  fail_unless(memcmp(degree, d, sizeof(degree)) == 0,
      "`degree` not computed correctly");
  free(d);
  d = NULL;
}
END_TEST

/* Test for
 *
 * A = [1 1 1 0 0;
 *      1 1 0 0 0;
 *      1 0 1 0 0;
 *      0 0 0 1 1;
 *      0 0 0 1 1]
 */
START_TEST (test_calc_degree02) {
  int xadj[] = {1, 3, 4, 5, 6, 7};
  int *d =calc_degree(5, xadj);
  int degree[] = {2, 1, 1, 1, 1};
  fail_unless(memcmp(degree, d, sizeof(degree)) == 0,
      "`degree` not computed correctly");
  free(d);
  d = NULL;
}
END_TEST

/* Test for
 *
 * A = [1 0 1 0 0;
 *      0 1 0 1 1;
 *      1 0 1 0 0;
 *      0 1 0 1 0;
 *      0 1 0 0 1]
 */
START_TEST (test_calc_degree03) {
  int xadj[] = {1, 2, 4, 5, 6, 7};
  int *d =calc_degree(5, xadj);
  int degree[] = {1, 2, 1, 1, 1};
  fail_unless(memcmp(degree, d, sizeof(degree)) == 0,
      "`degree` not computed correctly");
  free(d);
  d = NULL;
}
END_TEST

Suite *calc_degree_suite(void) {
  Suite *s = suite_create("calc_degree");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_calc_degree00);
  tcase_add_test(tc_core, test_calc_degree01);
  tcase_add_test(tc_core, test_calc_degree02);
  tcase_add_test(tc_core, test_calc_degree03);
  suite_add_tcase(s, tc_core);

  return s;
}

/* Test for
 *
 * A = [1 0 0 0 0;
 *      0 1 0 0 0;
 *      0 0 1 0 0;
 *      0 0 0 1 0;
 *      0 0 0 0 1]
 */
START_TEST (test_bubble_sort00) {
  int *adjncy = NULL;
  int degree[] = {0, 0, 0, 0, 0};
  bubble_sort(adjncy, degree, 0, 5);
  fail_unless(adjncy == NULL,
      "bubble sort not work correctly");
}
END_TEST

/* Test for
 *
 * A = [1 1 1 1 1;
 *      1 1 0 0 0;
 *      1 0 1 0 0;
 *      1 0 0 1 0;
 *      1 0 0 0 1]
 */
START_TEST (test_bubble_sort01) {
  int adjncy[] = {2, 3, 4, 5, 1, 1, 1, 1};
  int degree[] = {4, 1, 1, 1, 1};
  bubble_sort(adjncy, degree, 0, 4);
  int adjncy_exp[] = {2, 3, 4, 5, 1, 1, 1, 1};
  fail_unless(memcmp(adjncy_exp, adjncy, sizeof(adjncy_exp)) == 0,
      "bubble sort not work correctly");
}
END_TEST

/* Test for
 *
 * A = [1 1 1 0 0;
 *      1 1 0 0 0;
 *      1 0 1 0 0;
 *      0 0 0 1 1;
 *      0 0 0 1 1]
 */
START_TEST (test_bubble_sort02) {
  int adjncy[] = {2, 3, 1, 1, 5, 4};
  int degree[] = {2, 1, 1, 1, 1};
  bubble_sort(adjncy, degree, 0, 2);
  int adjncy_exp[] = {2, 3, 1, 1, 5, 4};
  fail_unless(memcmp(adjncy_exp, adjncy, sizeof(adjncy_exp)) == 0,
      "bubble sort not work correctly");
}
END_TEST

/* Test for
 *
 * A = [1 0 1 0 0;
 *      0 1 0 1 1;
 *      1 0 1 0 0;
 *      0 1 0 1 0;
 *      0 1 0 0 1]
 */
START_TEST (test_bubble_sort03) {
  int adjncy[] = {3, 4, 5, 1, 2, 2};
  int degree[] = {1, 2, 1, 1, 1};
  bubble_sort(adjncy, degree, 0, 1);
  int adjncy_exp[] = {3, 4, 5, 1, 2, 2};
  fail_unless(memcmp(adjncy_exp, adjncy, sizeof(adjncy_exp)) == 0,
      "bubble sort not work correctly");
}
END_TEST

/* Test for
 *
 * A = [1 1 1 1 1;
 *      1 1 1 0 0;
 *      1 1 1 0 1;
 *      1 0 0 1 0;
 *      1 0 1 0 1]
 */
START_TEST (test_bubble_sort04) {
  int adjncy[] = {2, 3, 4, 5, 1, 3, 1, 2, 5, 1, 1, 3};
  int degree[] = {4, 2, 3, 1, 2};
  bubble_sort(adjncy, degree, 0, 4);
  int adjncy_exp[] = {4, 2, 5, 3, 1, 3, 1, 2, 5, 1, 1, 3};
  fail_unless(memcmp(adjncy_exp, adjncy, sizeof(adjncy_exp)) == 0,
      "bubble sort not work correctly");
}
END_TEST

Suite *bubble_sort_suite(void) {
  Suite *s = suite_create("bubble_sort");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_bubble_sort00);
  tcase_add_test(tc_core, test_bubble_sort01);
  tcase_add_test(tc_core, test_bubble_sort02);
  tcase_add_test(tc_core, test_bubble_sort03);
  tcase_add_test(tc_core, test_bubble_sort04);
  suite_add_tcase(s, tc_core);

  return s;
}

int main(void) {
  int number_failed;
  SRunner *sr = srunner_create(CC_suite());
  srunner_add_suite(sr, fcc_suite());
  srunner_add_suite(sr, calc_degree_suite());
  srunner_run_all(sr, CK_NORMAL);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
