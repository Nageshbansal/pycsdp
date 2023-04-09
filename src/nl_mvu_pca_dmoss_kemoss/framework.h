/****************************************************************************
**
** Copyright (C) 2017
** Indian Institute of Technology, Roorkee (www.iitr.ac.in)
** The University of Sheffield (www.sheffield.ac.uk)
**
** This file is part of NL-MVU-PCA.
**
** GNU Lesser General Public License Usage
** This file may be used under the terms of the GNU Lesser General
** Public License version 2.1 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPL included in the
** packaging of this file.  Please review the following information to
** ensure the GNU Lesser General Public License version 2.1 requirements
** will be met: http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html.
**
****************************************************************************/

#ifndef FRAMEWORK_H_
#define FRAMEWORK_H_

#include <math.h>
#include <float.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>

#include "csdp_library/declarations.h"

/* framework.c */
gsl_vector* framework(const gsl_matrix* objectives, int nobj, int popsize, int k_neighbours, double theta);

/* pca.c */
gsl_vector* norm_errors(const gsl_vector* errors, int nobj);
int myround(double number);
gsl_matrix* compute_correlation(const gsl_matrix* objectives, int M, int popsize);
gsl_vector* select_obs(const gsl_matrix* objectives_picked, const size_t* imp_index, int Nv, int nobj);
gsl_matrix* objs_per_pca(const gsl_matrix* V, int nobj);
int variance_threshold(const gsl_vector* percent_eig, const size_t* imp_index, double theta, int nobj);
size_t* index_sort_eigen(const gsl_vector* eval, int nobj);
gsl_vector* percentage_eig(const gsl_vector* eval, int nobj);
void compute_eigen_vectors_values(const gsl_matrix* cormat, gsl_matrix* V, gsl_vector* eval, int nobj);
gsl_matrix* compute_correlation(const gsl_matrix* objectives, int M, int popsize);
gsl_matrix* rcm_analysis_signs(const gsl_matrix* original_correlation,const gsl_vector* fe, int nobj);
double compute_tcor(const gsl_vector* percent_eig, const size_t* imp_index, int m2, int nobj);
gsl_matrix* rcm_analysis_signs_threshold(const gsl_matrix* original_correlation, const gsl_vector* fe, double tcor, int nobj);
gsl_vector* objective_selection_score(const gsl_vector* percent_eig, const gsl_matrix* eval, const size_t* imp_index, int Nv, int nobj);
gsl_vector* final_reduction_rcm(const gsl_matrix* correlated_threshold_rcm, const gsl_vector* score_obj, const gsl_vector* fe, int nobj);
gsl_vector* norm_objective_selection_score(const gsl_vector* percent_eig, const size_t* imp_index, const gsl_matrix* V, int nobj);
double compute_error(const gsl_vector* c, const gsl_matrix* original_correlation, const gsl_matrix* correlated_threshold_rcm, const gsl_vector* fs, int nobj);
double compute_error_individual(const gsl_vector* c, const gsl_matrix* original_correlation, const gsl_matrix* correlated_threshold_rcm, const gsl_vector* fs, int nobj, int req_obj);
gsl_vector* compute_errors(const gsl_vector* fs, const gsl_vector* c, const gsl_matrix* original_correlation, const gsl_matrix* correlated_threshold_rcm, int nobj);
gsl_vector* cumulative_error(const gsl_vector* errors, const size_t* serror, int nobj);
int* get_deltamoss(const gsl_vector* fs, const gsl_vector* cerror, int nobj);
double* get_kemoss(const gsl_vector* cerror, int nobj);
gsl_matrix* get_dmoss_obj_counter(const gsl_vector* fs, const gsl_vector* cerror, const size_t* serror, int nobj);
gsl_matrix* get_kemoss_obj_counter(const size_t* serror, int nobj);

/* kernel.c */
gsl_matrix* kernel_build(const gsl_matrix* objectives, int nobj, int popsize, int k_neighbours);
gsl_matrix* call_csdp(gsl_vector* c, gsl_matrix* At, gsl_vector* b_in, int nobj, int n_variables, int n_constraints);
gsl_matrix* transpose_matrix(const gsl_matrix* A, int rows, int cols);
gsl_matrix* centering_matrix(const gsl_matrix* A, int rows, int cols);
gsl_matrix* gramian_matrix(const gsl_matrix* A, int rows, int cols);
gsl_matrix* k_nearest_neighbour(const gsl_matrix* A, int rows, int cols, const int k);
gsl_vector* constraints_right_hand_side(const gsl_matrix* a, int n_constraints, const gsl_matrix* G, int nobj);
gsl_matrix* constraints_left_hand_side(const gsl_matrix* a, int n_constraints, int n_variables, int nobj);
gsl_vector* function_coefficients(const gsl_matrix* a, int n_variables, int nobj);

/* report.c */
void print_data(const gsl_matrix* original_correlation,
                const gsl_matrix* cormat,
                const gsl_vector* eval,
                const gsl_matrix* V,
                const size_t* imp_index,
                const gsl_vector* percent_eig,
                int Nv,
                double theta,
                const gsl_matrix* objectives_picked,
                const gsl_vector* fe,
                int nobj,
                const gsl_matrix* correlated_rcm,
                double tcor,
                const gsl_matrix* correlated_threshold_rcm,
                const gsl_vector* score_obj,
                const gsl_vector* fs,
                const gsl_vector* c,
                double error,
                int m2,
                double m2sigma,
                const size_t* sci,
                const gsl_vector* errors,
                const gsl_vector* nerrors,
                const size_t* serror,
                const gsl_vector* cerror,
                const int* deltamoss,
                const gsl_matrix* dmoss_obj_counter,
                const double* kemoss,
                const gsl_matrix* kemoss_obj_counter
);

void print_kernel_data(const gsl_matrix* pop,
                       const gsl_matrix* popt,
                       const gsl_matrix* poptc,
                       const gsl_matrix* G,
                       const gsl_matrix* a,
                       const gsl_vector* b,
                       const gsl_matrix* At,
                       const gsl_vector* c,
                       const gsl_matrix* K,
                       int k_neighbours,
                       int popsize,
                       int nobj
);

#endif /* FRAMEWORK_H_ */
