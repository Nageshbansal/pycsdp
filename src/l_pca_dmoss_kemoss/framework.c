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

#include "framework.h"

gsl_vector* framework(const gsl_matrix* objectives, int nobj, int popsize, int k_neighbours, double theta) {

    int Nv;
    gsl_matrix* original_correlation;
    gsl_vector* eval;
    gsl_matrix* V;
    size_t* imp_index;
    gsl_vector* percent_eig;
    gsl_matrix* objectives_picked;
    gsl_vector* fe;
    gsl_matrix* correlated_rcm;
    double tcor;
    double m2sigma;
    gsl_matrix* correlated_threshold_rcm;
    gsl_vector* score_obj;
    int m2;
    gsl_vector* fs;
    gsl_vector* norm_score_obj;
    double error;
    size_t* sci;
    gsl_vector* errors;
    gsl_vector* nerrors;
    size_t* serror;
    gsl_vector* cerror;
    int* deltamoss;
    double* kemoss;
    gsl_matrix* dmoss_obj_counter;
    gsl_matrix* kemoss_obj_counter;

    /* 1) Compute Correlation Matrix (R) */
    original_correlation = compute_correlation(objectives, nobj, popsize);

    /* 2) Compute Eigenvalues and Eigenvectors */
    eval=gsl_vector_alloc(nobj);   /* eigenvalues */
    V=gsl_matrix_alloc(nobj,nobj); /* eigenvectors */
    compute_eigen_vectors_values(original_correlation, V, eval, nobj);

    /* 3) Compute sorted eigenvalues index */
    imp_index = index_sort_eigen(eval, nobj);

    /* 4) Compute percentage (ei) for each eigenvalue */
    percent_eig = percentage_eig(eval, nobj);

    /* 5) Compute the number of principal components
     * necessary to meet the variance threshold theta */
    Nv = variance_threshold(percent_eig, imp_index, theta, nobj);

    /* 6) Objectives selected from each principal component */
    objectives_picked = objs_per_pca(V, nobj);
    fe=select_obs(objectives_picked,imp_index,Nv,nobj);

    /* 7) Objectives Correlated based just on their signs */
    correlated_rcm = rcm_analysis_signs(original_correlation, fe, nobj);

    /* 8) Correlation Threshold */
    m2sigma=0.954;
    m2=variance_threshold(percent_eig, imp_index, m2sigma, nobj);
    tcor=compute_tcor(percent_eig, imp_index, m2, nobj);

    /* 9) Objectives Correlated based on signs and correlation threshold */
    correlated_threshold_rcm = rcm_analysis_signs_threshold(original_correlation, fe,tcor, nobj);

    /* 10) Objective selection score */
    score_obj=objective_selection_score(percent_eig, V, imp_index, Nv, nobj);

    /* 11) Final Reduction based on RCM */
    fs=final_reduction_rcm(correlated_threshold_rcm, score_obj, fe, nobj);

    /* 12) Normalised objective selection score and sorted */
    norm_score_obj=norm_objective_selection_score(percent_eig,imp_index,V,nobj);
    sci=(size_t*) malloc(nobj*sizeof(size_t));
    gsl_sort_vector_largest_index(sci,nobj,norm_score_obj);

    /* 13) Compute error incurred in the reduction */
    error=compute_error(norm_score_obj, original_correlation, correlated_threshold_rcm, fs, nobj);

    /* 14) Compute error associated with the removal of each objective */
    errors=compute_errors(fs, norm_score_obj, original_correlation, correlated_threshold_rcm, nobj);

    /* 15) Sorted error associated with the removal of each objective */
    serror=(size_t*) malloc(nobj*sizeof(size_t));
    gsl_sort_vector_smallest_index(serror,nobj,errors);

    /* 16) Normalised error associated with the removal of each objective */
    nerrors=norm_errors(errors,nobj);

    /* 17) Compute error cumulative */
    cerror=cumulative_error(nerrors,serror,nobj);

    /* 18) delta-MOSS */
    deltamoss=get_deltamoss(fs,cerror,nobj);
    dmoss_obj_counter=get_dmoss_obj_counter(fs,cerror,serror,nobj);

    /* 19) k-EMOSS */
    kemoss=get_kemoss(cerror,nobj);
    kemoss_obj_counter=get_kemoss_obj_counter(serror,nobj);

    /* Print all the computed data */
    print_data(original_correlation, eval, V, imp_index, percent_eig, Nv, theta, objectives_picked, fe, nobj, correlated_rcm, tcor, correlated_threshold_rcm, score_obj, fs, norm_score_obj, error, m2, m2sigma, sci, errors, nerrors, serror, cerror, deltamoss, dmoss_obj_counter, kemoss, kemoss_obj_counter);

    gsl_matrix_free(dmoss_obj_counter);
    gsl_matrix_free(kemoss_obj_counter);
    free(serror);
    free(kemoss);
    free(deltamoss);
    gsl_vector_free(cerror);
    gsl_vector_free(errors);
    gsl_vector_free(nerrors);
    free(sci);
    gsl_vector_free(norm_score_obj);
    /*gsl_vector_free(fs);*/
    gsl_vector_free(score_obj);
    gsl_matrix_free(correlated_rcm);
    gsl_matrix_free(correlated_threshold_rcm);
    gsl_vector_free(fe);
    gsl_matrix_free(objectives_picked);
    free(imp_index);
    gsl_vector_free(percent_eig);
    gsl_vector_free(eval);
    gsl_matrix_free(V);
    gsl_matrix_free(original_correlation);

    return fs;
}
