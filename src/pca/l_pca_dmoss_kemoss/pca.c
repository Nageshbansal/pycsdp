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

gsl_matrix* compute_correlation(const gsl_matrix* objectives, int M, int popsize) {

    int i, j, k;

    double m1;
    double var1;

    double* variance1;
    int cnt;

    /*************************/
    int *corindex;
    double** corindex_pca;
    /*************************/
    size_t *imp_index;
    size_t *vrealsort_index;
    size_t *vabssort_index;

    int* S;
    int* RS;
    int* RS_final;
    double* abs_E_mat;
    double** objmat;
    gsl_matrix* original_correlation;

    original_correlation = gsl_matrix_alloc(M,M);
    objmat = (double **)malloc(M*sizeof(double *));
    for(i=0;i<M;i++)
        objmat[i]=(double*)malloc(popsize*sizeof(double));

    variance1 = (double*) malloc(M * sizeof(double));
    corindex = (int *) calloc(M * M, sizeof(int));
    imp_index = (size_t *) calloc(M, sizeof(size_t));
    vrealsort_index = (size_t *) calloc(M, sizeof(size_t));
    vabssort_index = (size_t *) calloc(M, sizeof(size_t));

    S = (int*) calloc(M, sizeof(int));
    RS = (int*) calloc(M, sizeof(int));
    RS_final = (int*) calloc(M, sizeof(int));
    abs_E_mat = (double*) calloc(M, sizeof(double));

    corindex_pca = (double **) calloc(M, sizeof(double*));
    for (i = 0; i < M; i++)
        corindex_pca[i] = (double*) calloc(M, sizeof(double));

    cnt = 0;
    for (i = 0; i < popsize; i++) {
        for (j = 0; j < M; j++)
            objmat[j][cnt] = gsl_matrix_get(objectives,i,j);
        cnt++;
    }

    for (j = 0; j < M; j++) {
        m1 = 0.0;
        for (i = 0; i < cnt; i++)
            m1 += objmat[j][i];
        m1 = m1 / cnt;

        var1 = 0.0;
        for (i = 0; i < cnt; i++) {
            objmat[j][i] -= m1;
            var1 += objmat[j][i] * objmat[j][i];
        }
        variance1[j] = sqrt(var1);
    }

    for (j = 0; j < M; j++) {
        for (k = 0; k < M; k++) {
            m1 = 0.0;
            for (i = 0; i < cnt; i++)
                m1 += objmat[j][i] * objmat[k][i];
            gsl_matrix_set(original_correlation,j,k,m1 / (variance1[j] * variance1[k]));
        }
    }

    for(i=0;i<M;i++)
        free(objmat[i]);
    free(objmat);
    free(variance1);
    free(corindex);
    free(imp_index);
    free(vrealsort_index);
    free(vabssort_index);
    free(S);
    free(RS);
    free(RS_final);
    free(abs_E_mat);
    for (i = 0; i < M; i++)
        free(corindex_pca[i]);
    free(corindex_pca);

    return original_correlation;
}

void compute_eigen_vectors_values(const gsl_matrix* cormat, gsl_matrix* V, gsl_vector* eval, int nobj) {

    gsl_eigen_symmv_workspace* ws;
    gsl_matrix* cormat_temp;

    ws = gsl_eigen_symmv_alloc(4*nobj);
    cormat_temp = gsl_matrix_alloc(nobj,nobj);

    gsl_matrix_memcpy(cormat_temp, cormat);
    gsl_eigen_symmv(cormat_temp,eval,V,ws);

    gsl_matrix_free(cormat_temp);
    gsl_eigen_symmv_free(ws);

    return;
}

size_t* index_sort_eigen(const gsl_vector* eval, int nobj) {

    size_t* imp_index;

    imp_index = (size_t *) malloc(nobj*sizeof(size_t));

    /* getting indices from sorting */
    gsl_sort_vector_largest_index(imp_index,nobj,eval);

    return imp_index;
}

gsl_vector* percentage_eig(const gsl_vector* eval, int nobj) {

    int i;
    gsl_vector* percent_eig;
    double eig_sum;

    percent_eig=gsl_vector_alloc(nobj);

    eig_sum=0;
    for(i=0;i<nobj;i++) {
        eig_sum += gsl_vector_get(eval,i);
    }

    for(i=0;i<nobj;i++) {
        gsl_vector_set(percent_eig,i,gsl_vector_get(eval,i)/eig_sum);
    }

    return percent_eig;
}

int variance_threshold(const gsl_vector* percent_eig, const size_t* imp_index, double theta, int nobj) {

    int i, Nv;
    double value;

    value=0.0;
    for(i=0; i<nobj; i++) {
        value+=gsl_vector_get(percent_eig,imp_index[i]);
        if(value>theta)
            break;
    }

    Nv=i+1;

    return Nv;
}

gsl_matrix* objs_per_pca(const gsl_matrix* V, int nobj) {

    int pcas; /* along the principal components */
    int objs; /* along the objectives */
    int flag;
    int index, index2;
    double dmax;
    double value;
    gsl_matrix* pca_objs;

    pca_objs=gsl_matrix_alloc(nobj,nobj);

    for(pcas=0;pcas<nobj;pcas++) { /* along the principal components */

        /* check if signs are all the same */
        flag=0;
        value = gsl_matrix_get(V,0,pcas);
        for(objs=1;objs<nobj;objs++) { /* along the objectives */
            if(value>=0.0) { /* positive */
                if( gsl_matrix_get(V,objs,pcas)<0.0 ) {
                    flag=1;
                }
            } else { /* negative */
                if(gsl_matrix_get(V,objs,pcas)>=0.0  ) {
                    flag=1;
                }
            }
        }

        /* the elements of this principal component have all the same sign */
        if(flag==0) {

            /* get the first maximum */
            dmax=DBL_MIN;
            index=0;
            for(objs=0;objs<nobj;objs++) { /* along the objectives */
                if( fabs(gsl_matrix_get(V,objs,pcas))>dmax ) {
                    index=objs;
                    dmax=fabs(gsl_matrix_get(V,objs,pcas));
                }
            }

            /* get the second maximum */
            dmax=DBL_MIN;
            index2=0;
            for(objs=0;objs<nobj;objs++) {
                if(objs!=index) {
                    if(fabs(gsl_matrix_get(V,objs,pcas))>dmax) {
                        index2=objs;
                        dmax=fabs(gsl_matrix_get(V,objs,pcas));
                    }
                }
            }

            /* fill the matrix with two previously picked maximums */
            for(objs=0;objs<nobj;objs++) { /* along the objectives */

                if( objs==index ) {
                    gsl_matrix_set(pca_objs,objs,pcas,1.0);
                } else {
                    if( objs==index2 ) {
                        gsl_matrix_set(pca_objs,objs,pcas,1.0);
                    } else {
                        gsl_matrix_set(pca_objs,objs,pcas,0.0);
                    }
                }
            }

        }
        /* the elements of this principal component have different signs */
        else {

            /* get the absolute maximum */
            dmax=DBL_MIN;
            index=0;
            for(objs=0;objs<nobj;objs++) { /* along the objectives */
                if( fabs(gsl_matrix_get(V,objs,pcas))>dmax ) {
                    index=objs;
                    dmax=fabs(gsl_matrix_get(V,objs,pcas));
                }
            }

            for(objs=0;objs<nobj;objs++) { /* along the objectives */

                if( objs==index ) {/* select the max absolute value */
                    gsl_matrix_set(pca_objs,objs,pcas,1.0);
                } else {
                    /* select opposite sign values to maximum
                     * previously picked absolute value */
                    if( gsl_matrix_get(V,index,pcas)>0.0 ) { /* positive */
                        if(gsl_matrix_get(V,objs,pcas)<0.0) { /* pick all negative */
                            gsl_matrix_set(pca_objs,objs,pcas,1.0);
                        } else {
                            gsl_matrix_set(pca_objs,objs,pcas,0.0);
                        }
                    }

                    if( gsl_matrix_get(V,index,pcas)<0.0 ) { /* negative */
                        if(gsl_matrix_get(V,objs,pcas)>0.0) { /* pick all positive */
                            gsl_matrix_set(pca_objs,objs,pcas,1.0);
                        } else {
                            gsl_matrix_set(pca_objs,objs,pcas,0.0);
                        }
                    }
                }
            }
        }
    }

    return pca_objs;
}

gsl_vector* select_obs(const gsl_matrix* objectives_picked, const size_t* imp_index, int Nv, int nobj) {

    int pcas;
    int objs;
    double value;
    gsl_vector* fe;

    fe=gsl_vector_alloc(nobj);
    for(objs=0;objs<nobj;objs++) /* for each objective */
        gsl_vector_set(fe,objs,0.0);

    for(pcas=0;pcas<Nv;pcas++) { /* for each principal component */
        for(objs=0;objs<nobj;objs++) {
            value=gsl_matrix_get(objectives_picked,objs,imp_index[pcas]);
            if(value==1.0) {
                gsl_vector_set(fe,objs,1.0);
            }
        }
    }

    return fe;
}


gsl_matrix* rcm_analysis_signs(const gsl_matrix* original_correlation,const gsl_vector* fe, int nobj) {

    int i,j,k;
    gsl_matrix* rcm_matrix;
    rcm_matrix=gsl_matrix_alloc(nobj,nobj);
    double value1, value2;
    int flag;

    for(i=0;i<nobj;i++) {
        for(j=0;j<nobj;j++) {
            gsl_matrix_set(rcm_matrix,i,j,0.0);
        }
    }

    for(i=0;i<nobj;i++) {
        if(gsl_vector_get(fe,i)==1.0) {
            for(j=0;j<nobj;j++) {
                if(gsl_vector_get(fe,j)==1.0) {
                    if(i!=j) {
                        flag=1;
                        for(k=0;k<nobj;k++) {
                            if(gsl_vector_get(fe,k)==1.0) {
                                value1=gsl_matrix_get(original_correlation,i,k);
                                value2=gsl_matrix_get(original_correlation,j,k);
                                if(value1>=0.0) {/* positive */
                                    if(value2<0.0)
                                        flag=0;
                                } else { /* negative */
                                    if(value2>0.0)
                                        flag=0;
                                }
                            }
                            if(flag==1)
                                gsl_matrix_set(rcm_matrix,i,j,1.0);
                            else
                                gsl_matrix_set(rcm_matrix,i,j,0.0);
                        }
                    }
                }
            }
        }
    }

    return rcm_matrix;
}

double compute_tcor(const gsl_vector* percent_eig, const size_t* imp_index, int m2, int nobj) {

    double value;
    value=1.0-gsl_vector_get(percent_eig,imp_index[0])*(1.0-((double)m2/(double)nobj) );
    return value;
}

gsl_matrix* rcm_analysis_signs_threshold(const gsl_matrix* original_correlation, const gsl_vector* fe, double tcor, int nobj) {
    int i,j,k;
    gsl_matrix* rcm_matrix;
    rcm_matrix=gsl_matrix_alloc(nobj,nobj);
    double value1, value2;
    int flag;

    for(i=0;i<nobj;i++) {
        for(j=0;j<nobj;j++) {
            gsl_matrix_set(rcm_matrix,i,j,0.0);
        }
    }

    for(i=0;i<nobj;i++) {
        if(gsl_vector_get(fe,i)==1.0) {
            for(j=0;j<nobj;j++) {
                if(gsl_vector_get(fe,j)==1.0) {
                    if(i!=j) {
                        flag=1;
                        for(k=0;k<nobj;k++) {
                            if(gsl_vector_get(fe,k)==1.0) {
                                value1=gsl_matrix_get(original_correlation,i,k);
                                value2=gsl_matrix_get(original_correlation,j,k);
                                if(value1>=0.0) {/* positive */
                                    if(value2<0.0)
                                        flag=0;
                                } else { /* negative */
                                    if(value2>0.0)
                                        flag=0;
                                }
                            }
                            if(flag==1 && (gsl_matrix_get(original_correlation,i,j)>=tcor) )
                                gsl_matrix_set(rcm_matrix,i,j,1.0);
                            else
                                gsl_matrix_set(rcm_matrix,i,j,0.0);
                        }
                    }
                }
            }
        }
    }

    return rcm_matrix;
}

gsl_vector* objective_selection_score(const gsl_vector* percent_eig, const gsl_matrix* V, const size_t* imp_index, int Nv, int nobj) {

    int objs,pcas;
    gsl_vector* obj_score;
    double value;

    obj_score=gsl_vector_alloc(nobj);

    for(objs=0;objs<nobj;objs++) {

        value=0.0;
        for(pcas=0;pcas<Nv;pcas++) {
            value += gsl_vector_get(percent_eig,imp_index[pcas])*
                     fabs( gsl_matrix_get(V,objs,imp_index[pcas]));
        }
        gsl_vector_set(obj_score,objs,value);
    }

    return obj_score;
}

gsl_vector* final_reduction_rcm(const gsl_matrix* correlated_threshold_rcm, const gsl_vector* score_obj, const gsl_vector* fe, int nobj) {

    int i,j;
    gsl_vector* fs;
    double dmax;
    int index;

    fs=gsl_vector_alloc(nobj);

    for(i=0;i<nobj;i++)
        gsl_vector_set(fs,i,0.0);

    for(i=0;i<nobj;i++) {
        if(gsl_vector_get(fe,i)==1.0) {
            dmax=gsl_vector_get(score_obj,i);
            index=i;
            for(j=0;j<nobj;j++) {
                if(gsl_vector_get(fe,j)==1.0) {
                    if(i!=j) {
                        if(gsl_matrix_get(correlated_threshold_rcm,i,j)==1.0) {
                            if( gsl_vector_get(score_obj,j)>dmax ) {
                                dmax=gsl_vector_get(score_obj,j);
                                index=j;
                            }
                        }
                    }
                }
            }
            gsl_vector_set(fs,index,1.0);
        }
    }

    return fs;
}

gsl_vector* norm_objective_selection_score(const gsl_vector* percent_eig, const size_t* imp_index, const gsl_matrix* V, int nobj) {

    int i,pcas;
    double value;
    gsl_vector* c;

    c=gsl_vector_alloc(nobj);

    for(i=0;i<nobj;i++) {
        value=0.0;
        for(pcas=0;pcas<nobj;pcas++) {
            value+=gsl_vector_get(percent_eig,imp_index[pcas])*pow(gsl_matrix_get(V,i,imp_index[pcas]),2.0);
        }
        gsl_vector_set(c,i,value);
    }

    return c;
}

double compute_error(const gsl_vector* c, const gsl_matrix* original_correlation, const gsl_matrix* correlated_threshold_rcm, const gsl_vector* fs, int nobj) {

    int i,j;
    double error,value;
    double value1, value2;

    error=0.0;
    for(i=0;i<nobj;i++) {
        if(gsl_vector_get(fs,i)==0.0) { /* only redundant objectives */
            value=DBL_MIN;
            for(j=0;j<nobj;j++) {
                if(gsl_vector_get(fs,j)==1.0) { /* only for important objectives */
                    value1=gsl_matrix_get(correlated_threshold_rcm,i,j);
                    value2=gsl_matrix_get(original_correlation,i,j);
                    if( (value1*value2)>value ) {
                        value=value1*value2;
                    }
                }
            }
            error += gsl_vector_get(c,i)*(1.0-value);
        }
    }

    return error;
}

double compute_error_individual(const gsl_vector* c, const gsl_matrix* original_correlation, const gsl_matrix* correlated_threshold_rcm, const gsl_vector* fs, int nobj, int req_obj) {

    int j;
    double error,value;
    double value1, value2;

    error=0.0;
    value=DBL_MIN;
    for(j=0;j<nobj;j++) {
        if(gsl_vector_get(fs,j)==1.0) { /* only for important objectives */
            value1=gsl_matrix_get(correlated_threshold_rcm,req_obj,j);
            value2=gsl_matrix_get(original_correlation,req_obj,j);
            if( (value1*value2)>value ) {
                value=value1*value2;
            }
        }
    }
    error = gsl_vector_get(c,req_obj)*(1.0-value);

    return error;

}

gsl_vector* compute_errors(const gsl_vector* fs, const gsl_vector* c, const gsl_matrix* original_correlation, const gsl_matrix* correlated_threshold_rcm, int nobj) {

    int i;
    double value;
    gsl_vector* errors;

    errors=gsl_vector_alloc(nobj);

    for(i=0;i<nobj;i++) {
        if( gsl_vector_get(fs,i)==1.0 ) {
            value=gsl_vector_get(c,i);
        } else {
            value=compute_error_individual(c, original_correlation, correlated_threshold_rcm, fs, nobj, i);
        }
        gsl_vector_set(errors,i,value);
    }

    return errors;
}

gsl_vector* norm_errors(const gsl_vector* errors, int nobj) {

    int i;
    double sum, value;
    gsl_vector* nerrors;

    nerrors=gsl_vector_alloc(nobj);

    sum=0.0;
    for(i=0;i<nobj;i++) {
        sum+=gsl_vector_get(errors,i);
    }
    for(i=0;i<nobj;i++) {
        value=gsl_vector_get(errors,i);
        gsl_vector_set(nerrors,i,value/sum);
    }

    return nerrors;
}

gsl_vector* cumulative_error(const gsl_vector* nerrors, const size_t* serror, int nobj) {

    int i;
    double value;
    gsl_vector* cumulative;

    cumulative=gsl_vector_alloc(nobj);

    value=gsl_vector_get(nerrors,serror[0]);
    gsl_vector_set(cumulative,0,value);
    for (i = 0; i < (nobj - 1); i++) {
        value=gsl_vector_get(cumulative,i)+gsl_vector_get(nerrors,serror[i+1]);
        gsl_vector_set(cumulative,i+1,value);
    }

    return cumulative;
}

int myround(double number) {
    return (number >= 0) ? (int)(number + 0.5) : (int)(number - 0.5);
}

int* get_deltamoss(const gsl_vector* fs, const gsl_vector* cerror, int nobj) {

    int i,j;
    int* deltamoss;
    int d,size;
    double error;
    int counter;

    size=10;
    deltamoss = (int*) malloc(size*sizeof(int));

    counter=0;
    for(i=0;i<nobj;i++) {
        if(gsl_vector_get(fs,i)==1.0) {
            counter++;
        }
    }
    deltamoss[0]=counter;	/* 0% delta equals to fs */

    /* for remaining errors 10%-90% */
    error=0.10;
    for(i=0; i<(size-1);i++) {
        d=0;
        for(j=0; j<nobj; j++) {
            if(gsl_vector_get(cerror,j) >= error) {
                d++;
            }
        }
        deltamoss[i+1]=d;
        error += 0.1;
    }

    return deltamoss;
}

gsl_matrix* get_dmoss_obj_counter(const gsl_vector* fs, const gsl_vector* cerror, const size_t* serror, int nobj) {

    int i,j,size;
    double error;
    int* dmoss_obj_counter;
    gsl_matrix* m_dmoss_obj_counter;

    dmoss_obj_counter = (int*) calloc(nobj,sizeof(int));
    size=10;
    m_dmoss_obj_counter=gsl_matrix_alloc(size,nobj);

    /* delta=0% */
    for(i=0;i<nobj;i++) {
        gsl_matrix_set(m_dmoss_obj_counter,0,i,gsl_vector_get(fs,i));
    }

    error=0.10;
    for(i=0; i<(size-1);i++) {
        for(j=0; j<nobj; j++) {
            if(gsl_vector_get(cerror,j) >= error) {
                dmoss_obj_counter[serror[j]]=1;
            }
        }
        for(j=0; j<nobj; j++) {
            gsl_matrix_set(m_dmoss_obj_counter,i+1,j,dmoss_obj_counter[j]);
            dmoss_obj_counter[j]=0;
        }
        error += 0.1;
    }

    free(dmoss_obj_counter);

    return m_dmoss_obj_counter;
}

double* get_kemoss(const gsl_vector* cerror, int nobj) {

    int i;
    double* kemoss;

    kemoss=(double*) malloc(nobj*sizeof(double));

    for(i=0;i<nobj-1;i++) {
        kemoss[i]=gsl_vector_get(cerror,nobj-i-2);
    }
    kemoss[nobj-1]=0.0;

    return kemoss;
}

gsl_matrix* get_kemoss_obj_counter(const size_t* serror, int nobj) {

    int i,j;
    int* kemoss_obj_counter;
    gsl_matrix* mkemoss_obj_counter;

    kemoss_obj_counter=(int*) calloc(nobj,sizeof(int));
    mkemoss_obj_counter = gsl_matrix_alloc(nobj,nobj);

    for(i=0;i<nobj-1;i++) {
        for (j = nobj-i-1; j < nobj; j++) {
            kemoss_obj_counter[serror[j]]=1;
        }
        for(j=0;j<nobj;j++) {
            gsl_matrix_set(mkemoss_obj_counter,i,j,kemoss_obj_counter[j]);
        }
    }
    for(j=0;j<nobj;j++) {
        gsl_matrix_set(mkemoss_obj_counter,nobj-1,j,1.0);
    }

    free(kemoss_obj_counter);

    return mkemoss_obj_counter;
}
