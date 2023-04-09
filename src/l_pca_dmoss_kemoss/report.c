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

void print_data(const gsl_matrix* original_correlation,
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
               ) {

    int i,j;
    FILE *report;
    report = fopen("data_lpca.log","w");

    fprintf(report,"# Dimensionality Reduction using Maximum Variance Unfolding (MVU)\n");

    /* Original Correlation */
    fprintf(report,"\n# Original correlation matrix with values\n");
    for(i=0;i<nobj;i++) {
        for(j=0;j<nobj;j++) {
            fprintf(report,"%f\t",gsl_matrix_get(original_correlation,i,j));
        }
        fprintf(report,"\n");
    }
    fflush(report);
/*
    fprintf(report,"R14 = %f\n",gsl_matrix_get(original_correlation,0,3));
    fprintf(report,"R24 = %f\n",gsl_matrix_get(original_correlation,1,3));
    fprintf(report,"R34 = %f\n",gsl_matrix_get(original_correlation,2,3));
*/

    /* Original Correlation with only signs */
    fprintf(report,"\n# Original correlation matrix with only signs\n");
    for(i=0;i<nobj;i++) {
        for(j=0;j<nobj;j++) {
            fprintf(report,"%d ",i+1);
            if(gsl_matrix_get(original_correlation,i,j)>=0.0) {
                fprintf(report,"+\t");
            } else {
                fprintf(report,"-\t");
            }
        }
        fprintf(report,"\n");
    }
    fflush(report);

    /* Sorted Eigenvalues */
    /*
    fprintf(report,"\n# Sorted Eigenvalues\n");
    for(i=0;i<nobj;i++) {
        fprintf(report,"e%d %f\n",(int)imp_index[i]+1,gsl_vector_get(eval,imp_index[i]));
}
    fflush(report);
    */

    /* Sorted Normalised Eigenvalues */
    fprintf(report,"\n# Sorted Normalised Eigenvalues (ei) of Kernel matrix\n");
    for(i=0;i<nobj;i++) {
        fprintf(report,"e%d %f\n",i+1,gsl_vector_get(percent_eig,imp_index[i]));
    }
    fflush(report);

    /* Number of Principal Components that account for defined variance */
    fprintf(report,"\n# Number of Principal Components (Nv) that account for (theta)\n");
    fprintf(report,"Nv=%d\n",Nv);
    fprintf(report,"theta=%f\n",theta);
    fflush(report);

    /* Eigenvectors */
    /*
    fprintf(report,"\n# Eigenvectors\n");
    for(i=0;i<nobj;i++) {
        for(j=0;j<nobj;j++) {
            fprintf(report,"%f\t",gsl_matrix_get(V,i,j));
        }
        fprintf(report,"\n");
}
    fflush(report);
    */

    /* Sorted Eigenvectors */
    fprintf(report,"\n# Sorted Eigenvectors (Vj) of Kernel matrix\n");
    for(i=0;i<nobj;i++) {
        fprintf(report,"V%d(PCA%d)\t",i+1,i+1);
    }
    fprintf(report,"\n");
    for(i=0;i<nobj;i++) {
        for(j=0;j<nobj;j++) {
            fprintf(report,"%f\t",gsl_matrix_get(V,i,imp_index[j]));
        }
        fprintf(report,"\n");
    }
    fflush(report);

    /* Sorted Objectives picked by each PCA */
    fprintf(report,"\n# Objectives picked by each sorted PCA denoted by 1 else by 0\n");
    for(i=0;i<nobj;i++) {
        fprintf(report,"V%d(PCA%d)\t",i+1,i+1);
    }
    fprintf(report,"\n");
    for(i=0;i<nobj;i++) {
        for(j=0;j<nobj;j++) {
            fprintf(report,"f%d %d\t\t",i+1,(int)gsl_matrix_get(objectives_picked,i,imp_index[j]));
        }
        fprintf(report,"\n");
    }
    fflush(report);

    /* Objectives reduced by Eigenvalue Analysis (Fe) */
    fprintf(report,"\n# Objectives set after reduction by Eigenvalue Analysis\nFe = { ");
    for(i=0;i<nobj;i++) {
        if(gsl_vector_get(fe,i)==1.0)
            fprintf(report,"f%d ",i+1);
    }
    fprintf(report,"}\n");
    fflush(report);

    /* Objectives Correlated based just on their signs */
    fprintf(report,"\n# Matrix of objectives correlated based just on their signs (Equation 3.1)\n");
    for(i=0;i<nobj;i++) {
        fprintf(report,"  %d\t",i+1);
    }
    fprintf(report,"\n");
    for(i=0;i<nobj;i++) {
        for(j=0;j<nobj;j++) {
            fprintf(report,"%d %d\t",i+1,(int)gsl_matrix_get(correlated_rcm,i,j));
        }
        fprintf(report,"\n");
    }
    fflush(report);

    /* Correlation Threshold */
    fprintf(report,"\n# Computation of Tcor (Equation 4)\n"
    		"Tcor = 1.0-e1(1.0-M2sigma/M) = 1.0-%f(1.0-%d/%d) = %f",
    		gsl_vector_get(percent_eig,imp_index[0]),m2,nobj,tcor);
    fprintf(report,"\n");
    fflush(report);

    /* Objectives Correlated based on signs and correlation threshold */
    /*
    fprintf(report,"\n# Matrix of objectives correlated based on signs and correlation threshold\n");
    for(i=0;i<nobj;i++) {
        fprintf(report,"  %d\t",i+1);
    }
    fprintf(report,"\n");
    for(i=0;i<nobj;i++) {
        for(j=0;j<nobj;j++) {
            fprintf(report,"%d %d\t",i+1,(int)gsl_matrix_get(correlated_threshold_rcm,i,j));
        }
        fprintf(report,"\n");
    }
    fflush(report);
    */

    /* Without Threshold Based:Correlation within critical-objectives based Original correlation Matrix */
    fprintf(report,"\n# Without Threshold: Identically Correlated set of objectives in Fe\n");
    for(i=0;i<nobj;i++) {
    	fprintf(report,"S%d={ f%d ",i+1,i+1);
        for(j=0;j<nobj;j++) {
        	if( (gsl_matrix_get(correlated_rcm,i,j)==1.0) ) {
        		fprintf(report,"f%d ",j+1);
        	}
        }
        fprintf(report,"}\n");
    }
    fflush(report);

    /* Threshold Based:Correlation within critical-objectives based Original correlation Matrix */
    fprintf(report,"\n# Threshold Based: Identically Correlated set of objectives in Fe\n");
    for(i=0;i<nobj;i++) {
    	fprintf(report,"S%d={ f%d ",i+1,i+1);
        for(j=0;j<nobj;j++) {
        	if( (gsl_matrix_get(correlated_threshold_rcm,i,j)==1.0) ) {
        		fprintf(report,"f%d ",j+1);
        	}
        }
        fprintf(report,"}\n");
    }
    fflush(report);

    /* Objective selection score (ci) */
    fprintf(report,"\n# Objective selection score (ci=sum(ei*|fij|))\n");
    for(i=0;i<nobj;i++) {
        fprintf(report,"f%d %f\n",i+1,gsl_vector_get(score_obj,i));
    }
    fflush(report);

    /* Normalised objective selection score (ci) */
    /*
    fprintf(report,"\n# Sorted variance accounted by each objective over all Principal Components (ciM=sum(ei*fij^2))\n");
    for(i=0;i<nobj;i++) {
        fprintf(report,"f%d %f\n",(int)sci[i]+1,gsl_vector_get(c,sci[i]));
    }
    fflush(report);
    */

    /* Normalised objective selection score (ci) */
    fprintf(report,"\n# Variance accounted by each objective over all Principal Components (ciM=sum(ei*fij^2))\n");
    for(i=0;i<nobj;i++) {
        fprintf(report,"f%d %f\n",i+1,gsl_vector_get(c,i));
    }
    fflush(report);

    /* Sorted normalised objective selection score (ci) */
    fprintf(report,"\n# Sorted variance accounted by each objective over all Principal Components (ciM=sum(ei*fij^2))\n");
    for(i=0;i<nobj;i++) {
        fprintf(report,"f%d %f\n",(int)sci[i]+1,gsl_vector_get(c,(int)sci[i]));
    }
    fflush(report);

    /* Final Set */
    fprintf(report,"\n# Final Set\n(Fs) = ");
    j=0;
    for(i=0;i<nobj;i++) {
        if(gsl_vector_get(fs,i)==1.0) {
            fprintf(report,"%d ",i+1);
            j++;
        }
    }
    fprintf(report,"\n");
    fprintf(report,"Size = %d\n",j);
    fflush(report);

    fprintf(report,"\n########################\n");
    fprintf(report,"# delta-MOSS and k-EMOSS\n");
    fprintf(report,"########################\n");
    fflush(report);

    /* Error */
    fprintf(report,"\n# Error incurred in this reduction (Et)\n");
    fprintf(report,"Et = %e\n",error);
    fflush(report);

    /* Error per objective */
    fprintf(report,"\n# Error incurred per objective\n");
    for(i=0;i<nobj;i++) {
    	fprintf(report,"f%d %e\n",i+1,gsl_vector_get(errors,i));
    }
    fflush(report);

    /* Sorted error per objective */
    fprintf(report,"\n# Sorted error incurred per objective\n");
    for(i=0;i<nobj;i++) {
    	fprintf(report,"f%d %e\n",(int)serror[i]+1,gsl_vector_get(errors,(int)serror[i]));
    }
    fflush(report);

    /* Normalised sorted error per objective */
    fprintf(report,"\n# Normalised sorted error incurred per objective\n");
    for(i=0;i<nobj;i++) {
    	fprintf(report,"sf%d %e\n",(int)serror[i]+1,gsl_vector_get(nerrors,(int)serror[i]));
    }
    fflush(report);

    /* Error cumulative */
    fprintf(report,"\n# Error cumulative\n");
    for(i=0;i<nobj;i++) {
    	fprintf(report,"f%d %e\n",(int)serror[i]+1, gsl_vector_get(cerror,i));
    }
    fflush(report);

    /* delta-MOSS */
    fprintf(report,"\n# delta-MOSS\n");
    for(i=0;i<10;i++) {
    	fprintf(report,"delta-MOSS(%02d) = %d\n",i*10,deltamoss[i]);
    }
    fprintf(report,"\n");
    for(i=0;i<10;i++) {
    	fprintf(report,"dmoss_objselected(%02d):",i*10);
    	for(j=0;j<nobj;j++) {
    		fprintf(report," %d",(int)gsl_matrix_get(dmoss_obj_counter,i,j));
    	}
    	fprintf(report,"\n");
    }
    fflush(report);

    /* k-EMOSS */
    fprintf(report,"\n# k-EMOSS\n");
    for(i=0;i<nobj-1;i++) {
    	fprintf(report,"k-EMOSS(%02d) = %f\n",i+1, kemoss[i]);
    }
    fprintf(report,"k-EMOSS(%02d) = %f\n",nobj, kemoss[nobj-1]);
    fprintf(report,"\n");
    for(i=0;i<nobj;i++) {
    	fprintf(report,"kemoss_objselected(%02d):",i+1);
    	for(j=0;j<nobj;j++) {
    		fprintf(report," %d",(int)gsl_matrix_get(kemoss_obj_counter,i,j));
    	}
    	fprintf(report,"\n");
    }
    fflush(report);

    fclose(report);

    return;
}

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
                      ) {

    int i,j;
    double value;

    int n_constraints, n_variables;
    int isometry_constraints, centering_constraint;

    n_variables = nobj * nobj;
    isometry_constraints = nobj * k_neighbours;
    centering_constraint = 1;
    n_constraints = isometry_constraints + centering_constraint;

    printf("\nPopulation\n");
    for (i = 0; i < popsize; i++) {
        for (j = 0; j < nobj; j++) {
            value = gsl_matrix_get(pop, i, j);
            printf("%f ", value);
        }
        printf("\n");
    }

    printf("\nTranspose:\n");
    for (i = 0; i < nobj; i++) {
        for (j = 0; j < popsize; j++) {
            value = gsl_matrix_get(popt, i, j);
            printf("%f ", value);
        }
        printf("\n");
    }

    printf("\nCentering:\n");
    for (i = 0; i < nobj; i++) {
        for (j = 0; j < popsize; j++) {
            value = gsl_matrix_get(poptc, i, j);
            printf("%f ", value);
        }
        printf("\n");
    }

    printf("\nGramian Matrix\n");
    for (i = 0; i < nobj; i++) {
        for (j = 0; j < nobj; j++) {
            value = gsl_matrix_get(G, i, j);
            printf("%f ", value);
        }
        printf("\n");
    }
    printf("\nk-nearest Neighbours\n");
    for (i = 0; i < (nobj * k_neighbours); i++) {
        for (j = 0; j < 2; j++) {
            value = gsl_matrix_get(a, i, j);
            printf("%f ", value);
        }
        printf("\n");
    }

    printf("\nConstraints (b)\n");
    for (i = 0; i < n_constraints; i++) {
        value = gsl_vector_get(b, i);
        printf("%f ", value);
    }
    printf("\n");

    printf("\nConstraints (At)\n");
    for (i = 0; i < n_constraints; i++) {
        for (j = 0; j < n_variables; j++) {
            value = gsl_matrix_get(At, i, j);
            printf("%d\t", (int) value);
        }
        printf("\n");
    }

    printf("\nCoefficients (c)\n");
    for (i = 0; i < n_variables; i++) {
        value = gsl_vector_get(c, i);
        printf("%f ", value);
    }
    printf("\n");

    printf("\nLearnt Kernel Matrix:\n");
    for(i=0;i<nobj;i++) {
        for(j=0;j<nobj;j++) {
            printf("%.3f ",gsl_matrix_get(K,i,j));
        }
        printf("\n");
    }
    printf("\n");

    return;
}
