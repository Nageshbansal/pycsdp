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

gsl_matrix* kernel_build(const gsl_matrix* pop, int nobj, int popsize, int k_neighbours) {

    /*int i,j;*/
    int n_constraints, n_variables;
    int isometry_constraints, centering_constraint;

    gsl_matrix* popt; /* (nobj,popsize) */
    gsl_matrix* poptc; /* (nobj,popsize) */
    gsl_matrix* G; /* (nobj,nobj)    */
    gsl_matrix* a; /* (nobj*k_neighbours,2) */
    gsl_vector* b; /* (n_constraints)*/
    gsl_matrix* At; /* (n_constraints,n_variables) */
    gsl_vector* c; /* (n_variables) */
    gsl_matrix* K; /* (nobj,nobj) */

    n_variables = nobj * nobj;
    isometry_constraints = nobj * k_neighbours;
    centering_constraint = 1;
    n_constraints = isometry_constraints + centering_constraint;

    /* 1) Transpose the matrix */
    popt = transpose_matrix(pop, popsize, nobj);

    /* 2) Mean-centre the matrix */
    poptc = centering_matrix(popt, nobj, popsize);

    /* 3) Compute the Gramian matrix (or also called Correlation Matrix) */
    G = gramian_matrix(poptc, nobj, popsize);

    /* 4) Compute k-nearest neighbours */
    a = k_nearest_neighbour(poptc, nobj, popsize, k_neighbours);

    /* 5) Compute the right hand side constraints (b) */
    b = constraints_right_hand_side(a, n_constraints, G, nobj);

    /* 6) Compute the left hand side constraints (At) */
    At = constraints_left_hand_side(a, n_constraints, n_variables, nobj);

    /* 7) Compute coefficients of the function to maximise (Tr(K)) */
    c = function_coefficients(a, n_variables, nobj);

    /* 8) Compute the Kernel Matrix using 'At' 'b' and 'c' */
    K = call_csdp(c, At, b, nobj, n_variables, n_constraints);

    /* Uncomment for debug purposes*/
    /*print_kernel_data(pop, popt, poptc, G, a, b, At, c, K, k_neighbours, popsize, nobj);*/

    gsl_vector_free(c);
    gsl_matrix_free(At);
    gsl_vector_free(b);
    gsl_matrix_free(popt);
    gsl_matrix_free(poptc);
    gsl_matrix_free(G);
    gsl_matrix_free(a);

    return K;
}

gsl_matrix* transpose_matrix(const gsl_matrix* A, int rows, int cols) {

    int i, j;
    double value;
    gsl_matrix* At = gsl_matrix_alloc(cols, rows);

    for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++) {
            value = gsl_matrix_get(A, i, j);
            gsl_matrix_set(At, j, i, value);
        }

    return At;
}

gsl_matrix* centering_matrix(const gsl_matrix* A, int rows, int cols) {

    int i, j;
    double temp, temp2;

    gsl_matrix* Ac = gsl_matrix_alloc(rows, cols);
    gsl_vector* mean = gsl_vector_alloc(rows);

    /*1. build the means */
    for (i = 0; i < rows; i++) {
        temp = 0.0;
        for (j = 0; j < cols; j++) {
            temp += gsl_matrix_get(A, i, j);
        }
        gsl_vector_set(mean, i, temp / cols);
    }

    /*2. subtract the means and build Ac*/
    for (i = 0; i < rows; i++) {
        temp2 = gsl_vector_get(mean, i);
        for (j = 0; j < cols; j++) {
            temp = gsl_matrix_get(A, i, j);
            gsl_matrix_set(Ac, i, j, temp - temp2);
        }
    }

    gsl_vector_free(mean);

    return Ac;
}

gsl_matrix* gramian_matrix(const gsl_matrix* A, int rows, int cols) {

    int i, j, k;
    double temp1, temp2;
    double value;
    gsl_matrix* G = gsl_matrix_alloc(rows, rows);

    for (i = 0; i < rows; i++) {

        for (j = 0; j < rows; j++) {
            value = 0.0;
            for (k = 0; k < cols; k++) {
                temp1 = gsl_matrix_get(A, i, k);
                temp2 = gsl_matrix_get(A, j, k);
                value += temp1 * temp2;
            }
            gsl_matrix_set(G, i, j, value);
        }
    }

    return G;
}

gsl_matrix* k_nearest_neighbour(const gsl_matrix* A, int rows, int cols,
                                const int k) {

    int i, j, l, dd, xh, aa, xi, xjj;
    double d;
    double temp1, temp2;
    gsl_matrix* a;
    gsl_vector* D;
    int* xj;
    gsl_permutation * p;

    a = gsl_matrix_alloc(rows * k, 2);
    D = gsl_vector_alloc(rows - 1);
    p = gsl_permutation_alloc(rows - 1);
    xj = (int*) malloc((rows - 1) * sizeof(int));

    /* Fix one point */
    aa = 0;
    for (i = 0; i < rows; i++) {
        xi = i;
        dd = 0;
        xjj = 0;
        /* Determine the distances */
        for (j = 0; j < rows; j++) {
            if (i != j) {
                xj[xjj] = j;
                xjj++;
                d = 0.0;
                for (l = 0; l < cols; l++) {
                    temp1 = gsl_matrix_get(A, i, l);
                    temp2 = gsl_matrix_get(A, j, l);
                    d += temp1 * temp2;
                    d = fabs(d);
                }
                gsl_vector_set(D, dd, d);
                dd++;
            }
        }
        /* Sort the obtained distances (get only indices) increasing order */
        gsl_sort_vector_index(p, D);

        for (j = 0; j < k; j++) {
            xh = gsl_permutation_get(p, j);
            gsl_matrix_set(a, aa, 0, xi);
            gsl_matrix_set(a, aa, 1, xj[xh]);
            aa++;
        }
    }

    free(xj);
    gsl_vector_free(D);
    gsl_permutation_free(p);

    return a;
}

gsl_vector* constraints_right_hand_side(const gsl_matrix* a, int n_constraints,
                                        const gsl_matrix* G, int nobj) {

    int i, i1, j1;
    double value;

    gsl_vector* b;
    b = gsl_vector_alloc(n_constraints);

    /* Set first constraint (Sum(Kij)==0) */
    gsl_vector_set(b, 0, 0.0);

    /* Set remaining constraints accordingly to G and a*/
    for (i = 1; i < n_constraints; i++) {

        i1 = gsl_matrix_get(a, i - 1, 0);
        j1 = gsl_matrix_get(a, i - 1, 1);

        value = gsl_matrix_get(G, i1, i1) + gsl_matrix_get(G, j1, j1)
                - gsl_matrix_get(G, i1, j1) - gsl_matrix_get(G, j1, i1);

        gsl_vector_set(b, i, value);
    }

    return b;
}

gsl_matrix* constraints_left_hand_side(const gsl_matrix* a, int n_constraints,
                                       int n_variables, int nobj) {

    int i;
    int i1, j1;
    int c1, c2, c3, c4;
    gsl_matrix* At;
    At = gsl_matrix_alloc(n_constraints, n_variables);

    gsl_matrix_set_zero(At); /* Set all elements to zero */

    /* Set first constraint (Sum(Kij)==0) */
    for (i = 0; i < n_variables; i++)
        gsl_matrix_set(At, 0, i, 1.0);

    /* Set remaining constraints according with Nearest-Neighbours */
    for (i = 1; i < n_constraints; i++) {

        i1 = gsl_matrix_get(a, i - 1, 0); /* 1 - 0 */
        j1 = gsl_matrix_get(a, i - 1, 1); /* 2 - 1 */

        c1 = i1 * nobj + i1; /* -> 1 - 0 */
        c2 = j1 * nobj + j1; /* -> 5 - 4 */
        c3 = i1 * nobj + j1; /* -> 2 - 1 */
        c4 = j1 * nobj + i1; /* -> 4 - 3 */

        /*printf("c1 %d c2 %d c3 %d c4 %d\n", c1, c2, c3, c4);*/

        gsl_matrix_set(At, i, c1, 1.0);
        gsl_matrix_set(At, i, c2, 1.0);
        gsl_matrix_set(At, i, c3, -1.0);
        gsl_matrix_set(At, i, c4, -1.0);
    }

    return At;
}

gsl_vector* function_coefficients(const gsl_matrix* a, int n_variables,
                                  int nobj) {

    int i, j, k;
    double value;
    gsl_vector* c;
    c = gsl_vector_alloc(n_variables);

    k = 0;
    for (i = 0; i < nobj; i++) {
        for (j = 0; j < nobj; j++) {
            if (i == j)
                value = -2.0 / nobj;
            else
                value = 1.0 / nobj;
            gsl_vector_set(c, k, value);
            k++;
        }
    }

    return c;
}

gsl_matrix* call_csdp(gsl_vector* c, gsl_matrix* At, gsl_vector* b_in,
                      int nobj, int n_variables, int n_constraints) {

    int i,j,constr,n_entries,counter;
    gsl_matrix* K;

    struct blockmatrix C;
    double* b;
    struct constraintmatrix *constraints;
    struct blockmatrix X, Z;
    double* y;
    double pobj, dobj;
    struct sparseblock* blockptr;
    int ret;

    /***********************************/
    /* Objective Function Coefficients */
    /***********************************/

    C.nblocks = 1;
    C.blocks = (struct blockrec *) malloc(2 * sizeof(struct blockrec));
    if (C.blocks == NULL) {
        printf("Couldn't allocate storage for C!\n");
        exit(1);
    };

    /* fill the block */
    C.blocks[1].blockcategory = MATRIX;
    C.blocks[1].blocksize = nobj;
    C.blocks[1].data.mat = (double *) malloc(nobj * nobj * sizeof(double));
    if (C.blocks[1].data.mat == NULL) {
        printf("Couldn't allocate storage for C!\n");
        exit(1);
    }
    for (i = 0; i < nobj; i++) {
        for (j = 0; j < nobj; j++) {
            C.blocks[1].data.mat[ijtok(i+1,j+1,nobj)] =
                gsl_vector_get(c,(nobj * i) + j);
        }
    }

    /*******************************/
    /* Constraints right hand side */
    /*******************************/

    b = (double*) malloc((n_constraints + 1) * sizeof(double));
    if (b == NULL) {
        printf("Could not allocate storage for b!\n");
        exit(-1);
    }
    for (constr = 1; constr <= n_constraints; constr++) {
        b[constr] = gsl_vector_get(b_in,constr - 1);
    }

    /******************************/
    /* Constraints left hand side */
    /******************************/

    constraints = (struct constraintmatrix*) malloc(
                      (n_constraints + 1) * sizeof(struct constraintmatrix));
    if (constraints == NULL) {
        printf("Failed to allocate storage for constraints");
    }

    for (constr = 1; constr <= n_constraints; constr++) {

        constraints[constr].blocks = NULL;

        blockptr = (struct sparseblock *) malloc(sizeof(struct sparseblock));
        if (blockptr == NULL) {
            printf("Allocation of constraint block failed!\n");
            exit(1);
        }

        /* find out the number of non-zero entries */
        n_entries = 0;
        for (i = 0; i < nobj; i++) {
            for (j = i; j < nobj; j++) {
                if (gsl_matrix_get(At,constr - 1,(nobj * i) + j) != 0.0)
                    n_entries++;
            }
        }

        blockptr->blocknum = 1;
        blockptr->blocksize = nobj;
        blockptr->constraintnum = constr;
        blockptr->next = NULL;
        blockptr->nextbyblock = NULL;
        blockptr->entries = (double *) malloc((n_entries + 1) * sizeof(double));
        if (blockptr->entries == NULL) {
            printf("Allocation of constraint block failed!\n");
            exit(1);
        }
        blockptr->iindices = (int *) malloc((n_entries + 1) * sizeof(int));
        if (blockptr->iindices == NULL) {
            printf("Allocation of constraint block failed!\n");
            exit(1);
        }
        blockptr->jindices = (int *) malloc((n_entries + 1) * sizeof(int));
        if (blockptr->jindices == NULL) {
            printf("Allocation of constraint block failed!\n");
            exit(1);
        }

        blockptr->numentries = n_entries;

        counter = 1;
        for (i = 0; i < nobj; i++) {
            for (j = i; j < nobj; j++) {
                if (gsl_matrix_get(At,constr - 1,nobj * i + j) != 0.0) {
                    blockptr->iindices[counter] = i + 1;
                    blockptr->jindices[counter] = j + 1;
                    blockptr->entries[counter] =
                        gsl_matrix_get(At,constr - 1,(nobj * i) + j);
                    counter++;
                }
            }
        }

        blockptr->next = constraints[constr].blocks;
        constraints[constr].blocks = blockptr;
    }

    /*
    write_prob("prob.dat-s", nobj, n_constraints, C, b, constraints);
    ret = read_prob("prob.dat-s", &nobj, &counter, &C, &b, &constraints, 1);
    if (ret != 0) {
    	printf("Giving up.\n");
    	exit(10);
}
    */

    initsoln(nobj, n_constraints, C, b, constraints, &X, &y, &Z);
    ret = easy_sdp(nobj, n_constraints, C, b, constraints, 0.0, &X, &y, &Z,
                   &pobj, &dobj);
    /*
    if (ret == 0)
    	printf("The objective value is %.7e\n", (dobj + pobj) / 2);
    else
    	printf("SDP failed.\n");
    write_sol("prob.sol", nobj, n_constraints, X, y, Z);*/

    K = gsl_matrix_alloc(nobj, nobj);

    for(i=0; i<nobj; i++) {
        for(j=0;j<nobj; j++) {
            gsl_matrix_set(K,i,j,X.blocks[1].data.mat[ijtok(i+1,j+1,nobj)]);
        }
    }

    free_prob(nobj, n_constraints, C, b, constraints, X, y, Z);

    return K;
}

