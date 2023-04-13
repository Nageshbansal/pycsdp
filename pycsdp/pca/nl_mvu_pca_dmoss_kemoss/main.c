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

#include "global.h"

int main(double *data, int ncols, int nrows)
{

    FILE *infile;
    int rows;
    int objset_size, i, j;
    char* read_set=NULL;
    char* pread_set=NULL;
    int* objset;
    gsl_matrix* objectives;
    int nobj = 0, nruns = 0, popsize;
    int *cumsizes = NULL;
    gsl_vector* fs;
    int new_nobj;

    int opt; /* it's actually going to hold a char */
    int longopt_index;

    /* framework parameters */
    int k_neighbours;
    double theta;

    /* see the man page for getopt_long for an explanation of these fields */
    static struct option long_options[] = {
    		{"help",		no_argument,		NULL,	'h'},
    		{"version",		no_argument,		NULL,	'V'},
    		{NULL, 0, NULL, 0} /* marks end of list */
    };

    /*while (0 < (opt = getopt_long(argc, argv, "s:Vvh",
                                  long_options, &longopt_index))) {

        switch (opt) {

        case 's':
        	read_set=(char*)malloc(500*sizeof(char));
        	pread_set=read_set;
        	strcpy(read_set,optarg);
        	break;

        case 'V':
        case 'v':/* --version *//*
            version();
            exit(EXIT_SUCCESS);

        case '?':
            /* getopt prints an error message right here *//*
            fprintf(stderr, "Try `%s --help' for more information.\n", program_invocation_short_name);
            exit(EXIT_FAILURE);
        case 'h':
            usage();
            exit(EXIT_SUCCESS);

        default: /* should never happen *//*
            abort();
        }
    }*/

    /*if (optind < argc)
        for (; optind < argc; optind++)
            if (strncmp(argv[optind],"-", strlen("-")+1)) {
                infile = fopen(argv[optind],"r");
                if (infile == NULL) {
                    fprintf (stderr, "%s: %s: %s\n", program_invocation_short_name, argv[optind], strerror(errno));
                    exit(EXIT_FAILURE);
                }
                read_input(infile, argv[optind], &data, &nobj, &cumsizes, &nruns);
                fclose(infile);
            } else
                read_input(stdin, "stdin", &data, &nobj, &cumsizes, &nruns);
    else
        read_input(stdin, "stdin", &data, &nobj, &cumsizes, &nruns);*/

    nobj = nrows;
    cumsizes = ncols;
    objset=(int*) malloc(nobj*sizeof(int));

    if(read_set) {
    	objset_size=0;
    	while( sscanf(read_set,"%d%n",&i,&j)==1 ) {
    		objset[objset_size]=i-1;
    		read_set+=j;
    		objset_size++;
    	}
    	if(objset_size != nobj) {
    		printf("The number of objectives in the Set must be equal to number of provided objectives\n");
    		exit(-1);
    	}
    }
    else {
    	for(i=0;i<nobj;i++)
    		objset[i]=i;
    }

    popsize=cumsizes;

    objectives = perform_conversion(data, nobj, popsize);

    /**********************/
    /* Call the Framework */
    /**********************/
    k_neighbours=nobj-1;
    theta=0.997;
    fs=framework(objectives, nobj, popsize, k_neighbours, theta);

    new_nobj=0;
    for(i=0;i<nobj;i++) {
    	if( gsl_vector_get(fs,i)==1.0 ) {
    		new_nobj++;
    	}
    }

    printf("(Fs) =");
    for(i=0;i<nobj;i++) {
    	if( gsl_vector_get(fs,i)==1.0 ) {
    		printf(" %d",objset[i]+1);
    	}
    }
    printf("\n");
    printf("Size = %d\n",new_nobj);

    gsl_vector_free(fs);
    if(pread_set)
    	free(pread_set);
    free(objset);
    gsl_matrix_free(objectives);
    free(data);

    exit(EXIT_SUCCESS);
}

gsl_matrix* perform_conversion(const double* data, int nobj, int rows) {

    gsl_matrix* objectives;
    int n,r,k;
    k=0;
    objectives=gsl_matrix_alloc(rows,nobj);
    for(r=0; r<rows; r++) {
        for(n=0;n<nobj; n++, k++)
            gsl_matrix_set(objectives,r,n,data[k]);
    }
    return objectives;
}

void usage(void) {
    printf("\nUsage: %s [OPTIONS] [FILE...]\n\n", program_invocation_short_name);
    printf("This population analyser applied different types of dimensionality\n"
           "reduction techniques to a population of solutions.				  \n"
           "Nonlinear Component Analysis Based on Maximum Variance Unfolding  \n\n");

    printf(	"Options:\n"
            " -h, --help        print this summary and exit.                  \n"
            " --version     	print version number and exit.                \n\n");
    return;
}

void version(void) {

    printf("%s version %s \n\n", program_invocation_short_name, VERSION);
    printf(	"Copyright (C) 2010  \nJoao Antonio Duro (j.a.duro@cranfield.ac.uk) and "
            "\nDhish Saxena (d.saxena@cranfield.ac.uk)\n\n"
            "This is free software, and you are welcome to redistribute it under certain\n"
            "conditions.  See the GNU General Public License for details. There is NO   \n"
            "warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n\n");
}
