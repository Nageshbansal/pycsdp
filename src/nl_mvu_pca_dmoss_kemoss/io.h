/*************************************************************************

 I/O functions

 ---------------------------------------------------------------------

                       Copyright (c) 2005, 2006
                  Carlos Fonseca <cmfonsec@ualg.pt>
             Manuel Lopez-Ibanez <m.lopez-ibanez@napier.ac.uk>
                    Luis Paquete <lpaquete@ualg.pt>

 This program is free software (software libre); you can redistribute
 it and/or modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2 of the
 License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, you can obtain a copy of the GNU
 General Public License at:
                 http://www.gnu.org/copyleft/gpl.html
 or by writing to:
           Free Software Foundation, Inc., 59 Temple Place,
                 Suite 330, Boston, MA 02111-1307 USA

 ----------------------------------------------------------------------

*************************************************************************/

#ifndef _HV_IO_H_
#define _HV_IO_H_

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>


/* If we're not using GNU C, elide __attribute__ */
#ifndef __GNUC__
#  define  __attribute__(x)  /* NOTHING */
#endif

void 
errprintf(const char * template,...) 
/* enables the compiler to check the format string against the
   parameters */  __attribute__ ((format(printf, 1, 2)));

void warnprintf(const char *template,...)
/* enables the compiler to check the format string against the
   parameters */  __attribute__ ((format(printf, 1, 2)));

int 
read_input(FILE *infile, const char *filename, 
           double **datap, int *ncolsp, int **cumsizesp, int *nrunsp);

#endif
