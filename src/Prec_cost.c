/*!\page LICENSE LICENSE

Copyright (C) 2003 by the Board of Trustees of Massachusetts Institute of
Technology, hereafter designated as the Copyright Owners.

License to use, copy, modify, sell and/or distribute this software and
its documentation for any purpose is hereby granted without royalty,
subject to the following terms and conditions:

1.  The above copyright notice and this permission notice must
appear in all copies of the software and related documentation.

2.  The names of the Copyright Owners may not be used in advertising or
publicity pertaining to distribution of the software without the specific,
prior written permission of the Copyright Owners.

3.  THE SOFTWARE IS PROVIDED "AS-IS" AND THE COPYRIGHT OWNERS MAKE NO
REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED, BY WAY OF EXAMPLE, BUT NOT
LIMITATION.  THE COPYRIGHT OWNERS MAKE NO REPRESENTATIONS OR WARRANTIES OF
MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT THE USE OF THE
SOFTWARE WILL NOT INFRINGE ANY PATENTS, COPYRIGHTS TRADEMARKS OR OTHER
RIGHTS. THE COPYRIGHT OWNERS SHALL NOT BE LIABLE FOR ANY LIABILITY OR DAMAGES
WITH RESPECT TO ANY CLAIM BY LICENSEE OR ANY THIRD PARTY ON ACCOUNT OF, OR
ARISING FROM THE LICENSE, OR ANY SUBLICENSE OR USE OF THE SOFTWARE OR ANY
SERVICE OR SUPPORT.

LICENSEE shall indemnify, hold harmless and defend the Copyright Owners and
their trustees, officers, employees, students and agents against any and all
claims arising out of the exercise of any rights under this Agreement,
including, without limiting the generality of the foregoing, against any
damages, losses or liabilities whatsoever with respect to death or injury to
person or damage to property arising from or out of the possession, use, or
operation of Software or Licensed Program(s) by LICENSEE or its customers.

*/



/* these are functions used in mulSetup.c in determining whether to go down */
/* further in partitioning levels */

#include "induct.h"

/* SRW */
double OneCubeCost(cube*****, int, int, int, int, int, double*);
double ratio_of_divided_segs(double, charge*, SYS*);
int is_gp_charge(charge*);
void add_to_counts(cube*, int, int*****, int*****);
/*
 void dump_evalcnts(ssystem*);
 void initCounters(ssystem*);
 int *****make_ints_for_cubes(ssystem*);
 */

/* this function estimates the size of the matrix which will be inverted
 for this cube for the preconditioner.  Since M has not been formed, and
 the size of the preconditioner is based on meshes, not filaments/charges,
 it will estimate that there is one mesh per normal filament and one mesh for
 every two ground plane filaments.  Thus the number of meshes in a cube
 is cube->multisize/2.
 */

double OneCubeCost(cube *****cubes, int i, int j, int k, int l, int side,
		double *dir_cost) {
	int m, n, p;
	double total, dir_total, this_size;

	this_size = cubes[i][j][k][l]->upnumeles[0];

	total = dir_total = 0;
	for (m = j - 2; m <= j + 2; m++)
		for (n = k - 2; n <= k + 2; n++)
			for (p = l - 2; p <= l + 2; p++)
				if ((m >= 0) && (n >= 0) && (p >= 0) && (m < side) && (n < side)
						&& (p < side) && (cubes[i][m][n][p] != NULL)) {
					dir_total += cubes[i][m][n][p]->upnumeles[0];
					if (abs(m - j) < 2 && abs(n - k) < 2 && abs(p - l) < 2)
						total += cubes[i][m][n][p]->multisize / 2.0;
				}

	*dir_cost += this_size * dir_total;

	return total * total * total;
}

double ratio_of_divided_segs(double length, charge *charges, SYS *indsys) {
	SEGMENT *seg;
	int totalfils = 0, broken = 0;
	double rat;

	if (indsys->opts->auto_refine == OFF)
		return 0.0;

	for (seg = indsys->segment; seg != NULL; seg = seg->next) {
		totalfils += seg->num_fils;
		if (seg->length > length) {
			broken += seg->num_fils;
		}
	}

	rat = (double) broken / (double) totalfils;
	if (indsys->opts->debug == ON)
		fprintf(stdout, "To be broken ratio: %lg\n", rat);
	return rat;
}

int is_gp_charge(charge *chg) {
	if (chg->fil->segm->node[0]->gp == NULL)
		return FALSE;
	else
		return TRUE;
}

void add_to_counts(cube *nc, int cols, int *****evals, int *****cnts) {
	cube *na;

	for (na = nc; na != NULL; na = na->parent)
		cnts[na->level][na->j][na->k][na->l] += cols * nc->upnumeles[0];

	evals[nc->level][nc->j][nc->k][nc->l] += 1;
	na = nc->parent;
	if (na != NULL)
		evals[na->level][na->j][na->k][na->l] += 1;
}

#if 1 == 0
void dump_evalcnts(ssystem *sys)
{
	cube *****cubes = sys->cubes;
	int i,j,k,m,side;
	cube *nc, *na;

	printf("      cube    parent          Q2P              L2P               M2P\n");
	printf(" lvl  j,k,l    j,k,l   cubes size_mats    cubes size_mats   cubes size_mats\n");

	for(side = 1, i=0; i <= sys->depth; side *= 2, i++) {
		for(j=0; j < side; j++) {
			for(k=0; k < side; k++) {
				for(m = 0; m < side; m++) {
					nc = cubes[i][j][k][m];
					if (nc != NULL) {
						na = nc->parent;
						if (na == NULL)
						na = nc;
						printf("%3i  %3i %3i %3i   %3i %3i %3i  %3d %10d  %3d %10d  %3d %10d\n",
								i, nc->j, nc->k, nc->l, na->j, na->k, na->l,
								sys->evalQ2Ps[i][j][k][m], sys->cntQ2Ps[i][j][k][m],
								sys->evalL2Ps[i][j][k][m], sys->cntL2Ps[i][j][k][m],
								sys->evalM2Ps[i][j][k][m], sys->cntM2Ps[i][j][k][m]);
					}
				}
			}
		}
	}
}

void initCounters(ssystem *sys)
{

	sys->evalQ2Ps = make_ints_for_cubes(sys);
	sys->evalL2Ps = make_ints_for_cubes(sys);
	sys->evalM2Ps = make_ints_for_cubes(sys);
	sys->cntQ2Ps = make_ints_for_cubes(sys);
	sys->cntL2Ps = make_ints_for_cubes(sys);
	sys->cntM2Ps = make_ints_for_cubes(sys);
}

int *****make_ints_for_cubes(ssystem *sys)
{
	int *****cubes;
	int i,j,k,m,side;

	CALLOC(cubes, sys->depth+1, int****, ON, AMSC);

	/* allocate for levels 0, 1, and 2 (always used) */
	for(side = 1, i=0; i <= sys->depth; side *= 2, i++) {
		CALLOC(cubes[i], side, int***, ON, AMSC);
		for(j=0; j < side; j++) {
			CALLOC(cubes[i][j], side, int**, ON, AMSC);
			for(k=0; k < side; k++) {
				CALLOC(cubes[i][j][k], side, int*, ON, AMSC);
				for(m = 0; m < side; m++)
				cubes[i][j][k][m] = 0;
			}
		}
	}

	return cubes;
}
#endif
