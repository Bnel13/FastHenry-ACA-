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


/* This contains functions needed for the new preconditioner */

#include "induct.h"
#include "spMatrix.h"

static char outfname[80];

/* used in shift_mutual */
static double radius_factor;

/* SRW */
void choose_and_setup_precond(SYS*);
void get_selfs(SYS*);
void fill_spPre(ssystem*, SYS*, double);
void create_sparMatrix(SYS*);
void fill_bySegment(ssystem*, SYS*, double);
void fill_diagL(ssystem*, SYS*, double);
void fill_diagR(SYS*);
double shift_mutual(FILAMENT*, FILAMENT*, ssystem*);

void choose_and_setup_precond(SYS *indsys) {
	int precond_choice = indsys->opts->precond;

	indsys->precond_type = precond_choice;
	indsys->precond_subtype = OFF; /* will be changed below if necessary */

	/* default */
	if (precond_choice == ON || precond_choice == SPARSE) {
		indsys->precond_type = SPARSE;
		indsys->precond_subtype = CUBEL;
	} else if (precond_choice == LOC) {
		indsys->precond_type = LOC;
		indsys->precond_subtype = OVERLAP;
	} else if (precond_choice == POSDEF_LOC) {
		indsys->precond_type = LOC;
		indsys->precond_subtype = POSDEF_LOC;
		if (indsys->opts->mat_vect_prod == DIRECT)
			fprintf(stderr,
					"Warning: -pposdef not supported with direct matrix\n \
vector products.  POSDEF option ignored.\n");
	} else if (precond_choice == SEGML) {
		indsys->precond_type = SPARSE;
		indsys->precond_subtype = SEGML;
	} else if (precond_choice == DIAGL) {
		indsys->precond_type = SPARSE;
		indsys->precond_subtype = DIAGL;
	} else if (precond_choice == SHELLS) {
		indsys->precond_type = SPARSE;
		indsys->precond_subtype = SHELLS;
		radius_factor = indsys->opts->shell_r0;
	}

	if (indsys->precond_type == SPARSE) {
		/*
		 indsys->diagL = (double *)MattAlloc(indsys->num_fils,sizeof(double));
		 get_selfs(indsys->diagL, indsys);
		 */
		create_sparMatrix(indsys);
		/*    get_selfs(indsys);*/
	}
}

void get_selfs(SYS *indsys) {
	int j, filnum_j;
	FILAMENT *fil_j;
	SEGMENT *seg1;
	double *diagL = indsys->diagL;

	for (seg1 = indsys->segment; seg1 != NULL; seg1 = seg1->next) {
		for (j = 0; j < seg1->num_fils; j++) {
			fil_j = &(seg1->filaments[j]);
			filnum_j = fil_j->filnumber;
			diagL[filnum_j] = selfterm(fil_j);
#if SUPERCON == ON
			/* note that this routine is not used */
			if (seg1->lambda != 0.0)
				diagL[filnum_j] += seg1->r2 * fil_j->length / fil_j->area;
#endif
		}
	}
}

void fill_spPre(ssystem *sys, SYS *indsys, double w)
/* double w;  frequency */
{

	int i, j, err, used;
	char *Matrix = indsys->sparMatrix;
	MELEMENT *melem, *melem2, **Mlist = indsys->Mlist;
	CX *elem1, *elem2;
	int num_mesh = indsys->num_mesh;
	double *diagL = indsys->diagL;
	double *R = indsys->R;
	double rtempsum, itempsum;

	spClear(Matrix);

	if (indsys->precond_subtype == DIAGL) {
		fill_diagL(sys, indsys, w);

#if 1==0
		for(i = 0; i < num_mesh; i++) {
			for(j = i; j < num_mesh; j++) {
				rtempsum = 0;
				itempsum = 0;
				used = 0;
				/* brute force search for elements */
				for(melem = Mlist[i]; melem != NULL; melem = melem->mnext) {
					for(melem2 = Mlist[j]; melem2 != NULL; melem2 = melem2->mnext) {
						if (melem2->filindex == melem->filindex) {
							used = 1;
							rtempsum += melem->sign*R[melem->filindex]*melem2->sign;
							itempsum += melem->sign*diagL[melem->filindex]*melem2->sign;
						}
					}
				}
				if (used == 1) {
					(elem1 = (CX *)spGetElement(Matrix,i+1,j+1))->real
					= (elem2 = (CX *)spGetElement(Matrix,j+1,i+1))->real = rtempsum;
					elem1->imag = elem2->imag = w*itempsum;
				}

			}
		}
#endif
	} else if (indsys->precond_subtype == CUBEL
			|| indsys->precond_subtype == SHELLS)
		indPrecond(sys, indsys, w);
	else if (indsys->precond_subtype == SEGML)
		fill_bySegment(sys, indsys, w);
	else {
		printf("Unknown Sparse Preconditioner\n");
		exit(1);
	}

	if (indsys->opts->dumpMats & PRE) {
		concat4(outfname, "Pre", indsys->opts->suffix, ".mat");
		if (spFileMatrix(Matrix, outfname, "Precond", 0, 1, 1) == 0)
			fprintf(stderr, "saving sparse matrix failed\n");
	}

}

/* this is called by choose_and_setup_precond() and also by main() if 
 dont_form_Z is true (fmin = 0) */
void create_sparMatrix(SYS *indsys) {
	int err;

	/*  I thought using real matrices might help, but not really in the end.
	 if (!indsys->dont_form_Z)
	 indsys->sparMatrix = (char *)spCreate(indsys->num_mesh, 1, &err);
	 else
	 indsys->sparMatrix = (char *)spCreate(indsys->num_mesh, 0, &err);
	 */

	indsys->sparMatrix = (char *) spCreate(indsys->num_mesh, 1, &err);

	if (err != 0) {
		fprintf(stderr, "Couldn't create sparse matrix, err %d\n", err);
		exit(1);
	}
}

void fill_bySegment(ssystem *sys, SYS *indsys, double w) {
	SEGMENT *seg;
	int max_fils;

	typedef struct _filcube {
		charge *pchg;
		cube *fc;
		int cindex;
	} Filcube;

	Filcube *fillist;
	double **mat;
	double minx = sys->minx;
	double miny = sys->miny;
	double minz = sys->minz;
	double length = sys->length;
	int depth = sys->depth;
	char *Matrix = indsys->sparMatrix;
	double *R = indsys->R;
	double **L = indsys->Z;
	int fili, filj, i, j, k, num_fils, nsize, cindex;
	int xindex, yindex, zindex;
	cube *nc;
	charge *pchg, **nc_pc;
	MELEMENT *mtranj, *mtrani;
	MELEMENT **Mtrans = indsys->Mtrans;
	int debug = 0;
	int ismulti;

	max_fils = 0;
	for (seg = indsys->segment; seg != NULL; seg = seg->next)
		if (seg->num_fils > max_fils)
			max_fils = seg->num_fils;

	fillist = (Filcube *) malloc(max_fils * sizeof(Filcube));
	mat = (double **) malloc(max_fils * sizeof(double *));
	for (i = 0; i < max_fils; i++)
		mat[i] = (double *) malloc(max_fils * sizeof(double));

	ismulti = (indsys->opts->mat_vect_prod == MULTIPOLE);

	for (seg = indsys->segment; seg != NULL; seg = seg->next) {
		num_fils = seg->num_fils;
		if (ismulti) {
			for (i = 0; i < num_fils; i++) {
				pchg = fillist[i].pchg = seg->filaments[i].pchg;
				xindex = (pchg->x - minx) / length;
				yindex = (pchg->y - miny) / length;
				zindex = (pchg->z - minz) / length;
				nc = fillist[i].fc = sys->cubes[depth][xindex][yindex][zindex];
				nsize = nc->directnumeles[0];
				nc_pc = nc->chgs;
				j = 0;
				while (j < nsize && nc_pc[j] != pchg)
					j++;
				if (j == nsize) {
					fprintf(stderr,
							"fill_bySegment: Hey, charge isn't in cube!\n");
					exit(1);
				}
				fillist[i].cindex = j;
			}
			/* let's go and find all the direct terms for this segment */
			for (i = 0; i < num_fils; i++) {
				nc = fillist[i].fc;
				cindex = fillist[i].cindex;
				for (j = i; j < num_fils; j++) {
					mat[i][j] = mat[j][i] = 0;
					if (nc == fillist[j].fc)
						/* both fils are in the same cube */
						mat[i][j] = mat[j][i] =
								nc->directmats[0][cindex][fillist[j].cindex];
					else {
						k = 0;
						/* lets find which cube the other fil is in */
						while (k < nc->numnbrs && nc->nbrs[k] != fillist[j].fc)
							k++;
						if (k < nc->numnbrs) {
							mat[i][j] =
									mat[j][i] =
											(nc->directmats[k])[cindex][fillist[j].cindex];
						} else {
							if (debug == 1)
								printf("Can't find other fil\n");
						}
					}
				}
			}
		} else {
			for (i = 0; i < num_fils; i++) {
				fili = seg->filaments[i].filnumber;
				for (j = i; j < num_fils; j++) {
					filj = seg->filaments[j].filnumber;
					mat[i][j] = mat[j][i] = L[fili][filj];
				}
			}
		}

		for (i = 0; i < num_fils; i++) {
			if (ismulti)
				fili = fillist[i].pchg->fil->filnumber;
			else
				fili = seg->filaments[i].filnumber;
			for (j = 0; j < num_fils; j++) {
				if (ismulti)
					filj = fillist[j].pchg->fil->filnumber;
				else
					filj = seg->filaments[j].filnumber;
				for (mtranj = Mtrans[filj]; mtranj != NULL;
						mtranj = mtranj->mnext) {
					for (mtrani = Mtrans[fili]; mtrani != NULL;
							mtrani = mtrani->mnext) {
						((CX *) spGetElement(Matrix, mtrani->filindex + 1,
								mtranj->filindex + 1))->imag += w * mtrani->sign
								* mat[i][j] * mtranj->sign;
						if (i == j)
							((CX *) spGetElement(Matrix, mtrani->filindex + 1,
									mtranj->filindex + 1))->real += mtrani->sign
									* R[fili] * mtranj->sign;
					}
				}
			}
		}

	}

	free(fillist);
	for (i = 0; i < max_fils; i++)
		free(mat[i]);
	free(mat);

}

/* this is mostly called by fill_spPre() to form the precondtioner */
/* it used be called from main() if we are solving by LU decomposition
 and fmin=0 (the MZMt matrix will be a sparse MRMt).  but
 now has it's own function, fill_diagR()
 */
void fill_diagL(ssystem *sys, SYS *indsys, double w) {
	SEGMENT *seg;
	double val, valR;
	double minx = sys->minx;
	double miny = sys->miny;
	double minz = sys->minz;
	double length = sys->length;
	int depth = sys->depth;
	char *Matrix = indsys->sparMatrix;
	double *R = indsys->R;
	double **L = indsys->Z;
	int fili, filj, i, j, k, num_fils, nsize, cindex;
	int xindex, yindex, zindex;
	cube *nc;
	charge *pchg, **nc_pc;
	MELEMENT *mtranj, *mtrani;
	MELEMENT **Mtrans = indsys->Mtrans;
	int dont_form_Z = indsys->dont_form_Z;
	CX *elem;
	int ismulti, filnum;

	ismulti = (indsys->opts->mat_vect_prod == MULTIPOLE
			&& indsys->opts->soln_technique != LUDECOMP);

	for (seg = indsys->segment; seg != NULL; seg = seg->next) {
		num_fils = seg->num_fils;
		for (i = 0; i < num_fils; i++) {
			filnum = seg->filaments[i].filnumber;
			if (ismulti) {
				pchg = seg->filaments[i].pchg;
				xindex = (pchg->x - minx) / length;
				yindex = (pchg->y - miny) / length;
				zindex = (pchg->z - minz) / length;
				nc = sys->cubes[depth][xindex][yindex][zindex];
				nsize = nc->directnumeles[0];
				nc_pc = nc->chgs;
				j = 0;
				while (j < nsize && nc_pc[j] != pchg)
					j++;
				if (j == nsize) {
					fprintf(stderr,
							"fill_bySegment: Hey, charge isn't in cube!\n");
					exit(1);
				}
				cindex = j;

				val = nc->directmats[0][cindex][cindex];
			} else {
				if (!dont_form_Z)
					val = L[filnum][filnum];
				else
					val = 0;
			}
			valR = R[filnum];
			for (mtranj = Mtrans[filnum]; mtranj != NULL; mtranj =
					mtranj->mnext) {
				for (mtrani = Mtrans[filnum]; mtrani != NULL;
						mtrani = mtrani->mnext) {
					(elem = (CX *) spGetElement(Matrix, mtrani->filindex + 1,
							mtranj->filindex + 1))->imag += w * mtrani->sign
							* val * mtranj->sign;
					elem->real += mtrani->sign * valR * mtranj->sign;
				}
			}
		}
	}
}

/*
 Called from main() if we are solving by LU decomposition
 and fmin=0 (the MZMt matrix will be a sparse MRMt).

 It puts the resistance matrix into a sparse matrix structure for a quick
 solve.  I thought that separating it out and using real matrices
 for this part would increase speed, but it didn't seem to.
 The reordering and fillin manipulation must be dominating and
 since we only factor this once, there is no real benefit.
 */
void fill_diagR(SYS *indsys) {
	SEGMENT *seg;
	double val, valR;
	char *Matrix = indsys->sparMatrix;
	double *R = indsys->R;
	int fili, filj, i, j, k, num_fils, nsize, cindex;
	MELEMENT *mtranj, *mtrani;
	MELEMENT **Mtrans = indsys->Mtrans;
	CX *elem;
	int ismulti, filnum;

	for (seg = indsys->segment; seg != NULL; seg = seg->next) {
		num_fils = seg->num_fils;
		for (i = 0; i < num_fils; i++) {
			filnum = seg->filaments[i].filnumber;
			valR = R[filnum];
			for (mtranj = Mtrans[filnum]; mtranj != NULL; mtranj =
					mtranj->mnext) {
				for (mtrani = Mtrans[filnum]; mtrani != NULL;
						mtrani = mtrani->mnext) {
					*(spGetElement(Matrix, mtrani->filindex + 1,
							mtranj->filindex + 1)) += mtrani->sign * valR
							* mtranj->sign;
				}
			}
		}
	}
}

/* computes  mu/(4 pi r0) * (l_i,l_j) where (l_i,l_j) is the dot product
 of the vectors from one end of each fil to the other end and r0 is
 the multipole cube side length times a factor (radius_factor) which
 is specified on the command line or defaults to 0.87 = sqrt(3)*0.5
 so that the sphere encompasses the whole cube
 */
double shift_mutual(FILAMENT *fil_i, FILAMENT *fil_j, ssystem *sys) {
	/* a factor to scale r0 */
	/* static double factor = 0.87; *//*(sqrt( 3) / 2*/

	return dotprod(fil_i, fil_j) * MUOVER4PI / (sys->length * radius_factor);
}
