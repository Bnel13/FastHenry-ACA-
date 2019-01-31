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


/* # ***** sort to /src/direct
 # ***** */
//#include "mulGlobal.h" //BAPN ed out
#include "induct.h"//BAPN change

/* SRW */
double **Q2PDiag(charge**, int, int*, int);
double **Q2P(charge**, int, int*, charge**, int, int);
double **Q2Pfull(cube*, int);
double **ludecomp(double**, int, int);
void solve(double**, double*, double*, int);
void invert(double**, int, int*);
int compressMat(double**, int, int*, int);
void expandMat(double**, int, int, int*, int);
void matcheck(double**, int, int);
void matlabDump(double**, int, char*);

double **Q2PDiag(charge **chgs, int numchgs, int *is_dummy, int calc) {
	double **mat = NULL;
	int i, j;
	//start BAPN
	/*
	int one = 1;
	int two = 0;
	int one_c = 1;
	double mel = 0;
	int two_c = 0;
	charge *my_one = chgs[0];
	charge *myy = chgs[0];
	FILAMENT *my_onef = my_one->fil;
	my_onef->lenvect;
	charge **group_one; 
	charge **group_two; 
	group_one = malloc(10 * sizeof(charge));
	group_two = malloc(10 * sizeof(charge));
	group_one[0] = myy;
	for (int j = 0; j < numchgs - 1; j++)
	{
		charge *myC = chgs[j];
		FILAMENT *myCC = myC->fil;

		charge *myB = chgs[j + 1];
		FILAMENT *myBB = myB->fil;
		if (fabs(vdotp(myCC->lenvect, myBB->lenvect))
			/ (myCC->length * myBB->length) < EPS)
		{
			if (one == 1)
			{
				one = 0;
				group_two[two_c] = myB;
				two_c++;
				two = 1;
			}
			else
			{
				two = 0;
				group_one[one_c] = myB;
				one_c++;
				one = 1;
			}
		}
		else
		{
			if (one == 0)
			{
				group_two[two_c] = myB;
				two_c++;
			}
			else
			{
				group_one[one_c] = myB;
				one_c++;
			}

		}
	}
	for (int i = 0; i < one_c-1; i++)
	{
		for (int j = 0; j < two_c - 1; j++)
		{
			mel = calcp(group_one[i], group_two[j], NULL);
		}
	}*/
	//end BAPN
	/* Allocate storage for the potential coefficients. */
	CALLOC(mat, numchgs, double*, ON, AQ2PD);
	for (i = 0; i < numchgs; i++)
		CALLOC(mat[i], numchgs, double, ON, AQ2PD);
	
	if (calc) {
		/* Compute the potential coeffs. */
		/* - exclude dummy panels when they would need to carry charge
		 - exclude dielec i/f panels when they would lead to evals at their
		 centers (only if using two-point flux-den-diff evaluations) */
		for (i = 0; i < numchgs; i++) {
#if NUMDPT == 2
			if (chgs[i]->dummy)
				; /* don't check surface of a dummy */
			else if (chgs[i]->surf->type == DIELEC
					|| chgs[i]->surf->type == BOTH)
				continue;
#endif
			for (j = 0; j < numchgs; j++) { /* need to have charge on them */
#if SKIPQD == ON
				if(chgs[j]->pos_dummy == chgs[i] || chgs[j]->neg_dummy == chgs[i])
				continue;
#endif
				if (!is_dummy[j])
					mat[i][j] = calcp(chgs[j], chgs[i], NULL);
				
			}
		}
	}

#if DSQ2PD == ON
	dispQ2PDiag(mat, chgs, numchgs, is_dummy);
#endif

	return (mat);
}

double **Q2P(charge **qchgs, int numqchgs, int *is_dummy, charge **pchgs,
		int numpchgs, int calc) {
	double **mat = NULL;
	int i, j;

	/* Allocate storage for the potential coefficients. P rows by Q cols. */
	CALLOC(mat, numpchgs, double*, ON, AQ2P);
	for (i = 0; i < numpchgs; i++) {
		CALLOC(mat[i], numqchgs, double, ON, AQ2P);
		if (calc) {
			/* exclude:
			 - dummy panels in the charge list
			 - dielectric i/f panels in the eval list (if doing 2-point E's)*/
#if NUMDPT == 2
			if (pchgs[i]->dummy)
				; /* don't check the surface of a dummy */
			else if (pchgs[i]->surf->type == DIELEC
					|| pchgs[i]->surf->type == BOTH)
				continue;
#endif
			for (j = 0; j < numqchgs; j++) { /* only dummy panels in the charge list */
				if (!is_dummy[j]) /* (not the eval list) are excluded */
					mat[i][j] = calcp(qchgs[j], pchgs[i], NULL);
				/* old: pchgs[i],qchgs[j] */
			}
		}
	}

#if DISQ2P == ON
	dispQ2P(mat, qchgs, numqchgs, is_dummy, pchgs, numpchgs);
#endif

	return (mat);
}

/*
 used only in conjunction with DMPMAT == ON  and DIRSOL == ON
 to make 1st directlist mat = full P mat
 */
double **Q2Pfull(cube *directlist, int numchgs) {
	int i, j, fromp, fromq, top, toq;
	double **mat = NULL;
	cube *pq, *pp;
	charge **pchgs, **qchgs, *eval;

	/* allocate room for full P matrix */
	MALLOC(mat, numchgs, double*, ON, AQ2P);
	for (i = 0; i < numchgs; i++)
		MALLOC(mat[i], numchgs, double, ON, AQ2P);

	/* load the matrix in the style of Q2P() - no attempt to exploit symmetry */
	/* calculate interaction between every direct list entry and entire dlist */
	for (pp = directlist; pp != NULL; pp = pp->dnext) {
		pchgs = pp->chgs;
		fromp = pchgs[0]->index - 1; /* row index range */
		top = fromp + pp->upnumeles[0];
		for (pq = directlist; pq != NULL; pq = pq->dnext) {
			qchgs = pq->chgs;
			fromq = qchgs[0]->index - 1; /* column index range */
			toq = fromq + pq->upnumeles[0];

			for (i = fromp; i < top; i++) {
				for (j = fromq; j < toq; j++) {
					eval = qchgs[j - fromq];
					mat[i][j] = calcp(pchgs[i - fromp], eval, NULL);
				}
			}

		}
	}
	return (mat);
}

/*
 - returned matrix has L below the diagonal, U above (GVL1 pg 58)
 - if allocate == TRUE ends up storing P and LU (could be a lot)
 */
double **ludecomp(double **matin, int size, int allocate) {
	extern int fulldirops;
	double factor, **mat = NULL;
	int i, j, k;

	if (allocate == TRUE) {
		/* allocate for LU matrix and copy A */
		MALLOC(mat, size, double*, ON, AMSC);
		for (i = 0; i < size; i++) {
			MALLOC(mat[i], size, double, ON, AMSC);
			for (j = 0; j < size; j++)
				mat[i][j] = matin[i][j];
		}
	} else
		mat = matin;

	for (k = 0; k < size - 1; k++) { /* loop on rows */
		if (mat[k][k] == 0.0) {
			fprintf(stderr, "ludecomp: zero piovt\n");
			exit(1);
		}
		for (i = k + 1; i < size; i++) { /* loop on remaining rows */
			factor = (mat[i][k] /= mat[k][k]);
			fulldirops++;
			for (j = k + 1; j < size; j++) { /* loop on remaining columns */
				mat[i][j] -= (factor * mat[k][j]);
				fulldirops++;
			}
		}
	}
	return (mat);
}

/*
 For direct solution of Pq = psi, used if DIRSOL == ON or if preconditioning.
 */
void solve(double **mat, double *x, double *b, int size) {
	extern int fulldirops;
	int i, j;

	/* copy rhs */
	if (x != b)
		for (i = 0; i < size; i++)
			x[i] = b[i];

	/* forward elimination */
	for (i = 0; i < size; i++) { /* loop on pivot row */
		for (j = i + 1; j < size; j++) { /* loop on elimnation row */
			x[j] -= mat[j][i] * x[i];
			fulldirops++;
		}
	}

	/* back substitution */
	for (i--; i > -1; i--) { /* loop on rows */
		for (j = i + 1; j < size; j++) { /* loop on columns */
			x[i] -= mat[i][j] * x[j];
			fulldirops++;
		}
		x[i] /= mat[i][i];
		fulldirops++;
	}
}

/* 
 In-place inverts a matrix using guass-jordan.
 - is_dummy[i] = 0 => ignore row/col i
 */
void invert(double **mat, int size, int *reorder) {
	int i, j, k, best;
	double normal, multiplier, bestval, nextbest;
	/*
	 matlabDump(mat,size,"p");
	 */
	for (i = 0; i < size; i++) {
		best = i;
		bestval = ABS(mat[i][i]);
		for (j = i + 1; j < size; j++) {
			nextbest = ABS(mat[i][j]);
			if (nextbest > bestval) {
				best = j;
				bestval = nextbest;
			}
		}

		/* If reordering, find the best pivot. */
		if (reorder != NULL) {
			reorder[i] = best;
			if (best != i) {
				for (k = 0; k < size; k++) {
					bestval = mat[k][best];
					mat[k][best] = mat[k][i];
					mat[k][i] = bestval;
				}
			}
		}

		/* First i^{th} column of A. */
		normal = 1.0 / mat[i][i];
		for (j = 0; j < size; j++) {
			mat[j][i] *= normal;
		}
		mat[i][i] = normal;

		/* Fix the backward columns. */
		for (j = 0; j < size; j++) {
			if (j != i) {
				multiplier = -mat[i][j];
				for (k = 0; k < size; k++) {
					if (k != i)
						mat[k][j] += mat[k][i] * multiplier;
					else
						mat[k][j] = mat[k][i] * multiplier;
				}
			}
		}
	}

	/* Unravel the reordering, starting with the last column. */
	if (reorder != NULL) {
		for (i = size - 2; i >= 0; i--) {
			if (reorder[i] != i) {
				for (k = 0; k < size; k++) {
					bestval = mat[k][i];
					mat[k][reorder[i]] = mat[k][i];
					mat[k][i] = bestval;
				}
			}
		}
	}
	/*
	 matlabDump(mat,size,"c");
	 */

}

/*
 Used in conjuction with invert() to remove dummy row/columns
 comp_rows = TRUE => remove rows corresponding to ones in is_dummy[]
 comp_rows = FALSE => remove cols corresponding to ones in is_dummy[]
 comp_rows = BOTH => remove both rows and columns
 returns number of rows/cols in compressed matrix
 */
int compressMat(double **mat, int size, int *is_dummy, int comp_rows) {
	static int *cur_order;
	static int cur_order_array_size = 0;
	int cur_order_size, i, j, k;

	if (cur_order_array_size < size) {
		CALLOC(cur_order, size, int, ON, AMSC);
	}

	/* figure the new order vector (cur_order[i] = index of ith row/col) */
	for (i = cur_order_size = 0; i < size; i++) {
		if (!is_dummy[i])
			cur_order[cur_order_size++] = i;
	}

	if (comp_rows == TRUE || comp_rows == BOTH) {
		/* compress by removing rows from the matrix */
		for (i = 0; i < cur_order_size; i++) {
			if ((k = cur_order[i]) == i)
				continue; /* if not reordered */
			for (j = 0; j < size; j++) { /* copy the row to its new location */
				mat[i][j] = mat[k][j];
			}
		}
	}
	if (comp_rows == FALSE || comp_rows == BOTH) {
		/* compress by removing columns from the matrix */
		for (j = 0; j < cur_order_size; j++) {
			if ((k = cur_order[j]) == j)
				continue; /* if not reordered */
			for (i = 0; i < size; i++) { /* copy the col to its new location */
				mat[i][j] = mat[i][k];
			}
		}
	}
	return (cur_order_size);
}

/*
 Used in conjuction with invert() to add dummy row/columns
 exp_rows = TRUE => add rows corresponding to ones in is_dummy[]
 exp_rows = FALSE => add cols corresponding to ones in is_dummy[]
 exp_rows = BOTH => add rows and columns
 */
void expandMat(double **mat, int size, int comp_size, int *is_dummy,
		int exp_rows) {
	int i, j, k, next_rc;

	if (exp_rows == TRUE || exp_rows == BOTH) {
		next_rc = comp_size - 1;
		/* add rows to the matrix starting from the bottom */
		for (i = size - 1; i >= 0; i--) {
			if (is_dummy[i]) { /* zero out dummy row */
				for (j = 0; j < size; j++)
					mat[i][j] = 0.0;
			} else { /* copy the row from its compressed location */
				for (j = 0; j < size; j++)
					mat[i][j] = mat[next_rc][j];
				next_rc--;
			}
		}
	}
	if (exp_rows == FALSE || exp_rows == BOTH) {
		next_rc = comp_size - 1;
		/* add columns to the matrix starting from the right */
		for (j = size - 1; j >= 0; j--) {
			if (is_dummy[j]) { /* zero out dummy column */
				for (i = 0; i < size; i++)
					mat[i][j] = 0.0;
			} else { /* copy the col from its compressed location */
				for (i = 0; i < size; i++)
					mat[i][j] = mat[i][next_rc];
				next_rc--;
			}
		}
	}
}

/* 
 Checks to see if the matrix has the M-matrix sign pattern and if
 it is diagonally dominant.
 */
void matcheck(double **mat, int rows, int size) {
	double rowsum;
	int i, j;

	for (i = rows - 1; i >= 0; i--) {
		for (rowsum = 0.0, j = size - 1; j >= 0; j--) {
			if ((i != j) && (mat[i][j] > 0.0)) {
				printf("violation mat[%d][%d] =%g\n", i, j, mat[i][j]);
			}
			if (i != j)
				rowsum += ABS(mat[i][j]);
		}
		printf("row %d diag=%g rowsum=%g\n", i, mat[i][i], rowsum);
		if (rowsum > mat[i][i]) {
			for (j = size - 1; j >= 0; j--) {
				printf("col%d = %g ", j, mat[i][j]);
			}
			printf("\n");
		}
	}
}

void matlabDump(double **mat, int size, char *name) {
	FILE *foo;
	int i, j;
	char fname[100];

	sprintf(fname, "%s.m", name);
	/* SRW -- this is ascii data */
	foo = fopen(fname, "w");
	fprintf(foo, "%s = [\n", name);
	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			fprintf(foo, "%.10e  ", mat[i][j]);
		}
		fprintf(foo, "\n");
	}
	fprintf(foo, "]\n");
}

