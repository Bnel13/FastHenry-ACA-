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



/* this sets up the vectors to call Fastcap's ComputePsi */
/* It will be called twice for each coordinate direction.  Once for real
and once for imaginary */

#include "induct.h"
#include "spMatrix.h"
#include <time.h>//BAP include
#define MLACA 1;//BAP change
#define MLACA2 0;//BAP change

/* SRW */
void SetupComputePsi(CX*, ssystem*, CX*, int, charge*, double, double*, SYS*);
void realmatCXvec(CX*, double**, CX*, int);
void fixEvalDirect(charge**, int, int*, charge**, int, double**);


/*The preconditioner is still constructed with the same inputs as the original FastHenry.*/

/* Vs will contain the result,  Im is the 'q',  Size is the size of vectors. */
/* This will alter Im.  Im = Precond*Im */
void SetupComputePsi(CX *Vs, ssystem *sys, CX *Im, int size, charge *chglist,
	double w, double *R, SYS *indsys)
	/* double w;  radian frequency */
	/* double *R; resistance vector */
{

	extern double dirtime;
	double *q, *p;
	static CX *Ib = NULL, *Vb = NULL, *Vdirect = NULL, *ctemp;
	int branches;
	CX temp;
	MELEMENT *mtemp;
	charge *chg;
	int i, j;
	double rtemp;
	MELEMENT **Mtrans, **Mlist;
	double maxdiff, pdiff;
	int maxindx;
	int ind_opcnt_mult = 0, ind_opcnt_real = 0;

	branches = indsys->num_fils;
	Mtrans = indsys->Mtrans;
	Mlist = indsys->Mlist;

	if (Ib == NULL) {
		Ib = (CX *)MattAlloc(branches, sizeof(CX));
		Vb = (CX *)MattAlloc(branches, sizeof(CX));
		ctemp = (CX *)MattAlloc(size, sizeof(CX));
#ifndef NODEBUG
		Vdirect = (CX *)MattAlloc(branches, sizeof(CX));
#endif
	}

	for (i = 0; i < branches; i++)
		Vb[i] = CXZERO;

	q = sys->q;
	p = sys->p;
	ASSERT(size == indsys->num_mesh);

	if (indsys->precond_type == LOC) {
		multPrecond(indsys->Precond, Im, ctemp, size);
		for (i = 0; i < size; i++)
			Im[i] = ctemp[i];
	}
	else if (indsys->precond_type == SPARSE)
		spSolve(indsys->sparMatrix, (spREAL*)Im, (spREAL*)Im);

	/* do  Ib = Mtrans*Im */
	for (i = 0; i < branches; i++) {
		Ib[i] = CXZERO;
		for (mtemp = Mtrans[i]; mtemp != NULL; mtemp = mtemp->mnext) {
			if (mtemp->sign == 1)
				cx_add(Ib[i], Ib[i], Im[mtemp->filindex]);
			else
				cx_sub(Ib[i], Ib[i], Im[mtemp->filindex]);
		}
	}

	/* Evaluate M*L*Mt*Im = M*L*Ib using the multipole algorithm */

	/* Do all of the non-direct parts first */
	sys->DirectEval = FALSE;

	//Start BAP change

	int *p_index = (int *)calloc(indsys->num_fils, sizeof(int));
	for (chg = chglist; chg != NULL; chg = chg->next) {
		p_index[chg->fil->filnumber] = chg->index;
	}
	/*
	for (i = 0; i < branches; i++)
	p[i] = 0;
	int track = 0;
	for (int a = 0; a < sys->depth - 1; a++)
	{
	for (int b = 0; b < sys->number_level_groups[a + 1]; b++)
	{

	for (int k = 0; k < sys->group_interaction_number[a][b]; k++)
	{
	//	FILE *matrix = fopen("matrix.txt", "w");
	for (int i = 0; i < sys->kid_num[a + 1][b]; i++)
	{

	for (int j = 0; j < sys->kid_num[a + 1][sys->interactions[a][b][k]]; j++)
	{
	double sum = 0;
	for (int nn = 0; nn < sys->k[track]; nn++)
	{
	sum += sys->U[track][i][nn] * sys->V[track][nn][j];
	}
	double or = (sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->lenvect[0] / sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->length)*(sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->lenvect[0] / sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->length);
	or += (sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->lenvect[1] / sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->length)*(sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->lenvect[1] / sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->length);
	or += (sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->lenvect[2] / sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->length)*(sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->lenvect[2] / sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->length);
	p[p_index[sys->child_index[a + 1][b][i] - 1]] += sum* or;
	p[p_index[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]] += sum* or;
	}
	}
	++track;
	}

	}


	}
	*/
	//End BAP change
#if MLACA2 
	for (i = 0; i < 3; i++) { /* for each of the coordinate directions */

							  /* do the real part */
		for (chg = chglist; chg != NULL; chg = chg->next) {
			/* fill the pseudo-charge vector */
			q[chg->index] = Ib[chg->fil->filnumber].real * chg->fil->lenvect[i];
#if OPCNT == ON
			ind_opcnt_mult++;
#endif

		}
		computePsi(sys, q, p, branches, chglist); //Bap comment out

												  //Start BAP change
#endif
#if MLACA //start comment bap 

		int track = 0;
		double *VlI_real = (double *)malloc(sys->k_max * sizeof(double));
		double *VlI_imag = (double *)malloc(sys->k_max * sizeof(double));
		double *UlI_real_s = (double *)malloc(sys->k_max * sizeof(double));
		double *UlI_imag_s = (double *)malloc(sys->k_max * sizeof(double));

		double *UVlI_real = (double *)malloc(sys->el_max * sizeof(double));
		double *UVlI_imag = (double *)malloc(sys->el_max * sizeof(double));

		double *VUlI_real_s = (double *)malloc(sys->el_max * sizeof(double));
		double *VUlI_imag_s = (double *)malloc(sys->el_max * sizeof(double));

		double mul_lI_real;
		double mul_lI_imag;

		double mul_lIt_real;
		double mul_lIt_imag;

		int to_count;
		int *to_do = (int *)malloc(sys->el_max * sizeof(int));

		clock_t t;
		t = clock();
		for (int a = 0; a < sys->depth; a++)
		{
			for (int b = 0; b < sys->number_level_groups[a]; b++)
			{
				int kid_numm = sys->kid_num[a][b];
				for (int k = 0; k < sys->group_interaction_number[a][b]; k++)
				{
					//	FILE *matrix = fopen("matrix.txt", "w");




					int interact = sys->interactions[a][b][k];
					int kid_numm_i = sys->kid_num[a][interact];
					for (int c = 0; c < 3; c++)
					{

						memset(VlI_real, 0, sys->k[track] * sizeof(double));
						memset(VlI_imag, 0, sys->k[track] * sizeof(double));
						memset(UlI_real_s, 0, sys->k[track] * sizeof(double));
						memset(UlI_imag_s, 0, sys->k[track] * sizeof(double));

						memset(UVlI_real, 0, kid_numm * sizeof(double));
						memset(UVlI_imag, 0, kid_numm * sizeof(double));

						memset(VUlI_real_s, 0, kid_numm_i * sizeof(double));
						memset(VUlI_imag_s, 0, kid_numm_i * sizeof(double));


						to_count = 0;

						for (int j = 0; j < kid_numm_i; j++)
						{
							/*specific l x specific I*/

							//sum = (sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->lenvect[c] / sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->length)*Ib[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1].real;
							//sum = (sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->lenvect[c] / sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->length)*Ib[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1].imag;
							//sum = (sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->lenvect[c] / sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->length)*Ib[sys->child_index[a + 1][b][i] - 1].real;
							//sum = (sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->lenvect[c] / sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->length)*Ib[sys->child_index[a + 1][b][i] - 1].imag;
							if (sys->fill[sys->child_index[a][interact][j] - 1]->fil->lenvect[c] != 0)
							{
								to_do[to_count] = j;
								++to_count;

								mul_lI_real = sys->fill[sys->child_index[a][interact][j] - 1]->fil->lenvect[c] * Ib[sys->child_index[a][interact][j] - 1].real;
								mul_lI_imag = sys->fill[sys->child_index[a][interact][j] - 1]->fil->lenvect[c] * Ib[sys->child_index[a][interact][j] - 1].imag;

								for (int nn = 0; nn < sys->k[track]; nn++)
								{
									//sum += sys->U[track][i][nn] * sys->V[track][nn][j];

									VlI_real[nn] += sys->V[track][nn][j] * mul_lI_real;/*(sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->lenvect[c]) /* sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->length)*Ib[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1].real;*/
									VlI_imag[nn] += sys->V[track][nn][j] * mul_lI_imag; /*(sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->lenvect[c]) /* sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->length)*Ib[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1].imag;*/

																						// place here V*(l*I) 
								}


								//double or = (sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->lenvect[0] / sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->length)*(sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->lenvect[0] / sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->length);
								//or += (sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->lenvect[1] / sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->length)*(sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->lenvect[1] / sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->length);
								//or += (sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->lenvect[2] / sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->length)*(sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->lenvect[2] / sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->length);
								//p[p_index[sys->child_index[a + 1][b][i] -1]] += sum * or * Ib[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1].real;
								//p[p_index[sys->child_index[a + 1][sys->interactions[a][b][k]][j] -1]] += sum * or * Ib[sys->child_index[a + 1][b][i] - 1].real;
								//Vb[sys->child_index[a + 1][b][i] - 1].real += sum * or *Ib[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1].real;
								//Vb[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1].real += sum * or *Ib[sys->child_index[a + 1][b][i] - 1].real;
								//Vb[sys->child_index[a + 1][b][i] - 1].imag += sum * or *Ib[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1].imag;
								//Vb[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1].imag += sum * or *Ib[sys->child_index[a + 1][b][i] - 1].imag;
							}
						}
						//double *UVlI_real = (double *)calloc(sys->kid_num[a + 1][b], sizeof(double));
						//double *UVlI_imag = (double *)calloc(sys->kid_num[a + 1][b], sizeof(double));
						for (int i = 0; i < kid_numm; i++)
						{
							/* after n loop is complete it becomes possible */
							if (sys->fill[sys->child_index[a][b][i] - 1]->fil->lenvect[c] != 0)
							{
								mul_lIt_real = (sys->fill[sys->child_index[a][b][i] - 1]->fil->lenvect[c])*Ib[sys->child_index[a][b][i] - 1].real;//second
								mul_lIt_imag = (sys->fill[sys->child_index[a][b][i] - 1]->fil->lenvect[c])*Ib[sys->child_index[a][b][i] - 1].imag;//second
								for (int nn = 0; nn < sys->k[track]; nn++)
								{

									//start second
									UlI_real_s[nn] += sys->U[track][nn][i] * mul_lIt_real;/*(sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->lenvect[c])/* / sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->length)*Ib[sys->child_index[a + 1][b][i] - 1].real;*/
									UlI_imag_s[nn] += sys->U[track][nn][i] * mul_lIt_imag;/* (sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->lenvect[c]) /* sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->length)*Ib[sys->child_index[a + 1][b][i] - 1].imag;*/
																						  //end second


									UVlI_real[i] += sys->U[track][nn][i] * VlI_real[nn];
									UVlI_imag[i] += sys->U[track][nn][i] * VlI_imag[nn];



									//place U x (V x ln x In)
									//sys->U[track][i][nn]
								}

								//lm x (Um x (V x ln x In)/ sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->length) * UVlI_real[sys->child_index[a + 1][b][i] - 1];

								Vb[sys->child_index[a][b][i] - 1].imag += (sys->fill[sys->child_index[a][b][i] - 1]->fil->lenvect[c]) /* sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->length)*/ * UVlI_imag[i] * MUOVER4PI;
								Vb[sys->child_index[a][b][i] - 1].real += (sys->fill[sys->child_index[a][b][i] - 1]->fil->lenvect[c]) /* sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->length)*/ * UVlI_real[i] * MUOVER4PI;




							}
						}

						//free(UVlI_real);
						//free(UVlI_imag);

						//for (int i = 0; i < sys->kid_num[a + 1][b]; i++)
						//{
						//if (sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->lenvect[c] != 0)
						//{

						//}
						//mul_Il_real = (sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->lenvect[c])*Ib[sys->child_index[a + 1][b][i] - 1].real;
						//mul_Il_imag = (sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->lenvect[c])*Ib[sys->child_index[a + 1][b][i] - 1].imag;
						//for (int nn = 0; nn < sys->k[track]; nn++)
						//{
						//	IlU_real_s[nn] += sys->U[track][i][nn] * mul_Il_real;/*(sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->lenvect[c])/* / sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->length)*Ib[sys->child_index[a + 1][b][i] - 1].real;*/
						//IlU_imag_s[nn] += sys->U[track][i][nn] * mul_Il_imag;/* (sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->lenvect[c]) /* sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->length)*Ib[sys->child_index[a + 1][b][i] - 1].imag;*/
						//}
						//}
						//}
						/*double *IlUV_real_s = (double *)calloc(sys->kid_num[a + 1][sys->interactions[a][b][k]], sizeof(double));*/
						//double *IlUV_imag_s = (double *)calloc(sys->kid_num[a + 1][sys->interactions[a][b][k]], sizeof(double));
						int j;
						//int interacti = sys->interactions[a][b][k];
						//for (int w = 0; w < to_count; w++)
						for (int j = 0; j < sys->kid_num[a][interact]; j++)
						{

							//j = to_do[w];
							if (sys->fill[sys->child_index[a][sys->interactions[a][b][k]][j] - 1]->fil->lenvect[c] != 0)
							{
								//
								for (int nn = 0; nn < sys->k[track]; nn++)
								{
									VUlI_real_s[j] += UlI_real_s[nn] * sys->V[track][nn][j];
									VUlI_imag_s[j] += UlI_imag_s[nn] * sys->V[track][nn][j];
								}

								Vb[sys->child_index[a][interact][j] - 1].real += (sys->fill[sys->child_index[a][interact][j] - 1]->fil->lenvect[c]) /* sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->length)*/* VUlI_real_s[j] * MUOVER4PI;
								Vb[sys->child_index[a][interact][j] - 1].imag += (sys->fill[sys->child_index[a][interact][j] - 1]->fil->lenvect[c]) /* sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->length)*/* VUlI_imag_s[j] * MUOVER4PI;
							}
						}
						//free(IlUV_real_s);
						//free(IlUV_imag_s);




					}

					++track;


				}

			}


		}
		t = clock() - t;
		double time_taken = ((double)t) / CLOCKS_PER_SEC; // in seconds

		free(VlI_real);
		free(VlI_imag);
		free(UlI_real_s);
		free(UlI_imag_s);

		free(UVlI_real);
		free(UVlI_imag);

		free(VUlI_real_s);
		free(VUlI_imag_s);

		free(to_do);



		/*
		//Start BAP fast
		int ab = 0;
		int V_count = 0;
		int U_count = 0;
		int tr = 0;
		for (int a = 0; a < sys->depth; a++)
		{
		for (int b = 0; b < sys->number_level_groups[a]; b++)
		{
		++ab;


		for (int k = 0; k < sys->group_interaction_number[a][b]; k++)
		{
		for (int i = 0; i < sys->kid_num[a][sys->interactions[a][b][k]]; i++)
		{
		for (int j = 0; j < sys->k[tr]; j++)
		{
		++V_count;
		}
		}
		for (int i = 0; i < sys->kid_num[a][b]; i++)
		{
		for (int j = 0; j < sys->k[tr]; j++)
		{
		++U_count;
		}
		}
		++tr;
		}
		}
		}



		double *VV = (double *)malloc(V_count * sizeof(double));
		double *UU = (double *)malloc(U_count * sizeof(double));

		int V_count_new = 0;
		int U_count_new = 0;


		int *K_new = (int *)malloc(track * sizeof(int));
		int *kid_num_new = (int *)malloc(ab*sizeof(int));
		int *interactions_new = (int *)malloc(track *sizeof(int));
		int ab_new = 0;
		int track_new = 0;
		double **fill_len_new = (double **)malloc(indsys->num_fils*sizeof(double *));
		for (int i = 0; i < indsys->num_fils; i++)
		{
		fill_len_new[i] = (double *)malloc(3 * sizeof(double));
		for (int j = 0; j < 3; j++)
		{
		fill_len_new[i][j] = sys->fill[i]->fil->lenvect[j];
		}
		}



		for (int a = 0; a < sys->depth; a++)
		{
		for (int b = 0; b < sys->number_level_groups[a]; b++)
		{
		kid_num_new[ab_new] = sys->kid_num[a][b];
		++ab_new;
		for (int k = 0; k < sys->group_interaction_number[a][b]; k++)
		{
		interactions_new[track_new] = sys->interactions[a][b][k];
		K_new[track_new] = sys->k[track_new];


		for (int i = 0; i < sys->kid_num[a][sys->interactions[a][b][k]]; i++)
		{
		for (int j = 0; j < K_new[track_new]; j++)
		{
		VV[V_count_new] = sys->V[track_new][j][i];
		++V_count_new;
		}

		}
		for (int i = 0; i < sys->kid_num[a][b]; i++)
		{
		for (int j = 0; j < K_new[track_new]; j++)
		{
		UU[U_count_new] = sys->U[track_new][j][i];
		++U_count_new;
		}

		}
		++track_new;
		}
		}
		}
		//End BAP fast

		//start speed test BAP
		int track_speed = 0;
		int V_count_speed_1 = 0;
		int V_count_speed_2 = 0;
		int U_count_speed = 0;
		double *VlI_real_speed = (double *)malloc(sys->k_max * sizeof(double));
		double *VlI_imag_speed = (double *)malloc(sys->k_max * sizeof(double));
		double *UlI_real_s_speed = (double *)malloc(sys->k_max * sizeof(double));
		double *UlI_imag_s_speed = (double *)malloc(sys->k_max * sizeof(double));

		double *UVlI_real_speed = (double *)malloc(sys->el_max * sizeof(double));
		double *UVlI_imag_speed = (double *)malloc(sys->el_max * sizeof(double));

		double *VUlI_real_s_speed = (double *)malloc(sys->el_max * sizeof(double));
		double *VUlI_imag_s_speed = (double *)malloc(sys->el_max * sizeof(double));

		double mul_lI_real_speed;
		double mul_lI_imag_speed;

		double mul_lIt_real_speed;
		double mul_lIt_imag_speed;

		int to_count_speed;
		int *to_do_speed = (int *)malloc(sys->el_max * sizeof(int));


		for (int a = 0; a < sys->depth; a++)
		{
		for (int b = 0; b < sys->number_level_groups[a]; b++)
		{
		int kid_numm = sys->kid_num[a][b];
		for (int k = 0; k < sys->group_interaction_number[a][b]; k++)
		{
		//	FILE *matrix = fopen("matrix.txt", "w");




		int interact = sys->interactions[a][b][k];
		int kid_numm_i = sys->kid_num[a][interact];
		int VVV = V_count_speed_1;
		int UUU = U_count_speed;
		for (int c = 0; c < 3; c++)
		{
		memset(VlI_real_speed, 0, K_new[track_speed] * sizeof(double));
		memset(VlI_imag_speed, 0, K_new[track_speed] * sizeof(double));
		memset(UlI_real_s_speed, 0, K_new[track_speed] * sizeof(double));
		memset(UlI_imag_s_speed, 0, K_new[track_speed] * sizeof(double));

		memset(UVlI_real_speed, 0, kid_numm * sizeof(double));
		memset(UVlI_imag_speed, 0, kid_numm * sizeof(double));

		memset(VUlI_real_s_speed, 0, kid_numm_i * sizeof(double));
		memset(VUlI_imag_s_speed, 0, kid_numm_i * sizeof(double));

		to_count_speed = 0;

		if (c != 2)
		{
		V_count_speed_1 = VVV;
		V_count_speed_2 = VVV;
		U_count_speed = UUU;
		}



		for (int j = 0; j < kid_numm_i; j++)
		{

		if (fill_len_new[sys->child_index[a][interact][j] - 1][c] != 0)
		{
		to_do_speed[to_count] = j;
		++to_count_speed;

		mul_lI_real_speed = (fill_len_new[sys->child_index[a][interact][j] - 1][c])*Ib[sys->child_index[a][interact][j] - 1].real;
		mul_lI_imag_speed = (fill_len_new[sys->child_index[a][interact][j] - 1][c])*Ib[sys->child_index[a][interact][j] - 1].imag;

		for (int nn = 0; nn < K_new[track_speed]; nn++)
		{

		VlI_real_speed[nn] += VV[V_count_speed_1] * mul_lI_real_speed;
		VlI_imag_speed[nn] += VV[V_count_speed_1] * mul_lI_imag_speed;



		}
		}
		++V_count_speed_1;
		}



		for (int i = 0; i < kid_numm; i++)
		{

		if (fill_len_new[sys->child_index[a][b][i] - 1][c] != 0)
		{
		mul_lIt_real_speed = (fill_len_new[sys->child_index[a][b][i] - 1][c])*Ib[sys->child_index[a][b][i] - 1].real;
		mul_lIt_imag_speed = (fill_len_new[sys->child_index[a][b][i] - 1][c])*Ib[sys->child_index[a][b][i] - 1].imag;
		for (int nn = 0; nn < K_new[track_speed]; nn++)
		{


		UlI_real_s_speed[nn] += UU[U_count_speed] * mul_lIt_real_speed;
		UlI_imag_s_speed[nn] += UU[U_count_speed] * mul_lIt_imag_speed;



		UVlI_real_speed[i] += UU[U_count_speed] * VlI_real_speed[nn];
		UVlI_imag_speed[i] += UU[U_count_speed] * VlI_imag_speed[nn];




		}



		Vb[sys->child_index[a][b][i] - 1].imag += (fill_len_new[sys->child_index[a][b][i] - 1][c])  * UVlI_imag_speed[i];
		Vb[sys->child_index[a][b][i] - 1].real += (fill_len_new[sys->child_index[a][b][i] - 1][c])  * UVlI_real_speed[i];




		}
		++U_count_speed;
		}




		int j;

		for (int j = 0; j < kid_numm_i; j++)
		{

		if (fill_len_new[sys->child_index[a][interact][j] - 1][c] != 0)
		{

		V_count_speed_2 = j + VVV;
		for (int nn = 0; nn < K_new[track_speed]; nn++)
		{
		VUlI_real_s_speed[j] += UlI_real_s_speed[nn] * VV[V_count_speed_2];
		VUlI_imag_s_speed[j] += UlI_imag_s_speed[nn] * VV[V_count_speed_2];
		}

		Vb[sys->child_index[a][interact][j] - 1].real += (fill_len_new[sys->child_index[a][interact][j] - 1][c]) * VUlI_real_s_speed[j];
		Vb[sys->child_index[a][interact][j] - 1].imag += (fill_len_new[sys->child_index[a][interact][j] - 1][c]) * VUlI_imag_s_speed[j];
		}
		++V_count_speed_2;
		}





		}

		++track_speed;


		}

		}


		}


		free(VlI_real_speed);
		free(VlI_imag_speed);
		free(UlI_real_s_speed);
		free(UlI_imag_s_speed);

		free(UVlI_real_speed);
		free(UVlI_imag_speed);

		free(VUlI_real_s_speed);
		free(VUlI_imag_s_speed);

		free(to_do_speed);

		//end speed test BAP

		*/



		//End BAP change

#if 0
		int track = 0;
		double *VlI_real = (double *)malloc(sys->k_max * sizeof(double));
		double *VlI_imag = (double *)malloc(sys->k_max * sizeof(double));
		double *UlI_real_s = (double *)malloc(sys->k_max * sizeof(double));
		double *UlI_imag_s = (double *)malloc(sys->k_max * sizeof(double));

		double *UVlI_real = (double *)malloc(sys->el_max * sizeof(double));
		double *UVlI_imag = (double *)malloc(sys->el_max * sizeof(double));

		double *VUlI_real_s = (double *)malloc(sys->el_max * sizeof(double));
		double *VUlI_imag_s = (double *)malloc(sys->el_max * sizeof(double));

		double mul_lI_real;
		double mul_lI_imag;

		double mul_lIt_real;
		double mul_lIt_imag;

		int to_count;
		int *to_do = (int *)malloc(sys->el_max * sizeof(int));

		clock_t t;
		t = clock();

		double *q_real = (double *)malloc(indsys->num_fils * sizeof(double));
		double *q_imag = (double *)malloc(indsys->num_fils * sizeof(double));

		for (int c = 0; c < 3; c++)
		{
			track = 0;

			for (int i = 0; i < indsys->num_fils; i++)
			{
				q_real[i] = sys->fill[i]->fil->lenvect[c] * Ib[i].real;
				q_imag[i] = sys->fill[i]->fil->lenvect[c] * Ib[i].imag;
			}



			for (int a = 0; a < sys->depth; a++)
			{
				for (int b = 0; b < sys->number_level_groups[a]; b++)
				{
					int kid_numm = sys->kid_num[a][b];
					for (int k = 0; k < sys->group_interaction_number[a][b]; k++)
					{
						//	FILE *matrix = fopen("matrix.txt", "w");




						int interact = sys->interactions[a][b][k];
						int kid_numm_i = sys->kid_num[a][interact];


						memset(VlI_real, 0, sys->k[track] * sizeof(double));
						memset(VlI_imag, 0, sys->k[track] * sizeof(double));
						memset(UlI_real_s, 0, sys->k[track] * sizeof(double));
						memset(UlI_imag_s, 0, sys->k[track] * sizeof(double));

						memset(UVlI_real, 0, kid_numm * sizeof(double));
						memset(UVlI_imag, 0, kid_numm * sizeof(double));

						memset(VUlI_real_s, 0, kid_numm_i * sizeof(double));
						memset(VUlI_imag_s, 0, kid_numm_i * sizeof(double));


						to_count = 0;

						for (int j = 0; j < kid_numm_i; j++)
						{
							/*specific l x specific I*/

							//sum = (sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->lenvect[c] / sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->length)*Ib[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1].real;
							//sum = (sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->lenvect[c] / sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->length)*Ib[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1].imag;
							//sum = (sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->lenvect[c] / sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->length)*Ib[sys->child_index[a + 1][b][i] - 1].real;
							//sum = (sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->lenvect[c] / sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->length)*Ib[sys->child_index[a + 1][b][i] - 1].imag;
							if (sys->fill[sys->child_index[a][interact][j] - 1]->fil->lenvect[c] != 0)
							{
								to_do[to_count] = j;
								++to_count;

								mul_lI_real = q_real[sys->child_index[a][interact][j] - 1];
								mul_lI_imag = q_imag[sys->child_index[a][interact][j] - 1];

								for (int nn = 0; nn < sys->k[track]; nn++)
								{
									//sum += sys->U[track][i][nn] * sys->V[track][nn][j];

									VlI_real[nn] += sys->V[track][nn][j] * mul_lI_real;/*(sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->lenvect[c]) /* sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->length)*Ib[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1].real;*/
									VlI_imag[nn] += sys->V[track][nn][j] * mul_lI_imag; /*(sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->lenvect[c]) /* sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->length)*Ib[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1].imag;*/

																						// place here V*(l*I) 
								}


								//double or = (sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->lenvect[0] / sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->length)*(sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->lenvect[0] / sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->length);
								//or += (sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->lenvect[1] / sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->length)*(sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->lenvect[1] / sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->length);
								//or += (sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->lenvect[2] / sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->length)*(sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->lenvect[2] / sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->length);
								//p[p_index[sys->child_index[a + 1][b][i] -1]] += sum * or * Ib[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1].real;
								//p[p_index[sys->child_index[a + 1][sys->interactions[a][b][k]][j] -1]] += sum * or * Ib[sys->child_index[a + 1][b][i] - 1].real;
								//Vb[sys->child_index[a + 1][b][i] - 1].real += sum * or *Ib[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1].real;
								//Vb[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1].real += sum * or *Ib[sys->child_index[a + 1][b][i] - 1].real;
								//Vb[sys->child_index[a + 1][b][i] - 1].imag += sum * or *Ib[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1].imag;
								//Vb[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1].imag += sum * or *Ib[sys->child_index[a + 1][b][i] - 1].imag;
							}
						}
						//double *UVlI_real = (double *)calloc(sys->kid_num[a + 1][b], sizeof(double));
						//double *UVlI_imag = (double *)calloc(sys->kid_num[a + 1][b], sizeof(double));
						for (int i = 0; i < kid_numm; i++)
						{
							/* after n loop is complete it becomes possible */
							if (sys->fill[sys->child_index[a][b][i] - 1]->fil->lenvect[c] != 0)
							{
								mul_lIt_real = q_real[sys->child_index[a][b][i] - 1];//second
								mul_lIt_imag = q_imag[sys->child_index[a][b][i] - 1];//second
								for (int nn = 0; nn < sys->k[track]; nn++)
								{

									//start second
									UlI_real_s[nn] += sys->U[track][nn][i] * mul_lIt_real;/*(sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->lenvect[c])/* / sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->length)*Ib[sys->child_index[a + 1][b][i] - 1].real;*/
									UlI_imag_s[nn] += sys->U[track][nn][i] * mul_lIt_imag;/* (sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->lenvect[c]) /* sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->length)*Ib[sys->child_index[a + 1][b][i] - 1].imag;*/
																						  //end second


									UVlI_real[i] += sys->U[track][nn][i] * VlI_real[nn];
									UVlI_imag[i] += sys->U[track][nn][i] * VlI_imag[nn];



									//place U x (V x ln x In)
									//sys->U[track][i][nn]
								}

								//lm x (Um x (V x ln x In)/ sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->length) * UVlI_real[sys->child_index[a + 1][b][i] - 1];

								Vb[sys->child_index[a][b][i] - 1].imag += (sys->fill[sys->child_index[a][b][i] - 1]->fil->lenvect[c]) /* sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->length)*/ * UVlI_imag[i] * MUOVER4PI;
								Vb[sys->child_index[a][b][i] - 1].real += (sys->fill[sys->child_index[a][b][i] - 1]->fil->lenvect[c]) /* sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->length)*/ * UVlI_real[i] * MUOVER4PI;




							}
						}

						//free(UVlI_real);
						//free(UVlI_imag);

						//for (int i = 0; i < sys->kid_num[a + 1][b]; i++)
						//{
						//if (sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->lenvect[c] != 0)
						//{

						//}
						//mul_Il_real = (sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->lenvect[c])*Ib[sys->child_index[a + 1][b][i] - 1].real;
						//mul_Il_imag = (sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->lenvect[c])*Ib[sys->child_index[a + 1][b][i] - 1].imag;
						//for (int nn = 0; nn < sys->k[track]; nn++)
						//{
						//	IlU_real_s[nn] += sys->U[track][i][nn] * mul_Il_real;/*(sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->lenvect[c])/* / sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->length)*Ib[sys->child_index[a + 1][b][i] - 1].real;*/
						//IlU_imag_s[nn] += sys->U[track][i][nn] * mul_Il_imag;/* (sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->lenvect[c]) /* sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->length)*Ib[sys->child_index[a + 1][b][i] - 1].imag;*/
						//}
						//}
						//}
						/*double *IlUV_real_s = (double *)calloc(sys->kid_num[a + 1][sys->interactions[a][b][k]], sizeof(double));*/
						//double *IlUV_imag_s = (double *)calloc(sys->kid_num[a + 1][sys->interactions[a][b][k]], sizeof(double));
						//int j;
						//int interacti = sys->interactions[a][b][k];
						//for (int w = 0; w < to_count; w++)
						for (int j = 0; j < sys->kid_num[a][interact]; j++)
						{

							//j = to_do[w];
							if (sys->fill[sys->child_index[a][sys->interactions[a][b][k]][j] - 1]->fil->lenvect[c] != 0)
							{
								//
								for (int nn = 0; nn < sys->k[track]; nn++)
								{
									VUlI_real_s[j] += UlI_real_s[nn] * sys->V[track][nn][j];
									VUlI_imag_s[j] += UlI_imag_s[nn] * sys->V[track][nn][j];
								}

								Vb[sys->child_index[a][interact][j] - 1].real += (sys->fill[sys->child_index[a][interact][j] - 1]->fil->lenvect[c]) /* sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->length)*/* VUlI_real_s[j] * MUOVER4PI;
								Vb[sys->child_index[a][interact][j] - 1].imag += (sys->fill[sys->child_index[a][interact][j] - 1]->fil->lenvect[c]) /* sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->length)*/* VUlI_imag_s[j] * MUOVER4PI;
							}
						}
						//free(IlUV_real_s);
						//free(IlUV_imag_s);



						++track;
					}




				}

			}


		}
		t = clock() - t;
		double time_taken = ((double)t) / CLOCKS_PER_SEC; // in seconds

		free(VlI_real);
		free(VlI_imag);
		free(UlI_real_s);
		free(UlI_imag_s);

		free(UVlI_real);
		free(UVlI_imag);

		free(VUlI_real_s);
		free(VUlI_imag_s);

		free(to_do);
#endif
#endif //end comment bap 
#if MLACA2
		for (chg = chglist; chg != NULL; chg = chg->next) {
			/* add potential due to i direction */
			//Start BAP change
			//Vb[chg->fil->filnumber].real = p[chg->index];
			//End BAP change
			Vb[chg->fil->filnumber].real += p[chg->index] * chg->fil->lenvect[i]
				* MUOVER4PI;
#if OPCNT == ON
			ind_opcnt_mult++;
#endif
		}

		/* do the imaginary part */
		for (chg = chglist; chg != NULL; chg = chg->next) {
			/* fill the pseudo-charge vector */
			q[chg->index] = Ib[chg->fil->filnumber].imag * chg->fil->lenvect[i];
#if OPCNT == ON
			ind_opcnt_mult++;
#endif
		}
		computePsi(sys, q, p, branches, chglist);
#endif
		//Start BAP change
		/*
		int how_many = 0;
		for (i = 0; i < branches; i++)
		p[i] = 0;
		track = 0;
		for (int a = 0; a < sys->depth - 1; a++)
		{
		for (int b = 0; b < sys->number_level_groups[a + 1]; b++)
		{

		for (int k = 0; k < sys->group_interaction_number[a][b]; k++)
		{
		//	FILE *matrix = fopen("matrix.txt", "w");
		for (int i = 0; i < sys->kid_num[a + 1][b]; i++)
		{

		for (int j = 0; j < sys->kid_num[a + 1][sys->interactions[a][b][k]]; j++)
		{
		double sum = 0;
		for (int nn = 0; nn < sys->k[track]; nn++)
		{
		sum += sys->U[track][i][nn] * sys->V[track][nn][j];
		}
		double or = (sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->lenvect[0] / sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->length)*(sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->lenvect[0] / sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->length);
		or += (sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->lenvect[1] / sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->length)*(sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->lenvect[1] / sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->length);
		or += (sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->lenvect[2] / sys->fill[sys->child_index[a + 1][b][i] - 1]->fil->length)*(sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->lenvect[2] / sys->fill[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1]->fil->length);
		//p[p_index[sys->child_index[a + 1][b][i] -1]] += sum * or * Ib[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1].imag;
		//p[p_index[sys->child_index[a + 1][sys->interactions[a][b][k]][j] -1]] += sum * or * Ib[sys->child_index[a + 1][b][i] - 1].imag;
		//Vb[sys->child_index[a + 1][b][i] - 1].imag += sum * or * Ib[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1].imag;
		//Vb[sys->child_index[a + 1][sys->interactions[a][b][k]][j] - 1].imag += sum * or * Ib[sys->child_index[a + 1][b][i] - 1].imag;
		how_many += 2;
		}
		}
		++track;
		}

		}


		*/
		//End BAP change
		//start direct BAP change


#if MLACA //start comment bap 

		int how_many = 0;
		int tally = 0;
		for (int i = 0; i < sys->group_num; i++)//number of groups on lowest level 
		{
			for (int j = 0; j < sys->kid_num[sys->depth - 1][i]; j++)//number of filaments in observation group 
			{
				for (int k = 0; k < sys->kid_num[sys->depth - 1][i]; k++)//number of filaments in observation group 
				{

					Vb[sys->child_index[sys->depth - 1][i][j] - 1].real += sys->direct_inductance[tally] * Ib[sys->child_index[sys->depth - 1][i][k] - 1].real;
					Vb[sys->child_index[sys->depth - 1][i][j] - 1].imag += sys->direct_inductance[tally] * Ib[sys->child_index[sys->depth - 1][i][k] - 1].imag;
					++tally;
					++how_many;
				}
			}
		}


		for (int i = 0; i < sys->group_num; i++)//number of groups on lowest level 
		{
			for (int j = 0; j < sys->number_next_iteration[i]; j++)//near interaction groups
			{
				for (int k = 0; k < sys->kid_num[sys->depth - 1][i]; k++)//number of filaments in observation group 
				{
					for (int l = 0; l < sys->kid_num[sys->depth - 1][sys->next_iteration[i][j]]; l++)//number of filaments in source group
					{

						Vb[sys->child_index[sys->depth - 1][i][k] - 1].imag += sys->direct_inductance[tally] * Ib[sys->child_index[sys->depth - 1][sys->next_iteration[i][j]][l] - 1].imag;// = calcp(sys->fill[sys->child_index[sys->depth][i][k]], sys->fill[sys->child_index[sys->depth][sys->next_iteration[i][j]][l]], NULL);
						Vb[sys->child_index[sys->depth - 1][sys->next_iteration[i][j]][l] - 1].imag += sys->direct_inductance[tally] * Ib[sys->child_index[sys->depth - 1][i][k] - 1].imag;
						Vb[sys->child_index[sys->depth - 1][i][k] - 1].real += sys->direct_inductance[tally] * Ib[sys->child_index[sys->depth - 1][sys->next_iteration[i][j]][l] - 1].real;// = calcp(sys->fill[sys->child_index[sys->depth][i][k]], sys->fill[sys->child_index[sys->depth][sys->next_iteration[i][j]][l]], NULL);
						Vb[sys->child_index[sys->depth - 1][sys->next_iteration[i][j]][l] - 1].real += sys->direct_inductance[tally] * Ib[sys->child_index[sys->depth - 1][i][k] - 1].real;
						++tally;
						how_many += 2;
					}

				}
			}
		}








		//end direct BAP change
#endif //end comment bap 

#if MLACA2
		for (chg = chglist; chg != NULL; chg = chg->next) {
			/* add potential due to i direction */
			//Vb[chg->fil->filnumber].imag = p[chg->index];
			Vb[chg->fil->filnumber].imag += p[chg->index] * chg->fil->lenvect[i]
				* MUOVER4PI;
#if OPCNT == ON
			ind_opcnt_mult++;
#endif
		}

	}//link 3
#endif
#if MLACA2
	 /* do the direct parts */
	sys->DirectEval = TRUE;
	/* do the real part of the Direct part */
	for (i = 1; i <= branches; i++)
		p[i] = 0;


	//Start BAP change

	//End BAP change


	for (chg = chglist; chg != NULL; chg = chg->next)
		/* fill the pseudo-charge vector */
		q[chg->index] = Ib[chg->fil->filnumber].real;

	/* starttimer; */
	mulDirect(sys);
	mulEval(sys);
	/* stoptimer; */
	dirtime += dtime;

	for (chg = chglist; chg != NULL; chg = chg->next) {
		/* add potential due to i direction */
		Vb[chg->fil->filnumber].real += p[chg->index];//Bap comment out 
	}

	/* do the imaginary part of the Direct part */
	for (i = 1; i <= branches; i++)
		p[i] = 0;
	for (chg = chglist; chg != NULL; chg = chg->next)
		/* fill the pseudo-charge vector */
		q[chg->index] = Ib[chg->fil->filnumber].imag;

	/* starttimer; */
	mulDirect(sys);
	mulEval(sys);
	/* stoptimer; */
	dirtime += dtime;

	for (chg = chglist; chg != NULL; chg = chg->next) {
		/* add potential due to i direction */
		Vb[chg->fil->filnumber].imag += p[chg->index];//Bap comment out 
	}
#endif
	/* do Vs = M*Vb*jw */
	for (i = 0; i < size; i++) {
		Vs[i] = CXZERO;
		for (mtemp = Mlist[i]; mtemp != NULL; mtemp = mtemp->mnext)
			if (mtemp->sign == 1)
				cx_add(Vs[i], Vs[i], Vb[mtemp->filindex]);
			else
				cx_sub(Vs[i], Vs[i], Vb[mtemp->filindex]);

		/* multiply by jw */
		rtemp = -Vs[i].imag * w;
		Vs[i].imag = Vs[i].real * w;
		Vs[i].real = rtemp;
	}

	/* add in M*R*Mt*Im = M*R*Ib */
	for (i = 0; i < size; i++) {
		for (mtemp = Mlist[i]; mtemp != NULL; mtemp = mtemp->mnext) {
			cx_scalar_mult(temp, mtemp->sign * R[mtemp->filindex],
				Ib[mtemp->filindex]);
			cx_add(Vs[i], Vs[i], temp);
#if OPCNT == ON
			ind_opcnt_mult += 2;
			ind_opcnt_real += 2;
#endif
		}
	}

#if OPCNT == ON
	printf("Inductance (mesh to branch) mults: %d\n", ind_opcnt_mult);
	printf("Just doing MRMtIm: %d\n", ind_opcnt_real);
	printops();
	exit(0);
#endif

#ifdef NODEBUG
	/* for debugging, compare to direct Vb = ZM Ib */
	realmatCXvec(Vdirect, indsys->Z, Ib, branches);
	maxdiff = 0;
	maxindx = 0;
	for (i = 0; i < branches; i++) {
		if (cx_abs(Vb[i]) > 1e-23) {
			cx_sub(temp, Vdirect[i], Vb[i]);
			pdiff = cx_abs(temp) / cx_abs(Vb[i]);
		}
		else
			pdiff = cx_abs(Vb[i]);

		if (pdiff > maxdiff) {
			maxdiff = pdiff;
			maxindx = i;
		}
	}
	if (maxdiff < .3)
		printf("maxdiff: %g  Vb[%d]=%g  Vdirect[%d]=%g\n",
			maxdiff, maxindx, cx_abs(Vb[maxindx]), maxindx, cx_abs(Vdirect[maxindx]));
	else
		printf("***maxdiff: %g  Vb[%d]=%g  Vdirect[%d]=%g***\n",
			maxdiff, maxindx, cx_abs(Vb[maxindx]), maxindx, cx_abs(Vdirect[maxindx]));

#endif
}

void realmatCXvec(CX *y, double **A, CX *x, int size) {
	int i, j;
	CX temp;

	for (i = 0; i < size; i++) {
		y[i] = CXZERO;
		for (j = 0; j < size; j++) {
			cx_scalar_mult(temp, A[i][j], x[j]);
			cx_add(y[i], y[i], temp);
		}
	}
}

/* this function fixes Eval matrices which are computed directly */
/* This is necessary since direct mutual terms are not componentwise,
but the multipole routines are called once for each component direction.
Basically, componentwise multiplication will cause the elements
to be multiplied by the dot product of the fil->lenvect vectors of
the two filaments.  This will divide that product out.  Also, MUOVER4PI
must also be divided out
*/

void fixEvalDirect(charge **qchgs, int numqchgs, int *is_dummy, charge **pchgs,
	int numpchgs, double **mat) {
	int i, j, k;
	double dotprod, magi, magj;
	double *lenvecti, *lenvectj;
	static double eps = EPS;

	for (i = 0; i < numpchgs; i++) {
		lenvecti = pchgs[i]->fil->lenvect;
		magi = 0;
		for (k = 0; k < 3; k++)
			magi += lenvecti[k] * lenvecti[k];
		for (j = 0; j < numqchgs; j++) {
			lenvectj = qchgs[j]->fil->lenvect;
			magj = dotprod = 0;
			for (k = 0; k < 3; k++) {
				magj += lenvectj[k] * lenvectj[k];
				dotprod += lenvecti[k] * lenvectj[k];
			}
			if (fabs(dotprod) / sqrt(magi * magj) > EPS) /* filaments aren't perpendicular */
				mat[i][j] = mat[i][j] / (dotprod * MUOVER4PI);
			else { /* if they are, mat[i][j] == 0.0, hopefully */
				if (mat[i][j] != 0.0)
					printf(
						"Warning: dot product = %lg < EPS, but mat[i][j] = %lg\n",
						dotprod, mat[i][j]);
			}
		}
	}
}
