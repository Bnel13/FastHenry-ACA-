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



/* # ***** sort to /src/main
 # ***** */
//#include "mulGlobal.h"
#include "induct.h"//BAPN change

int *localcnt, *multicnt, *evalcnt; /* counts of builds done by level */

int **Q2Mcnt, **Q2Lcnt, **Q2Pcnt, **L2Lcnt; /* counts of xformation mats */
int **M2Mcnt, **M2Lcnt, **M2Pcnt, **L2Pcnt, **Q2PDcnt;

/* SRW */
void mulMatDirect(ssystem*);
void bdmulMatPrecond(ssystem*);
void olmulMatPrecond(ssystem*);
void find_flux_density_row(double**, double**, int, int, int, int, int,
		charge**, charge**, int*, int*);
void mulMatUp(ssystem*);
void mulMatEval(ssystem*);
void mulMatDown(ssystem*);


long double mutualfill33(FILAMENT *fil1, FILAMENT *fil2) {
	long double inductance = 0;
	long double volume1 = 0;
	long double volume2 = 0;
	long double distance = 0;
	long double disx = 0.0;
	long double disy = 0;
	long double disz = 0;
	long double permeability = 0;
	long double combo = 0;
	disx = (fil1->x[0] + fil1->x[1]) / 2 - (fil2->x[0] + fil2->x[1]) / 2;// ^ 2;
	disy = (fil1->y[0] + fil1->y[1]) / 2 - (fil2->y[0] + fil2->y[1]) / 2;// ^ 2;
	disz = (fil1->z[0] + fil1->z[1]) / 2 - (fil2->z[0] + fil2->z[1]) / 2;// ^ 2);
	combo = disx * disx + disy * disy + disz * disz;
	long double a1 = 0;
	long double a2 = 0;
	a1 = fil1->area;
	a2 = fil2->area;
	volume1 = fil1->area*fil1->length;
	volume2 = fil2->area*fil2->length;

	distance = sqrt(combo);
	inductance = (MU0 / (4 * PI *a1*a2))*((volume1*volume2) / (distance));
	return inductance;
}



void MLACA(int depth, int *number_of_groups_level, int **number_of_interactions_group, int ***interactions, int **num_elements, int ***elements, charge **filaments, double **total,ssystem *ssys)
{
	//ssys->number_groups_level = (int *)calloc(depth,sizeof(int));
	int hshshshshss = 0;
	//int nummm = 0;
	/*
	double **total = (double **)malloc(1643*sizeof(double*));
	for (int i = 0; i < 1643; i++)
	{
	total[i] = (double *)malloc(1643*sizeof(double));
	}
	*/
	int count = 0;
	for (int i = 0; i < depth - 1; i++)
	{
		for (int j = 0; j < number_of_groups_level[i + 1]; j++)
		{
			for (int k = 0; k < number_of_interactions_group[i][j]; k++)
			{
				++count;
			}
		}
	}
	//ssys->k = (int*)malloc(count * sizeof(int));
	long double ***UU = (long double ***)calloc(count, sizeof(long double **));
	//ssys->U = (long double ***)calloc(count, sizeof(long double **));
	long double ***VV = (long double ***)calloc(count, sizeof(long double **));
	//ssys->V = (long double ***)calloc(count, sizeof(long double **));
	/*
	int count1 = 0;
	for (int i = 0; i < depth - 1; i++)
	{
		for (int j = 0; j < number_of_groups_level[i + 1]; j++)
		{
			for (int k = 0; k < number_of_interactions_group[i][j]; k++)
			{
				UU[count1] = (long double **)calloc(num_elements[i + 1][j], sizeof(long double*));
				VV[count1] = (long double **)calloc(num_elements[i + 1][interactions[i][j][k]], sizeof(long double*));//k

				for (int q = 0; q < num_elements[i + 1][j]; q++)
				{
					UU[count1][q] = (long double *)calloc(num_elements[i + 1][j], sizeof(long double));//k
				}
				for (int p = 0; p < num_elements[i + 1][interactions[i][j][k]]; p++)//k
				{
					VV[count1][p] = (long double *)calloc(num_elements[i + 1][interactions[i][j][k]], sizeof(long double));
				}
				++count1;
			}
		}
	}
	*///End working version
	//for (int i = 0; i<1; i++)
	//U[i] = (double *)malloc( sizeof(double));
	//for (int i = 0; i<1; i++)
	//V[i] = (double *)malloc(sizeof(double));
	int *KK = (int *)malloc(count * sizeof(int));
	int count2 = 0;
	//unsigned long int memory = 0;
	int memory = 0;
	//ssys->number_of_interactions_groups = (int **)malloc(depth * sizeof(int*));
	//ssys->interact = (int***)calloc(depth , sizeof(int**));
	for (int i = 0; i < depth - 1; i++)
	{
		//ssys->interact[i] = (int**)calloc(number_of_groups_level[i + 1], sizeof(int*));
		//ssys->number_of_interactions_groups[i] = (int *)malloc(number_of_groups_level[i]* sizeof(int));
		for (int j = 0; j < number_of_groups_level[i + 1]; j++)
		{
			//ssys->interact[i][j] = (int**)calloc(number_of_groups_level[i + 1], sizeof(int*));
			//++ssys->number_groups_level[i + 1];
			//ssys->number_of_interactions_groups[i][j] = number_of_interactions_group[i][j];
			for (int k = 0; k < number_of_interactions_group[i][j]; k++)
			{
				//ssys->interact[i][j][k] = interactions[i][j][k];
				//Start memory allocation change
				UU[count2] = (long double **)calloc(num_elements[i + 1][j], sizeof(long double*));
				//ssys->U[count2] = (long double **)calloc(num_elements[i + 1][j], sizeof(long double*));
				VV[count2] = (long double **)calloc(num_elements[i + 1][interactions[i][j][k]], sizeof(long double*));//k
				//ssys->V[count2] = (long double **)calloc(num_elements[i + 1][interactions[i][j][k]], sizeof(long double*));//k
				for (int q = 0; q < num_elements[i + 1][j]; q++)
				{
					UU[count2][q] = (long double *)calloc(num_elements[i + 1][j], sizeof(long double));//k
					//ssys->U[count2][q] = (long double *)calloc(num_elements[i + 1][j], sizeof(long double));//k
				}
				for (int p = 0; p < num_elements[i + 1][interactions[i][j][k]]; p++)//k
				{
					VV[count2][p] = (long double *)calloc(num_elements[i + 1][interactions[i][j][k]], sizeof(long double));
					//ssys->V[count2][p] = (long double *)calloc(num_elements[i + 1][interactions[i][j][k]], sizeof(long double));
				}
				//end memory allocation change


				//fprintf(matrix, "%d\t", memory);
				/*
				double **U = (double **)malloc(num_elements[i + 1][j] * sizeof(double *));
				for (int i = 0; i<num_elements[i + 1][j]; i++)
				U[i] = (double *)malloc(num_elements[i + 1][interactions[i][j][k]] / 2 * sizeof(double));

				double **V = (double **)malloc(num_elements[i + 1][j] / 2 * sizeof(double *));
				for (int i = 0; i<num_elements[i + 1][j] / 2; i++)
				V[i] = (double *)malloc(num_elements[i + 1][interactions[i][j][k]] * sizeof(double));
				*/
				
				
				//memory = ACA_new(num_elements[i + 1][j], num_elements[i + 1][interactions[i][j][k]], elements[i + 1][j], elements[i + 1][interactions[i][j][k]], filaments,UU[count2],VV[count2]);
				KK[count2] = memory;
				//ssys->k[count2] = memory;
				//total[elements[a + 1][b][i] - 1][elements[a + 1][interactions[a][b][k]][j] - 1]
				/*
				FILE *matrix = fopen("matrix.txt", "w");
				*/
				/*
				for (int q = 0; q < num_elements[i + 1][j]; q++)
				{
					for (int p = 0; p < num_elements[i + 1][interactions[i][j][k]]; p++)
					{
						double sum = 0;
						for (int i = 0; i < memory; i++)
						{
							sum += UU[count2][q][i] * VV[count2][i][p];
						}
						//total[elements[i + 1][j][q] - 1][elements[i + 1][interactions[i][j][k]][p] - 1] = sum;
						//fprintf(matrix, "%e\t", calcp(filaments[elements[a + 1][b][i] - 1], filaments[elements[a + 1][interactions[a][b][k]][j] - 1], NULL) - sum);
					}
					//fprintf(matrix, "\n");
				}
				*/
				/*
				fclose(matrix);
				*/


				for (int q = 0; q < num_elements[i + 1][j]; q++)
				{
					UU[count2][q] = (double *)realloc(UU[count2][q], memory * sizeof(double));
					//ssys->U[count2][q] = (double *)realloc(UU[count2][q], memory * sizeof(double));
				}
				for (int p = memory; p < num_elements[i + 1][interactions[i][j][k]]; p++)//k
				{
					free(VV[count2][p]);
					//free(ssys->V[count2][p]);
				}

				++count2;

				/*
				UU = (double ***)realloc(UU, ka * sizeof(double**));


				UU[ka - 1] = (double **)malloc(num_elements[i + 1][j] * sizeof(double *));

				for (int j = 0; j < num_elements[i + 1][j]; j++) {

				UU[ka - 1][j] = (double *)malloc(kk * sizeof(double));
				memcpy(&(UU[ka - 1][j][0]), &(U[j][0]), kk * sizeof(double));
				}

				VV = (double ***)realloc(VV, ka * sizeof(double**));
				VV[ka - 1] = (double **)malloc(kk * sizeof(double *));

				for (int j = 0; j < kk; j++) {

				VV[ka - 1][j] = (double *)malloc(num_elements[i + 1][interactions[i][j][k]] * sizeof(double));
				memcpy(&(VV[ka - 1][j][0]), &(V[j][0]), num_elements[i + 1][interactions[i][j][k]] * sizeof(double));

				}
				*/
				//memory += num_elements[i+1][j] * num_elements[i+1][interactions[i][j][k]];
			}

		}

	}



	//for (a = 0, nextc = sys->directlist; nextc != NULL; nextc =
	//nextc->dnext) 
	int see = 0;
	int track = 0;
	for (int a = 0; a < depth - 1; a++)
	{
		for (int b = 0; b < number_of_groups_level[a + 1]; b++)
		{
			for (int k = 0; k < number_of_interactions_group[a][b]; k++)
			{
				//	FILE *matrix = fopen("matrix.txt", "w");
				for (int i = 0; i < num_elements[a + 1][b]; i++)
				{

					for (int j = 0; j < num_elements[a + 1][interactions[a][b][k]]; j++)
					{
						long double sum = 0;
						for (int nn = 0; nn < KK[track]; nn++)
						{
							sum += UU[track][i][nn] * VV[track][nn][j];
						}
						long double or = (filaments[elements[a + 1][b][i] - 1]->fil->lenvect[0] / filaments[elements[a + 1][b][i] - 1]->fil->length)*(filaments[elements[a + 1][interactions[a][b][k]][j] - 1]->fil->lenvect[0] / filaments[elements[a + 1][interactions[a][b][k]][j] - 1]->fil->length);
						or += (filaments[elements[a + 1][b][i] - 1]->fil->lenvect[1] / filaments[elements[a + 1][b][i] - 1]->fil->length)*(filaments[elements[a + 1][interactions[a][b][k]][j] - 1]->fil->lenvect[1] / filaments[elements[a + 1][interactions[a][b][k]][j] - 1]->fil->length);
						or += (filaments[elements[a + 1][b][i] - 1]->fil->lenvect[2] / filaments[elements[a + 1][b][i] - 1]->fil->length)*(filaments[elements[a + 1][interactions[a][b][k]][j] - 1]->fil->lenvect[2] / filaments[elements[a + 1][interactions[a][b][k]][j] - 1]->fil->length);
						//if(((calcp(filaments[elements[a + 1][b][i] - 1], filaments[elements[a + 1][interactions[a][b][k]][j] - 1], NULL) - sum* or)/ calcp(filaments[elements[a + 1][b][i] - 1], filaments[elements[a + 1][interactions[a][b][k]][j] - 1], NULL))>0.1)
						total[elements[a + 1][b][i] - 1][elements[a + 1][interactions[a][b][k]][j] - 1] = total[elements[a + 1][interactions[a][b][k]][j] - 1][elements[a + 1][b][i] - 1] = 
							sum * or;
							//(mutualfill33(filaments[elements[a + 1][b][i] - 1]->fil, filaments[elements[a + 1][interactions[a][b][k]][j] - 1]->fil))* or;//	calcp(filaments[elements[a + 1][b][i] - 1], filaments[elements[a + 1][interactions[a][b][k]][j] - 1], NULL);//sum* or ;
						//mutualfill3(filaments[elements[a + 1][b][i] - 1], filaments[elements[a + 1][interactions[a][b][k]][j] - 1]);
						//total[elements[a + 1][interactions[a][b][k]][j] - 1][elements[a + 1][b][i] - 1] = //sum* or ;
						//fprintf(matrix, "%e\t", calcp(filaments[elements[a + 1][b][i] - 1], filaments[elements[a + 1][interactions[a][b][k]][j] - 1], NULL) - sum*or);
						see += 2;
						//++nummm;
						if (((mutualfill33(filaments[elements[a + 1][b][i] - 1]->fil, filaments[elements[a + 1][interactions[a][b][k]][j] - 1]->fil)) - sum)/ mutualfill33(filaments[elements[a + 1][b][i] - 1]->fil, filaments[elements[a + 1][interactions[a][b][k]][j] - 1]->fil) > 0.0001)
							++hshshshshss;
					}
					//fprintf(matrix, "\n");
				}
				//	fclose(matrix);
				++track;
			}

		}


	}
	//ssys->count = track;
	/*
	FILE *matrix = fopen("matrix.txt", "w");
	for (int i = 0; i < 1643; i++)
	{
	for (int j = 0; j < 1643; j++)
	{
	fprintf(matrix, "%e\t", total[i][j]);
	}
	fprintf(matrix, "\n");
	}
	fclose(matrix);
	*/
}
//End BAP

//Start BAP
void branch(int ***Group_kids, double ***Center, int *paren_index, int **kids_index, int level, int group_number, double distance,
	int previous_unknowns, int *unknown_back, int ***interaction, int parent_index, int depth, int **interaction_number, int *memoryyy, int **num_elements, int **still_to_do, int *number_still_to_do)
{
	int unknowns = 0;
	int *unknown = (int *)malloc(paren_index[level] * sizeof(int));
	for (int i = 0; i < paren_index[level]; i++)
	{
		unknown[i] = 0;
	}
	int unknown_nexts = 0;
	int *unknown_next = (int *)malloc(paren_index[level + 1] * sizeof(int));
	for (int i = 0; i < paren_index[level + 1]; i++)
	{
		unknown_next[i] = 0;
	}
	int index = 0;
	//write code to place the unknowns from the same group that are still able to exploit symetry
	//int max_kids = kids_index[level - 1][parent_index];
	//int j = 0;
	//while (group_number != Group_kids[level - 1][parent_index][max_kids - j - 1])//goes down symetric branches 
	//{
	//	unknown[unknowns] = Group_kids[level - 1][parent_index][max_kids - j - 1];
	//	++j;
	//	++unknowns;
	//}
	interaction_number[level - 1][group_number] = 0;
	for (int i = 0; i < previous_unknowns; i++)//Go through all possible interactions on the level
	{
		int source = unknown_back[i];
		//for (int j = 0; j < kids_index[level - 1][uncover_unknown]; j++)
		//{
		//int source = Group_kids[level - 1][uncover_unknown][j];
		double distance_x = fabs(Center[level][group_number][0] - Center[level][source][0]);
		double distance_y = fabs(Center[level][group_number][1] - Center[level][source][1]);
		double distance_z = fabs(Center[level][group_number][2] - Center[level][source][2]);
		if (distance_x  >  distance * 2 || distance_y > distance * 2 || distance_z > distance * 2)
		{
			//this is where the interactions are determained and stored in no particular order
			interaction[level - 1][group_number][index] = source;
			++index;
			++interaction_number[level - 1][group_number];
		}
		else if (depth == level + 1)
		{
			still_to_do[group_number][number_still_to_do[group_number]] = source;
			++number_still_to_do[group_number];
		}
		else
		{
			//This is where all the pices that need to still be investigated on the next level should be stored.
			unknown[unknowns] = source;
			++unknowns;
		}

		//}
	}
	//function compairing start_observation with all the possible interactions.
	for (int i = 0; i < unknowns; i++)
	{
		for (int j = 0; j < kids_index[level][unknown[i]]; j++)
		{
			unknown_next[unknown_nexts] = Group_kids[level][unknown[i]][j];
			++unknown_nexts;
		}
	}

	for (int i = 0; i < kids_index[level][group_number] - 1; i++)//call each time for every child go down all posible brances
	{
		unknown_next[unknown_nexts] = Group_kids[level][group_number][kids_index[level][group_number] - i - 1];;
		++unknown_nexts;
	}

	if (depth == level + 1)
	{
		for (int i = 0; i < unknown_nexts; i++)
		{
			//memoryyy[1] += num_elements[5][group_number] * num_elements[5][unknown_next[i]];
		}


		//memoryyy[0] += num_elements[5][group_number] * num_elements[5][group_number];

	}

	if (depth > level - 1)
	{
		for (int i = 0; i < kids_index[level][group_number]; i++)//call each time for every child go down all posible brances
		{
			branch(Group_kids, Center, paren_index, kids_index, level + 1, Group_kids[level][group_number][i], distance / 2,
				unknown_nexts - i, unknown_next, interaction, group_number, depth, interaction_number, memoryyy, num_elements, still_to_do, number_still_to_do);
		}

	}

}
//End Bap



//Start Bap
void top_levels(int ***Group_kids, int *paren_index, double ***Center, int **kids_index, double distance, int depth, charge **fill, int **number_els, int ***els, int number_groups, int number_fils, double **total, ssystem *syys)
{
	int **group_interaction_number;
	int mem = 0;
	int mem2 = 0;
	int *memoryy = (int*)malloc(2 * sizeof(int));
	memoryy[0] = 0;
	memoryy[1] = 0;
	group_interaction_number = (int **)malloc(depth * sizeof(int *));

	int **extra = (int **)malloc(number_groups * sizeof(int *));
	for (int i = 0; i < number_groups; i++)
	{
		extra[i] = (int *)malloc(100 * sizeof(int));
	}

	int *number_extra = (int *)malloc(number_groups * sizeof(int));
	for (int i = 0; i < number_groups; i++)
	{
		number_extra[i] = 0;
	}

	for (int i = 0; i < depth; i++)
	{
		group_interaction_number[i] = (int *)malloc(number_groups * sizeof(int));
	}

	int observation;
	distance *= (depth - 1);
	int unknowns = 0;
	int x = 0;
	int y = 0;
	int *unknown = (int *)malloc(paren_index[1] * sizeof(int));
	int ***interaction = (int ***)malloc(depth * sizeof(int**));
	for (int i = 0; i < depth - 1; i++)
	{
		interaction[i] = (int **)malloc(paren_index[i + 1] * sizeof(int *));
		for (int j = 0; j < paren_index[i + 1]; j++)
		{
			interaction[i][j] = (int *)malloc((paren_index[i + 1] - 1) * sizeof(int));
		}

	}

	for (int i = 0; i < paren_index[1]; i++)
	{
		unknowns = 0;
		if (y == kids_index[0][x])
		{
			++x;
			y = 0;
		}
		observation = Group_kids[0][x][y];
		++y;
		for (int j = observation + 1; j < paren_index[1]; j++)
		{
			unknown[unknowns] = j;
			++unknowns;
			mem += number_els[1][observation] * number_els[1][j];
		}
		mem2 += number_els[1][observation] * number_els[1][observation];
		branch(Group_kids, Center, paren_index, kids_index, 1, observation, distance*2, unknowns, unknown, interaction, x, depth, group_interaction_number, memoryy, number_els, extra, number_extra);
	}

	int track = 0;
	for (int i = 0; i < number_groups; i++)
	{
		for (int j = 0; j < number_extra[i]; j++)
		{
			//memoryy[1] += number_els[6][i] * number_els[6][extra[i][j]];
			for (int y = 0; y < number_els[depth - 1][i]; y++)
			{
				for (int t = 0; t < number_els[depth - 1][extra[i][j]]; t++)
				{
					//total[elements[a + 1][b][i] - 1][elements[a + 1][interactions[a][b][k]][j] - 1]
					total[els[depth - 1][i][y] - 1][els[depth - 1][extra[i][j]][t] - 1] = calcp(fill[els[depth - 1][i][y] - 1], fill[els[depth - 1][extra[i][j]][t] - 1], NULL);
					total[els[depth - 1][extra[i][j]][t] - 1][els[depth - 1][i][y] - 1] = calcp(fill[els[depth - 1][i][y] - 1], fill[els[depth - 1][extra[i][j]][t] - 1], NULL);
					track += 2;
				}
			}

		}

	}
	
	for (int i = 0; i < number_groups; i++)
	{
		for (int y = 0; y < number_els[depth - 1][i]; y++)
		{
			for (int o = 0; o < number_els[depth - 1][i]; o++)
			{
				total[els[depth - 1][i][y] - 1][els[depth - 1][i][o] - 1] = calcp(fill[els[depth - 1][i][y] - 1], fill[els[depth - 1][i][o] - 1], NULL);
				++track;
			}
		}
	}
	/*
	for (int i = 0; i < 1643; i++)
	{
		for (int y = 0; y < 1643; y++)
		{	
				total[i][y] = calcp(fill[i], fill[y], NULL);
		}
	}
	*/
	//FILE *matrix = fopen("matrix.txt", "w");
	//int a = 0, b = 0;
	//fprintf(matrix, "%d\t", level_count[4]);
	//fprintf(matrix, "%d\t", level_count[4]);

	//for (a = 0, nextc = sys->directlist; nextc != NULL; nextc =
	//nextc->dnext) 

	//for (a = 0; a < 300; a++)
	//{

	//	for (b = 0; b < 2; b++)
	//	{
	//	fprintf(matrix, "%d\t", memoryy[b]);
	//	fprintf(matrix, "\n");
	//	}
	//	fprintf(matrix, "\n");
	//	}

	//fclose(matrix);
	
	MLACA(depth, paren_index, group_interaction_number, interaction, number_els, els, fill,total,syys);
	/*
	FILE *matrix = fopen("matrix.txt", "w");
	for (int i = 0; i < 1643; i++)
	{
	for (int j = 0; j < 1643; j++)
	{
	fprintf(matrix, "%e\t", total[i][j]- calcp(fill[i], fill[j], NULL));
	}
	fprintf(matrix, "\n");
	}
	fclose(matrix);
	*/
}
//end Bap




/*
 MulMatDirect creates the matrices for the piece of the problem that is done
 directly exactly.
 */
void mulMatDirect(ssystem *sys) {
	//sys->matrix = (double**)malloc(23226*sizeof(double*));
	//for (int i = 0; i < 23226; i++)
	//{
		//sys->matrix[i] = (double*)malloc(23226 * sizeof(double));
	//}
	int memoryyy = 0;
	int memory_off = 0;
	cube *nextc, *nextnbr;



	







	
	int i, nummats, **temp = NULL;
	extern double lutime, dirtime;

#if DIRSOL == ON || EXPGCR == ON
	extern double *trimat, *sqrmat; /* flattened triangular, square matrices */
	extern int up_size, eval_size;
	extern int *real_index; /* for map btwn condensed/expanded vectors */
#endif

	/* First count the number of matrices to be done directly. */
	for (nextc = sys->directlist; nextc != NULL; nextc = nextc->dnext) {
		for (nummats = 1, i = 0; i < nextc->numnbrs; i++) {
			nextnbr = nextc->nbrs[i];
			ASSERT(nextnbr->upnumvects > 0);
			nummats++;
		}

		/* Allocate space for the vects and mats. */
		nextc->directnumvects = nummats;
		if (nummats > 0) {
			CALLOC(nextc->directq, nummats, double*, ON, AMSC);
			CALLOC(temp, nummats, int*, ON, AMSC);
			CALLOC(nextc->directnumeles, nummats, int, ON, AMSC);
			CALLOC(nextc->directmats, nummats, double**, ON, AMSC);
			/*      CALLOC(nextc->precondmats, nummats, double**, ON, AMSC); */
		}

		/* initialize the pointer from this cube to its part of dummy vector
		 - save the self part found in indexkid() */
		temp[0] = nextc->nbr_is_dummy[0];
		nextc->nbr_is_dummy = temp;
	}

	/* Now place in the matrices. */
	for (nextc = sys->directlist; nextc != NULL; nextc = nextc->dnext) {
		nextc->directq[0] = nextc->upvects[0];
		nextc->directnumeles[0] = nextc->upnumeles[0];

		/* starttimer; */
#if DIRSOL == ON || EXPGCR == ON
		if(nextc == sys->directlist) {
			if(eval_size < MAXSIZ) {
				fprintf(stderr,
						"mulMatDirect: non-block direct methods not supported\n");
				exit(1);
				/* if this is going to work, need a special, condensing Q2P
				 as well as some way to use it in the framework of the GCR loop */
				nextc->directmats[0] = Q2P(nextc->chgs, eval_size,
						nextc->nbr_is_dummy[0], nextc->chgs,
						eval_size, TRUE);
			}
			else blkQ2Pfull(sys->directlist, up_size, eval_size,
					&trimat, &sqrmat, &real_index, sys->is_dummy);
		}
		else nextc->directmats[0]
		= Q2PDiag(nextc->chgs, nextc->upnumeles[0], nextc->nbr_is_dummy[0],
				TRUE);
#else
		nextc->directmats[0] = Q2PDiag(nextc->chgs, nextc->upnumeles[0],
				nextc->nbr_is_dummy[0],
				TRUE);
		memoryyy += nextc->upnumeles[0] * nextc->upnumeles[0];
		/*
		 nextc->precondmats[0]
		 = Q2PDiag(nextc->chgs, nextc->upnumeles[0], nextc->nbr_is_dummy[0],
		 FALSE);
		 */
		/*dumpMatCor(nextc->directmats[0], (double *)NULL, nextc->upnumeles[0]);*/

#endif
		/* stoptimer; */
		dirtime += dtime;

#if DSQ2PD == ON
		dumpQ2PDiag(nextc);
#endif

#if DMTCNT == ON
		Q2PDcnt[nextc->level][nextc->level]++;
#endif

#if DIRSOL == ON
		/* transform A into LU */
		if(eval_size > MAXSIZ) {
			blkLUdecomp(sqrmat, trimat, up_size);
		}
		else if(nextc == sys->directlist) {
			/* starttimer; */
			nextc->directlu = ludecomp(nextc->directmats[0], eval_size, TRUE);
			/* stoptimer; */
			lutime += dtime;
		}
#endif

		/* starttimer; */
		for (nummats = 1, i = 0; i < nextc->numnbrs; i++) {
			nextnbr = nextc->nbrs[i];
			ASSERT(nextnbr->upnumvects > 0);
			nextc->directq[nummats] = nextnbr->upvects[0];
			nextc->nbr_is_dummy[nummats] = nextnbr->nbr_is_dummy[0];
			nextc->directnumeles[nummats] = nextnbr->upnumeles[0];
			nextc->directmats[nummats] = Q2P(nextnbr->chgs,
					nextnbr->upnumeles[0], nextnbr->nbr_is_dummy[0],
					nextc->chgs, nextc->upnumeles[0],
					TRUE);
			nummats++;
			memory_off += nextnbr->upnumeles[0] * nextc->upnumeles[0];
			/*
			 nextc->precondmats[nummats++] = Q2P(nextnbr->chgs,
			 nextnbr->upnumeles[0],
			 nextnbr->nbr_is_dummy[0],
			 nextc->chgs, nextc->upnumeles[0],
			 FALSE);
			 */
#if DMTCNT == ON
			Q2Pcnt[nextc->level][nextnbr->level]++;
#endif
		}
		/* stoptimer; */
		dirtime += dtime;
	}
	/*
	FILE * fpp;

	fpp = fopen("file_mlfma_direct.txt", "w+");
	fprintf(fpp, "  MLFMA direct memory: %d\n", memoryyy * 8);
	fclose(fpp);
	*/
	/*
		FILE * fppp;

	fppp = fopen("file_mlfma_direct_off.txt", "w+");
	fprintf(fppp, "  MLFMA direct memory: %d\n", memory_off * 8);
	fclose(fppp);
	*/
}

/*
 MulMatPrecond creates the preconditioner matrix
 */
void bdmulMatPrecond(ssystem *sys) {
	cube *nc, *kid, *kidnbr;
	double **mat, **nbrmat;
	int i, j, k, l, kidi;
	int kidsize, nbrsize, size, row, col, first, offset;
	double factor;
	charge *pc;
	surface *surf;

	printf("This Preconditioner is not used in FastHenry\n");
	exit(1);

	for (nc = sys->precondlist; nc != NULL; nc = nc->pnext) {

		/* find total number of charges in cube to dimension P. */
		for (size = 0, i = 0; i < nc->numkids; i++) {
			kid = nc->kids[i];
			if (kid != NULL) {
				ASSERT(kid->level == sys->depth);
				size += kid->directnumeles[0]; /* Equals number of charges. */
			}
		}

		/* allocate and zero a preconditioner matrix. */
		MALLOC(mat, size, double*, ON, AMSC);
		for (i = 0; i < size; i++) {
			MALLOC(mat[i], size, double, ON, AMSC);
		}
		for (i = 0; i < size; i++) {
			for (j = 0; j < size; j++) {
				mat[i][j] = 0.0;
			}
		}

		/* Chase through the kids to place in potential coeffs. */
		for (first = TRUE, row = 0, kidi = 0; kidi < nc->numkids; kidi++) {
			kid = nc->kids[kidi];
			if (kid != NULL) {
				/* Exploit the hierarchical charge numbering to get precond vector. */
				if (first == TRUE) {
					first = FALSE;
					nc->prevectq = kid->directq[0];
					nc->prevectp = kid->eval;
				}
				/* Get the diagonal block of P^{-1}. */
				kidsize = kid->directnumeles[0];
				for (k = kidsize - 1; k >= 0; k--) {
					for (l = kidsize - 1; l >= 0; l--) {
						mat[row + k][row + l] = kid->directmats[0][k][l];
					}
				}
				/* Get the off-diagonals of P^{-1}. */
				for (col = 0, i = 0; i < nc->numkids; i++) {
					kidnbr = nc->kids[i];
					if (kidnbr != NULL) {
						if (kidnbr != kid) {
							/* Chase thru list of nbrs to get matrix associated with this
							 kidnbr.  Note, this is because the kid list and nbr matrix
							 list are in different orders, could be fixed. */
							for (j = kid->numnbrs - 1; j >= 0; j--) {
								if (kidnbr == kid->nbrs[j]) {
									nbrmat = kid->directmats[j + 1];
									nbrsize = kidnbr->directnumeles[0];
									for (k = kidsize - 1; k >= 0; k--) {
										for (l = nbrsize - 1; l >= 0; l--) {
											mat[row + k][col + l] =
													nbrmat[k][l];
										}
									}
									break;
								}
							}
						}
						col += kidnbr->directnumeles[0];
					}
				}
				ASSERT(col == size);
				row += kidsize;
			}
		}
		ASSERT(row == size);

		nc->precond = ludecomp(mat, size, FALSE);
		nc->presize = size;
	}
}

/* This near picks up only the hamming distance one cubes. */
#define HNEAR(nbr, nj, nk, nl) \
((ABS((nbr)->j - (nj)) + ABS((nbr)->k - (nk)) + ABS((nbr)->l - (nl))) <= 1)

/* This near picks up all 27 neighboring cubes. */
#define NEAR(nbr, nj, nk, nl) \
((ABS((nbr)->j - (nj)) <= 1) && \
 (ABS((nbr)->k - (nk)) <= 1) && \
 (ABS((nbr)->l - (nl)) <= 1))

/* This near picks only the diagonal, for testing. */
#define DNEAR(nbr, nj, nk, nl) \
(((nbr)->j == (nj)) && \
 ((nbr)->k == (nk)) && \
 ((nbr)->l == (nl)) )

void olmulMatPrecond(ssystem *sys) {
	cube *nc, *nnbr, *nnnbr;
	double **mat, **nmat;
	int i, j, k, l, m;
	int maxsize, nsize, nnsize, nnnsize, *reorder;
	int nj, nk, nl, offset, noffset;
	int dindex, *nc_dummy, *nnbr_dummy, *nnnbr_dummy;
	static int *is_dummy; /* local dummy flag vector, stays around */
	static int big_mat_size = 0; /* size of previous mat */
	charge **nnnbr_pc, **nnbr_pc, **nc_pc, **mpc, *dp;
	surface *surf;
	double factor;

	printf("This Preconditioner is not used in FastHenry\n");
	exit(1);

	/* Figure out the max number of elements in any set of near cubes. */
	for (maxsize = 0, nc = sys->directlist; nc != NULL; nc = nc->dnext) {
		nsize = nc->directnumeles[0];
		nj = nc->j;
		nk = nc->k;
		nl = nc->l;
		for (i = 0; i < nc->numnbrs; i++) {
			nnbr = nc->nbrs[i];
			if (NEAR(nnbr, nj, nk, nl))
				nsize += nnbr->directnumeles[0];
		}
		maxsize = MAX(nsize, maxsize);
	}

	/* Allocate a matrix big enough for any set of 7. */
#if JACDBG == ON
	printf("max direct size =%d\n", maxsize);
#endif
	MALLOC(reorder, maxsize, int, ON, AMSC);
	MALLOC(mat, maxsize, double*, ON, AMSC);
	for (i = 0; i < maxsize; i++) {
		MALLOC(mat[i], maxsize, double, ON, AMSC);
	}

	/* Now go fill-in a matrix. */
	for (maxsize = 0, nc = sys->directlist; nc != NULL; nc = nc->dnext) {
		nsize = nc->directnumeles[0];
		nc_dummy = nc->nbr_is_dummy[0];
		nc_pc = nc->chgs;
#if CHKDUM == ON
		chkDummyList(nc_pc, nc_dummy, nsize);
#endif
		nj = nc->j;
		nk = nc->k;
		nl = nc->l;
		for (i = nsize - 1; i >= 0; i--) {
			if (nc_dummy[i])
				continue; /* dummy rows copied only in divided diff */
			if (nc_pc[i]->surf->type != DIELEC) {
				for (j = nsize - 1; j >= 0; j--) {
					mat[i][j] = nc->directmats[0][i][j];
				}
			} else {
#if DPCOMP == ON
				fprintf(stdout, "Source mat, nc to nc\n");
				dumpMat(nc->directmats[0], nsize, nsize);
#endif
				find_flux_density_row(mat, nc->directmats[0], i, nsize, nsize,
						0, 0, nc_pc, nc_pc, nc_dummy, nc_dummy);
			}
		}
		offset = nsize;
		for (k = 0; k < nc->numnbrs; k++) { /* loop on neighbors of nc */
			nnbr = nc->nbrs[k];
			if (NEAR(nnbr, nj, nk, nl)) {
				nnsize = nc->directnumeles[k + 1];
				nmat = nc->directmats[k + 1];
				ASSERT(nc->directnumeles[k + 1] == nnbr->directnumeles[0]);
				nnbr_dummy = nnbr->nbr_is_dummy[0];
				nnbr_pc = nnbr->chgs;
#if CHKDUM == ON
				chkDummyList(nnbr_pc, nnbr_dummy, nnsize);
#endif
				for (i = nsize - 1; i >= 0; i--) {
					if (nc_dummy[i])
						continue;
					if (nc_pc[i]->surf->type != DIELEC) {
						for (j = nnsize - 1; j >= 0; j--) {
							mat[i][offset + j] = nmat[i][j];
						}
					} else {
#if DPCOMP == ON
						fprintf(stdout, "Source mat, nnbr to nc\n");
						dumpMat(nmat, nsize, nnsize);
#endif
						find_flux_density_row(mat, nmat, i, nnsize, nsize, 0,
								offset, nc_pc, nnbr_pc, nc_dummy, nnbr_dummy);
					}
				}
				/* Get the row of the big matrix associated with this nnbr. */
				for (noffset = 0, l = -1; l < nc->numnbrs; l++) { /* lp on nc's nbrs */
					if (l < 0)
						nnnbr = nc;
					else
						nnnbr = nc->nbrs[l];
					if (NEAR(nnnbr, nj, nk, nl)) { /* Note, near to nc!! */
						if (nnbr == nnnbr)
							m = -1;
						else { /* Find this nnnbr's position in nnbr's list */
							for (m = 0; m < nnbr->numnbrs; m++) {
								if (nnbr->nbrs[m] == nnnbr)
									break;
							}
							ASSERT(m < nnbr->numnbrs);
						}
						nnnsize = nnbr->directnumeles[m + 1];
						nmat = nnbr->directmats[m + 1];
						ASSERT(
								nnbr->directnumeles[m + 1]
										== nnnbr->directnumeles[0]);
						nnnbr_pc = nnnbr->chgs; /* panels in nnnbr */
						nnnbr_dummy = nnnbr->nbr_is_dummy[0];
#if CHKDUM == ON
						chkDummyList(nnnbr_pc, nnnbr_dummy, nnnsize);
#endif
						for (i = nnsize - 1; i >= 0; i--) { /* loop on panels in nnbr */
							if (nnbr_dummy[i])
								continue;
							if (nnbr_pc[i]->surf->type != DIELEC) {
								for (j = nnnsize - 1; j >= 0; j--) {
									mat[offset + i][noffset + j] = nmat[i][j];
								}
							} else {
#if DPCOMP == ON
								fprintf(stdout, "Source mat, nnnbr to nnbr\n");
								dumpMat(nmat, nnsize, nnnsize);
#endif
								find_flux_density_row(mat, nmat, i, nnnsize,
										nnsize, offset, noffset, nnbr_pc,
										nnnbr_pc, nnbr_dummy, nnnbr_dummy);
							}
						}
						noffset += nnnsize;
					}
				}
				offset += nnsize;
			}
		}

		/* set up the local is_dummy vector for the rows/cols of mat */
		/* THIS COULD BE AVOIDED BY USING CUBE is_dummy's INSIDE invert() */
		if (big_mat_size < offset) { /* allocate only if larger array needed */
			CALLOC(is_dummy, offset, int, ON, AMSC);
		}
		/* dump sections of the dummy vector in order cubes appear in nbr lst */
		/* (use fragment of Jacob's loop above) */
		nnnsize = noffset = nc->directnumeles[0];
		nc_dummy = nc->nbr_is_dummy[0];
		for (i = nnnsize - 1; i >= 0; i--) {
			is_dummy[i] = nc_dummy[i];
		}
		for (l = 0; l < nc->numnbrs; l++) {
			nnnbr = nc->nbrs[l];
			if (NEAR(nnnbr, nj, nk, nl)) {
				nnnsize = nnnbr->directnumeles[0];
				nc_dummy = nnnbr->nbr_is_dummy[0];
				for (i = nnnsize - 1; i >= 0; i--) {
					is_dummy[i + noffset] = nc_dummy[i];
				}
				noffset += nnnsize;
			}
		}

		/* The big Matrix is filled in, invert it and get the preconditioner. */
#if DPCOMP == ON
		fprintf(stdout, "Before compression\n");
		dumpMat(mat, offset, offset);
#endif
		nnnsize = compressMat(mat, offset, is_dummy, BOTH);
#if DPCOMP == ON
		fprintf(stdout, "After compression\n");
		dumpMat(mat, nnnsize, nnnsize);
#endif
		invert(mat, nnnsize, NULL);
		expandMat(mat, offset, nnnsize, is_dummy, BOTH);
#if DPCOMP == ON
		fprintf(stdout, "After expansion\n");
		dumpMat(mat, offset, offset);
#endif

		/* Copy out the preconditioner to the saved matrices. */
		for (i = nsize - 1; i >= 0; i--) {
			for (j = nsize - 1; j >= 0; j--) {
				nc->precondmats[0][i][j] = mat[i][j];
			}
		}
		offset = nsize;
		for (k = 0; k < nc->numnbrs; k++) {
			nnbr = nc->nbrs[k];
			if (NEAR(nnbr, nj, nk, nl)) {
				nnsize = nc->directnumeles[k + 1];
				nmat = nc->precondmats[k + 1];
				for (i = nsize - 1; i >= 0; i--) {
					for (j = nnsize - 1; j >= 0; j--) {
						nmat[i][j] = mat[i][offset + j];
					}
				}
				offset += nnsize;
			} else
				nc->precondmats[k + 1] = NULL;
		}
	}
}

/*
 finds a row of flux density coeffs from three potential coeff rows
 - to_mat[eval_row][] is the destination row; from_mat[eval_row][]
 initially contains the potential coefficients for evals at the
 center of eval_panels[eval_row] (unless NUMDPT == 2, is garbage then)
 - the eval panels are scaned until eval_panels[eval_row]'s
 dummies are found and the corresponding two rows are identified
 - the divided differences built with entries in the same columns in
 these three rows replace the to_mat[eval_row][] entries:
 to_mat[eval_row][j] = a1*from_mat[eval_row][j]
 + a2*from_mat[pos_dum_row][j] + a3*from_mat[neg_dum_row][j]
 - if a dummy panel is not found in the panel list, its row is generated
 using explicit calcp() calls (shouldn't happen much)
 - global flags used here
 NUMDPT = number of divided diff points, 2 or 3
 SKIPDQ = ON=>don't do cancellation-prone add-subtract of identical
 influence of DIELEC/BOTH panels' charges on dummy panel pot. evals
 */
void find_flux_density_row(double **to_mat, double **from_mat, int eval_row,
		int n_chg, int n_eval, int row_offset, int col_offset,
		charge **eval_panels, charge **chg_panels, int *eval_is_dummy,
		int *chg_is_dummy) {
	int dindex, j;
	double factor;
	charge *dp;
	surface *surf = eval_panels[eval_row]->surf;

	/* do divided difference w/ three rows to get dielectric row */
#if NUMDPT == 3
	/* - dielectric panel row first */
	factor = -(surf->outer_perm + surf->inner_perm)/
	(eval_panels[eval_row]->pos_dummy->area);
#if DPDDIF == ON
	fprintf(stdout, "Center row, factor = %g\n", factor);
#endif
	for(j = n_chg - 1; j >= 0; j--) { /* loop on columns */
		if(!chg_is_dummy[j])
		to_mat[row_offset + eval_row][col_offset + j]
		= from_mat[eval_row][j]*factor;
#if DPDDIF == ON
		fprintf(stdout, " %.16e", from_mat[eval_row][j]);
#endif				/* #if DPDDIF == ON */
	}
#endif				/* #if NUMDPT == 3 */
	/* - do positive dummy row */
	/*   first find the dummy row */
	dindex = -1;
	dp = eval_panels[eval_row]->pos_dummy; /* get dummy panel from eval panel */
	for (j = n_eval - 1; j >= 0; j--) {
		if (!eval_is_dummy[j])
			continue;
		if (dp == eval_panels[j]) {
			dindex = j;
			break;
		}
	}
	if (dindex != -1) { /* dummy row found */
#if NUMDPT == 3
		factor = surf->outer_perm/eval_panels[dindex]->area;
#else
		/* this is the only factor required for two dummy rows in two point case */
		factor = (surf->inner_perm - surf->outer_perm)
				/ (eval_panels[eval_row]->neg_dummy->area
						+ eval_panels[eval_row]->pos_dummy->area);
#endif
#if DPDDIF == ON
		fprintf(stdout, "\nPos dummy row, factor = %g\n", factor);
#endif
		for (j = n_chg - 1; j >= 0; j--) {
#if SKIPDQ == ON
			if(chg_panels[j]->index == eval_panels[eval_row]->index) {
				to_mat[row_offset + eval_row][col_offset + j] = 0.0;
				continue;
			}
#endif
			if (!chg_is_dummy[j])
#if NUMDPT == 3
				to_mat[row_offset + eval_row][col_offset + j]
				+= from_mat[dindex][j]*factor;
#else				/* make sure to overwrite possible garbage */
				to_mat[row_offset + eval_row][col_offset + j] =
						-from_mat[dindex][j] * factor;
#endif
#if DPDDIF == ON
			fprintf(stdout, " %.16e (%d)", from_mat[dindex][j],chg_panels[j]->index);
#endif
		}
	} else { /* dummy row out of cube => build it w/calcp */
#if NUMDPT == 3
		factor = surf->outer_perm/dp->area;
#else
		/* this is the only factor required for two dummy rows in two point case */
		factor = (surf->inner_perm - surf->outer_perm)
				/ (eval_panels[eval_row]->neg_dummy->area
						+ eval_panels[eval_row]->pos_dummy->area);
#endif
#if DPDDIF == ON
		fprintf(stdout, "\nPos dummy calcp row, factor = %g\n", factor);
#else
		fprintf(stderr, "\nolmulMatPrecond: building pos. dummy row\n");
#endif
		for (j = n_chg - 1; j >= 0; j--) {
#if SKIPQD == ON
			if(chg_panels[j]->index == eval_panels[eval_row]->index) {
				to_mat[row_offset + eval_row][col_offset + j] = 0.0;
				continue;
			}
#endif
			if (!chg_is_dummy[j]) {
#if NUMDPT == 3
				to_mat[row_offset + eval_row][col_offset + j]
				+= calcp(chg_panels[j], dp, NULL)*factor;
#else
				to_mat[row_offset + eval_row][col_offset + j] = -calcp(
						chg_panels[j], dp, NULL) * factor;
#endif
#if DPDDIF == ON
				fprintf(stdout, " %.16e (%d)",
						calcp(chg_panels[j], dp, NULL),
						chg_panels[j]->index);
			}
			else {
				fprintf(stdout, " dummy");
#endif
			}
		}
	}
	/* - do negative dummy row */
	/*   first find the dummy row */
	dindex = -1;
	dp = eval_panels[eval_row]->neg_dummy; /* get dummy panel from eval panel */
	for (j = n_eval - 1; j >= 0; j--) {
		if (!eval_is_dummy[j])
			continue;
		if (dp == eval_panels[j]) {
			dindex = j;
			break;
		}
	}
	if (dindex != -1) { /* dummy row found */
#if NUMDPT == 3
		factor = surf->inner_perm/eval_panels[dindex]->area;
#endif
#if DPDDIF == ON
		fprintf(stdout, "\nNeg dummy row, factor = %g\n", factor);
#endif
		for (j = n_chg - 1; j >= 0; j--) {
#if SKIPQD == ON
			if(chg_panels[j]->index == eval_panels[eval_row]->index) continue;
#endif
			if (!chg_is_dummy[j])
				to_mat[row_offset + eval_row][col_offset + j] +=
						from_mat[dindex][j] * factor;
#if DPDDIF == ON
			fprintf(stdout, " %.16e (%d)", from_mat[dindex][j],chg_panels[j]->index);
#endif
		}
	} else { /* dummy row out of cube => build it w/calcp */
		factor = surf->inner_perm / dp->area;
#if DPDDIF == ON
		fprintf(stdout, "\nNeg dummy calcp row, factor = %g\n", factor);
#else
		fprintf(stderr, "olmulMatPrecond: building neg. dummy row\n");
#endif
		for (j = n_chg - 1; j >= 0; j--) {
#if SKIPQD == ON
			if(chg_panels[j]->index == eval_panels[eval_row]->index) continue;
#endif
			if (!chg_is_dummy[j]) {
				to_mat[row_offset + eval_row][col_offset + j] += calcp(
						chg_panels[j], dp, NULL) * factor;
#if DPDDIF == ON
				fprintf(stdout, " %.16e (%d)",
						calcp(chg_panels[j], dp, NULL),
						chg_panels[j]->index);
			}
			else {
				fprintf(stdout, " dummy");
#endif
			}
		}
	}
#if NUMDPT == 2
	/* - do row entry due to panel contribution
	 - entry only necessary if eval panel is in chg panel list */
	/*   search for the eval panel in the charge panel list */
	dp = NULL;
	for (j = n_chg - 1; j >= 0; j--) {
		if (!chg_is_dummy[j]) {
			if (eval_panels[eval_row] == chg_panels[j]) {
				dp = eval_panels[eval_row];
				break;
			}
		}
	}
	/*   set entry if eval panel found in chg panel list
	 - this is an overwrite; contributions of other rows should cancel */
	if (dp != NULL) {
		to_mat[row_offset + eval_row][col_offset + j] = -(2 * M_PI
				* (surf->inner_perm + surf->outer_perm)
				/ eval_panels[eval_row]->area);
	}
#endif

#if DPDDIF == ON
	fprintf(stdout, "\nDivided difference row (%d)\n",
			eval_panels[eval_row]->index);
	for(j = n_chg - 1; j >= 0; j--) {
		fprintf(stdout, " %.16e (%d)",
				to_mat[row_offset + eval_row][col_offset + j],
				chg_panels[j]->index);
	}
	fprintf(stdout, "\n\n");
#endif
}

/* 
 MulMatUp computes the multipole to multipole or charge to
 multipole matrices that map to a parent's multipole coeffs from its
 children's multipoles or charges. Note that only one set of
 multipole to multipole matrices is computed per level by exploiting the
 uniform break-up of three-space (ie many shifts have similar geometries).
 */
void mulMatUp(ssystem *sys) {
	cube *nextc, *kid;
	int i, j, numterms, depth, order = sys->order;
	double **multimats[8];

#if OFF == ON			/* OFF == OFF produces the M2M error!! */
	for(i=0; i < 8; i++) multimats[i] = NULL;
#endif

	numterms = multerms(order);

	if (sys->depth < 2) {
		/* fprintf(stdout, "\nWarning: no multipole acceleration\n");*/
		return; /* return if upward pass not possible */
	}

	/* Handle the lowest level cubes first (set up Q2M's). */
	for (nextc = sys->multilist[sys->depth]; nextc != NULL; nextc =
			nextc->mnext) {
		nextc->multisize = numterms;
		CALLOC(nextc->multi, numterms, double, ON, AMSC);
		CALLOC(nextc->upmats, 1, double**, ON, AMSC);
		nextc->upmats[0] = mulQ2Multi(nextc->chgs, nextc->nbr_is_dummy[0],
				nextc->upnumeles[0], nextc->x, nextc->y, nextc->z, order);

#if DISSYN == ON
		multicnt[nextc->level]++;
#endif

#if DMTCNT == ON
		Q2Mcnt[nextc->level][nextc->level]++;
#endif

	}

	if (sys->locallist[sys->depth] == NULL && sys->multilist[sys->depth] == NULL) {
		fprintf(stdout, "No expansions at level %d (lowest)\n", sys->depth);
		/*if(sys->depth < 3)
		 fprintf(stdout, " (Warning: no multipole acceleration)\n");*/
	} else if (sys->locallist[sys->depth] == NULL) {
		fprintf(stdout, "No local expansions at level %d (lowest)\n",
				sys->depth);
	} else if (sys->multilist[sys->depth] == NULL) {
		fprintf(stdout, "No multipole expansions at level %d (lowest)\n",
				sys->depth);
		/*if(sys->depth < 3)
		 fprintf(stdout, " (Warning: no multipole acceleration)\n");*/
	}

	/* Allocate the vectors and matrices for the cubes. */
	/* no multipoles over root cube or its kids (would not be used if made) */
	for (depth = (sys->depth - 1); depth > 1; depth--) {
		/* set up M2M's and Q2M's to compute the multipoles needed for this level */
		if (sys->locallist[depth] == NULL && sys->multilist[depth] == NULL) {
			fprintf(stdout, "No expansions at level %d\n", depth);
			/*if(depth < 3) fprintf(stdout, " (Warning: no multipole acceleration)\n");
			 else fprintf(stdout, "\n");*/
		} else if (sys->locallist[depth] == NULL) {
			fprintf(stdout, "No local expansions at level %d\n", depth);
		} else if (sys->multilist[depth] == NULL) {
			fprintf(stdout, "No multipole expansions at level %d\n", depth);
			/*if(depth < 3) fprintf(stdout, " (Warning: no multipole acceleration)\n");
			 else fprintf(stdout, "\n");*/
		}

#if ON == ON			/* ON == OFF produces the M2M error!! */
		/* NULL out pointers to same-geometry M2M mats for this level */
		for (i = 0; i < 8; i++)
			multimats[i] = NULL;
#endif

		/* Hit nonempty cubes at this level assigning ptrs to precomputed   */
		/* M2M mats (for this lev), or if a kid is exact, computing Q2M matrices. */
		for (nextc = sys->multilist[depth]; nextc != NULL; nextc =
				nextc->mnext) {

#if DISSYN == ON
			multicnt[nextc->level]++;
#endif

			/* Save space for upvector sizes, upvect ptrs, and upmats. */
			nextc->multisize = numterms;
			if (numterms > 0) {
				CALLOC(nextc->multi, numterms, double, ON, AMSC);
			}
			if (nextc->upnumvects) {
				CALLOC(nextc->upnumeles, nextc->upnumvects, int, ON, AMSC);
				CALLOC(nextc->upvects, nextc->upnumvects, double*, ON, AMSC);
				CALLOC(nextc->upmats, nextc->upnumvects, double**, ON, AMSC);
			}
			/* Go through nonempty kids and fill in upvectors and upmats. */
			for (i = 0, j = 0; j < nextc->numkids; j++) {
				if ((kid = nextc->kids[j]) != NULL) { /* NULL => empty kid cube */
					if (kid->mul_exact == FALSE) { /* if kid has a multi */
						nextc->upvects[i] = kid->multi;
						nextc->upnumeles[i] = kid->multisize;
						if (multimats[j] == NULL) { /* Build the needed matrix only once. */
							multimats[j] = mulMulti2Multi(kid->x, kid->y,
									kid->z, nextc->x, nextc->y, nextc->z,
									order);
						}
						nextc->upmats[i] = multimats[j];

#if DMTCNT == ON
						M2Mcnt[kid->level][nextc->level]++; /* cnts usage, ~computation */
#endif				/* only at most 8 mats really built/level */

					} else { /* if kid is exact, has no multi */
						nextc->upvects[i] = kid->upvects[0];
						nextc->upnumeles[i] = kid->upnumeles[0];
						nextc->upmats[i] = mulQ2Multi(kid->chgs,
								kid->nbr_is_dummy[0], kid->upnumeles[0],
								nextc->x, nextc->y, nextc->z, order);

#if DMTCNT == ON
						Q2Mcnt[kid->level][nextc->level]++;
#endif

					}
					i++; /* only increments if kid is not empty */
				}
			}
		}
	}
}

/*
 builds the transformation matrices for the final evaluation pass (M2P, L2P)
 for all the direct list (generated by linkcubes(), non-empty
 lowest level cubes) cubes:

 for each cube A in the direct list:
 1) if the cube is not exact (always the case if ADAPT = OFF)
 a) and if DNTYPE = GRENGD build an L2P matrix from A to A
 b) and if DNTYPE = NOSHFT build an L2P matrix from each of A's ancestors
 with level > 1 (including A) to A
 c) and if DNTYPE = NOLOCL build an M2P matrix from each of A's fake
 ilist entries to A (same action as 2b)
 2) if the cube is exact, find the 1st ancestor of A, cube B,
 which either is not exact and is at level 2,3,4... or is at level 1
 a) if B is at level 2,3,4...
 i) if DNTYPE = GRENGD, construct an L2P from B to A and M2P's
 from the cubes in the true interaction lists of A and all its
 ancestors up to and including B (a partial fake interaction list)
 j) if DNTYPE = NOSHFT, find cube C, the ancestor of B at level 1;
 construct L2P's from the ancestors of B (including B but not C)
 to A and Q- or M2P's from the cubes in the true interaction lists
 of A and all its ancestors up to and including B (a partial fake
 interaction list)
 k) if DNTYPE = NOLOCL, do 2b
 b) if B is at level 1 construct M2P's for all the cubes in A's
 fake interaction list

 true interaction list - RADINTER = OFF, those sibling
 (same level) cubes of a given cube who are children of the neighbors
 of the given cube's parent and are not neighbors of the given cube
 - ie those cubes required to cover charges well separated from the given
 cube but not accounted for in the parent's local expansion
 - the flag NNBRS is the number of sibling cube "shells" taken as neighbors

 fake interaction list - RADINTER = OFF, the combined true interaction lists
 of a given cube and all its ancestors at levels 2,3,4...

 if RADINTER = ON, any 8 siblings of the given cube which form a well
 separated cube one level up are included in the lists as a single higher
 level cube

 if ADAPT = OFF, no cube is exact so step 2 is never done

 this routine is used alone if compiled with DNTYPE = NOLOCL or after
 mulMatDown, which produces M2L and L2L matrices (DNTYPE = GRENGD) or
 just M2L matrices (DNTYPE = NOSHFT) --  DNTYPE = GRENGD does the full
 Greengard hiearchical downward pass

 */
void mulMatEval(ssystem *sys) {
	int i, j, k, ttlvects, vects;
	cube *na, *nc, *nexti;

	if (sys->depth < 2)
		return; /* ret if upward pass not possible/worth it */

	for (nc = sys->directlist; nc != NULL; nc = nc->dnext) {

		ASSERT(nc->level == sys->depth);
		ASSERT(nc->upnumvects > 0);

		/* allocate space for evaluation pass vectors; check nc's ancestors */
		/* First count the number of transformations to do. */
		for (na = nc, ttlvects = 0; na->level > 1; na = na->parent) {
			if (na->loc_exact == FALSE && DNTYPE != NOLOCL) {
				ttlvects++; /* allow for na to na local expansion (L2P) */
				if (DNTYPE == GRENGD)
					break; /* Only one local expansion if shifting. */
			} else {
				ttlvects += na->interSize; /* room for Q2P and M2P xformations */
			}
		}
		nc->evalnumvects = ttlvects; /* save ttl # of transformations to do */
		if (ttlvects > 0) {
			CALLOC(nc->evalvects, ttlvects, double*, ON, AMSC);
			CALLOC(nc->evalnumeles, ttlvects, int, ON, AMSC);
			CALLOC(nc->evalmats, ttlvects, double**, ON, AMSC);
			CALLOC(nc->eval_isQ2P, ttlvects, int, ON, AMSC);
		}

#if DILIST == ON
		fprintf(stdout, "\nInteraction list (%d entries) for ", ttlvects);
		disExParsimpcube(nc);
#endif

		/* set up exp/charge vectors and L2P, Q2P and/or M2P matrices as req'd */
		for (j = 0, na = nc, ttlvects = 0; na->level > 1; na = na->parent) {
			if (na->loc_exact == FALSE && DNTYPE != NOLOCL) {
				/* build matrices for local expansion evaluation */
				nc->evalmats[j] = mulLocal2P(na->x, na->y, na->z, nc->chgs,
						nc->upnumeles[0], sys->order);
				nc->evalnumeles[j] = na->localsize;
				nc->evalvects[j] = na->local;
				/*	add_to_counts(nc, na->localsize, sys->evalL2Ps, sys->cntL2Ps); */
				nc->eval_isQ2P[j] = FALSE;
				j++;

#if DMTCNT == ON
				L2Pcnt[na->level][nc->level]++;
#endif

#if DILIST == ON
				fprintf(stdout, "L2P: ");
				disExtrasimpcube(na);
#endif
				if (DNTYPE == GRENGD)
					break; /* Only one local expansion if shifting. */
			} else { /* build matrices for ancestor's (or cube's if 1st time) ilist */
				for (i = 0; i < na->interSize; i++) {
					nexti = na->interList[i];
					if (nexti->mul_exact == TRUE) {
						nc->evalvects[j] = nexti->upvects[0];
						nc->evalmats[j] = Q2P(nexti->chgs, nexti->upnumeles[0],
								nexti->nbr_is_dummy[0], nc->chgs,
								nc->upnumeles[0], TRUE);
						/* this is a hack to fix the fact that direct stuff is */
						/* don't done componentwise */

						/* obsolete as of 12/92  due to eval_isQ2P stuff
						 fixEvalDirect(nexti->chgs, nexti->upnumeles[0],
						 nexti->nbr_is_dummy[0], nc->chgs,
						 nc->upnumeles[0], nc->evalmats[j]);
						 */

						nc->evalnumeles[j] = nexti->upnumeles[0];
						nc->eval_isQ2P[j] = TRUE;
						/*	    add_to_counts(nc, nexti->upnumeles[0], sys->evalQ2Ps, sys->cntQ2Ps); */
						j++;

#if DMTCNT == ON
						Q2Pcnt[nexti->level][nc->level]++;
#endif

#if DILIST == ON
						fprintf(stdout, "Q2P: ");
						disExtrasimpcube(nexti);
#endif
					} else {
						nc->evalvects[j] = nexti->multi;
						nc->evalmats[j] = mulMulti2P(nexti->x, nexti->y,
								nexti->z, nc->chgs, nc->upnumeles[0],
								sys->order);
						nc->evalnumeles[j] = nexti->multisize;
						nc->eval_isQ2P[j] = FALSE;
						/*	    add_to_counts(nc, nexti->multisize, sys->evalM2Ps, sys->cntM2Ps);*/

						j++;

#if DMTCNT == ON
						M2Pcnt[nexti->level][nc->level]++;
#endif

#if DILIST == ON
						fprintf(stdout, "M2P: ");
						disExtrasimpcube(nexti);
#endif
					}
				}
			}
		}
	}
}

/* 
 sets up matrices for the downward pass
 For each cube in local list (parents always in list before kids):
 1) parent's local to child's local unless DNTYPE = NOSHFT or no parent local
 2) multipoles for (Parent+parent's nbrs - child nbrs) to child's local
 -eval is sum of ancestral local evals for each lowest lev cube if NOSHFT
 otherwise only lowest level local is evaluated (see mulMatEval)
 -with ADAPT = OFF no cube is exact so local list is all non-empty cube lev>1
 -mats that give potentials (M2P, L2P, Q2P) are calculated in mulMatEval()
 -this routine makes only L2L, M2L and Q2L matrices
 */
void mulMatDown(ssystem *sys) {
	int i, j, vects;
	cube *nc, *parent, *ni;
	int depth;

	ASSERT(DNTYPE != NOLOCL); /* use mulMatEval() alone if NOLOCL */

	for (depth = 2; depth <= sys->depth; depth++) { /* no locals before level 2 */
		for (nc = sys->locallist[depth]; nc != NULL; nc = nc->lnext) {

			/* Allocate for interaction list, include one for parent if needed */
			if ((depth <= 2) || (DNTYPE == NOSHFT))
				vects = nc->interSize;
			else
				vects = nc->interSize + 1;
			nc->downnumvects = vects;
			if (vects > 0) {
				CALLOC(nc->downvects, vects, double*, ON, AMSC);
				CALLOC(nc->downnumeles, vects, int, ON, AMSC);
				CALLOC(nc->downmats, vects, double**, ON, AMSC);
			}

			parent = nc->parent;
			ASSERT(parent->loc_exact == FALSE); /* has >= #evals of any of its kids*/

#if DISSYN == ON
			localcnt[nc->level]++;
#endif

			if ((depth <= 2) || (DNTYPE == NOSHFT))
				i = 0; /* No parent local. */
			else { /* Create the mapping matrix for the parent to kid. */
				i = 1;

				nc->downmats[0] = mulLocal2Local(parent->x, parent->y,
						parent->z, nc->x, nc->y, nc->z, sys->order);
				nc->downnumeles[0] = parent->localsize;
				nc->downvects[0] = parent->local;

#if DMTCNT == ON
				L2Lcnt[parent->level][nc->level]++;
#endif
			}

			/* Go through the interaction list and create mapping matrices. */
			for (j = 0; j < nc->interSize; j++, i++) {
				ni = nc->interList[j];
				if (ni->mul_exact == TRUE) { /* ex->ex (Q2P) xforms in mulMatEval */
					nc->downvects[i] = ni->upvects[0];
					nc->downmats[i] = mulQ2Local(ni->chgs, ni->upnumeles[0],
							ni->nbr_is_dummy[0], nc->x, nc->y, nc->z,
							sys->order);
					nc->downnumeles[i] = ni->upnumeles[0];
#if DMTCNT == ON
					Q2Lcnt[ni->level][nc->level]++;
#endif
				} else {
					nc->downvects[i] = ni->multi;
					nc->downmats[i] = mulMulti2Local(ni->x, ni->y, ni->z, nc->x,
							nc->y, nc->z, sys->order);
					nc->downnumeles[i] = ni->multisize;
#if DMTCNT == ON
					M2Lcnt[ni->level][nc->level]++;
#endif
				}
			}
		}
	}
}

