extern "C"
{
#include "induct.h"//BAPN change
}



/*This calculates all U and V matrices by calling ACA_new to apply ACA and then SVD_QR
to apply QR decomposition and SVD recompression resulting in the final U and V.*/
void MLACA_cal(int depth, int *number_of_groups_level, int **number_of_interactions_group, int ***interactions, int **num_elements, int ***elements, charge **filaments, ssystem *sys, double tol)
{
	clock_t t;
	t = clock();
	sys->distance *= 4;
	int svd_mem = 0;
	int hshshshshss = 0;
	sys->k_max = 0;
	sys->el_max = 0;
	sys->count = 0;
	for (int i = 0; i < depth; i++)
	{
		for (int j = 0; j < number_of_groups_level[i]; j++)
		{
			for (int k = 0; k < number_of_interactions_group[i][j]; k++)
			{
				++sys->count;
			}
		}
	}

	sys->U = (double ***)malloc(sys->count * sizeof(double **));

	sys->V = (double ***)malloc(sys->count * sizeof(double **));


	sys->k = (int *)malloc(sys->count * sizeof(int));
	int count2 = 0;
	int memory = 0;
	long long int ACA_memory = 0;

	int ASC_count = 0; //ASC grouping
	int memory_ASC = 0;
	for (int i = 0; i < depth; i++)
	{
		//ssys->interact[i] = (int**)calloc(number_of_groups_level[i + 1], sizeof(int*));
		//ssys->number_of_interactions_groups[i] = (int *)malloc(number_of_groups_level[i] * sizeof(int));
		for (int j = 0; j < number_of_groups_level[i]; j++)
		{
			//ssys->interact[i][j] = (int**)calloc(number_of_groups_level[i + 1], sizeof(int*));
			//++ssys->number_groups_level[i + 1];
			//ssys->number_of_interactions_groups[i][j] = number_of_interactions_group[i][j];
			if (sys->el_max < num_elements[i][j]) sys->el_max = num_elements[i][j];

			//Grouping ASC
			ASC_count = 0;
			//end ASC grouping
			for (int k = 0; k < number_of_interactions_group[i][j]; k++)
			{



				double **UU = (double **)malloc(num_elements[i][j] * sizeof(double*));

				//ssys->U[count2] = (long double **)calloc(num_elements[i + 1][j], sizeof(long double*));

				//sys->V[count2] = (double **)malloc(num_elements[i][interactions[i][j][k]]* sizeof(double*));//k
				//sys->V[count2] = (double **)malloc(num_elements[i][interactions[i][j][k]] * sizeof(double*));//k
				double **VV = (double **)malloc(num_elements[i][interactions[i][j][k]] * sizeof(double*));//k


																										  //sys->V[count2] = (double **)malloc(number_of_interactions_group[i][j] * sizeof(double*));//k
																										  //ssys->V[count2] = (long double **)calloc(num_elements[i + 1][interactions[i][j][k]], sizeof(long double*));//k
																										  //	for (int q = 0; q < num_elements[i][j]; q++)
																										  //{
																										  //sys->U[count2][q] = (double *)calloc(num_elements[i + 1][j], sizeof(double));//k
																										  //sys->U[count2][0] = (double *)malloc(num_elements[i][j] * sizeof(double));
				UU[0] = (double *)malloc(num_elements[i][j] * sizeof(double));


				//ssys->U[count2][q] = (long double *)calloc(num_elements[i + 1][j], sizeof(long double));//k
				//}
				//for (int p = 0; p < num_elements[i][interactions[i][j][k]]; p++)//k
				//{
				//sys->V[count2][p] = (double *)calloc(num_elements[i + 1][interactions[i][j][k]], sizeof(double));
				//sys->V[count2][0] = (double *)malloc(num_elements[i][interactions[i][j][k]] * sizeof(double));
				VV[0] = (double *)malloc(num_elements[i][interactions[i][j][k]] * sizeof(double));




				//ssys->V[count2][p] = (long double *)calloc(num_elements[i + 1][interactions[i][j][k]], sizeof(long double));
				//}



				double distance_x = fabs(sys->Center[i][j][0] - sys->Center[i][k][0]);
				double distance_y = fabs(sys->Center[i][j][1] - sys->Center[i][k][1]);
				double distance_z = fabs(sys->Center[i][j][2] - sys->Center[i][k][2]);


				memory = ACA_new(num_elements[i][j], num_elements[i][interactions[i][j][k]], elements[i][j], elements[i][interactions[i][j][k]], filaments, UU, VV, tol);

				if (memory >sys->k_max) sys->k_max = memory;

				ACA_memory += memory*num_elements[i][j] + memory*num_elements[i][interactions[i][j][k]];




				int rank;
				rank = SVD_QR(UU, VV, num_elements[i][j], num_elements[i][interactions[i][j][k]], memory, tol*10);

				svd_mem += rank*num_elements[i][j] + rank*num_elements[i][interactions[i][j][k]];
				sys->k[count2] = rank;



				sys->U[count2] = (double **)malloc(rank * sizeof(double*));
				sys->V[count2] = (double **)malloc(rank * sizeof(double*));//k
				for (int p = 0; p < rank; p++)
				{
					sys->U[count2][p] = (double *)malloc(num_elements[i][j] * sizeof(double));
					for (int q = 0; q < num_elements[i][j]; q++)
					{
						sys->U[count2][p][q] = UU[p][q];
					}
					sys->V[count2][p] = (double *)malloc(num_elements[i][interactions[i][j][k]] * sizeof(double));
					for (int q = 0; q < num_elements[i][interactions[i][j][k]]; q++)
					{
						sys->V[count2][p][q] = VV[p][q];

					}
				}

				for (int p = 0; p < memory; p++)
				{
					free(UU[p]);
				}
				free(UU);

				for (int p = 0; p < memory; p++)//k
				{
					free(VV[p]);
				}
				free(VV);

				++count2;

			}
		}
	}
	t = clock() - t;

	if (tol == 0)
	{

		FILE * fpp;

		fpp = fopen("file_mlaca.txt", "w+");
		fprintf(fpp, "  MLACA memory: %lld\n", ACA_memory);
		fclose(fpp);
		double time_taken = ((double)t) / CLOCKS_PER_SEC; // in seconds

	}
	sys->memory = svd_mem;
	sys->memory_aa = ACA_memory;
}
//End BAP


//Start BAP
/*This recursively calculates all group interactions ACA should be applied to 
as well as storing the remaining group interactions that need to be calculated directly.*/
void branch(int ***Group_kids, double ***Center, int *paren_index, int **kids_index, int level, int group_number, double distance,
	int previous_unknowns, int *unknown_back, int ***interaction, int depth, int **interaction_number, int **still_to_do, int *number_still_to_do, int ffffff)
{


	interaction[level][group_number] = (int *)malloc(previous_unknowns * sizeof(int));




	int unknowns = 0;

	int *unknown = (int *)malloc(previous_unknowns * sizeof(int));

	for (int i = 0; i < previous_unknowns; i++)
	{
		unknown[i] = 0;
	}

	int unknown_nexts = 0;

	int index = 0;

	int source;
	interaction_number[level][group_number] = 0;
	double distance_x;
	double distance_y;
	double distance_z;
	for (int i = 0; i < previous_unknowns; i++)//Go through all possible interactions on the level
	{
		source = unknown_back[i];

		distance_x = fabs(Center[level][group_number][0] - Center[level][source][0]);
		distance_y = fabs(Center[level][group_number][1] - Center[level][source][1]);
		distance_z = fabs(Center[level][group_number][2] - Center[level][source][2]);
		if ((distance_x  >  distance * 2 || distance_y > distance * 2 || distance_z > distance * 2 || ffffff == 3))
		{
			//this is where the interactions are determained and stored in no particular order
			interaction[level][group_number][index] = source;
			++index;
			++interaction_number[level][group_number];
		}
		else if (depth == level + 1)
		{
			still_to_do[group_number][number_still_to_do[group_number]] = source;
			//memoryyy[0] += num_elements[level][group_number] * num_elements[level][still_to_do[group_number][number_still_to_do[group_number]]];
			++number_still_to_do[group_number];

		}
		else
		{
			//This is where all the pices that need to still be investigated on the next level should be stored.
			unknown[unknowns] = source;
			++unknowns;
		}


	}

	//function compairing start_observation with all the possible interactions.
	for (int i = 0; i < unknowns; i++)
	{
		for (int j = 0; j < kids_index[level][unknown[i]]; j++)
		{
			++unknown_nexts;
		}
	}

	unknown_nexts += kids_index[level][group_number] - 1;


	int *unknown_next = (int *)malloc(unknown_nexts * sizeof(int));
	unknown_nexts = 0;

	for (int i = 0; i < unknowns; i++)
	{
		for (int j = 0; j < kids_index[level][unknown[i]]; j++)
		{
			unknown_next[unknown_nexts] = Group_kids[level][unknown[i]][j];

			++unknown_nexts;
		}
	}
	free(unknown);

	for (int i = 0; i < kids_index[level][group_number] - 1; i++)//call each time for every child go down all posible brances
	{
		unknown_next[unknown_nexts] = Group_kids[level][group_number][kids_index[level][group_number] - i - 1];

		++unknown_nexts;
	}




	if (depth == level + 1)
	{
		for (int i = 0; i < number_still_to_do[group_number]; i++)
		{
			//memoryyy[0] += num_elements[level][group_number] * num_elements[level][still_to_do[group_number][i]];
		}


		//memoryyy[0] += num_elements[5][group_number] * num_elements[5][group_number];

	}
	//free(unknown_back);
	if (depth > level - 1)
	{
		for (int i = 0; i < kids_index[level][group_number]; i++)//call each time for every child go down all posible brances
		{
			branch(Group_kids, Center, paren_index, kids_index, level + 1, Group_kids[level][group_number][i], distance / 2,
				unknown_nexts - i, unknown_next, interaction, depth, interaction_number, still_to_do, number_still_to_do, ffffff);
		}
		free(unknown_next);
	}





}

//End BAP


//Start Bap
/*Initiates Branch to recursively compute interactions.*/
void top_level(ssystem *sys, int ddd)
{
	//sys->group_interaction_number;
	int mem = 0;
	int mem2 = 0;
	int *memoryy = (int*)malloc(2 * sizeof(int));
	memoryy[0] = 0;
	memoryy[1] = 0;
	sys->group_interaction_number = (int **)malloc(sys->depth * sizeof(int *));

	sys->next_iteration = (int **)malloc(sys->group_num * sizeof(int *));
	for (int i = 0; i < sys->group_num; i++)
	{
		sys->next_iteration[i] = (int *)malloc(sys->group_num * sizeof(int));
	}

	sys->number_next_iteration = (int *)malloc((sys->group_num) * sizeof(int));
	for (int i = 0; i < sys->group_num; i++)
	{
		sys->number_next_iteration[i] = 0;
	}

	for (int i = 0; i < sys->depth; i++)
	{
		//sys->group_interaction_number[i] = (int *)malloc(sys->group_num * sizeof(int));
		sys->group_interaction_number[i] = (int *)malloc(sys->group_num * sizeof(int));
	}

	sys->observation;
	sys->distance *= pow(2, sys->depth - ddd)*1.000001;//1:3rd,2:2nd
	sys->unknown_num = 0;
	sys->parent_index = 0;
	int y = 0;
	int *unknown_indices = (int *)malloc(sys->number_level_groups[0] * sizeof(int));
	//int *unknown_indices = (int *)malloc(sys->unknown_num * sizeof(int));
	sys->interactions = (int ***)malloc(sys->depth * sizeof(int**));
	for (int i = 0; i < sys->depth; i++)
	{
		sys->interactions[i] = (int **)malloc(sys->number_level_groups[i] * sizeof(int *));
		for (int j = 0; j < sys->number_level_groups[i]; j++)
		{
			//sys->interactions[i][j] = (int *)malloc((sys->number_level_groups[i] - 1) * sizeof(int));
		}

	}

	for (int i = 0; i < sys->number_level_groups[0]; i++)
	{
		sys->unknown_num = 0;
		if (y == sys->kids_index[0][sys->parent_index])
		{
			++sys->parent_index;
			y = 0;
		}
		sys->observation = sys->Parent_kids[0][sys->parent_index][y];

		++y;
		for (int j = i + 1; j < sys->number_level_groups[0]; j++)
		{
			unknown_indices[sys->unknown_num] = j;
			++sys->unknown_num;
			mem += sys->kid_num[1][sys->observation] * sys->kid_num[1][j];
		}
		mem2 += sys->kid_num[1][sys->observation] * sys->kid_num[1][sys->observation];
		branch(sys->Parent_kids, sys->Center, sys->number_level_groups, sys->kids_index, 0, i, sys->distance, sys->unknown_num, unknown_indices, sys->interactions
			, sys->depth, sys->group_interaction_number, sys->next_iteration, sys->number_next_iteration, ddd);
	}
	free(unknown_indices);

	//start free interactions

	//end free interactions




	int whatwhat = 0;
	for (int i = 0; i < sys->group_num; i++)
	{
		//track += sys->kid_num[sys->depth - 1][i] * sys->kid_num[sys->depth - 1][i];
		for (int j = 0; j < sys->number_next_iteration[i]; j++)
		{
			//memoryy[1] += number_els[6][i] * number_els[6][extra[i][j]];
			whatwhat += sys->kid_num[sys->depth - 1][i] * sys->kid_num[sys->depth - 1][sys->next_iteration[i][j]];
			for (int y = 0; y < sys->kid_num[sys->depth - 1][i]; y++)
			{
				for (int t = 0; t < sys->kid_num[sys->depth - 1][sys->next_iteration[i][j]]; t++)
				{
					//total[elements[a + 1][b][i] - 1][elements[a + 1][interactions[a][b][k]][j] - 1]
					//total[els[depth - 1][i][y] - 1][els[depth - 1][extra[i][j]][t] - 1] = calcp(fill[els[depth - 1][i][y] - 1], fill[els[depth - 1][extra[i][j]][t] - 1], NULL);
					//total[els[depth - 1][extra[i][j]][t] - 1][els[depth - 1][i][y] - 1] = calcp(fill[els[depth - 1][i][y] - 1], fill[els[depth - 1][extra[i][j]][t] - 1], NULL);
					//track += 2;
				}
			}

		}

	}

	for (int i = 0; i < sys->group_num; i++)
	{
		for (int y = 0; y < sys->kid_num[sys->depth - 1][i]; y++)
		{
			for (int o = 0; o < sys->kid_num[sys->depth - 1][i]; o++)
			{
				//total[els[depth - 1][i][y] - 1][els[depth - 1][i][o] - 1] = calcp(fill[els[depth - 1][i][y] - 1], fill[els[depth - 1][i][o] - 1], NULL);
				//track+=2;
			}
		}
	}

	//MLACA_cal(sys->depth, sys->number_level_groups, sys->group_interaction_number, sys->interactions, sys->kid_num, sys->child_index, sys->fill);

}
//end Bap




/*Places group information in a format suitable 
for use by the MLACA - SVD algorithm.*/
void Set_up_MLACA(ssystem *sys, int number_fills) {

	cube *nextc;

	/////
	cube *newc;
	////

	///
	sys->distance = 0;
	////
	sys->group_num = 0;

	for (nextc = sys->directlist; nextc != NULL; nextc =
		nextc->dnext)
	{
		++sys->group_num;

	}

	///////
	for (nextc = sys->directlist; nextc != NULL; nextc =
		nextc->dnext) {
		newc = nextc;
		newc->parent_num = -1;
		newc->kid_num = 0;

		for (int i = 0; i < sys->depth; i++)
		{
			newc = newc->parent;
			newc->parent_num = -1;
			newc->kid_num = 0;
			newc->path_open = 0;
		}
	}
	////
	sys->fill = (charge **)malloc(number_fills * sizeof(charge*));
	////

	sys->Center = (double ***)malloc(sys->depth * sizeof(double**));

	for (int i = 0; i < sys->depth; i++) {

		sys->Center[i] = (double **)malloc(sys->group_num * sizeof(double *));


		for (int j = 0; j < sys->group_num; j++) {

			sys->Center[i][j] = (double *)malloc(3 * sizeof(double));

		}

	}


	sys->Parent_kids = (int ***)malloc(sys->depth * sizeof(int**));
	for (int i = 0; i < sys->depth; i++) {

		sys->Parent_kids[i] = (int **)malloc(sys->group_num * sizeof(int *));

		for (int j = 0; j < sys->group_num; j++) {

			sys->Parent_kids[i][j] = (int *)malloc(20 * sizeof(int));
		}

	}


	sys->kids_index = (int **)malloc(sys->depth * sizeof(int*));
	for (int i = 0; i < sys->depth; i++)
	{
		sys->kids_index[i] = (int *)malloc(sys->group_num * sizeof(int));
	}
	for (int i = 0; i < sys->depth; i++)
	{
		for (int j = 0; j < sys->group_num; j++)
		{
			sys->kids_index[i][j] = 0;
		}
	}

	sys->number_level_groups = (int *)malloc(sys->depth * sizeof(int *));
	for (int i = 0; i < sys->depth; i++)
	{

		sys->number_level_groups[i] = 0;
	}
	sys->kid_num = (int **)malloc(sys->depth * sizeof(int *));
	for (int i = 0; i < sys->depth; i++)
	{
		sys->kid_num[i] = (int *)malloc(number_fills * sizeof(int));
	}



	int number_of_filaments = 0;
	int *test_test = (int *)malloc(number_fills * sizeof(int));
	for (nextc = sys->directlist; nextc != NULL; nextc =
		nextc->dnext)
	{
		for (int j = 0; j < nextc->upnumeles[0]; j++)
		{
			charge *temp = nextc->chgs[j];
			sys->fill[nextc->chgs[j]->fil->filnumber] = temp;
			++number_of_filaments;
			test_test[nextc->chgs[j]->fil->filnumber] = number_of_filaments;
		}
	}
	///////
	for (nextc = sys->directlist; nextc != NULL; nextc =
		nextc->dnext) {



		sys->distance = (sqrt(pow((nextc->parent->x - nextc->x), 2) + pow((nextc->parent->y - nextc->y), 2) + pow((nextc->parent->z - nextc->z), 2)) / sqrt(3))*2.001;
		newc = nextc;
		int back_kid = 0;
		int stop = 0;
		for (int i = sys->depth - 1; -1 < i; i--)
		{
			if (newc->parent_num == -1)
			{
				++sys->number_level_groups[i];
				newc->parent_num = sys->number_level_groups[i] - 1;
				sys->Center[i][newc->parent_num][0] = newc->x;
				sys->Center[i][newc->parent_num][1] = newc->y;
				sys->Center[i][newc->parent_num][2] = newc->z;
			}
			if (i < sys->depth - 1 && stop == 0)
			{
				if (newc->path_open == 1) stop = 1;
				sys->Parent_kids[i][newc->parent_num][sys->kids_index[i][newc->parent_num]] = back_kid;
				++sys->kids_index[i][newc->parent_num];
				newc->path_open = 1;
			}
			back_kid = newc->parent_num;

			for (int j = 0; j < nextc->upnumeles[0]; j++)
			{
				//Level_grouping[i][newc->parent_num][j + newc->kid_num] = nextc->chgs[j]->index;
			}
			newc->kid_num = newc->kid_num + nextc->upnumeles[0];
			sys->kid_num[i][newc->parent_num] = newc->kid_num;
			newc = newc->parent;
		}
	}

	////////

	sys->child_index = (int***)malloc(sys->depth * sizeof(int**));
	for (int i = 0; i < sys->depth; i++)
	{
		sys->child_index[i] = (int **)malloc(sys->number_level_groups[i] * sizeof(int*));
		for (int j = 0; j < sys->number_level_groups[i]; j++)
		{
			sys->child_index[i][j] = (int*)malloc(sys->kid_num[i][j] * sizeof(int));
		}
	}




	for (nextc = sys->directlist; nextc != NULL; nextc =
		nextc->dnext) {
		newc = nextc;
		newc->kid_num = 0;

		for (int i = 0; i < sys->depth; i++)
		{
			newc = newc->parent;
			newc->kid_num = 0;
		}
	}

	for (nextc = sys->directlist; nextc != NULL; nextc =
		nextc->dnext) {
		newc = nextc;
		for (int i = sys->depth - 1; -1 < i; i--)
		{

			for (int j = 0; j < nextc->upnumeles[0]; j++)
			{
				sys->child_index[i][newc->parent_num][j + newc->kid_num] = nextc->chgs[j]->fil->filnumber + 1;
			}
			newc->kid_num = newc->kid_num + nextc->upnumeles[0];
			newc = newc->parent;
		}
	}

	//End Bap
}