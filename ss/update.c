#include "ss.h"
#include <string.h>

/*
	The refinment procedure, based on the routine describes in Banga,
	it only replace the refSet members with better candidates by considering the fitness
	and the diversity (equality) to avoid duplicate. The refinement (performing local search)
	will be executed seperately.
 */
void update_ref_set(SSType *ssParams){

	// Sort the candidate set...
	quick_sort_set(ssParams, ssParams->candidates_set, ssParams->candidates_set_size, 'c');
	
	int i = 0;
	if (ssParams->candidates_set->members[i].cost < ssParams->ref_set->members[i].cost ){

		copy_ind(ssParams, &(ssParams->ref_set->members[i]), &(ssParams->candidates_set->members[i]));
		i++;
	}
	
	// The loop continues until there is not better result in candidate set
	while (ssParams->candidates_set->members[i].cost < ssParams->ref_set->members[ ssParams->ref_set_size - 1].cost )
	{
		/////////////////////
		// Diversity Check //
		/////////////////////
		int duplicate_index = is_exist(ssParams, ssParams->ref_set, ssParams->ref_set_size, &(ssParams->candidates_set->members[i]));
		// int duplicate_index = -1;
		if (-1 == duplicate_index)
		{
		   if ( ssParams->perform_flatzone_detection )
		   {
			   	if ( !is_in_flatzone(ssParams, ssParams->ref_set , ssParams->ref_set_size, &(ssParams->candidates_set->members[i])) )
				// if ( ssParams->candidates_set->members[i].cost < ssParams->ref_set->members[ssParams->ref_set_size -  1].cost * ( 1 - ssParams->fitness_epsilon) )
			   	{
					replace(ssParams, &(ssParams->ref_set->members[ssParams->ref_set_size - 1]), &(ssParams->candidates_set->members[i]), 'i');
			   	}
		   }
		   else
		   {
				replace(ssParams, &(ssParams->ref_set->members[ssParams->ref_set_size - 1]), &(ssParams->candidates_set->members[i]), 'i');
		   }
		}else
		{	
			ssParams->n_duplicates++;
		}
			/* Check if the candidate has better fitness than the duplicated member; or at least it has reasonable 
				differences with it (passing the flatzone filter). If yes, then we could replace it.
			*/
				// IMP: This makes problem since it replace duplicated ind with candidate although they are not actually the same
				// 
				// if ( ssParams->candidates_set->members[i].cost < ssParams->ref_set->members[duplicate_index].cost * ( 1 - ssParams->fitness_epsilon) )
				// {
				// 	copy_ind(ssParams, &(ssParams->ref_set->members[duplicate_index]), &(ssParams->candidates_set->members[i]));
				// 	ssParams->n_ref_set_update++;
				// 	// TODO: BubbleSort or QuickSort needed! BubbleSort is perfect since the list is already sorted
				// 	// and one run of BubbleSort will put the replaced member in the right place.
				// 	// bubble_sort(ssParams, ssParams->ref_set, duplicate_index, 'c');
				// 	printf("0\n");
					
				// 	#ifdef STATS
				// 		update_frequency_matrix(ssParams, &(ssParams->candidates_set->members[i]));
				// 	#endif	
				// }

			// }
		// }
		i++;
	}

}

void replace(SSType *ssParams, individual *dest, individual *src, char sort_to_perform){

	copy_ind(ssParams, dest, src);
	insertion_sort(ssParams, ssParams->ref_set, ssParams->ref_set_size, 'c');
	ssParams->n_ref_set_update++;

	#ifdef STATS
		update_frequency_matrix(ssParams, dest);
	#endif		
}

/**
 * [re_gen_ref_set description]
 * @param ssParams [description]
 * @param set      [description]
 * @param set_size [description]
 * @param type     Either `s` or `n`, `s` indicate sorted scatter_set and `n` indidates unsorted scatter_set
 */
void re_gen_ref_set(SSType *ssParams, Set *set, int set_size, char type, Input *inp, ScoreOutput *out){
	int n = ssParams->nreal;
	int b = ssParams->ref_set_size;
	int g = ssParams->max_elite;
	int h = b - g;


	// print_set(ssParams, ssParams->ref_set, ssParams->ref_set_size, ssParams->nreal);
	// refSet:
	//////////// h /////////////
	////////////////////////////
	//	h = b - g 	|		g //
	////////////////////////////
// printf("%d\n", ssParams->scatter_set_size);

	double **M = (double **)malloc(n * sizeof(double*));
	for (int i = 0; i < n; ++i){
		M[i] = (double *)malloc( (b) * sizeof(double));
	}
	// compute_Mt(ssParams, ssParams->ref_set, ssParams->ref_set_size, M, n, h - 2);
	int m = ssParams->scatter_set_size;
	init_scatter_set(ssParams, ssParams->scatter_set);
	if (type == 's'){
		evaluate_set(ssParams, ssParams->scatter_set, ssParams->scatter_set_size, inp, out);
		quick_sort_set(ssParams, ssParams->scatter_set, ssParams->scatter_set_size, 'c');
	}

	double **P = (double **)malloc(m * sizeof(double *));
	for (int i = 0; i < m; ++i){
		P[i] = (double *)malloc( b * sizeof(double *) );
	}

	double **tmp_row = (double **)malloc(1 * sizeof(double *));
			 tmp_row[0] = (double *)malloc(n * sizeof(double));
	double *msp = (double *)malloc(m * sizeof(double));

	int max_index;
	int min_index;
	for (int k = h; k < b; ++k, --m)
	{
		compute_Mt(ssParams, ssParams->ref_set, ssParams->ref_set_size, M, n, k);
		// print_double_matrix(ssParams, M, n, k);
		for (int i = 0; i < m; ++i)
		{
			for (int j = 0; j < n; ++j)
			{
				tmp_row[0][j] = ssParams->ref_set->members[0].params[j] - ssParams->scatter_set->members[i].params[j];
			}
			matrix_product(ssParams, tmp_row, 1, n, M, n, k, &(P[i]), 1, k);
			msp[i] = max(P[i], k, &max_index);
			// printf("%d\n", max_index);

		}
		min(msp, m, &min_index);
// printf("-----%d", min_index);
		evaluate_ind(ssParams, &(ssParams->scatter_set->members[min_index]), inp, out);
// printf("hi\n");
		copy_ind(ssParams, &(ssParams->ref_set->members[k]), &(ssParams->scatter_set->members[min_index]));
		delete_and_shift(ssParams, ssParams->scatter_set, ssParams->scatter_set_size, min_index);

		#ifdef STATS
			update_frequency_matrix(ssParams, &(ssParams->ref_set->members[k]));
		#endif
	}

	// print_set(ssParams, ssParams->ref_set, ssParams->ref_set_size, ssParams->nreal);
	// TODO: Free arrays

}

void compute_Mt(SSType *ssParams, Set *set, int set_size, double **M, int m_row, int m_col){
	for (int i = 0; i < m_row; ++i){		
		for (int j = 0; j < m_col; ++j)
		{
			M[i][j] = ssParams->ref_set->members[0].params[i] - ssParams->ref_set->members[j + 1].params[i];
		}
	}
}

/**
 * Generate set of subsets for recombination process using the members of set
 */
void select_subsets_list(SSType *ssParams, Set *set, int set_size){

	int k = 0;
	for (int i = 0; i < ssParams->ref_set_size; ++i)
	{
		for (int j = i + 1; j < ssParams->ref_set_size; ++j)
		{

			if( !is_equal(ssParams, &(set->members[i]), &(set->members[j])) 
				 && !is_exist_in_subsets_list(ssParams, &(set->members[j]), &(set->members[i])) ){
				// TODO: The candidateSet generation could be formed here to save some computation time.
				// 		 but the code will lose its clarity.
				
				copy_ind(ssParams, &(ssParams->subsets_list[k].members[0]), &(set->members[i]));

				copy_ind(ssParams, &(ssParams->subsets_list[k].members[1]), &(set->members[j]));

				k++;
				ssParams->subsets_list_size = k;
			}
		}
	}
}

