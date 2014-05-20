#include "ss.h"


/*
	Allocate memory for ssParams
 */
void init_ssParams(SSType *ssParams){


	// Estimating the Reference Set Size
	ssParams->subsets_list_size = (ssParams->ref_set_size * ssParams->ref_set_size);		// Maximum possible size of subsets_list
	
	ssParams->n_refinement        = 0;
	ssParams->n_ref_set_update    = 0;
	ssParams->n_duplicates        = 0;
	ssParams->n_flatzone_detected = 0;
	ssParams->n_function_evals    = 0;

	ssParams->freqs_matrix = (int **)malloc(ssParams->nreal * sizeof(int *));
	ssParams->probs_matrix = (double **)malloc(ssParams->nreal * sizeof(double *));

	// Generating the sub regions matrixes
	ssParams->min_boundary_matrix = (double **)malloc( ssParams->nreal * sizeof(double *));
	ssParams->max_boundary_matrix = (double **)malloc( ssParams->nreal * sizeof(double *));

	for (int i = 0; i < ssParams->nreal; ++i)
	{
		ssParams->min_boundary_matrix[i] = (double *)malloc((ssParams->p) * sizeof(double));
		ssParams->max_boundary_matrix[i] = (double *)malloc((ssParams->p) * sizeof(double));

		ssParams->freqs_matrix[i] = (int *)malloc(ssParams->p * sizeof(int));
		ssParams->probs_matrix[i] = (double *)malloc(ssParams->p * sizeof(double));

		for (int j = 1; j <= ssParams->p; ++j)
		{
			/* Building the sub region matrixes */	// One matrix would be enough but for clarity I use two.
			ssParams->min_boundary_matrix[i][j - 1] = ssParams->min_real_var[i] + ((ssParams->max_real_var[i] - ssParams->min_real_var[i]) / ssParams->p) * (j - 1);
			ssParams->max_boundary_matrix[i][j - 1] = ssParams->min_real_var[i] + ((ssParams->max_real_var[i] - ssParams->min_real_var[i]) / ssParams->p) * j; 

			/* Building the freqs_matrix: frequently_matrix */
			ssParams->freqs_matrix[i][j - 1] = 1;
			ssParams->probs_matrix[i][j - 1] = 1./(ssParams->p);		// Since all values are 1.
		}
	}

	#ifdef STATS
		write_int_matrix(ssParams, ssParams->freqs_matrix, ssParams->nreal, ssParams->p, freqs_matrix_file, 0, 'w');
	#endif

	// print_double_matrix(ssParams, ssParams->min_boundary_matrix, ssParams->nreal, ssParams->p);
	// print_double_matrix(ssParams, ssParams->max_boundary_matrix, ssParams->nreal, ssParams->p);

}

/*
	Generate Scatter Set based on the Banga's instruction.
 */
void init_scatter_set(SSType *ssParams, Set *set){

	// Generating Scatter Set with the size of p
	for (int k = 0; k < ssParams->p; ++k)
	{
		for (int i = 0; i < ssParams->nreal; ++i)
		{
			set->members[k].params[i] = rndreal(ssParams->min_boundary_matrix[i][k], ssParams->max_boundary_matrix[i][k]);
		}
	}

	// Expanding the Scatter Set to the size of 10*ref_set_size

	int a = 0;					// a indicated the index of selected sub-regions that the new value should be generated from.
	double probs_sum = 0;
	double rnd = 0;

	for (int k = 0; k < ssParams->diverse_set_size - ssParams->p; ++k)
	{	
		/* Generate new item to add to scatter_set
		   set->members[ssParams->p + k] is going to be generated */
		
		for (int i = 0; i < ssParams->nreal; ++i)
		{	
			rnd = rndreal(0, 1);

			for (int j = 0; j < ssParams->p; ++j)
			{	
				/* Compute `probs` and find the index */

				probs_sum += ssParams->probs_matrix[i][j];

				if ( rnd <= probs_sum ){
					a = j;
					
					// Updating the f matrix, and prob_matrix
					ssParams->freqs_matrix[i][a]++;
					
					// Updating the ssParams->probs_matrix
					// IMPROVE: This could be rewritten much nicer
					double f_col_prob = 0;
					for (int t = 0; t < ssParams->p; ++t){
						f_col_prob += (1. / ssParams->freqs_matrix[i][t]);
					}
					ssParams->probs_matrix[i][a] = (1. / ssParams->freqs_matrix[i][a] ) / (f_col_prob);

					probs_sum = 0;
					break;
				}
			}

			// sub-region index is selected now: `a`
			set->members[ssParams->p + k].params[i] = rndreal(ssParams->min_boundary_matrix[i][a], ssParams->max_boundary_matrix[i][a] );

		}

	}	// Scatter Set is now extented to 10*b (b: ref_set_size)

	#ifdef STATS
		write_int_matrix(ssParams, ssParams->freqs_matrix, ssParams->nreal, ssParams->p, freqs_matrix_file, 0, 'w');
	#endif

	//TODO: Free stuffs
}


/*
	Form the ref_set using the procedure describes in Banga paper
	It assumes that the ref_set is initialize and the proper memory allocated for it.
 */
void init_ref_set(SSType *ssParams){

	printf("Forming the refSet...\n");

	int m = ssParams->diverse_set_size;
	int b = ssParams->ref_set_size;
	int h = b / 2;	// max_elite

	// Sorting the diverse_set
	quick_sort_set(ssParams, ssParams->diverse_set, ssParams->diverse_set_size, 'c');

	for (int i = 0; i < h; ++i)
	{
		/* Copying first `h` best items of the diverse_set with respect to the cost */
		copy_ind(ssParams, &(ssParams->ref_set->members[i]), &(ssParams->diverse_set->members[i]));
	}

	// Expanding the ref_set by applying the diversity distance
	
	// dist_matrix will store the distances between the rest of diverse_set and the ref_set
	double **dist_matrix;
	dist_matrix = (double **)malloc( (m - h) * sizeof(double *));
	for (int i = 0; i < (m - h); ++i){
		dist_matrix[i] = (double *)malloc( (b) * sizeof(double));
	}

	double *min_dists = (double *)malloc((m - h) * sizeof(double));


	int min_index, max_index;
	min_index = 0;
	max_index = 0;

	int i = 0;

	for (int k = h; k < b; ++k, --m)		// Iterate through the ref_set for expanding
	{										// --m is reduce the depth that the second `for` should go through it
											// since the array is shifted after each iteration

		for (i = h; i < m; ++i)		// Iterate through the rest of diverse_set
		{							 

			for (int j = 0; j < k; ++j)		// Iterate trhough the selected members of ref_set
			{								
				/* Computing the distance */
				dist_matrix[i - (b / 2)][j] = euclidean_distance(ssParams, &(ssParams->diverse_set->members[i]), &(ssParams->ref_set->members[j]));
			}

			// Find the minimum of each row into a new variable, then finding it's max
			min_dists[i - (b / 2)] = min(dist_matrix[i - (b / 2)], k, &min_index);
		}
		
		// Now, the min_dists of vectors to all the members of R is computed,
		// we could compute the maximum of them and also it's index
		
		// Finding the place of the max/min distance
		max(min_dists, i - (b / 2), &max_index); 

		// Adding the vector in max_index to the ref_set
		copy_ind(ssParams, &(ssParams->ref_set->members[k]), &(ssParams->diverse_set->members[h + max_index]));

		// After assigment, the rest of should be shift into the left form the position that the 
		// vector is selected
		
		// It could be done more complicated!
		delete_and_shift(ssParams, ssParams->diverse_set, ssParams->diverse_set_size, (b / 2) + max_index);

	}

	#ifdef STATS
		write_int_matrix(ssParams, ssParams->freqs_matrix, ssParams->nreal, ssParams->p, freqs_matrix_file, 0, 'w');
	#endif

	free(min_dists);
	for (int i = 0; i < (ssParams->diverse_set_size - b/2); ++i){
		free(dist_matrix[i]);
	}
	free(dist_matrix);

}


