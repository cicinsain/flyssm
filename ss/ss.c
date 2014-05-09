#include "ss.h"

FILE *ref_set_history_file;
FILE *best_sols_history_file;
FILE *freqs_matrix_file;
FILE *freq_mat_final_file;
FILE *prob_mat_final_file;
FILE *ref_set_final_file;

ScoreOutput out;


/*
	Initialize the Scatter Search process
 */
void InitSS(Input *inp, SSType *ssParams, char *inname){


    out.score          = 1e38;           // start with a very large number
    out.penalty        = 0;
    out.size_resid_arr = 0;
    out.jacobian       = NULL;
    out.residuals      = NULL;

	// TODO: To be decide based on the parameters input file
	// char *ref_set_history_file_name;
	// char *best_sols_history_file_name;

	printf("\nInitialize Scatter Search...\n");

	ref_set_history_file = fopen("ref_set_history_file.out", "w");
	best_sols_history_file = fopen("best_sols_history_file.out", "w");
	freqs_matrix_file = fopen("freqs_matrix_history.out", "w");


	// FIXME: I am hardcoded!

	// Preparing the output files
	init_report_files(ssParams);

	// it generates sample ssParams for test functions
	// init_sample_params(ssParams);		

	// Allocate memory of ssParams variables
	init_ssParams(ssParams);

	if ( !ssParams->perform_warm_start )
	{
		// Initialize the Diverse Set
		ssParams->diverse_set = (Set *)malloc(sizeof(Set));
		allocate_set_memory(ssParams, ssParams->diverse_set, ssParams->diverse_set_size, ssParams->nreal);

		init_scatter_set(ssParams, ssParams->diverse_set);				// Banga version
		evaluate_set(ssParams, ssParams->diverse_set, ssParams->diverse_set_size, inp, &out);

	// print_set(ssParams, ssParams->diverse_set, ssParams->diverse_set_size, ssParams->nreal);

		// Initialize the Reference Set
		ssParams->ref_set = (Set *)malloc(sizeof(Set));
		allocate_set_memory(ssParams, ssParams->ref_set, ssParams->ref_set_size, ssParams->nreal);

		/* Expanding the ref_set by appending solution with least cost
			and then adding solutions with least distance from max_elite
		 	selected solutions in ref_set */
		// It should be sorted since the item added based on the distance are not in the fitness order
		init_ref_set(ssParams);		// Banga version
		quick_sort_set(ssParams, ssParams->ref_set, ssParams->ref_set_size, 'c');		
															
		// Writing the sorted ref_set to the file
		
		// TODO: Now, the diverse set could be freed

		// Initialize the best solution
		ssParams->best = (individual *)malloc(sizeof(individual));
		// allocate_ind_memory(ssParams, ssParams->best, ssParams->nreal);	

		// Assigning the best solutions
		ssParams->best = &(ssParams->ref_set->members[0]);				// The first members of ref_set is always the best

	}else{	// WARM START
		warm_start(ssParams);
	}

#ifdef LOG
	write_set(ssParams, ssParams->ref_set, ssParams->ref_set_size, ssParams->nreal, ref_set_history_file, 0, 'w');
	write_ind(ssParams, ssParams->best, ssParams->nreal, best_sols_history_file, 0, 'w');
#endif

}


/*
	Run the Scatter Search process
 */
void RunSS(Input *inp, SSType *ssParams, char *inname){

	bool wasChanged = false;

	// Initialize the SubSets List
	ssParams->subsets_list = (Set *)malloc( ssParams->subsets_list_size * sizeof(Set) );
	for (int i = 0; i < ssParams->subsets_list_size; ++i)
	{
		/* 	allocating memory for each subset, with the size equal to `pair_size` and lenth of `nreal`. */
		allocate_set_memory(ssParams, &(ssParams->subsets_list[i]), ssParams->pair_size, ssParams->nreal);
	}

	// Declare the candidate set with the maximum size possible!
	ssParams->candidates_set = (Set *)malloc(sizeof(Set));
	allocate_set_memory(ssParams, ssParams->candidates_set, ssParams->ref_set_size * ssParams->ref_set_size * 6, ssParams->nreal);


	// TODO: Add new criteria for the stop, maybe minimum distance to the minumum result
	printf("Starting the optimization procedure...\n");
	for (int iter = 1; iter < ssParams->max_iter; ++iter)
	{

		// Selecting the SubSets List
		select_subsets_list(ssParams, ssParams->ref_set, ssParams->ref_set_size);


		/* --------------- Banga Specific ------------ */
		// Generate new candidates
		generate_candiates(ssParams);
		evaluate_set(ssParams, ssParams->candidates_set, ssParams->candidates_set_size, inp, &out);
		quick_sort_set(ssParams, ssParams->candidates_set, ssParams->candidates_set_size, 'c');

	
		// Update refSet by replacing new cadidates
		update_ref_set(ssParams);

		// Perform the local_search
		if (ssParams->perform_local_search ){
			refine_set(ssParams, ssParams->ref_set, ssParams->ref_set_size, 's', inp, &out);
		}
		quick_sort_set(ssParams, ssParams->ref_set, ssParams->ref_set_size, 'c');

#ifdef LOG
		// Append the ref_set to the file
		write_set(ssParams, ssParams->ref_set, ssParams->ref_set_size, ssParams->nreal, ref_set_history_file, iter, 'w');
		fflush(ref_set_history_file);

		// Append the best solution to the sol history file
		write_ind(ssParams, ssParams->best, ssParams->nreal, best_sols_history_file, iter, 'w');
		fflush(best_sols_history_file);
#endif

		// print_set(ssParams, ssParams->ref_set, ssParams->ref_set_size, ssParams->nreal);
		loadBar(iter, ssParams->max_iter, 50, 50);

		#ifdef STATS
			write_int_matrix(ssParams, ssParams->freqs_matrix, ssParams->nreal, ssParams->p, freqs_matrix_file, iter, 'w');
		#endif		

		printf("\nBest Solution:\n");
		// print_ind(ssParams, ssParams->best, ssParams->nreal);	
		printf("rms:%lf\n", ssParams->best->cost);
		printf("Number of substitution in Reference Set: %d\n", ssParams->n_ref_set_update);
		printf("Number of Local Search Performed: %d\n", ssParams->n_refinement);

	}

	write_params_to_fly_output_standard(ssParams, inp, inname);

	// printf("Reference Set:\n");
	// print_set(ssParams, ssParams->ref_set, ssParams->ref_set_size, ssParams->nreal);

	printf("\n====================================\n");
	printf("Best Solution:\n");
	print_ind(ssParams, ssParams->best, ssParams->nreal);

	printf("====================================\n");
	printf("\nStatistics: \n");
	printf("Number of substitution in Reference Set: %d\n", ssParams->n_ref_set_update);
	printf("Number of Local Search Performed: %d\n", ssParams->n_refinement);


	printf("\nExporting ref_set_final.csv, freq_mat_final.csv and prob_mat_final.csv for warm start...\n");

	freq_mat_final_file = fopen("freq_mat_final.csv", "w");
	prob_mat_final_file = fopen("prob_mat_final.csv", "w");
	ref_set_final_file = fopen("ref_set_final.csv", "w");

	write_set(ssParams, ssParams->ref_set, ssParams->ref_set_size, ssParams->nreal, ref_set_final_file, -1, 'w');
	write_int_matrix(ssParams, ssParams->freqs_matrix, ssParams->nreal, ssParams->p, freq_mat_final_file, -1, 'w');
	write_double_matrix(ssParams, ssParams->probs_matrix, ssParams->nreal, ssParams->p, prob_mat_final_file, -1, 'w');


	deallocate_ssParam(ssParams);

	fclose(ref_set_history_file);
	fclose(best_sols_history_file);


}
