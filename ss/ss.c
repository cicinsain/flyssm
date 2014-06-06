#include "ss.h"

FILE *ref_set_history_file;
FILE *best_sols_history_file;
FILE *freqs_matrix_file;
FILE *freq_mat_final_file;
FILE *prob_mat_final_file;
FILE *ref_set_final_file;
FILE *stats_file;


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

	printf("\nInitializing Scatter Search...\n");

	// Preparing the output files
	init_report_files(ssParams);

	// it generates sample ssParams for test functions, some problem specific parameters should be initialize here.
	// init_sample_params(ssParams);		

	// Allocate memory of ssParams variables, and initialize some parameters
	init_ssParams(ssParams);

	if ( !ssParams->perform_warm_start ){

		init_scatter_set(ssParams, ssParams->scatter_set);
		evaluate_set(ssParams, ssParams->scatter_set, ssParams->scatter_set_size, inp, &out);

		init_ref_set(ssParams);
		quick_sort_set(ssParams, ssParams->ref_set, ssParams->ref_set_size, 'c');		
															
		// Initialize the best solution
		ssParams->best = (individual *)malloc(sizeof(individual));

		// Assigning the best solutions
		
		ssParams->best = &(ssParams->ref_set->members[0]);
	}
	else{
		init_warm_start(ssParams);
		// FIX: Sometimes crashes...
		// re_gen_ref_set(ssParams, ssParams->ref_set, ssParams->ref_set_size, 's');
		// ssParams->n_regen++;

	}	

	write_set(ssParams, ssParams->ref_set, ssParams->ref_set_size, ssParams->nreal, ref_set_history_file, 0, 'w');
	write_ind(ssParams, ssParams->best, ssParams->nreal, best_sols_history_file, 0, 'w');


}


/*
	Run the Scatter Search process
 */
void RunSS(Input *inp, SSType *ssParams, char *inname){

	int iter;
	int n_ref_set_update    = 0;
	int n_refinement        = 0;
	int n_duplicates        = 0;
	int n_function_evals    = 0;
	int n_flatzone_detected = 0;
	// bool wasChanged = false;

	// TODO: Add new criteria for the stop, maybe minimum distance to the minumum result
	printf("Starting the optimization procedure...\n");
	for (iter = 1; iter < ssParams->max_iter; ++iter)
	{
		// printf("hi\n");
		// Selecting the SubSets List
		select_subsets_list(ssParams, ssParams->ref_set, ssParams->ref_set_size);
		// Generate new candidates
		generate_candiates(ssParams);
		evaluate_set(ssParams, ssParams->candidates_set, ssParams->candidates_set_size, inp, &out);
		// if (ssParams->perform_local_search && (iter % ssParams->local_search_freq == 0)  ){
		// 	printf("%s", KGRN);
		// 	printf("-");
		// 	printf("%s", KNRM);
		// 	refine_set(ssParams, ssParams->candidates_set, ssParams->candidates_set_size, 's', inp, &out);
		// }
		// Update refSet by replacing new cadidates
		update_ref_set(ssParams);
		// print_set(ssParams, ssParams->ref_set, ssParams->ref_set_size, ssParams->nreal);
		// print_set(ssParams, ssParams->candidates_set, ssParams->candidates_set_size, ssParams->nreal);

		// Perform the local_search
		if (ssParams->perform_local_search && (iter % ssParams->local_search_freq == 0)  ){
			refine_set(ssParams, ssParams->ref_set, ssParams->ref_set_size, 's', inp, &out);
		}
		quick_sort_set(ssParams, ssParams->ref_set, ssParams->ref_set_size, 'c');
		
		// Append the ref_set to the file
		write_set(ssParams, ssParams->ref_set, ssParams->ref_set_size, ssParams->nreal, ref_set_history_file, iter, 'w');
		fflush(ref_set_history_file);

		// Append the best solution to the sol history file
		write_ind(ssParams, ssParams->best, ssParams->nreal, best_sols_history_file, iter, 'w');
		fflush(best_sols_history_file);

		// loadBar(iter, ssParams->max_iter, 50, 50);

		if (ssParams->perform_stop_criteria)
			if (fabs(ssParams->ref_set->members[0].cost - ssParams->ref_set->members[ssParams->ref_set_size - 1].cost) < ssParams->stop_criteria){
				printf("\n%s   Stop by difference criteria!\n   The difference between the best and worst memebers of refSet is smaller than %lf\n\n%s", KRED, ssParams->stop_criteria, KNRM);
				break;
			}

		if ((ssParams->perform_ref_set_regen && iter != 1) &&
			(
			 (double)(ssParams->n_duplicates - n_duplicates) / (double)ssParams->candidates_set_size > 0.7
			 // || (double)(ssParams->n_flatzone_detected - n_flatzone_detected) / (double)ssParams->candidates_set_size > 0.7
			 || iter % ssParams->ref_set_regen_freq == 0)
			)
		{
			// if (rndreal(0, 1) < 0.5)
			// 	re_gen_ref_set(ssParams, ssParams->ref_set, ssParams->ref_set_size, 's', inp, &out);
			// else
				re_gen_ref_set(ssParams, ssParams->ref_set, ssParams->ref_set_size, 'n', inp, &out);

			// print_set(ssParams, ssParams->ref_set, ssParams->ref_set_size, ssParams->nreal);

			ssParams->n_regen++;
			quick_sort_set(ssParams, ssParams->ref_set, ssParams->ref_set_size, 'c');
			printf("%s", KBLU);
			printf("+");
			printf("%s", KNRM);

		}

		// printf("hi\n");

		#ifdef STATS
			write_int_matrix(ssParams, ssParams->freqs_matrix, ssParams->nreal, ssParams->p, freqs_matrix_file, iter, 'w');
		#endif		

		ssParams->n_iter++;
		fprintf(stats_file, "%d\t%d\t%d\t%d\t%d\t%d\t%d\n", iter, ssParams->n_ref_set_update -  n_ref_set_update, ssParams->n_flatzone_detected -  n_flatzone_detected, ssParams->n_refinement -  n_refinement, ssParams->n_duplicates -  n_duplicates, ssParams->n_function_evals -  n_function_evals, ssParams->candidates_set_size);

		#ifdef DEBUG
			printf("\nStats - (%d):\n", iter);
			// print_ind(ssParams, ssParams->best, ssParams->nreal);	
			printf("\t\tBest RMS: %lf\n", ssParams->best->cost);
			printf("\t\t# Replacement in Reference Set: %d\n", ssParams->n_ref_set_update - n_ref_set_update);
			printf("\t\t# of Local Search Performed: %d\n", ssParams->n_refinement - n_refinement);
			printf("\t\t# Duplicates: %d\n", ssParams->n_duplicates - n_duplicates);
			printf("\t\t# Flatzone: %d\n", ssParams->n_flatzone_detected - n_flatzone_detected);
			printf("\t\t================= candidateSetSize: %d\n", ssParams->candidates_set_size);
		#endif	

		n_ref_set_update    = ssParams->n_ref_set_update;
		n_refinement        = ssParams->n_refinement;
		n_duplicates        = ssParams->n_duplicates;
		n_function_evals    = ssParams->n_function_evals;
		n_flatzone_detected = ssParams->n_flatzone_detected;

	}

	printf("%s", KYEL);
	printf("\nReference Set:\n");
	print_set(ssParams, ssParams->ref_set, ssParams->ref_set_size, ssParams->nreal);
	printf("%s", KNRM);
	
	printf("%s", KCYN);
	printf("\n====================================\n");
	printf("Best Solution:\n");
	print_ind(ssParams, ssParams->best, ssParams->nreal);
	printf("====================================\n");
	printf("%s", KNRM);

	printf("%s", KGRN);
	printf("\nStatistics: \n");
	printf("# of iterations: %d\n", iter);
	printf("# Replacement in Reference Set: %d\n", ssParams->n_ref_set_update);
	printf("# Flatzone detected: %d\n", ssParams->n_flatzone_detected);
	printf("# Local Search Performed: %d\n", ssParams->n_refinement);
	printf("# Duplicates: %d\n", ssParams->n_duplicates);
	printf("# Function Evalation: %d\n", ssParams->n_function_evals);
	printf("# refSet regen: %d\n", ssParams->n_regen);
	printf("%s", KNRM);

	freq_mat_final_file = fopen("freq_mat_final.csv", "w");
	prob_mat_final_file = fopen("prob_mat_final.csv", "w");
	ref_set_final_file = fopen("ref_set_final.csv", "w");


	printf("\nExporting ref_set_final.csv, freq_mat_final.csv and prob_mat_final.csv for warm start...\n");
	write_set(ssParams, ssParams->ref_set, ssParams->ref_set_size, ssParams->nreal, ref_set_final_file, -1, 'w');
	write_int_matrix(ssParams, ssParams->freqs_matrix, ssParams->nreal, ssParams->p, freq_mat_final_file, -1, 'w');
	write_double_matrix(ssParams, ssParams->probs_matrix, ssParams->nreal, ssParams->p, prob_mat_final_file, -1, 'w');


	deallocate_ssParam(ssParams);

	fclose(ref_set_history_file);
	fclose(best_sols_history_file);

}
