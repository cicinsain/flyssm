#include "ss.h"
#include <string.h>


void refine_subsets_list(SSType *ssParams, char method, Input *inp, ScoreOutput *out){

	for (int i = 0; i < ssParams->subsets_list_size; ++i)
	{
		refine_set(ssParams, &(ssParams->subsets_list[i]), ssParams->pair_size, method, inp, out);
	}

	// deallocate_ind_memory(ssParams, new_candidate);
}

void refine_set(SSType *ssParams, Set *set, int set_size, char method, Input *inp, ScoreOutput *out){

	for (int i = set_size - 1; i >= 0; --i)
	{
		/* The local search won't apply on members with bad fitness; larger than ssParams->local_search_f1_criteria * sol */
		// FIXME: There is a bug here that cause `nan` value as a cost! for some functions!
		// printf("%lf < %lf\n", fabs(set->members[i].cost - ssParams->sol), ssParams->local_search_f1_criteria);
		if (ssParams->local_search_1_filter){
			if ( !( ( fabs(set->members[i].cost - ssParams->sol)) > ssParams->local_search_f1_criteria) )
		// if(0)
			{
				// printf("0\n");
			// If the first filter pass, then we check to see if the area already discovered with
			// local search procedure or not, by locating the clsest members in refSet with the 
			// selected memeber, if their fitness don't vary a lot then we perform the local search.
				second_filter:
				if (ssParams->local_search_2_filter){
					if ( ! ((fabs(set->members[i].cost - 
						set->members[ closest_member(ssParams, set, set_size, &(set->members[i]), i) ].cost)) 
						< ssParams->local_search_f2_criteria) )
					{
						local_search:
						ssParams->n_refinement++;
						refine_individual(ssParams, set, set_size, &(set->members[i]), method, inp, out);
						// print_ind(ssParams, &(set->members[i]), ssParams->nreal);
					}
				}else{
					goto local_search;
				}
			}
		}
		else{
			goto second_filter;
		}
	}
}


void refine_individual(SSType *ssParams, Set *set, int set_size, individual *ind, char method, Input *inp, ScoreOutput *out){

	// TODO: Make these temporary variables global to improve the performance
 	// new_candidate = (individual *)malloc(sizeof(individual));
	// allocate_ind_memory(ssParams, &(new_candidate), ssParams->nreal);

individual new_candidate;


	double *new_params = (double *)malloc( ssParams->nreal * sizeof(double));

	for (int i = 0; i < ssParams->max_no_improve; ++i)
	{
		/* Run the local optimization procedure */
		switch (method){
			case 't':		// Stochastic Hill Climbing 
				{
					ssParams->n_refinement++;
					take_step(ssParams, ind->params, new_params);
					new_candidate.params = new_params;
					evaluate_ind(ssParams, &( new_candidate), inp, out);
			
					if (new_candidate.cost < ind->cost)
					{
						/* Replace ind->params with newly generated params */
						copy_ind(ssParams, ind, &( new_candidate));
					}
					break;
				}
			case 's':	// SMART!
				{

					take_step(ssParams, ind->params, new_params);
					new_candidate.params = new_params;
					evaluate_ind(ssParams, &( new_candidate), inp, out);
			
					if (new_candidate.cost < ind->cost)
					{
						/* Replace ind->params with newly generated params */
						copy_ind(ssParams, ind, &( new_candidate));
						// TODO: the matrix shouldn't update every iteration, should define a flag here
						//		 and check if it perform then at the end of the routine update the matrix
						#ifdef STATS
							update_frequency_matrix(ssParams, &(new_candidate));
						#endif
					}
					break;	
				}
			}
	}

	free(new_params);

	// deallocate_ind_memory(ssParams, new_candidate);
}











