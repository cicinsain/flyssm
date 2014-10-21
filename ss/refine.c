#include "ss.h"
#include <string.h>


void refine_subsets_list(SSType *ssParams, char method, Input *inp, ScoreOutput *out){

	for (int i = 0; i < ssParams->subsets_list_size; ++i)
	{
		refine_set(ssParams, &(ssParams->subsets_list[i]), ssParams->pair_size, method, inp, out);
	}
}

void refine_set(SSType *ssParams, Set *set, int set_size, char method, Input *inp, ScoreOutput *out){

	int closest_member_index;
	for (int i = 0; i < set_size; ++i)
	{
		if (ssParams->local_search_1_filter)
		{
			if ( ( fabs(set->members[i].cost - ssParams->sol)) < ssParams->local_search_f1_criteria) 
			{
				goto second_filter;

				second_filter:
				{
					if (ssParams->local_search_2_filter)
					{
						closest_member_index = closest_member(ssParams, set, set_size, &(set->members[i]), i);
						/**
						 * The idea was to check if the closest member is actually close enough to be a reason to stop the local search
						 * but I guess finding the reasonable value for that (which is '5') might be a bit hard!
						 */
						// printf("++%lf\n", euclidean_distance(ssParams, &(set->members[i]), &(set->members[closest_member_index])) );
						if ( !(euclidean_distance(ssParams, &(set->members[i]), &(set->members[closest_member_index])) < ssParams->min_distance_for_local_search) )
						{	
							// if ( ! ((fabs(set->members[i].cost - set->members[ closest_member_index ].cost)) < ssParams->local_search_f2_criteria) )
							if ( !( set->members[i].cost < set->members[closest_member_index].cost + (ssParams->local_search_f2_criteria * set->members[closest_member_index].cost )
									&& set->members[i].cost > set->members[closest_member_index].cost - (ssParams->local_search_f2_criteria * set->members[closest_member_index].cost ) ) )
								// Will check if the new candidate is in flatzone with radius of local_search_f2_criteria
							{
								// goto local_search;

								// local_search:
								// {
									refine_individual(ssParams, set, set_size, &(set->members[i]), method, inp, out);
								 	//ssParams->n_refinement++;
								// }
							}
						}
					}
					else
					{
						refine_individual(ssParams, set, set_size, &(set->members[i]), method, inp, out);
					}
				}
			}
		}
		else
		{
			goto second_filter;
		}
	}
}


void refine_individual(SSType *ssParams, Set *set, int set_size, individual *ind, char method, Input *inp, ScoreOutput *out){

	// TODO: Make these temporary variables global to improve the performance
 	// new_candidate = (individual *)malloc(sizeof(individual));
	// allocate_ind_memory(ssParams, &(new_candidate), ssParams->nreal);

individual new_candidate;
printf("refining individual\n");

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

	ssParams->n_refinement++;
        printf("nrefinement: %d\n", ssParams->n_refinement);
	free(new_params);

	// deallocate_ind_memory(ssParams, new_candidate);
}











