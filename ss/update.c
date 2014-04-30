#include "ss.h"
#include <string.h>




/*
	The refinment procedure, based on the routine describes in Banga,
	it only replace the refSet members with better candidates by considering the fitness
	and the diversity (equality) to avoid duplicate. The refinement (performing local search)
	will be executed seperately.
 */
// TODO: This could be moved to the select.c and refine_* function from local_search could be moved here.
void update_ref_set(SSType *ssParams){

	// Candidates set is sorted!
	// Case 1: One candidate beats the best solutions already found
	
	int i = 0;
	if (ssParams->candidates_set->members[i].cost < ssParams->ref_set->members[i].cost )
	{
			copy_ind(ssParams, &(ssParams->ref_set->members[i]), &(ssParams->candidates_set->members[i]));
			i++;
	}


	// The loop continues until there is not better result in candidate set
	while (ssParams->candidates_set->members[i].cost < ssParams->ref_set->members[ ssParams->ref_set_size - 1].cost )
	{
			/* 
				Check the existance of the candidate in the ref_set
				This avoid the presence of duplicate solutions, as a consequence there will be 
				only one global minima in the solution set, and the points around it will be
				in the ref_set as well..
			 */
			// We send the `ref_set_size - 1` since we want to avoid the last memeber.
			int duplicate_index = is_exist(ssParams, ssParams->ref_set, ssParams->ref_set_size - 1, &(ssParams->candidates_set->members[i]));
			if (-1 == duplicate_index)
			{
				/* Check if the new candidates is located in the flatzone, but it's bit tricky,  
				   avoid many duplicate in the scatter search; consequently, less risk to stuck in
				   a local minima somewhere
				   */
				// if ( !is_in_flatzone(ssParams, ssParams->ref_set , ssParams->ref_set_size - 1, &(ssParams->candidates_set->members[i])) )
				// if ( ssParams->candidates_set->members[i].cost < ssParams->ref_set->members[ssParams->ref_set_size -  1].cost * ( 1 - ssParams->fitness_epsilon) )
				{
					copy_ind(ssParams, &(ssParams->ref_set->members[ssParams->ref_set_size - 1]), &(ssParams->candidates_set->members[i]));
					insertion_sort(ssParams, ssParams->ref_set, ssParams->ref_set_size, 'c');

					ssParams->n_ref_set_update++;

					#ifdef STATS
						update_frequency_matrix(ssParams, &(ssParams->candidates_set->members[i]));
					#endif					
				}
			}
			else{
				/* Check if the candidate has better fitness than the duplicated member; or at least it has reasonable 
					differences with it (passing the flatzone filter). If yes, then we could replace it.
				*/
					if ( ssParams->candidates_set->members[i].cost < ssParams->ref_set->members[duplicate_index].cost * ( 1 - ssParams->fitness_epsilon) )
					{
						copy_ind(ssParams, &(ssParams->ref_set->members[duplicate_index]), &(ssParams->candidates_set->members[i]));
						ssParams->n_ref_set_update++;
						// TODO: BubbleSort or QuickSort needed! BubbleSort is perfect since the list is already sorted
						// and one run of BubbleSort will put the replaced member in the right place.
						
						#ifdef STATS
							update_frequency_matrix(ssParams, &(ssParams->candidates_set->members[i]));
						#endif	
					}

				// }
			}
		i++;
	}

}

/*
	Generate set of subsets for recombination process using the members of set
 */
// TODO: Need to be improved
void select_subsets_list(SSType *ssParams, Set *set, int set_size){

	// printf("%d\n", ssParams->subsets_list_size);

	int k = 0;
	// printf("%d\n", ssParams->ref_set_size);
	for (int i = 0; i < ssParams->ref_set_size; ++i)
	{
		for (int j = i + 1; j < ssParams->ref_set_size; ++j)
		{

			if( !is_equal(ssParams, &(set->members[i]), &(set->members[j])) 
				/* && is_exist(ssParams, &(set->members[j]), &(set->members[i])  */ ){
				/*TODO*/
				// TODO: The candidateSet generation could be formed here to save some computation time.
				// 		 but the code will lose its clarity.
				
				copy_ind(ssParams, &(ssParams->subsets_list[k].members[0]), &(set->members[i]));

				copy_ind(ssParams, &(ssParams->subsets_list[k].members[1]), &(set->members[j]));

				k++;
			}
		}
	}
	// printf("%d\n", k);
	ssParams->subsets_list_size = k;

}

