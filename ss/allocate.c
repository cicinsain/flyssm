#include "ss.h"

/*
	Allocate memory for an individual type
 */
void allocate_ind_memory(SSType *ssParams, individual *ind, int member_length){

	ind->params = (double *)malloc(member_length * sizeof(double));
	ind->cost = 0;
	ind->distance = 0;
}

void allocate_subset_memory(SSType *ssParams, individual *pair){

	
}

/*
	Allocate for set of individuals
 */
void allocate_set_memory(SSType *ssParams, Set *set, int set_size, int member_length){

	set->members = (individual *)malloc( set_size * sizeof(individual) );
	for (int i = 0; i < set_size; ++i)
	{
		allocate_ind_memory(ssParams, &(set->members[i]), member_length);
	}

}

/*
	Free memory of an individual
 */
void deallocate_ind_memory(SSType *ssParams, individual *ind){

	free(ind->params);
}

/*
	Free memory of set of individuals
 */
void deallocate_set_memory(SSType *ssParams, Set *set, int set_size){

	for (int i = 0; i < set_size; ++i)
	{
		deallocate_ind_memory(ssParams, &(set->members[i]));
	}
	free(set->members);

}

/*
	Free the memory allocated with sub_sets_list
 */
void deallocate_subsets_list_memory(SSType *ssParams){
						// max size of sub_sets_list
	for (int i = 0; i < ssParams->ref_set_size * ssParams->ref_set_size; ++i)
	{
		deallocate_set_memory(ssParams, &(ssParams->subsets_list[i]), ssParams->pair_size);
	}
}


/*
	Free all the memory allocated witb ssParams;
 */
void deallocate_ssParam(SSType *ssParams){

	free(ssParams->min_real_var);
	free(ssParams->max_real_var);

	for (int i = 0; i < ssParams->nreal; ++i){
		free(ssParams->freqs_matrix[i]);
		free(ssParams->probs_matrix[i]);
		free(ssParams->min_boundary_matrix[i]);
		free(ssParams->max_boundary_matrix[i]);
	}
	free(ssParams->freqs_matrix);
	free(ssParams->probs_matrix);
	free(ssParams->min_boundary_matrix);
	free(ssParams->max_boundary_matrix);

	deallocate_set_memory(ssParams, ssParams->ref_set, ssParams->ref_set_size);
	free(ssParams->ref_set);
	deallocate_set_memory(ssParams, ssParams->candidates_set, ssParams->ref_set_size * ssParams->ref_set_size * 4);
	free(ssParams->candidates_set);
	if ( !ssParams->perform_warm_start )
	{
		deallocate_set_memory(ssParams, ssParams->diverse_set, ssParams->diverse_set_size);
		free(ssParams->diverse_set);
	}
	deallocate_subsets_list_memory(ssParams);
	free(ssParams->subsets_list);


}