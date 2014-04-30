#include "ss.h"


/*
	Use the sub_sets_list members to produce candidates,
	the recombined_set will be generated.
 */
void generate_candiates(SSType *ssParams){

	int candidates_count = 0;
	double *steps = (double *)malloc( ssParams->nreal * sizeof(double));
	double mid_cost = ssParams->ref_set->members[ssParams->max_elite].cost;


	for (int i = 0; i < ssParams->subsets_list_size; ++i)
	{
		for (int j = 0; j < ssParams->nreal; ++j)
		{
			/* Step sizes in every direction */
			// TOOD: Check the formula
			steps[j] = rndreal(0, 1) * ( ssParams->subsets_list[i].members[0].params[j] -  ssParams->subsets_list[i].members[1].params[j] ) / 2;
		}

		/* Generating new candidates */
		if ( ssParams->subsets_list[i].members[0].cost < mid_cost && ssParams->subsets_list[i].members[1].cost < mid_cost )		// type 1
		{
			generate_ind_candidate(ssParams, &(ssParams->subsets_list[i].members[0]), &(ssParams->candidates_set->members[ candidates_count++ ]), steps, '1');
			
			generate_ind_candidate(ssParams, &(ssParams->subsets_list[i].members[0]), &(ssParams->candidates_set->members[ candidates_count++ ]), steps, '2');

			generate_ind_candidate(ssParams, &(ssParams->subsets_list[i].members[1]), &(ssParams->candidates_set->members[ candidates_count++ ]), steps, '3');
			
			generate_ind_candidate(ssParams, &(ssParams->subsets_list[i].members[1]), &(ssParams->candidates_set->members[ candidates_count++ ]), steps, '3');

		}else
		if ( ssParams->subsets_list[i].members[0].cost < mid_cost && ssParams->subsets_list[i].members[1].cost >= mid_cost )	// type 2
		{
			
			generate_ind_candidate(ssParams, &(ssParams->subsets_list[i].members[0]), &(ssParams->candidates_set->members[ candidates_count++ ]), steps, '1');

			generate_ind_candidate(ssParams, &(ssParams->subsets_list[i].members[0]), &(ssParams->candidates_set->members[ candidates_count++ ]), steps, '2');

			generate_ind_candidate(ssParams, &(ssParams->subsets_list[i].members[1]), &(ssParams->candidates_set->members[ candidates_count++ ]), steps, '3');

		}else
		if ( ssParams->subsets_list[i].members[0].cost >= mid_cost && ssParams->subsets_list[i].members[1].cost >= mid_cost )	// type 3
		{
			
			if (rndreal(0, 1) < 0.5)
			{
				generate_ind_candidate(ssParams, &(ssParams->subsets_list[i].members[0]), &(ssParams->candidates_set->members[ candidates_count++ ]), steps, '1');
			}else
			{
				generate_ind_candidate(ssParams, &(ssParams->subsets_list[i].members[1]), &(ssParams->candidates_set->members[ candidates_count++ ]), steps, '3');
			}

			generate_ind_candidate(ssParams, &(ssParams->subsets_list[i].members[0]), &(ssParams->candidates_set->members[ candidates_count++ ]), steps, '2');

		}


	}

	ssParams->candidates_set_size = candidates_count;

	free(steps);

}


void generate_ind_candidate(SSType *ssParams, individual *base, individual *candidate,  double *steps, char type){
	// FIXME: It happens that the variables go outside the boundaries, I could be solved by finding the MIN and MAX between new values and min and max of variables.

	int i;
	double new_value;
	switch (type){
		case '1':
			for (i = 0; i < ssParams->nreal; ++i)
			{
				new_value = base->params[i] - steps[i];
				new_value = MIN(new_value, ssParams->max_real_var[i]);
				new_value = MAX(new_value, ssParams->min_real_var[i]);
				candidate->params[i] = new_value;
			}
			break;
		case '2':
			for (i = 0; i < ssParams->nreal; ++i)
			{
				new_value = base->params[i] + steps[i];
				new_value = MIN(new_value, ssParams->max_real_var[i]);
				new_value = MAX(new_value, ssParams->min_real_var[i]);
				candidate->params[i] = new_value;
			}
			break;
		case '3':
			for (i = 0; i < ssParams->nreal; ++i)
			{
				new_value = base->params[i] + steps[i];
				new_value = MIN(new_value, ssParams->max_real_var[i]);
				new_value = MAX(new_value, ssParams->min_real_var[i]);
				candidate->params[i] = new_value;
			}
			break;
	}

}