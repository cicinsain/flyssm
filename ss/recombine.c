#include "ss.h"


/*
	Use the sub_sets_list members to produce candidates,
	the recombined_set will be generated.
 */
void generate_candiates(SSType *ssParams){

	int candidates_count = 0;
	double *dists = (double *)malloc( ssParams->nreal * sizeof(double));
	double mid_cost = ssParams->ref_set->members[ssParams->max_elite].cost;
	double diff;

	for (int i = 0; i < ssParams->subsets_list_size; ++i)
	{
		for (int j = 0; j < ssParams->nreal; ++j)
		{
			/* Step sizes in every direction */
			// TOOD: Check the formula
			diff = ( ssParams->subsets_list[i].members[0].params[j] -  ssParams->subsets_list[i].members[1].params[j] );
			// dists[j] = (diff * diff) / 2;
			dists[j] = diff / 2;
		}

		/* Generating new candidates */
		if ( ssParams->subsets_list[i].members[0].cost < mid_cost && ssParams->subsets_list[i].members[1].cost < mid_cost )		// type 1
		{
			// generate_ind_candidate(ssParams, &(ssParams->subsets_list[i].members[0]), &(ssParams->candidates_set->members[ candidates_count++ ]), dists, '0');

			// generate_ind_candidate(ssParams, &(ssParams->subsets_list[i].members[1]), &(ssParams->candidates_set->members[ candidates_count++ ]), dists, '0');

			generate_ind_candidate(ssParams, &(ssParams->subsets_list[i].members[0]), &(ssParams->candidates_set->members[ candidates_count++ ]), dists, '1');
			
			generate_ind_candidate(ssParams, &(ssParams->subsets_list[i].members[0]), &(ssParams->candidates_set->members[ candidates_count++ ]), dists, '2');

			generate_ind_candidate(ssParams, &(ssParams->subsets_list[i].members[1]), &(ssParams->candidates_set->members[ candidates_count++ ]), dists, '3');
			
			generate_ind_candidate(ssParams, &(ssParams->subsets_list[i].members[1]), &(ssParams->candidates_set->members[ candidates_count++ ]), dists, '3');

		}else
		if ( ssParams->subsets_list[i].members[0].cost < mid_cost && ssParams->subsets_list[i].members[1].cost >= mid_cost )	// type 2
		{
			// generate_ind_candidate(ssParams, &(ssParams->subsets_list[i].members[0]), &(ssParams->candidates_set->members[ candidates_count++ ]), dists, '0');
			
			generate_ind_candidate(ssParams, &(ssParams->subsets_list[i].members[0]), &(ssParams->candidates_set->members[ candidates_count++ ]), dists, '1');

			generate_ind_candidate(ssParams, &(ssParams->subsets_list[i].members[0]), &(ssParams->candidates_set->members[ candidates_count++ ]), dists, '2');

			generate_ind_candidate(ssParams, &(ssParams->subsets_list[i].members[1]), &(ssParams->candidates_set->members[ candidates_count++ ]), dists, '3');

		}else
		if ( ssParams->subsets_list[i].members[0].cost >= mid_cost && ssParams->subsets_list[i].members[1].cost >= mid_cost )	// type 3
		{
			
			if (rndreal(0, 1) < 0.5)
			{
				generate_ind_candidate(ssParams, &(ssParams->subsets_list[i].members[0]), &(ssParams->candidates_set->members[ candidates_count++ ]), dists, '1');
			}else
			{
				generate_ind_candidate(ssParams, &(ssParams->subsets_list[i].members[1]), &(ssParams->candidates_set->members[ candidates_count++ ]), dists, '3');
			}

			generate_ind_candidate(ssParams, &(ssParams->subsets_list[i].members[0]), &(ssParams->candidates_set->members[ candidates_count++ ]), dists, '2');

		}


	}

	ssParams->candidates_set_size = candidates_count;
	// printf("canSet size: %d\n", candidates_count);

	free(dists);

}


void generate_ind_candidate(SSType *ssParams, individual *base, individual *candidate, double *dists, char type){
	// FIXME: It happens that the variables go outside the boundaries, I could be solved by finding the MIN and MAX between new values and min and max of variables.

	int i;
	double new_value;
	double rnd;
	switch (type){
		case '0':
			rnd = rndreal(0, 1);
			for (i = 0; i < ssParams->nreal; ++i)
			{
				new_value = base->params[i] - (dists[i] * rnd);
				// new_value = MIN(new_value, ssParams->max_real_var[i]);
				// new_value = MAX(new_value, ssParams->min_real_var[i]);
				candidate->params[i] = new_value;
			}
			break;
		case '1':
			for (i = 0; i < ssParams->nreal; ++i)
			{
				new_value = base->params[i] - (dists[i] * rndreal(0,1));
				// new_value = MIN(new_value, ssParams->max_real_var[i]);
				// new_value = MAX(new_value, ssParams->min_real_var[i]);
				candidate->params[i] = new_value;
			}
			break;
		case '2':
			for (i = 0; i < ssParams->nreal; ++i)
			{
				new_value = base->params[i] + (dists[i] * rndreal(0,1));
				// new_value = MIN(new_value, ssParams->max_real_var[i]);
				// new_value = MAX(new_value, ssParams->min_real_var[i]);
				candidate->params[i] = new_value;
			}
			break;
		case '3':
			for (i = 0; i < ssParams->nreal; ++i)
			{
				new_value = base->params[i] + (dists[i] * rndreal(0,1));
				// new_value = MIN(new_value, ssParams->max_real_var[i]);
				// new_value = MAX(new_value, ssParams->min_real_var[i]);
				candidate->params[i] = new_value;
			}
			break;
	}

}
