#include "ss.h"

/*
	Evaluate the cost of an individual
 */
void evaluate_ind(SSType *ssParams, individual *ind, Input *inp, ScoreOutput *out){
	
	// sample_obj_func(ssParams, ind);
	ind->cost = objective_function(ind->params, ssParams, inp, out);
	// printf("%lf\n", objective_function(ind->params, ssParams, inp, out));
}

/*
	Evaluate the cost of every individual in a set
 */
void evaluate_set(SSType *ssParams, Set *set, int set_size, Input *inp, ScoreOutput *out){

	for (int i = 0; i < set_size; ++i)
	{
		evaluate_ind(ssParams, &(set->members[i]), inp, out);
	}

}

