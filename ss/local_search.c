#include "ss.h"
// #include "string.h"

/*
        Nelder-Mead
*/
void nelder_mead(SSType *ssParams, individual *ind){

}


/*
	Stochastic Hill Climbing
*/
void take_step(SSType *ssParams, double *params, double *new_params){

	double max_v, min_v;
	for (int i = 0; i < ssParams->nreal; ++i)
	{
		/* Determining the lower_bound and upper_bound for drawing a random number in neighborhood */
		min_v = MAX(ssParams->min_real_var[i], params[i] - ssParams->step_size);
		max_v = MIN(ssParams->max_real_var[i], params[i] + ssParams->step_size);
		new_params[i] = rndreal(min_v, max_v);
	}

}