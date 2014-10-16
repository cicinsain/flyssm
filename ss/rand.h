#include "ss.h"

double rndreal (double low, double high) 
{
    //printf("%lg, %lg, %lg %lg %lg %lg\n", low, high, (high-low), (double)rand(), ((double)rand() / RAND_MAX), (low + (high-low) * ((double)rand() / RAND_MAX) ));
    return (low + (high-low) * ((double)rand() / RAND_MAX) );
}

/*
	Assign a random real vector to every member of the set
 */
void generate_random_set(SSType *ssParams, Set *set, int set_size){
	for (int i = 0; i < set_size; ++i)
	{
		random_ind(ssParams, &(set->members[i]));
	}
}

/*
	Generate random values with respect to the min_real_var and max_real_var for every individual in a set
 */
//TODO: Could be generalize by adding two more variables `min`, `max`
void random_ind(SSType *ssParams, individual *ind){

	for (int i = 0; i < ssParams->nreal; ++i)
	{
		ind->params[i] = rndreal(ssParams->min_real_var[i], ssParams->max_real_var[i]);
	}
}