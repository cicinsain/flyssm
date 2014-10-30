#include "ss.h"
#include <string.h>
#include <gsl/gsl_vector_double.h>

SSType *ssParamsLocal;
Input *inpLocal;
ScoreOutput *outLocal;
individual indLocal;

void refine_subsets_list(SSType *ssParams, char method, Input *inp, ScoreOutput *out){

	for (int i = 0; i < ssParams->subsets_list_size; ++i)
	{
		refine_set(ssParams, &(ssParams->subsets_list[i]), ssParams->pair_size, method, inp, out);
	}
}

void refine_set(SSType *ssParams, Set *set, int set_size, char method, Input *inp, ScoreOutput *out){

	int closest_member_index;
        
        //DAMJAN TEST - REMOVE
        //set_size = 1;
        //REMOVE
        
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
								 	ssParams->n_refinement++;
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
        
        outLocal = out;
        inpLocal = inp;
        ssParamsLocal = ssParams;
    
        individual new_candidate;

	double *new_params = (double *)malloc( ssParams->nreal * sizeof(double));

        //Damjan
        if (method == 'n') {
            //printf("running Nelder Mead cost before = %lg\n", ind->cost);
            int debug = 0;
            const gsl_multimin_fminimizer_type *T;
            gsl_multimin_fminimizer *s;
            int status;
            unsigned int iter = 0;
            const size_t p = ssParams->nreal;
            gsl_multimin_function f;
            // print_Ind(eSSParams, ind);
            // print_Ind(eSSParams, eSSParams->best);

            gsl_vector_view x = gsl_vector_view_array (ind->params, p);
            // gsl_vector_fprintf(stdout, &x.vector, "%f");
            // printf("---\n");
            const gsl_rng_type * type;
            gsl_rng * r;

            gsl_vector *ss = gsl_vector_alloc(p);
            gsl_vector_set_all(ss, 0.1);

            gsl_rng_env_setup();

            type = gsl_rng_default;
            r = gsl_rng_alloc (type);

            f.f = &nelder_objfn;                                    
            f.n = p;		
            f.params = &x.vector;		

            double size;

            T = gsl_multimin_fminimizer_nmsimplex2;
            s = gsl_multimin_fminimizer_alloc (T, p);
            gsl_multimin_fminimizer_set(s, &f, &x.vector, ss);

            do
            {
                iter++;
                status = gsl_multimin_fminimizer_iterate(s);

                if (debug == 1)
                {
                    for (int i = 0; i < s->x->size; ++i){
                        printf("%lf, ", gsl_vector_get(s->x, i));
                    }
                    printf("--> %lf\n ", s->fval);
                }

                if (status)
                    break;

                size = gsl_multimin_fminimizer_size(s);
                status = gsl_multimin_test_size( size, 1e-3);
                //printf("iter=%d\n", iter);
            }            
            while (status == GSL_CONTINUE && iter < ssParams->max_no_improve);
            //while (status == GSL_CONTINUE && iter < 100);
                                    
            if (iter != 0){
            // printf("%d\n", iter);
                ind->params = s->x->data;
                ind->cost = s->fval;
                //eSSParams->stats->n_successful_localSearch++;
                //eSSParams->stats->n_local_search_iterations += iter;
            }
            //printf("Nelder finished - free mem\n");
            //gsl_multimin_fminimizer_free(s); //crashing here
            gsl_vector_free(ss);            
            gsl_rng_free (r);
            return;
        }
        //Damjan
        
	for (int i = 0; i < ssParams->max_no_improve; ++i)
	{
            //printf("running localsearch max interations %d\n", ssParams->max_no_improve);
		/* Run the local optimization procedure */
		switch (method){
			case 't':		// Stochastic Hill Climbing 
				{
					//ssParams->n_refinement++;
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
	free(new_params);

	// deallocate_ind_memory(ssParams, new_candidate);
}


double nelder_objfn(const gsl_vector *x, void *data){
    indLocal.params = x->data;
    //memcpy(new_candidate.params[i], x->data, ssParamsLocal->nreal*sizeof(double));            
    int i;
    
    for (i = 0; i < x->size; i++) {
        //printf ("%lg ", indLocal.params[i]);
    }
    //printf("\n");
    //printf ("%lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg\n",indLocal.params[0],indLocal.params[1],indLocal.params[2],indLocal.params[3],indLocal.params[4],indLocal.params[5],indLocal.params[6],indLocal.params[7]);
    evaluate_ind(ssParamsLocal, &(indLocal), inpLocal, outLocal);
    //printf("cost = %lg\n", indLocal.cost);
    return indLocal.cost;
}






