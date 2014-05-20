#include "ss.h"
// #include "test_problems/Ackleys.h"
// #include "test_problems/AluffiPentini.h"		// Accurate!
// #include "test_problems/BeckerLago.h"		// Accurate!
// #include "test_problems/Bohachevsky1.h"		// Accurate!
// #include "test_problems/Bohachevsky2.h"		// Accurate!
// #include "test_problems/Branin.h"			// Accurate!
// #include "test_problems/Camel3.h"			// Accurate!
// #include "test_problems/Camel6.h"			// Accurate!
// #include "test_problems/CosMix2.h"			// Accurate!
// #include "test_problems/CosMix4.h"			// Accurate!
// #include "test_problems/DekkersAarts.h"		// Accurate!
// #include "test_problems/Easom.h"				// Accurate!
// #include "test_problems/EMichalewicz.h"		// Accurate!
// #include "test_problems/Expo.h"				// Accurate!
// #include "test_problems/GoldPrice.h"			// Accurate!
// #include "test_problems/Hartman3.h"			// Accurate!
// #include "test_problems/Hartman6.h"			// Accurate!
// #include "test_problems/Kowalik.h"			// Accurate!
// #include "test_problems/LM1.h"				// Accurate!
// #include "test_problems/LM2n10.h"			// Accurate!
// #include "test_problems/LM2n5.h"				// Accurate!
// #include "test_problems/MeyerRoth.h"			// Accurate!
// #include "test_problems/MieleCantrell.h"		// Accurate!
// #include "test_problems/ModRosenbrock.h"		// Accurate!
// #include "test_problems/MultiGauss.h"		// Accurate!
// #include "test_problems/Neumaier2.h"			// Accurate!
// #include "test_problems/Neumaier3.h"			// Accurate!
// #include "test_problems/Paviani.h"			// Accurate!
// #include "test_problems/Periodic.h"			// Accurate!
// #include "test_problems/PowellQ.h"			// Accurate!
// #include "test_problems/Schaffer1.h"			// Accurate!
// #include "test_problems/Schaffer2.h"			// Accurate!
// #include "test_problems/Schubert.h"			// Accurate!
// #include "test_problems/Shekel10.h"			// Accurate!
// #include "test_problems/Shekel5.h"			// Accurate!
// #include "test_problems/Shekel7.h"			// Accurate!
// #include "test_problems/Shekelfox10.h"		// Accurate!
// #include "test_problems/Shekelfox5.h"		// Accurate!
// #include "test_problems/Wood.h"				// Accurate!
// #include "test_problems/Zeldasine10.h"		// Accurate!
// #include "test_problems/Zeldasine20.h"		// Accurate!
// #include "test_problems/Gulf.h"				// Accurate!
// #include "test_problems/Rosenbrock.h"		// Accurate!
// #include "test_problems/Beale.h"				// Accurate!

// #include "test_problems/Helical.h"			// Not accurate!	// not many subsituation
// #include "test_problems/Modlangerman.h"		// Not accurate!
// #include "test_problems/Oddsquare.h"			// Not accurate!
// #include "test_problems/PriceTransistor.h"	// Not accurate!
// #include "test_problems/Rastrigin.h"			// Not accurate!
// #include "test_problems/Salomon.h"			// Not accurate!	// not many subsituation
// #include "test_problems/Schwefel.h"			// Not accurate!
// #include "test_problems/STChebychev9.h"		// Not accurate!
// #include "test_problems/STChebychev17.h"		// Not accurate!
// #include "test_problems/Griewank.h"			// Not accurate!

// #include "test_problems/Hosaki.h"			// totally off!
// #include "test_problems/McCormic.h"			// totally off!



// void init_sample_params(SSType *ssParams){

// 	ssParams->ref_set_size = 20;		// TODO: Check the (ref_set_size < diverse_set_size); otherwise the segmentation will happen
// 	ssParams->max_iter = 500;
// 	ssParams->step_size = 0.0005;	// TODO: Should be computed using the boundaries of variables
// 	ssParams->max_no_improve = 100;	
// 	// ssParams->diverse_set_size = 200;
// 	ssParams->diverse_set_size = 10 * ssParams->ref_set_size;
// 	ssParams->max_elite = ssParams->ref_set_size / 2;
// 	ssParams->pair_size = 2;
// 	ssParams->subsets_list_size = 0;
// 	// ssParams->p = ssParams->ref_set_size / 2;		// (ref_set_size / 2) or (ref_set_size)
// 	ssParams->p = 8;		// (ref_set_size / 2) or (ref_set_size)
// 							// NOTE: p should not be large otherwise it breaks the balance of the frequency matrix!
// 	ssParams->local_search_criteria = 5;

// 	ssParams->dist_epsilon = 0.0001;
// 	ssParams->fitness_epsilon = .00001;

// 	ssParams->sol = (double)SOL;

// #ifdef TEST_PROBLEM
// 	ssParams->nreal = N;
// 	printf("\n--------------------------------");
// 	printf("\n Analytic Solution: %lf \n ", (double)SOL);
// 	printf("---------------------------------\n\n");
// #else
// 	ssParams->nreal = 2;
// #endif

// 	ssParams->min_real_var = (double *)malloc(ssParams->nreal * sizeof(double));
// 	ssParams->max_real_var = (double *)malloc(ssParams->nreal * sizeof(double));


// #ifdef TEST_PROBLEM
// 	bounds(ssParams->min_real_var, ssParams->max_real_var);
// #else
// 	for (int i = 0; i < ssParams->nreal; ++i)
// 	{
// 		ssParams->min_real_var[i] = -5;	// TODO: It should be -100
// 		ssParams->max_real_var[i] = 5;
// 	}
// #endif

// 	// Determining the step_size

// 	// What? This is what implemented in pyCleverAlgorithm but I don't like it.
// 	// ssParams->step_size = (ssParams->max_real_var[0] - ssParams->min_real_var[0]) * 0.005;
// 	// The step should be small
// 	ssParams->step_size = 0.0005;


// }

// void sample_obj_func(SSType *ssParams, individual *ind){
// 	double cost = 0;
// 	double x, y;
// 	x = ind->params[0];
// 	y = ind->params[1];

// #if !defined(TEST_PROBLEM)
// 	// nreal = 3
// 	// ---------------
// 	for (int i = 0; i < ssParams->nreal; ++i)
// 	{
// 		cost += ind->params[i] * ind->params[i];
// 	}; 	ind->cost = cost;
	
// 	// nreal = 2
// 	// -------------------
// 	// ind->cost =  -20 * exp(-0.2 * sqrt(.5 * (x*x + y*y))) - exp(0.5 * (cos(2 * pi * x) + cos(2 * pi * y))) + 20 + eul;
// 	// ind->cost = -cos(x) * cos(y) * exp(-((x-pi)*(x-pi) + (y-pi)*(y-pi)));
// 	// ind->cost = -0.0001 * pow(fabs( sin(x) * sin(y) * exp( fabs( 100 - sqrt( x*x + y*y )/pi ) ) + 1), 0.1);
// 	// ind->cost = -(y + 4) * sin( sqrt( fabs( y + x/2 + 47 ))) - x * sin( sqrt( fabs(x - (y + 47))));
// 	// ind->cost = 0.5 + ( pow( sin( x*x - y*y ), 2) - 0.5 ) / pow( ( 1 + 0.001*( x*x + y*y) ), 2);
// 	// ind->cost = 10 * 2 + (x − 10 * cos(2* pi * x)) + (y − 10 * cos(2 * pi * y));

// #else
// 	ind->cost = objfn(ind->params);
// #endif

// }



double objective_function(double *s, SSType *ssParams, Input *inp, ScoreOutput *out){


    for (int i = 0; i < inp->tra.size; ++i){
        *( inp->tra.array[i].param  ) = s[i];
        // printf("%lf, ", s[i]);
    }
    // printf("\n\n");
    
    // inp->lparm = CopyParm( inp->zyg.parm, &( inp->zyg.defs ) );
    // Computing the out using the inp configuration
    Score( inp, out, 0 );
    ssParams->n_function_evals++;
    
    // FIXME: Many of the run end up in this situation =>  Fixed by shrinking the range of T, E, and M from (-1, 1) to (-0.1, 0.1)
    // TODO: I can track the FORBIDDEN_MOVE and do something after getting several FORBIDDEN_MOVE in a row like increasing the mutation rate or something else.
    // if (out.score != FORBIDDEN_MOVE){
    //     printf("Got it!\n");
    // }
    // printf("%lf\n", out->score);
    return sqrt(out->score/1970);
    // return out->score;

    // Replacing the objectives value with scores from Score()
    // for (g = 0; g < inp->zyg.defs.ngenes; ++g){
    //     // TODO: Shall I perform some scaling here
    //     amosaParams->d_eval[g] = out->nsga2_outs[g]->score;
    //     // printf("%lf,", out->nsga2_outs[g]->score);
    // }
}

