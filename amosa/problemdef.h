/* Test problem definitions */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "amosa.h"
// # include "rand.h"

# include "score.h"
# include "zygotic.h"

/*# define ctp1 */


// ParamList *ptab;                /* array of pointers to parameters and ranges */



# define GeneBased
/*# define GeneTimeBased*/
/*# define GeneNucleiBased*/

# define MAX_DIFF 3771450

// ScoreOutput out;
/*
 Gene based problem
 # of real variable = depends on the number of genes (ngenes^2 + 4*ngenes)
 # of bin variable  = 0
 # of objectives    = inp->zyg.defs.ngenes
 # of constraints   = 0 (Actually there are 3 box constraints that are handling with the
 NSGA2 algorithm.)
 */
#ifdef GeneBased


void objective_function(double *s, AMOSAType *amosaParams, Input *inp, ScoreOutput *out){
    int i, g;
    
    // out.score          = 1e38;           // start with a very large number
    // out.penalty        = 0;
    // out.size_resid_arr = 0;
    // out.jacobian       = NULL;
    // out.residuals      = NULL;
    
    
    // ptab = inp->tra.array;
    // Updating the 'inp->zyg.parm' variables with parameter set of each individual.
    for (i = 0; i < inp->tra.size; ++i){
        // *( ptab[i].param  ) = xreal[i];
        *( inp->tra.array[i].param  ) = s[i];
    }
    
    // inp->lparm = CopyParm( inp->zyg.parm, &( inp->zyg.defs ) );
    // Computing the out using the inp configuration
    Score( inp, out, 0 );
    
    // FIXME: Many of the run end up in this situation =>  Fixed by shrinking the range of T, E, and M from (-1, 1) to (-0.1, 0.1)
    // TODO: I can track the FORBIDDEN_MOVE and do something after getting several FORBIDDEN_MOVE in a row like increasing the mutation rate or something else.
    // if (out.score != FORBIDDEN_MOVE){
    //     printf("Got it!\n");
    // }
    
    // Replacing the objectives value with scores from Score()
    for (g = 0; g < inp->zyg.defs.ngenes; ++g){
        // TODO: Shall I perform some scaling here
        amosaParams->d_eval[g] = out->nsga2_outs[g]->score;
        // printf("%lf,", out->nsga2_outs[g]->score);
    }
    // printf("\n");
    // free(&out);
    // free(ptab);
    
}

#endif


