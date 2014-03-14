/* Test problem definitions */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "nsga2.h"
# include "rand.h"

# include "score.h"
# include "zygotic.h"

/*# define ctp1 */


// ParamList *ptab;                /* array of pointers to parameters and ranges */



# define GeneBased
/*# define GeneTimeBased*/
/*# define GeneNucleiBased*/

# define MAX_DIFF 3771450

ScoreOutput out;
/*
    Gene based problem
    # of real variable = depends on the number of genes (ngenes^2 + 4*ngenes)
    # of bin variable  = 0
    # of objectives    = inp->zyg.defs.ngenes
    # of constraints   = 0 (Actually there are 3 box constraints that are handling with the 
                                NSGA2 algorithm.)
*/
#ifdef GeneBased

void objective_function(Input *inp, ScoreOutput *out, double *xreal, double *xbin, int **gene, double *obj, double *constr){
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
        *( inp->tra.array[i].param  ) = xreal[i];
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
        obj[g] = out->nsga2_outs[g]->score;
    }
    // free(&out);
    // free(ptab);

}

#endif






/*  Test problem CTP1
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 2
    */
#ifdef ctp1
void objective_function (Input *inp, double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    double g;
    g = 1.0 + xreal[1];
    obj[0] = xreal[0];
    obj[1] = g*exp(-obj[0]/g);
    constr[0] = obj[1]/(0.858*exp(-0.541*obj[0]))-1.0;
    constr[1] = obj[1]/(0.728*exp(-0.295*obj[0]))-1.0;
    return;
}
#endif





    // printf("GeneBased Objective Function\n");
/*================================= [Aug 29, 2013]
    void objective_function(Input *inp, double *xreal, double *xbin, int **gene, double *obj, double *constr)
        // Amir: With this definition, I should:
                    1. I can just use the PArrPtr which is the list of parameters to twaek!
                        1.1. inp->tra.array.param (PArrPtr)
                        1.2. inp->tra.array.param_range (Range)
                    2. Send the inp to the Score(&inp, out, jacobian ) to compute the score
                        2.1. Modify the score function to return an array of ScoreOutput instead of one!
                    3. Take the output, set the obj[]
                    4. free up the memory.
                    - By this method, I don't need to change the 'static Input inp', just pass it from one to another.
*/
    /*====================== [Aug 27, 2013]
        Amir:
            1. Call the Score function with xreal variable
            2. Getting the answers in a form of output struct
                1. Report the exact values for the objective from Eval functions in zygnote.c (easier)
                2. Or send the answer here and compute the objectives RMS here. (not recommended).
            3. Set the value to the obj pointer and DONE!
    */ 
    // Amir: I should call the MoveSA() function, and return the out->score for each genes as *obj
