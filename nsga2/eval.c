/* Routine for evaluating population members  */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "nsga2.h"
# include "rand.h"

/* Routine to evaluate objective function values and constraints for a population */
void evaluate_pop (Input *inp, ScoreOutput *out, NSGA2Type *nsga2Params, population *pop)
{
    int i;
    /*======
        Amir: TODO: Possible parallelization can be done here but I need to take care of
                    'inp' that carry variables, I may need to define several out for each
                    node. That shouldn't be hard!*/
    for (i=0; i<nsga2Params->popsize; i++)
    {
        evaluate_ind (inp, out, nsga2Params, &(pop->ind[i]));
    }
    return;
}

/* Routine to evaluate objective function values and constraints for an individual */
void evaluate_ind (Input *inp, ScoreOutput *out, NSGA2Type *nsga2Params, individual *ind)
{
    int j;
    objective_function (inp, out, ind->xreal, ind->xbin, ind->gene, ind->obj, ind->constr);

    if (nsga2Params->ncon==0)
    {
        ind->constr_violation = 0.0;
    }
    else
    {
        ind->constr_violation = 0.0;
        for (j=0; j<nsga2Params->ncon; j++)
        {
            if (ind->constr[j]<0.0)
            {
                ind->constr_violation += ind->constr[j];
            }
        }
    }
    return;
}
