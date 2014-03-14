/* Routines for storing population data into files */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "nsga2.h"
# include "rand.h"
# include "fly_io.h"

/* Function to print the information of a population in a file */
void report_pop (NSGA2Type *nsga2Params,  population *pop, FILE *fpt)
{
    int i, j, k;
    for (i=0; i<nsga2Params->popsize; i++)
    {
        for (j=0; j<nsga2Params->nobj; j++)
        {
            fprintf(fpt,"%e\t",pop->ind[i].obj[j]);
        }
        if (nsga2Params->ncon!=0)
        {
            for (j=0; j<nsga2Params->ncon; j++)
            {
                fprintf(fpt,"%e\t",pop->ind[i].constr[j]);
            }
        }
        if (nsga2Params->nreal!=0)
        {
            for (j=0; j<nsga2Params->nreal; j++)
            {
                fprintf(fpt,"%13.9f\t",pop->ind[i].xreal[j]);
            }
        }
        if (nsga2Params->nbin!=0)
        {
            for (j=0; j<nsga2Params->nbin; j++)
            {
                for (k=0; k<nsga2Params->nbits[j]; k++)
                {
                    fprintf(fpt,"%d\t",pop->ind[i].gene[j][k]);
                }
            }
        }
        fprintf(fpt,"%e\t",pop->ind[i].constr_violation);
        fprintf(fpt,"%d\t",pop->ind[i].rank);
        fprintf(fpt,"%e\n",pop->ind[i].crowd_dist);
    }
    return;
}

/* Function to print the information of feasible and non-dominated population in a file */
void report_feasible (NSGA2Type *nsga2Params,  population *pop, FILE *fpt)
{
    int i, j, k;
    for (i=0; i<nsga2Params->popsize; i++)
    {
        if (pop->ind[i].constr_violation == 0.0 && pop->ind[i].rank==1)
        {
            for (j=0; j<nsga2Params->nobj; j++)
            {
                fprintf(fpt,"%e\t",pop->ind[i].obj[j]);
            }
            if (nsga2Params->ncon!=0)
            {
                for (j=0; j<nsga2Params->ncon; j++)
                {
                    fprintf(fpt,"%e\t",pop->ind[i].constr[j]);
                }
            }
            if (nsga2Params->nreal!=0)
            {
                for (j=0; j<nsga2Params->nreal; j++)
                {
                    fprintf(fpt,"%13.9f\t",pop->ind[i].xreal[j]);
                }
            }
            if (nsga2Params->nbin!=0)
            {
                for (j=0; j<nsga2Params->nbin; j++)
                {
                    for (k=0; k<nsga2Params->nbits[j]; k++)
                    {
                        fprintf(fpt,"%d\t",pop->ind[i].gene[j][k]);
                    }
                }
            }
            fprintf(fpt,"%e\t",pop->ind[i].constr_violation);
            fprintf(fpt,"%d\t",pop->ind[i].rank);
            fprintf(fpt,"%e\n",pop->ind[i].crowd_dist);
        }
    }
    return;
}


void write_params_to_fly_output_standard(NSGA2Type *nsga2Params, Input *inp, population *pop, char *inname){
    int i, j;
    ParamList *ptab;

    // printf("**********************\n%s\n", inname);

    printf("\nCreating output files for each parameters set.\n");

    ptab = inp->tra.array;      // Get the pointer to the 'inp' array.
    for (i = 0; i < nsga2Params->popsize; ++i){
        /* code */
        for (j = 0; j < inp->tra.size; ++j){
            /* Modifying the local parameters of 'inp' struct */
            *( ptab[j].param  ) = pop->ind[i].xreal[j];
        }
        // Each run will create an output file in form of 'inname_parm_XXXXXX' in which 'XXXXXXX'
        // will replace by random string.
        WriteParameters(inname, &(inp->zyg.parm), "eqparms", 9, inp->zyg.defs);
    }
    printf("Done.\n");
}