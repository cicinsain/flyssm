/* Data initializtion routines */

# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <errno.h>
# include <math.h>

# include "nsga2.h"
# include "rand.h"

/* Function to initialize a population randomly */
void initialize_pop (NSGA2Type *nsga2Params, population *pop)
{
    int i;
    for (i=0; i<nsga2Params->popsize; i++)
    {
        initialize_ind (nsga2Params, &(pop->ind[i]));
    }
    return;
}

/* Function to initialize an individual randomly */
void initialize_ind (NSGA2Type *nsga2Params, individual *ind)
{
    int j, k;
    if (nsga2Params->nreal!=0)
    {
        for (j=0; j<nsga2Params->nreal; j++)
        {
            // The min_realvar and max_realvar are computed using scramble.c so the randomly generated number will fit into the boundaries
            ind->xreal[j] = rndreal (nsga2Params->min_realvar[j], nsga2Params->max_realvar[j]);
            // printf("%f\n", ind->xreal[j]);
        }
        // printf("\n");
    }
    if (nsga2Params->nbin!=0)
    {
        for (j=0; j<nsga2Params->nbin; j++)
        {
            for (k=0; k<nsga2Params->nbits[j]; k++)
            {
                if (randomperc() <= 0.5)
                {
                    ind->gene[j][k] = 0;
                }
                else
                {
                    ind->gene[j][k] = 1;
                }
            }
        }
    }
    return;
}

void initialize_pop_from_file(NSGA2Type *nsga2Params, population *pop, char *inname){
    int i = 0;
    // char init_pop_fname[256]="/home/amir/Dropbox/Study/master/fly_nsga2/tllg58c13-v5.55_scrambled_init_pop";
    char init_pop_fname[256]="/home/amir/Dropbox/Study/master/fly_nsga2/output";
    // sprintf(init_pop_fname, "%s_%s", inname, "scrambled_init_pop");
    FILE* init_pop_file = fopen(init_pop_fname, "r");


    char line[1024];
    // printf("%d\n", nsga2Params->popsize);
    while(fgets(line, 1024, init_pop_file) && (i < nsga2Params->popsize)){
        // printf("%d", i);
        char* tmp = strdup(line);
        initialize_ind_from_file (nsga2Params, &(pop->ind[i++]), tmp);
        free(tmp);
    }

    // for (i=0; i<nsga2Params->popsize; i++){
    //     initialize_ind (nsga2Params, &(pop->ind[i]));
    // }
}

void initialize_ind_from_file(NSGA2Type *nsga2Params, individual *ind, char* line)
{
    int j = 0;
    const char* tok;
    for (tok = strtok(line, ",");tok && *tok && (j < nsga2Params->nreal);tok = strtok(NULL, ",\n")){
        ind->xreal[j++] = atof(tok);
        printf("%lf\n", ind->xreal[j-1]);
    }
    printf("\n");
    // printf("-%d,", j);
}

