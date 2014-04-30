/* NSGA-II routine (implementation of the 'main' function) */

# include <string.h>
# include <time.h>

# include "nsga2.h"
# include "rand.h"
// # include "fly_io.h"
// # include "maternal.h"

// # include "global.h"


FILE *fpt1;
FILE *fpt2;
FILE *fpt3;
FILE *fpt4;
FILE *fpt5;
FILE *gp;
population *parent_pop;
population *child_pop;
population *mixed_pop;
ScoreOutput out;



void InitNSGA2(Input *inp, NSGA2Type *nsga2Params, char *inname)
{
    int i, g;
    char *b_name;

    out.score          = 1e38;           // start with a very large number
    out.penalty        = 0;
    out.size_resid_arr = 0;
    out.jacobian       = NULL;
    out.residuals      = NULL;
    out.nsga2_outs     = (ScoreOutput **) malloc( inp->zyg.defs.ngenes * sizeof(ScoreOutput));
    for (g = 0; g < inp->zyg.defs.ngenes; ++g){
        out.nsga2_outs[g] = (ScoreOutput *) malloc ( sizeof(ScoreOutput) );
    }

    printf("Initialize NSGA2\n");
    
    char *fpt1_name;    fpt1_name = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    char *fpt2_name;    fpt2_name = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    char *fpt3_name;    fpt3_name = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    char *fpt4_name;    fpt4_name = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    char *fpt5_name;    fpt5_name = ( char * ) calloc( MAX_RECORD, sizeof( char ) );

    b_name = strrchr(inname, '/');
    sprintf(fpt1_name, "initial_pop_%s.out", &b_name[1]);
    sprintf(fpt2_name, "final_pop_%s.out", &b_name[1]);
    sprintf(fpt3_name, "best_pop_%s.out", &b_name[1]);
    sprintf(fpt4_name, "all_pop_%s.out", &b_name[1]);
    sprintf(fpt5_name, "params_%s.out", &b_name[1]);

    fpt1 = fopen(fpt1_name, "w"); free(fpt1_name);
    fpt2 = fopen(fpt2_name, "w"); free(fpt2_name);
    fpt3 = fopen(fpt3_name, "w"); free(fpt3_name);
    fpt4 = fopen(fpt4_name, "w"); free(fpt4_name);
    fpt5 = fopen(fpt5_name, "w"); free(fpt5_name);

    fprintf(fpt1,"# This file contains the data of initial population\n");
    fprintf(fpt2,"# This file contains the data of final population\n");
    fprintf(fpt3,"# This file contains the data of final feasible population (if found)\n");
    fprintf(fpt4,"# This file contains the data of all generations\n");
    fprintf(fpt5,"# This file contains information about inputs as read by the program\n");

    printf("\nInput data successfully entered, now performing initialization \n");
    fprintf(fpt5,"\n Population size = %d",nsga2Params->popsize);
    fprintf(fpt5,"\n Number of generations = %d",nsga2Params->ngen);
    fprintf(fpt5,"\n Number of objective functions = %d",nsga2Params->nobj);
    fprintf(fpt5,"\n Number of constraints = %d",nsga2Params->ncon);
    fprintf(fpt5,"\n Number of real variables = %d",nsga2Params->nreal);
    if (nsga2Params->nreal!=0)
    {
        for (i=0; i<nsga2Params->nreal; i++)
        {
            fprintf(fpt5,"\n Lower limit of real variable %d = %e",i+1,nsga2Params->min_realvar[i]);
            fprintf(fpt5,"\n Upper limit of real variable %d = %e",i+1,nsga2Params->max_realvar[i]);
        }
        fprintf(fpt5,"\n Probability of crossover of real variable = %e",nsga2Params->pcross_real);
        fprintf(fpt5,"\n Probability of mutation of real variable = %e",nsga2Params->pmut_real);
        fprintf(fpt5,"\n Distribution index for crossover = %e",nsga2Params->eta_c);
        fprintf(fpt5,"\n Distribution index for mutation = %e",nsga2Params->eta_m);
    }
    fprintf(fpt5,"\n Number of binary variables = %d",nsga2Params->nbin);
    if (nsga2Params->nbin!=0)
    {
        for (i=0; i<nsga2Params->nbin; i++)
        {
            fprintf(fpt5,"\n Number of bits for binary variable %d = %d",i+1,nsga2Params->nbits[i]);
            fprintf(fpt5,"\n Lower limit of binary variable %d = %e",i+1,nsga2Params->min_binvar[i]);
            fprintf(fpt5,"\n Upper limit of binary variable %d = %e",i+1,nsga2Params->max_binvar[i]);
        }
        fprintf(fpt5,"\n Probability of crossover of binary variable = %e",nsga2Params->pcross_bin);
        fprintf(fpt5,"\n Probability of mutation of binary variable = %e",nsga2Params->pmut_bin);
    }
    fprintf(fpt5,"\n Seed for random number generator = %e",nsga2Params->seed);
    fprintf(fpt5,"\n Number of data points = %d", inp->zyg.ndp);
    nsga2Params->bitlength = 0;
    if (nsga2Params->nbin!=0)
    {
        for (i=0; i<nsga2Params->nbin; i++)
        {
            nsga2Params->bitlength += nsga2Params->nbits[i];
        }
    }
    fprintf(fpt1,"# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n",nsga2Params->nobj,nsga2Params->ncon,nsga2Params->nreal,nsga2Params->bitlength);
    fprintf(fpt2,"# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n",nsga2Params->nobj,nsga2Params->ncon,nsga2Params->nreal,nsga2Params->bitlength);
    fprintf(fpt3,"# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n",nsga2Params->nobj,nsga2Params->ncon,nsga2Params->nreal,nsga2Params->bitlength);
    fprintf(fpt4,"# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n",nsga2Params->nobj,nsga2Params->ncon,nsga2Params->nreal,nsga2Params->bitlength);

    parent_pop = (population *)malloc(sizeof(population));
    child_pop = (population *)malloc(sizeof(population));
    mixed_pop = (population *)malloc(sizeof(population));
    allocate_memory_pop (nsga2Params, parent_pop, nsga2Params->popsize);
    allocate_memory_pop (nsga2Params, child_pop, nsga2Params->popsize);
    allocate_memory_pop (nsga2Params, mixed_pop, 2*nsga2Params->popsize);

    randomize(nsga2Params->seed);
    initialize_pop (nsga2Params,  parent_pop);
    // initialize_pop_from_file (nsga2Params,  parent_pop, inname);
    printf("\nInitialization done, now performing first generation\n");
    decode_pop(nsga2Params, parent_pop);
    evaluate_pop (inp, &out, nsga2Params, parent_pop);
    assign_rank_and_crowding_distance (nsga2Params, parent_pop);
    report_pop (nsga2Params, parent_pop, fpt1);
    fprintf(fpt4,"# gen = 1\n");
    report_pop(nsga2Params, parent_pop,fpt4);
    printf("\n gen = 1");
    fflush(stdout);
    // if (nsga2Params->choice!=0)    onthefly_display (nsga2Params, parent_pop,gp,1);
    fflush(fpt1);
    fflush(fpt2);
    fflush(fpt3);
    fflush(fpt4);
    fflush(fpt5);
    sleep(1);
}




int RunNSGA2(Input *inp, NSGA2Type *nsga2Params, char *inname)
{
    int i;


    /*=======================
        Amir: Running the algorithm
    	*/
    // printf("Start the generations...\n");

    for (i=2; i<=nsga2Params->ngen; i++)
    {
        selection (nsga2Params,  parent_pop, child_pop);
        mutation_pop (nsga2Params,  child_pop);
        decode_pop(nsga2Params,  child_pop);
        evaluate_pop(inp, &out, nsga2Params,  child_pop);
        merge (nsga2Params,  parent_pop, child_pop, mixed_pop);
        fill_nondominated_sort (nsga2Params,  mixed_pop, parent_pop);
        /* Comment following four lines if information for all
        generations is not desired, it will speed up the execution */
        fprintf(fpt4,"# gen = %d\n",i);
        report_pop(nsga2Params, parent_pop,fpt4);
        fflush(fpt4);
        // if (nsga2Params->choice!=0)    onthefly_display (nsga2Params, parent_pop,gp,i);
        printf("\n gen = %d",i);
    }
    printf("\n\nGenerations finished, now reporting solutions");

    report_pop(nsga2Params,  parent_pop,fpt2);
    report_feasible(nsga2Params,  parent_pop,fpt3);
    write_params_to_fly_output_standard(nsga2Params, inp, parent_pop, inname);


    if (nsga2Params->nreal!=0)
    {
        fprintf(fpt5,"\n Number of crossover of real variable = %d",nsga2Params->nrealcross);
        fprintf(fpt5,"\n Number of mutation of real variable = %d",nsga2Params->nrealmut);
    }
    if (nsga2Params->nbin!=0)
    {
        fprintf(fpt5,"\n Number of crossover of binary variable = %d",nsga2Params->nbincross);
        fprintf(fpt5,"\n Number of mutation of binary variable = %d",nsga2Params->nbinmut);
    }

    fflush(stdout);
    fflush(fpt1);
    fflush(fpt2);
    fflush(fpt3);
    fflush(fpt4);
    fflush(fpt5);
    fclose(fpt1);
    fclose(fpt2);
    fclose(fpt3);
    fclose(fpt4);
    fclose(fpt5);
    // if (nsga2Params->choice!=0)
    // {
    //     pclose(gp);
    // }
    if (nsga2Params->nreal!=0)
    {
        free (nsga2Params->min_realvar);
        free (nsga2Params->max_realvar);
    }
    if (nsga2Params->nbin!=0)
    {
        free (nsga2Params->min_binvar);
        free (nsga2Params->max_binvar);
        free (nsga2Params->nbits);
    }
    deallocate_memory_pop (nsga2Params,  parent_pop, nsga2Params->popsize);
    deallocate_memory_pop (nsga2Params,  child_pop, nsga2Params->popsize);
    deallocate_memory_pop (nsga2Params,  mixed_pop, 2*nsga2Params->popsize);
    free (parent_pop);
    free (child_pop);
    free (mixed_pop);


    free(out.nsga2_outs);
    free(out.residuals);


    printf("\n Routine successfully exited \n");
    return (0);
}


void print_nsga2Params(NSGA2Type *nsga2Params){
    int i;

    printf("NSGA2 Parameters:\n");
    printf("\n seed number is : %lf", nsga2Params->seed);
    printf("\n population size read is : %d",nsga2Params->popsize);
    printf("\n number of generations read is : %d",nsga2Params->ngen);
    printf("\n number of objectives entered is : %d",nsga2Params->nobj);
    printf("\n number of constraints entered is : %d",nsga2Params->ncon);
    printf("\n number of real variables entered is : %d",nsga2Params->nreal);
    printf("\n variables bounds: ");
    for (i=0; i<nsga2Params->nreal; i++){
            printf("[%lf", nsga2Params->min_realvar[i]);
            printf(" %lf], ", nsga2Params->max_realvar[i]);
    }
    printf("\n Probability of crossover entered is : %e",nsga2Params->pcross_real);
    printf("\n Probability of mutation entered is : %e",nsga2Params->pmut_real);
    printf("\n The eta_c entered is : %e",nsga2Params->eta_c);
    printf("\n The eta_m entered is : %e",nsga2Params->eta_m);
    // if (nsga2Params.nbin != 0){
    //     printf ("\n number of binary variables entered is : %d",nsga2Params.nbin);
    //     printf("\n Probability of crossover entered is : %e",nsga2Params.pcross_bin);
    //     printf("\n Probability of mutation entered is : %e",nsga2Params.pmut_bin);
    // }
}



