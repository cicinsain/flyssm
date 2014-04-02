#include "amosa.h"
#include "fly_io.h"

void write_params_to_fly_output_standard(AMOSAType *amosaParams, Input *inp, char *inname){
    int i, j;
    ParamList *ptab;
    
    // printf("**********************\n%s\n", inname);
    
    printf("\nCreating output files for each parameters set.\n");
    
    ptab = inp->tra.array;      // Get the pointer to the 'inp' array.
    
    // printf("%d\n", amosaParams->i_archivesize);
    for(int i=0;i<amosaParams->i_archivesize;i++)
    {
        // printf("%d\n", amosaParams->i_totalno_var);
        for(int h=0;h<amosaParams->i_totalno_var;h++){
        	*( ptab[h].param  ) = amosaParams->d_archive[i][h];
        	// printf("%lf,", amosaParams->d_archive[i][h]);
            // fprintf(fp, "%f\t", amosaParams->d_func_archive[i][h]);
        }
        // printf("\n");
        WriteParameters(inname, &(inp->zyg.parm), "eqparms", 9, inp->zyg.defs);
        
    }//End of for loop
    
    
    // for (i = 0; i < amosaParams->popsize; ++i){
    //     /* code */
    //     for (j = 0; j < inp->tra.size; ++j){
    //         /* Modifying the local parameters of 'inp' struct */
    //         *( ptab[j].param  ) = pop->ind[i].xreal[j];
    //     }
    //     // Each run will create an output file in form of 'inname_parm_XXXXXX' in which 'XXXXXXX'
    //     // will replace by random string.
    //     WriteParameters(inname, &(inp->zyg.parm), "eqparms", 9, inp->zyg.defs);
    // }
    printf("Done.\n");
}