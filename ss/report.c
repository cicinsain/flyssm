#include "ss.h"
#include "fly_io.h"



/*
	Write the set to the file
 */
void write_set(SSType *ssParams, Set *set, int set_size, int member_length, FILE *fpt, int iter, char mode){

	for (int i = 0; i < set_size; ++i)
	{
		write_ind(ssParams, &(set->members[i]), member_length, fpt, iter, mode);
	}
	// fflush(fpt);

}

void write_ind(SSType *ssParams, individual *ind, int member_length, FILE *fpt, int iter, char mode){

	// fprintf(fpt, "%d\t", iter);
	if (iter != -1)
		fprintf(fpt, "%d\t", iter);
	
	for (int i = 0; i < member_length; ++i)
	{
		fprintf(fpt, "%.5lf\t", ind->params[i]);
	}
	fprintf(fpt, "%lf\n", ind->cost);
	// fprintf(fpt, "\n");


}


void print_set(SSType *ssParams, Set *set, int set_size, int member_length){
	printf("-----------------------------------\n");
	for (int i = 0; i < set_size; ++i)
	{
		printf("%d: ", i); print_ind(ssParams, &(set->members[i]), member_length);
	}
	printf("\n");

}

void print_ind(SSType *ssParams, individual *ind, int member_length){

	for (int i = 0; i < member_length; ++i)
	{
            printf("%.5lf, ", ind->params[i]);
	}	
	//printf("\t (cost: %lf)\t(distance: %lf)\n", ind->cost, ind->distance);
	printf("\t (cost: %lf)\n", ind->cost/*, ind->distance*/);

}

/*
	Print a subset list
*/
void print_subsets_list(SSType *ssParams){

	for (int i = 0; i < ssParams->subsets_list_size; ++i)
	{
		printf("[i: %d]\n", i);
		print_set(ssParams, &(ssParams->subsets_list[i]), ssParams->pair_size, ssParams->nreal);
	}
}

void print_double_matrix(SSType *ssParams, double **matrix, int row, int col){
	for (int r=0; r<row; r++)
	{
		printf("%d: ", r);
	    for(int c=0; c<col; c++)
	         printf("%.4lf     ", matrix[r][c]);
	    printf("\n");
	 }
	 printf("===========================================\n");

}

void print_int_matrix(SSType *ssParams, int **matrix, int row, int col){
	for (int r=0; r<row; r++)
	{
		printf("%d: ", r);
	    for(int c=0; c<col; c++)
	         printf("%d\t", matrix[r][c]);
	    printf("\n");
	 }
	 printf("===========================================\n");

}

// Process has done i out of n rounds,
// and we want a bar of width w and resolution r.
void loadBar(int x, int n, int r, int w)
{
    // Only update r times.
    if ( x % (n/r +1) != 0 ) return;
 
    // Calculuate the ratio of complete-to-incomplete.
    float ratio = x/(float)n;
    int   c     = ratio * w;
 
    // Show the percentage complete.
    printf("%3d%% [", (int)(ratio*100) );
 
    // Show the load bar.
    for (int x=0; x<c; x++)
       printf("=");
 
    for (int x=c; x<w; x++)
       printf(" ");
 
    // ANSI Control codes to go back to the
    // previous line and clear it.
    printf("]\n\033[F\033[J");
}

void write_int_matrix(SSType *ssParams, int **matrix, int row, int col, FILE *fpt, int iter, char mode){

	for (int r=0; r<row; r++)
	{
		if (iter != -1)
			fprintf(fpt, "%d", iter);
		
	    for(int c=0; c<col; c++)
	         fprintf(fpt, "\t%d", matrix[r][c]);
	    fprintf(fpt, "\n");
	 }
}

void write_double_matrix(SSType *ssParams, double **matrix, int row, int col, FILE *fpt, int iter, char mode){

	for (int r=0; r<row; r++)
	{
		if (iter != -1)
			fprintf(fpt, "%d", iter);
		
	    for(int c=0; c<col; c++)
	         fprintf(fpt, "\t%lf", matrix[r][c]);
	    fprintf(fpt, "\n");
	 }
}

void write_params_to_fly_output_standard(SSType *ssParams, Input *inp, char *inname){
    int i, j;
    ParamList *ptab;
    
    printf("**********************\n%s\n", inname);
    
    printf("\nCreating output files for each parameters set.\n");
    
    ptab = inp->tra.array;      // Get the pointer to the 'inp' array.
    
    // printf("%d\n", ssParams->i_archivesize);
    // printf("%d\n", ssParams->i_hardl);
    for(int i=0;i<ssParams->ref_set_size;i++)
    {
        // printf("%d\n", ssParams->i_totalno_var);
        for(int h=0;h<ssParams->nreal;h++){
        	*( ptab[h].param  ) = ssParams->ref_set->members[i].params[h];
            // printf("%lf,", ssParams->d_archive[i][h]);
         //    printf("---%lf", ssParams->d_func_archive[i][h]);
        	// printf("---%lf\t", ssParams->d_solution[i][h]);

            // fprintf(fp, "%f\t", ssParams->d_func_archive[i][h]);
        }
        // printf("\n");
        WriteParameters(inname, &(inp->zyg.parm), "eqparms", 9, inp->zyg.defs);
        
    }//End of for loop
    
    
    // for (i = 0; i < ssParams->popsize; ++i){
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