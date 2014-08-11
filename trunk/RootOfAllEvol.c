/** 
 * @file:   RootOfAllEvol.c
 * @author: Damjan Cicin-Sain
 * @email: damjan.cicin@crg.es
 *
 * Created on June 30, 2010, 12:10 PM
 */

#include "global.h"
/*#include "fly_io.h"*/
#include "mex.h"
#include "mathLib.h"            /* Trunc function is in there */    
#include "fly_sa.h"             /* MoveX function is in here */
#include <string.h>
#include "RootOfAllEvol.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /*printf("Mex Function Entered Successfully\n");*/
    static double *x, *y, *z, *m;
    static int *mask;
    static int i, size, rowLen, masksize;
    static int nParm = 0;
    static int init = 1;    /* init = 1 means 'initialization loop' */
    static int jacobian = 0;    /* set this if you need the Jacobian matrix */
    static Files files;
    static char *inputfile = NULL;

    static ScoreOutput out;
    
    out.score = 1e38;       /*start with a very large number*/
    out.penalty = 0;
    out.size_resid_arr = 0;
    out.jacobian = NULL;
    out.residuals = NULL;
    
    
    /* allocate memory for static file names */

    files.inputfile = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    files.statefile = ( char * ) calloc( MAX_RECORD, sizeof( char ) );

    inputfile = mxArrayToString(prhs[2]);
    /*printf("%s\n", inputfile);*/
        
    rowLen = mxGetN(prhs[0]);
    if (rowLen < nParm) {   /*for some strange and/or paranormal reason rowLen is 1 whenever local search is called, despite that the input array length is still 56 */
        rowLen = nParm;
    } else {
        nParm = rowLen;
    }
    
    masksize = mxGetN(prhs[1]);
    mask = (int *) calloc( masksize, sizeof(int) );

    x = mxGetPr(prhs[0]);
    m = mxGetPr(prhs[1]);

    for (i=0; i<masksize; i++) {
        mask[i] = (int) m[i];        
    }
    for (i=0; i<nParm; i++) {
        x[i] = Trunc(x[i], 5);      /*this is needed to avoid errors due to lack of precision*/
        /*if (x[i] > 100)
            printf("Warning - %lf > 100 for i = %d!!!\n", x[i], i);*/
    }
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    y = mxGetPr(plhs[0]);
       
    strcpy( files.inputfile, inputfile);
    
    sprintf( files.statefile, "%s.state", files.inputfile );
    /* check if a state file exists (access() is in unistd.h) */
    /*if( 0 == access( files.statefile, F_OK ) )
        stateflag = 1; //use this for restore when implemented*/
    
    MoveX(x, mask, &out, &files, init, jacobian, 1); /*0 Rkck, 1 Direct-Band */
    /*printf("SCORE = %lf, PENALTY = %lf, RETURNED %lf\n", out.score, out.penalty, out.score + out.penalty);*/
    /*printf("RETURNED %lf\n", out.score + out.penalty);*/
    if (out.score < 0) {    /* maybe eliminate later and deal only with 1e38 */
        printf("OUT_OF_BOUND - setting score to 0 and penalty to 1e38\n");
        out.score = 0;
        /*out.score = 1e38; /*testing (without penalties) purposes*/
        /*out.penalty = 0; /*testing (without penalties) purposes*/
        out.penalty = 1e38;
    }
    
    size = out.size_resid_arr;
    plhs[1] = mxCreateDoubleMatrix(1,size,mxREAL);
    z = mxGetPr(plhs[1]);
    for (i=0; i<size; i++) {
        z[i] = out.residuals[i];
    }
    y[0] = out.score + out.penalty; /*comment this to test without penalties*/
    
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    init = 0;
    mxFree(inputfile);
    free(files.inputfile);
    free(files.statefile);
    free(out.jacobian);
    free(out.residuals);
    free(mask);
    return;
}

