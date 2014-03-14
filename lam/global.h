/**
 *                                                               
 *   @file global.h                                              
 *                                                               
 *****************************************************************
 *                                                               
 *   written by JR, modified by Yoginho                          
 *                                                               
 *****************************************************************
 *                                                               
 *   THIS HEADER CONTAINS *ONLY* STUFF THAT IS NEEDED IN BOTH    
 *   LSA.C AND COST FUNCTION-SPECIFIC FILES.                     
 *                                                               
 *****************************************************************
 *                                                               
 *   PLEASE THINK TWICE BEFORE PUTTING ANYTHING IN HERE!!!!      
 *                                                               
 *****************************************************************
 *                                                               
 *   global.h gets included automatically by other header files  
 *   if needed and is not included explicitly in .c files.       
 *                                                               
 */


/*==========
	Amir
	*/
    // #define NSGA2


#ifndef GLOBAL_INCLUDED
#define GLOBAL_INCLUDED

#include <float.h>

// #ifdef AMOSA
//     #include "amosa.h"
// #elif NSGA2
//     #include "nsga2.h"
// #endif


/*** A global in global for debugging **************************************/

int debug;                      /* debugging flag */
int proc_id;

/*** Constants *************************************************************/

/** IMPORTANT NOTE: the following used to be only 80 to be backward 
 *                compatible with punch cards that were of that maximum length.
*                 We have decided to abandon punch-card compatibility,    
*                 since few people are actually using them anymore...     
*                 It was a tough decision though.      JR & JJ, July 2001 
*/ 
#define MAX_RECORD 1024          /* max. length of lines read from file */

/** The following defines the maximum float precision that is supported by  
 * the code.
 */
extern const int MAX_PRECISION;

/** the following constant as a score tells the annealer to reject a move,  
* no matter what. It had better not be a number that could actually be a  
* score.
*/
extern const double FORBIDDEN_MOVE;     /* the biggest possible score, ever */

/** out of bound control in score */
extern const int OUT_OF_BOUND;  

/** The whole output */
typedef struct ScoreOutput {    
    int size_resid_arr;
    double score;
    double penalty;
    double *residuals;
    double *jacobian;
    /*============
    	Amir: adding array of scores for multi-objective optimization purpose.
    	*/
    // #if defined(AMOSA) || defined(NSGA2)
    //      TODO: Found the problem here: If I use the #if then the code won't work. It seems that the NSGA2 is not defined here...
    	struct ScoreOutput **nsga2_outs;
    // #endif
} ScoreOutput;


#endif
