/**
 * @file moves.h                                               
 * @authors JR, modified by Yoginho                          
 *                                                 
 * @copyright Copyright (C) 1989-2003 John Reinitz, 2009-2013 Damjan Cicin-Sain, 
 * Anton Crombach and Yogi Jaeger
 *
 * @brief Problem-specific stuff for the Lam annealer.
 *
 * There are functions for move generation and translation of    
 * model data structures into annealing parameter arrays. Moreover 
 * there are a couple of functions that read annealing and  
 * tune parameters and the tweak table from the data file.        
 */

#ifndef MOVES_INCLUDED
#define MOVES_INCLUDED

/* this def needed for func. defs that refer to (* FILE) */
#include <stdio.h>
/* following for structures & consts used thruout */
#include <global.h>
/* following for StopStyle enum style */
#include <sa.h>
/* ... and this for the Range struct */
#include <maternal.h>


/*** CONSTANTS *************************************************************/

/** minimum value of theta_bar (move size) */
extern const double THETA_MIN;  
/** initial value for all theta_bar (move size) */
extern const double THETA_INIT; 

/** following two are for drand to */
extern const int LOWBITS;       
/** erand conversion */
extern const int BYTESIZE;      


/*** TONS OF STRUCTS *******************************************************/

/** @brief Annealing parameters that are not specific to the Lam algorithm.
 *
 * In general the parameters should be used in moves.c or fly_sa.c but    
 * \em not in lsa.c. In the data file, the members of the struct labeled RO  
 * are read from the $annealing_input section. They are used for initial   
 * conditions of the annealer and do not change during a run. Members      
 * labeled OUT are written to the $annealing_output section upon completion 
 * of a run.                                                          
 */
typedef struct {
    /** seed for random number generator RO */
    long seed;                  
    /** the initial equilibration temperature RO */
    double start_tempr;         
    /** gain for proportional control of move size RO */
    double gain;                
    /** the final energy of the answer OUT */
    double stop_energy;         
    /** total number of iterations OUT */
    int max_count;              
    /** number of sweeps between updating theta_bar RO */
    int interval;               
    /*int distribution;    1 - uniform; 2 - exp; 3 - normal; 4 - lorentz RO */
} AParms;

/** @brief Acceptance statistics are used for keeping the acceptance 
 * ratio for each parameter as close as possible to 0.44.
 * 
 * There is one such struct for each parameter to be tweaked.           
 */
typedef struct {
    /** acceptance ratio for parameter */
    double acc_ratio;           
    /** theta bar is proportional to move size */
    double theta_bar;           
    /** number of moves since last call to \c UpdateControl() */
    int hits;                   
    /** number of these moves that were accepted */
    int success;                
} AccStats;


/* Some annealers need to know about the search space to make the right    
 * moves. For penalty type ranges, the penalty needs to be converted into  
 * explicit range limits using Penalty2Limits (in score.c). Probably, this 
 * should be called from UpdateControl and give the ranges corresponding   
 * to a defined penalty level                                               
 * BUT: PASSING PARAM LIMITS TO THE ANNEALER IS CURRENTLY NOT SUPPORTED    */
/*
typedef struct {
  double    *param;                // pointers to parameters to be tweaked
  Range     *param_range;        // pointers to corresponding range limits
} ParamList;


typedef struct {
  int       size;                           // size of the ParamList array
  ParamList *array;            // points to 1st element of ParamList array
  double    *pen_vec;   // penalty vector: see score.h, struct SearchSpace
} PArrPtr;
 */

/** @brief Contains copies of the static variables of moves.c and the values 
 * of parameters undergoing annealing.            
 */
typedef struct {
    /** Used during a save to point to annealed-on parameters */
    ParamList *pt;              
    /** points to current acceptance stats */
    AccStats *acc_tab_ptr;      
    /** points to array of annealed-on doubles for restore */
    double *newval;             
    /** energy before the last move */
    double old_energy;          
    /** # of parameters to be tweaked */
    int nparams;                
    /** index of parameter to be tweaked during a sweep */
    int index;                  
    /** number of moves already made */
    int nhits;                  
    /** number of completed sweeps */
    int nsweeps;                
} MoveState;


/* Tweak struct is for tweaking individual parameters or not; each pointer 
 * below points to an array of ints which represent each paramter          
 * 1 means tweak, 0 means leave it alone                                   
 * see ../doc/dataformatX.X for further details on the tweak section       */
/*
typedef struct {
 int       *Rtweak;                             // which Rs to be tweaked
 int       *Ttweak;                             // which Ts to be tweaked
 int       *Etweak;                             // which Es to be tweaked
 int       *mtweak;                             // which ms to be tweaked
 int       *htweak;                             // which hs to be tweaked
 int       *dtweak;                             // which ds to be tweaked
 int       *lambdatweak;                   // which lambdas to be tweaked
 int       *tautweak;                   // which taus to be tweaked
} Tweak;
 */

/** @brief Saving command line options in savestate.c 
 */
typedef struct {
    /** filename of input file */
    char *inname;               
    /** filename of output file */
    char *outname;              
    /** original command line */
    char *argv;                 
    /** derivative function */
    char *derivfunc;            
    /** solver function */
    char *solver;               
    /** stop criterion */
    StopStyle stop_flag;        
    /** prolix flag */
    int prolix_flag;            
    /** landscape flag */
    int landscape_flag;         
    /** flag for timing code */
    int time_flag;              
    /** log display flag */
    int log_flag;               
    /** frequency for writing state files */
    long state_write;           
    /** frequency for printing status to stdout */
    long print_freq;            
    /** opt for printing captions */
    long captions;              
    /** division style flag */
    int olddivstyle;            
    /** output precision */
    int precision;              
    /** solver step size */
    double stepsize;            
    /** flag for quenchit mode (T=0 immediately) */
    int quenchit;               
    /** flag for equilibration mode (T = const) */
    int equil;                  

#ifdef MPI
    /** flag for tuning mode */
    int tuning;                 
    /** index for sample interval (=covar_index * tau) */
    int covar_index;            
    /** how many times to write tune statistics */
    int write_tune_stat;        
    /** auto-stop tuning runs? */
    int auto_stop_tune;         
#endif
} Opts;




/*** FUNCTION PROTOTYPES ***************************************************/

/* fly_sa.c: I/O functions for miscellaneous stuff */

/*** ReadTune: reads the tune_parameters section in a data file and ********
 *             turns a SAType structure to the caller                      *
 ***************************************************************************/

//SAType ReadTune(FILE *fp);  //it's in fly_io.h now

/** ReadAParameters: reads the AParm struct from an annealing_input sec- 
 *                    tion; these are the annealing parameters that are    
 *                    not Lam-specific (and should NOT go into lsa.c       
 */
AParms ReadAParameters( FILE * fp );

/*** WriteAParameters: writes the aparm struct into a new section in the ***
 *                     file specified by filename; the new 'annealing_out- *
 *                     put section is inserted right after the 'tune_para- *
 *                     meters' section; to achieve this, we need to write  *
 *                     to a temporary file which is then renamed to the    *
 *                     output file name                                    *
 ***************************************************************************/

//void WriteAParameters(char *filename, AParms aparm);

/*** PrintAParameters: writes an 'annnealing_output' section with 'title' **
 *                     to the stream specified by fp                       *
 ***************************************************************************/

//void PrintAParameters(FILE *fp, AParms aparm, char *title);

/** WriteEquil: writes the equilibrate_variance section to the data file 
 *               right after the $equilibrate section                      
 */
void WriteEquil( char *filename, double *equil_var );


/** PrintEquil: writes an 'equilibrate_variance' section with 'title' 
 *               to the stream specified by fp                             
 */
void PrintEquil( FILE * fp, double *equil_var, char *title );


/** PrintTimes: writes two (parallel: three) times sections */
void PrintTimes( FILE * fp, double *delta );



/* functions that communicate with savestate.c */

/**  GetOptions: returns command line options to savestate.c 
 *               for the detailed meaning of all these options see Parse-  
 *               CommandLine() above); Opts struct defined in moves.h      
 */
Opts *GetOptions( void );

/**  RestoreOptions: restores the values of the command line opt variables 
 *                   from the Opts struct (used for restoring a run)       
 */
void RestoreOptions( Opts * options );


/**  Translate: creates an array of pointers that point to the parameters 
 *              to be tweaked; the PArrPtr that is returned also includes  
 *              pointers to the corresponding parameter ranges for each    
 *              parameter, although this feature is not yet used anywhere  
 *              in the annealing code                                      
 *     CAUTION: InitZygote and InitScoring have to be called first!        
 */
PArrPtr Translate( Input * inp );


/* moves.c: functions for move generation */

/* Initializing and restoring functions */

/** InitMoves: initializes the following moves.c-specific stuff: 
 *              - static annealing parameter struct (ap)                   
 *              - static tweak struct (tweak) in translate.c               
 *              - initializes random number generator in lsa.c             
 *              - receives parameter list from Translate, stores nparams   
 *              - initializes acc_tab for acceptance statistics            
 *              - set mixing interval in lsa.c (parallel code only)        
 *                                                                         
 *              it then returns the initial temperature to the caller      
 */
double InitMoves( FILE * fp, Input * inp );

/** RestoreMoves: restores move generator from state file 
 *           NOTE: InitMoves will be called before this function during    
 *                 a restore                                               
 */
void RestoreMoves( MoveState * MovePtr );

/** RestoreProlix: restores prolix file after a run has been interrupted */
void RestoreProlix( char *prolixfile );

/* a function for finalizing a run */

/** GetFinalInfo: collects stop energy and final count for output to the 
 *                 data file                                               
 */
AParms GetFinalInfo( void );

/* move generation functions that are used in moves.c, but not lsa.c */

/** Move: tweaks the parameter-to-be-tweaked according to the current 
 *         value of index; also calls UpdateControl if necessary           
 */
void Move( Files * files, DistParms * DistP );

/** UpdateControl: each interval number of steps, acceptance stats are 
 *                  dated here; this function also prints prolix stuff, if 
 *                  required.                                              
 */
void UpdateControl( Files * files );

/* functions that communicate with other source files */

/** SetProlix: sets flag for printing prolix output on acceptance stats 
 *              and initializes the prolix file if required                
 */
char *SetProlix( int value, char *file, int init_flag );

/** MoveSave: returns a MoveState struct in which the current state of 
 *             moves is saved; use for writing state file                  
 */
MoveState *MoveSave( void );

/* savestate.c */

/* funcs that read/write/erase the state file */

/**  StateRead: reads Lam statistics, move state and erand state from a 
 *              state file and restores the annealer's state to the same   
 *              state it was in before it got interrupted                  
 *     CAUTION: InitMoves must be called before calling StateRead!         
 */
void StateRead( char *statefile, Opts * options, MoveState * move_ptr, double *stats, char *rand, double *times );

/**  StateRm: removes the state file after the run has been completed; 
 *            unless we're tuning in parallel, only the root node needs to  
 *            delete a state file                                          
 */
void StateRm( void );

#endif
