/**
 * @file moves.c                                               
 * @authors JR, modified by Yoginho
 *
 * @copyright Copyright (C) 1989-2003 John Reinitz, 2009-2013 Damjan Cicin-Sain, 
 * Anton Crombach and Yogi Jaeger
 *                                                               
 * @brief Function implementation for the move generation of 
 * the Lam-schedule annealer.
 * 
 * All this stuff is problem-dependent, but you shouldn't have to worry 
 * about the actual structure of the data in the model to be scored. That
 * is why moves.c communicates with model-specific stuff via translate.c. 
 * Funcs in this file also shouldn't use any Lam-specific stuff like 
 * \c SAType structs and such.
 */

#ifdef ICC
#include <mathimf.h>
#else
#include <math.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <error.h>
#include <moves.h>
#include <random.h>
//#include <distributions.h>               /* DistP.variables and prototypes */
//#include <sa.h>                /* *ONLY* for random number funcs and flags */
#include <fly_io.h>

#ifdef MPI                      
#include <mpi.h>                /* this is the official MPI interface */
#include <MPI.h>                /* our own structs and such only needed by parallel code */
#endif


/*** CONSTANTS *************************************************************/

const double THETA_MIN = 0;     /* minimum value of theta_bar (move size) */
const double THETA_INIT = 5.0;  /* initial value for all theta_bar (move size) */

const int LOWBITS = 0x330E;     /* following two are for drand to */
const int BYTESIZE = 8;         /* erand conversion */



/*** STATIC VARIABLES ******************************************************/

AParms ap;                      /* static copy of annealing parameters */

ParamList *ptab;                /* array of pointers to parameters and ranges */
AccStats *acc_tab;              /* struct to accumulate acceptance statistics */

int nparams;                    /* number of parameters to be tweaked */
int idx;                        /* index in ptab of thing-to-be-tweaked */

int nhits;                      /* number of moves since start of execution */
int nsweeps;                    /* number of sweeps since start of execution */

double pretweak;                /* used to restore param. value after rejection */

double old_energy = -999.;      /* two static vars used by Generate- */
double new_energy;              /* Move to calculate delta_e */

#ifdef MPI
static long *hits;              /* used for pooling number of moves */
static long *success;           /* used for pooling number of successes */
static long *tmp;               /* temp array for MPI_Allreduce sendbuf */
#endif

ScoreOutput out_local;

/*** FUNCTIONS *************************************************************/

/*** INITIALIZING AND RESTORING FUNCTIONS **********************************/

/**  InitMoves: initializes the following moves.c-specific stuff: 
 *              - static annealing parameter struct (ap)                   
 *              - static tweak struct (tweak) in translate.c               
 *              - initializes random number generator in lsa.c             
 *              - receives parameter list from Translate, stores nparams   
 *              - initializes acc_tab for acceptance statistics            
 *              - set mixing interval in lsa.c (parallel code only)        
 *                                                                         
 *              it then returns the initial temperature to the caller      
 */
double
InitMoves( FILE * fp, Input * inp ) {
    int i;                      /* local loop counter */

    /* following is used to initialize random number generator() */

    long seedval;               /* contains random number generator seed */

    /*
    int left;
    unsigned short left16;
    unsigned short middle16;
    unsigned short *xsubj;

    xsubj = ( unsigned short * ) calloc( 3, sizeof( unsigned short ) );
*/
    /* read annealing paramters and parameters-to-be-tweaked */

    ap = ReadAParameters( fp ); /* ap: static annealing params */
    ap.max_count = 0;

    if( equil == 1 ) {          /* read equilibration params and put them into lsa.c */
        SetEquilibrate( InitEquilibrate( fp ) );
    }

    /* initialze the random number generator, now dSFMT() */
#ifdef MPI
    seedval = ap.seed + myid;   /* each processor gets a different seed */
#else
    seedval = ap.seed;
#endif

   /* xsubj[0] = LOWBITS;

    middle16 = ( unsigned short ) seedval;
    xsubj[1] = middle16;

    left = seedval >> ( BYTESIZE * sizeof( unsigned short ) );
    left16 = ( unsigned short ) left;
    xsubj[2] = left16;*/

    InitRand( seedval );         /* makes the xsubj array static to lsa.c */

    /* Set up data structure for tweaking */
    nparams = inp->tra.size;    /* nparams is static to moves.c */
    ptab = inp->tra.array;      /* make parameter array static to moves.c */

    /* acc_tab is for statistics like acceptance ratio etc. */

    acc_tab = ( AccStats * ) calloc( nparams, sizeof( AccStats ) );

    for( i = 0; i < nparams; i++ ) {
        acc_tab[i].acc_ratio = 0;
        acc_tab[i].theta_bar = THETA_INIT;
        acc_tab[i].hits = 0;
        acc_tab[i].success = 0;
    }

#ifdef MPI

    /* allocate static arrays for parallel code */

    hits = ( long * ) calloc( nparams, sizeof( long ) );
    success = ( long * ) calloc( nparams, sizeof( long ) );
    tmp = ( long * ) calloc( nparams, sizeof( long ) );

#endif

    /* Finally, return the start temperature. */
    return ap.start_tempr;
}

/**  RestoreMoves: restores move generator from state file 
 *           NOTE: InitMoves will be called before this function during    
 *                 a restore                                               
 */
void
RestoreMoves( MoveState * MovePtr ) {
    int i;                      /* local loop counter */

    nparams = MovePtr->nparams; /* restore static stuff */
    idx = MovePtr->index;
    nhits = MovePtr->nhits;
    nsweeps = MovePtr->nsweeps;
    old_energy = MovePtr->old_energy;
    for( i = 0; i < nparams; i++ ) {    /* restore parameters */
        *( ptab[i].param ) = MovePtr->newval[i];
    }
    free( MovePtr->newval );
    free( acc_tab );            /* InitMoves has already been called. */
    acc_tab = MovePtr->acc_tab_ptr;     /* restore acceptance stats  */
    free( MovePtr );
}


/*** A FUNCTION FOR FINALIZING A RUN ***************************************/

/**  GetFinalInfo: collects stop energy and final count for output to the 
 *                 data file                                               
 */
AParms
GetFinalInfo( void ) {
    ap.stop_energy = old_energy;
    ap.max_count = nhits;

    return ap;
}




/*** MOVE GENERATION *******************************************************/

/** GenerateMove: evaluates the old energy, changes a parameter, then eval- 
 *               uates the new energy; returns the difference between old  
 *               and new energy to the caller                              
 */
double
GenerateMove( Files * files, DistParms * distp, ScoreOutput * out ) {
    double delta_e;             /* energy difference before and after move */
    /* for first call: check for valid parameters */

    if( old_energy == -999. ) {
        MoveSA( NULL, distp, out, NULL, 0, 0 );
        old_energy = out->score + out->penalty;
        if( old_energy == FORBIDDEN_MOVE )
            error( "GenerateMove: 1st call gave forbidden move" );
    }

    /* make a move, score and return either FORBIDDEN_MOVE or delta_e */

    Move( files, distp );

    acc_tab[idx].hits++;
    //new_energy = Score();

    MoveSA( NULL, distp, out, NULL, 0, 0 );

    new_energy = out->score + out->penalty;
    if( new_energy >= FORBIDDEN_MOVE ) {
        return ( FORBIDDEN_MOVE );
    } else {
        delta_e = new_energy - old_energy;

        //printf("%d) %lg - %lg = %lg ENERGYCHANGE\n", idx, new_energy, old_energy, delta_e);
        return ( delta_e );
    }
}

/** AcceptMove: sets new energy as the old energy for the next step and 
 *               keeps track of the number of successful moves          
 */
void
AcceptMove( void ) {
    old_energy = new_energy;
    acc_tab[idx].success++;
}

/**  RejectMove: simply resets the tweaked parameter to the pretweak value */
void
RejectMove( void ) {
    *( ptab[idx].param ) = pretweak;
}





/*** MOVE GENERATION - PART 2: FUNCS NEEDED IN MOVES.C (BUT NOT LSA.C) *****/

/** Move: tweaks the parameter-to-be-tweaked according to the current 
 *         value of index; also calls UpdateControl if necessary           
 */
void
Move( Files * files, DistParms * DistP ) {
    double tweakee;
    double theta;
    double sign;

    /* update counters */

    idx++;
    nhits++;

    idx = idx % nparams;

#ifdef MPI
    nsweeps = ( nhits / nparams ) * nnodes;
#else
    nsweeps = ( nhits / nparams );
#endif

    /* update statistics if interval passed & at least one sweep completed */
    if( !( nsweeps % ap.interval ) && !( idx ) && ( nsweeps ) ) {
        UpdateControl( files ); /* see comments in moves.h */
    }
    //printf("idx=%d param=%lg lower=%lg upper=%lg\n", idx, *(ptab[idx].param), ptab[idx].param_range->lower, ptab[idx].param_range->upper);
    tweakee = *( ptab[idx].param );
    pretweak = tweakee;


    /* this section does the tweaking */
    theta = generate_dev( acc_tab[idx].theta_bar, DistP );
    /* the sign stuff is needed for exponential distribution because always positive */
    /* may (but not likely) need in future if wish to evaluate poisson or pareto distributions for fly  */

    if( DistP->distribution == 1 ) {    /* need  positive + negative values for theta */
        sign = RandomReal(  ) - 0.5;
        if( sign <= 0 )
            theta = -theta;
    }
    //printf("idx %d theta %lg\n", idx, theta);

    tweakee += theta;

    *( ptab[idx].param ) = tweakee;     /* original eqparms in score.c tweaked */
    //printf("param = %lg\n", tweakee);
}

/** UpdateControl: each interval number of steps, acceptance stats are 
 *                  dated here; this function also prints prolix stuff, if 
 *                  required.                                              
 */
void
UpdateControl( Files * files ) {
    int i;                      /* local loop counter */
    double x;                   /* temp variable to manipulate theta_bar */

#ifdef MPI
    
    /* if parallel, pool the accpetance statistics */
    for( i = 0; i < nparams; i++ ) {
        hits[i] = ( long ) acc_tab[i].hits;
        success[i] = ( long ) acc_tab[i].success;
    }

    for( i = 0; i < nparams; i++ ) {
        tmp[i] = hits[i];
    }
    MPI_Allreduce( tmp, hits, nparams, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );

    for( i = 0; i < nparams; i++ ) {
        tmp[i] = success[i];
    }
    MPI_Allreduce( tmp, success, nparams, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );

    for( i = 0; i < nparams; i++ ) {
        acc_tab[i].hits = ( int ) hits[i];
        acc_tab[i].success = ( int ) success[i];
    }
#endif

    for( i = 0; i < nparams; i++ ) {

        acc_tab[i].acc_ratio = ( ( double ) acc_tab[i].success ) / ( ( double ) acc_tab[i].hits );
        x = log( acc_tab[i].theta_bar );
        x += ap.gain * ( acc_tab[i].acc_ratio - 0.44 );
        acc_tab[i].theta_bar = exp( x );
        if( acc_tab[i].theta_bar < THETA_MIN )
            acc_tab[i].theta_bar = THETA_MIN;

      
        acc_tab[i].hits = 0;
        acc_tab[i].success = 0;

    }

}


/*** FUNCTIONS THAT COMMUNICATE WITH OTHER SOURCE FILES ********************/

/** MoveSave: returns a MoveState struct in which the current state of 
 *             moves is saved; use for writing state file              
 */
MoveState *
MoveSave( void ) {
    MoveState *move_stuff;

    move_stuff = ( MoveState * ) malloc( sizeof( MoveState ) );

    move_stuff->old_energy = old_energy;        /* pretty straightforward, no? */
    move_stuff->pt = ptab;
    move_stuff->acc_tab_ptr = acc_tab;
    move_stuff->newval = NULL;
    move_stuff->index = idx;
    move_stuff->nhits = nhits;
    move_stuff->nparams = nparams;
    move_stuff->nsweeps = nsweeps;

    return move_stuff;
}



#ifdef MPI

/*** FUNCTIONS FOR SENDING MOVE STATES UPON MIXING (PARALLEL CODE ONLY) ****
 *   prototypes for these are in MPI.h; there's an extensive comment on    *
 *   how move state communication should be done at the beginning of lsa.c */

/** MakeStateMsg: function to prepare a move state message which is then 
 *                 passed to other nodes via MPI; lsa.c doesn't know about 
 *                 the structs we use for acceptance statistics in         
 *                 move(s).c, but we can safely assume that what we have   
 *                 to send can be communicated as longs and doubles; thus, 
 *                 we split the move state message in two arrays, one for  
 *                 the longs and one for the doubles; then we return the   
 *                 arrays and their sizes to lsa.c                         
 */
void
MakeStateMsg( long **longbuf, int *lsize, double **doublebuf, int *dsize ) {
    int i, n_l, n_d;



    /* calculate buffer size, compare with move parameters below */

    *lsize = 2 * nparams + 3;
    *dsize = 3 * nparams + 1;

    /* allocate buffer */

    *longbuf = ( long * ) calloc( *lsize, sizeof( long ) );
    *doublebuf = ( double * ) calloc( *dsize, sizeof( double ) );

    /* pack longs into their buffer */

    ( *doublebuf )[0] = old_energy;

    ( *longbuf )[0] = ( long ) idx;
    ( *longbuf )[1] = ( long ) nhits;
    ( *longbuf )[2] = ( long ) nsweeps;

    for( i = 0; i < nparams; i++ ) {

        n_l = i * 2 + 3;
        n_d = i * 3 + 1;

        ( *longbuf )[n_l] = ( long ) acc_tab[i].hits;
        ( *longbuf )[n_l + 1] = ( long ) acc_tab[i].success;

        ( *doublebuf )[n_d] = *( ptab[i].param );
        ( *doublebuf )[n_d + 1] = acc_tab[i].acc_ratio;
        ( *doublebuf )[n_d + 2] = acc_tab[i].theta_bar;

    }
}

/** AcceptMsg: gets the move state message from lsa.c and reinstalls acc- 
 *              eptance statistics into move.c; see the comment for Make-  
 *              StateMsg above for the rationale behind the two arrays     
 *              that are passed                                            
 */
void
AcceptStateMsg( long *longbuf, double *doublebuf ) {
    int i, n_l, n_d;

    idx = ( int ) longbuf[0];
    nhits = ( int ) longbuf[1];
    nsweeps = ( int ) longbuf[2];

    old_energy = doublebuf[0];

    for( i = 0; i < nparams; i++ ) {

        n_l = i * 2 + 3;
        n_d = i * 3 + 1;

        acc_tab[i].hits = ( int ) longbuf[n_l];
        acc_tab[i].success = ( int ) longbuf[n_l + 1];

        *( ptab[i].param ) = doublebuf[n_d];
        acc_tab[i].acc_ratio = doublebuf[n_d + 1];
        acc_tab[i].theta_bar = doublebuf[n_d + 2];

    }

    free( longbuf );
    free( doublebuf );

}
#endif
