/**
 * @file fly_nsga2.c                                                
 * @author JR, modified by Yoginho,                            
 *   -D landscape option by Lorraine Greenwald in Oct 2002,            
 *   -g option by Yousong Wang in Feb 2002,        
 *   -a option by Marcel Wolf in Apr 2002                          
 *
 * @copyright Copyright (C) 1989-2003 John Reinitz, 2009-2013 Damjan Cicin-Sain, 
 * Anton Crombach and Yogi Jaeger
 * 
 * @brief Although main() is in lsa.c, this is the file that 'defines'    
 * the fly_nsga2 program.
 *
 * It contains most of its problem-      
 * specific code (except for move generation -> moves.c, saving    
 * of intermediate state files -> savestate.c and communication    
 * with the specific cost function that is used -> translate.c).   
 *                                                                 
 * After I've told you all that's NOT in this file, here's what    
 * the funcs below actually do: parsing fly_nsga2 command line opts   
 * is one of its jobs; there are funcs that make the first and     
 * last moves and funcs that read and write Lam and Lam-independent 
 * annealing parameters to the problem-specific data file.    
 */


#include <float.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>             /* for command line option stuff */

/*======
    Utils*/
#include "error.h"              /* error handling funcs */
#include "distributions.h"      /* DistP.variables and prototypes */
#include "integrate.h"
#include "random.h"             /* for InitRand() */
/*=====
    SA*/
#include "sa.h"                 /* problem-independent annealing funcs */
 
/*======
    Solver*/
#include "maternal.h"           /* for olddivstyle and such */
#include "moves.h"              /* problem-specific annealing funcs */
#include "score.h"              /* for init and Score funcs */
#include "solvers.h"            /* for name of solver funcs */
#include "zygotic.h"            /* for init, mutators and derivative funcs */
#include "fly_io.h"

/*=================
    NSGA2
    */
# include "amosa.h"



#ifdef MPI
#include <mpi.h>                /* this is the official MPI interface */
#include <MPI.h>                /* our own structs and such only needed by parallel code */
#endif

/*========
    defining the NSGA2Type paramters...
    */
#ifdef NSGA2
NSGA2Type nsga2Params;
#endif

AMOSAType amosaParams;

/*** Constants *************************************************************/

/* command line option string */

/* #define  OPTS       ":a:b:Bc:C:d:De:Ef:g:hi:lLnopQr:s:StTvw:W:y:" */

const char *OPTS = ":a:b:Bc:C:De:Ef:g:hi:lLm:nNopQr:s:StTvw:W:y:";
/* command line option string */
/* D will be debug, like scramble, score */
/* must start with :, option with argument must have a : following */


/*** STATIC VARIABLES ******************************************************/

/* Help, usage and version messages */

#ifdef MPI
static const char usage[] =
    "Usage: fly_nsga2.mpi [-b <bkup_freq>] [-B] [-C <covar_ind>] \n"
    "                  [-D] [-e <freeze_crit>][-E] [-f <param_prec>] [-g <g(u)>]\n"
    "                  [-h] [-i <stepsize>] [-l] [-L] [-n] [-N] [-p] [-s <solver>]\n"
    "                  [-S] [-t] [-T] [-v] [-w <out_file>]\n" "                  [-W <tune_stat>] [-y <log_freq>]\n" "                  <datafile>\n";
#else
static const char usage[] =
    "Usage: fly_nsga2 [-a <accuracy>] [-b <bkup_freq>] [-B] [-e <freeze_crit>] [-E]\n"
    "              [-f <param_prec>] [-g <g(u)>] [-h] [-i <stepsize>] [-l] [-L] \n"
    "              [-m <score_method>] [-n] [-N] [-p] [-Q] [-s <solver>] [-t] [-v]\n"
    "              [-w <out_file>] [-y <log_freq>]\n" "              <datafile>\n";
#endif

static const char help[] =
#ifdef MPI
    "Usage: fly_nsga2.mpi [options] <datafile>\n\n"
#else
    "Usage: fly_nsga2 [options] <datafile>\n\n"
#endif
    "Argument:\n"
    "  <datafile>          input data file\n\n"
    "Options:\n"
    "  -a <accuracy>       solver accuracy for adaptive stepsize ODE solvers\n"
    "  -b <bkup_freq>      write state file every <bkup_freq> * tau moves\n" "  -B                  run in benchmark mode (only do fixed initial steps)\n"
#ifdef MPI
    "  -C <covar_ind>      set covar sample interval to <covar_ind> * tau\n"
#endif
    "  -D                  debugging mode, prints all kinds of debugging info\n"
    "  -e <freeze_crit>    set annealing freeze criterion to <freeze_crit>\n"
    "  -E                  run in equilibration mode\n"
    "  -f <param_prec>     float precision of parameters is <param_prec>\n"
    "  -g <g(u)>           chooses g(u): e = exp, h = hvs, s = sqrt, t = tanh\n"
    "  -h                  prints this help message\n"
    "  -i <stepsize>       sets ODE solver stepsize (in minutes)\n" "  -l                  echo log to the terminal\n"
#ifdef MPI
    "  -L                  write local logs (llog files)\n"
#endif
    "  -m <score_method>   w = wls, o=ols score calculation method\n"
    "  -n                  nofile: don't print .log or .state files\n"
    "  -N                  generates landscape to .landscape file in equilibrate mode \n"
    "  -o                  use oldstyle cell division times (3 div only)\n" "  -p                  prints move acceptance stats to .prolix file\n"
#ifndef MPI
    "  -Q                  quenchit mode, T is lowered immediately to zero\n"
#endif
    "  -s <solver>         choose ODE solver\n"
#ifdef MPI
    "  -S                  disable tuning stop flag\n"
#endif
    "  -t                  write timing information to .times file\n"
#ifdef MPI
    "  -T                  run in tuning mode\n"
#endif
    "  -v                  print version and compilation date\n" "  -w <out_file>       write output to <out_file> instead of <datafile>\n"
#ifdef MPI
    "  -W <tune_stat>      tuning stats written <tune_stat> times per interval\n"
#endif
    "  -y <log_freq>       write log every <log_freq> * tau moves\n\n" "Please report bugs to <yoginho@usa.net>. Thank you!\n";

static char version[MAX_RECORD];        /* version gets set below */

static const char verstring[] =
    "compiled by:      %s\n" "         on:      %s\n" "      using:      %s\n" "      flags:      %s\n" "       date:      %s at %s\n";

/* other static variables */

static char *inname;            /* filename of input file */
static char *outname;           /* filename of output file */

static char *argvsave;          /* static string for saving command line */

static double stepsize = 1.;    /* stepsize for solver */
static double accuracy = 0.001; /* accuracy for solver (not used yet) */
static int precision = 8;       /* precision for eqparms */
static int prolix_flag = 0;     /* to prolix or not to prolix */
static int landscape_flag = 0;  /* generate energy landscape data */
static int method = 0;          /* 0 for wls, 1 for ols */
/* set the landscape flag (and the landscape filename) in lsa.c */

/* the following lines define a pointers to:                               */
/*            - pd:    dvdt function, Dvdt_sqrt or DvdtOrig in zygotic.c   */
/*            - pj:    Jacobian function, in zygotic.c                     */
/* This pointer needs to be static since both InitialMove and Restore-     */
/* State need it; it gets set in ParseCommandLine()                        */
/*                                                                         */
/* NOTE: ps (solver) is declared as global in integrate.h                  */

static Input inp;               // The whole input - this is static in order to not to read data from file at every loop

void ( *pd ) ( double *, double, double *, int, SolverInput *, Input * );
void ( *pj ) ( double, double *, double *, double **, int, SolverInput *, Input * );






/*** FUNCTIONS *************************************************************/

/*** COMMAND LINE OPTS ARE ALL PARSED HERE *********************************/

/**  ParseCommandLine: well, parses the command line and returns an index  
 *                     to the 1st argument after the command line options  
 */
int
ParseCommandLine( int argc, char **argv ) {






    int c, i;                   /* used to parse command line options */

    /* external declarations for command line option parsing (unistd.h) */
    extern char *optarg;        /* command line option argument */
    extern int optind;          /* pointer to current element of argv */
    extern int optopt;          /* contain option character upon error */
    /* set the version string */



#ifdef MPI
    sprintf( version, "fly_nsga2 version %s parallel", VERS );
#else
#ifdef ALPHA_DU
    sprintf( version, "fly_nsga2 version %s serial dec-kcc-dxml", VERS );
#else
//    sprintf(version, "fly_nsga2 version %s serial", VERS);
#endif
#endif

    /* following part sets default values for command line options */

    pd = DvdtOrig;              /* old default derivative function */
    //pd = Dvdt_sqrt; /* default derivative function */
    pj = JacobnOrig;            /* default Jacobian function */
    ps = Rkck;                  /* default solver */
    dd = DvdtDelay;             /* delayed derivative fnuction */

    captions = 100000000;       /* default freq for writing captions (off) */
    print_freq = 100;           /* default freq for writing to log file */
    state_write = 1000;         /* default freq for writing to state file */

    stop_flag = absolute_freeze;        /* type of stop criterion */
    time_flag = 0;              /* flag for timing */
    log_flag = 0;               /* flag for writing logs to the screen */
    nofile_flag = 0;            /* flog for not writing .log and .state files */

#ifdef MPI
    tuning = 0;                 /* tuning mode is off by default */
    covar_index = 1;            /* covariance sample will be covar_index * tau */
    write_tune_stat = 1;        /* how many times do we write tuning stats? */
    auto_stop_tune = 1;         /* auto stop tuning runs? default: on */
    write_llog = 0;             /* write local llog files when tuning; default: off */
#endif

    /* following part parses command line for options and their arguments      */
    optarg = NULL;
    while( ( c = getopt( argc, argv, OPTS ) ) != -1 ) {
        switch ( c ) {
        case 'a':
            accuracy = atof( optarg );
            if( accuracy <= 0 )
                error( "fly_nsga2: accuracy (%g) is too small", accuracy );
            break;
        case 'b':              /* -b sets backup frequency (to write state file) */
            state_write = strtol( optarg, NULL, 0 );
            if( state_write < 1 )
                error( "fly_nsga2: max. backup frequency is every tau steps i.e. -b 1" );
            if( state_write == LONG_MAX )
                error( "fly_nsga2: argument for -b too large" );
            break;
        case 'B':              /* -B runs in benchmark mode -> quit after intial steps */
            bench = 1;
            time_flag = 1;
            break;
        case 'c':              /* -c sets the frequency for printing captions */
            error( "fly_nsga2: -c is not supported anymore, captions off for good" );
            /* if you want to be able to insert captions into .log files, uncomment    *
             * the following lines; CAUTION: make sure that RestoreLog() works proper- *
             * ly with captions before you do this!                                    */
            /*      captions = strtol(optarg, NULL, 0);
               if ( captions < 1 )
               error("fly_nsga2: can't print captions more than every line (-c 1)");
               if ( captions == LONG_MAX )
               error("fly_nsga2: argument for -c too large");                        */
            break;
        case 'C':              /* parallel code: -C sets the covar sample index */
#ifdef MPI
            /* -C is currently disabled since I (Yogi) could not get it to work and I  *
             * got tired of debugging and I didn't really need it at that time and so  *
             * on and so forth; try running on a smaller number of nodes or get it to  *
             * work yourself; that's what I say. Grmbl.                                */
            /*    covar_index = atoi(optarg); 
               if ( covar_index < 1 )
               error("fly_nsga2: covariation sample index must be >= 1");            */
            error( "fly_nsga2: -C does not work yet, try to run on fewer nodes" );
#else
            error( "fly_nsga2: can't use -C in serial, tuning only in parallel" );
#endif
            break;
        case 'D':
            debug = 1;
            break;
        case 'e':              /* -e sets the stopping criterion */
            if( !( strcmp( optarg, "pfreeze" ) ) )
                stop_flag = proportional_freeze;
            else if( !( strcmp( optarg, "afreeze" ) ) )
                stop_flag = absolute_freeze;
            else if( !( strcmp( optarg, "abs" ) ) )
                stop_flag = absolute_energy;
            else
                error( "fly_nsga2: valid stopping criteria are pfreeze, afreeze, abs" );
            break;
        case 'E':              /* -E does equilibration runs */
            equil = 1;
            break;
        case 'f':
            precision = atoi( optarg ); /* -f determines float precision */
            if( precision < 0 )
                error( "fly_nsga2: what exactly would a negative precision be???" );
            if( precision > MAX_PRECISION )
                error( "fly_nsga2: max. float precision is %d!", MAX_PRECISION );
            break;
        case 'g':              /* -g choose g(u) function */
            pd = DvdtOrig;      //
            if( !( strcmp( optarg, "s" ) ) )
                gofu = Sqrt;
            else if( !( strcmp( optarg, "t" ) ) )
                gofu = Tanh;
            else if( !( strcmp( optarg, "e" ) ) )
                gofu = Exp;
            else if( !( strcmp( optarg, "h" ) ) )
                gofu = Hvs;
            else if( !( strcmp( optarg, "k" ) ) ) {
                gofu = Kolja;
            } else
                error( "fly_nsga2: %s is an invalid g(u), should be e, h, s or t", optarg );
            break;
        case 'h':              /* -h help option */
            PrintMsg( help, 0 );
            break;
        case 'i':              /* -i sets the stepsize */
            stepsize = atof( optarg );
            if( stepsize < 0 )
                error( "fly_nsga2: going backwards? (hint: check your -i)" );
            if( stepsize == 0 )
                error( "fly_nsga2: going nowhere? (hint: check your -i)" );
            if( stepsize > MAX_STEPSIZE )
                error( "fly_nsga2: stepsize %g too large (max. is %g)", stepsize, MAX_STEPSIZE );
            break;
        case 'l':              /* -l displays the log to the screen */
            log_flag = 1;
            break;
        case 'L':              /* -L writes local .llog files when tuning */
#ifdef MPI
            write_llog = 1;
#else
            error( "fly_nsga2: can't use -L in serial, tuning only in parallel" );
#endif
            break;
        case 'm':              /* -m sets the score method: w for wls, o for ols */
            if( !( strcmp( optarg, "w" ) ) )
                method = 0;
            else if( !( strcmp( optarg, "o" ) ) )
                method = 1;
            break;
        case 'n':              /* -n suppresses .state and .log files */
            nofile_flag = 1;
            break;
        case 'N':              /* -N sets laNdscape flag and Equilibrate mode */
            equil = 1;
            landscape_flag = 1;
            break;
        case 'o':              /* -o sets old division style (ndivs = 3 only! ) */
            olddivstyle = 1;
            break;
        case 'p':              /* -p sets prolix printing */
            prolix_flag = 1;
            break;
        case 'Q':              /* -Q sets quenchit mode (serial code only) */
#ifdef MPI
            error( "fly_nsga2: can't use -Q in parallel, quenchit only in serial" );
#else
            quenchit = 1;
#endif
            break;
        case 'r':
            error( "fly_nsga2: -r is currently not supported, use -g instead" );
            break;
        case 's':              /* -s sets solver to be used */
            if( !( strcmp( optarg, "a" ) ) )
                ps = Adams;
            else if( !( strcmp( optarg, "bd" ) ) )
                ps = BaDe;
            else if( !( strcmp( optarg, "bs" ) ) )
                ps = BuSt;
            else if( !( strcmp( optarg, "e" ) ) )
                ps = Euler;
            else if( !( strcmp( optarg, "h" ) ) )
                ps = Heun;
            else if( !( strcmp( optarg, "mi" ) ) || !( strcmp( optarg, "m" ) ) )
                ps = Milne;
            else if( !( strcmp( optarg, "me" ) ) )
                ps = Meuler;
            else if( !( strcmp( optarg, "r4" ) ) || !( strcmp( optarg, "r" ) ) )
                ps = Rk4;
            else if( !( strcmp( optarg, "r2" ) ) )
                ps = Rk2;
            else if( !( strcmp( optarg, "rck" ) ) )
                ps = Rkck;
            else if( !( strcmp( optarg, "rf" ) ) )
                ps = Rkf;
            else if( !( strcmp( optarg, "sd" ) ) )
                ps = SoDe;
            else if( !( strcmp( optarg, "kr" ) ) )
                ps = Krylov;
            /* else if (!(strcmp(optarg, "bnd")))
               ps = Band; */
            else
                error( "fly_nsga2: invalid solver (%s), use: a,bs,e,h,kr,mi,me,r{2,4,ck,f}", optarg );
            break;
        case 'S':              /* -S unsets the auto_stop_tune flag */
#ifdef MPI
            auto_stop_tune = 0;
#else
            error( "fly_nsga2: can't use -S in serial, tuning only in parallel" );
#endif
            break;
        case 't':              /* -t saves times in data and .times files */
            time_flag = 1;
            break;
        case 'T':              /* -T does tuning (parallel code only) */
#ifdef MPI
            tuning = 1;
#else
            error( "fly_nsga2: can't use -T in serial, tuning only in parallel" );
#endif
            break;
        case 'v':              /* -v prints version message */
            fprintf( stderr, "%s\n", version );
            //fprintf(stderr, verstring, USR, MACHINE, COMPILER, FLAGS, __DATE__, __TIME__);
            exit( 0 );
        case 'w':              /* -w sets output file */
            outname = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
            outname = strcpy( outname, optarg );
            SetOutname( outname );      /* communicates outname to lsa.c */
            break;
        case 'W':              /* -W determines the frequency of writing tuning stats */
#ifdef MPI
            write_tune_stat = atoi( optarg );
            if( write_tune_stat < 1 )
                error( "fly_nsga2: frequency of writing tune stats must be >= 1" );
#else
            error( "fly_nsga2: can't use -W in serial, tuning only in parallel" );
#endif
            break;
        case 'y':              /* -y set frequency to print log */
            print_freq = strtol( optarg, NULL, 0 );
            if( print_freq < 1 )
                error( "fly_nsga2: can't print status more than every tau steps (-y 1)" );
            if( print_freq == LONG_MAX )
                error( "fly_nsga2: argument for -y too large" );
            break;
        case ':':
            error( "fly_nsga2: need an argument for option -%c", optopt );
            break;
        case '?':
        default:
            error( "fly_nsga2: unrecognized option -%c", optopt );
        }
    }

    /* error checking here */
#ifdef MPI
    if( ( tuning == 1 ) && ( equil == 1 ) )
        error( "fly_nsga2: can't combine -E with -T" );
    if( write_llog && !tuning )
        error( "fly_nsga2: -L only makes sense when tuning" );
#else
    if( ( quenchit == 1 ) && ( equil == 1 ) )
        error( "fly_nsga2: can't combine -E with -Q" );
#endif
    if( ( ( argc - ( optind - 1 ) ) != 2 ) )
        PrintMsg( usage, 1 );

    argvsave = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    for( i = 0; i < argc; i++ ) {
        if( i > 0 )
            argvsave = strcat( argvsave, " " );
        argvsave = strcat( argvsave, argv[i] );
    }
    argvsave = strcat( argvsave, "\n" );

    return optind;
}

/** MoveSA: This function actually does almost everything.
 * First it creates a static Input structure 'inp', where it puts all the 
 * information from the input file. This part is executed only once (when init == 1).
 * Then it creates a ScoreOutput structure 'out' where the score, penalty 
 * and residual vectors will be stored.
 * At the end it runs the score function, where all the calculation is done.
 */
double
MoveSA( NucStatePtr state_ptr, DistParms * distp, ScoreOutput * out, Files * files, int init, int jacobian ) {
    //printf("MoveSA Function Entered Successfully\n");
    char *p;

    double i_temp;
    i_temp = 0.0;
    
    if( init == 1 ) {           //initialization loop - this part is done only the first time, and it reads information from the input file, and it saves
        //all that info in various subsections of the structure 'inp'
        FILE *infile;
        FILE *slogfile = NULL;  //solver log file where Blastoderm function writes it's output (only in debugging mode)


        SAType in_tune;         // temporary struct to read in tune parameters

        /* save some things in static storage for FinalMove() and StateWrite() */

        inname = ( char * ) calloc( MAX_RECORD + 1, sizeof( char ) );   // input file
        inname = strcpy( inname, files->inputfile );

        if( !outname ) {        // in case there's no -w, outname = inname
            outname = ( char * ) calloc( MAX_RECORD + 1, sizeof( char ) );
            outname = strcpy( outname, inname );
        }

        /* set the prolix flag (and the .prolix filename) in moves.c, if necessary */

        if( prolix_flag ) {
            files->prolixfile = SetProlix( prolix_flag, outname, 1 );   /* 1 means: new .prolix file */
        }

        if( landscape_flag )    /* 1 means gen landscape files */
            InitLandscape( landscape_flag, outname );   /*lives in lsa.c */
        /* initialize some Lam/Greening structures */
        p = state_ptr->tune.progname;   /* tune.progname contains program name */
        p = strcpy( p, version );

        infile = fopen( inname, "r" );
        if( !infile )
            file_error( "fly_nsga2 error opening input file" );
        if( debug ) {
            slogfile = fopen( strcat( inname, ".slog" ), "w" );
            if( !slogfile )
                file_error( "fly_*** error opening slog file" );        //check if this message is ok for this kind of error*/
        }
        inp.zyg = InitZygote( infile, pd, pj, &inp, "input" );
        inp.sco = InitScoring( infile, method, &inp );
        inp.his = InitHistory( infile, &inp );  //It fills the polations vector
        inp.ext = InitExternalInputs( infile, &inp );
        inp.ste = InitStepsize( stepsize, accuracy, slogfile, inname );
        // read the list of parameters to be tweaked
        inp.twe = InitTweak( infile, NULL, inp.zyg.defs );
        inp.tra = Translate( &inp );
        in_tune = ReadSATune( infile ); /* read tune_parameter section */

        // Reads all the nsga2 parameters from the input file, and translate the variables boundaries 
        // from SearchSpace to the nsga2.max/min
        // #ifdef NSGA2
        // nsga2Params = ReadNSGA2Parameters(infile, &inp);    
        // #endif

        amosaParams = ReadAMOSAParameters(infile, &inp);

        i_temp = InitMoves( infile, &inp );     /* set initial temperature and initialize */
        // initialize distribution stuff
        inp.dis = InitDistribution( infile );
        memcpy( distp, &( inp.dis ), sizeof( DistParms ) );
        fclose( infile );
        /* initialize Lam parameters (see sa.h for further detail) */

        state_ptr->tune.lambda = in_tune.lambda;
        state_ptr->tune.lambda_mem_length_u = in_tune.lambda_mem_length_u;
        state_ptr->tune.lambda_mem_length_v = in_tune.lambda_mem_length_v;
        state_ptr->tune.control = in_tune.control;
        state_ptr->tune.initial_moves = in_tune.initial_moves;
        state_ptr->tune.tau = in_tune.tau;
        state_ptr->tune.freeze_count = in_tune.freeze_count;
        state_ptr->tune.update_S_skip = in_tune.update_S_skip;
        state_ptr->tune.criterion = in_tune.criterion;
#ifdef MPI
        state_ptr->tune.mix_interval = in_tune.mix_interval;
#endif

    }

    inp.lparm = CopyParm( inp.zyg.parm, &( inp.zyg.defs ) );

    // #ifdef NSGA2
    //     InitNSGA2(&inp, &nsga2Params, inname);
    //     RunNSGA2(&inp, &nsga2Params, inname);
    // #endif
        // ReadAMOSAParams(&amosaParams);
        InitAMOSA(&inp, &amosaParams);
        RunAMOSA(&inp, &amosaParams, inname);


    // Ignore the Score function in order to avoid running the SA
#ifndef NSGA2
    //In this function all the calculations are made
    // Score( &inp, out, jacobian );
#endif
    FreeMutant( inp.lparm );


    //        exit(1);
    
    return ( i_temp );          /* initial temperature */

}

/** RestoreState: called when an interrupted run is restored; does the 
 *                 following:                                              
 *                 - stores various static things for filenames and such   
 *                 - initializes Lam parameters for lsa.c                  
 *                 - initializes model and scoring funcs & solver stepsize 
 *                 - initializes move generation in moves.c                
 *                 - restores state at which previous run was interrupted  
 *                   according to state file                               
 *                                                                         
 * Comment by JR:  RestoreState was originally going to be implemented with
 * branches in InitialMove. That won't work because when when this func.   
 * returns, control should go right to Loop(), skipping all the initiali-  
 * zation stuff in Initialize. Hence most of the code in InitialMove is    
 * just repeated here.                                                     
 */
void
RestoreState( NucStatePtr state_ptr, DistParms * distp, Files * files ) {
    FILE *infile, *slogfile = NULL;     /* pointer to original input file */
    char *p;                    /* temporary string */
    double i_temp;
    SAType in_tune;             /* temporary struct for tune parameters */

    Opts *options;              /* used to restore command line options */
    MoveState *move_ptr;        /* used to restore moves */
    double *stats;              /* used to restore Lam stats */
    char *rand;                 /* used to restore ERand48 */
    double delta[2];            /* used to restore times */

    char *statefile = files->statefile;

    /* allocate memory for structures that will be returned (stats gets allo-  *
     * cated in StateRead(), since we need to know if we're tuning or not)     */

    options = ( Opts * ) malloc( sizeof( Opts ) );
    options->inname = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    options->outname = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    options->argv = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    options->derivfunc = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    options->solver = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    options->landscape_flag = 0;

    rand = ( char * ) calloc( INT_MAX, sizeof( char ) );
    stats = ( double * ) calloc( 31, sizeof( double ) );
    move_ptr = ( MoveState * ) malloc( sizeof( MoveState ) );
    //rand = ( unsigned short * ) calloc( 3, sizeof( unsigned short ) );

    StateRead( statefile, options, move_ptr, stats, rand, delta );
    /* restore options in fly_nsga2.c (and some in lsa.c) */

    RestoreOptions( options );

    /* initialize some Lam/Greening structures */
    p = state_ptr->tune.progname;       /* tune.progname contains program name */
    p = strcpy( p, version );

    p = strcpy( p, "The Other One" );   /* Grateful Dead tune, what else? */
    state_ptr->tune.tunefile = NULL;

    /* read data file and initialize different things for scoring and move     */
    /* generation (in zygotic.c, score.c and moves.c                           */
    inname = files->inputfile;
    infile = fopen( inname, "r" );
    if( !infile )
        file_error( "fly_nsga2" );

    if( debug ) {
        slogfile = fopen( strcat( inname, ".slog" ), "a" );
        if( !slogfile )
            file_error( "fly_*** error opening slog file" );
    }

    inp.zyg = InitZygote( infile, pd, pj, &inp, "input" );
    inp.sco = InitScoring( infile, method, &inp );
    inp.his = InitHistory( infile, &inp );      //It fills the polations vector
    inp.ext = InitExternalInputs( infile, &inp );
    inp.ste = InitStepsize( stepsize, accuracy, slogfile, inname );
    //Check if we need the following in Restore
    inp.twe = InitTweak( infile, NULL, inp.zyg.defs );
    inp.tra = Translate( &inp );
    in_tune = ReadSATune( infile );     /* read tune_parameter section */
    i_temp = InitMoves( infile, &inp ); /* set initial temperature and initialize */

    // initialize distribution stuff
    //printf("...ok!\nInitDist...");
    inp.dis = InitDistribution( infile );
    memcpy( distp, &( inp.dis ), sizeof( DistParms ) );

    fclose( infile );

    /* initialize Lam parameters (see sa.h for further detail) */

    state_ptr->tune.lambda = in_tune.lambda;
    state_ptr->tune.lambda_mem_length_u = in_tune.lambda_mem_length_u;
    state_ptr->tune.lambda_mem_length_v = in_tune.lambda_mem_length_v;
    state_ptr->tune.control = in_tune.control;
    state_ptr->tune.initial_moves = in_tune.initial_moves;
    state_ptr->tune.tau = in_tune.tau;
    state_ptr->tune.freeze_count = in_tune.freeze_count;
    state_ptr->tune.update_S_skip = in_tune.update_S_skip;
    state_ptr->tune.criterion = in_tune.criterion;
#ifdef MPI
    state_ptr->tune.mix_interval = in_tune.mix_interval;
#endif
    RestoreMoves( move_ptr );
    RestoreLamstats( stats );
    if( time_flag )
        RestoreTimes( delta );
    
    RestoreRand( rand );
    
    if( prolix_flag )
        RestoreProlix( files->prolixfile );
    free(rand);
}





/*** THE FINAL MOVE FUNCTION ***********************************************/

/**  FinalMove: determines the final energy and move count and then prints 
 *              those to wherever they need to be printed to; also should  
 *              do the cleaning up, i.e freeing stuff and such after a run 
 */
void
FinalMove( void ) {
    // int i;
#ifdef MPI
    int i;                      /* loop counter */
#endif

    AParms ap;

    double equil_var[2];        /* array for results of equilibration run */

    char *shell_cmd;            /* shell command to be executed by system */

#ifdef MPI
    int winner = 0;             /* id of the winning node */
    double minyet = DBL_MAX;    /* minimum score, used to find winner */

    double *final_e;            /* array of final energies of all nodes */

    final_e = ( double * ) calloc( nnodes, sizeof( double ) );
#endif

    /* get the answer and some additional info */

    ap = GetFinalInfo(  );      /* reads final energy and move count */

    if( equil )
        GetEquil( equil_var );  /* get the final equilibration results */

    
#ifdef MPI

    /* parallel code: find the node with the lowest score */

    for( i = 0; i < nnodes; i++ )       /* initialize the energy array */
        final_e[i] = 0;
    /* collect the final scores from all nodes */
    MPI_Allgather( &ap.stop_energy, 1, MPI_DOUBLE, final_e, 1, MPI_DOUBLE, MPI_COMM_WORLD );

    for( i = 0; i < nnodes; i++ ) {     /* evaluate the node with the lowest score */
        if( final_e[i] <= minyet ) {
            minyet = final_e[i];
            winner = i;
        }
    }

    /* write the answer */

    if( myid == winner ) {
#endif /* outfile different from infile? */
        if( strcmp( inname, outname ) ) {       /* -> create copy of infile */
            shell_cmd = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
            sprintf( shell_cmd, "cp -f %s %s", inname, outname );
            if( -1 == system( shell_cmd ) )
                error( "FinalMove: error creating output file %s", outname );
            free( shell_cmd );
        }

        /* all the funcs below write a section at its appropriate position in the  */
        /* data file; to achieve this, they create a temporary file which is then  */
        /* renamed to the final output file name                                   */
        /* printf("Residuals:\n");
           for (i=0; i<out->size_resid_arr; i++) {
           printf("%lg\t", out->residuals[i]);
           }
           printf("\n"); */

        WriteVersion( outname, version, argvsave );
        if( equil )
            WriteEquil( outname, equil_var );
        else {
            printf( "WARNING We are entering the writing of final output\n" );
            fflush( stdout );
            WriteParameters( outname, &( inp.zyg.parm ), "eqparms", precision, inp.zyg.defs );
            WriteAParameters( outname, ap );
            printf( "WARNING And we have done so successfully\n" );
            fflush( stdout );
        }

#ifdef MPI
    }
#endif

    /* clean up the state file and free memory */

    if( !equil && !nofile_flag )
#ifdef MPI
        if( !tuning )
#endif
            StateRm(  );
    FreeZygote(  );
    FreeHistory( inp.zyg.nalleles, inp.his );
    FreeExternalInputs( inp.zyg.nalleles, inp.ext );
}

/** WriteEquil: writes the equilibrate_variance section to the data file 
 *               right after the $equilibrate section                      
 */
void
WriteEquil( char *filename, double *equil_var ) {
    char *temp;                 /* temporary file name */
    char *record;               /* record to be read and written */
    char *record_ptr;           /* pointer used to remember record for 'free' */
    char *saverec;              /* used to save following section title */
    char *shell_cmd;            /* used by 'system' below */

    FILE *outfile;              /* name of output file */
    FILE *tmpfile;              /* name of temporary file */


    temp = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    record = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    saverec = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    shell_cmd = ( char * ) calloc( MAX_RECORD, sizeof( char ) );

    record_ptr = record;        /* this is to remember record for 'free' */

    /* open output and temporary file */

    outfile = fopen( filename, "r" );   /* open outfile for reading */
    if( !outfile )              /* sorry for the confusion! */
        error( "WriteEquil: error opening %s", filename );

    temp = strcpy( temp, "equilXXXXXX" );       /* required by mkstemp() */
    if( mkstemp( temp ) == -1 ) /* get unique name for temp file */
        error( "WriteEquil: error creating temporary file name" );

    tmpfile = fopen( temp, "w" );       /* ... and open it for writing */
    if( !tmpfile )
        error( "WriteEquil: error opening temporary file %s", temp );

    if( FindSection( outfile, "equilibrate_variance" ) ) {
        fclose( outfile );      /* erase section if already present */
        KillSection( filename, "equilibrate_variance" );
        outfile = fopen( filename, "r" );
    }
    rewind( outfile );

    /* the follwoing three loops look for the appropriate file position to     */
    /* write the equilibrate_variance section                                  */

    while( strncmp( record = fgets( record, MAX_RECORD, outfile ), "$equilibrate", 12 ) )
        fputs( record, tmpfile );
    fputs( record, tmpfile );

    while( strncmp( record = fgets( record, MAX_RECORD, outfile ), "$$", 2 ) )
        fputs( record, tmpfile );
    fputs( record, tmpfile );

    do {
        record = fgets( record, MAX_RECORD, outfile );
        if( !record )
            break;
    } while( strncmp( record, "$", 1 ) );

    fputs( "\n", tmpfile );

    if( record )
        saverec = strcpy( saverec, record );

    /* now we write the eqparms section into the tmpfile */

    PrintEquil( tmpfile, equil_var, "equilibrate_variance" );

    fprintf( tmpfile, "\n" );

    /* ... and then write all the rest, if there is any */

    if( record )
        fputs( saverec, tmpfile );

    while( ( record = fgets( record, MAX_RECORD, outfile ) ) )
        fputs( record, tmpfile );

    fclose( outfile );
    fclose( tmpfile );

    /* rename tmpfile into new file */

    sprintf( shell_cmd, "cp -f %s %s", temp, filename );

    if( -1 == system( shell_cmd ) )
        error( "WriteEquil: error renaming temp file %s", temp );

    if( remove( temp ) )
        warning( "WriteEquil: temp file %s could not be deleted", temp );

    /* clean up */

    free( record_ptr );
    free( saverec );
    free( temp );
    free( shell_cmd );
}

/** WriteTimes: writes the timing information to wherever it needs to be 
 *               written to at the end of a run                            
 */
void
WriteTimes( double *times ) {
    char *timefile;             /* name of the .times file */
    FILE *timeptr;              /* file pointer for .times file */

    /* create time file name by appending .times to input file name */
    timefile = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    timefile = strcpy( timefile, outname );
    timefile = strcat( timefile, ".times" );

    timeptr = fopen( timefile, "w" );
    if( !timeptr )
        file_error( "main" );

    PrintTimes( timeptr, times );       /* write times to .times file */

    fclose( timeptr );          /* clean up */
    free( timefile );
}

/** PrintTimes: writes two (parallel: three) times sections */
void
PrintTimes( FILE * fp, double *times ) {
    fprintf( fp, "wallclock: %.3f\n", times[0] );
    fprintf( fp, "user:      %.3f\n", times[1] );
}





/*** FUNCTIONS THAT COMMUNICATE WITH SAVESTATE.C ***************************/

/**  GetOptions: returns command line options to savestate.c 
 *               for the detailed meaning of all these options see Parse-  
 *               CommandLine() above); Opts struct defined in moves.h      
 */
Opts *
GetOptions( void ) {
    Opts *options;

    options = ( Opts * ) malloc( sizeof( Opts ) );

    options->inname = inname;
    options->outname = outname;
    options->argv = argvsave;

    options->derivfunc = ( char * ) calloc( MAX_RECORD, sizeof( char ) );

    if( pd == DvdtOrig )
        options->derivfunc = strcpy( options->derivfunc, "DvdtOrig" );
    else if( pd == Dvdt_sqrt )
        options->derivfunc = strcpy( options->derivfunc, "Dvdt_sqrt" );
    else
        error( "GetOptions: unknown derivative function" );

    options->solver = ( char * ) calloc( MAX_RECORD, sizeof( char ) );

    if( ps == Adams )
        options->solver = strcpy( options->solver, "Adams" );
    else if( ps == BaDe )
        options->solver = strcpy( options->solver, "BaDe" );
    else if( ps == BuSt )
        options->solver = strcpy( options->solver, "BuSt" );
    else if( ps == Euler )
        options->solver = strcpy( options->solver, "Euler" );
    else if( ps == Heun )
        options->solver = strcpy( options->solver, "Heun" );
    else if( ps == Milne )
        options->solver = strcpy( options->solver, "Milne" );
    else if( ps == Meuler )
        options->solver = strcpy( options->solver, "Meuler" );
    else if( ps == Rk2 )
        options->solver = strcpy( options->solver, "Rk2" );
    else if( ps == Rk4 )
        options->solver = strcpy( options->solver, "Rk4" );
    else if( ps == Rkck )
        options->solver = strcpy( options->solver, "Rkck" );
    else if( ps == Rkf )
        options->solver = strcpy( options->solver, "Rkf" );
    else if( ps == SoDe )
        options->solver = strcpy( options->solver, "SoDe" );
    else if( ps == Krylov )
        options->solver = strcpy( options->solver, "Krylov" );
    /*else if (ps == Band)
       options->solver = strcpy(options->solver, "Band");
     */
    else
        error( "GetOptions: unknown solver function" );

    options->stop_flag = stop_flag;
    options->prolix_flag = prolix_flag;
    options->landscape_flag = landscape_flag;
    options->log_flag = log_flag;
    options->time_flag = time_flag;
    options->state_write = state_write;
    options->print_freq = print_freq;
    options->captions = captions;
    options->olddivstyle = olddivstyle;
    options->precision = precision;
    options->stepsize = stepsize;
    options->quenchit = quenchit;
    options->equil = equil;

#ifdef MPI
    options->tuning = tuning;
    options->covar_index = covar_index;
    options->write_tune_stat = write_tune_stat;
    options->auto_stop_tune = auto_stop_tune;
#endif

    return options;
}

/**  RestoreOptions: restores the values of the command line opt variables 
 *                   from the Opts struct (used for restoring a run)       
 */
void
RestoreOptions( Opts * options ) {

    /* restore input/output file names and the full command line string; note  *
     * that the output file name needs to be communicated to lsa.c, since we   * 
     * need it there for setting up log and tuning file names                  */
    inname = options->inname;
    outname = options->outname;
    SetOutname( outname );
    argvsave = options->argv;

    /* all the other options */
    if( !strcmp( options->derivfunc, "DvdtOrig" ) )
        pd = DvdtOrig;
    else if( !strcmp( options->derivfunc, "Dvdt_sqrt" ) )
        pd = Dvdt_sqrt;
    else
        error( "RestoreOptions: unknown deriv function %s", options->derivfunc );

    if( !strcmp( options->solver, "Adams" ) )
        ps = Adams;
    else if( !strcmp( options->solver, "BaDe" ) )
        ps = BuSt;
    else if( !strcmp( options->solver, "BuSt" ) )
        ps = BuSt;
    else if( !strcmp( options->solver, "Euler" ) )
        ps = Euler;
    else if( !strcmp( options->solver, "Heun" ) )
        ps = Heun;
    else if( !strcmp( options->solver, "Milne" ) )
        ps = Milne;
    else if( !strcmp( options->solver, "Meuler" ) )
        ps = Meuler;
    else if( !strcmp( options->solver, "Rk2" ) )
        ps = Rk2;
    else if( !strcmp( options->solver, "Rk4" ) )
        ps = Rk4;
    else if( !strcmp( options->solver, "Rkck" ) )
        ps = Rkck;
    else if( !strcmp( options->solver, "Rkf" ) )
        ps = Rkf;
    else if( !strcmp( options->solver, "SoDe" ) )
        ps = SoDe;
    else if( !strcmp( options->solver, "Krylov" ) )
        ps = Krylov;
    /*
       else if (!strcmp(options->solver, "Band"))
       ps = Band;
     */
    else
        error( "RestoreOptions: unknown solver %s", options->solver );

    stop_flag = options->stop_flag;
    prolix_flag = options->prolix_flag;
    if( prolix_flag )
        SetProlix( prolix_flag, outname, 0 );
    landscape_flag = options->landscape_flag;
    if( landscape_flag )
        error( "RestoreOptions: cannot restore an equilibration (lanDscape) run" );

    log_flag = options->log_flag;
    time_flag = options->time_flag;

    state_write = options->state_write;
    print_freq = options->print_freq;
    captions = options->captions;

    olddivstyle = options->olddivstyle;
    precision = options->precision;
    stepsize = options->stepsize;
    quenchit = options->quenchit;
    equil = options->equil;
    if( equil )
        error( "RestoreOptions: cannot restore an equilibration run" );

#ifdef MPI
    tuning = options->tuning;
    if( tuning )
        error( "RestoreOptions: cannot restore a tuning run" );
    covar_index = options->covar_index;
    write_tune_stat = options->write_tune_stat;
    auto_stop_tune = options->auto_stop_tune;
#endif
    free( options->derivfunc );
    free( options->solver );
    free( options );
}
