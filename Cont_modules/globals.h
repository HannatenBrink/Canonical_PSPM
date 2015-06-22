/***
   NAME
     globals
   DESCRIPTION
     Header file with global definitions, such as error codes, return
     values and function prototypes.

   Last modification: AMdR - Mar 08, 2014
***/
#ifndef GLOBALS
#define GLOBALS
#include "ctype.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "stdint.h"
#include "stdarg.h"
#include "string.h"
#include "sys/param.h"
#include <sys/stat.h>
#include <signal.h>
#include <setjmp.h>
#include <unistd.h>
#include "sys/ioctl.h"

#if defined (__APPLE__)
#include <Accelerate/Accelerate.h>
#else
#include "cblas.h"
#endif



/*
 *===========================================================================
 *	Import the population and environment dimension settings
 *===========================================================================
 */
#ifndef		CONTINUATION
#define		CONTINUATION	0
#endif

#if defined(PROBLEMHEADER)
#include PROBLEMHEADER				// Include header file
#else
#if (!defined(CURVE) && !defined(COHORT) && !defined(IO) && !defined(ODE_DOPRI) && !defined(UTILS))
#warning No header file defined!
#endif
#endif


/*
 *===========================================================================
 *	Some numerical settings
 *===========================================================================
 */

#ifndef 	MAX_PNTDIM
#define         MAX_PNTDIM      	100
#endif
#ifndef 	MAX_OUTPUTDIM
#define         MAX_OUTPUTDIM      	100
#endif
#define		MAX_PNTDIM_SQ		(MAX_PNTDIM*MAX_PNTDIM)
#ifndef		MAXITER
#define 	MAXITER			40
#endif
#ifndef		MIN_STEP
#define 	MIN_STEP		1.0E-6
#endif
#ifndef		MINPARSTEPVAL
#define		MINPARSTEPVAL		1.0E-8
#endif
#ifndef		UPDATE_JAC
#define 	UPDATE_JAC		5	// Frequency of update
#endif
#ifndef		JACOBIAN_STEP
#define 	JACOBIAN_STEP		1.0E-4	// % change for jacobian
#endif
#ifndef		DYTOL
#define 	DYTOL			1.0E-7
#endif
#ifndef		RHSTOL
#define 	RHSTOL			1.0E-8
#endif
#ifndef		RHSMAX
#define 	RHSMAX			0.99
#endif
#ifndef		ALLOWNEGATIVE
#define		ALLOWNEGATIVE		0	// Allow negative solutions?
#endif
#ifndef		MOOREPENROSE			// Use Moore-Penrose (1) or
#define		MOOREPENROSE		0	// pseudo-arclength (0) continuation
#endif
#ifndef		PARSTEPSCALING
#define 	PARSTEPSCALING		1
#endif
#ifndef		CENTRALDIFF
#define		CENTRALDIFF		1	// Compute derivative by central difference
#else
#if (CENTRALDIFF != 1)
#undef 		CENTRALDIFF
#define		CENTRALDIFF		0	// Otherwise by forward difference
#endif
#endif
#ifndef		CVFFILE
#define		CVFFILE			1	// Use CVF file for input?
#endif

#define		STEP_DOUBLE		4
#define		STEP_HALF		(MAXITER/2)
#define		MAX_STEPREDUCE		2048
#define 	MAX_EXP			50.0

#define		MILLI			1.0E-3
#define		MICRO			1.0E-6
#define		NANO			1.0E-9
#define		PICO			1.0E-12
#define		FEMTO			1.0E-15
#define		ATTO			1.0E-18

#define BIFPARONE	((int)(parameter[14]+0.5))
#define BIFPARTWO	((int)(parameter[15]+0.5))
#define BIFPARTHREE ((int)(parameter[16]+0.5))
#define BIFPARFOUR  ((int)(parameter[17]+0.5))
/*
 *===========================================================================
 *	Numerical settings used in odesolve.c
 *===========================================================================
 */
#define 	INIT_H			0.001
#define 	LARGEST_STEP		2.0

#define 	ABS_ERR			FEMTO
#define 	ACCURACY		PICO
#define 	FUNCTOL			1.0E-8



/*
 *===========================================================================
 *	Macros used in curve.c
 *===========================================================================
 */
#define SUCCES			0
#define FAILURE			1

#define MEM_ALLOC_FAIL		100
#define MEM_COPY_FAIL		101
#define ZERO_DIVISION		102
#define SINGULARITY		103
#define NORM_OVERFLOW		104
#define NO_CONVERGENCE		105
#define ILLEGAL_INPUT		106
#define FAILED_EVALUATION	107
#define ILLEGAL_BOUNDS		108

#define MACH_PREC		1.0E-16
#define MACH_MAX		1.0E100

#define	FORWARD			0		// Computational methods for derivatives
#define CENTRAL			1
#define RICHARDSON		2


/*
 *===========================================================================
 *	Various other macro definitions
 *===========================================================================
 */
#define MAX_STR_LEN		1024
#if defined(DBL_MAX)				// Stub value for missing
#define MISSING_VALUE		(DBL_MAX)	// data point
#elif defined(MAXDOUBLE)
#define MISSING_VALUE		(MAXDOUBLE)
#else
#define MISSING_VALUE		(1.23456789e+307)
#endif

#define max(a, b)       	(((a) > (b)) ? (a) : (b))
#define min(a, b)       	(((a) < (b)) ? (a) : (b))
#define sign(a)       		(((a) < 0.0) ? (-1.0) : (1.0))
#define isodd(a)		(((a)%2))
#define iseven(a)		(!((a)%2))
#define issane(a)		((fpclassify(a)==FP_ZERO) || (fpclassify(a)==FP_NORMAL)) 

/*
 *===========================================================================
 *		Type definitions
 *===========================================================================
 */
typedef double	*cohort;



/*
 *===========================================================================
 *	Global variable definitions
 *===========================================================================
 */

#if POPULATION_NR
#ifndef MAX_COHORTS
#define MAX_COHORTS	1000
#endif
#ifndef COHORT_SIZE
#define   COHORT_SIZE    (1+I_STATE_DIM+I_CONST_DIM)	// Dimension of cohorts
#endif
#endif

#if (defined(CURVE) || defined(COHORT) || defined(IO) || defined(ODE_DOPRI) || defined(UTILS))
#undef	EXTERN
#define EXTERN	extern
#else
#undef	EXTERN
#define EXTERN
#endif

EXTERN double	pnt_scale[MAX_PNTDIM];
#if ENVIRON_DIM
EXTERN double	env[ENVIRON_DIM];
#endif
#if POPULATION_NR
EXTERN cohort	pop[POPULATION_NR][MAX_COHORTS];
EXTERN int	cohort_no[POPULATION_NR];
EXTERN double	logDead[POPULATION_NR];
#endif
#if PARAMETER_NR
EXTERN double	parameter[PARAMETER_NR];
#endif
EXTERN double	equal2zero;
EXTERN char	runname[MAX_STR_LEN];
EXTERN FILE	*errfile, *outfile;
EXTERN double	parstep, Maxparstep;
EXTERN jmp_buf  JumpBuffer;
EXTERN int	Stepchange;
EXTERN int	Stepreduce;
EXTERN int	NewStateComputed, Report2Terminal;
EXTERN double	RhsNorm;


/*
 *===========================================================================
 *	Single/double precision and BLAS/LAPACK function mapping
 *===========================================================================
 */

// BLAS Level 1 functions

#define ASUM		cblas_dasum
#define AXPY		cblas_daxpy
#define COPY		cblas_dcopy
#define DOT		cblas_ddot
#define NRM2		cblas_dnrm2
#define SCAL		cblas_dscal

// BLAS Level 2 functions

#define GEMV		cblas_dgemv
#define GER		cblas_dger

// Lapack functions
#define	GETRF		dgetrf_
#define	GETRS		dgetrs_
#define	GERFS		dgerfs_
#define	GECON		dgecon_
#define	GELS		dgels_

#if !defined (__APPLE__)
void 	GETRF(int *,  int *, double *, int *, int *, int *);
void 	GETRS(char *, int *, int *, double *, int *, int *, double *, int *, int *);
void 	GERFS(char *, int *, int *, double *, int *, double *, int *, int *,
	      double *, int *, double *, int *, double *, double *,
	      double *, int *, int *);
void 	GECON(char *, int *, double *, int *, double *, double *,
	      double *, int *, int *);
//void 	GEMV(const char *, int *, int *, double *, void *, int *, void *, int *, double *, void *, int *);
#endif


/*
 *===========================================================================
 *	Function prototyping: cohort.c
 *===========================================================================
 */

void 	SetCohortSize(const int);
void 	ReadIsfFile(const char *, const int);
void 	LabelState(int, const char *, ...);
#if POPULATION_NR
void 	WriteBinStateToFile(const char *,
			    void (*)(double e[ENVIRON_DIM],
				     cohort p[POPULATION_NR][MAX_COHORTS]));
#endif
void	WriteEsfFile(void);

/*
 *===========================================================================
 *	Function prototyping: curve.c
 *===========================================================================
 */

int			FindPoint(const int, double *, double *, double,
				  double, const int, int (*)(double *, double *));
int			FindPointGlobal(const int, double *, double *, double,
				        double, const int, int (*)(double *, double *));
int			TangentVec(const int, double *, double *,
				   int (*)(double *, double *));
int			Jacobian(const int pntdim, double *pnt, const int fncdim, double *jac,
				 int (*fnc)(double *, double *), int method);
int			CurveFuncDeriv(const int,  double *, double *, const int, double *, int (*)(double *, double *), int);
double			LPcondition(const int, double *, int (*)(double *, double *), int);
int             Evogradient(const int pntdim, const int evodim, double *pnt, int (*fnc)(double *, double *), int *evopars, double *k1);
double			SelectionGradient(const int pntdim, double *y, int (*fnc)(double *, double *), int varindex, int R0index);
double			SelectionGradient2(const int pntdim, double *y, int (*fnc)(double *, double *), int parindex, int R0index);
bool            Noevocheck(double *k1);


/*
 *===========================================================================
 *	Function prototyping: odesolve.c
 *===========================================================================
 */

int			odesolve(double *, int, double *, double,
				 void (*)(double, double *, double *),
				 double (*)(double, double *));

/*
 *===========================================================================
 *	Function prototyping: io.c
 *===========================================================================
 */

void			ReportMsg(int, const char *, ...);
void			ErrorMsg(const char *, const int, const char *, ...);
void	  		ErrorExit(const char *, const int, const char *, ...);
void			NumProcError(const char *, const int, const int);
void 			PrettyPrint(FILE *, const int, ...);
void 			PrettyPrintArray(FILE *, const int, double *);
void 			PrettyPrintHeader(FILE *, const int, ...);
char			*ReadDouble(double *, char *);
int			ReadDoublesOnLine(FILE *, double *, int, const int);
void	  		ReadCvfFile(const char *);
void 			ReportNote(const char *fmt, ...);
void	  		WriteReport(const char *, int, char **);
void  			SIGINThandler(int);
void 			installSIGINThandler(void);
void 			checkESCpressed(void);

/*
 *===========================================================================
 *	Function prototyping: utils.c
 *===========================================================================
 */

double 			anorm(int, int, double *);
inline double 		bexp(double);
inline int 		imax(int, int);
inline int 		imin(int, int);
inline double 		dmax(double, double);
inline double 		dmin(double, double);
int			SetScalesSimple(double *point, int basedim);
int			SetScales(double *point, int basedim, int dim);
double			determinant(double *, int);
int			ComputeEigenvector(const int, double *, double, double *);

/*
 *===========================================================================
 *	Function prototypeing: Problem-specific user file
 *===========================================================================
 */

int 			SetPntDim(void);
void			UserInit(int argc, char *argv[]);
int			Equation(double *argument, double *result);
int			DefineOutput(double *argument, double *output);
#if POPULATION_NR
void 			UpdateCohorts(double *env, cohort pop[POPULATION_NR][MAX_COHORTS]);
#endif

/*
 *===========================================================================
 *	Implementation of the generic main routine.
 *===========================================================================
 */

#if !(defined(CURVE) || defined(COHORT) || defined(IO) || defined(ODE_DOPRI) || defined(UTILS))

// Global variables for use in user-routines and main()

static double		point[MAX_PNTDIM];
#if (CVFFILE)
static int		restart = 0;
#endif

#if (CONTINUATION == 0)

int main(int argc, char **argv)

{
  int			i, j;
  int			maxind = 0, pntdim;
  double		result[MAX_PNTDIM];
  double		maxdif = 0.0;

  // Initialize some variables
  errfile = NULL;	outfile = NULL;
  NewStateComputed = 0; Report2Terminal = 1;

#if (CVFFILE)
  char			*ch;

  if (argc < (BASE_PNTDIM+2))
    ErrorExit(__FILE__, __LINE__, "Usage: %s [<Name of CVF, ISF or ESF file>] %s", argv[0], PROBLEMVARS);

  (void)strcpy(runname, argv[1]);

  // Strip possible CVF, ISF or ESF extensions
  if ((ch = strstr(runname, ".cvf"))) *ch = '\0';
  else if ((ch = strstr(runname, ".isf"))) *ch = '\0';
  else if ((ch = strstr(runname, ".esf"))) { *ch = '\0'; restart = 1; }
#else
  if (argc < (BASE_PNTDIM+1))
    ErrorExit(__FILE__, __LINE__, "Usage: %s %s", argv[0], PROBLEMVARS);
#endif

  // Map all command-line variables into the argument vector
  memset((void *)point,  0, MAX_PNTDIM*sizeof(double));
  for (i=0, j=CVFFILE+1; i<BASE_PNTDIM; i++, j++) point[i] = atof(argv[j]);

#if (POPULATION_NR)
  SetCohortSize(COHORT_SIZE);
#endif

#if (CVFFILE)
  ReadCvfFile(runname);
#endif

  UserInit(argc, argv);

#if (CVFFILE)
  if (!restart) WriteReport(runname, argc, argv);
#endif

  pntdim = SetPntDim();
  if (pntdim > MAX_PNTDIM) ErrorExit(__FILE__, __LINE__, "Dimension MAX_PNTDIM exceeded!");

  (void)SetScales(point, BASE_PNTDIM, pntdim);

  Equation(point, result);

  PrettyPrint(stdout, 1, point[0]);
  for (i=1; i<BASE_PNTDIM; i++)
    {
      PrettyPrint(stdout, 2, point[i], result[i-1]);
      if (fabs(result[i-1]) > maxdif)
	{
	  maxdif = fabs(result[i-1]);
	  maxind = i-1;
	}
    }

  fprintf(stderr, "\n\n\n");
  for (i=0; i<BASE_PNTDIM; i++)
    fprintf(stderr, "Scale variable %2d: %.1E\n", i, pnt_scale[i]);

  fprintf(stderr, "\n\nMaximum difference in component %d: %f\n",
	  maxind+1, result[maxind]);

#if POPULATION_NR
  Equation(point, NULL);
  WriteBinStateToFile(runname, UpdateCohorts);
#endif

  return 0;
}




/*==============================================================================*/

#elif (CANONICAL == 1)

static int		CurveEnd = 0;

static void Usage(char *progname)

{
#if (CVFFILE)
  fprintf(stderr, "\n\nUsage: %s [-v <0, 1, 2>] <Name of CVF, ISF or ESF file> %s %s\n\n",
      progname, PROBLEMVARS, "<step> [<min. step>] <max. step> <min. par.> <max. par.>");
#else
  fprintf(stderr, "\n\nUsage: %s [-v <0, 1, 2>] %s %s\n\n",
      progname, PROBLEMVARS, "<step> [<min. step>] <max. step> <min. par.> <max. par.>");
#endif
  exit(1);

  return;
}

/*==============================================================================*/
int main(int argc, char **argv)

{
    register int		i, j, outmax;
    double		oldpoint[MAX_PNTDIM], oldpointselec[MAX_PNTDIM], oldpar;
    double		tanvec[MAX_PNTDIM], selvec[MAX_PNTDIM];
    int			pntnr = 0, pntdim, oldpntdim, nargs;
    int			cycles, last = 1, retval, setpoint, fp;
    char			fname[MAX_STR_LEN];
    double		Output[MAX_OUTPUTDIM];
    double		minbound, maxbound;
    double		maxparstep, Minparstep;
    int			my_argc;
    char		**argpnt1 = NULL, **argpnt2 = NULL, **my_argv = NULL;
    double      h;
    int         evopar[EVODIM];
    double      k1[EVODIM];

    
    
    // Initialize some variables
    errfile = NULL;	outfile = NULL;
    parstep = 1.0;	Maxparstep = -0.001; Minparstep = MINPARSTEPVAL;
    Stepchange = 0;	Stepreduce = 1;
    NewStateComputed = 0; Report2Terminal = 1;
    
    my_argv = (char **)calloc((size_t)argc, sizeof(char *));
    
    /*
     * Process all intrinsic command line parameters. Possibilities:
     *
     *	-v 		 	: Level of reporting information, written to
     *			   	  the terminal
     *
     *	-?   | --help		: Print usage message
     */
    argpnt1 = argv;
    argpnt2 = my_argv;
    my_argc = 0;
    while (*argpnt1)
    {
        if (!strcmp(*argpnt1, "-?") || !strcmp(*argpnt1, "--help"))
        {
            Usage(argv[0]);
        }
        else if (!strcmp(*argpnt1, "-v"))
        {
            argpnt1++;
            if (!*argpnt1)
            {
                fprintf(stderr, "\nNo reporting level specified!\n");
                Usage(argv[0]);
            }
            switch (atoi(*argpnt1))
            {
                case 0:
                case 1:
                case 2:
                    Report2Terminal = atoi(*argpnt1);
                    break;
                default:
                    fprintf(stderr, "\nWrong reporting level specifier: %s\n", *argpnt1);
                    Usage(argv[0]);
                    break;
            }
        }
        else if ((!strncmp(*argpnt1, "--", 2)))
        {
            fprintf(stderr, "\nUnknown command line option: %s\n", *argpnt1);
            Usage(argv[0]);
        }
        else if ((!strncmp(*argpnt1, "-", 1)) && isalpha(*(*argpnt1+1)))
        {
            fprintf(stderr, "\nUnknown command line option: %s\n", *argpnt1);
            Usage(argv[0]);
        }
        else
        {
            *argpnt2 = *argpnt1;
            my_argc++;
            argpnt2++;
        }
        argpnt1++;
    }
    
#if (CVFFILE)
    char			*ch;
    
    nargs = BASE_PNTDIM+6;
    if ((my_argc != nargs) && (my_argc != (nargs+1))) Usage(argv[0]);
    (void)strcpy(runname, my_argv[1]);
    
    // Strip possible CVF, ISF or ESF extensions
    if ((ch = strstr(runname, ".cvf"))) *ch = '\0';
    else if ((ch = strstr(runname, ".isf"))) *ch = '\0';
    else if ((ch = strstr(runname, ".esf"))) { *ch = '\0'; restart = 1; }
#else
    nargs = BASE_PNTDIM+5;
    if ((my_argc != nargs) && (my_argc != (nargs+1))) Usage(argv[0]);
    (void)strcpy(runname, my_argv[0]);
#endif
    
    // Map all command-line variables into the argument vector
    memset((void *)point,  0, MAX_PNTDIM*sizeof(double)); //copieert allemaal nullen naar point
    for (i=0, j=CVFFILE+1; i<BASE_PNTDIM; i++, j++) point[i] = atof(my_argv[j]); //init point
    
    parstep    = atof(my_argv[j++]);
    if (my_argc == nargs)
    {
        Maxparstep = atof(my_argv[j++]);
        if (Maxparstep > 0.0)
            Minparstep = min(MINPARSTEPVAL, MILLI*Maxparstep);
        else
            Minparstep = MINPARSTEPVAL;
    }
    else if (my_argc == (nargs+1))
    {
        Minparstep = atof(my_argv[j++]);
        Maxparstep = atof(my_argv[j++]);
    }
    
    minbound   = atof(my_argv[j++]);
    maxbound   = atof(my_argv[j++]);
    
#if (POPULATION_NR)
    SetCohortSize(COHORT_SIZE);
#endif
    
#if (CVFFILE)
    ReadCvfFile(runname);
#endif
    
    UserInit(my_argc, my_argv);
    
#if (CVFFILE)
    if (!restart) WriteReport(runname, my_argc, my_argv);
#endif
    
    pntdim = SetPntDim();
    if (pntdim > MAX_PNTDIM) ErrorExit(__FILE__, __LINE__, "Dimension MAX_PNTDIM exceeded!");
    
    (void)SetScales(point, BASE_PNTDIM, pntdim);
    
    (void)strcpy(fname, runname);
    (void)strcat(fname, ".err");
    errfile=fopen(fname, "a");
    
    (void)strcpy(fname, runname);
    (void)strcat(fname, ".out");
    outfile=fopen(fname, "a");
    
    // Use Ctrl-C to change step size
    installSIGINThandler();
    
    // Find the first point on the curve
    memset((void *)tanvec, 0, MAX_PNTDIM*sizeof(double));
    tanvec[0] = 1.0;
    
      retval = FindPoint(pntdim, point, tanvec, DYTOL, RHSTOL, 2*MAXITER, Equation);
      if (retval != SUCCES)
        {
          NumProcError(__FILE__, __LINE__, retval);
          ErrorMsg(__FILE__, __LINE__, "Failure to locate first point");
          return 1;
        }
   
    //INITIALIZATION OF EVOPAR //Niet echt mooi zo, maar werkt wel...
    
    evopar[0] = BIFPARONE;
#if EVODIM == 2
    evopar[1] = BIFPARTWO;
#endif
#if EVODIM == 3
    evopar[1] = BIFPARTWO;
    evopar[2] = BIFPARTHREE;
#endif
#if EVODIM == 4
    evopar[1] = BIFPARTWO;
    evopar[2] = BIFPARTHREE;
    evopar[3] = BIFPARFOUR;
#endif
    
 
    ///=====NIEUWE CODE CANONICAL=====///
    

    
    // Continue the curve
    while (1)
    {
        
        cycles = 0;
        
        //**** The following is new and important!!!
        memset((void *)tanvec, 0, pntdim*sizeof(double));
        tanvec[0] = 1.0;
        
        retval = FindPoint(pntdim, point, tanvec, DYTOL, RHSTOL, MAXITER, Equation);
        

        if (retval != SUCCES)
        {
            NumProcError(__FILE__, __LINE__, retval);
            ErrorMsg(__FILE__, __LINE__, "Failed to locate a solution point");
            break;
        }
        
        // Scale the point vector anew if necessary and redo the current point
        if ((retval = SetScales(point, BASE_PNTDIM, pntdim)))
            ReportMsg(1, "\n\nVariable %d rescaled!\n", retval);
        // Otherwise generate output and predict new point on the curve
        else
            if (Report2Terminal)
            {
                ReportMsg(1, "Point found:\t");
                for (i=0; i<BASE_PNTDIM; i++)
                ReportMsg(1, "%10.5E\t", point[i]*pnt_scale[i]);
     }
        
        
        
        
        for (i=0; (pntnr > 20) && (i<pntdim); i++)
        {
            

            
            CurveEnd = CurveEnd || ((point[0]*pnt_scale[0]) < minbound) || ((point[0]*pnt_scale[0]) > maxbound) || Noevocheck(k1);
#if (!ALLOWNEGATIVE)
            CurveEnd = CurveEnd || (point[i] < -DYTOL);
#endif
            
        }
        
        if (outfile)
        {
            outmax = DefineOutput(point, Output);
            if (outmax)
            {
                Output[outmax++] = RhsNorm;
                PrettyPrintArray(outfile, outmax, Output);
                fflush(outfile);
            }
        }
        if (CurveEnd) return 0;
        
        pntnr++;
        
        
      
        Evogradient(pntdim, EVODIM, point, Equation, evopar, k1); //Calculate selectiongradient
        
       point[0] = point[0]+DELTAS*k1[0];  //update bifpar
        
        for (i=1;i<EVODIM;i++)
        {
            parameter[evopar[i]] = k1[i]*DELTAS+parameter[evopar[i]]; //update all other parameters
        }

        if (errfile) fflush(errfile);
        
    }
    
    
    return 0;
}
/*==========================================================================*/
#else

static int		CurveEnd = 0;

static void Usage(char *progname)

{
#if (CVFFILE)
    fprintf(stderr, "\n\nUsage: %s [-v <0, 1, 2>] <Name of CVF, ISF or ESF file> %s %s\n\n",
            progname, PROBLEMVARS, "<step> [<min. step>] <max. step> <min. par.> <max. par.>");
#else
    fprintf(stderr, "\n\nUsage: %s [-v <0, 1, 2>] %s %s\n\n",
            progname, PROBLEMVARS, "<step> [<min. step>] <max. step> <min. par.> <max. par.>");
#endif
    exit(1);
    
    return;
}

/*==========================================================================*/
int main(int argc, char **argv)

{
  register int		i, j, outmax;
  double		oldpoint[MAX_PNTDIM];
  double		tanvec[MAX_PNTDIM];
  int			pntnr = 0, pntdim, oldpntdim, nargs;
  int			cycles, last = 1, retval, setpoint;
  char			fname[MAX_STR_LEN];
  double		Output[MAX_OUTPUTDIM];
  double		minbound, maxbound;
  double		maxparstep, Minparstep;
  int			my_argc;
  char			**argpnt1 = NULL, **argpnt2 = NULL, **my_argv = NULL;

  // Initialize some variables
  errfile = NULL;	outfile = NULL;
  parstep = 1.0;	Maxparstep = -0.001; Minparstep = MINPARSTEPVAL;
  Stepchange = 0;	Stepreduce = 1;
  NewStateComputed = 0; Report2Terminal = 1;

  my_argv = (char **)calloc((size_t)argc, sizeof(char *));

  /*
   * Process all intrinsic command line parameters. Possibilities:
   *
   *	-v 		 	: Level of reporting information, written to
   *			   	  the terminal
   *
   *	-?   | --help		: Print usage message
   */
  argpnt1 = argv;
  argpnt2 = my_argv;
  my_argc = 0;
  while (*argpnt1)
    {
      if (!strcmp(*argpnt1, "-?") || !strcmp(*argpnt1, "--help"))
	{
	  Usage(argv[0]);
	}
      else if (!strcmp(*argpnt1, "-v"))
	{
	  argpnt1++;
	  if (!*argpnt1)
	    {
	      fprintf(stderr, "\nNo reporting level specified!\n");
	      Usage(argv[0]);
	    }
	  switch (atoi(*argpnt1))
	    {
	    case 0:
	    case 1:
	    case 2:
	      Report2Terminal = atoi(*argpnt1);
	      break;
	    default:
	      fprintf(stderr, "\nWrong reporting level specifier: %s\n", *argpnt1);
	      Usage(argv[0]);
	      break;
	    }
	}
      else if ((!strncmp(*argpnt1, "--", 2)))
	{
	  fprintf(stderr, "\nUnknown command line option: %s\n", *argpnt1);
	  Usage(argv[0]);
	}
      else if ((!strncmp(*argpnt1, "-", 1)) && isalpha(*(*argpnt1+1)))
	{
	  fprintf(stderr, "\nUnknown command line option: %s\n", *argpnt1);
	  Usage(argv[0]);
	}
      else
	{
	  *argpnt2 = *argpnt1;
	  my_argc++;
	  argpnt2++;
	}
      argpnt1++;
    }

#if (CVFFILE)
  char			*ch;

  nargs = BASE_PNTDIM+6;
  if ((my_argc != nargs) && (my_argc != (nargs+1))) Usage(argv[0]);
  (void)strcpy(runname, my_argv[1]);

  // Strip possible CVF, ISF or ESF extensions
  if ((ch = strstr(runname, ".cvf"))) *ch = '\0';
  else if ((ch = strstr(runname, ".isf"))) *ch = '\0';
  else if ((ch = strstr(runname, ".esf"))) { *ch = '\0'; restart = 1; }
#else
  nargs = BASE_PNTDIM+5;
  if ((my_argc != nargs) && (my_argc != (nargs+1))) Usage(argv[0]);
  (void)strcpy(runname, my_argv[0]);
#endif

  // Map all command-line variables into the argument vector
  memset((void *)point,  0, MAX_PNTDIM*sizeof(double));
  for (i=0, j=CVFFILE+1; i<BASE_PNTDIM; i++, j++) point[i] = atof(my_argv[j]);

  parstep    = atof(my_argv[j++]);
  if (my_argc == nargs)
    {
      Maxparstep = atof(my_argv[j++]);
      if (Maxparstep > 0.0)
	Minparstep = min(MINPARSTEPVAL, MILLI*Maxparstep);
      else
	Minparstep = MINPARSTEPVAL;
    }
  else if (my_argc == (nargs+1))
    {
      Minparstep = atof(my_argv[j++]);
      Maxparstep = atof(my_argv[j++]);
    }

  minbound   = atof(my_argv[j++]);
  maxbound   = atof(my_argv[j++]);

#if (POPULATION_NR)
  SetCohortSize(COHORT_SIZE);
#endif

#if (CVFFILE)
  ReadCvfFile(runname);
#endif

  UserInit(my_argc, my_argv);

#if (CVFFILE)
  if (!restart) WriteReport(runname, my_argc, my_argv);
#endif

  pntdim = SetPntDim();
  if (pntdim > MAX_PNTDIM) ErrorExit(__FILE__, __LINE__, "Dimension MAX_PNTDIM exceeded!");

  (void)SetScales(point, BASE_PNTDIM, pntdim);

  (void)strcpy(fname, runname);
  (void)strcat(fname, ".err");
  errfile=fopen(fname, "a");

  (void)strcpy(fname, runname);
  (void)strcat(fname, ".out");
  outfile=fopen(fname, "a");

  // Use Ctrl-C to change step size
  installSIGINThandler();

  // Find the first point on the curve
  memset((void *)tanvec, 0, MAX_PNTDIM*sizeof(double));
  tanvec[0] = 1.0;

//  retval = FindPoint(pntdim, point, tanvec, DYTOL, RHSTOL, 2*MAXITER, Equation);
//  if (retval != SUCCES)
//    {
//      NumProcError(__FILE__, __LINE__, retval);
//      ErrorMsg(__FILE__, __LINE__, "Failure to locate first point");
//      return 1;
//    }

  // Continue the curve
  while (1)
    {
      // Use Ctrl-C to change step size and restart here
      setpoint = setjmp(JumpBuffer);
      if (setpoint)
	{
	  // Here we return after changing the step size. We should only get
	  // here if the following operations are valid!
	  Stepreduce = 1;
	  COPY(pntdim, oldpoint, 1, point, 1);
	  AXPY(pntdim, (parstep/pnt_scale[0]), tanvec, 1, point, 1);

	  ReportMsg(1, "\n\nPrediction :\t");
	  for (i=0; i<BASE_PNTDIM; i++)
	    ReportMsg(1, "%10.5E\t", point[i]*pnt_scale[i]);
	  ReportMsg(1, "\n");
	}

      // Compute fixed point with new varied parameter
      cycles = 0;
      while(Stepreduce <= MAX_STEPREDUCE)
	{
	  retval = FindPoint(pntdim, point, tanvec, DYTOL, RHSTOL, MAXITER, Equation);

	  // If not successful with standard Newton method retry with error-based damped Newton algorithm
	  if (retval != SUCCES)
	    retval = FindPointGlobal(pntdim, point, tanvec, DYTOL, RHSTOL, MAXITER, Equation);

	  if (retval == SUCCES)
	    {
	      if (Report2Terminal)
		{
		  ReportMsg(1, "Point found:\t");
		  for (i=0; i<BASE_PNTDIM; i++)
		    ReportMsg(1, "%10.5E\t", point[i]*pnt_scale[i]);
		  ReportMsg(1, "\n");
		}
	      else
		{
		  for (i=0; i<BASE_PNTDIM; i++)
		    printf("%15.10E\t", point[i]*pnt_scale[i]);
		  printf("\n");
		}
	      break;
	    }

	  NumProcError(__FILE__, __LINE__, retval);

	  // If unsuccesfull while repeating point exit
	  if (!Stepchange) ErrorExit(__FILE__, __LINE__, "Failed to locate solution point!");

	  // Generate prediction of solution point with smaller step size
	  cycles++;
	  Stepreduce *= 8;
	  COPY(pntdim, oldpoint, 1, point, 1);

	  AXPY(pntdim, (parstep/(Stepreduce*pnt_scale[0])), tanvec, 1, point, 1);

	  ReportMsg(1, "\n\nPrediction :\t");
	  for (i=0; i<BASE_PNTDIM; i++)
	    ReportMsg(1, "%10.5E\t", point[i]*pnt_scale[i]);
	  ReportMsg(1, "\n");
	}

      // If unsuccesfull exit
      if (retval != SUCCES) ErrorExit(__FILE__, __LINE__, NULL);

      // Increase step when both this and previous point were located with
      // the current step size
      if ((last && !cycles) && (Stepreduce > 1)) Stepreduce /= 2;
      last = (!cycles);

      oldpntdim = pntdim;
      pntdim = SetPntDim();
      if (pntdim > MAX_PNTDIM) ErrorExit(__FILE__, __LINE__, "Dimension MAX_PNTDIM exceeded!");

      if (oldpntdim != pntdim)
	{
	  (void)ReportMsg(1, "\n\nProblem dimension changed!\n");

	  // For added problem dimensions extrapolate both the point and the
	  // tangent vector elements
	  for (i=oldpntdim; i<pntdim; i++)
	    {
	      point[i]  = max(0.0, point[oldpntdim-1] + (i-oldpntdim+1)*(point[oldpntdim-1] - point[oldpntdim-2]));
	      tanvec[i] = tanvec[oldpntdim-1] + (i-oldpntdim+1)*(tanvec[oldpntdim-1] - tanvec[oldpntdim-2]);
	    }

	  Stepchange = 0;
	}
      // Scale the point vector anew if necessary and redo the current point
      else if ((retval = SetScales(point, BASE_PNTDIM, pntdim)))
	{
	  ReportMsg(1, "\n\nVariable %d rescaled!\n", retval);
	  // I think I have to compute the tangent vector anew here, because of the change in scaling
	  retval   = TangentVec(pntdim, point, tanvec, Equation);

	  Stepchange = 0;
	}
      // Otherwise generate output and predict new point on the curve
      else
	{
	  // The following computation of the tangent vector is not necessary in case of MoorePenrose continuation
#if (MOOREPENROSE != 1)
	  retval   = TangentVec(pntdim, point, tanvec, Equation);
#endif
	  // Signal curve stop if one of the components has become negative or parameter is out of bounds
	  // and this is not the curve beginning
	  for (i=0; (pntnr > 20) && (i<pntdim); i++)
	    {
	      CurveEnd = CurveEnd || ((point[0]*pnt_scale[0]) < minbound) || ((point[0]*pnt_scale[0]) > maxbound);
#if (!ALLOWNEGATIVE)
	      CurveEnd = CurveEnd || (point[i] < -DYTOL);
#endif
	    }

	  // Generate output: invoked after setting CurveEnd to allow for output of last solution point on branch
	  if (outfile)
	    {
	      outmax = DefineOutput(point, Output);
	      if (outmax)
		{
		  Output[outmax++] = RhsNorm;
		  PrettyPrintArray(outfile, outmax, Output);
		  fflush(outfile);
		}
	    }

	  // Exit at end of curve after generation of last output point
	  if (CurveEnd) return 0;

	  pntnr++;

	  // Store info on current solution point
	  COPY(pntdim, point, 1, oldpoint, 1);
	  Stepchange = 1;

	  // Determine appropriate stepsize
	  ReportMsg(2, "Tangent vector in component   0: %.8G\n", tanvec[0]);
	  ReportMsg(2, "Targeted  step in component   0: %.8G\n", parstep*tanvec[0]/Stepreduce);
#if (PARSTEPSCALING == 1)
	  maxparstep = Maxparstep*pnt_scale[0];
	  if (Maxparstep > 0.0)
	    {
	      if ((Stepreduce == 1) && (fabs(parstep*tanvec[0]) > Maxparstep))
		parstep = sign(parstep)*fabs(Maxparstep/tanvec[0]);
	    }
	  else
	    {
	      if ((Stepreduce == 1) && (fabs(parstep*tanvec[0]) > fabs(maxparstep*point[0])))
		parstep = sign(parstep)*fabs(maxparstep*point[0]/tanvec[0]);
	    }
	  if (fabs(parstep*tanvec[0]) < Minparstep)
	    parstep = sign(parstep)*fabs(Minparstep/tanvec[0]);
	  AXPY(pntdim, (parstep/(Stepreduce*pnt_scale[0])), tanvec, 1, point, 1);
#else
	  if ((Stepreduce == 1) && (fabs(parstep) > Maxparstep))
	    parstep = sign(parstep)*fabs(Maxparstep);
	  if (fabs(parstep) < Minparstep)
	    parstep = sign(parstep)*fabs(Minparstep);
	  AXPY(pntdim, (parstep/Stepreduce), tanvec, 1, point, 1);
#endif
	  ReportMsg(2, "Realized  step in component   0: %.8G\n", parstep*tanvec[0]/Stepreduce);

	  // Prediction
	  ReportMsg(1, "\n\nPrediction :\t");
	  for (i=0; i<BASE_PNTDIM; i++)
	    ReportMsg(1, "%10.5E\t", point[i]*pnt_scale[i]);
	  ReportMsg(1, "\n");
	}

      if (errfile) fflush(errfile);
    }

  return 0;
}

#endif
#endif // !(defined(CURVE) || defined(COHORT) || defined(IO) || defined(ODE_DOPRI) || defined(UTILS))

/*==============================================================================*/

inline double	bexp(double x)

{
    double		pw = x;
    
    pw = max(pw, -MAX_EXP);
    pw = min(pw,  MAX_EXP);
    
    return exp(pw);
}


/*==============================================================================*/

inline int	  imin(int a, int b)

{
    return (a < b) ? a : b;
}


/*==============================================================================*/

inline int	  imax(int a, int b)

{
    return (a > b) ? a : b;
}


/*==============================================================================*/

inline double dmax(double a, double b)
{
    double		dmaxa = a, dmaxb = b;
    return (dmaxa > dmaxb) ? dmaxa:dmaxb;
}


/*==============================================================================*/

inline double dmin(double a, double b)
{
    double		dmina = a, dminb = b;
    return (dmina < dminb) ? dmina:dminb;
}


/*==============================================================================*/

/*==============================================================================*/

#endif // GLOBALS
