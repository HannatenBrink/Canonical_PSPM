/*============================================================================
 *
 * Main interface routine nleq_res() originally defined in nleq_res.c
 *
 *============================================================================
 *
    NLEQ_RES - RESidual-based Damped Newton algorithm
    
 *  Written by        L. Weimann 
 *  Purpose           Solution of systems of nonlinear equations
 *  Method            RESidual-based Damped Newton algorithm
                      (see reference below)
 *  Category          F2a. - Systems of nonlinear equations
 *  Keywords          Nonlinear equations, Newton methods
 *  Version           1.1.1
 *  Revision          May 2006
 *  Latest Change     June 2006
 *  Library           NewtonLib
 *  Code              C, Double Precision
 *  Environment       Standard C environment on PC's,
                      workstations and hosts.
 *  Copyright     (c) Konrad-Zuse-Zentrum fuer
                      Informationstechnik Berlin (ZIB)
                      Takustrasse 7, D-14195 Berlin-Dahlem
                      phone : + 49/30/84185-0
                      fax   : + 49/30/84185-125
 *  Contact           Lutz Weimann
                      ZIB, Division Scientific Computing, 
                           Department Numerical Analysis and Modelling
                      phone : + 49/30/84185-185
                      fax   : + 49/30/84185-107
                      e-mail: weimann@zib.de
 
 *    References:
 
      /1/ P. Deuflhard:
          Newton Methods for Nonlinear Problems. -
          Affine Invariance and Adaptive Algorithms.
          Series Computational Mathematics 35, Springer (2004)
 
   ---------------------------------------------------------------
 
 * Licence
     You may use or modify this code for your own non commercial
     purposes for an unlimited time. 
     In any case you should not deliver this code without a special 
     permission of ZIB.
     In case you intend to use the code commercially, we oblige you
     to sign an according licence agreement with ZIB.
 
 * Warranty 
     This code has been tested up to a certain level. Defects and
     weaknesses, which may be included in the code, do not establish
     any warranties by ZIB. ZIB does not take over any liabilities
     which may follow from acquisition or application of this code.
 
 * Software status 
     This code is under care of ZIB and belongs to ZIB software class 2.
 
      ------------------------------------------------------------
 
 *    Parameters description
      ======================
 
      The calling interface looks as follows:

      extern void nleq_res(struct NLEQ_FUN funs, int n, double *x,
                           struct NLEQ_OPT *opt, 
                           struct NLEQ_INFO *info)

      The structures used within the parameter list are defined
      as follows:
      ---
      struct NLEQ_FUN
      {
        NLEQ_FFUN *fun;
        NLEQ_JFUN *jac;
      };
      
      where the types used within this structure are defined by
      typedef void NLEQ_FFUN(int*,double*,double*,int*);
      and
      typedef void NLEQ_JFUN(int*,int*,double*,double*,int*);
      ---
      struct NLEQ_OPT
      {
         double tol;
         int maxiter;
         LOGICAL restricted, nleqcalled; 
         PRINT_LEVEL errorlevel, monitorlevel, datalevel;
         FILE *errorfile, *monitorfile, *datafile;
         PROBLEM_TYPE nonlin;
         double *scale;
      };
      
      where the types used within this structure are defined by
      typedef enum {None=0, Minimum=1, Verbose=2, Debug=3} PRINT_LEVEL;
      typedef enum {False=0, True=1} LOGICAL ;
      typedef enum {Mildly_Nonlinear=2, Highly_Nonlinear=3,
                    Extremely_Nonlinear=4} PROBLEM_TYPE ;
      ---
      struct NLEQ_INFO
      {
         double precision, normdx;
         double *fx;
         int iter, rcode, subcode, nofunevals, nojacevals;
      };
      ---
      
      A detailed description of the parameters follows: 
      
      struct NLEQ_FUN funs :
      
      The field funs.fun must contain a pointer to the problem function fun -
      The required parameters interface of fun is described in detail below
      
      The field funs.jac must either contain a pointer to the Jacobian function
      jac or a NULL pointer. If a NULL pointer is supplied, then the Jacobian
      will be approximately computed by an internal function of the nleq_res
      package.
      
      int n :
      The number of equations and unknown variables of the nonlinear system.
      
      double *x :
      A pointer to an array of double values of size n.
      The array must contain on input an initial guess of the problems solution,
      which is used as the start-vector of the damped Newton iteration.
      On output, the pointed array contains an approximate solution vector x*,
      which fits the small residual condition
      || funs.fun(x*) || <= opt->tol,
      where ||...|| is a scaled Euclidian norm.
      
      struct NLEQ_OPT *opt:
      A pointer to an options structure. The pointed fields of the structure
      contain input options to nleq_res.
      
      opt->tol is of type double and must contain the residuum threshold
      which funs.fun(x*) must fit for the solution vector x*.
      
      opt->maxiter is of type int and must contain the maximum number of allowed
      iterations. if a zero or negative value is supplied, then the maximum
      iteration count will be set to 50.
      
      opt->nonlin is of type PROBLEM_TYPE and must classify the problem to
      be solved. The following classifications may be used:
      Mildly_Nonlinear: The problem is considered to be mildly nonlinear and
                        nleq_res starts up with dampingfactor=1.
      Highly_Nonlinear: The problem is considered to be highly nonlinear and
                        nleq_res starts up with dampingfactor=1.0e-4.
      Extremely_Nonlinear: The problem is considered to be extremely nonlinear
                        and nleq_res starts up with dampingfactor=1.0e-6.
                        Moreover, opt->restricted is set automatically to True.

      opt->restricted is of type LOGICAL.
      If set to True, then the restricted monotonicity test will be applied for
      determination whether the next iterate (and the associate damping factor
      lambda) will be accepted. This means, with
      theta = ||F(x(k+1))|| / ||F(x(k))||, the condition
      theta <= 1.0 - lambda/4 must be fit
      If set to False, then the standard monotonicity test will be applied, i.e.
      the following condition must be fit:
      theta < 1.0.
      
      opt->nleqcalled is of type LOGICAL and only used internally. This field
      should always be set to False.
      
      opt->errorlevel is of type PRINT_LEVEL. If it is set to level None,
      then no error message will be printed if an error condition occurs.
      If it is set to level Minimum or any higher level, then error messages
      will be printed, if appropriate.
      
      opt->monitorlevel is of type PRINT_LEVEL. If it is set to level None,
      then no monitor output will be printed.
      If it is set to level Minimum, a few infomation will be printed.
      If set to level Verbose, then some infomation about each Newton iteration
      step, fitting into a single line, will be printed. The higher level Debug
      is reserved for future additional information output.
      
      opt->datalevel is of type PRINT_LEVEL. If it is set to level None,
      then no data output will be printed.
      If it is set to level Minimum, then the values of the initial iteration
      vector x and the final vector x will be printed.
      If set to level Verbose, then the iteration vector x will be printed for
      each Newton step. The higher level Debug is reserved for future additional
      information output.
      
      opt->errorfile is of type pointer to FILE, as defined by the <stdio.h>
      header file. It must be set either to a NULL pointer, stderr, stdout,
      or to another file pointer which has been initialized by a fopen call.
      If it is set to NULL, opt->errorfile will be set to stdout. The error 
      messages will be printed to opt->errorfile.
      
      opt->monitorfile is of type pointer to FILE, as defined by the <stdio.h>
      header file. It must be set either to a NULL pointer, stderr, stdout,
      or to another file pointer which has been initialized by a fopen call.
      If it is set to NULL, opt->monitorfile will be set to stdout. The monitor 
      output will be printed to opt->monitorfile.
      
      opt->datafile is of type pointer to FILE, as defined by the <stdio.h>
      header file. It must be set either to a NULL pointer, stderr, stdout,
      or to another file pointer which has been initialized by a fopen call.
      If it is set to NULL, a file named "nleq_res.data" will be opened by a
      fopen call and opt->datafile will be set to the filepointer which the
      fopen returns. The data output will be printed to opt->datafile.
      
      opt->iterfile is of type pointer to FILE, as defined by the <stdio.h>
      header file. It must be set either to a NULL pointer or to file pointer
      which has been initialized by a fopen call. The iteration number and
      the iteration vector will be written out to the associated file, for
      each Newton iteration step. If opt->iterfile is set to NULL, no such 
      data will be written out.
      
      opt->resfile is of type pointer to FILE, as defined by the <stdio.h>
      header file. It must be set either to a NULL pointer or to file pointer
      which has been initialized by a fopen call. The iteration number and
      the residuum vector will be written out to the associated file, for
      each Newton iteration step. If opt->resfile is set to NULL, no such 
      data will be written out.
      
      opt->miscfile is of type pointer to FILE, as defined by the <stdio.h>
      header file. It must be set either to a NULL pointer or to file pointer
      which has been initialized by a fopen call. The iteration number, an
      identification number of the calling code (1 for NLEQ_RES), the norm
      of the residuum, the norm of the Newton correction, a zero value 
      as a dummy placeholder value, the accepted damping factor, and another
      zero value as a dummy placeholder value will be written out, for
      each Newton iteration step. If opt->miscfile is set to NULL, no such 
      data will be written out. For additional information on the file output,
      refer to the description of this option in the QNRES documentation.
     
      Note: The output to the files opt->iterfile, opt->resfile and
            opt->miscfile is written as a single for each iteration step.
            Such, the data in these files are suitable as input to the
            graphics utility GNUPLOT.
      
      opt->scale is of type pointer to a double array of size n. 
      This array must, if present, contain positive scaling values, which are
      used in computations of scaled norms and Jacobian scaling, as follows:
      || f || = squareroot(sum(1 to n) ( f_i/scale_i )^2)
      The pointer may be initialized with a NULL pointer. In this case, all
      scaling values are internally set to 1.
      
      opt->scaleopt is of type enum{...}. This option is only meaningful, if
      the user has not supplied a scaling vector via opt->scale.
      In this case, if opt->scaleopt is set to StandardScale, all scaling
      vector components are set to 1. If it is set to StartValueScale,
      then the setting is fscale_i = max(1,abs(f0_i)) for i=0,...,n-1,
      where f0=(f0_0,...,f0_(n-1))=problem_function(x0), x0=start vector.
      
      struct NLEQ_INFO *info:
      A pointer to an info structure. The pointed fields of the structure
      are set output info of nleq_res.
      
      info->precision is of type double and is set to the achieved scaled norm
      of the residuum at the final iterate.
      
      info->normdx is of type double and is set to the unscaled norm of the
      last Newton correction.
      
      info->fx is a pointer to a double array of size n, which contains the
      final residuum vector.
      
      info->iter is set to number of Newton iteration steps done.
      
      info->nofunevals is set to the number of done calls to the problem
      function funs.fun.
      
      info->nojacevals is set to the number of done calls to the Jacobian
      function funs.jac.
      
      info->rcode is set to the return-code of nleq_res. A return-code 0
      means that nleq_res has terminated sucessfully. For the meaning of a
      nonzero return-code, see the error messages list below.
      
      info->subcode is set for certain failure conditions to the error code
      which has been returned by a routine called from nleq_res.


      Parameter definitions of the required problem routine funs.fun
      and the optional Jacobian routine funs.jac
      --------------------------------------------------------------
      
      void fun(int *n, double *x, double *f, int *fail);
        int    *n     input  Number of vector components.
        double *x     input  Vector of unknowns, of size *n .
        double *f     output Vector of function values.
        int    *fail  output fun evaluation-failure indicator.
                      On input:  undefined.
                      On output: Indicates failure of fun evaluation, 
                                 if having a value <= 2.
                      If <0 or >2: nleq_res will be terminated with
                                   error code = 82, and *fail will be stored to
                                   info->subcode.
                      If =1: A new trial Newton iterate will
                             computed, with the damping factor
                             reduced to it's half.
                      If =2: A new trial Newton iterate will computed, with the
                             damping factor reduced by a reduction factor, which
                             must be output through f[0] by fun, and it's value
                             must be >0 and <1.
                      Note, that if IFAIL = 1 or 2, additional conditions 
                      concerning the damping factor, e.g. the minimum damping
                      factor may also influence the value of the reduced 
                      damping factor.
 
      void jac(int *n, int *ldjac, double *x, double *dfdx, int *fail);
                        Ext    Jacobian matrix subroutine
        int    *n     input  Number of vector components.
        int    *ldjac input  Leading dimension of Jacobian array, i.e.
                             the total row length for C-style two-dimensional
                             arrays, or the total column length for 
                             Fortran-style two-dimensional arrays.
                             See Note below!
        double *x     input  Vector of unknowns, of size *n .
        double *dfdx  output dfdx[i][k]: partial derivative of i-th component
                             of output parameter *f from fun with respect 
                             to x[k].
        int    *fail  output fun evaluation-failure indicator.
                      On input:  undefined.
                      On output: Indicates failure of jac evaluation
                      and causes termination of nleq_res, f set to a nonzero
                      value on output.
                                 
      Note: The calling interfaces of the user routines fun and jac has
            been designed to be compatible with routines programmed for
            use with the Fortran codes NLEQ1 and NLEQ2. However, note
            that the Fortran matrix storage mode is columnwise while
            the C matrix storage mode is rowwise. If you intend to link 
            a Jacobian routine, which has been programmed in Fortran for
            use with NLEQ1 or NLEQ2, you must either transpose the Jacobian,
            or you must compile the nleq_res package for use with Fortran
            matrix storage mode, by setting the C preprocessor flag FMAT,
            i.e. setting the gcc compiler option -DFMAT, when compiling the
            file jacobian_and_linalg.c .


      The following error conditions may occur: (returned via info->rcode)
      --------------------------------------------------------------------
      
      -999 routine nleq_fwalloc failed to allocate double memory via malloc.
      -998 routine nleq_iwalloc failed to allocate int memory via malloc.
      -997 routine nleq_pfwalloc failed to allocate double pointer memory
           via malloc.
      -995 Internal i/o control block could not be allocated via malloc.
      -994 Internally used data structure could not be allocated via malloc.
      -989 Default data-output file could not be opened via fopen call.
       -99 NULL pointer obtained from funs.fun field - the problem function
           must be defined!
         1 Singular Jacobian matrix (detected by routine nleq_linfact),
           nleq_res cannot proceed the iteration.
         2 Maximum number of Newton iteration (as set by opt->maxiter) exceeded.
         3 No convergence of Newton iteration, damping factor became too small.
        20 Nonpositive input for dimensional parameter n.
        21 Nonpositive value for opt->tol supplied.
        22 Negative scaling value for some component of vector opt->scale
           supplied.
        80 nleq_linfact returned with an error other than singular Jacobian.
           Check info->subcode for the nleq_linfact failure code.
        81 nleq_linsol returned with an error.
           Check info->subcode for the nleq_linsol failure code.
        82 The user defined problem function funs.fun returned a nonzero code
           other than 1 or 2. 
           Check info->subcode for the user-function failure code.
        83 The user defined Jacobian function funs.jac returned a nonzero code.
           Check info->subcode for the Jacobian-function failure code.
         
      Summary of changes:
      -------------------
      
      Version  Date        Changes
      1.1.1    2006/06/06  Missing int return code in function 
                           nleqres_initscale, bug fixed.
      1.1      2006/06/02  Added the output data files iterfile, resfile and
                           miscfile, where optional data output is provided,
                           a single line, starting with the iteration number,
                           for each iteration step. 
      1.0      2006/05/30  Initial release.
      
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "nleq.h"

#undef  THETA_MAX
#define THETA_MAX				0.25
#undef  LAMBDA_START_DEFAULT
#define LAMBDA_START_DEFAULT			1.0e-2
#undef  LAMBDA_START_EXTREMELY_DEFAULT
#define LAMBDA_START_EXTREMELY_DEFAULT		1.0e-4
#undef  LAMBDA_MIN_DEFAULT
#define LAMBDA_MIN_DEFAULT			1.0e-4
#undef  LAMBDA_MIN_EXTREMELY_DEFAULT
#define LAMBDA_MIN_EXTREMELY_DEFAULT		1.0e-8

#define MAX_ITER_DEFAULT			50
#define NONLIN					opt->nonlin

struct NLEQ_IO					*nleq_ioctl;

// Static routines originally defined in nleq_res.c
static int	nleqres_initscale(int n, double **scale, struct NLEQ_OPT *opt);
static void	nleqres_monitor(int k, int n, double normdx, double normf, double lambda);

// Static functions originally defined in qnres.c
static void 	qnres(struct NLEQ_FUN fun, int n, double *x, struct NLEQ_OPT *opt, struct NLEQ_INFO *info);
static int	qnr_initscale(int n, double **scale);
static void	qnr_monitor(int k, int n, double normdx, double normf, char qnres_id[9], char fail_reason[7]);

// Static routines originally defined in jacobian_and_linalg.c
static int	nleq_jacalloc(int n, double **jac, int *ldjac, struct NLEQ_OPT *opt);
static void	nleq_jacrow_scale(int n, double *jac, double *scale, struct NLEQ_OPT *opt);
static int	nleq_linfact(int n, double *jac, struct NLEQ_OPT *opt);
static void	nleq_linfree(void);
static int	nleq_linsol(int n, double *b, struct NLEQ_OPT *opt);
static int	nleq_numjac(NLEQ_FFUN *f, int n, double *x, double *fx, double *scale, double *jac, int *nfcn, struct NLEQ_OPT *opt);

// Static routines originally defined in utils.c
static void	nleq_dataout(int k, int n, double *x, struct NLEQ_DATA *data);
static int	nleq_iwalloc(int size, int **ptr, char vname[]);
static double	nleq_norm2(int n, double *v);
static int	nleq_parcheck_and_print(int n, struct NLEQ_OPT *opt, struct NLEQ_FUN fun, int nleq_code);
static int	nleq_pfwalloc(int size, double ***ptr, char vname[]);
static void	nleq_scale(int n, double *v1, double *v2, double *scale);
static double	nleq_scaled_norm2(int n, double *v, double *scale);
static double	nleq_scaled_sprod(int n, double *v1, double *v2, double *scale);


/*==========================================================================*/

void nleq_res(struct NLEQ_FUN fun, int n, double *x, struct NLEQ_OPT *opt, struct NLEQ_INFO *info)
{
  double		ftol=opt->tol;
  double		lambda, lambda_new, mue, normfk, normfkm1, normfkp1, reduction_factor, normdx, theta, s;
  double		lambda_min;
  int			i, k=0, fail=0, ldjac, max_iter=opt->maxiter;
  LOGICAL		qnres_iter=False, saved_nleqcalled=opt->nleqcalled, restricted=opt->restricted, io_allocated=False, reducted;
  PRINT_LEVEL		error_level;
  double		*dx, *fxk, *fxkp1, *fscale, *xkp1, *w;
  double		*jac;
  int			scale_allocated=0, nfcn=0, njac=0;
  NLEQ_FFUN		*f = fun.fun;
  NLEQ_JFUN		*jc = fun.jac;
  struct NLEQ_DATA	*data=malloc(sizeof(struct NLEQ_DATA));

  if (!nleq_ioctl) nleq_ioctl=malloc(sizeof(struct NLEQ_IO));
  if (!nleq_ioctl)
    {
      fprintf(stderr,"\n could not allocate output controlblock\n");
      RCODE=-995; return;
    }
  else
    io_allocated = True;
  if (!data)
    {
      fprintf(stderr,"\n could not allocate struct data\n");
      RCODE=-994; return;
    }

  data->codeid		= NLEQ_RES;
  data->normdxbar	= 0.0;
  data->theta		= 0.0;
  data->mode		= Initial;
  ERRORLEVEL		= opt->errorlevel;
  MONITORLEVEL		= opt->monitorlevel;
  DATALEVEL		= opt->datalevel;
  error_level		= opt->errorlevel;
  ERROR			= opt->errorfile;
  MONITOR		= opt->monitorfile;
  DATA			= opt->datafile;
  FITER			= opt->iterfile;
  FRES			= opt->resfile;
  FMISC			= opt->miscfile;

  if ( !ERROR && ERRORLEVEL>0 )
    ERROR   = stdout;
  if ( !MONITOR && MONITORLEVEL>0 )
    MONITOR = stdout;
  if ( !DATA && DATALEVEL>0 )
    {
      DATA=fopen("nleq_res.data","w");
      if (!DATA && ERRORLEVEL>0)
	{
	  fprintf(ERROR,"\n fopen of file nleq_res.data failed\n");
	  RCODE=-989; return;
	}
    }
  opt->errorfile	= ERROR;
  opt->monitorfile	= MONITOR;
  opt->datafile		= DATA;

  if ( MONITORLEVEL > 0 )
    fprintf(MONITOR,"\n NLEQ_RES - Version 1.1\n");
  RCODE = nleq_parcheck_and_print(n,opt,fun,1);
  if ( RCODE !=0 )
    {
      if (io_allocated)
	{
	  free(nleq_ioctl);
	  nleq_ioctl=NULL;
	}
      if (data) free(data);
      return;
    }

  opt->nleqcalled = True;
  if ( max_iter <= 0 )
    max_iter = MAX_ITER_DEFAULT;
  if ( NONLIN==Mildly_Nonlinear )
    {
      lambda     = 1.0;
      lambda_min = LAMBDA_MIN_DEFAULT;
    }
  else if ( NONLIN==Highly_Nonlinear )
    {
      lambda     = LAMBDA_START_DEFAULT;
      lambda_min = LAMBDA_MIN_DEFAULT;
    }
  else if ( NONLIN==Extremely_Nonlinear )
    {
      lambda     = LAMBDA_START_EXTREMELY_DEFAULT;
      lambda_min = LAMBDA_MIN_EXTREMELY_DEFAULT;
      restricted = True;
    } ;

  RCODE = nleq_jacalloc(n,&jac,&ldjac,opt);
  if ( RCODE !=0 ) return;
  RCODE = nleq_fwalloc(n,&dx,"dx");
  if ( RCODE !=0 ) return;
  RCODE = nleq_fwalloc(n,&xkp1,"xkp1");
  if ( RCODE !=0 ) return;
  RCODE = nleq_fwalloc(n,&fxk,"fxk");
  if ( RCODE !=0 ) return;
  RCODE = nleq_fwalloc(n,&fxkp1,"fxkp1");
  if ( RCODE !=0 ) return;
  RCODE = nleq_fwalloc(n,&w,"w");
  if ( RCODE !=0 ) return;

  data->fx = fxk;
  data->dx = dx;
  f(&n,x,fxk,&fail);
  nfcn++;
  if (fail != 0)
    {
      RCODE=82;
      goto errorexit;
    }

  fscale = opt->scale;
  if ( fscale == NULL )
    {
      RCODE=nleqres_initscale(n,&fscale,opt);
      if (RCODE !=0) goto errorexit;
      scale_allocated = 1;
      opt->scale = fscale;
    }

  if ( opt->scaleopt == StartValueScale && scale_allocated == 1 )
    for (i=0;i<n;i++) fscale[i] = MAX(fabs(fxk[i]),1.0);

  if ( MONITORLEVEL > 1 )
    fprintf(MONITOR,"\n iter      norm(dx)  norm_scl(fk)    lambda \n\n");

  normfk = nleq_scaled_norm2(n,fxk,fscale);
  normdx = 0.0;
  RCODE = 2;

  do
    {
      if ( normfk <= ftol )
	{
	  RCODE=0;
	  break;
	} 						// stop, x contains the solution!

      if (*jc != NULL)
	{
	  jc(&n,&ldjac,x,jac,&fail);  njac++;
	  if (fail !=0)
	    {
	      RCODE = 83;
	      goto errorexit;
	    }
	  }
      else
	{
	  fail=nleq_numjac(f,n,x,fxk,NULL,jac,&nfcn,opt);
	  if (fail !=0)
	    {
	      RCODE = 82;
	      goto errorexit;
	    }
	}

      // scale jacobian and change sign of it
      nleq_jacrow_scale(n,jac,fscale,opt);
      fail = nleq_linfact(n,jac,opt);  			// compute LU-factorization of Jacobian

      if ( fail < 0 )
	{
	  RCODE=80;
	  goto errorexit;
	}
      else if ( fail > 0 )
	{
	  RCODE=1;
	  goto errorexit;
	}

      // compute newton correction
      nleq_scale(n,fxk,dx,fscale);
      fail = nleq_linsol(n,dx,opt);

      if ( fail != 0 )
	{
	  RCODE=81;
	  goto errorexit;
	}
      normdx = nleq_norm2(n,dx);
      if ( k>0 )
	{
	  mue = (normfkm1/normfk)*mue;
	  lambda = MIN(1.0,mue);
	}
      if ( MONITORLEVEL > 1 )
	normfkp1 = normfk; // set normfkp1 for print only
      reducted = False;
      checkregularity:

      if ( MONITORLEVEL > 1 && fail==0 )
	nleqres_monitor(k,n,normdx,normfkp1,lambda);
      if ( lambda < lambda_min )
	{
	  RCODE=3;
	  break;
	}  // stop, convergence failure!

      for (i=0;i<n;i++)
	xkp1[i]=x[i]+lambda*dx[i]; 			// new trial iterate

      f(&n,xkp1,fxkp1,&fail);
      nfcn++;
      if ( fail<0 || fail>2 )
	{
	  RCODE=82;
	  goto errorexit;
	}
      else if ( fail==1 || fail== 2 )
	{
	  if ( fail==1 )
	    reduction_factor = 0.5;
	  else
	    reduction_factor = fxkp1[0];
	  if ( reduction_factor <= 0.0 || reduction_factor >= 1.0 )
	    {
	      RCODE=82;
	      goto errorexit;
	    }
	  if ( MONITORLEVEL>1 )
	    fprintf(MONITOR," %4i  FUN could not be evaluated  %7f\n", k,lambda);
	  if ( lambda > lambda_min )
	    lambda = MAX(lambda*reduction_factor,lambda_min);
	  else
	    lambda = lambda*reduction_factor;
	  reducted = True;
	  goto checkregularity;
	}

      normfkp1 = nleq_scaled_norm2(n,fxkp1,fscale);
      theta = normfkp1/normfk;
      s = 1.0-lambda;
      for (i=0;i<n;i++)
	w[i] = fxkp1[i]-s*fxk[i];
      mue = (0.5*normfk*lambda*lambda) / nleq_scaled_norm2(n,w,fscale);

      if ( ( !restricted && theta >= 1.0 ) ||
	  (  restricted && theta > 1.0-lambda/4.0) )
	{
	  lambda_new = MIN(mue,0.5*lambda);
	  if ( lambda <= lambda_min )
	    lambda = lambda_new;
	  else
	    lambda = MAX(lambda_new,lambda_min);
	  reducted = True;
	  goto checkregularity;
	}
      lambda_new = MIN(1.0,mue);
      if ( lambda==1.0 && lambda_new==1.0 && theta < THETA_MAX )
	qnres_iter = True;
      else
	{
	  if( lambda_new >= 4.0*lambda && !reducted )
	    {
	      lambda=lambda_new;
	      goto checkregularity;
	    }
	}
      data->normf  = normfk;
      data->normdx = normdx;
      data->lambda = lambda;

      nleq_dataout(k,n,x,data);
      data->mode = Intermediate;
      for (i=0;i<n;i++)
	x[i]=xkp1[i];					// accept new iterate

      // next step
      k++;
      normfkm1 = normfk;
      normfk   = normfkp1;

      // start qnres only, if residuum not yet small enough
      qnres_iter = qnres_iter && normfk > ftol;
      // perform qnres steps, if choosen
      if (qnres_iter)
	{
	  info->fx        = fxk;
	  info->iter      = k;
	  info->normdx    = normdx;
	  opt->maxiter    = max_iter-k+1;
	  opt->errorlevel = 0;

	  qnres(fun,n,x,opt,info);
	  nfcn += info->nofunevals;
	  k = info->iter;
	  normfkp1=info->precision;
	  opt->errorlevel = error_level;
	  ERRORLEVEL      = error_level;

	  // if QNRES failed, try to continue NLEQ_RES
	  if ( RCODE != 0 )
	    {
	      RCODE = 2;
	      qnres_iter=False;
	    }
	}
      else
	for (i=0;i<n;i++)
	  fxk[i]=fxkp1[i];
    }
  while ( k <= max_iter && RCODE == 2 );

  if ( !qnres_iter )
    {
      if ( MONITORLEVEL > 1 )
	nleqres_monitor(k,n,normdx,normfk,lambda);
      data->normf  = normfk;
      data->normdx = normdx;
      data->lambda = lambda;
      data->mode   = ( RCODE==0 ? Solution : Final );
      nleq_dataout(k,n,x,data);
    }

errorexit:
  if ( ERRORLEVEL > 0 && RCODE != 0 )
    {
      switch ( RCODE )
      {
	case     1:
	  fprintf(ERROR,"\n Error return from nleq_linfact: Singular Jacobian\n");
	  break;
	case     2:
	  fprintf(ERROR,"\n Error - Maximum allowed number of iterations exceeded\n");
	  break;
	case     3:
	  fprintf(ERROR,"\n Error - no convergence, damping factor became too small\n");
	  break;
	case    80:
	  fprintf(ERROR,"\n Error return from nleq_linfact: fail=%i\n",fail);
	  break;
	case    81:
	  fprintf(ERROR,"\n Error return from nleq_linsol: fail=%i\n",fail);
	  break;
	case    82:
	  fprintf(ERROR,"\n Error return from problem function: fail=%i\n",
	      fail);
	  break;
	case    83:
	  fprintf(ERROR,"\n Error return from Jacobian function: fail=%i\n",
	      fail);
	  break;
	default   :
	  fprintf(ERROR,"\n Error, code=%i,  subcode=%i\n",RCODE,fail);
      }
    }

  info->subcode = fail;
  if (io_allocated)
    {
      free(nleq_ioctl);
      nleq_ioctl=NULL;
    }

  free(data);
  free(dx);
  free(xkp1);
  if (qnres_iter) free(fxk);
  free(fxkp1);
  free(w);
  if (scale_allocated)
    {
      free(fscale);
      opt->scale = NULL;
    }

  nleq_linfree();
  if (!qnres_iter)
    {
      info->precision = normfk;
      info->normdx    = normdx;
    }

  // restore original values
  opt->maxiter     = max_iter;
  opt->nleqcalled  = saved_nleqcalled;
  info->iter       = k;
  info->nofunevals = nfcn;
  info->nojacevals = njac;
  if (!qnres_iter)
    info->fx       = fxk;
}


/*==========================================================================*/

static int nleqres_initscale(int n, double **scale, struct NLEQ_OPT *opt)
{
  int				i, rcode;
  rcode = nleq_fwalloc(n,scale,"scale");
  if ( rcode != 0 ) return rcode;
  for (i=0;i<n;i++)
    (*scale)[i]=1.0;
  return 0;
}


/*==========================================================================*/

static void nleqres_monitor(int k, int n, double normdx, double normf, double lambda)
{
  fprintf(MONITOR," %4i  %12e  %12e  %7f\n",k,normdx,normf,lambda);
  return;
}


/*============================================================================
 *
 *  QNRES - RESidual-based Quasi-Newton algorithm
 *
 *  Routines originally defined in qnres.c
 *
 *============================================================================
 *
 *  Written by        L. Weimann
 *  Purpose           Solution of systems of nonlinear equations
 *  Method            RESidual-based Quasi-Newton algorithm
                      (see reference below)
 *  Category          F2a. - Systems of nonlinear equations
 *  Keywords          Nonlinear equations, Newton methods
 *  Version           1.1.1
 *  Revision          May 2006
 *  Latest Change     June 2006
 *  Library           NewtonLib
 *  Code              C, Double Precision
 *  Environment       Standard C environment on PC's,
                      workstations and hosts.
 *  Copyright     (c) Konrad-Zuse-Zentrum fuer
                      Informationstechnik Berlin (ZIB)
                      Takustrasse 7, D-14195 Berlin-Dahlem
                      phone : + 49/30/84185-0
                      fax   : + 49/30/84185-125
 *  Contact           Lutz Weimann
                      ZIB, Division Scientific Computing,
                           Department Numerical Analysis and Modelling
                      phone : + 49/30/84185-185
                      fax   : + 49/30/84185-107
                      e-mail: weimann@zib.de

 *    References:

      /1/ P. Deuflhard:
          Newton Methods for Nonlinear Problems. -
          Affine Invariance and Adaptive Algorithms.
          Series Computational Mathematics 35, Springer (2004)

   ---------------------------------------------------------------

 * Licence
     You may use or modify this code for your own non commercial
     purposes for an unlimited time.
     In any case you should not deliver this code without a special
     permission of ZIB.
     In case you intend to use the code commercially, we oblige you
     to sign an according licence agreement with ZIB.

 * Warranty
     This code has been tested up to a certain level. Defects and
     weaknesses, which may be included in the code, do not establish
     any warranties by ZIB. ZIB does not take over any liabilities
     which may follow from acquisition or application of this code.

 * Software status
     This code is under care of ZIB and belongs to ZIB software class 2.

      ------------------------------------------------------------

 *    Parameters description
      ======================

      The calling interface looks as follows:

      extern void qnres(struct NLEQ_FUN funs, int n, double *x,
                        struct NLEQ_OPT *opt,
                        struct NLEQ_INFO *info)

      The structures used within the parameter list are defined
      as follows:
      ---
      struct NLEQ_FUN
      {
        NLEQ_FFUN *fun;
        NLEQ_JFUN *jac;
      };

      where the types used within this structure are defined by
      typedef void NLEQ_FFUN(int*,double*,double*,int*);
      and
      typedef void NLEQ_JFUN(int*,int*,double*,double*,int*);
      ---
      struct NLEQ_OPT
      {
         double rtol;
         int maxiter;
         LOGICAL restricted, nleqcalled;
         PRINT_LEVEL errorlevel, monitorlevel, datalevel;
         FILE *errorfile, *monitorfile, *datafile;
         PROBLEM_TYPE nonlin;
         double *scale;
      };

      where the types used within this structure are defined by
      typedef enum {None=0, Minimum=1, Verbose=2, Debug=3} PRINT_LEVEL;
      typedef enum {False=0, True=1} LOGICAL ;
      typedef enum {Mildly_Nonlinear=2, Highly_Nonlinear=3,
                    Extremely_Nonlinear=4} PROBLEM_TYPE ;
      ---
      struct NLEQ_INFO
      {
         double precision, normdx;
         double *fx;
         int iter, rcode, subcode, nofunevals, nojacevals;
      };
      ---

      A detailed description of the parameters follows:

      struct NLEQ_FUN funs :

      The field funs.fun must contain a pointer to the problem function fun -
      The required parameters interface of fun is described in detail below

      The field funs.jac must either contain a pointer to the Jacobian function
      jac or a NULL pointer. If a NULL pointer is supplied, then the Jacobian
      will be approximately computed by an internal function of the qnres
      package.

      int n :
      The number of equations and unknown variables of the nonlinear system.

      double *x :
      A pointer to an array of double values of size n.
      The array must contain on input an initial guess of the problems solution,
      which is used as the start-vector of the Quasi-Newton iteration.
      On output, the pointed array contains an approximate solution vector x*,
      which fits the small residual condition
      || funs.fun(x*) || <= opt->tol,
      where ||...|| is a scaled Euclidian norm.

      struct NLEQ_OPT *opt:
      A pointer to an options structure. The pointed fields of the structure
      contain input options to qnres.

      opt->tol is of type double and must contain the residuum threshold
      which funs.fun(x*) must fit for the solution vector x*.

      opt->maxiter is of type int and must contain the maximum number of allowed
      iterations. if a zero or negative value is supplied, then the maximum
      iteration count will be set to 50.

      opt->nleqcalled is of type LOGICAL and only used internally. This field
      must always be set to False.

      opt->errorlevel is of type PRINT_LEVEL. If it is set to level None,
      then no error message will be printed if an error condition occurs.
      If it is set to level Minimum or any higher level, then error messages
      will be printed, if appropriate.

      opt->monitorlevel is of type PRINT_LEVEL. If it is set to level None,
      then no monitor output will be printed.
      If it is set to level Minimum, a few infomation will be printed.
      If set to level Verbose, then some infomation about each Quasi-Newton
      iteration step, fitting into a single line, will be printed. The higher
      level Debug is reserved for future additional information output.

      opt->datalevel is of type PRINT_LEVEL. If it is set to level None,
      then no data output will be printed.
      If it is set to level Minimum, then the values of the initial iteration
      vector x and the final vector x will be printed.
      If set to level Verbose, then the iteration vector x will be printed for
      each Quasi-Newton step. The higher level Debug is reserved for future
      additional information output.

      opt->errorfile is of type pointer to FILE, as defined by the <stdio.h>
      header file. It must be set either to a NULL pointer, stderr, stdout,
      or to another file pointer which has been initialized by a fopen call.
      If it is set to NULL, opt->errorfile will be set to stdout. The error
      messages will be printed to opt->errorfile.

      opt->monitorfile is of type pointer to FILE, as defined by the <stdio.h>
      header file. It must be set either to a NULL pointer, stderr, stdout,
      or to another file pointer which has been initialized by a fopen call.
      If it is set to NULL, opt->monitorfile will be set to stdout. The monitor
      output will be printed to opt->monitorfile.

      opt->datafile is of type pointer to FILE, as defined by the <stdio.h>
      header file. It must be set either to a NULL pointer, stderr, stdout,
      or to another file pointer which has been initialized by a fopen call.
      If it is set to NULL, a file named "qnres.data" will be opened by a
      fopen call and opt->datafile will be set to the filepointer which the
      fopen returns. The data output will be printed to opt->datafile.

      opt->iterfile is of type pointer to FILE, as defined by the <stdio.h>
      header file. It must be set either to a NULL pointer or to file pointer
      which has been initialized by a fopen call. The iteration number and
      the iteration vector will be written out to the associated file, for
      each Quasi-Newton iteration step. If opt->iterfile is set to NULL,
      no such data will be written out.

      opt->resfile is of type pointer to FILE, as defined by the <stdio.h>
      header file. It must be set either to a NULL pointer or to file pointer
      which has been initialized by a fopen call. The iteration number and
      the residuum vector will be written out to the associated file, for
      each Quasi-Newton iteration step. If opt->resfile is set to NULL,
      no such data will be written out.

      opt->miscfile is of type pointer to FILE, as defined by the <stdio.h>
      header file. It must be set either to a NULL pointer or to file pointer
      which has been initialized by a fopen call. The iteration number, an
      identification number of the calling code (0 for QNRES), the norm
      of the residuum, the norm of the Quasi-Newton correction, a zero value
      as a dummy placeholder value, a 1.0 as damping factor value, and
      the control parameter theta will be written out, for each Quasi-Newton
      iteration step. If opt->miscfile is set to NULL, no such data will be
      written out.

      Note: The output to the files opt->iterfile, opt->resfile and
            opt->miscfile is written as a single for each iteration step.
            Such, the data in these files are suitable as input to the
            graphics utility GNUPLOT.

      opt->scale is of type pointer to a double array of size n.
      This array must, if present, contain positive scaling values, which are
      used in computations of scaled norms and Jacobian scaling, as follows:
      || f || = squareroot(sum(1 to n) ( f_i/scale_i )^2)
      The pointer may be initialized with a NULL pointer. In this case, all
      scaling values are internally set to 1.

      opt->scaleopt is of type enum{...}. This option is only meaningful, if
      the user has not supplied a scaling vector via opt->scale.
      In this case, if opt->scaleopt is set to StandardScale, all scaling
      vector components are set to 1. If it is set to StartValueScale,
      then the setting is fscale_i = max(1,abs(f0_i)) for i=0,...,n-1,
      where f0=(f0_0,...,f0_(n-1))=problem_function(x0), x0=start vector.

      struct NLEQ_INFO *info:
      A pointer to an info structure. The pointed fields of the structure
      are set output info of qnres.

      info->precision is of type double and is set to the achieved scaled norm
      of the residuum at the final iterate.

      info->normdx is of type double and is set to the unscaled norm of the
      last Quasi-Newton correction.

      info->fx is a pointer to a double array of size n, which contains the
      final residuum vector.

      info->iter is set to number of Quasi-Newton iteration steps done.

      info->nofunevals is set to the number of done calls to the problem
      function funs.fun.

      info->nojacevals is set to the number of done calls to the Jacobian
      function funs.jac.

      info->rcode is set to the return-code of qnres. A return-code 0
      means that qnres has terminated sucessfully. For the meaning of a
      nonzero return-code, see the error messages list below.

      info->subcode is set for certain failure conditions to the error code
      which has been returned by a routine called from qnres.


      Parameter definitions of the required problem routine funs.fun
      and the optional Jacobian routine funs.jac
      --------------------------------------------------------------

      void fun(int *n, double *x, double *f, int *fail);
        int    *n     input  Number of vector components.
        double *x     input  Vector of unknowns, of size *n .
        double *f     output Vector of function values.
        int    *fail  output fun evaluation-failure indicator.
                      On input:  undefined.
                      On output: Indicates failure of fun evaluation,
                                 if set to a nonzero value. In this case,
                                 qnres will be terminated with error code=82,
                                 and *fail will be stored to info->subcode.

      void jac(int *n, int *ldjac, double *x, double *dfdx, int *fail);
                        Ext    Jacobian matrix subroutine
        int    *n     input  Number of vector components.
        int    *ldjac input  Leading dimension of Jacobian array, i.e.
                             the total row length for C-style two-dimensional
                             arrays, or the total column length for
                             Fortran-style two-dimensional arrays.
                             See Note below!
        double *x     input  Vector of unknowns, of size *n .
        double *dfdx  output dfdx[i][k]: partial derivative of i-th component
                             of output parameter *f from fun with respect
                             to x[k].
        int    *fail  output fun evaluation-failure indicator.
                      On input:  undefined.
                      On output: Indicates failure of jac evaluation,
                                 if set to a nonzero value. In this case,
                                 qnres will be terminated with error code=83,
                                 and *fail will be stored to info->subcode.

      Note: The calling interfaces of the user routines fun and jac has
            been designed to be compatible with routines programmed for
            use with the Fortran codes NLEQ1 and NLEQ2. However, note
            that the Fortran matrix storage mode is columnwise while
            the C matrix storage mode is rowwise. If you intend to link
            a Jacobian routine, which has been programmed in Fortran for
            use with NLEQ1 or NLEQ2, you must either transpose the Jacobian,
            or you must compile the qnres package for use with Fortran
            matrix storage mode, by setting the C preprocessor flag FMAT,
            i.e. setting the gcc compiler option -DFMAT, when compiling the
            file jacobian_and_linalg.c .


      The following error conditions may occur: (returned via info->rcode)
      --------------------------------------------------------------------

      -999 routine nleq_fwalloc failed to allocate double memory via malloc.
      -998 routine nleq_iwalloc failed to allocate int memory via malloc.
      -997 routine nleq_pfwalloc failed to allocate double pointer memory
           via malloc.
      -995 Internal i/o control block could not be allocated via malloc.
      -994 Internally used data structure could not be allocated via malloc.
      -989 Default data-output file could not be opened via fopen call.
       -99 NULL pointer obtained from funs.fun field - the problem function
           must be defined!
         1 Singular Jacobian matrix (detected by routine nleq_linfact),
           qnres cannot proceed the iteration.
         2 Maximum number of Quasi-Newton iteration (as set by opt->maxiter)
           exceeded.
         3 No convergence of Quasi-Newton iteration, theta became too large.
         6 Ill conditioned update (kappa became too large)
        80 nleq_linfact returned with an error other than singular Jacobian.
           Check info->subcode for the nleq_linfact failure code.
        81 nleq_linsol returned with an error.
           Check info->subcode for the nleq_linsol failure code.
        82 The user defined problem function funs.fun returned a nonzero code.
           Check info->subcode for the user-function failure code.
        83 The user defined Jacobian function funs.jac returned a nonzero code.
           Check info->subcode for the Jacobian-function failure code.

      Summary of changes:
      -------------------

      Version  Date        Changes
      1.1.1    2006/06/06  Missing int return code in function
                           qnr_initscale, bug fixed.
      1.1      2006/06/02  Added the output data files iterfile, resfile and
                           miscfile, where optional data output is provided,
                           a single line, starting with the iteration number,
                           for each iteration step.
      1.0      2006/05/30  Initial release.

*/

#define THETA_MAX				0.25
#undef  KAPPA_MAX
#define KAPPA_MAX				1.0e5
#define MAX_ITER_DEFAULT			50


/*==========================================================================*/

static void qnres(struct NLEQ_FUN fun, int n, double *x, struct NLEQ_OPT *opt, struct NLEQ_INFO *info)
{
  int			inc_one=1;
  double		one=1.0;
  double		sigmak, sigmakp1, kappa, thetak, beta, ftol=opt->tol;
  double		s, normfk, normdx;
  int			i, j, k=0, k1, kprint, fail=0, ldjac, max_iter=opt->maxiter;
  LOGICAL		nleqcalled = opt->nleqcalled, skipstep, io_allocated = False;
  double		*dfxkp1, *fxk, *fxkp1, *v, *w;
  double		*dx, *fscale, *fgamma;
  double		**fx, **dfx;
  double		*jac;
  int			scale_allocated=0, nfcn=0, njac=0;
  char			qnres_id[9]="        \0", fail_reason[7]="      \0";
  NLEQ_FFUN		*f = fun.fun;
  NLEQ_JFUN		*jc = fun.jac;
  struct NLEQ_DATA	*data=malloc(sizeof(struct NLEQ_DATA));

  if (!nleq_ioctl)
    nleq_ioctl=malloc(sizeof(struct NLEQ_IO));
  if (!nleq_ioctl)
    {
      fprintf(stderr,"\n could not allocate output controlblock\n");
      RCODE=-995;
      return;
    }
  else
    io_allocated = True;
  if (!data)
    {
      fprintf(stderr,"\n could not allocate struct data\n");
      RCODE=-994;
      return;
    }

  data->codeid    = QNRES;
  data->normdxbar = 0.0;
  data->lambda    = 1.0;
  ERRORLEVEL      = opt->errorlevel;
  MONITORLEVEL    = opt->monitorlevel;
  DATALEVEL       = opt->datalevel;
  ERROR           = opt->errorfile;
  MONITOR  	  = opt->monitorfile;
  DATA     	  = opt->datafile;
  FITER    	  = opt->iterfile;
  FRES     	  = opt->resfile;
  FMISC    	  = opt->miscfile;

  if ( !ERROR && ERRORLEVEL>0 )
    ERROR   = stdout;
  if ( !MONITOR && MONITORLEVEL>0 )
    MONITOR = stdout;
  if ( !DATA && DATALEVEL>0 )
    {
      DATA=fopen("qnres.data","w");
      if (!DATA && ERRORLEVEL>0)
	{
	  fprintf(ERROR,"\n fopen of file qnres.data failed\n");
	  RCODE=-989;
	  return;
	}
    }
  opt->errorfile   = ERROR;
  opt->monitorfile = MONITOR;
  opt->datafile    = DATA;

  if ( !nleqcalled )
    {
      if ( MONITORLEVEL > 0 )
	fprintf(MONITOR,"\n QNRES - Version 1.1\n");
      RCODE = nleq_parcheck_and_print(n,opt,fun,0);
      if ( RCODE !=0 )
	{
	  if (io_allocated)
	    {
	      free(nleq_ioctl);
	      nleq_ioctl=NULL;
	    }
	  if (data) free(data);
	  return;
	}
    }

  if ( max_iter <= 0 )
    max_iter = MAX_ITER_DEFAULT;

  RCODE = nleq_pfwalloc(max_iter+2,&fx,"fx");
  if ( RCODE !=0 ) return;
  RCODE = nleq_pfwalloc(max_iter+2,&dfx,"dfx");
  if ( RCODE !=0 ) return;
  RCODE = nleq_fwalloc(max_iter+2,&fgamma,"fgamma");
  if ( RCODE !=0 ) return;
  RCODE = nleq_fwalloc(n,&dx,"dx");
  if ( RCODE !=0 ) return;

  data->dx = dx;
  fscale = opt->scale;
  if ( fscale == NULL )
    {
      RCODE=qnr_initscale(n,&fscale);
      if (RCODE !=0) goto errorexit;
      scale_allocated = 1;
    }
  skipstep=nleqcalled;
  if(!skipstep)
    {
      RCODE = nleq_jacalloc(n,&jac,&ldjac,opt);
      if ( RCODE !=0 ) return;
      RCODE = nleq_fwalloc(n,&fx[0],"fx[0]");
      if ( RCODE !=0 ) return;
      f(&n,x,fx[0],&fail);
      nfcn++;
      if (fail != 0)
	{
	  RCODE=82;
	  goto errorexit;
	}

      fxk = fx[0];
      if ( opt->scaleopt == StartValueScale && scale_allocated == 1 )
	for (i=0;i<n;i++)
	  fscale[i] = MAX(fabs(fxk[i]),1.0);
      if (*jc != NULL)
	{
	  jc(&n,&ldjac,x,jac,&fail);
	  njac++;
	  if (fail !=0)
	    {
	      RCODE = 83;
	      goto errorexit;
	    }
	}
      else
	{
	  fail = nleq_numjac(f,n,x,fx[0],NULL,jac,&nfcn,opt);
	  if (fail !=0)
	    {
	      RCODE = 82;
	      goto errorexit;
	    }
	}
      sigmak = nleq_scaled_sprod(n,fx[0],fx[0],fscale);
      normfk = sqrt(sigmak);

      // scale jacobian and change sign of it
      nleq_jacrow_scale(n,jac,fscale,opt);
      fail = nleq_linfact(n,jac,opt);				// compute LU-factorization of Jacobian
      if ( fail < 0 )
	{
	  RCODE=80;
	  goto errorexit;
	}
      else if ( fail > 0 )
	{
	  RCODE=1;
	  goto errorexit;
	}

      nleq_scale(n,fx[0],dx,fscale);
      fail = nleq_linsol(n,dx,opt);				// compute first quasi newton correction
      if ( fail != 0 )
	{
	  RCODE=81;
	  goto errorexit;
	}
      normdx = nleq_norm2(n,dx);
      data->normf  = normfk;
      data->normdx = normdx;
      data->theta  = 0.0;
      data->mode   = Initial;
      nleq_dataout(k,n,x,data);
      data->mode   = Intermediate;

      if ( MONITORLEVEL > 1 )
	{
	  fprintf(MONITOR,"\n iter      norm(dx)  norm_scl(fk)\n\n");
	  qnr_monitor(k,n,normdx,normfk,qnres_id,fail_reason);
	}
    }
  else
    {
      strcpy(qnres_id,"   QNRES\0");
      fx[0] = info->fx;
      sigmak = nleq_scaled_sprod(n,fx[0],fx[0],fscale);
      data->mode = Intermediate;
      normdx = info->normdx;
    }
  kappa = 1.0;

  do
    {
      if (!skipstep)
	cblas_daxpy(n, one, dx, inc_one, x, inc_one);		// x(k+1)=x(k)+dx(k)
      skipstep = False;
      RCODE = nleq_fwalloc(n,&fx[k+1],"fx[k+1]");
      if ( RCODE !=0 ) return;
      RCODE = nleq_fwalloc(n,&dfx[k+1],"dfx[k+1]");
      if ( RCODE !=0 ) return;
      RCODE = 2;

      f(&n,x,fx[k+1],&fail);
      nfcn++;
      if (fail != 0)
	{
	  RCODE=82;
	  break;
	}
      dfxkp1   = dfx[k+1];
      fxkp1    = fx[k+1];
      fxk      = fx[k];
      data->fx = fxk;

      for (i=0;i<n;i++)
	dfxkp1[i] = fxkp1[i]-fxk[i];

      sigmakp1 = nleq_scaled_sprod(n,fxkp1,fxkp1,fscale);
      normfk = sqrt(sigmakp1);

      if ( sigmakp1 <= ftol*ftol )
	{
	  k++;
	  RCODE=0;
	  break;}						// stop, x contains the solution!

      thetak = sqrt(sigmakp1/sigmak);
      if ( thetak >= THETA_MAX )
	{
	  k++;
	  RCODE=3;
	  break;
	}							// stop, no convergence!

      w=dfxkp1;
      v=dx;
      fgamma[k] = nleq_scaled_sprod(n,w,w,fscale);
      kappa = kappa/(1.0-2.0*thetak);
      if ( kappa >= KAPPA_MAX )
	{
	  k++;
	  RCODE=6;
	  break;
	}							// stop, ill conditioned update!

      s = 1.0 - nleq_scaled_sprod(n,w,fxkp1,fscale) / fgamma[k];
      for (i=0;i<n;i++)
	v[i] = s*fxkp1[i];
      for (j=k-1;j>=0;j--)
	{
	  beta = - nleq_scaled_sprod(n,dfx[j+1],v,fscale)/fgamma[j];
	  cblas_daxpy(n, beta, fx[j+1], inc_one, v, inc_one);		// v=v+beta*fx[j+1]
	}
      nleq_scale(n,v,v,fscale);

      // compute new quasi newton correction
      fail = nleq_linsol(n,v,opt);
      if ( fail != 0 )
	{
	  RCODE=81;
	  break;
	}
      normdx = nleq_norm2(n,v);
      k++;
      kprint = ( nleqcalled ? (info->iter)+k-1 : k );
      data->normf  = normfk;
      data->normdx = normdx;
      data->theta  = thetak;
      nleq_dataout(kprint,n,x,data);

      if ( MONITORLEVEL > 1 )
	qnr_monitor(kprint,n,normdx,normfk,qnres_id,fail_reason);
      sigmak = sigmakp1;
    }
  while ( k <= max_iter && RCODE == 2 );

  kprint = ( nleqcalled ? (info->iter)+k-1 : k );
  if ( RCODE != 2 )
    { if ( MONITORLEVEL > 1 )
      { if (nleqcalled)
	{  if ( RCODE==3 ) strcpy(fail_reason,"THETA!\0");
	else if ( RCODE==6 ) strcpy(fail_reason,"KAPPA!\0");
	}
      qnr_monitor(kprint,n,normdx,normfk,qnres_id,fail_reason);
      }
    }

  data->normf = normfk;  data->normdx = normdx;  data->theta = thetak;
  data->mode = ( RCODE==0 ? Solution : Final );
  if ( RCODE==0 || !nleqcalled ) nleq_dataout(kprint,n,x,data);

  errorexit:
  if ( ERRORLEVEL > 0 && RCODE != 0 )
    {
      switch ( RCODE )
      {
	case     1:
	  fprintf(ERROR,"\n Error return from nleq_linfact: Singular Jacobian\n");
	  break;
	case     2:
	  fprintf(ERROR,"\n Error - Maximum allowed number of iterations exceeded\n");
	  break;
	case     3:
	  fprintf(ERROR,"\n Error - QNRES failed - theta = %e\n",thetak);
	  break;
	case     6:
	  fprintf(ERROR,"\n Error - QNRES failed - kappa = %e\n",kappa);
	  break;
	case    80:
	  fprintf(ERROR,"\n Error return from nleq_linfact: fail=%i\n",fail);
	  break;
	case    81:
	  fprintf(ERROR,"\n Error return from nleq_linsol: fail=%i\n",fail);
	  break;
	case    82:
	  fprintf(ERROR,"\n Error return from problem function\n");
	  break;
	case    83:
	  fprintf(ERROR,"\n Error return from Jacobian evaluation function\n");
	  break;
	default   :
	  fprintf(ERROR,"\n Error, code=%i,  subcode=%i\n",RCODE,fail);
      }
    }
  info->subcode = fail;
  if ( !nleqcalled && io_allocated ) { free(nleq_ioctl); nleq_ioctl=NULL; }
  free(data);  free(fgamma);
  if (scale_allocated) free(fscale);
  k1 = ( RCODE == 2 ? k : k+1 );
  info->fx         = fx[k1];
  for (i=1;i<=k;i++) {free(fx[i]);  free(dfx[i]);}
  if ( k1==k+1 ) free(dfx[k1]);
  if (!nleqcalled) {free(fx[0]); nleq_linfree();}
  free(dx);  free(fx);  free(dfx);
  info->precision  = sqrt(sigmakp1);
  info->normdx     = normdx;
  info->iter       = (nleqcalled ? (info->iter)+k-1 : k );
  info->nofunevals = nfcn;
  info->nojacevals = njac;
}


/*==========================================================================*/

static int qnr_initscale(int n, double **scale)
{
  int			i, rcode;
  rcode = nleq_fwalloc(n,scale,"scale");
  if ( rcode !=0 ) return rcode;
  for (i=0;i<n;i++)
    (*scale)[i]=1.0;
  return 0;
}


/*==========================================================================*/

static void qnr_monitor(int k, int n, double normdx, double normf, char qnres_id[9], char fail_reason[7])
{
  fprintf(MONITOR," %4i  %12e  %12e  %s  %s\n",k,normdx,normf,qnres_id,fail_reason);
  return;
}


/*============================================================================
 *
 *  Jacobian related routines for the NewtonLib package.
 *
 *  Routines originally defined in jacobian_and_linalg.c
 *
 *============================================================================
 *
 *  Written by        L. Weimann
 *  Purpose           Jacobian and linear algebra related routines
 *  Category          ???. - Utilities
 *  Keywords          Linear system solution, Numerical differentiation,
                      Jacobian scaling
 *  Version           1.1
 *  Revision          May 2006
 *  Latest Change     May 2006
 *  Library           NewtonLib
 *  Code              C, Double Precision
 *  Environment       Standard C environment on PC's,
                      workstations and hosts.
 *  Copyright     (c) Konrad-Zuse-Zentrum fuer
                      Informationstechnik Berlin (ZIB)
                      Takustrasse 7, D-14195 Berlin-Dahlem
                      phone : + 49/30/84185-0
                      fax   : + 49/30/84185-125
 *  Contact           Lutz Weimann
                      ZIB, Division Scientific Computing,
                           Department Numerical Analysis and Modelling
                      phone : + 49/30/84185-185
                      fax   : + 49/30/84185-107
                      e-mail: weimann@zib.de

   ---------------------------------------------------------------

 * Licence
     You may use or modify this code for your own non commercial
     purposes for an unlimited time.
     In any case you should not deliver this code without a special
     permission of ZIB.
     In case you intend to use the code commercially, we oblige you
     to sign an according licence agreement with ZIB.

 * Warranty
     This code has been tested up to a certain level. Defects and
     weaknesses, which may be included in the code, do not establish
     any warranties by ZIB. ZIB does not take over any liabilities
     which may follow from acquisition or application of this code.

 * Software status
     This code is under care of ZIB and belongs to ZIB software class 2.

   ------------------------------------------------------------

   The implementation of the routines within this file depend on the
   chosen Jacobian format and the selected linear system solver.
   The current implementation supports Jacobians stored in full
   storage mode, and uses the Fortran LAPACK solver routines DGETRF
   for the LU-decomposition and DGETRS for the solution (i.e.
   backsubstitution).

   The routines and functions are the following:
   ---------------------------------------------

   int nleq_jacalloc(int n, double **jac, int *ldjac,
                     struct NLEQ_OPT *opt)
   This function must allocate memory for the Jacobian, using either
   directly the malloc function or the functions nleq_fwalloc,
   nleq_iwalloc, nleq_pfwalloc (see file "utils.c"), as suitable.
   The parameters are:
   int      n        (input)   The number of equations and unknowns
   double   **jac    (output)  A pointer of type (double *) returned
                               by the memory allocation function, and
                               pointing to the reserved memory.
   int     *ldjac    (output)  A positive int number specifying the
                               leading dimension size of the Jacobian
                               (asuming the Jacobian to be a two dimensional
                                array of (double)).
   struct NLEQ_OPT *opt (input) The input options of the main routine.

   The int return code of nleq_jacalloc must be 0, if everything worked fine
   with the memory allocation, and otherwise a number from -999 to -996.
   ---
   Note: You may expand the struct NLEQ_OPT, which is defined in the "nleq.h"
         header file, by adding fields which may describe characteristics
         of your Jacobian, for example, such as lower and upper bandwidth
         of a Jacobian which will be stored in band mode.
   ---
   int nleq_linfact(int n, double *jac, struct NLEQ_OPT *opt)
   This function must call the linear algebra routine for computation of
   the (LU-)decomposition.
   The parameters are:
   int      n        (input)   The number of equations and unknowns
   double   *jac     (in/out)  On input:  The Jacobian
                               On output: The decomposition of the Jacobian.
   struct NLEQ_OPT *opt (input) The input options of the main routine.

   The int return code must be 0, if the matrix decomposition was successful.
   A singular matrix must be indicated by a positive return code, any other
   error by a negative return code.
   ---
   Note: To keep all needed information of the matrix decomposition, it
         may be necessary to allocate additional memory on the first call
         of nleq_linfact, i.e. to store pivot indices or other information.
   ---
   int nleq_linsol(int n, double *b, struct NLEQ_OPT *opt)
   This function must call the linear algebra routine for the solution
   (i.e. backsubstitution) of the linear system, after the Jacobian has
   been factorized through nleq_linfact.
   The parameters are:
   int      n        (input)   The number of equations and unknowns
   double   *b       (in/out)  On input: The right hand side vector of the
                                         linear equations system.
                               On output: The solution of the linear system.
   struct NLEQ_OPT *opt (input) The input options of the main routine.

   A successful completion must be indicated by a return code 0, a failure
   condition by any nonzero return code.
   ---
   void nleq_linfree()
   This routine must free any memory, which has been allocated either by
   the routine nleq_jacalloc or the routine nleq_linfact.
   ---
   int nleq_numjac(NLEQ_FFUN *f, int n, double *x, double *fx,
                   double *scale, double *jac, int *nfcn,
                   struct NLEQ_OPT *opt)
   This function should compute the Jacobian by numerical difference
   approximation. The main codes will also work if this function just
   does nothing more than returning a nonzero value, but for this case
   the Jacobian routine must be always supplied by the user.
   The parameters are:
   NLEQ_FFUN *f       (input)   A pointer to the user problem function f.
   int       n        (input)   The number of equations and unknowns
                                and size of any vectors.
   double    *x       (input)   The input vector where to compute the Jacobian.
                                The content of x must be preserved on output!
   double    *fx      (input)   The vector f(x).
                                The content of fx must be preserved on output!
   double    *scale   (input)   The scaling vector.
                                The content of scale must be preserved on
                                output!
   double   *jac      (output)  A pointer to the Jacobian storage.
   int      *nfcn     (in/out)  The value of *nfcn must be incremented by one
                                for each call of the user problem function.
   struct NLEQ_OPT *opt (input) The input options of the main routine.

   A successful completion must be indicated by a return code 0, a failure
   condition by any nonzero return code.
   ---
   void nleq_jaccolumn_scale(int n, double *jac, double *scale,
                             struct NLEQ_OPT *opt)
   This routine must scale the columns of the Jacobian and change the sign
   of all Jacobian elements. The funtion is used by the main codes
   nleq_err and qnerr.
   The scaling must accomplish the following transformation:

   jac[i][j] = - jac[i][j]*scale[j] for j=0, ..., n-1  and i=0, ..., n-1.

   The parameters are:
   int       n        (input)   The number of equations and unknowns
                                and size of any vectors.
   double   *jac      (in/out)  A pointer to the Jacobian storage.
   double    *scale   (input)   The scaling vector.
                                The content of scale must be preserved on
                                output!
   struct NLEQ_OPT *opt (input) The input options of the main routine.
 */

#ifdef FMAT
#define JMODE		'N'
#else
#define JMODE		'T'
#endif

#define AJDEL		1.0e-8
#define AJMIN		1.0e-4

static int 		*pivot=NULL;
static double	 	*mat=NULL;
static double		*v=NULL;


/*==========================================================================*/

static int nleq_jacalloc(int n, double **jac, int *ldjac, struct NLEQ_OPT *opt)
{
  *ldjac = n;
  return nleq_fwalloc(n*n, jac, "jac");
}


/*==========================================================================*/

static void nleq_jacrow_scale(int n, double *jac, double *scale, struct NLEQ_OPT *opt)
{
  int		i,j;
  double	s;

  // scale jacobian and change sign of it
  for (i=0;i<n;i++)
    {
      s = -scale[i];
#ifndef FMAT
      for (j=i*n;j<(i+1)*n;j++)
	jac[j] /= s;
#else
      for (j=i;j<n*n;j+=n)
	jac[j] /= s;
#endif
    }

  return;
}


/*==========================================================================*/

static int nleq_linfact(int n, double *jac, struct NLEQ_OPT *opt)
{
  int		fail=0;

  if (!pivot)
    {
      fail=nleq_iwalloc(n, &pivot, "pivot");
      if ( fail !=0 ) return fail;
    }
  mat=jac;
  dgetrf_(&n, &n, mat, &n, pivot, &fail);
  return fail;
}


/*==========================================================================*/

static void nleq_linfree()
{
  if(pivot) free(pivot);
  pivot = NULL;
  if (v)    free(v);
  v     = NULL;
  if (mat)  free(mat);
  mat   = NULL;
  return;
}


/*==========================================================================*/

static int nleq_linsol(int n, double *b, struct NLEQ_OPT *opt)
{
  int		fail=0, nrhs=1;
  char		mode=JMODE;

  if (mat)
    dgetrs_(&mode, &n, &nrhs, mat, &n, pivot, b, &n, &fail);
  else
    fail = -980;

  return fail;
}


/*==========================================================================*/

static int nleq_numjac(NLEQ_FFUN *f, int n, double *x, double *fx, double *scale, double *jac, int *nfcn, struct NLEQ_OPT *opt)
{
  int		i, k;
  int		fail=0;
  double	w, u;

  if (!v) fail=nleq_fwalloc(n, &v, "v");
  if (fail!=0) return -992;
  for (k=0;k<n;k++)
    {
      w = x[k];
      if (scale)
	u = MAX( fabs(x[k]) , scale[k] );
      else
	u = fabs(x[k]);
      u = MAX(u, AJMIN)*AJDEL*SIGN(x[k]);
      x[k] = w+u;
      f(&n, x, v, &fail);
      (*nfcn)++;
      if ( fail != 0 ) return fail;
      x[k] = w;
#ifndef FMAT
      for (i=0;i<n;i++) jac[i*n+k] = (v[i]-fx[i]) / u ;
#else
      kn = k*n;
      for (i=0;i<n;i++) jac[kn+i] = (v[i]-fx[i]) / u ;
#endif
    }
  return fail;
}


/*============================================================================
 *
 *  Common utility routines for the NewtonLib package.
 *
 *  Routines originally defined in utils.c
 *
 *============================================================================
 *
 *  Written by        L. Weimann
 *  Purpose           Performing certain common tasks of NewtonLib codes
 *  Category          ???. - Utilities
 *  Keywords          Memory allocation, scaled norm, scaled scalarproduct,
                      data output
 *  Version           1.1
 *  Revision          May 2006
 *  Latest Change     May 2006
 *  Library           NewtonLib
 *  Code              C, Double Precision
 *  Environment       Standard C environment on PC's,
                      workstations and hosts.
 *  Copyright     (c) Konrad-Zuse-Zentrum fuer
                      Informationstechnik Berlin (ZIB)
                      Takustrasse 7, D-14195 Berlin-Dahlem
                      phone : + 49/30/84185-0
                      fax   : + 49/30/84185-125
 *  Contact           Lutz Weimann
                      ZIB, Division Scientific Computing,
                           Department Numerical Analysis and Modelling
                      phone : + 49/30/84185-185
                      fax   : + 49/30/84185-107
                      e-mail: weimann@zib.de

   ---------------------------------------------------------------

 * Licence
     You may use or modify this code for your own non commercial
     purposes for an unlimited time.
     In any case you should not deliver this code without a special
     permission of ZIB.
     In case you intend to use the code commercially, we oblige you
     to sign an according licence agreement with ZIB.

 * Warranty
     This code has been tested up to a certain level. Defects and
     weaknesses, which may be included in the code, do not establish
     any warranties by ZIB. ZIB does not take over any liabilities
     which may follow from acquisition or application of this code.

 * Software status
     This code is under care of ZIB and belongs to ZIB software class 2.

      ------------------------------------------------------------


    This file contains the following routines and functions:
    --------------------------------------------------------

    nleq_fwalloc - allocate memory for a double precision array
    nleq_iwalloc - allocate memory for an integer array
    nleq_pfwalloc - allocate memory for an array of pointers to
                    double precision arrays
    nleq_scaled_norm2 - compute the scaled Euclidian-norm of a vector
    nleq_scaled_sprod - compute a scaled scalar product
    nleq_norm2 - compute the unscaled Euclidian-norm of a vector
    nleq_descale - compute the unscaled vector from a scaled vector
    nleq_dataout - write data output
    nleq_parcheck_and_print - check and print parameter and options settings

    The calls of the routines/functions are as follows:
    ---------------------------------------------------

    int nleq_fwalloc(int size, double **ptr, char vname[])

    This function allocates memory for a double array via the malloc function.
    The parameters are:
    int    size    (input)  The number of elements of type double to be allocated.
    double **ptr   (output) The pointer to memory, of type (double *), which has
                            been returned by the malloc function.
    char   vname[] (input)  An identifying character string of the memory portion
                            to be allocated, used for print within a possible error
                            message.
    The int return code is 0, if the memory allocation was successful, or -999
    if the allocation failed.
    ---
    int nleq_iwalloc(int size, int **ptr, char vname[])

    This function allocates memory for a int array via the malloc function.
    The parameters are:
    int    size    (input)  The number of elements of type int to be allocated.
    int    **ptr   (output) The pointer to memory, of type (int *), which has
                            been returned by the malloc function.
    char   vname[] (input)  An identifying character string of the memory portion
                            to be allocated, used for print within a possible error
                            message.
    The int return code is 0, if the memory allocation was successful, or -998
    if the allocation failed.
    ---
    int nleq_pfwalloc(int size, double ***ptr, char vname[])

    This function allocates memory for a pointer array via the malloc function.
    The parameters are:
    int    size    (input)  The number of elements of type pointer to be allocated.
    int    ***ptr  (output) The pointer to memory, of type (double **) (i.e.
                            to pointers which are pointing to memory allocated
                            for double arrays), which has been returned by the
                            malloc function.
    char   vname[] (input)  An identifying character string of the memory portion
                            to be allocated, used for print within a possible error
                            message.
    The int return code is 0, if the memory allocation was successful, or -997
    if the allocation failed.
    ---
    Note: In order to activate some debug output on dynamic memory allocation
          set the C preprocessor flag DEBUG, i.e. set the option -DDEBUG if
          you use the GNU C-compiler (gcc).
    ---
    double nleq_scaled_norm2(int n, double *v, double *scale)

    This function computes the scaled norm of the (double) vector v of
    size n, using the (double) vector scale for scaling, as below:
    result := Sqrt ( ( Sum(i=0 to n-1) (v[i]/scale[i])^2 ) / n ) .
    ---
    double nleq_scaled_sprod(int n, double *v1, double *v2, double *scale)

    This function computes the scaled scalar product of the (double) vectors
    v1 and v2 of size n, using the (double) vector scale for scaling, as below:
    result := ( Sum(i=0 to n-1) (v1[i]/scale[i])*(v2[i]/scale[i]) ) / n .
    ---
    double nleq_norm2(int n, double *v)

    This function computes the (ordinary) norm of the (double) vector v of
    size n, as below:
    result := Sqrt ( ( Sum(i=0 to n-1) v[i]^2 ) / n ) .
    ---
    void nleq_descale(int n, double *v1, double *v2, double *scale)

    This routine computes the descaled vector of the vector v1 of size n and
    stores the result to the vector v2, as below:
    v2[i] = v1[i]*scale[i]  for i=0, ..., n-1 .
    ---
    void nleq_dataout(int k, int n, double *x, struct NLEQ_DATA *data)

    This routine is designed to print out computed data within each iteration
    step. It may be replaced or extended for special purposes.
    The parameters are (all input arguments):

    int    k    The current iteration step number.
    int    n    The size of any double arrays mentioned below.
    double *x   An array which holds the current iterate x^k .

    The fields of the struct NLEQ_DATA are:

    double* data->fx       An array which holds fun(x^k).
    double* data->dx       An array which holds the latest Newton correction.
    double  data->lambda   The current damping factor.
    double  data->normdx   The unscaled norm of the latest  Newton correction.
    double  data->normf    The scaled norm of fun(x^k).
    enum    data->mode     The mode of the current iterate:
                           Initial (=1): This is the first call of nleq_dataout.
                           Intermediate (=2): This is an intermediate call of
                                              nleq_dataout.
                           Solution (=3): This is a final call of nleq_dataout,
                                          and *x holds an approximate solution.
                           Final (=4) : This is a final call of nleq_dataout,
                                        but *x does not hold a solution.
    ---
    int nleq_parcheck_and_print(int n, struct NLEQ_OPT *opt,
                                struct NLEQ_FUN fun, int nleq_code)

    This function checks and prints parameter and options settings.
    The parameters are:

    int             n     The number of equations and unknowns.
    struct NLEQ_OPT *opt  A pointer to an options data structure. See routines
                          qnres, nleq_res, qnerr and nleq_err for details on it.
    struct NLEQ_FUN fun   A structure, holding pointers to the problem function
                          routine and (optionally) the Jacobian routine.
                          For more details, see the description in one main
                          routine e.g. nleq_res.
    int        nleq_code  The identification number of the calling code.
                          Valid identifications are:
                          0 : qnres
                          1 : nleq_res
                          2 : qnerr
                          3 : nleq_err

 */

/*==========================================================================*/

static void nleq_dataout(int k, int n, double *x, struct NLEQ_DATA *data)
{
  int			i;
  double		*fx = data->fx;
  if (FITER)
    {
      fprintf(FITER, "%5i", k);
      for (i=0;i<n;i++)
	fprintf(FITER, "  % 14.10e", x[i]);
      fprintf(FITER, "\n");
    }
  if (FRES)
    {
      fprintf(FRES, "%5i", k);
      for (i=0;i<n;i++)
	fprintf(FRES, "  % 14.10e", fx[i]);
      fprintf(FRES, "\n");
    }
  if (FMISC)
    fprintf(FMISC, "%5i  %1i  % 14.10e  % 14.10e  % 14.10e  %10.8f  %9.3e\n",
	k, data->codeid, data->normf, data->normdx, data->normdxbar, data->lambda, data->theta);

  if ( DATALEVEL==0 ) return;
  if ( k==0 )
    {
      fprintf(DATA, "  Start data:\n  N = %i\n\n", n);
      fprintf(DATA, "  Format: iteration-number, (x(i), i=1, ...N), Normf , Normdx\n\n");
      fprintf(DATA, "    Initial data:\n\n");
    }
  else if ( data->mode==Solution )
    fprintf(DATA, "\n    Solution data:\n\n");
  else if ( data->mode==Final )
    fprintf(DATA, "\n    Final data:\n\n");
  else if ( k==1 && DATALEVEL>1 )
    fprintf(DATA, "\n    Intermediate data:\n\n");
  if ( k==0 || data->mode==Solution || data->mode==Final || DATALEVEL > 1 )
    {
      fprintf(DATA, "  %4i\n", k);
      for (i=0;i<n-2;i+=3)
	fprintf(DATA, "             % 14.10e  % 14.10e  % 14.10e\n", x[i], x[i+1], x[i+2]);
      if (i<n)
	{
	  fprintf(DATA, "             % 14.10e", x[i]); i++;
	  if (i<n)
	    fprintf(DATA, "  % 14.10e", x[i]);
	  fprintf(DATA, "\n");
	}
      fprintf(DATA, "             % 14.10e  % 14.10e\n", data->normf, data->normdx);
    }
  return;
}


/*==========================================================================*/

int nleq_fwalloc(int size, double **ptr, char vname[])
{
  int			i;
#ifdef DEBUG
  fprintf(stdout, "\n allocation of %s follows, size=%i\n", vname, size);
#endif
  *ptr = (double *) malloc((int) size*sizeof(double)) ;
#ifdef DEBUG
  fprintf(stdout, "\n allocation of %s done\n", vname);
#endif
  if (*ptr == NULL)
    {
      if (ERRORLEVEL>0 && ERROR)
	fprintf(ERROR, "\n allocation of %s failed!\n", vname);
      return -999;
    }
  else
    {
      for(i=0;i<size;i++)
	(*ptr)[i]=0.0;
      return 0;
    }
}


/*==========================================================================*/

static int nleq_iwalloc(int size, int **ptr, char vname[])
{
  int			i;
#ifdef DEBUG
  fprintf(stdout, "\n allocation of %s follows, size=%i\n", vname, size);
#endif
  *ptr = (int *) malloc((int) size*sizeof(int)) ;
#ifdef DEBUG
  fprintf(stdout, "\n allocation of %s done\n", vname);
#endif
  if (*ptr == NULL)
    {
      if (ERRORLEVEL>0 && ERROR)
	fprintf(ERROR, "\n allocation of %s failed!\n", vname);
      return -998;
    }
  else
    {
      for(i=0;i<size;i++) (*ptr)[i]=0;
      return 0;
    }
}


/*==========================================================================*/

static double nleq_norm2(int n, double *v)
{
  int			i;
  double		rval = 0.0;

  for (i=0;i<n;i++)
    rval += v[i]*v[i];

  return sqrt( rval / (double)n );
}


/*==========================================================================*/

#define TOLMIN 1.0e-15
#define TOLMAX 1.0e-1

static int nleq_parcheck_and_print(int n, struct NLEQ_OPT *opt, struct NLEQ_FUN fun, int nleq_code)
{
  int			i;
  LOGICAL		restricted = opt->restricted;
  double		large = 1.0/SMALL, default_scale;

  if (!fun.fun)
    {
      if ( ERRORLEVEL>0 )
	fprintf(ERROR, "\n Error - Problem function fun.fun has not been defined\n");
      return -99;
    }
  if ( n<=0 )
    {
      if ( ERRORLEVEL>0 )
	fprintf(ERROR, "\n Error - Number of Equations/unknowns must be >0\n");
      return 20;
    }
  if ( opt->tol <= 0 )
    {
      if ( ERRORLEVEL>0 )
	fprintf(ERROR, "\n Error - opt->tol must be positive\n");
      return 21;
    }
  else
    {
      if ( opt->tol > TOLMAX )
	{
	  opt->tol = TOLMAX;
	  if ( ERRORLEVEL>1 )
	    fprintf(ERROR, "\n User prescribed RTOL decreased to reasonable largest value RTOL=%e\n", opt->tol);
	}
      else if ( opt->tol < TOLMIN )
	{
	  opt->tol = TOLMIN;
	  if ( ERRORLEVEL>1 )
	    fprintf(ERROR, "\n User prescribed RTOL increased to reasonable smallest value RTOL=%e\n", opt->tol);
	}
    }
  if ( nleq_code==3 )
    {
      if ( opt->nonlin >= Highly_Nonlinear )
	default_scale = opt->tol;
      else
	default_scale = 1.0;
    }
  else
    default_scale = 1.0;
  if (opt->scale)
    {
      for (i=0;i<n;i++)
	{
	  if ( (opt->scale)[i] < 0.0 )
	    {
	      if ( ERRORLEVEL>0 )
		fprintf(ERROR, "\n Error - negative value (opt->scale)[%i] supplied\n", i);
	      return 22;
	    }
	  else if ( (opt->scale)[i] == 0.0 )
	    (opt->scale)[i] = default_scale;
	  else if ( (opt->scale)[i] < SMALL )
	    {
	      if ( ERRORLEVEL>1 )
		fprintf(ERROR, "\n Warning (opt->scale)[%i] too small - increased to %e\n", i, SMALL);
	      (opt->scale)[i] = SMALL;
	    }
	  else if ( (opt->scale)[i] > large )
	    {
	      if ( ERRORLEVEL>1 )
		fprintf(ERROR, "\n Warning (opt->scale)[%i] too large - decreased to %e\n", i, large);
	      (opt->scale)[i] = large;
	    }
	}
    }

  if ( MONITORLEVEL==0 ) return 0;

  fprintf(MONITOR, "\n Problem dimension: n = %i\n", n);
  if ( nleq_code <= 1 )
    fprintf(MONITOR, "\n Prescribed residuum threshold: %e\n", opt->tol);
  else if ( nleq_code <= 3 )
    fprintf(MONITOR, "\n Prescribed relative precision: %e\n", opt->tol);
  fprintf(MONITOR, "\n The Jacobian is supplied by ");
  if (fun.jac)
    fprintf(MONITOR, "a user routine\n\n");
  else
    fprintf(MONITOR, "numerical differentiation\n\n");
  if ( nleq_code == 1 || nleq_code == 3 )
    {
      fprintf(MONITOR, " The problem is specified as being ");
      switch ( opt->nonlin )
      {
	case Mildly_Nonlinear:
	  fprintf(MONITOR, "mildly nonlinear\n");
	  break;
	case Highly_Nonlinear:
	  fprintf(MONITOR, "highly nonlinear\n");
	  break;
	case Extremely_Nonlinear:
	  fprintf(MONITOR, "extremely nonlinear\n");
	  restricted = True;
      }
      if ( restricted )
	fprintf(MONITOR, " The restricted monotonicity test will be applied\n");
      else
	fprintf(MONITOR, " The standard monotonicity test will be applied\n");
    }
  fprintf(MONITOR, " The maximum permitted number of iteration steps is: %i\n", opt->maxiter);

  return 0;
}


/*==========================================================================*/

static int nleq_pfwalloc(int size, double ***ptr, char vname[])
{
  int i;
#ifdef DEBUG
  fprintf(stdout, "\n allocation of %s follows, size=%i\n", vname, size);
#endif
  *ptr = (double **) malloc((int) size*sizeof(double*)) ;
#ifdef DEBUG
  fprintf(stdout, "\n allocation of %s done\n", vname);
#endif
  if (*ptr == NULL)
    {
      if (ERRORLEVEL>0 && ERROR)
	fprintf(ERROR, "\n allocation of %s failed!\n", vname);
      return -997;
    }
  else
    {
      for(i=0;i<size;i++) (*ptr)[i]=NULL;
      return 0;
    }
}


/*==========================================================================*/

static void nleq_scale(int n, double *v1, double *v2, double *scale)
{
  int			i;

  for (i=0;i<n;i++)
    v2[i]=v1[i]/scale[i];

  return;
}


/*==========================================================================*/

static double nleq_scaled_norm2(int n, double *v, double *scale)
{
  int			i;
  double		t, rval = 0.0;

  for (i=0;i<n;i++)
    {
      t=v[i]/scale[i];
      rval += t*t;
    }

  return sqrt( rval / (double)n );
}


/*==========================================================================*/

static double nleq_scaled_sprod(int n, double *v1, double *v2, double *scale)
{
  int			i;
  double		t1, t2, rval = 0.0;
  for (i=0;i<n;i++)
    {
      t1=v1[i]/scale[i];
      t2=v2[i]/scale[i];
      rval += t1*t2;
    }

  return rval / (double)n ;
}


/*==========================================================================*/
