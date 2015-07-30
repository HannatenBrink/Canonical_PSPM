/***
   NAME
     cvode
   DESCRIPTION
     Interface to sundials' CVODE solver of an ODE system

   Last modification: AMdR - Jan 09, 2014
***/
#include "stdio.h"
#include "string.h"

// Includes for the sundials integrator
#include <cvode/cvode.h>             // prototypes for CVODE fcts., consts.
#include <nvector/nvector_serial.h>  // serial N_Vector types, fcts., macros

// Type definitions for the callback rooutines
typedef void		RHSFUN(double, double *, double *);
typedef double		STOPFUN(double, double *);

static RHSFUN		*UserRhs;
static STOPFUN		*UserStop;


// Default integration settings

#ifndef INIT_H
#define INIT_H		0.001
#endif
#ifndef LARGEST_STEP
#define LARGEST_STEP	1.0
#endif
#ifndef SMALLEST_STEP
#define SMALLEST_STEP	1.0E-8
#endif
#ifndef ABS_ERR
#define ABS_ERR		1.0E-13
#endif
#ifndef ACCURACY
#define ACCURACY	1.0E-9
#endif
#ifndef FUNCTOL
#define FUNCTOL		1.0E-8
#endif
#ifndef PROFILE
#define PROFILE		0
#endif
#ifndef MXSTEPS
#define MXSTEPS		10000
#endif

#ifndef SUCCES
#define SUCCES			0
#endif
#ifndef FAILURE
#define FAILURE			1
#endif

/*
 *===========================================================================
 *
 * STATIC CALLBACK FUNCTIONS USED INTERNALLY IN THE MODULE
 *
 *===========================================================================
 */

static int RhsCallBack(realtype x, N_Vector y, N_Vector ydot, void *user_data)
{
  double *ydata    = (double *)NV_DATA_S(y);
  double *ydotdata = (double *)NV_DATA_S(ydot);

  (*UserRhs)((double)x, ydata, ydotdata);

  return 0;
}


/*==============================================================================*/

static int 	StopCallBack(realtype x, N_Vector y, realtype *stop, void *user_data)
{
  double *ydata    = (double *)NV_DATA_S(y);

  *stop = (realtype)(*UserStop)((double)x, ydata);

  return 0;
}


/*
 *===========================================================================
 *
 * UTILITY ROUTINES FROM SUNDIALS EXAMPLE PROGRAMS
 *
 *===========================================================================
 */

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer
 */

static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;

  // Check if SUNDIALS function returned NULL pointer - no memory allocated
  if (opt == 0 && flagvalue == NULL)
    {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	  funcname);
      return(1);
    }

  // Check if flag < 0
  else if (opt == 1)
    {
      errflag = (int *) flagvalue;
      if (*errflag < 0)
	{
	  fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
	  return(1);
	}
    }

  // Check if function returned NULL pointer - no memory allocated
  else if (opt == 2 && flagvalue == NULL)
    {
      fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	  funcname);
      return(1);
    }

  return(0);
}


/*===========================================================================*/

static void PrintFinalStats(void *cvode_mem)
{
  long int nst, nfe, nsetups, nni, ncfn, netf, nge;
  int flag;

  flag = CVodeGetNumSteps(cvode_mem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

  flag = CVodeGetNumGEvals(cvode_mem, &nge);
  check_flag(&flag, "CVodeGetNumGEvals", 1);

  fprintf(stderr, "\nFinal Statistics:\n");
  fprintf(stderr, "nst = %-6ld nfe  = %-6ld nsetups = %-6ld nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n \n",
	  nst, nfe, nsetups, nni, ncfn, netf, nge);
}


/*
 *===========================================================================
 *
 * INTERFACE FUNCTIONS TO BE EXPOSED TO THE CALLING PROGRAM
 *
 *===========================================================================
 */


int		odesolve(double *yvec, int vecdim,
			 double *xinit, double xmax,
			 void (*rhs)(double, double *, double *),
			 double (*stop)(double, double *))

/*
 * odesolve - Solve the system of ODEs specified in the function (*rhs)() up to
 * 	      a maximum time xmax or until the function (*stop)() returns a 0
 *	      result
 */
  
{
  int			flag;
  N_Vector 		y = NULL;
  void			*cvode_mem = NULL;

  UserRhs  = rhs;
  UserStop = stop;
  // Call CVodeCreate to create the solver memory and specify the
  // standard settings for non-stiff problems (CV_ADAMS, CV_FUNCTIONAL)
  // See 4.5.1. in user-guide 2.7.0
  cvode_mem = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);
  if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return FAILURE;

  // Create serial vector of length vecdim for I.C.
  y = N_VMake_Serial(vecdim, yvec);
  if (check_flag((void *)y, "N_VMake_Serial", 0))
    {
      CVodeFree(&cvode_mem);
      return FAILURE;
    }

  // Call CVodeInit to initialize the integrator memory and specify the
  // user's right hand side function in y'=f(t,y), the initial time xval, and
  // the initial dependent variable vector y.
  flag = CVodeInit(cvode_mem, RhsCallBack, *xinit, y);
  if (check_flag(&flag, "CVodeInit", 1))
    {
      N_VDestroy_Serial(y);
      CVodeFree(&cvode_mem);
      return FAILURE;
    }

  // Call CVodeSStolerances to specify the scalar relative tolerance
  // and absolute tolerances
  flag = CVodeSStolerances(cvode_mem, ACCURACY, ABS_ERR);
  if (check_flag(&flag, "CVodeSStolerances", 1))
    {
      // Dispose of previous allocated memory
      N_VDestroy_Serial(y);
      CVodeFree(&cvode_mem);
      return FAILURE;
    }

  // Call CVodeRootInit to specify the root function StopCallBack with 1 component
  if (UserStop)
    {
      // Check only crossings from negative to positive
      int rootdir[] = {1};
      flag = CVodeRootInit(cvode_mem, 1, StopCallBack);
      if (check_flag(&flag, "CVodeRootInit", 1))
        {
          // Dispose of previous allocated memory
          N_VDestroy_Serial(y);
          CVodeFree(&cvode_mem);
          return FAILURE;
        }
      flag = CVodeSetRootDirection(cvode_mem, rootdir);
      if (check_flag(&flag, "CVodeSetRootDirection", 1))
        {
          // Dispose of previous allocated memory
          N_VDestroy_Serial(y);
          CVodeFree(&cvode_mem);
          return FAILURE;
        }
    }

  // Call CVodeSetInitStep to set the initial step size
  flag = CVodeSetInitStep(cvode_mem, INIT_H);
  if (check_flag(&flag, "CVodeSetInitStep", 1))
    {
      // Dispose of previous allocated memory
      N_VDestroy_Serial(y);
      CVodeFree(&cvode_mem);
      return FAILURE;
    }

  // Call CVodeSetMinStep to set the minimum step size
  flag = CVodeSetMinStep(cvode_mem, SMALLEST_STEP);
  if (check_flag(&flag, "CVodeSetMinStep", 1))
    {
      // Dispose of previous allocated memory
      N_VDestroy_Serial(y);
      CVodeFree(&cvode_mem);
      return FAILURE;
    }

  // Call CVodeSetMaxStep to set the maximum step size
  flag = CVodeSetMaxStep(cvode_mem, LARGEST_STEP);
  if (check_flag(&flag, "CVodeSetMaxStep", 1))
    {
      // Dispose of previous allocated memory
      N_VDestroy_Serial(y);
      CVodeFree(&cvode_mem);
      return FAILURE;
    }

  // Call CVodeSetMaxStep to set the maximum step size
  flag = CVodeSetMaxNumSteps(cvode_mem, MXSTEPS);
  if (check_flag(&flag, "CVodeSetMaxNumSteps", 1))
    {
      // Dispose of previous allocated memory
      N_VDestroy_Serial(y);
      CVodeFree(&cvode_mem);
      return FAILURE;
    }

  // Now do the actual integration
  
  flag = CVode(cvode_mem, xmax, y, xinit, CV_NORMAL);
  if (check_flag(&flag, "CVode", 1))
    {
      // Dispose of previous allocated memory
      N_VDestroy_Serial(y);
      CVodeFree(&cvode_mem);
      return FAILURE;
    }
  // Save the final dependent variable values
  (void)memcpy(yvec, NV_DATA_S(y), vecdim*sizeof(double));

  if (PROFILE) PrintFinalStats(cvode_mem);

  // Dispose of previous allocated memory
  N_VDestroy_Serial(y);
  CVodeFree(&cvode_mem);

  return SUCCES;
}

/*===========================================================================*/
