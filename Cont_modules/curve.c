/***
   NAME
     curve
   DESCRIPTION
     This module implements routines that are specifically used in locating
     points on the equilibrium branch of a structured population model.

     Last modification: AMdR - Mar 07, 2014
***/
#ifndef CURVE
#define CURVE
#endif
#include <globals.h>
#include "nleq.h"

int			ipiv[MAX_PNTDIM];
double			Jac[MAX_PNTDIM_SQ];
double			JacCopy[MAX_PNTDIM_SQ];
double			oldvec[MAX_PNTDIM];
double			work[4*MAX_PNTDIM];
int			lwork = 4*MAX_PNTDIM;



#if (MOOREPENROSE == 1)
/*===========================================================================*/

int	FindPoint(const int pntdim, double *guess, double *tanvec,
		  double ytol, double rhstol, const int max_iter, int (*fnc)(double *, double *))

  /*
   * FindPoint - Routine locates a point on a curve determined by a
   *		 system of non-linear, algebraic equations.
   *		 Moore-Penrose continuation is applied to locate the point
   *		 and compute at the same time the tangent vector to the curve.
   *		 For details, see page 15 in the file MoorePenroseContinuation.pdf,
   *		 entitled "Lecture 1: Continuation problems. Numerical continuation of
   *		 equilibria and limit cycles of ODEs. Yu.A. Kuznetsov (Utrecht University, NL)
   *		 May 4, 2009"
   *
   * Arguments - pntdim	   : The dimension of the solution point on the curve.
   *			     The dimension of the system of equations is
   *			     assumed to be exactly 1 less.
   *		 guess	   : Pointer to an array containing the initial point
   *			     to start the iteration from. The first element of
   *			     the vector is assumed to be non-adjustable parameter.
   *		 ytol	   : Tolerance determining when change in y equals zero.
   *		 rhstol	   : Tolerance determining when RHS equals zero.
   *		 max_iter  : Maximum stepnumber allowed in iteration.
   *		 fnc	   : Pointer to function specifying the system of
   *			     equations. The function must have a (double)
   *			     pointer as first argument, containing the point
   *			     in which to evaluate the system and a (double)
   *			     pointer as second argument, containing the
   *			     results after evaluation of the equations.
   */

{
  register int		iter, i;
  int			sysdim  = pntdim-1;
  int			pntdim1 = pntdim;
  int			pntdim2 = pntdim*pntdim;
  int			nrhs = 1;
  int			info;
  static int		fast_iters = 0, slow_iters = 0, first = 1;
  char			trans[1] = {"N"};
  double		ynorm, dynorm, rhsnorm, vnorm;
  double		y[pntdim];
  double		tv[pntdim];
  double		dy[pntdim];
  double		dv[pntdim];
  double		B[pntdim*pntdim];
  double		BLU[pntdim*pntdim];
  double		Q[pntdim];
  double		R[pntdim];

  if (first)
    {
      ReportMsg(2, "\n\nUsing LU decomposition for solving linear systems\n");
    }
  first = 0;

  COPY(pntdim, guess, 1, y, 1);  //pndtdim = array size, guess=input, 1=input stride, dy=output, 1=output stride. copy of guess naar vec y
  COPY(pntdim, tanvec, 1, tv, 1);  //copy of tanvec naar vec tv (tanvec bestaat in eerste instantie uit 1 en verder nullen
  memset((void *)dy,  0, pntdim*sizeof(double));  //vector dy wordt gevuld met nullen

  ReportMsg(2, "\nNew point:\tParameter: %10.5E\t", y[0]);  //eerste parameter wordt geprint
  ReportMsg(2, "Dimension: %4d\n", pntdim);

  // The iteration loop
  for(iter=0; iter<max_iter; iter++)
    {
      // Compute norm of Y and of dY
      ynorm  = NRM2(pntdim, y, 1);   //returns Euclidian norm of vec y (guess) (pntdim=dimensie vector, y=vector, 1=increments between elements of y (dus 1 is alles)
      dynorm = NRM2(pntdim, dy, 1);  //returns euclidian norm of vec dy (vol met nullen)
      dynorm = dynorm/(1.0+dynorm);
      if (!issane(ynorm) || !issane(dynorm))  //checkt of ze normaal zijn (dus niet nan, unknown, inf of subnormal(?))
	{
	  ErrorMsg(__FILE__, __LINE__, "Norm overflow in FindPoint");
	  return NORM_OVERFLOW;
	}

      checkESCpressed();
      memset((void *)Q,  0, pntdim*sizeof(double));  //vector Q wordt gevuld met nullen
      if ((*fnc)(y, Q) == FAILURE)     ///in de equation wordt de guess gezet, plus een lege vector Q die gevuld wordt met de resultaten
	{
	  ErrorMsg(__FILE__, __LINE__, "Right-hand side computation failed");
	  return FAILED_EVALUATION;
	}

      // Compute rhsnorm
      rhsnorm = NRM2(sysdim, Q, 1); ///de Euclidian norm of vec Q (res)
      rhsnorm = rhsnorm/(1.0+rhsnorm);
      ReportMsg(1, "\tdY  norm: %10.5E\tRHS norm: %10.5E\n", dynorm, rhsnorm); //dynorm and rhsnorm wordt geschreven in command line (dy is altijd eerst nul)

      // Return if converged or diverged
      if ((!issane(rhsnorm)) || (rhsnorm > RHSMAX))
	{
	  ErrorMsg(__FILE__, __LINE__, "Norm overflow in FindPoint");
	  return NORM_OVERFLOW;
	}
      else if ((rhsnorm < rhstol) && (dynorm < ytol)) //als de normen klein genoeg zijn stop de iteratie
      {
	  COPY(pntdim, y, 1, guess, 1);   ///copy of y naar guess
	  COPY(pntdim, tv, 1, tanvec, 1); ///copy of tv naar tanvec
	  if ((Stepchange) && (Stepreduce == 1))
	    {
	      if (iter < STEP_DOUBLE)
		{
		  fast_iters++;
		  slow_iters = 0;
		  if (fast_iters == 2)
		    {
		      parstep   *= 1.5;
		      fast_iters = 0;
		    }
		}
	      else if (iter > STEP_HALF)
		{
		  slow_iters++;
		  fast_iters = 0;
		  if (slow_iters == 2)
		    {
		      parstep   *= 0.5;
		      slow_iters = 0;
		    }
		}
	      else
		{
		  fast_iters = 0;
		  slow_iters = 0;
		}
	    }

	  return SUCCES;
      }
      // Compute Jacobian every UPDATE_JAC steps, otherwise the Jacobian is updated
      // via a Broyden update
      if (!(iter%UPDATE_JAC))
	{
	  ReportMsg(2, "%-s", "Computing jacobian");
	  memset((void *)Jac, 0, MAX_PNTDIM_SQ*sizeof(double));
	  Jacobian(pntdim, y, sysdim, Jac, fnc, CENTRALDIFF);		// Compute J = F_x(X^k)
	  ReportMsg(2, ".....Ok!\n");
	}
      else								// Broyden update of Jacobian
	{
	  dynorm = DOT(pntdim, dy, 1, dy, 1);
	  GER(CblasColMajor, pntdim, pntdim, 1.0/dynorm, Q, 1, dy, 1, Jac, pntdim);
	}
      //  		       | J     |
      // Create the matrix B = | V^k^T |
      for (i=0; i<pntdim; i++)
	{
	  COPY(sysdim, Jac+i*sysdim, 1, B+i*pntdim, 1);
	  *(B+i*pntdim+sysdim) = tv[i];
	}

      /*
       * Notice that the Jacobian is stored as
       *
       *				|dF1/dy1 ... dFn/dy1|
       *				|dF1/dy2 ... dFn/dy2|
       *			J =	|   .		.   |
       *				|   .		.   |
       *				|dF1/dyn ... dFn/dyn|
       *
       * Meaning that all coefficients pertaining to yi are to be found in ROW i (as opposed to column i).
       * The matrix is hence stored in column-wise (fortran) style.
       * Solving J.dy = -F(y) with dy = (dy1 ... dyn) and F(y) = (F1(y) ... Fn(y)) requires the variable
       * trans[1] to be defined as {"N"} (see programs/various/testlapack.c for details).
       */

      // LU decompose the matrix B
      memset((void *)BLU, 0, pntdim2*sizeof(double));
      COPY(pntdim2, B, 1, BLU, 1);
      GETRF(&pntdim1, &pntdim1, BLU, &pntdim1, ipiv, &info);
      if (info < 0)
	{
	  ErrorMsg(__FILE__, __LINE__,
	      "Illegal value for parameter %d in GETRF", abs(info));
	  return ILLEGAL_INPUT;
	}
      else if (info > 0)
	{
	  ErrorMsg(__FILE__, __LINE__,
	      "Singular Jacobian matrix found in GETRF\n");
	  return SINGULARITY;
	}

      // Solve B (X^(k+1) - X^k) = - Q for dy = X^(k+1) - X^k
      COPY(pntdim, Q, 1, dy, 1);
      SCAL(pntdim, -1.0, dy, 1);
      GETRS(trans, &pntdim1, &nrhs, BLU, &pntdim1, ipiv, dy, &pntdim1, &info);
      if (info < 0)
	{
	  ErrorMsg(__FILE__, __LINE__,
		   "Illegal value for parameter %d in GETRS", abs(info));
	  return ILLEGAL_INPUT;
	}

      AXPY(pntdim, 1.0, dy, 1, y, 1);

      // Create the vector R = (J V^k 0)^T
      memset((void *)R, 0, pntdim*sizeof(double));
      GEMV(CblasColMajor, CblasNoTrans, sysdim, pntdim, 1.0, Jac, sysdim, tv, 1, 1.0, R, 1);
      R[sysdim] = 0;
      // Solve B (W - V^k) = - R for dv = W - V^k
      COPY(pntdim, R, 1, dv, 1);
      SCAL(pntdim, -1.0, dv, 1);
      GETRS(trans, &pntdim1, &nrhs, BLU, &pntdim1, ipiv, dv, &pntdim1, &info);
      if (info < 0)
	{
	  ErrorMsg(__FILE__, __LINE__,
		   "Illegal value for parameter %d in GETRS", abs(info));
	  return ILLEGAL_INPUT;
	}

      AXPY(pntdim, 1.0, dv, 1, tv, 1);
      vnorm = NRM2(pntdim, tv, 1);
      SCAL(pntdim,  1.0/vnorm, tv, 1);
    }

  return NO_CONVERGENCE;
}



#else  // (MOOREPENROSE != 1)
/*===========================================================================*/

int	FindPoint(const int pntdim, double *guess, double *tanvec,
		  double ytol, double rhstol, const int max_iter, int (*fnc)(double *, double *))

  /*
   * FindPoint - Routine locates a point on a curve determined by a
   *		 system of non-linear, algebraic equations.
   *	         The iteration adjusts the vector-elements following a simple
   *		 Newton-Chord method with Broyden update (see Kuznetsov pg. 418).
   *		 Pseudo-arclength continuation is used to continue past curve folds.
   *
   * Arguments - pntdim	   : The dimension of the solution point on the curve.
   *			     The dimension of the system of equations is
   *			     assumed to be exactly 1 less. 
   *		 guess	   : Pointer to an array containing the initial point
   *			     to start the iteration from. The first element of
   *			     the vector is assumed to be non-adjustable parameter.
   *		 ytol	   : Tolerance determining when change in y equals zero.
   *		 rhstol	   : Tolerance determining when RHS equals zero.
   *		 max_iter  : Maximum stepnumber allowed in iteration.
   *		 fnc	   : Pointer to function specifying the system of
   *			     equations. The function must have a (double)
   *			     pointer as first argument, containing the point
   *			     in which to evaluate the system and a (double)
   *			     pointer as second argument, containing the
   *			     results after evaluation of the equations.
   */

{
  register int		iter, i;
  int			sysdim  = pntdim-1;
  int			pntdim1 = pntdim;
  int			pntdim2 = pntdim*pntdim;
  int			nrhs = 1;
  int			info;
  static int		fast_iters = 0, slow_iters = 0, first = 1;
  char			trans[1] = {"N"};
  double		ynorm, dynorm, rhsnorm;
  double		y[pntdim];
  double		tv[pntdim];
  double		dy[pntdim];
  double		rhs[pntdim];

  if (first)
    {
      ReportMsg(2, "\n\nUsing QR factorization for solving linear systems\n");
    }
  first = 0;

  COPY(pntdim, guess, 1, y, 1); //copieert guess naar vec y
  memset((void *)dy,  0, pntdim*sizeof(double));

  ReportMsg(2, "\nNew point:\tParameter: %10.5E\t", y[0]);
  ReportMsg(2, "Dimension: %4d\n", pntdim);

  // The iteration loop
  for(iter=0; iter<max_iter; iter++)
    {
      // Compute norm of Y and of dY
      ynorm  = NRM2(pntdim, y, 1);
      dynorm = NRM2(pntdim, dy, 1);
      dynorm = dynorm/(1.0+dynorm);
      if (!issane(ynorm) || !issane(dynorm))
	{
	  ErrorMsg(__FILE__, __LINE__, "Norm overflow in FindPoint");
	  return NORM_OVERFLOW;
	}

      checkESCpressed();
      memset((void *)rhs,  0, pntdim*sizeof(double));
      if ((*fnc)(y, rhs) == FAILURE)
	{
	  ErrorMsg(__FILE__, __LINE__, "Right-hand side computation failed");
	  return FAILED_EVALUATION;
	}

      // Compute rhsnorm
      rhsnorm = NRM2(pntdim, rhs, 1);
      rhsnorm = rhsnorm/(1.0+rhsnorm);

      ReportMsg(1, "\tdY  norm: %10.5E\tRHS norm: %10.5E\n", dynorm, rhsnorm);

      // Return if converged or diverged
      if ((!issane(rhsnorm)) || (rhsnorm > RHSMAX))
	{
	  ErrorMsg(__FILE__, __LINE__, "Norm overflow in FindPoint");
	  return NORM_OVERFLOW;
	}
      else if ((rhsnorm < rhstol) && (dynorm < ytol))
	{
	  COPY(pntdim, y, 1, guess, 1); //copieert vec y naar vec guess
	  RhsNorm = rhsnorm;
	  if ((Stepchange) && (Stepreduce == 1))
	    {
	      if (iter < STEP_DOUBLE)
		{
		  fast_iters++;
		  slow_iters = 0;
		  if (fast_iters == 2)
		    {
		      parstep   *= 1.5;
		      fast_iters = 0;
		    }
		}
	      else if (iter > STEP_HALF)
		{
		  slow_iters++;
		  fast_iters = 0;
		  if (slow_iters == 2)
		    {
		      parstep   *= 0.5;
		      slow_iters = 0;
		    }
		}
	      else
		{
		  fast_iters = 0;
		  slow_iters = 0;
		}
	    }

	  return SUCCES;
	}

      // When not converged find new point via pseudo-arclength continuation
      COPY(pntdim, y, 1, tv, 1);  //copieert y naar vec tv
      AXPY(pntdim, -1.0, guess, 1, tv, 1);  /// tv = -1 * guess+tv
      rhs[sysdim] = DOT(pntdim, tv, 1, tanvec, 1);

      // Compute Jacobian every UPDATE_JAC steps, otherwise the Jacobian is updated
      // via a Broyden update (see below)
      if (!(iter%UPDATE_JAC))
	{
	  ReportMsg(2, "%-s", "Computing jacobian");
	  memset((void *)JacCopy, 0, MAX_PNTDIM_SQ*sizeof(double));
	  Jacobian(pntdim, y, sysdim, JacCopy, fnc, CENTRALDIFF);		// Compute J = F_x(X^k)
	  for (i=0; i<pntdim; i++)
	    {
	      COPY(sysdim, JacCopy+i*sysdim, 1, Jac+i*pntdim, 1);		// Extract dF/dx
	      *(Jac+i*pntdim+sysdim) = tanvec[i];
	    }
	  ReportMsg(2, ".....Ok!\n");
	}
      else					      // Broyden update of Jacobian
	{
	  dynorm = DOT(pntdim, dy, 1, dy, 1);
	  GER(CblasColMajor, pntdim, pntdim, 1.0/dynorm, rhs, 1, dy, 1, Jac, pntdim);
	}


      /*
       * Notice that the Jacobian is stored as 
       *
       *				|dF1/dy1 ... dFn/dy1|
       *				|dF1/dy2 ... dFn/dy2|
       *			J =	|   .		.   |
       *				|   .		.   |
       *				|dF1/dyn ... dFn/dyn|
       *
       * Meaning that all coefficients pertaining to yi are to be found in ROW i (as opposed to column i). 
       * The matrix is hence stored in column-wise (fortran) style.
       * Solving J.dy = -F(y) with dy = (dy1 ... dyn) and F(y) = (F1(y) ... Fn(y)) requires the variable
       * trans[1] to be defined as {"N"} (see programs/various/testlapack.c for details).
       */

      // Use QR factorization to solve the linear system
      memset((void *)JacCopy, 0, MAX_PNTDIM_SQ*sizeof(double));
      COPY(pntdim2, Jac, 1, JacCopy, 1);
      COPY(pntdim,  rhs,    1, dy,     1);
      SCAL(pntdim, -1.0,   dy, 1);

      GELS(trans, &pntdim1, &pntdim1, &nrhs, JacCopy, &pntdim1, dy, &pntdim1, work, &lwork, &info );

      if (info < 0)
	{
	  ErrorMsg(__FILE__, __LINE__, "Illegal value for parameter %d in GELS", abs(info));
	  return ILLEGAL_INPUT;
	}
      else if (info > 0)
	{
	  ErrorMsg(__FILE__, __LINE__, "Singular Jacobian matrix found in GELS\n");
	  return SINGULARITY;
	}

      // Adjust point
      AXPY(pntdim, 1.0, dy, 1, y, 1); //daxpy functie uit BLAS (AXPY(pntdim, alhpa,x,incrx,y,incry) = y=alpha*x+y). y=1.0*dy+y Y(guess)=nu dus toegenomen.
    }

  return NO_CONVERGENCE;
}


#endif // (MOOREPENROSE == 1)

/*==============================================================================*/

int	TangentVec(const int pntdim, double *sol, double *tanvec,
		   int (*fnc)(double *, double *))
  
  /*
   * TangentVec - routine determines the direction of the curve defined by the
   *	          system of equations
   *
   *					F(y) = 0
   *
   *		  The point y is considered to have a dimension of exactly 1
   *		  larger than the number of equations (i.e. the dimension of
   *		  F(y)).
   *
   * Arguments - pntdim	   : The dimension of the solution point on the curve.
   *		 y	   : Pointer to an array containing the fixed point
   *		 tanvec	   : Pointer to return tangent vector
   *		 fnc	   : Pointer to function specifying the system of
   *			     equations. The function must have a (double)
   *			     pointer as first argument, containing the point
   *			     in which to evaluate the system and a (double)
   *			     pointer as second argument, containing the
   *			     results after evaluation of the equations.
   */

{
  register int		j;
  int			sysdim = pntdim-1, matdim = pntdim*pntdim;
  int			rows   = pntdim,   cols   = pntdim;
  int			nrhs = 1, info;
  char			trans[1] = {"N"};
  double		norm;
  static int		first = 1, lastdim;
  double		y[pntdim];
  
  // Initialize
  COPY(pntdim, sol, 1, y, 1);
  norm = NRM2(pntdim, y, 1);
  if (!issane(norm))
    {
      ErrorMsg(__FILE__, __LINE__, "Norm overflow in curvedir");
      return NORM_OVERFLOW;
    }

  // Determine the Jacobian of the extended system (variable plus parameter
  // dependence).
  ReportMsg(2, "\nComputing curve direction     ");
  Jacobian(pntdim, y, sysdim, JacCopy, fnc, CENTRALDIFF);
  ReportMsg(2, ".....Ok!\n");

  // Append the current tangent vector as the last row to the jacobian to
  // preserve direction. See the matcont manual at
  // http://www.matcont.ugent.be/manual.pdf, page 10 & 11
  // Notice, however, it is here added as the last COLUMN because of the
  // Fortran column-wise storage!

  for (j=0; j<pntdim; j++)
    {
      COPY(sysdim, JacCopy+j*sysdim, 1, Jac+j*pntdim, 1);		// Extract dF/dy
      *(Jac+j*pntdim+sysdim) = tanvec[j];
    }

  memset((void *)JacCopy, 0, MAX_PNTDIM_SQ*sizeof(double));
  COPY(matdim, Jac, 1, JacCopy, 1);
  memset((void *)tanvec, 0, pntdim*sizeof(double));
  tanvec[sysdim] = 1.0;
  Stepchange = 0;

  // Use QR factorization to solve the linear system
  GELS(trans, &rows, &cols, &nrhs, JacCopy, &rows, tanvec, &rows, work, &lwork, &info );

  if (info < 0)
    {
      ErrorMsg(__FILE__, __LINE__, "Illegal value for parameter %d in GELS", abs(info));
      memset((void *)tanvec, 0, pntdim*sizeof(double));
      tanvec[0] = 1.0;
      return ILLEGAL_INPUT;
    }
  else if (info > 0)
    {
      ErrorMsg(__FILE__, __LINE__, "Singular Jacobian matrix found in GELS\n");
      memset((void *)tanvec, 0, pntdim*sizeof(double));
      tanvec[0] = 1.0;
      return SINGULARITY;
    }

  norm = NRM2(pntdim, tanvec, 1);		/* Normalize and store      */
  SCAL(pntdim,  1.0/norm, tanvec, 1);
  
  if ((first && (tanvec[0] < 0.0)) ||
      ((lastdim != pntdim) && (DOT(min(pntdim, lastdim), tanvec, 1, oldvec, 1) < 0.0)))
    SCAL(pntdim, -1.0, tanvec, 1);
  first = 0;
  lastdim = pntdim;

  memset((void *)oldvec, 0, MAX_PNTDIM*sizeof(double));
  COPY(pntdim, tanvec, 1, oldvec, 1);

  return SUCCES;
}



/*==============================================================================*/
#define MAX_EXTRAPOLATIONS	10

#define DIFFDELTA       1.0E-9                  // error goal
#define DIFFTOL         1.0E-9                  // relative error goal

int	Jacobian(const int pntdim, double *pnt,
		 const int fncdim, double *jac,
		 int (*fnc)(double *, double *), int method)
  /*
   * Routine determines the Jacobian of the n-dimensional function F(y) w.r.t. the m-dimensional
   * variable y at the current point given by 'pnt'. The routine hence returns in 'jac' the
   * following matrix of partial derivatives:
   *
   *				|dF1/dy1 ... dFn/dy1|
   *				|   .		.   |
   *			Df =	|   .		.   |
   *				|   .		.   |
   *				|dF1/dym ... dFn/dym|
   *
   * Notice that all coefficients pertaining to yi are to be found in ROW i (as opposed to column i).
   * The matrix is hence stored in column-wise (fortran) style.
   */


{
  register int		i, j, k, m, alldone;
  int			jacDone[fncdim];
  double		Djac[MAX_EXTRAPOLATIONS+1][MAX_EXTRAPOLATIONS+1][pntdim*fncdim];
  double		y[pntdim];
  double		ydif, yfac, old;
  double		rhs[fncdim];
  double		err, rerr;

  // Initialize
  COPY(pntdim, pnt, 1, y, 1);
  if (method == FORWARD)
    {
      if ((*fnc)(y, rhs) == FAILURE)
	{
	  ErrorMsg(__FILE__, __LINE__, "Right-hand side computation failed");
	  return FAILURE;
	}
    }

  memset((void *)jac, 0, (pntdim*fncdim)*sizeof(double));
  for (j=0; j<pntdim; j++)
    {
      checkESCpressed();
      old  = y[j];
      switch (method)
	{
	case FORWARD:
	  ydif = max(fabs(JACOBIAN_STEP*y[j]), MIN_STEP);
	  yfac = 1.0/ydif;
	  break;
	case CENTRAL:
	  ydif = max(fabs(JACOBIAN_STEP*y[j]), MIN_STEP);
	  yfac = 1.0/(2.0*ydif);
	  break;
	default:							// Richardsons extrapolation
	  ydif = max(fabs(0.1*y[j]), 0.1);					// Problematic??
	  yfac = 1.0/(2.0*ydif);
	  break;
	}

      memset((void *)jacDone,  0, fncdim*sizeof(int));
      for (i = 0; i < MAX_EXTRAPOLATIONS; i++)
	{
	  y[j] = old + ydif;
	  if ((*fnc)(y, Djac[i][0]+j*fncdim) == FAILURE)
	    {
	      ErrorMsg(__FILE__, __LINE__, "Right-hand side computation failed");
	      return FAILURE;
	    }

	  if (method != FORWARD)
	    {
	      y[j]  = old - ydif;
	      if ((*fnc)(y, rhs) == FAILURE)
		{
		  ErrorMsg(__FILE__, __LINE__, "Right-hand side computation failed");
		  return FAILURE;
		}
	    }
	  AXPY(fncdim,  -1.0, rhs, 1, Djac[i][0]+j*fncdim, 1);
	  SCAL(fncdim,   yfac, Djac[i][0]+j*fncdim, 1);

	  if (method == RICHARDSON)
	    {
	      // Use Richardson extrapolation to refine the derivatives
	      unsigned powerof4 = 1;
	      for (k = 0; k < i; k++)
		{
		  powerof4 <<= 2; // this value is equivalent to pow(4,k+1)
		  for (m=0; m<fncdim; m++)
		    Djac[i][k+1][j*fncdim+m] = Djac[i][k][j*fncdim+m] + (Djac[i][k][j*fncdim+m] - Djac[i-1][k][j*fncdim+m]) / (powerof4 -1);
		}

	      if (i)
		{
		  // compute errors
		  alldone = 1;
		  for (m=0; m<fncdim; m++)
		    {
		      if (jacDone[m]) continue;
		      err  = fabs(Djac[i][i][j*fncdim+m]-Djac[i-1][i-1][j*fncdim+m]);
		      rerr = 2*err/(fabs(Djac[i][i][j*fncdim+m])+fabs(Djac[i-1][i-1][j*fncdim+m]) + MACH_PREC);
		      if ((rerr < DIFFTOL) || (err < DIFFDELTA))
			{
			  jac[j*fncdim+m] = Djac[i][i][j*fncdim+m];
			  jacDone[m] = 1;
			}
		      else alldone = 0;
		    }
		  if (alldone) break;
		}
	      ydif *= 0.5;
	      yfac  = 1.0/(2.0*ydif);
	    }
	  else	// In case of FORWARD and CENTRAL differencing finish right after i = 0 (single pass)
	    {
	      COPY(fncdim,  Djac[0][0]+j*fncdim, 1, jac+j*fncdim, 1);
	      break;
	    }
	}
      if (i == MAX_EXTRAPOLATIONS)
 	{
	  ErrorMsg(__FILE__, __LINE__, "Richardsons extrapolation failed");

	  // In case the Richardsons extrapolation failed, fall back on central differencing
 	  COPY(fncdim,  Djac[0][0]+j*fncdim,   1, jac+j*fncdim, 1);

#if (DEBUG == 1)
	  fprintf(stderr, "\n\ndF1/dy%d ... dFn/dy%d:\n", j, j);
	  fprintf(stderr, "====================\n");
 	  for (k=0; k<MAX_EXTRAPOLATIONS; k++)
 	    {
 	      for (m=0; m<fncdim; m++)
 		fprintf(stderr, "\t%14.8G", Djac[k][k][j*fncdim+m]);
 	      fprintf(stderr, "\n");
 	    }

 	  fprintf(stderr, "\nAbsolute error:\n===============\n");
 	  for (k=1; k<MAX_EXTRAPOLATIONS; k++)
 	    {
 	      for (m=0; m<fncdim; m++)
 		fprintf(stderr, "\t%14.8G", fabs(Djac[k][k][j*fncdim+m]-Djac[k-1][k-1][j*fncdim+m]));
 	      fprintf(stderr, "\n");
 	    }

 	  fprintf(stderr, "\nRelative error:\n===============\n");
 	  for (k=1; k<MAX_EXTRAPOLATIONS; k++)
 	    {
 	      for (m=0; m<fncdim; m++)
 		fprintf(stderr, "\t%14.8G", 2*fabs(Djac[k][k][j*fncdim+m]-Djac[k-1][k-1][j*fncdim+m])/(fabs(Djac[k][k][j*fncdim+m])+fabs(Djac[k-1][k-1][j*fncdim+m]) + MACH_PREC));
 	      fprintf(stderr, "\n");
 	    }
 	  fprintf(stderr, "\n\n");
#endif
 	}
      y[j] = old;
    }

  return SUCCES;
}



/*==============================================================================*/

int	CurveFuncDeriv(const int pntdim,  double *pnt, double *dy,
		       const int fncdim, double *dfuncdp,
		       int (*fnc)(double *, double *), int method)

  /*
   * CurveFuncDeriv  - routine determines the derivative of the variables y w.r.t.
   *                   to the curve parameter p along the curve defined by the
   *		       system of equations
   *
   *					F(y,p) = 0
   *
   *		       The vector y should have the same dimension as the number
   *		       of equations (i.e. the dimension of F(y,p)), which should
   *		       be equal to pntdim-1. The vector argument 'pnt' contains
   *		       the value of y and as last element the value of the parameter
   *		       p and hence has a dimension equal to 'pntdim'.
   *		       In addition, the routine determines the derivative w.r.t.
   *		       the curve parameter p of functions G(y, p) of the
   *		       variables y and the parameter p. The dimension of the
   *		       function vector G is determined by 'fncdim'
   *
   * Arguments - pntdim	   : The dimension of the argument vector 'pnt'.
   *		 pnt	   : Pointer to an array containing as first elements the
   *		 	     the state variables y and as the last element the
   *		 	     value of the parameter p. Together 'pnt' contains
   *		 	     the coordinates of a point on the curve determined
   *		 	     by F(y, p).
   *		 dy	   : Pointer to return vector with dy/dp values
   *		 fncdim    : The dimension of the function vector G(y, p)
   *		 dfuncdp   : Point to return vector with dG/dp values
   *		 fnc	   : Pointer to function specifying the system of
   *			     equations. The function must have a (double)
   *			     pointer as first argument, containing the point
   *			     in which to evaluate the system and a (double)
   *			     pointer as second argument, containing the
   *			     results after evaluation of the equations.
   *			     This array of results should contain as first elements
   *			     the values of F(y,p), following by the values of G(y, p).
   *		 method	   : FORWARD, CENTRAL or RICHARDSON, determining the manner to
   *		 	     compute the derivative
   */

{
  register int		j, k;
  int			sysdim  = pntdim-1;
  int			nrhs = 1, info;
  char			trans[1] = {"N"};
  double		J[pntdim*((pntdim-1)+fncdim)];
  double		jac[pntdim*(pntdim-1)];
  double		dfuncdy[pntdim*fncdim];
  double		work[4*pntdim];
  int			lwork = 4*pntdim;

  // Determine the Jacobian of the extended system (variable plus parameter
  // dependence). Use central differencing for the derivatives.
  // The derivatives w.r.t. to the parameter p are stored in the last row, because
  // the last vector element of the argument 'pnt' is the parameter.
  if (Jacobian(pntdim, pnt, sysdim+fncdim, J, fnc, method) == FAILURE)
    return FAILED_EVALUATION;

#if (DEBUG == 1)
  int m, n;
  fprintf(stderr, "\nJacobian:\n=========\n");
  for (m=0; m<pntdim; m++)
    {
      for (n=0; n<(sysdim+fncdim); n++)
	fprintf(stderr, "\t%14.8G",J[m*(sysdim+fncdim)+n]);
      fprintf(stderr, "\n");
    }
#endif

  // The value dG(y(p), p)/dp is computed by:
  //
  //    dG(y(p), p)/dp = del(G(y(p), p))/del(y) dy/dp + del(G(y(p), p))/del(p)
  //
  // where dy/dp is solved from:
  //
  //   del(F(y(p),p))/del(y) dy/dp + del(F(y(p), p))/del(p) = 0
  //
  for (j=0; j<pntdim; j++)
    {
      COPY(sysdim, J+j*(sysdim+fncdim), 1, jac+j*sysdim, 1);		// Extract dF/dy & dF/dp
      COPY(fncdim, J+j*(sysdim+fncdim)+sysdim, 1, dfuncdy+j*fncdim, 1);	// Extract dG/dy & dG/dp
    }
  COPY(sysdim, jac+(pntdim-1)*sysdim, 1, dy, 1);			// Store -dF/dp in dy
  SCAL(sysdim, -1.0, dy, 1);

  // dy/dp is solved from:
  //
  //   del(F(y(p),p))/del(y) dy/dp + del(F(y(p), p))/del(p) = 0
  //
  // Solving DF[1..n][1..n].dy = -DF[n+1][1..n] with dy = (dy1/dp ... dyn/dp) requires the variable
  // trans[1] to be defined as {"N"} (see programs/various/testlapack.c for details).
  GELS(trans, &sysdim, &sysdim, &nrhs, jac, &sysdim, dy, &sysdim, work, &lwork, &info );

  if (info < 0)
	{
	  ErrorMsg(__FILE__, __LINE__, "Illegal value for parameter %d in GELS", abs(info));
	  return ILLEGAL_INPUT;
	}
  else if (info > 0)
	{
	  ErrorMsg(__FILE__, __LINE__, "Singular Jacobian matrix found in GELS\n");
	  return SINGULARITY;
	}
  /*
   * Now dy contains the derivatives (dy1/dp, ..., dyn/dp) of the variables w.r.t. to the
   * curve parameter p. To compute dG/dp the matrix DG[1..n][1..m] is multiplied with the
   * vector dy and the last row of the matrix DG[n+1][1..m] is added to this result
   */
#if (DEBUG == 1)
  fprintf(stderr, "\ndY/dP:\n======\n");
  for (m=0; m<sysdim; m++)
    fprintf(stderr, "\t%14.8G",dy[m]);
  fprintf(stderr, "\n");
#endif

  // The value dG(y(p), p)/dp is computed by:
  //
  //    dG(y(p), p)/dp = del(G(y(p), p))/del(y) dy/dp + del(G(y(p), p))/del(p)
  //
  for (k=0; k<fncdim; k++)
    {
      dfuncdp[k] = 0.0;
      for(j=0; j<sysdim; j++)
	{
	  dfuncdp[k] += dfuncdy[j*fncdim+k]*dy[j];
	}
      dfuncdp[k] += dfuncdy[j*fncdim+k];	// j has value sysdim, the index of dGi/dp
    }

#if (DEBUG == 1)
  fprintf(stderr, "\ndG/dy:\n======\n");
  for (m=0; m<sysdim; m++)
    fprintf(stderr, "\t%14.8G",dfuncdy[m]);
  fprintf(stderr, "\n");

  fprintf(stderr, "\ndG/dP:\n======\n");
  for (m=0; m<fncdim; m++)
    fprintf(stderr, "\t%14.8G",dfuncdp[m]);
  fprintf(stderr, "\n\n");
#endif

  return SUCCES;
}


/*==============================================================================*/

double		LPcondition(const int pntdim, double *y, int (*fnc)(double *, double *), int method)

  /*
   * LPcondition - Routine computes the factor determining the location of a limit point, i.e.
   * 		   the parameter component of the tangent vector. This component always has
   * 		   index 0 in the vector of the solution point.
   *
   * Arguments - pntdim    : The dimension of the argument vector 'y'. Notice that this
   * 			     equals 2 (for the parameters) plus the dimension of the vector of
   * 			     state variables.
   *		 y	   : Pointer to an array containing as first element the value
   *		 	     of the parameter p and as subsequent elements the values of
   *		 	     the state variables y. The last element is assumed to be the
   *		 	     second parameter in the LP continuation
   *		 fnc	   : Pointer to function specifying the system of
   *			     equations. The function must have a (double)
   *			     pointer as first argument, containing the point
   *			     in which to evaluate the system and a (double)
   *			     pointer as second argument, containing the
   *			     results after evaluation of the equations.
   */
{
  static int		LPImmediateReturn = 0;
  register int		j;
  const int		lppntdim = pntdim-1;
  int			sysdim = lppntdim-1, matdim = lppntdim*lppntdim;
  int			rows   = lppntdim,   cols   = lppntdim;
  int			ldm    = lppntdim,   ldv    = lppntdim;
  int			ipiv[lppntdim],      nrhs = 1;
  int			info, iwork[lppntdim];
  int			maxind;
  char			trans[1] = {"N"}, whichnorm[1] = {"1"};
  double		rhs[pntdim-1];
  double		jacmat[pntdim*pntdim];
  double		ludmat[lppntdim*lppntdim];
  double		tanvec[lppntdim], work[4*lppntdim];
  double		norm, maxcond, cond;

  // Prevent recurrence
  if (LPImmediateReturn) return 0.0;
  LPImmediateReturn = 1;

  memset((void *)jacmat, 0, pntdim*pntdim*sizeof(double));
  // Determine the Jacobian of the extended system (variable plus parameter dependence).
  // Notice that we have to call Jacobian() with the full dimension (pntdim) of y to pass
  // the entire argument vector to Equation(). As a result, the Jacobian matrix will have
  // one row more than we really need, representing the derivatives w.r.t. to the 2nd parameter
  // of the LP continuation. This additional row, with index lppntdim = pntdim-1 will be ignored.
  // Also notice that the last column will be 0, i.e. the value of LPcondition() when
  // called from this routine.
  Jacobian(pntdim, y, pntdim-1, jacmat, fnc, method);
  for (j=0; j<lppntdim; j++) jacmat[j*lppntdim+sysdim] = 0.0;

  /*
    * The resulting Jacobian hence equals the following matrix of partial derivatives:
    *
    *				|dF1/dy1 ... dFn/dy1  0|
    *				|   .		.     0|
    *			Df =	|   .		.     0|
    *				|   .		.     0|
    *				|dF1/dym ... dFn/dym  0|
    *
    * In which m = n+1. Notice that all coefficients pertaining to yi are to be found
    * in ROW i (as opposed to column i). The matrix is hence stored in column-wise (fortran) style.
    */

  // Additional call to reset global variables
  (*fnc)(y, rhs);

  // Find the most non-singular matrix (largest inverse condition)
  for (j=0, maxcond=0.0, maxind=-1; j<lppntdim; j++)
    {
      // LU decompose the matrix and compute determinant
      COPY(matdim, jacmat, 1, ludmat, 1);
      ludmat[j*lppntdim+sysdim] = 1.0;
      norm = anorm(lppntdim, lppntdim, ludmat);
      GETRF(&rows, &cols, ludmat, &ldm, ipiv, &info);
      if (info < 0)
	{
	  ErrorMsg(__FILE__, __LINE__,
		   "Illegal value for parameter %d in DGETRF", abs(info));
	  exit(1);
	}
      else if (info == 0)
	{
	  GECON(whichnorm, &cols, ludmat, &ldm, &norm, &cond, work,
		iwork, &info);
	  if (info < 0)
	    {
	      ErrorMsg(__FILE__, __LINE__,
		       "Illegal value for parameter %d in DGECON", abs(info));
	      exit(1);
	    }
	  if (cond > maxcond+FEMTO)
	    {
	      maxcond = cond;
	      maxind  = j;
	    }
	}
    }

  if (maxind == -1)
    {
      ErrorMsg(__FILE__, __LINE__,
	       "No non-singular matrix found in LPcondition()");
      exit(1);
    }

  jacmat[maxind*lppntdim+sysdim] = 1.0;
  COPY(matdim, jacmat, 1, ludmat, 1);
  memset((void *)tanvec, 0, lppntdim*sizeof(double));
  tanvec[sysdim] = 1.0;

  // Solve for tangent vector
  GELS(trans, &rows, &cols, &nrhs, ludmat, &ldm, tanvec, &ldv, work, &lwork, &info );

  if (info < 0)
	{
	  ErrorMsg(__FILE__, __LINE__, "Illegal value for parameter %d in GELS", abs(info));
	  return ILLEGAL_INPUT;
	}
  else if (info > 0)
	{
	  ErrorMsg(__FILE__, __LINE__, "Singular Jacobian matrix found in GELS\n");
	  return SINGULARITY;
	}

  norm = NRM2(lppntdim, tanvec, 1);		/* Normalize and store      */
  SCAL(lppntdim,  1.0/norm, tanvec, 1);

  LPImmediateReturn = 0;

  return tanvec[0];
}

/*==============================================================================*/

int Evogradient(const int pntdim, const int evodim, double *pnt, int (*fnc)(double *, double *), int *evopars, double *selgrad, double *minbounds, double *maxbounds)
{
   
    double y[pntdim];
    double selec[evodim];
    
    int i;
    
    COPY(pntdim, pnt, 1, y, 1); //copieert pnt naar vec y
    selec[0] =  SelectionGradient(BASE_PNTDIM, y, Equation, 0, 0);
    if ((y[0] <=minbounds[0]&&selec[0]<0) || (y[0]>maxbounds[0]&&selec[0]>0))
    {
        selec[0] = 0;
    }
    for (i=1;i<evodim;i++)
    {
        selec[i] =  SelectionGradient2(BASE_PNTDIM, y, Equation, evopars[i], 0);
        if ((parameter[evopars[i]] <= minbounds[i] && selec[i]<0) || (parameter[evopars[i]] >= maxbounds[i] && selec[i]>0))
        {
            selec[i] = 0;
        }
    }
    
    
    COPY(evodim, selec, 1, selgrad, 1);
    return SUCCES;
}


/*==============================================================================*/
bool Noevocheck(double *selgrad)
{
    for (int i=0; i < EVODIM; i++)
    {
        if (selgrad[i] > CANTOL)
            return false;
    }
    return true;
}
/*==============================================================================*/


double		SelectionGradient(const int pntdim, double *pnt, int (*fnc)(double *, double *), int varindex, int R0index)

  /*
   * SelectionGradient - Routine computes the factor determining the location of an ESS, i.e.
   * 		         the component of the Jacobian that represents the derivative w.r.t.
   * 		         to the parameter, which should always has index 0 in the vector of
   * 		         the solution point.
   *
   * Arguments - pntdim    : The dimension of the argument vector 'y'. Notice that this
   * 			     equals 2 (for the parameters) plus the dimension of the vector of
   * 			     state variables.
   *		 y	   : Pointer to an array containing as first element the value
   *		 	     of the parameter p and as subsequent elements the values of
   *		 	     the state variables y. The last element is assumed to be the
   *		 	     second parameter in the LP continuation
   *		 fnc	   : Pointer to function specifying the system of
   *			     equations. The function must have a (double)
   *			     pointer as first argument, containing the point
   *			     in which to evaluate the system and a (double)
   *			     pointer as second argument, containing the
   *			     results after evaluation of the equations.
   *        varindex Index of y waar parameter staat.
   *		R0index	   : Index of the equilibrium condition (R0-1 = 0) in the result vector
   *			     returned by (*fnc)().
   */
{
  static int		SGImmediateReturn = 0;
  double		y[pntdim], old, ydif;
  double		rhs[pntdim-1], dR0dp;

  // Prevent recurrence
  if (SGImmediateReturn) return 0.0;
  SGImmediateReturn = 1;

  // Initialize
  COPY(pntdim, pnt, 1, y, 1); //copieert pnt naar y vec
  old  = y[varindex];   ///Waarde van de parameter waarna de R0 wordt gedif
  ydif = max(fabs(JACOBIAN_STEP*y[varindex]), MIN_STEP);  ///stap waarmee de parameeter wordt veranderd

  y[varindex] = old + ydif;   ///Nieuwe parameterwaarde, y hoger dan oude
  if ((*fnc)(y, rhs) == FAILURE)  ///Hier wordt het systeem opnieuw ge-evalueerd, dus met de nieuwe parameter en met oude environment
    {
      ErrorMsg(__FILE__, __LINE__, "Right-hand side computation failed");
      return FAILURE;
    }
  dR0dp = rhs[R0index];   ///dit is de nieuwe R0 na evaluatie van de nieuwe y (R0+)

  y[varindex]  = old - ydif;  //Nieuwe parameterwaarde, y lager dan oude
  if ((*fnc)(y, rhs) == FAILURE)
    {
      ErrorMsg(__FILE__, __LINE__, "Right-hand side computation failed");
      return FAILURE;
    }
  dR0dp -= rhs[R0index];  /// R0+ - R0-
  dR0dp /= 2.0*ydif;

  y[varindex] = old;
  // Additional call to reset global variables
  (*fnc)(y, rhs);

  SGImmediateReturn = 0;

  return dR0dp;					// Return the derivative w.r.t. the scaled or unscaled parameter?
//  return (dR0dp/pnt_scale[0]);			// par=scaledpar*scale => dR0/dscaledpar = dR0/d(par/scale) = scale*DR0/dpar
}



/*============================*/
 
double		SelectionGradient2(const int pntdim, double *pnt, int (*fnc)(double *, double *), int parindex, int R0index)
 
 /*
 * SelectionGradient - Routine computes the factor determining the location of an ESS, i.e.
 * 		         the component of the Jacobian that represents the derivative w.r.t.
 * 		         to the parameter, which should always has index 0 in the vector of
 * 		         the solution point.
 *
 * Arguments - pntdim    : The dimension of the argument vector 'y'. Notice that this
 * 			     equals 2 (for the parameters) plus the dimension of the vector of
 * 			     state variables.
 *		 y	   : Pointer to an array containing as first element the value
 *		 	     of the parameter p and as subsequent elements the values of
 *		 	     the state variables y. The last element is assumed to be the
 *		 	     second parameter in the LP continuation
 *		 fnc	   : Pointer to function specifying the system of
 *			     equations. The function must have a (double)
 *			     pointer as first argument, containing the point
 *			     in which to evaluate the system and a (double)
 *			     pointer as second argument, containing the
 *			     results after evaluation of the equations.
 *        parindex Parameterindex (anders dan in selectiongradient waar het de var index is, nu is het parameter index.
 *		R0index	   : Index of the equilibrium condition (R0-1 = 0) in the result vector
 *			     returned by (*fnc)().
 */
{
    static int	SGImmediateReturn = 0;
    double		y[pntdim], old, ydif;
    double		rhs[pntdim-1], dR0dp;
    
    // Prevent recurrence
    if (SGImmediateReturn) return 0.0;
    SGImmediateReturn = 1;
    
    // Initialize
    COPY(pntdim, pnt, 1, y, 1);  //punt(guess) wordt naar yvec gekop)
    old  = parameter[parindex];   ///Waarde van de parameter waarna de R0 wordt gedif
    ydif = max(fabs(JACOBIAN_STEP*parameter[parindex]), MIN_STEP);  ///stap waarmee de parameter wordt veranderd.
    parameter[parindex] = old + ydif;   ///Nieuwe parameterwaarde, y hoger dan oude
    if ((*fnc)(y, rhs) == FAILURE)  ///Hier wordt het systeem opnieuw ge-evalueerd, dus met de nieuwe parameter en oude environment
    {
        ErrorMsg(__FILE__, __LINE__, "Right-hand side computation failed");
        return FAILURE;
    }
    dR0dp = rhs[R0index];   ///dit is de nieuwe R0 na evaluatie van de nieuwe y (R0+)

    parameter[parindex]  = old - ydif;  //Nieuwe parameterwaarde, y lager dan oude
    if ((*fnc)(y, rhs) == FAILURE)
    {
        ErrorMsg(__FILE__, __LINE__, "Right-hand side computation failed");
        return FAILURE;
    }
    dR0dp -= rhs[R0index];  /// R0+ - R0-
    dR0dp /= 2.0*ydif;
 
    parameter[parindex] = old;
    // Additional call to reset global variables
    (*fnc)(y, rhs);
 
    SGImmediateReturn = 0;
 
    
    return dR0dp;					// Return the derivative w.r.t. the scaled or unscaled parameter?
    //  return (dR0dp/pnt_scale[0]);			// par=scaledpar*scale => dR0/dscaledpar = dR0/d(par/scale) = scale*DR0/dpar
 
}



/*==============================================================================*/

typedef	int		USERFUN(double *, double *);
static USERFUN		*UserFun = NULL;
static double		*tv, *yinit, *dy;

void nleq_callback(int *argn, double *y, double *fy, int *iflag)
{
  (*UserFun)(y, fy);

  COPY(*argn, y, 1, dy, 1);
  AXPY(*argn, -1.0, yinit, 1, dy, 1);
  fy[*argn-1] = DOT(*argn, dy, 1, tv, 1);

  *iflag=0;
  return;
}


/*==============================================================================*/

int	FindPointGlobal(const int pntdim, double *guess, double *tanvec,
		        double ytol, double rhstol, const int max_iter, int (*fnc)(double *, double *))

  /*
   * FindPointGlobal -	Routine locates a point on a curve determined by a
   *			system of non-linear, algebraic equations.
   *			This routine uses the ERRor-based Damped Newton algorithm (NLEQ_ERR)
   *  			written by L. Weimann, as presented in P. Deuflhard (2004), Newton
   *  			Methods for Nonlinear Problems. - Affine Invariance and Adaptive
   *  			Algorithms. Series Computational Mathematics 35, Springer
   *
   *			See file nleq_err.c for more details
   *
   * Arguments - pntdim	   : The dimension of the solution point on the curve.
   *			     The dimension of the system of equations is
   *			     assumed to be exactly 1 less.
   *		 guess	   : Pointer to an array containing the initial point
   *			     to start the iteration from. The first element of
   *			     the vector is assumed to be non-adjustable parameter.
   *		 ytol	   : Tolerance determining when change in y equals zero.
   *		 rhstol	   : Tolerance determining when RHS equals zero.
   *		 max_iter  : Maximum stepnumber allowed in iteration.
   *		 fnc	   : Pointer to function specifying the system of
   *			     equations. The function must have a (double)
   *			     pointer as first argument, containing the point
   *			     in which to evaluate the system and a (double)
   *			     pointer as second argument, containing the
   *			     results after evaluation of the equations.
   */
{
  double		y[pntdim];
  double		rhs[pntdim];
  double		ydif[pntdim];

  UserFun = fnc;
  tv      = tanvec;
  yinit   = guess;
  dy      = ydif;

  COPY(pntdim, guess, 1, y, 1);

  struct NLEQ_FUN 	problem;
  struct NLEQ_OPT 	opt;
  struct NLEQ_INFO 	info;

  problem.fun = nleq_callback;
  problem.jac = NULL;

  opt.tol		= rhstol;
  opt.maxiter		= max_iter;
  opt.errorlevel	= Minimum;
  opt.monitorlevel	= None;
  opt.datalevel		= Verbose;
  opt.errorfile		= errfile;
  opt.monitorfile	= NULL;
  opt.datafile		= errfile;
  opt.iterfile		= NULL;
  opt.resfile		= NULL;
  opt.miscfile		= NULL;
  opt.nonlin		= Highly_Nonlinear;
  opt.restricted	= False;
  opt.scale		= NULL;

  nleq_res(problem, pntdim, y, &opt, &info);

  if ( info.rcode == 0 )
    {
      COPY(pntdim, y, 1, guess, 1);
      memset((void *)rhs,  0, pntdim*sizeof(double));
      (*fnc)(y, rhs);

      // Compute rhsnorm
      RhsNorm = NRM2(pntdim, rhs, 1);
      RhsNorm = RhsNorm/(1.0+RhsNorm);
      return SUCCES;
    }
  else
    {
      switch ( info.rcode )
      {
	case -999:
	  ErrorMsg(__FILE__, __LINE__, "Routine nleq_fwalloc failed to allocate double memory via malloc");
	  break;
	case -998:
	  ErrorMsg(__FILE__, __LINE__, "Routine nleq_iwalloc failed to allocate double memory via malloc");
	  break;
	case -997:
	  ErrorMsg(__FILE__, __LINE__, "Routine nleq_pfwalloc failed to allocate double memory via malloc");
	  break;
	case -995:
	  ErrorMsg(__FILE__, __LINE__, "Internal i/o control block could not be allocated via malloc");
	  break;
	case -994:
	  ErrorMsg(__FILE__, __LINE__, "Internally used data structure could not be allocated via malloc");
	  break;
	case -989:
	  ErrorMsg(__FILE__, __LINE__, "Default data-output file could not be opened via fopen call");
	  break;
	case -99:
	  ErrorMsg(__FILE__, __LINE__, "NULL pointer obtained from funs.fun field - the problem function must be defined");
	  break;
	case 1:
	  ErrorMsg(__FILE__, __LINE__, "Singular Jacobian matrix (detected by routine nleq_linfact), nleq_res cannot proceed the iteration");
	  break;
	case 2:
	  ErrorMsg(__FILE__, __LINE__, "Maximum number of Newton iteration (as set by opt->maxiter) exceeded");
	  break;
	case 3:
	  ErrorMsg(__FILE__, __LINE__, "No convergence of Newton iteration, damping factor became too small");
	  break;
	case 20:
	  ErrorMsg(__FILE__, __LINE__, "Nonpositive input for dimensional parameter n");
	  break;
	case 21:
	  ErrorMsg(__FILE__, __LINE__, "Nonpositive value for opt->dxtol supplied");
	  break;
	case 22:
	  ErrorMsg(__FILE__, __LINE__, "Negative scaling value for some component of vector opt->scale supplied");
	  break;
	case 80:
	  ErrorMsg(__FILE__, __LINE__, "Routine nleq_linfact returned with an error other than singular Jacobian");
	  break;
	case 81:
	  ErrorMsg(__FILE__, __LINE__, "Routine nleq_linsol returned with an error");
	  break;
	case 82:
	  ErrorMsg(__FILE__, __LINE__, "The user defined problem function funs.fun returned a nonzero code other than 1 or 2");
	  break;
	case 83:
	  ErrorMsg(__FILE__, __LINE__, "The user defined Jacobian function funs.jac returned a nonzero code");
	  break;
      }
      return NO_CONVERGENCE;
    }
}


/*==============================================================================*/
