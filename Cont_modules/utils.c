/***
   NAME
     utils
   DESCRIPTION
     Utility functions
   Last modification: AMdR - Oct 03, 2013
***/
#ifndef UTILS
#define UTILS
#endif

#include "globals.h"


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

double anorm(int rows, int cols, double *a)
{
  register int		i;
  double		tmp, maxval = 0.0;
  
  for (i=0; i<rows; i++)
    {
      // AMdR - Sep 02, 2009: Error corrected: i is a row selector!
      // tmp = ASUM(cols, a+i*rows, 1);
      tmp = ASUM(cols, a+i*cols, 1);
      maxval = max(maxval, tmp);
    }

  return maxval;
}


/*==============================================================================*/


#define SAFETY		0.8
#define MIN_SCALE	1.0E-6
#define MAX_SCALE	1.0E6

static int		oldscale[MAX_PNTDIM];
const  double		upper = ((2.0-SAFETY)*10.0), lower = SAFETY;


int	SetScalesSimple(double *point, int basedim)

/*
 * This routine scales the vector of problem variables to within reasonable
 * bounds. This makes them more comparable, which is advantageous for the
 * computations. 
 *
 * This routine is only to be used for problems, in which all variables are
 * independent from each other. It is therefore not to be used to computation
 * of equilibria in structured population models with an infinite-dimensional
 * environment. In that case all variables can perhaps not be scaled
 * independently.
 *
 * For example in cannibalism problems, the cannibalistic pressure an
 * individual is exposed to at different body sizes should be scaled for all
 * victim body sizes identically! In these cases the user should supply a
 * routine that is largely similar to this one, but does not scale the elements
 * of the infinite-dimensional environment independently.
 */

{
  register int		i, scaleset = 0, newscale;
  static int		first = 1;
  double		tmp, val;

  if (first)
    {
      for (i=0; i<basedim; i++) pnt_scale[i] = 1.0;
      for (i=0; i<basedim; i++) oldscale[i]  = 0;
    }
  first = 0;
  
  for (i=0; i<basedim; i++)
    {
      val = fabs(point[i]);
      if ((val < lower) || (val > upper))
	{
	  tmp = max(val*pnt_scale[i], MIN_SCALE);
	  tmp = min(tmp, MAX_SCALE);
	  tmp = floor(log10(tmp)+FUNCTOL);
	  newscale = (int)tmp;

	  if (newscale != oldscale[i])
	    {
	      tmp           = pow(10.0, tmp);
	      point[i]     *= pnt_scale[i]/tmp;
	      pnt_scale[i]  = tmp;
	      oldscale[i]   = newscale;
	      scaleset	= i+1;
	    }
	}
    }

	return scaleset; 
}


/*==============================================================================*/

int	SetScales(double *point, int basedim, int dim)

{
  register int		i, scaleset = 0, newscale;
  static int		first = 1;
  double		tmp, maxicm;
#if CANONICAL == 1

      for (i=0; i<(basedim+1); i++) pnt_scale[i] = 1.0;
      for (i=0; i<(basedim+1); i++) oldscale[i]  = 1.0;

  return scaleset;
#else
  if (first)
    {
      for (i=0; i<(basedim+1); i++) pnt_scale[i] = 1.0;
      for (i=0; i<(basedim+1); i++) oldscale[i]  = 0;
    }
  first = 0;
  
  scaleset = SetScalesSimple(point, basedim);

  if (dim > basedim)
    {
      for(i=basedim, maxicm = 0.0; i<dim; i++) maxicm = max(maxicm, point[i]);
      if ((maxicm < lower) || (maxicm > upper))
        {
          tmp = max(maxicm*pnt_scale[basedim], MIN_SCALE);
          tmp = min(tmp, MAX_SCALE);
          tmp = floor(log10(tmp)+FUNCTOL);
          newscale = (int)tmp;

          if (newscale != oldscale[basedim])
	    {
	      tmp         = pow(10.0, tmp);
	      for(i=basedim; i<dim; i++)
	        point[i] *= pnt_scale[basedim]/tmp;
	      pnt_scale[basedim]  = tmp;
	      oldscale[basedim]   = newscale;
	      scaleset = basedim+1;
	    }
        }
    }

  return scaleset;
#endif
}


/*==============================================================================*/

double		determinant(double *mat, int matdim)

/*
 * determinant - Routine computes the determinant of the matrix mat using a LU
 *		 decomposition (or Gaussian elimination) with complete pivot
 *		 selection (see Stoer & Burlisch, 4.1 adn 4.5).
 *
 * Return - the value of the determinant.
 *
 */

{
  double	 *lftupcor, maxcol, divisor;
  double	 tmprow[matdim], tmpdbl, det1, ten;
  int		 row, tmpint, detsign = 1, det2;
  register int	 i, j, k;

  for(i=0; i<matdim; i++)
    {
      // Complete pivoting is chosen for stability
      // Search for the pivot element in the minor matrix
      lftupcor = (mat+i*(matdim+1));
      row = i; tmpint = i;
      maxcol = fabs(*lftupcor);
      for(j=i; j<matdim; j++)
	for(k=i; k<matdim; k++)
	  {
	    tmpdbl = fabs(*(mat+j*matdim+k));
	    if(tmpdbl > maxcol)
	      {
		maxcol = tmpdbl;
		row = j;
		tmpint = k;
	      }
	  }

      if (maxcol < ATTO) return 0.0;

      // Swap the row with largest element to ith row.
      if (row > i)
	{
	  detsign *= -1;
	  memcpy(tmprow,		  lftupcor,		   (matdim-i)*sizeof(double));
	  memcpy(lftupcor,		  lftupcor+(row-i)*matdim, (matdim-i)*sizeof(double));
	  memcpy(lftupcor+(row-i)*matdim, tmprow,		   (matdim-i)*sizeof(double));
	  /*
	  for(j=0; j<matdim-i; j++)
	    {
	      tmpdbl			 = lftupcor[j];
	      lftupcor[j]		 = lftupcor[(row-i)*matdim+j];
	      lftupcor[(row-i)*matdim+j] = tmpdbl;
	    }
	  */
	}

      // Swap the column with largest element to ith column.
      if (tmpint > i)
	{
	  detsign *= -1;
	  for(j=0; j<matdim; j++)
	    {
	      tmpdbl		   = mat[j*matdim+i];
	      mat[j*matdim+i]	   = mat[j*matdim+tmpint];
	      mat[j*matdim+tmpint] = tmpdbl;
	    }
	}

      // Create the next matrix (Stoer & Burlish, pg 163).
      for(j=1; j<matdim-i; j++)
	{
	  if (*lftupcor) *(lftupcor+j*matdim) /= *lftupcor;
	  else return MISSING_VALUE;

	  divisor = *(lftupcor+j*matdim);
	  for(k=1; k<matdim-i; k++)
	    *(lftupcor+j*matdim+k) -= divisor**(lftupcor+k);
	}
    }

  det1 = 1.0;
  det2 = 0;
  ten  = 10.0;
  for (i=0; i<matdim; i++)
    {
      det1 *= mat[i*matdim+i];
      if (det1 == 0.0) return 0.0;

      while (fabs(det1) < 1.0)
	{
	  det1 *= ten;
	  det2 -= 1;
	}
      while (fabs(det1) > 10.0)
	{
	  det1 /= ten;
	  det2 += 1;
	}
    }

  return (detsign*det1*pow(10.0, det2));
}


/*===========================================================================*/

int ComputeEigenvector(const int matdim, double *mat, double eigenval, double *eigenvec)

/*
 * Routine computes the normalized, right eigenvector of the matrix mat,
 * pertaining to the eigenvalue eigenval. The result is returned in
 * eigenvec. The eigenvector is normalized such that the sum of its elements
 * equals 1.
 */
{
  register int		i, j;
  int			maxind;
  int			rows   = matdim,   cols   = matdim;
  int			ldm    = matdim,   ldv    = matdim;
  double		cond, maxcond, norm, ferr[1], berr[1];
  double		localmat[matdim*matdim];
  double		ludmat[matdim*matdim];
  int			ipiv[matdim], nrhs = 1, info;
  char			trans[1] = {"N"}, whichnorm[1] = {"1"};
  int			iwork[matdim];
  double		work[4*matdim];
  double		tv0[matdim];
  double		tv1[matdim];

  if (Report2Terminal)
    (void)fprintf(stderr,  "%-30s", "Finding non-singular matrix");
  if (errfile)
    (void)fprintf(errfile, "%-30s", "Finding non-singular matrix");

  for (i=0, maxcond=0.0, maxind=-1; i<matdim; i++)
    {
      if (Report2Terminal) (void)fprintf(stderr,  ".");
      if (errfile)	   (void)fprintf(errfile, ".");

      // Make a working copy of the matrix for LU decomposition
      memcpy(localmat, mat, matdim*matdim*sizeof(double));

      // Substract the eigenvalue
      for (j=0; j<matdim; j++) localmat[j*matdim+j] -= eigenval;

      // Replace row i of the matrix with only 1's
      // NB: Fortran column-wise storage!!
      for (j=0; j<matdim; j++) localmat[j*matdim+i] = 1.0;

      // LU decompose the matrix and compute determinant
      GETRF(&rows, &cols, localmat, &ldm, ipiv, &info);
  
      // If invalid operation skip to the next row
      if (info == 0)
	{
	  // Compute 1-norm
	  norm = anorm(matdim, matdim, localmat);

	  GECON(whichnorm, &cols, localmat, &ldm, &norm, &cond, work, iwork, &info);

	  // If invalid operation skip to the next row
	  if ((info >= 0) && (cond > maxcond+FEMTO))
	    {
	      maxcond = cond;
	      maxind  = i;
	    }
	}
    }

  if (maxind == -1)
    {
      ErrorMsg(__FILE__, __LINE__, "No non-singular matrix found in ComputeEigenvector()");
      memset((void *)eigenvec, 0, matdim*sizeof(double));

      return FAILURE;
    }

  // Make a working copy of the matrix for LU decomposition
  memcpy(localmat, mat, matdim*matdim*sizeof(double));

  // Substract the eigenvalue
  for (j=0; j<matdim; j++) localmat[j*matdim+j] -= eigenval;

  // Replace row maxind of the matrix with only 1's
  // NB: Fortran column-wise storage!!
  for (j=0; j<matdim; j++) localmat[j*matdim+maxind] = 1.0;

  // Make a copy for LU-decomposition
  memcpy(ludmat, localmat, matdim*matdim*sizeof(double));

  GETRF(&rows, &cols, ludmat, &ldm, ipiv, &info);

  // Solve for the eigenvector
  memset((void *)eigenvec, 0, matdim*sizeof(double));
  memset((void *)tv0,      0, matdim*sizeof(double));
  tv0[maxind] = eigenvec[maxind] = 1.0;

  GETRS(trans, &rows, &nrhs, ludmat, &ldm, ipiv, eigenvec, &ldv, &info);
  if (info < 0)
    {
      ErrorMsg(__FILE__, __LINE__,
	       "Illegal value for parameter %d in GETRS", abs(info));
      memset((void *)eigenvec, 0, matdim*sizeof(double));
      return ILLEGAL_INPUT;
    }

  // Store the computed eigenvector
  memcpy(tv1, eigenvec, matdim*sizeof(double));

  // Refine the solution
  GERFS(trans, &rows, &nrhs, localmat, &ldm, ludmat, &ldm, ipiv, tv0, &ldv, 
	eigenvec, &ldv, ferr, berr, work, iwork, &info);
  if (info < 0)
    {
      ErrorMsg(__FILE__, __LINE__,
	       "Illegal value for parameter %d in GERFS", abs(info));

      // Return the unrefined result
      memcpy(eigenvec, tv1, matdim*sizeof(double));

      return SUCCES;
    }

  return SUCCES;
}


/*==============================================================================*/
