/***
   NAME
     cohort
   DESCRIPTION
     Implements routines to carry out population and cohort manipulations.

   Last modification: AMdR - Jan 29, 2010
***/
#ifndef COHORT
#define COHORT
#endif
#include <globals.h>



/*
 *===========================================================================
 * Some local type definitions
 *===========================================================================
 */

// Magic key of the type of CSB file written by this module

#define		CSB_MAGIC_KEY		20030509

typedef struct envdim {
			double		timeval;
			int		columns;
			int		data_offset;
			uint32_t	memory_used;
		       } Envdim;

typedef struct popdim {
			double		timeval;
			int		population;
			int		cohorts;
			int		columns;
			int		data_offset;
			int		lastpopdim;
		       } Popdim;

/*
 *===========================================================================
 * Some global variables
 *===========================================================================
 */

static double		*CohortSpace = NULL;
static int		CohortSize;
static size_t		DataMemAllocated = 0;
static char     	statelabels[POPULATION_NR][MAX_STR_LEN];




/*==========================================================================*/
/*
 * Start of function implementations.
 */

void SetCohortSize(const int newsize)

{
  register int		i, j;
  size_t	     	mem_req;

  CohortSize = newsize;
  mem_req    = (POPULATION_NR*MAX_COHORTS)*CohortSize;

  if (mem_req > DataMemAllocated)		/* Create memory for new    */
    {						/* cohort size if necessary */
      DataMemAllocated = mem_req;
      if (CohortSpace)
	CohortSpace = (double *)realloc((void *)CohortSpace,
					DataMemAllocated*sizeof(double));
      else
	CohortSpace = (double *)malloc(DataMemAllocated*sizeof(double));

      if (!CohortSpace)
	ErrorExit(__FILE__, __LINE__, "Memory allocation failure for cohort space!");
    }
  
  /*
   * Set up the cohort pointers
   */
  for (i=0; i<POPULATION_NR; i++)
    for (j=0; j<MAX_COHORTS; j++)
      pop[i][j] = CohortSpace + (i*MAX_COHORTS + j)*CohortSize;
  
  return;
}


 
/*==========================================================================*/

void ReadIsfFile(const char *fname, const int restart)

  /* 
   * ReadIsfFile - Routine reads all values from the named .isf file.
   */

{
  register int		i, j;
  char			*ch, input[MAX_STR_LEN];
  int			done, read_no, first, cohnr, warnics;
  double		tmpval;
  FILE			*infile;

  //* Add isf extension to file name if required and open the file
  strcpy(input, fname);
  if (restart) (void)strcat(input, ".esf");
  else (void)strcat(input, ".isf");
  infile=fopen(input, "r");
  if (!infile) ErrorExit(__FILE__, __LINE__, "Unable to read ISF or ESF file!");

  //* Read the environmental variables
  read_no=0;
  while((!feof(infile)) && (!ferror(infile)))
    {
      ch=fgets(input, MAX_STR_LEN, infile);
      if (!ch || (*ch == '#')) continue;

      while(ch)					//* Read double from string
	{					//* stop on line end
	  if((ch=ReadDouble(&tmpval, ch)) != NULL)
	    {
	      env[read_no++] = tmpval;
	    }
	  if(read_no == ENVIRON_DIM) break;	//* Ignore additional values
	}
      if (read_no) break;			//* Stop if read some data
    }

  //* Read the population data
  for(i=0, first=1; i<POPULATION_NR; i++)
    {
      done    = 0;
      cohnr   = 0;
      warnics = 1;
      while((!feof(infile)) && (!ferror(infile)) && (!done))
	{
	  ch=fgets(input, MAX_STR_LEN, infile);
	  if (!ch || (*ch == '#')) continue;

	  //* If first cohort, count number on line and set CohortSize
	  if (first)
	    {
	      read_no = 0;
	      while(ch)				/* Count doubles in string, */
		{				/* stop on line end	    */
		  if((ch=ReadDouble(&tmpval, ch)) != NULL) read_no++;
		}
	      if (read_no)
		{
		  SetCohortSize(read_no);
		  first = 0;
		  ch = input;
		}
	      else continue;			//* Read an empty start line
	    }
	  
	  //* Initialize all data to default: MISSING_VALUE
	  for(j=0; j<CohortSize; j++) pop[i][cohnr][j] = MISSING_VALUE;

	  /*
	   * Read complete cohort from one line.
	   * Stop on line end or full cohort.
	   */
	  for(j=0, read_no=0; (j<CohortSize) && (ch); j++)
	    if((ch=ReadDouble(pop[i][cohnr]+j, ch)) != NULL) read_no++;
					
	  //* End of population reached
          if ((read_no == 0) && cohnr) done = 1;
	  else if (read_no > 0)
	    {
	      //* Signal incomplete cohort
	      if ((read_no != CohortSize) && warnics)
		{
		  ErrorMsg(__FILE__, __LINE__,
			   "%-45s", "Incomplete cohort(s) in population %d", i);
		  warnics =0;
		}
	      cohnr++;
	    }
	}
      cohort_no[i] = cohnr;
    }

  for(i=0; i<POPULATION_NR; i++)
    if(!cohort_no[i])
      ErrorMsg(__FILE__, __LINE__, "Initial population %d is empty!", i);

  return;
}



/*==========================================================================*/

void LabelState(int popnr, const char *fmt, ...)

{
  va_list		argpnt;
  char			tmpbuf[MAX_STR_LEN];

  va_start(argpnt, fmt);
  vsprintf(tmpbuf, fmt, argpnt);
  va_end(argpnt);
  strncpy(statelabels[popnr], tmpbuf, (MAX_STR_LEN-1)*sizeof(char));
  statelabels[popnr][MAX_STR_LEN-1] = '\0';

  return;
}



/*==========================================================================*/

void WriteBinStateToFile(const char *fname,
			 void (*updatepop)(double e[ENVIRON_DIM],
					   cohort p[POPULATION_NR][MAX_COHORTS]))

/* 
 * WriteBinStateToFile - Routine writes the entire state of the populations 
 *		         and the environment in binary format to the named
 *			 file (with csb extension) after updating the
 *			 population state using the procedure (*updatepop)().
 */
  
{
  register int		i, j, k;
  size_t		hdrdbls, lbldbls;
  Envdim		cenv[2];
  Popdim		cpop[2];
  double		zero = 0.0;
  double		tmpval;
  uint32_t		tmpuint32;
  int			tmpint;
  char			input[MAX_STR_LEN];
  FILE			*fp;
  struct stat           st;

  //* Add csb extension to file name if required and open the file
  strcpy(input, fname);
  if (strcmp(input+strlen(input)-4, ".csb")) (void)strcat(input, ".csb");

  if ((stat(input, &st) != 0) || (st.st_size == 0))
    {
      // File does not exist yet. Write magic key and parameters
      fp=fopen(input, "wb");
      if (!fp) ErrorExit(__FILE__, __LINE__, "Unable to open CSB file!");

      tmpuint32 = CSB_MAGIC_KEY;
      fwrite((void *)(&tmpuint32), 1, sizeof(uint32_t), fp);
      tmpint  = PARAMETER_NR;
      fwrite((void *)(&tmpint), 1, sizeof(int), fp);
      fwrite((void *)(parameter), PARAMETER_NR, sizeof(double), fp);
    }
  else
    {
      fp=fopen(input, "ab");
      if (!fp) ErrorExit(__FILE__, __LINE__, "Unable to open CSB file!");
    }

  //* NB: The user is responsible for setting cohort_no[] appropriately!
  if (updatepop) (*updatepop)(env, pop);
  NewStateComputed = 1;

  hdrdbls  = (sizeof(Envdim)/sizeof(double))+1;
  (void)memset((void *)cenv, 0, hdrdbls*sizeof(double));

  cenv->timeval      = env[0];
  cenv->columns      = ENVIRON_DIM;
  cenv->data_offset  = hdrdbls;
  cenv->memory_used  = hdrdbls*sizeof(double);
  cenv->memory_used += ENVIRON_DIM*sizeof(double);
  for (i=0; i<POPULATION_NR; i++)
    {
      hdrdbls = (sizeof(Popdim)/sizeof(double))+1;
      lbldbls = (strlen(statelabels[i])*sizeof(char))/sizeof(double)+1;
      cenv->memory_used += (hdrdbls+lbldbls)*sizeof(double);
      cenv->memory_used += imax(cohort_no[i],1)*CohortSize*sizeof(double);
    }

  hdrdbls  = (sizeof(Envdim)/sizeof(double))+1;
  fwrite((void *)cenv, 1, hdrdbls*sizeof(double), fp);
  for (i=0; i<ENVIRON_DIM; i++)
    {
      tmpval = (double)env[i];
      fwrite((void *)(&tmpval), 1, sizeof(double), fp);
    }

  for (i=0; i<POPULATION_NR; i++)
    {
      hdrdbls = (sizeof(Popdim)/sizeof(double))+1;
      lbldbls = (strlen(statelabels[i])*sizeof(char))/sizeof(double)+1;
      (void)memset((void *)cpop, 0, hdrdbls*sizeof(double));
      cpop->timeval     = env[0];
      cpop->population  = i;
      cpop->columns     = CohortSize;
      cpop->cohorts     = imax(cohort_no[i], 1);
      cpop->data_offset = hdrdbls+lbldbls;
      cpop->lastpopdim  = (i == (POPULATION_NR-1));

      fwrite((void *)cpop, 1, hdrdbls*sizeof(double), fp);
      fwrite((void *)statelabels[i], 1, lbldbls*sizeof(double), fp);

      if (cohort_no[i])
	{
	  for(k=0; k<CohortSize; k++) 
	    for(j=0; j<cohort_no[i]; j++)
	      {
		tmpval = pop[i][j][k];
		fwrite((void *)(&tmpval), 1, sizeof(double), fp);
	      }
	}
      else
	{
	  for(k=0; k<CohortSize; k++) 
	    fwrite((void *)(&zero), 1, sizeof(double), fp);
	}
    }

  (void)fflush(fp);				//* Flush the file buffer
  (void)fclose(fp);
  
  return;
}



/*==========================================================================*/

void	  WriteEsfFile()

/* 
 * WriteEsfFile - Routine writes the entire state of the populations
 *		  and the environment to an ESF file
 */
  
{
  register int		i, j, k;
  register double	output;
  char			input[MAX_STR_LEN];
  FILE			*fp;

  //* Create ESF filename and open it
  sprintf(input, "state-%G.esf", env[0]);
  fp=fopen(input, "w");
  if (!fp) ErrorExit(__FILE__, __LINE__, "Unable to open ESF file!");

  for(i=0; i<ENVIRON_DIM; i++)			/* Write environment state  */
    {
      if (i > 0) (void)fprintf(fp, "\t");
      output = env[i];
      if (((fabs(output) <= 1.0E4) && (fabs(output) >= 1.0E-4)) ||
	  (output == 0))
	(void)fprintf(fp, "%f", output);
      else (void)fprintf(fp, "%E", output);
    }
  (void)fprintf(fp, "\n\n");
						/* Write population state   */
  for(i=0; i<POPULATION_NR; i++, (void)fprintf(fp, "\n"))
    {
      if (cohort_no[i])
	{
	  for(j=0; j<cohort_no[i]; j++, (void)fprintf(fp, "\n"))
	    {
	      for(k=0; k<CohortSize; k++) 
		{
		  if (k > 0) (void)fprintf(fp, "\t");
		  output = pop[i][j][k];
		  if (((fabs(output) <= 1.0E4) && (fabs(output) >= 1.0E-4)) ||
		      (output == 0)) (void)fprintf(fp, "%f", output);
		  else (void)fprintf(fp, "%E", output);
		}
	    }
	}
      else
	{
	  for(j=0; j<CohortSize; j++)
	    {
	      if (j) (void)fprintf(fp, "\t");
	      (void)fprintf(fp, "%f", 0.0);
	    }
	  (void)fprintf(fp, "\n");
	}
    }
  (void)fprintf(fp, "\n");

  return;
}




/*==========================================================================*/


