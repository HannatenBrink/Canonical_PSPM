/***
   NAME
     io
   DESCRIPTION
     Implements all low-level I/O routines

   Last modification: AMdR - Jan 23, 2014
***/
#ifndef IO
#define IO
#endif
#include <globals.h>

#define   DESCRIP_MAX    80     //* Maximum length of a quantity description
#define   REPORTNOTE_MAX 4096   //* Maximum length of a ReportNote
#define	  MAX_INPUT_LINE 2048	//* Maximum length a line of input


static struct	{
		 char identical_zero[DESCRIP_MAX];
		 char Dead[DESCRIP_MAX];
#if PARAMETER_NR
		 char parameter[PARAMETER_NR][DESCRIP_MAX];
#endif
		} description;

static char	**usernotes;			//* User defined notes


#define ECVF "Unexpected end/error while reading CVF file!"
#define EIZ  "Error during input of zero comparison value from CVF file!"
#define EMIN "Error during input of the cohort minima from the CVF file!"
#define MARN "Memory allocation failure in storing user defined report notes!"
#define REP  "Unable to open REP file for generating report file of the run!"

/*==========================================================================*/
/*
 * Start of function implementations.
 */

void ReportMsg(int reportlevel, const char *fmt, ...)

{
  va_list		argpnt;

  va_start(argpnt, fmt);
  if (Report2Terminal >= reportlevel)
    {
      vfprintf(stderr,  fmt, argpnt);
      (void)fflush(stderr);
    }
  va_end(argpnt);

  va_start(argpnt, fmt);
  if (errfile)
    {
      vfprintf(errfile, fmt, argpnt);
      (void)fflush(errfile);
    }
  va_end(argpnt);
  
  return;
}


 
/*==========================================================================*/

void ErrorMsg(const char *name, const int line, const char *fmt, ...)

{
  va_list		argpnt;

  if (Report2Terminal)
    {
      (void)fprintf(stderr, "\n** %-12s (%3d): ", name, line);
      va_start(argpnt, fmt);
      vfprintf(stderr, fmt, argpnt);
      va_end(argpnt);
      (void)fprintf(stderr, "\n");
      (void)fflush(stderr);
    }

  if (errfile)
    {
      (void)fprintf(errfile, "\n** %-12s (%3d): ", name, line);
      va_start(argpnt, fmt);
      vfprintf(errfile, fmt, argpnt);
      va_end(argpnt);
      (void)fprintf(errfile, "\n");
      (void)fflush(errfile);
    }
  
  return;
}


 
/*==========================================================================*/

void NumProcError(const char *name, const int line, const int NumProcErrorCode)

{
  switch (NumProcErrorCode)
    {
    case MEM_ALLOC_FAIL:
      ErrorMsg(name, line, "%-45s", "Memory allocation failure");
      break;
    case ZERO_DIVISION:
      ErrorMsg(name, line, "%-45s", "Division by 0");
      break;
    case SINGULARITY:
      ErrorMsg(name, line, "%-45s", "Singular matrix encountered");
      break;
    case NORM_OVERFLOW:
      ErrorMsg(name, line, "%-45s", "Norm overflow in Newton iteration");
      break;
    case NO_CONVERGENCE:
      ErrorMsg(name, line, "%-45s", "No convergence in Newton iteration");
      break;
    case ILLEGAL_BOUNDS:
      ErrorMsg(name, line, "%-45s", "Lower and upper value do not bound solution in curvepnt1D()");
      break;
    case ILLEGAL_INPUT:
    case FAILED_EVALUATION:
      break;
    default:
      ErrorMsg(name, line, "%-45s%d", "Unknown numerical error code: ", NumProcErrorCode);
      break;
    }
  
  return ;
}




/*==========================================================================*/

void ErrorExit(const char *name, const int line, const char *fmt, ...)

  /* 
   * ErrorExit - Routine issues an error message. The 
   *		 program subsequently terminates.
   */

{
  va_list		argpnt;

  if (fmt)					/* Error exit               */
    {
      if (Report2Terminal)
        {
	  (void)fprintf(stderr, "\nERROR in %-12s (%3d):\n\n", name, line);
	  (void)fprintf(stderr, "** ");

	  va_start(argpnt, fmt);
	  vfprintf(stderr, fmt, argpnt);
	  va_end(argpnt);

	  (void)fprintf(stderr, " **\n\n");
	  fprintf(stderr, "\a\a\a\a\a");
        }

      if (errfile)
	{
	  (void)fprintf(errfile, "\nERROR in %-12s (%3d):\n", name, line);
	  (void)fprintf(stderr, "** ");

	  va_start(argpnt, fmt);
	  vfprintf(stderr, fmt, argpnt);
	  va_end(argpnt);

	  (void)fprintf(stderr, " **\n\n\n");
	}
#if POPULATION_NR
      if (NewStateComputed) WriteEsfFile();
#endif
#if (DEBUG > 0)
      abort();
#else
      exit(1);
#endif
    }
  else						/* Normal exit              */
    {
#if POPULATION_NR
      if (NewStateComputed) WriteEsfFile();
#endif
#if (DEBUG > 0)
      abort();
#else
      exit(0);
#endif
    }

  return;
}



/*==========================================================================*/

void PrettyPrint(FILE *fp, const int dim, ...)
  
{	
  register int		i;
  double		tmp;
  va_list		pars;

  va_start(pars, dim);
  for (i=0; i<dim; i++)
    {
      tmp = (double)va_arg(pars, double);
      if (((fabs(tmp) <= 1.0E4) && (fabs(tmp) >= 1.0E-3)) || (tmp == 0.))
	{
	  if (i) (void)fprintf(fp, "%15.8f", tmp);
	  else   (void)fprintf(fp, "%13.8f", tmp);
	}
      else
	{
	  if (i) (void)fprintf(fp, "%15.6E", tmp);
	  else   (void)fprintf(fp, "%13.6E", tmp);
	}
    }
  (void)fprintf(fp, "\n");
  (void)fflush(fp);
  va_end(pars);
  
  return;
}




/*==========================================================================*/

void PrettyPrintArray(FILE *fp, const int dim, double *vec)
  
{	
  register int		i;
  double		tmp;

  for (i=0; i<dim; i++)
    {
      tmp = vec[i];
      if (((fabs(tmp) <= 1.0E4) && (fabs(tmp) >= 1.0E-3)) || (tmp == 0.))
	{
	  if (i) (void)fprintf(fp, "%15.8f", tmp);
	  else   (void)fprintf(fp, "%13.8f", tmp);
	}
      else
	{
	  if (i) (void)fprintf(fp, "%15.6E", tmp);
	  else   (void)fprintf(fp, "%13.6E", tmp);
	}
    }
  (void)fprintf(fp, "\n");
  (void)fflush(fp);
  
  return;
}




/*==========================================================================*/

void PrettyPrintHeader(FILE *fp, const int dim, ...)
  
{	
  register int		i, len;
  char			*s, buf[30];
  va_list		pars;

  va_start(pars, dim);
  for (i=0; i<dim; i++)
    {
      s = (char *)va_arg(pars, char *);
      len = strlen(s);
      if (i)
	{
	  strcpy(buf, "                         ");
	  strcpy(buf+15-(13-len)/2-len, s);
	  (void)fprintf(fp, "%-15s", buf);
	}
      else
	{
	  strcpy(buf, "#                        ");
	  strcpy(buf+13-(11-len)/2-len, s);
	  (void)fprintf(fp, "#\n%-13s", buf);
	}
    }
  (void)fprintf(fp, "\n#\n");
  (void)fflush(fp);
  va_end(pars);
  
  return;
}




/*==============================================================================*/

char	*ReadDouble(double *val, char *cpnt)

  /* 
   * ReadDouble - Routine reads a double value from the string pointed to
   *		  by "cpnt". Invalid characters are skipped. It returns a 
   *		  pointer to the rest of the string or NULL on error.
   */

{
  register char		*ch, *end=NULL;
  int			dot_start=0;

  ch=cpnt; while((!isdigit(*ch)) && *ch) ch++;	/* Skip non digits          */
  if(isdigit(*ch))
    {
      end=ch;
      if((ch!=cpnt) && (*(ch-1)=='.'))		/* Is previous a dot?	    */
	{
	  ch--; dot_start=1;
	}
      if((ch!=cpnt) && (*(ch-1)=='-')) ch--;	/* Is previous a minus?	    */

      while(isdigit(*end)) end++;		/* Skip digits		    */
      if((!dot_start) && (*end=='.'))		/* Dot starts mantissa	    */
	{
	  end++; while(isdigit(*end)) end++;
	}

      if((*end=='e') || (*end=='E'))		/* Possible exponent        */
	{
	  end++; if((*end=='+') || (*end=='-')) end++;
	  while(isdigit(*end)) end++;
	}
      *val = 0.0;
      (void)sscanf(ch, "%lg", val);
    }

  return end;
}



/*==========================================================================*/

int	  	ScanLine(FILE *infile, char *descr, double *value, int var_nr)

  /* 
   * ScanLine - Routine reads one line from the control variable file, 
   *		splits it up in the description and the value part and 
   *		returns the number of values read. "var_nr" is the number
   *		of variables to be read.
   */

{
  register char		*ch=NULL, *dsp;
  char			input[MAX_INPUT_LINE], tmp_str[DESCRIP_MAX];
  int			read_no=0;
						/* Input line		    */
  dsp=descr;					/* Skip empty lines	    */
  while((!feof(infile)) && (!ferror(infile)) && (!ch))
    {
      if ((ch=fgets(input, MAX_INPUT_LINE, infile)))
	{
	  while(*ch==' ') ch++;
	  if(*ch == '\n') ch=NULL;
	}
    }
  if(feof(infile) || ferror(infile)) ErrorExit(__FILE__, __LINE__, ECVF);

  if(*ch == '"')				/* Comment is in quotes	    */
    {
      ch++;					/* Skip the first quote	    */
      while(*ch != '"')				/* Search for the closing   */
	{					/* quotes, storing string   */
	  (*dsp)=(*ch); dsp++; ch++;		/* between them		    */
	}
      ch++;					/* Skip closing quotes	    */
    }
  else						/* Comment is not quoted    */
    {
      while((!isdigit(*ch)) && *ch)		/* Skip non-digits and put  */
	{					/* them in description	    */
	  (*dsp)=(*ch); dsp++; ch++;
	}					/* Remove trailing blanks   */
						/* and tabs                 */
      while((*(dsp-1)==' ') || (*(dsp-1)=='\t')) *(--dsp)='\0';
    }

  ch=strcpy(tmp_str, ch);			/* Copy rest to string      */

  while(ch)
      if((ch=ReadDouble(value+read_no, ch)) != NULL)
	  if((++read_no)==var_nr) return read_no;

  return read_no;
}



/*==========================================================================*/

int	ReadDoublesOnLine(FILE *fp, double *y, int n, const int skip)

  /* 
   * ReadDoublesOnLine - Reads n doubles from a single line
   */

{
  char			input[MAX_STR_LEN], *ch;
  int			read_no = 0;
  double 		tmpval;

  if (!fp) return FAILURE;

  while((!feof(fp)) && (!ferror(fp)))
    {						/* Input line, skip empty   */
      ch=fgets(input, MAX_STR_LEN, fp);		/* and lines starting with #*/
      if (skip && (*ch == '#')) continue;

      while(ch)					/* Read double from string, */
	{					/* stop on line end	    */
	  if((ch=ReadDouble(&tmpval, ch)) != NULL)
	    {
	      y[read_no++] = tmpval;
	      if (read_no == n) return SUCCES;
	    }
	}
      if (!skip) break;
    }

  return FAILURE;
}




/*==============================================================================*/

void  SIGINThandler(int sig)
{
  int			c;

  fprintf(stderr, "\n\n%s\n\n%s\n%s\n%s\n%s\n%s\n\n%s",
	  "Choose from the following options:",
	  "1. If you want to quit",
	  "2. If you want to change the step size",
	  "3. If you want to change the maximum step size",
	  "4. If you want to toggle terminal reporting on/off",
	  "5. If you want to continue without change",
	  "Press 1, 2, 3, 4 or 5 : ");
  while ((c < 1) || (c > 5)) scanf("%d", &c);
  fprintf(stderr, "\n");
  switch (c)
    {
    case 1:
      ErrorExit(__FILE__, __LINE__, NULL);
      break;
    case 2:
      if (!Stepchange)
	{
	  fprintf(stderr, "Step size can't be changed now!\n\n\n");
	  return;
	}
      fprintf(stderr, "Change step size from %.5E to : ", parstep);
      (void)scanf("%lg", &parstep);
      longjmp(JumpBuffer, 1);
      break;
    case 3:
      fprintf(stderr, "Change maximum step size from %.5E to : ", Maxparstep);
      (void)scanf("%lg", &Maxparstep);
      break;
    case 4:
      Report2Terminal = !Report2Terminal;
      break;
    default:
      break;
    }

  return;
}




/*==========================================================================*/

void installSIGINThandler()
  
{
  struct sigaction 	sa;
  int			result;

  sa.sa_handler = SIGINThandler;
  sa.sa_flags   = SA_NODEFER;
  result = sigaction(SIGINT, &sa, (struct sigaction *)NULL);
  
  if (!result)
    {
      if (Report2Terminal)
	fprintf(stderr, "Hit <RETURN> or Ctrl-C to change step size!\n");
    }
  else
    ErrorMsg(__FILE__, __LINE__, "Can't catch Ctrl-C interrupt!");

  return;
}




/*==============================================================================*/

void checkESCpressed()
  
{
  register int		j;
  int			success;
  long			bytes = 0L, nread;
  char			buf[MAX_STR_LEN];

  success = ioctl(STDIN_FILENO, FIONREAD, &bytes);
  if ((success != -1) && (bytes > 0))
    {
      nread = min(bytes, MAX_STR_LEN);
      if ((nread = read(STDIN_FILENO, buf, (unsigned)nread)) > 0)
	{
	  for (j=0; j < nread; j++)
	    {
	      if (buf[j] == 10)
		{
		  SIGINThandler(SIGINT);
		  return;
		}
	    }
	}
      else success = -1;
      bytes -= nread;
    }

  return;
}




/*==========================================================================*/

void	  ReadCvfFile(const char *fname)

  //* Routine reads parameter values from the .cvf file.

{
#if POPULATION_NR
  register int		i;
  double		tmp[POPULATION_NR];
#else
  double		tmp[1];
#endif
#if PARAMETER_NR
  register int		j;
  char			msg[MAX_STR_LEN];
#endif
  char			input[MAX_INPUT_LINE];
  int			read_no;
  FILE 			*infile;

  //* Add cvf extension to file name and open
  strcpy(input, fname);
  if (!strstr(input, ".cvf")) (void)strcat(input, ".cvf");
  infile=fopen(input, "r");
  if (!infile) ErrorExit(__FILE__, __LINE__, "Unable to read CVF file: %s !", input);

  read_no=ScanLine(infile, input, tmp, 1);		//* DISCARD accuracy
  read_no=ScanLine(infile, input, tmp, 1);		//* DISCARD cohort_limit

  read_no=ScanLine(infile, description.identical_zero, &equal2zero, 1);
  if(read_no != 1) ErrorExit(__FILE__, __LINE__, EIZ);

  read_no=ScanLine(infile, input, tmp, 1);		//* DISCARD max_time
  read_no=ScanLine(infile, input, tmp, 1);		//* DISCARD delt_out
  read_no=ScanLine(infile, input, tmp, 1);		//* DISCARD state_out

#if POPULATION_NR
  read_no=ScanLine(infile, description.Dead, tmp, POPULATION_NR);
  if(read_no != POPULATION_NR) ErrorExit(__FILE__, __LINE__, EMIN);
  for(i=0; i<POPULATION_NR; i++) logDead[i] = log(tmp[i]);
  
  for(i=0; i<I_STATE_DIM; i++)
    read_no=ScanLine(infile, input, tmp, 1);		//* DISCARD rel_tols[][i]
  for(i=0; i<I_STATE_DIM; i++)
    read_no=ScanLine(infile, input, tmp, 1);		//* DISCARD abs_tols[][i]
#endif

#if PARAMETER_NR
  for(j=0; j<PARAMETER_NR; j++)
    {
      read_no=ScanLine(infile, description.parameter[j], parameter+j, 1);
      if(read_no != 1)
	{
	  (void)sprintf(msg, "Error while reading parameter %d from CVF file!",
			j);
	  ErrorExit(__FILE__, __LINE__, msg);
	}
    }
#endif

  (void)fclose(infile);

  return;
}



/*==========================================================================*/

void ReportNote(const char *fmt, ...)

{
  va_list		argpnt;
  register int		i;
  static int		lines = 0, count = 0;
  
  if ((!usernotes) || (count >= (lines-1)))
    {
      lines += 10;
      usernotes = (char **)realloc(usernotes, (size_t)lines*sizeof(char *));
      if (!usernotes) ErrorExit(__FILE__, __LINE__, MARN);
      for (i=count; i<lines; i++) usernotes[i] = (char *)NULL;
    }
  usernotes[count] = (char *)calloc(REPORTNOTE_MAX, sizeof(char));
  if (!(usernotes[count])) ErrorExit(__FILE__, __LINE__, MARN);

  va_start(argpnt, fmt);
  vsnprintf(usernotes[count], REPORTNOTE_MAX, fmt, argpnt);
  va_end(argpnt);
  
  count++;

  return;
} /* ReportNotes */



/*==========================================================================*/

void	  WriteReport(const char *fname, int argc, char **argv)

{
  register int		i;
  char			filename[DESCRIP_MAX], *ch, **unp;
  FILE			*rep;

  //* Open REP file

  ch=strcpy(filename, fname);
  if (!strstr(filename, ".rep")) ch=strcat(filename, ".rep");
  rep=fopen(filename, "w");
  if(!rep) ErrorExit(__FILE__, __LINE__, REP);
  
  while(*ch) ch++;
  while(*(--ch)!='.') *ch=' '; *ch='\0';
  
  for(i=0; i<79; i++) (void)fprintf(rep, "*"); (void)fprintf(rep, "\n");
  (void)fprintf(rep, "\n%2s%-s\n", " ", "PROGRAM, RUN AND ARGUMENTS");
  (void)fprintf(rep, "%4sProgram   : %-s\n", " ", argv[0]);
  (void)fprintf(rep, "%4sRun       : %-s\n", " ", fname);
  (void)fprintf(rep, "%4sArguments : ", " ");
  for (i=1; i<argc; i++) (void)fprintf(rep, " %-s", argv[i]);
  (void)fprintf(rep, "\n");

  if (usernotes)
    {
      unp = usernotes;
      (void)fprintf(rep, "\n%2s%-s\n", " ", "MODEL SPECIFIC NOTES");
      while (*unp)
	{
	  (void)fprintf(rep, "%4s%-s\n", " ", *unp);
	  unp++;
	}
    }
  
  (void)fprintf(rep, "\n%2s%-s\n", " ", "USED VALUES FOR CONTROL VARIABLES");
  (void)fprintf(rep, "%4s%-65s%5s%-10.4G\n", " ", 
		description.identical_zero, "  :  ", equal2zero);

#if POPULATION_NR
  (void)fprintf(rep, "\n%2s%-s\n", " ", "USED VALUES OF TOLERANCES");

  (void)fprintf(rep, "%4s%-65s%5s", " ",
		description.Dead, "  :  ");
  for(i=0; i<POPULATION_NR; i++)
      (void)fprintf(rep, "%-7.4G", exp(logDead[i]));
  (void)fprintf(rep, "\n");
#endif

#if PARAMETER_NR
  (void)fprintf(rep, "\n%2s%-s\n", " ", "USED VALUES OF PARAMETERS");
  
  for(i=0; i<PARAMETER_NR; i++)
      (void)fprintf(rep, "%4s%-65s%5s%-10.4G\n", " ", 
		    description.parameter[i], "  :  ", parameter[i]);
#endif

  for(i=0; i<79; i++) (void)fprintf(rep, "*"); (void)fprintf(rep, "\n");
  
  (void)fclose(rep);
  
  return;
}



/*==========================================================================*/



      
