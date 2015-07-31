/***
 NAME
 EvoSym
 
 PURPOSE
 Module computes the internal equilibrium of a consumer population
 
 STARTING POINTS:
 
 Last modification: AMdR - Mar 07, 2014
 ***/

#include "globals.h"

/*
 *===========================================================================
 *		Definition of global variables and parameters
 *===========================================================================
 */

// These are the variables to solve for
static double		Requi, Bequi;
static int          Bifparone = 0;
static int          Bifpartwo = 0;
static int          Bifparthree = 0;
static int          Bifpartfour = 0;

// Dimensioning
#define VOLUME		1.0 // 1.0E+6

// These are other variables or parameters occurring in the equations
#define DELTA		parameter[ 0]	// Default: 0.1
#define RMAX		parameter[ 1]	// Default: 10

#define XB          parameter[ 2]	// Default: 0.1
#define XJ          parameter[ 3]	// Default: 1.0
#define XM          parameter[ 4]	// Default: 10.0

#define M           parameter[ 5]	// Default: 1.0
#define T           parameter[ 6]	// Default: 0.1
#define MH          parameter[ 7]	// Default: 0.333333

#define	SIGMA		parameter[ 8]	// Default: 0.5

#define MUC         parameter[ 9]	// Default: 0.01
#define MUJ         parameter[10]	// Default: 0.0
#define MUA         parameter[11]	// Default: 0.0
#define Q           parameter[12]	// Default: 1.0
#define P           parameter[13]	// Default: 1.0




/*
 *===========================================================================
 *	Implementation of user-defined routines
 *===========================================================================
 */

int             SetPntDim()

{
    // Set the dimension of the solution point to be solved for
    
    return BASE_PNTDIM;
}



/*===========================================================================*/

void		UserInit(int argc, char *argv[])

{
    Bifparone = BIFPARONE;
    parameter[Bifparone] = atof(argv[2]);
    
    Bifpartwo            = BIFPARTWO;
    Bifparthree          = BIFPARTHREE;
    Bifpartfour          = BIFPARFOUR;
#if ((SYSTEM_TYPE == SYSTEM_EVOSYM_BP) || (SYSTEM_TYPE == SYSTEM_EVOSYM_LP))
    parameter[Bifpartwo] = atof(argv[BASE_PNTDIM+1]);
#endif
    
#if (SYSTEM_TYPE == SYSTEM_EVOSYM_ESS)
    parameter[Bifpartwo] = atof(argv[BASE_PNTDIM]);
    //parameter[Bifparthree] = atof(argv[BASE_PNTDIM+1]);
#endif
    
    ReportNote("Consumer output per volume of: %G L", VOLUME);
    
    return;
}




/*==============================================================================*/

// The ODE system: definition of variables
#define	SIZE		x[0]
#define	SURVIVAL	x[1]
#define	CUMREPRO	x[2]
#define	CUMINTAKE	x[3]

#define	JUVENILE	0
#define	ADULT		1
#define STAGES		2

#define BASEODE_DIM     4
#define STAGEVAR_DIM	5
#define COHORT_DIM      0

#define ODEDIM		(BASEODE_DIM+STAGEVAR_DIM+COHORT_DIM)

// Global variables to hold variables shared among routines
static int          LifeStage, CurSysDim;
static double		Lbound[STAGES];
static double		SavedState[STAGES][BASEODE_DIM+STAGEVAR_DIM];

// The ODE system: define the right-hand side
void Lifehistory (double age, double *x, double *dxda)
{
    int             index;
    double          Ingest, netproduction, fecundity, growth, mort;
    double          L, kappa;
    double          Ingest_R;
    
    Ingest          = M*pow(SIZE, Q) * Requi / (Requi + 1);
    Ingest_R        = MH*pow(SIZE, Q) * Requi / (Requi + 1);
    netproduction   = SIGMA*Ingest - T*pow(SIZE, P);
    
    if (LifeStage == JUVENILE)
    {
        fecundity = 0.0;
        growth    = netproduction;
        mort      = MUC + MUJ;
    }
    else
    {
        L         = (SIZE - XJ)/(XM - XJ);
        kappa     = 1 - 3*L*L + 2*L*L*L;
        growth    = kappa*netproduction;
        fecundity = (1-kappa)*netproduction/XB;
        mort      = MUC + MUA;
    }
    
    
    if (!dxda) return;
    
    index = 0;
    dxda[index++]	= growth;
    dxda[index++]	= -mort*SURVIVAL;
    dxda[index++]	= fecundity*SURVIVAL;
    dxda[index++]	= Ingest_R*SURVIVAL;
    
    // Do not set the values below if integration is within iteration loop
    if (CurSysDim < ODEDIM) return;
    
    // Lifestage variables
    dxda[index++]	= SURVIVAL;                         // Number in size class
    dxda[index++]	= SIZE*SURVIVAL;					// Biomass in size class
    dxda[index++]	= growth*SURVIVAL;					// Growth in (structural) biomass
    dxda[index++]	= XB*fecundity*SURVIVAL;			// Increase in biomass through reproduction
    dxda[index++]	= mort*SIZE*SURVIVAL;				// Loss of biomass through mortality
    
    return;
}


// The ODE system: define possible stopping criteria

double Stop(double age, double *x)
{
    double	result = -1.0;
    
    result = SIZE - Lbound[LifeStage];
    result = max(result, logDead[0] - log(SURVIVAL));
    
    return result;
}


/*===========================================================================*/

int		Equation(double *argument, double *result)

{
    int			index, i, j;
    double		tval, tend, x[ODEDIM];
    
    index = 0;
    //******************************************************************************
    // Map current estimate of solution to global variables
    
#if (SYSTEM_TYPE == SYSTEM_EVOSYM_BP)
    parameter[Bifparone] = argument[index]*pnt_scale[index]; index++;
    parameter[Bifpartwo] = argument[index]*pnt_scale[index]; index++;
    Requi = RMAX;
#elif (SYSTEM_TYPE == SYSTEM_EVOSYM_EQ)
    parameter[Bifparone] = argument[index]*pnt_scale[index]; index++;
    Requi                = argument[index]*pnt_scale[index]; index++;
#elif (SYSTEM_TYPE == SYSTEM_EVOSYM_LP)
    parameter[Bifparone] = argument[index]*pnt_scale[index]; index++;
    Requi                = argument[index]*pnt_scale[index]; index++;
    parameter[Bifpartwo] = argument[index]*pnt_scale[index]; index++;
#elif (SYSTEM_TYPE == SYSTEM_EVOSYM_ESS)
    parameter[Bifparone]    = argument[index]*pnt_scale[index]; index++;
    Requi                   = argument[index]*pnt_scale[index]; index++;
    if (NOR0XX == 1) {
    parameter[Bifpartwo]    = argument[index]*pnt_scale[index]; index++;
    }
    //parameter[Bifparthree]  = argument[index]*pnt_scale[index]; index++;
#elif (SYSTEM_TYPE == SYSTEM_EVOSYM_CAN)
    parameter[Bifparone] = argument[index]*pnt_scale[index]; index++;
    Requi                = argument[index]*pnt_scale[index]; index++;
#endif
    
    Lbound[JUVENILE] = XJ;
    Lbound[ADULT]    = XM;
    
    //==============================================================================
    // Set the initial point for the ODEs
    
    for (i=0; i<ODEDIM; i++) x[i] = 0.0;
    SIZE      = XB;
    SURVIVAL  = 1.0;
    
    tval      = 0.0;
    tend	  = -logDead[0]/MUC;
    
    if (result)
    {
        // During iteration integrate reduced (base) ODE system only
        CurSysDim   = BASEODE_DIM;
    }
    else
    {
        // For output integrate the extended ODE system
        CurSysDim   = ODEDIM;
    }
    
    // Integration is carried out in 2 consecutive steps, for the juvenile and adult phase separately
    LifeStage   = JUVENILE;
    
    // The integration loop
    while (1)
    {
        if (odesolve(x, CurSysDim, &tval, tend, Lifehistory, Stop) == FAILURE)
        {
            ErrorMsg(__FILE__, __LINE__, "Integration failed!");
            return FAILURE;
        }
        
        if (fabs(SIZE-Lbound[LifeStage]) < MICRO)
        {
            if (!result)					// Only store values if routine is invoked for output
            {
                for (i=0;           i<(BASEODE_DIM+STAGEVAR_DIM); i++) SavedState[LifeStage][i] = x[i];
                for (i=BASEODE_DIM; i<(BASEODE_DIM+STAGEVAR_DIM); i++) x[i] = 0.0;
            }
            LifeStage++;
        }
        else if (SIZE >= (Lbound[LifeStage] - MICRO))
        {
            ErrorMsg(__FILE__, __LINE__, "Integration failed to stop at L = %.2f! (L = %.2f)", Lbound[LifeStage], SIZE);
            return FAILURE;
        }
        
        if ((-log(SURVIVAL) >= (-logDead[0]-MICRO)) || (tval >= (tend-MICRO))) break;
    }
    
    // Store data of the final stage & cohort
    if (!result)					// Only compute values if routine is invoked for output
    {
        for (i=0; i<(BASEODE_DIM+STAGEVAR_DIM); i++) SavedState[ADULT][i] = x[i];
        
        Bequi = VOLUME*DELTA*(RMAX - Requi)/CUMINTAKE;
        
        for (i=0; i<STAGES; i++)
            for (j=BASEODE_DIM; j<(BASEODE_DIM+STAGEVAR_DIM); j++)
                SavedState[i][j] *= Bequi;
        
        return SUCCES;
    }
    
    //==============================================================================
    // Compute the final values of the fixed point equation F(y)=0,
    
#if (SYSTEM_TYPE == SYSTEM_EVOSYM_BP)
    result[0]  = CUMREPRO - 1.0;
#elif (SYSTEM_TYPE == SYSTEM_EVOSYM_EQ)
    result[0]  = CUMREPRO - 1.0;
#elif (SYSTEM_TYPE == SYSTEM_EVOSYM_CAN)
    result[0]  = CUMREPRO - 1.0;
#elif (SYSTEM_TYPE == SYSTEM_EVOSYM_LP)
    result[0]  = CUMREPRO - 1.0;
    result[1]  = LPcondition(BASE_PNTDIM, argument, Equation, CENTRAL);
#elif (SYSTEM_TYPE == SYSTEM_EVOSYM_ESS)
    result[0]  = CUMREPRO - 1.0;
    // The last index in the following call indicates the index of the equation R0-1= 0 in the result[] vector
    if (NOR0XX == 1)
    {
    result[1]  = SelectionGradient(BASE_PNTDIM, argument, Equation, 0, 0);
    }
   
    
#endif
    
    return SUCCES;
}



/*
 *==========================================================================
 *
 *			UPDATE IDCARDS
 *
 *==========================================================================
 */

void UpdateCohorts(double *env, cohort pop[POPULATION_NR][MAX_COHORTS])

{
    return;
}


/*==============================================================================*/

int		DefineOutput(double *argument, double *output)

{
    int			outnr = 0;
    
    Equation(argument, NULL);
    
    // There are maximally 4 values in the point vector to solve for, which
    // occurs when continuing an LP for this system. In all cases the 4
    // values are written to the output file

    output[outnr++] = parameter[Bifparone];
    output[outnr++] = SelectionGradient(BASE_PNTDIM, argument, Equation, 0, 0);
    output[outnr++] = parameter[Bifpartwo];
    output[outnr++] = SelectionGradient2(BASE_PNTDIM, argument, Equation, 13, 0);
    output[outnr++] = Requi;
    output[outnr++] = Bequi;
#if (SYSTEM_TYPE == SYSTEM_EVOSYM_ESS)
    output[outnr++] = R0yy(BASE_PNTDIM,argument,Equation,0,0);
    output[outnr++] = R0xx(BASE_PNTDIM-1,argument,Equation,0,0);
#endif

 
    
    
    return outnr;
}



/*==============================================================================*/
