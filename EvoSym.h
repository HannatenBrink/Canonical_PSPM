// - header file -- AMdR - Jan 23, 2014

// System definition

#define SYSTEM_EVOSYM_BP		100
#define SYSTEM_EVOSYM_EQ		101
#define SYSTEM_EVOSYM_LP		102
#define SYSTEM_EVOSYM_ESS		103
#define SYSTEM_EVOSYM_CAN		104

#ifndef SYSTEM_TYPE
#define SYSTEM_TYPE		SYSTEM_EVOSYM_EQ
#endif

// Dimension of the solution point
#if (SYSTEM_TYPE == SYSTEM_EVOSYM_BP)
#define BASE_PNTDIM		2
#define	PROBLEMVARS		"Par.1 Par.2"
#elif   (SYSTEM_TYPE == SYSTEM_EVOSYM_EQ)
#define BASE_PNTDIM		2
#define	PROBLEMVARS		"Par.1 R"
#elif (SYSTEM_TYPE == SYSTEM_EVOSYM_LP)
#define BASE_PNTDIM		3
#define	PROBLEMVARS		"Par.1 R Par.2"
#elif (SYSTEM_TYPE == SYSTEM_EVOSYM_ESS)
#define BASE_PNTDIM		3
#define	PROBLEMVARS		"Par.1 R Par.2"
#elif   (SYSTEM_TYPE == SYSTEM_EVOSYM_CAN)
#define BASE_PNTDIM		2
#define	PROBLEMVARS		"Par.1 R"
#define     CANONICAL     	1
#endif // (SYSTEM_TYPE == SYSTEM_EVOSYM_BP)




#define PARSTEPSCALING		1

// Dimension of the structured population for output, parameter number and maximum output length
// NOT USED AT PRESENT
#define POPULATION_NR   	1
#define	ENVIRON_DIM         2
#define	I_STATE_DIM         2
#define	I_CONST_DIM         5

#define	PARAMETER_NR		18
#define	MAX_PNTDIM          6
#define MAX_OUTPUTDIM		100

// Numerical settings
#define DYTOL			1.0E-7		// Variable tolerance
#define RHSTOL			1.0E-7		// Function tolerance
#define CANTOL          1.0E-6      // Tolerance for selectiongradient=0
#define DELTAS           1.0E-4      // Factor to multiply with selection gradient
#define EVODIM              2           //number of evolving traits

#undef 	JACOBIAN_STEP
#define JACOBIAN_STEP   	1.0E-3  	// % change for jacobian
#undef 	UPDATE_JAC
#define	UPDATE_JAC          20
#define MOOREPENROSE		0
#define SMALLEST_STEP		1.0E-8
#define CENTRALDIFF         1
