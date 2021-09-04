/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_ETSTRUCTS__
#define __INCLUDE_ETSTRUCTS__

#include "myUtils.h"
#include "myLog.h"
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_spline.h>

#define DEFAULT_MASS 1e6
#define DEFAULT_SPIN 0.0
#define DEFAULT_ETA 1e-3

typedef struct tagCoreParams {
    REAL8 mass;                 /**< mass of companion 1 */
    REAL8 spin;
    REAL8 eta;
    REAL8 step_ratio;
    INT mode;
    char prefix[STR_COMM_SIZE];
} CoreParams;

typedef struct tagRK4Integrator
{
    gsl_odeiv_step    *step;
    gsl_odeiv_control *control;
    gsl_odeiv_evolve  *evolve;

    gsl_odeiv_system  *sys;

    int (* dydt) (double t, const double y[], double dydt[], void * params);
    int (* stop) (double t, const double y[], double dydt[], void * params);

    int retries;		/* retries with smaller step when derivatives encounter singularity */
    int stopontestonly;	/* stop only on test, use tend to size buffers only */

    int returncode;
} RK4Integrator;


#endif

