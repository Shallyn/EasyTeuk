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
#define DEFAULT_STEP_RATIO 0.9
#define DEFAULT_TIME_MAX 500.

typedef struct tagCoreParams {
    REAL8 mass;
    REAL8 spin;
    REAL8 eta;
    REAL8 step_ratio;
    REAL8 t_max;
    INT shot_step;
    INT source_type;
    INT mode;
    // dump
    REAL8 rt0;
    REAL8 theta0;
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

