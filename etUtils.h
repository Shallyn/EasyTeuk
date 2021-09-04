/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_ETUTILS__
#define __INCLUDE_ETUTILS__

#include "etStructs.h"
#include "myFileIO.h"
#define BEGINGSL \
        { \
          gsl_error_handler_t *__UNIQUE_ID(saveGSLErrorHandler); \
          __UNIQUE_ID(saveGSLErrorHandler) = gsl_set_error_handler_off();

#define ENDGSL \
          gsl_set_error_handler( __UNIQUE_ID(saveGSLErrorHandler) ); \
        }


INT get_DebugFlag();
void set_DebugFlag(INT flag);
#define IS_DEBUG (get_DebugFlag())
#define DEBUG_START \
do{set_DebugFlag(1);}while(0)
#define DEBUG_END \
do{set_DebugFlag(0);}while(0)

INT get_VersionFlag();
void set_VersionFlag(INT flag);
#define CODE_VERSION (get_VersionFlag())
#define SET_CODE_VERSION(ver) \
do{set_VersionFlag(ver);}while(0)

REAL8 numerical_discrete_delta(REAL8 x0, REAL8 index, REAL8 val, REAL8 dx);

RK4Integrator *CreateRK4Integrator(int dim, 
    int (*dydt) (double t, const double y[], double dydt[], void *params),   /* These are XLAL functions! */
    int (*stop) (double t, const double y[], double dydt[], void *params), 
    double eps_abs, double eps_rel);
void DestroyRK4Integrator(RK4Integrator *integrator);


#endif

