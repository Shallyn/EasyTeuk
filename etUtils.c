/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include "etUtils.h"

INT g_DebugFlag = 0;
INT get_DebugFlag()
{
    return g_DebugFlag;
}
void set_DebugFlag(INT flag)
{
    g_DebugFlag = flag;
    return;
}

INT g_VersionFlag = 1;
INT get_VersionFlag()
{
    return g_VersionFlag;
}
void set_VersionFlag(INT flag)
{
    g_VersionFlag = flag;
    return;
}

/**
 * @brief Numerical Discrete delta function
 * @param 
 *  REAL8 xi
 *  REAL8 val
 *  REAL8 dx
 */
REAL8 numerical_discrete_delta(REAL8 x0, REAL8 index, REAL8 val, REAL8 dx)
{
    INT idx_val = (val / dx);
    REAL8 h4 = pow(dx, 4.);
    if (index < idx_val - 1)
        return 0.0;
    if (index > idx_val + 2)
        return 0.0;
    REAL8 x_k = (idx_val-1) * dx;
    REAL8 x_kp1 = x_k + dx;
    REAL8 x_kp2 = x_kp1 + dx;
    REAL8 x_kp3 = x_kp2 + dx;
    if (index == idx_val - 1)
        return -(val - x_kp1) * (val - x_kp2) * (val - x_kp3) / 6. / h4;
    if (index == idx_val)
        return (val - x_k) * (val - x_kp2) * (val - x_kp3) / 2. / h4;
    if (index == idx_val+1)
        return -(val - x_k) * (val - x_kp1) * (val - x_kp3) / 2. / h4;
    if (index == idx_val+2)
        return (val - x_k) * (val - x_kp1) * (val - x_kp2) / 6. / h4;
    return 0.0;
}

RK4Integrator *CreateRK4Integrator(int dim, 
    int (*dydt) (double t, const double y[], double dydt[], void *params),   /* These are XLAL functions! */
    int (*stop) (double t, const double y[], double dydt[], void *params), 
    double eps_abs, double eps_rel)
{
    RK4Integrator *integrator;

    /* allocate our custom integrator structure */
    if (!(integrator = (RK4Integrator *) MYCalloc(1, sizeof(RK4Integrator)))) 
    {
        return NULL;
    }

    /* allocate the GSL ODE components */
    integrator->step = gsl_odeiv_step_alloc(gsl_odeiv_step_rkf45, dim);
    //XLAL_CALLGSL(integrator->step = gsl_odeiv_step_alloc(gsl_odeiv_step_rk8pd, dim));
    integrator->control = gsl_odeiv_control_y_new(eps_abs, eps_rel);
    integrator->evolve = gsl_odeiv_evolve_alloc(dim);

    /* allocate the GSL system (functions, etc.) */
    integrator->sys = (gsl_odeiv_system *) MYCalloc(1, sizeof(gsl_odeiv_system));

    /* if something failed to be allocated, bail out */
    if (!(integrator->step) || !(integrator->control) || !(integrator->evolve) || !(integrator->sys)) 
    {
        STRUCTFREE(integrator, RK4Integrator);
        return NULL;
    }

    integrator->dydt = dydt;
    integrator->stop = stop;

    integrator->sys->function = dydt;
    integrator->sys->jacobian = NULL;
    integrator->sys->dimension = dim;
    integrator->sys->params = NULL;

    integrator->retries = 6;
    integrator->stopontestonly = 0;

    return integrator;
}

void DestroyRK4Integrator(RK4Integrator *integrator)
{
    if (!integrator)
        return;
    if (integrator->evolve)
        gsl_odeiv_evolve_free(integrator->evolve);
    if (integrator->control)
        gsl_odeiv_control_free(integrator->control);
    if (integrator->step)
        gsl_odeiv_step_free(integrator->step);
    MYFree(integrator->sys);
    MYFree(integrator);
    return;
}
