/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 2021, Yueyang
**/

#include "etCore.h"
#define LOCAL_DEBUG 0

typedef struct tagTeukParams
{
    REAL8 delta;
    REAL8 sigma2;
    REAL8 sigma;
    REAL8 bb;
    REAL8 cc1;
    REAL8 cc2;
    REAL8 cc3;
    REAL8 m31;
    REAL8 m32;
    REAL8 a31;
    REAL8 a32;
    REAL8 a33;
    REAL8 a34;
    REAL8 lder2;
    REAL8 lder1;
}TeukParams;

typedef struct tagDiscreteGridParams
{
    REAL8 rt_min;
    REAL8 rt_max;
    REAL8 rt_step;
    UINT rt_num;

    REAL8 ang_min;
    REAL8 ang_max;
    REAL8 ang_step;
    UINT ang_num;

    REAL8 t_min;
    REAL8 t_max;
    REAL8 t_step;
    UINT t_num;
}DiscreteGridParams;

static INT compute_Teukolsky_mode(INT modeM, CoreParams *params, DiscreteGridParams *gp);

static INT compute_TeukParams(TeukParams *params, INT m, INT s,
                    REAL8 r, REAL8 ang, REAL8 a, REAL8 M);

static REAL8 calculate_rt(REAL8 r, REAL8 a2);

static REAL8 calculate_delta(REAL8 r, REAL8 a2);

static REAL8 calculate_drt_dr(REAL8 r, REAL8 a2);

static REAL8 calculate_radius_from_tortoise(REAL8 rt, REAL8 a, REAL8 eps);

static REAL8 calculate_small_radius_from_tortoise(REAL8 rt, REAL8 a, REAL8 eps);

static INT SetInitConditions_BellPulse(REAL8Array **init, CoreParams *cp, DiscreteGridParams *gp, 
                        INT modeL, INT modeM, INT type, REAL8 mag, REAL8 width, REAL8 r0, REAL8 v0);

/**
 * @brief Numerically solving Teukolsky equation in TD
 * 
 * @param params 
 * @return INT 
 */
INT Solve_Teukolsky_TD(CoreParams *params)
{
    DiscreteGridParams gp;
    INT status = CEV_SUCCESS;

    // set grid params
    gp.rt_min = -100.; // [M]
    gp.rt_max = 500.; // [M]
    gp.ang_min = 0.0;
    gp.ang_max = CST_PI;
    gp.t_min = 0.0;
    gp.t_max = 50.;

    gp.rt_num = 8001;
    gp.ang_num = 33;
    gp.rt_step = (gp.rt_max - gp.rt_min) / gp.rt_num;
    gp.ang_step = (gp.ang_max - gp.ang_min) / gp.ang_num;

    gp.t_step = params->step_ratio*GET_MIN(gp.rt_step, 5.*gp.ang_step);
    gp.t_num = (INT)((gp.t_max - gp.t_min) / gp.t_step);
    INT mm;
    mm = 2;
    PRINT_LOG_INFO(LOG_INFO, "Calculate mode m = %d", mm);
    status = compute_Teukolsky_mode(mm, params, &gp);
    // for (mm=-2; mm < 3; mm++)
    // {
    //     PRINT_LOG_INFO(LOG_INFO, "Calculate mode m = %d", mm);
    //     status = compute_Teukolsky_mode(mm, params, &gp);
    //     if (status != CEV_SUCCESS)
    //         goto QUIT;
    // }
QUIT:
    return status;
}

/**
 * @brief Calculate Teukolsky auxiliary parameters
 * 
 * @param params 
 * @param m 
 * @param s 
 * @param r 
 * @param ang 
 * @param a 
 * @param M 
 * @return INT 
 */
static INT compute_TeukParams(TeukParams *params, INT m, INT s,
                    REAL8 r, REAL8 ang, REAL8 a, REAL8 M)
{
    memset(params, 0, sizeof(TeukParams));
    REAL8 r2 = r*r;
    REAL8 a2 = a*a;
    REAL8 s_ang = GET_SIN(ang);
    REAL8 c_ang = GET_COS(ang);
    REAL8 s_ang2 = s_ang*s_ang;
    REAL8 c_ang2 = c_ang*c_ang;
    REAL8 prt_b;
    INT m2 = m*m;
    INT s2 = s*s;
    params->delta = r2 - 2.*M*r + a2;
    params->sigma2 = pow(r2+a2, 2.) - a2*params->delta*s_ang2;
    params->sigma = GET_SQRT(params->sigma2);
    prt_b = (params->delta/params->sigma) * ( 2.*r/(r2+a2) - (2.*r*(r2+a2) -a2*s_ang2*(r-M) )/params->sigma2 ); 
    params->bb = (r2+a2) / params->sigma2;
    params->cc1 = 2.*s*(-3.*M*r2+M*a2+r*(r2+a2)) / params->sigma2;
    params->cc2 = -2.*( r*params->delta*(1.+s)-(a2-r2)*M*s )/params->sigma2 - 
        6.*params->delta*params->bb/r/params->sigma;
    params->cc3 = 2.*a*(2.*r*M*m + params->delta*s*c_ang)/params->sigma2;
    params->a31 = params->delta*(m2 + 2.*c_ang*s*m + c_ang2*s2-s_ang2*s) / params->sigma2/s_ang2 - 
        6.*params->delta*(a2+r*(r*(s+2) - M*(s+3))) / r2 / params->sigma2;
    params->a32 = (4.*M*(r-1)*s*m*a*M + 6.*a*m*params->delta/r) / params->sigma2;
    params->a33 = params->cc1;
    params->a34 = -params->cc3;
    params->m31 = -params->bb*params->cc1+params->bb*prt_b + params->cc2;
    params->m32 = params->bb*params->cc3 + 2.*a*m*(r2+a2)/params->sigma2;
    params->lder2 = -params->delta / params->sigma2;
    params->lder1 = -(c_ang/s_ang) * params->delta / params->sigma2;
    return CEV_SUCCESS;
}

/**
 * @brief calculate single mode
 * 
 * @param modeM 
 * @param params 
 * @param gp 
 * @return INT 
 */
static INT compute_Teukolsky_mode(INT modeM, CoreParams *params, DiscreteGridParams *gp)
{
    register UINT it, jr, kf;
    CHAR fsave[STR_COMM_SIZE];
    TeukParams tp_rj;
    TeukParams tp_rjph;

    REAL8Array *rHalfCoord = NULL;
    REAL8Array *init = NULL;
    REAL8Array *buffer = NULL;
    REAL8Array *out = NULL;

    REAL8Array *rCoord = NULL;
    REAL8Array *rtCoord = NULL;
    REAL8Array *angCoord = NULL;
    REAL8Array *tCoord = NULL;

    UINT rt_half_num = 2*gp->rt_num - 1;
    // allocate memory
    PRINT_LOG_INFO(LOG_INFO, "[Mode %d]Allocate Memory", modeM);
    buffer = CreateREAL8Array(4, 4, 3, rt_half_num, gp->ang_num);
    out = CreateREAL8Array(4, 2, gp->t_num, gp->rt_num, gp->ang_num);
    memset(buffer->data, 0, buffer->size * sizeof(REAL8));
    memset(out->data, 0, out->size * sizeof(REAL8));
    UINT lu_rt = gp->ang_num;
    UINT lu_t = lu_rt * rt_half_num;
    UINT lu_u = lu_t * 3;
    REAL8 *uR, *uI, *vR, *vI, *UR, *UI;
    uR = buffer->data;
    uI = uR + lu_u;
    vR = uI + lu_u;
    vI = vR + lu_u;
    UINT LU_rt = gp->ang_num;
    UINT LU_t = LU_rt * gp->rt_num;
    UINT LU_u = LU_t * gp->t_num;
    UR = out->data;
    UI = UR + LU_u; 
    /**
     * @brief Set initial conditions
     * 
     */
    SetInitConditions_BellPulse(&init, params, gp, 2, modeM, 0, 0.1, 50, 75, 0.01);
    memcpy(uR, init->data, lu_t * sizeof(REAL8));
    memcpy(vR, init->data + lu_t, lu_t * sizeof(REAL8));
    // copy uR, uI to out array
    for (jr=0; jr < gp->rt_num; jr++)
    {
        for (kf=0; kf < gp->ang_num; kf++)
        {
            UR[kf + jr*LU_rt] = uR[kf + 2*jr*lu_rt];
            UI[kf + jr*LU_rt] = uI[kf + 2*jr*lu_rt];
        }
    }
    // sprintf(fsave, "%s/mode_%d.h5", params->prefix, modeM);
    // DumpREAL8ArrayTohdf5(fsave, "data", out, 1);
    // goto QUIT;
    /**
     * @brief Set r coordinates
     * 
     */
    PRINT_LOG_INFO(LOG_INFO, "[Mode %d]Set radius coordinates", modeM);
    rHalfCoord = CreateREAL8Array(1, rt_half_num);
    REAL8 r_j, rt_j;
    REAL8 r_jph, rt_jph;
    for (jr = 0; jr<rt_half_num; jr++ )
    {
        rt_j = gp->rt_min + jr * gp->rt_step/2.;
        rHalfCoord->data[jr] = calculate_small_radius_from_tortoise(rt_j, params->spin, 1e-3);
    }

    REAL8 ang_k;
    REAL8 Source, ForceT = 0.0, LderTerm;
    REAL8 duRdrt, duIdrt;
    PRINT_LOG_INFO(LOG_INFO, "[Mode %d]Start iteration", modeM);
    for (it = 0; it < gp->t_num-1; it++)
    {
        PRINT_LOG_INFO(LOG_INFO, "PROC:%d/%d", it, gp->t_num-1);
        for (jr = 0; jr < rt_half_num-2; jr++)
        {
            // radius
            rt_jph = gp->rt_min + (1.+jr)*gp->rt_step/2.;
            r_jph = rHalfCoord->data[jr+1];
            // First step
            for (kf = 1; kf < gp->ang_num-1; kf++)
            {
                // theta_k
                ang_k = gp->ang_min + kf * gp->ang_step;

                // compute aux parameters
                compute_TeukParams(&tp_rjph, modeM, -2, r_jph, ang_k, params->spin, 1.);
                uR[kf+(jr+1)*lu_rt+lu_t] = (uR[kf+(jr+2)*lu_rt]+uR[kf+(jr)*lu_rt])/2. - 
                    (gp->t_step/2.) * ( tp_rjph.bb*(uR[kf+(jr+2)*lu_rt] - uR[kf+(jr)*lu_rt])/gp->rt_step -
                    vR[kf+(jr+1)*lu_rt]);
                uI[kf+(jr+1)*lu_rt+lu_t] = (uI[kf+(jr+2)*lu_rt]+uI[kf+(jr)*lu_rt])/2. - 
                    (gp->t_step/2.) * ( tp_rjph.bb*(uI[kf+(jr+2)*lu_rt] - uI[kf+(jr)*lu_rt])/gp->rt_step -
                    vI[kf+(jr+1)*lu_rt]);
                // debug
#if LOCAL_DEBUG
                if (isnan(uR[kf+(jr+1)*lu_rt+lu_t]) || isnan(uI[kf+(jr+1)*lu_rt+lu_t]))
                {
                    print_debug("time[%d+1/2], ang[%d+1], r[%d]\n", it, jr, kf);
                    print_debug("uR = %e\n", uR[kf+(jr+1)*lu_rt+lu_t]);
                    print_debug("uI = %e\n", uI[kf+(jr+1)*lu_rt+lu_t]);
                    print_debug("vR = %e\n", vR[kf+(jr+1)*lu_rt+lu_t]);
                    print_debug("vI = %e\n", vI[kf+(jr+1)*lu_rt+lu_t]);
                    print_debug("b = %e\n", tp_rj.bb);
                    print_debug("last state:(t[%d+0], ang[%d], r[%d]):\n", it, jr, kf); 
                    print_debug("\t%e, %e, %e, %e\n", uR[kf+(jr)*lu_rt], uI[kf+(jr)*lu_rt], vR[kf+(jr)*lu_rt], vI[kf+(jr)*lu_rt]);
                    print_debug("last state:(t[%d+0], ang[%d], r[%d+1]):\n", it, jr, kf); 
                    print_debug("\t%e, %e, %e, %e\n", uR[kf+(jr+1)*lu_rt], uI[kf+(jr+1)*lu_rt], vR[kf+(jr+1)*lu_rt], vI[kf+(jr+1)*lu_rt]);
                    print_debug("last state:(t[%d+0], ang[%d], r[%d+2]):\n", it, jr, kf); 
                    print_debug("\t%e, %e, %e, %e\n", uR[kf+(jr+2)*lu_rt], uI[kf+(jr+2)*lu_rt], vR[kf+(jr+2)*lu_rt], vI[kf+(jr+2)*lu_rt]);
                    goto QUIT;
                }
#endif
                // l31(uR)
                LderTerm = tp_rjph.lder2 * (uR[kf+1+(jr+1)*lu_rt] + uR[kf-1+(jr+1)*lu_rt] - 2.*uR[kf+(jr+1)*lu_rt])/gp->ang_step/gp->ang_step + 
                    tp_rjph.lder1 * (uR[kf+1+(jr+1)*lu_rt] - uR[kf-1+(jr+1)*lu_rt])/2./gp->ang_step;
                // duR/drt
                duRdrt = (uR[kf+(jr+2)*lu_rt]-uR[kf+(jr)*lu_rt])/gp->rt_step;
                // duI/drt
                duIdrt = (uI[kf+(jr+2)*lu_rt]-uI[kf+(jr)*lu_rt])/gp->rt_step;
                Source = ForceT - tp_rjph.m31 * duRdrt - tp_rjph.m32*duIdrt - LderTerm - 
                    (tp_rjph.a31*uR[kf+(jr+1)*lu_rt] +tp_rjph.a32*uI[kf+(jr+1)*lu_rt] + tp_rjph.a33*vR[kf+(jr+1)*lu_rt]+tp_rjph.a34*vI[kf+(jr+1)*lu_rt]);
                vR[kf+(jr+1)*lu_rt+lu_t] = (vR[kf+(jr+2)*lu_rt]+vR[kf+(jr)*lu_rt])/2. - 
                    (gp->t_step/2.) * ( (-tp_rjph.bb)*(vR[kf+(jr+2)*lu_rt] - vR[kf+(jr)*lu_rt])/gp->rt_step -
                    Source);
                
                // l31*uI
                LderTerm = tp_rjph.lder2 * (uI[kf+1+(jr+1)*lu_rt] + uI[kf-1+(jr+1)*lu_rt] - 2.*uI[kf+(jr+1)*lu_rt])/gp->ang_step/gp->ang_step + 
                    tp_rjph.lder1 * (uI[kf+1+(jr+1)*lu_rt] - uI[kf-1+(jr+1)*lu_rt])/2./gp->ang_step;
                Source = ForceT + tp_rjph.m32 * duRdrt - tp_rjph.m31*duIdrt - LderTerm - 
                    (-tp_rjph.a32*uR[kf+(jr+1)*lu_rt] +tp_rjph.a31*uI[kf+(jr+1)*lu_rt] - tp_rjph.a34*vR[kf+(jr+1)*lu_rt]+tp_rjph.a33*vI[kf+(jr+1)*lu_rt]);
                vI[kf+(jr+1)*lu_rt+lu_t] = (vI[kf+(jr+2)*lu_rt]+vI[kf+(jr)*lu_rt])/2. - 
                    (gp->t_step/2.) * ( (-tp_rjph.bb)*(vI[kf+(jr+2)*lu_rt] - vI[kf+(jr)*lu_rt])/gp->rt_step -
                    Source);
#if LOCAL_DEBUG
                // debug
                if (isnan(vR[kf+(jr+1)*lu_rt+lu_t]) || isnan(vI[kf+(jr+1)*lu_rt+lu_t]))
                {
                    print_debug("time[%d+1/2], ang[%d+1], r[%d]\n", it, jr, kf);
                    print_debug("uR = %e\n", uR[kf+(jr+1)*lu_rt+lu_t]);
                    print_debug("uI = %e\n", uI[kf+(jr+1)*lu_rt+lu_t]);
                    print_debug("vR = %e\n", vR[kf+(jr+1)*lu_rt+lu_t]);
                    print_debug("vI = %e\n", vI[kf+(jr+1)*lu_rt+lu_t]);
                    print_debug("b = %e\n", tp_rj.bb);
                    print_debug("last state:(t[%d+0], ang[%d], r[%d]):\n", it, jr, kf); 
                    print_debug("\t%e, %e, %e, %e\n", uR[kf+(jr)*lu_rt], uI[kf+(jr)*lu_rt], vR[kf+(jr)*lu_rt], vI[kf+(jr)*lu_rt]);
                    print_debug("last state:(t[%d+0], ang[%d], r[%d+1]):\n", it, jr, kf); 
                    print_debug("\t%e, %e, %e, %e\n", uR[kf+(jr+1)*lu_rt], uI[kf+(jr+1)*lu_rt], vR[kf+(jr+1)*lu_rt], vI[kf+(jr+1)*lu_rt]);
                    print_debug("last state:(t[%d+0], ang[%d], r[%d+2]):\n", it, jr, kf); 
                    print_debug("\t%e, %e, %e, %e\n", uR[kf+(jr+2)*lu_rt], uI[kf+(jr+2)*lu_rt], vR[kf+(jr+2)*lu_rt], vI[kf+(jr+2)*lu_rt]);
                    print_debug("duRdrt = %e, duIdrt = %e\n", duRdrt, duIdrt);
                    print_debug("LderTerm = %e\n", LderTerm);
                    print_debug("Source = %e\n", Source);
                    print_debug("m31 = %e, m32 = %e\n", tp_rjph.m31, tp_rjph.m32);
                    print_debug("a31 = %e, a32 = %e, a33 = %e, a34 = %3\n", tp_rjph.a31, tp_rjph.a32, tp_rjph.a33, tp_rjph.a34);
                    print_debug("c1 = %e, c2 = %e, c3 = %e\n", tp_rjph.cc1, tp_rjph.cc2, tp_rjph.cc3);
                    goto QUIT;
                }
#endif

            }
            // fix angular boundary conditions
            // we dont need to care about v
            if (abs(modeM)%2)
            {
                // even mode, dudf = 0
                uR[(jr+1)*lu_rt+lu_t] = uR[1+(jr+1)*lu_rt+lu_t];
                uI[(jr+1)*lu_rt+lu_t] = uI[1+(jr+1)*lu_rt+lu_t];
                uR[gp->ang_num-1+(jr+1)*lu_rt+lu_t] = uR[gp->ang_num-2+(jr+1)*lu_rt+lu_t];
                uI[gp->ang_num-1+(jr+1)*lu_rt+lu_t] = uI[gp->ang_num-2+(jr+1)*lu_rt+lu_t];
            } else {
                // odd mode, u = 0
                uR[(jr+1)*lu_rt+lu_t] = 0.0;
                uI[(jr+1)*lu_rt+lu_t] = 0.0;
                uR[gp->ang_num-1+(jr+1)*lu_rt+lu_t] = 0.0;
                uI[gp->ang_num-1+(jr+1)*lu_rt+lu_t] = 0.0;
            }
        }

        // fix radius boundary conditions
        for (kf = 0; kf < gp->ang_num; kf++)
        {
            uR[kf + lu_t] = 0.0;
            uI[kf + (rt_half_num-1)*lu_rt + lu_t] = 0.0;
            vR[kf + lu_t] = 0.0;
            vI[kf + (rt_half_num-1)*lu_rt + lu_t] = 0.0;
        }

        // Second iteration
        for (jr = 1; jr < rt_half_num-1; jr++)
        {
            rt_j = gp->rt_min + jr * gp->rt_step/2.;
            r_j = rHalfCoord->data[jr];
            for (kf = 1; kf < gp->ang_num-1; kf++)
            {
                ang_k = gp->ang_min + kf * gp->ang_step;
                // compute aux parameters
                compute_TeukParams(&tp_rj, modeM, -2, r_j, ang_k, params->spin, 1.);
                uR[kf+(jr)*lu_rt+2*lu_t] = uR[kf+jr*lu_rt+lu_t] - 
                    (gp->t_step) * ( tp_rj.bb*(uR[kf+(jr+1)*lu_rt+lu_t] - uR[kf+(jr-1)*lu_rt+lu_t])/gp->rt_step -
                    vR[kf+(jr)*lu_rt+lu_t]);
                uI[kf+(jr)*lu_rt+2*lu_t] = uI[kf+jr*lu_rt+lu_t] - 
                    (gp->t_step) * ( tp_rj.bb*(uI[kf+(jr+1)*lu_rt+lu_t] - uI[kf+(jr-1)*lu_rt+lu_t])/gp->rt_step -
                    vI[kf+(jr)*lu_rt+lu_t]);
#if LOCAL_DEBUG

                // debug
                if (isnan(uR[kf+(jr)*lu_rt+2*lu_t]) || isnan(uI[kf+(jr)*lu_rt+2*lu_t]))
                {
                    print_debug("time[%d+1], ang[%d], r[%d]\n", it, jr, kf);
                    print_debug("uR = %e\n", uR[kf+(jr)*lu_rt+2*lu_t]);
                    print_debug("uI = %e\n", uI[kf+(jr)*lu_rt+2*lu_t]);
                    print_debug("vR = %e\n", vR[kf+(jr)*lu_rt+2*lu_t]);
                    print_debug("vI = %e\n", vI[kf+(jr)*lu_rt+2*lu_t]);
                    print_debug("b = %e\n", tp_rj.bb);
                    print_debug("last state:(t[%d+1/2], ang[%d], r[%d]):\n", it, jr, kf); 
                    print_debug("\t%e, %e, %e, %e\n", uR[kf+(jr)*lu_rt+lu_t], uI[kf+(jr)*lu_rt+lu_t], vR[kf+(jr)*lu_rt+lu_t], vI[kf+(jr)*lu_rt+lu_t]);
                    print_debug("last state:(t[%d+1/2], ang[%d], r[%d-1]):\n", it, jr, kf); 
                    print_debug("\t%e, %e, %e, %e\n", uR[kf+(jr-1)*lu_rt+lu_t], uI[kf+(jr-1)*lu_rt+lu_t], vR[kf+(jr-1)*lu_rt+lu_t], vI[kf+(jr-1)*lu_rt+lu_t]);
                    print_debug("last state:(t[%d+1/2], ang[%d], r[%d+1]):\n", it, jr, kf); 
                    print_debug("\t%e, %e, %e, %e\n", uR[kf+(jr+1)*lu_rt+lu_t], uI[kf+(jr+1)*lu_rt+lu_t], vR[kf+(jr+1)*lu_rt+lu_t], vI[kf+(jr+1)*lu_rt+lu_t]);
                    goto QUIT;
                }
#endif
                // l31(uR)
                LderTerm = tp_rj.lder2 * (uR[kf+1+(jr)*lu_rt+lu_t] + uR[kf-1+(jr)*lu_rt] - 2.*uR[kf+(jr)*lu_rt+lu_t])/gp->ang_step/gp->ang_step + 
                    tp_rj.lder1 * (uR[kf+1+(jr)*lu_rt+lu_t] - uR[kf-1+(jr)*lu_rt+lu_t])/2./gp->ang_step;
                // duR/drt
                duRdrt = (uR[kf+(jr+1)*lu_rt+lu_t]-uR[kf+(jr-1)*lu_rt+lu_t])/gp->rt_step;
                // duI/drt
                duIdrt = (uI[kf+(jr+1)*lu_rt+lu_t]-uI[kf+(jr-1)*lu_rt+lu_t])/gp->rt_step;

                Source = ForceT - tp_rj.m31 * duRdrt - tp_rj.m32*duIdrt - LderTerm - 
                    (tp_rj.a31*uR[kf+(jr)*lu_rt+lu_t] +tp_rj.a32*uI[kf+(jr)*lu_rt+lu_t] + tp_rj.a33*vR[kf+(jr)*lu_rt+lu_t]+tp_rj.a34*vI[kf+(jr)*lu_rt+lu_t]);
                vR[kf+(jr)*lu_rt+2*lu_t] = vR[kf+(jr)*lu_rt+lu_t] - 
                    (gp->t_step) * ( (-tp_rj.bb)*(vR[kf+(jr+1)*lu_rt+lu_t] - vR[kf+(jr-1)*lu_rt+lu_t])/gp->rt_step -
                    Source);

                // l31*uI
                LderTerm = tp_rj.lder2 * (uI[kf+1+(jr)*lu_rt+lu_t] + uI[kf-1+(jr)*lu_rt+lu_t] - 2.*uI[kf+(jr)*lu_rt+lu_t])/gp->ang_step/gp->ang_step + 
                    tp_rj.lder1 * (uI[kf+1+(jr)*lu_rt+lu_t] - uI[kf-1+(jr)*lu_rt+lu_t])/2./gp->ang_step;
                Source = ForceT + tp_rj.m32 * duRdrt - tp_rj.m31*duIdrt - LderTerm - 
                    (-tp_rj.a32*uR[kf+(jr)*lu_rt+lu_t] +tp_rj.a31*uI[kf+(jr)*lu_rt+lu_t] - tp_rj.a34*vR[kf+(jr)*lu_rt+lu_t]+tp_rj.a33*vI[kf+(jr)*lu_rt+lu_t]);
                vI[kf+(jr)*lu_rt+2*lu_t] = vI[kf+(jr)*lu_rt+lu_t] - 
                    (gp->t_step) * ( (-tp_rj.bb)*(vI[kf+(jr+1)*lu_rt+lu_t] - vI[kf+(jr-1)*lu_rt+lu_t])/gp->rt_step -
                    Source);
            }

            // fix angular boundary conditions
            // we dont need to care about v
            if (abs(modeM)%2)
            {
                // even mode, dudf = 0
                uR[(jr)*lu_rt+2*lu_t] = uR[1+(jr)*lu_rt+2*lu_t];
                uI[(jr)*lu_rt+2*lu_t] = uI[1+(jr)*lu_rt+2*lu_t];
                uR[gp->ang_num-1+(jr)*lu_rt+2*lu_t] = uR[gp->ang_num-2+(jr)*lu_rt+2*lu_t];
                uI[gp->ang_num-1+(jr)*lu_rt+2*lu_t] = uI[gp->ang_num-2+(jr)*lu_rt+2*lu_t];
            } else {
                // odd mode, u = 0
                uR[(jr)*lu_rt+2*lu_t] = 0.0;
                uI[(jr)*lu_rt+2*lu_t] = 0.0;
                uR[gp->ang_num-1+(jr)*lu_rt+2*lu_t] = 0.0;
                uI[gp->ang_num-1+(jr)*lu_rt+2*lu_t] = 0.0;
            }
        }
        // fix radius boundary conditions
        for (kf = 0; kf < gp->ang_num; kf++)
        {
            uR[kf + 2*lu_t] = 0.0;
            uI[kf + (rt_half_num-1)*lu_rt + 2*lu_t] = 0.0;
            vR[kf + 2*lu_t] = 0.0;
            vI[kf + (rt_half_num-1)*lu_rt + 2*lu_t] = 0.0;
        }

        // copy uR, uI to out array
        for (jr=0; jr < gp->rt_num; jr++)
        {
            for (kf=0; kf < gp->ang_num; kf++)
            {
                UR[kf + jr*LU_rt + (it+1)*LU_t] = uR[kf + 2*jr*lu_rt + 2*lu_t];
                UI[kf + jr*LU_rt + (it+1)*LU_t] = uI[kf + 2*jr*lu_rt + 2*lu_t];
            }
        }
        memcpy(uR, uR + 2*lu_t, lu_t * sizeof(REAL8));
        memcpy(uI, uI + 2*lu_t, lu_t * sizeof(REAL8));
        memcpy(vR, vR + 2*lu_t, lu_t * sizeof(REAL8));
        memcpy(vI, vI + 2*lu_t, lu_t * sizeof(REAL8));
    }

    tCoord = CreateREAL8Array(1, gp->t_num);
    rCoord = CreateREAL8Array(1, gp->rt_num);
    rtCoord = CreateREAL8Array(1, gp->rt_num);
    angCoord = CreateREAL8Array(1, gp->ang_num);
    for (it=0;it<gp->t_num;it++)
        tCoord->data[it] = gp->t_min + it * gp->t_step;
    for (it=0;it<gp->rt_num;it++)
    {
        rtCoord->data[it] = gp->rt_min + it * gp->rt_step;
        rCoord->data[it] = rHalfCoord->data[2*it];
    }
    for (it=0;it<gp->ang_num;it++)
        angCoord->data[it] = gp->ang_min + it * gp->ang_step;

    sprintf(fsave, "%s/mode_%d.h5", params->prefix, modeM);
    DumpREAL8ArrayTohdf5(fsave, "data", out, 1);
    DumpREAL8ArrayTohdf5(fsave, "time", tCoord, 0);
    DumpREAL8ArrayTohdf5(fsave, "radius", rCoord, 0);
    DumpREAL8ArrayTohdf5(fsave, "tortoise_radius", rtCoord, 0);
    DumpREAL8ArrayTohdf5(fsave, "azimuth", angCoord, 0);
QUIT:
    STRUCTFREE(out, REAL8Array);
    STRUCTFREE(buffer, REAL8Array);
    STRUCTFREE(rHalfCoord, REAL8Array);
    STRUCTFREE(rCoord, REAL8Array);
    STRUCTFREE(rtCoord, REAL8Array);
    STRUCTFREE(tCoord, REAL8Array);
    STRUCTFREE(angCoord, REAL8Array);
    STRUCTFREE(init, REAL8Array);
    return CEV_SUCCESS;
}

static REAL8 calculate_rt(REAL8 r, REAL8 a2)
{
    REAL8 rp = 1.+GET_SQRT(1-a2);
    REAL8 rm = 1.-GET_SQRT(1-a2);
    REAL8 rpm = rp-rm;
    return r + 2.*rp*log((r-rp)/2.)/rpm - 2.*rm*log((r-rm)/2.)/rpm;
}

static REAL8 calculate_delta(REAL8 r, REAL8 a2)
{
    return r*r - 2.*r + a2;
}

static REAL8 calculate_drt_dr(REAL8 r, REAL8 a2)
{
    return (r*r + a2) / calculate_delta(r, a2);
}

static REAL8 calculate_radius_from_tortoise(REAL8 rt, REAL8 a, REAL8 eps)
{
    REAL8 a2 = a*a;
    REAL8 rp = 1.+GET_SQRT(1-a2);
    REAL8 rm = 1.-GET_SQRT(1-a2);
    REAL8 rpm = rp - rm;
#define FUNC_CALC_RT(RR, RP, RM, RPM) \
    (RR + 2.*RP*log((RR-RP)/2.)/RPM - 2.*RM*log((RR-RM)/2.)/RPM)
#define FUNC_CALC_DRTDR(RR, A2) \
    ( (RR*RR+A2) / (RR*RR - 2.*RR + A2) )
    REAL8 r_0 = GET_MAX(rt, rp+0.1);
    REAL8 rt_0 = FUNC_CALC_RT(r_0, rp, rm, rpm);
    REAL8 r_1 = r_0 - (rt_0-rt) / FUNC_CALC_DRTDR(r_0, a2);
    REAL8 rt_1 = FUNC_CALC_RT(r_1, rp, rm, rpm);
    while (fabs(rt_1 - rt) > eps )
    {
        r_0 = r_1;
        rt_0 = FUNC_CALC_RT(r_0, rp, rm, rpm);
        r_1 = r_0 - (rt_0-rt) / FUNC_CALC_DRTDR(r_0, a2);
        rt_1 = FUNC_CALC_RT(r_1, rp, rm, rpm);
    }
#undef FUNC_CALC_RT
#undef FUNC_CALC_DRTDR
    return r_1;
}

/**
 * @brief Calculate radius from tortoise r by Newtonian iteration
 * 
 */
static REAL8 calculate_small_radius_from_tortoise(REAL8 rt, REAL8 a, REAL8 eps)
{
    REAL8 a2 = a*a;
    REAL8 rp = 1.+GET_SQRT(1-a2);
    REAL8 rm = 1.-GET_SQRT(1-a2);
    REAL8 rpm = rp - rm;
    REAL8 fp = 2.*rp / rpm;
    REAL8 fm = 2.*rm / rpm;
#define FUNC_CALC_RT_SMALL(PP, RX, RP, RM, FP, FM, RPM) \
    (FP*(-PP-CST_LN2) + RP + RX - FM*log((RX+RPM)/2.))
#define FUNC_CALC_DRDRT_SMALL(A2, RX, RX2, RP, RPM) \
    ( (RX2 + RX*RPM) / (RX2 + 2.*RX*RP + RP*RP + A2) )
    REAL8 p0 = rt > 0.0 ? -log(rt) : 1.0;
    REAL8 r0 = exp(-p0);
    REAL8 r02 = exp(-2.*p0);
    REAL8 rt0 = FUNC_CALC_RT_SMALL(p0, r0, rp, rm, fp, fm, rpm);
    REAL8 p1 = p0 - (rt0 - rt) * FUNC_CALC_DRDRT_SMALL(a2, r0, r02, rp, rpm) / (-r0);
    REAL8 r1 = exp(-p1);
    REAL8 r12 = exp(-2.*p1);
    REAL8 rt1 = FUNC_CALC_RT_SMALL(p1, r1, rp, rm, fp, fm, rpm);
    while (fabs(rt1 - rt) > eps)
    {
        p0 = p1;
        r0 = exp(-p0);
        r02 = exp(-2.*p0);
        rt0 = FUNC_CALC_RT_SMALL(p0, r0, rp, rm, fp, fm, rpm);
        p1 = p0 - (rt0 - rt) * FUNC_CALC_DRDRT_SMALL(a2, r0, r02, rp, rpm) / (-r0);
        r1 = exp(-p1);
        r12 = exp(-2.*p1);
        rt1 = FUNC_CALC_RT_SMALL(p1, r1, rp, rm, fp, fm, rpm);
    }
#undef FUNC_CALC_RT_SMALL
#undef FUNC_CALC_DRDRT_SMALL
    return r1 + rp;
}


/**
 * @brief Set the Init Conditions object
 * 
 * @param init 
 * @param gp 
 * @param modeM 
 * @param type 
 * @return INT 
 */
static INT SetInitConditions_BellPulse(REAL8Array **init, CoreParams *cp, DiscreteGridParams *gp, 
                        INT modeL, INT modeM, INT type, REAL8 mag, REAL8 width, REAL8 r0, REAL8 v0)
{
    INT i, j;
    UINT rt_num = 2*gp->rt_num - 1;
    UINT ang_num = gp->ang_num;
    UINT lu_rt = ang_num;
    UINT lu_u = rt_num * lu_rt;
    REAL8Array *out = NULL;
    out = CreateREAL8Array(3, 2, rt_num, ang_num);
    REAL8 *u = out->data;
    REAL8 *v = u + lu_u;
    REAL8 r, rt, ang;
    REAL8 sig2, r2, s_ang2, bb, delta, a2 = cp->spin*cp->spin;
    REAL8 Y;
    COMPLEX16 swY;
    for (i=0; i<rt_num; i++)
    {
        rt = gp->rt_min + 0.5*i*gp->rt_step;
        r = calculate_small_radius_from_tortoise(rt, cp->spin, 1e-3);
        r2 = r*r;
        delta = r2 - 2.*r + a2;
        for(j=0; j<ang_num; j++)
        {
            ang = gp->ang_min + j*gp->ang_step;
            s_ang2 = pow(GET_SIN(ang), 2.);
            SpinWeightedSphericalHarmonic(ang, 0.0, -2, modeL, modeM, &swY);
            Y = C_REAL(swY);
            sig2 = pow(r2+a2, 2.) - a2*delta*s_ang2;
            bb = (r2+a2) / sig2;
            u[j + i*lu_rt] = mag*exp(-pow((rt-r0)/width,2.)) * Y;
            v[j + i*lu_rt] = mag*2*u[j + i*lu_rt]*(r0-rt)*(bb-v0)*Y/width/width;
        }
    }
    *init = out;
    return CEV_SUCCESS;
}