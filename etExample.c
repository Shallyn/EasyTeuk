/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include "etExample.h"
#define LOCAL_USE_DEBUG 0

/**
 *
 * Example: 1-d string vibration excitation
 * 
 */

static REAL8 Example_compute_string_excitation(INT it, INT ir, REAL8 tStep, REAL8 rStep)
{
    // return 0.0;
    REAL8 t = it * tStep;
    REAL8 r = ir * rStep;
    REAL8 mag = 0.1 / (1 + t);
    REAL8 rt = 10.0;
    // REAL8 retP1, retP2;
    REAL8 delta_r, delta_t;
    delta_r = numerical_discrete_delta(0.0, ir, rt, rStep);
    // delta_t = numerical_discrete_delta(0.0, it, 10.0, tStep);
    // if (retP1 > 0.0)
    //     print_debug("(t,r,phi,valr,valphi) = (%f,%f,%f,%f,%f)\n", t, r, phi,rt,phit);
    return mag * delta_r;
}


INT Example_String_Vibration(CoreParams *params)
{
    INT i,j,k;
    INT failed = 0; /* used throughout */
    CHAR fsave[STR_COMM_SIZE];
    REAL8Array *tCoord = NULL;
    REAL8Array *rCoord = NULL;
    REAL8 tStep = 0.08;
    REAL8 rStep = 0.1;

    REAL8 tHalfStep = tStep / 2.;
    REAL8 rHalfStep = rStep / 2.;

    INT numr = 250;
    INT numhR = 2*numr-1;
    INT numt = 1000;
    INT numhT = 2*numt-1;

    /**
     * Here we use Lax-Wendroff iteration, which consists two steps
     * and we need half step.
     */
    REAL8Array *buffer = NULL;
    REAL8Array *arr = NULL;

    INT idx_t = numhR, idx_t_itp = numr;
    INT idx_u = idx_t * numhT, idx_u_itp = idx_t_itp * numt;
    buffer = CreateREAL8Array(3, 2, numhT, numhR);
    memset(buffer->data, 0, sizeof(REAL8) * buffer->size);
    REAL8 *u = buffer->data;
    REAL8 *v = buffer->data + idx_u;
    REAL8 cc = 0.8;
    REAL8 force = 0.0;
    
    for (i=0; i<numt-1; i++)
    {
        for (j=0; j<numhR-1; j++)
        {
            // t_i = 2i*hdt, r_j = j*hdr, f_k = k*df
            // t_i+1/2 = (2i+1)*hdt, r_j+1/2 = (j+1)*hdr
            u[j+1 + (2*i+1)*idx_t] = (u[j+2+(2*i)*idx_t]+u[j+(2*i)*idx_t])/2. - 
                tHalfStep*( (u[j+2+(2*i)*idx_t]-u[j+(2*i)*idx_t])/rStep - v[j+1+(2*i)*idx_t] );
            force = Example_compute_string_excitation(2*i, j+1, tHalfStep, rHalfStep);
            v[j+1 + (2*i+1)*idx_t] = (v[j+2+(2*i)*idx_t]+v[j+(2*i)*idx_t])/2. - 
                tHalfStep*( (-cc)*(v[j+2+(2*i)*idx_t]-v[j+(2*i)*idx_t])/rStep - force );
        }

        // fix boundary conditions
        u[(2*i+1)*idx_t] = 0.0;
        u[numhR-1 + (2*i+1)*idx_t] = 0.0;
        v[(2*i+1)*idx_t] = 0.0;
        v[numhR-1 + (2*i+1)*idx_t] = 0.0;

        for (j=1; j<numhR-1; j++)
        {
            u[j+(2*i+2)*idx_t] = u[j+(2*i)*idx_t] - 
                tStep * ( (u[j+1+(2*i+1)*idx_t] - u[j-1+(2*i+1)*idx_t])/rStep - v[j+(2*i+1)*idx_t]);
            force = Example_compute_string_excitation(2*i+1, j, tHalfStep, rHalfStep);
            v[j+(2*i+2)*idx_t] = v[j+(2*i)*idx_t] - 
                tStep * ( (-cc)*(v[j+1+(2*i+1)*idx_t] - v[j-1+(2*i+1)*idx_t])/rStep - force);
        }

        // fix boundary conditions
        u[(2*i+2)*idx_t] = 0.0;
        u[numhR-1 + (2*i+2)*idx_t] = 0.0;
        v[(2*i+2)*idx_t] = 0.0;
        v[numhR-1 + (2*i+2)*idx_t] = 0.0;
    }

    /**
     *
     * Interpolation
     * 
     */
    arr = CreateREAL8Array(3, 2, numt, numr);
    for (i=1; i<numt; i++)
    {
        for (j=0; j<numr; j++)
        {
            arr->data[j + i*idx_t_itp] = buffer->data[2*j + 2*i*idx_t];
            arr->data[j + i*idx_t_itp + idx_u_itp] = buffer->data[2*j + 2*i*idx_t + idx_u];
        }
    }

    // dump data
    tCoord = CreateREAL8Array(1, numt);
    rCoord = CreateREAL8Array(1, numr);
    for (i=0;i<numt;i++)
        tCoord->data[i] = i * tStep;
    for (i=0;i<numr;i++)
        rCoord->data[i] = i * rStep;
    sprintf(fsave, "%s/%s", params->prefix, "string_osci.h5");
    DumpREAL8ArrayTohdf5(fsave, "data", arr, 1);
    DumpREAL8ArrayTohdf5(fsave, "tCoord", tCoord, 0);
    DumpREAL8ArrayTohdf5(fsave, "rCoord", rCoord, 0);

QUIT:
    STRUCTFREE(buffer, REAL8Array);
    STRUCTFREE(arr, REAL8Array);
    STRUCTFREE(tCoord, REAL8Array);
    STRUCTFREE(rCoord, REAL8Array);
    return CEV_SUCCESS;
}

/**
 * 
 * Example: 2-d plane vibration excitation
 * 
 */
static REAL8 Example_compute_source(INT it, INT ir, INT iphi, REAL8 tStep, REAL8 rStep, REAL8 phiStep)
{
    // return 0.0;
    REAL8 t = it * tStep;
    REAL8 r = ir * rStep;
    REAL8 phi = iphi * phiStep;
    REAL8 mag = 0.01;
    REAL8 rt = 10.0;
    REAL8 omega = 0.2;
    REAL8 phit = omega * t;
    // REAL8 retP1, retP2;
    while (phit > CST_2PI)
        phit -= CST_2PI;
    REAL8 delta_r, delta_phi, delta_t;
    delta_r = numerical_discrete_delta(0.0, ir, rt, rStep);
    delta_phi = numerical_discrete_delta(0.0, iphi, phit, phiStep);
    delta_t = numerical_discrete_delta(0.0, it, 10.0, tStep);
    // if (retP1 > 0.0)
    //     print_debug("(t,r,phi,valr,valphi) = (%f,%f,%f,%f,%f)\n", t, r, phi,rt,phit);
    return mag * delta_t * delta_r * delta_phi;
}


/**
 * @brief Example: 2-d plane vibration
 * 
 * @param params
 *   CoreParams *params
 * @return INT SUCCESS
 */
INT Example_Plane_Vibration(CoreParams *params)
{
    /**
     * 
     * Step 1. parameters and variables declaration
     * 
     */
    INT i,j,k;
    INT failed = 0; /* used throughout */
    CHAR fsave[STR_COMM_SIZE];

    /* set grid step, time step should not greater than space step */
    REAL8Array *tCoord = NULL;
    REAL8Array *rCoord = NULL;
    REAL8Array *fCoord = NULL;
    REAL8 tStep = 0.08;
    REAL8 rStep = 0.1;
    REAL8 fStepDegree = 3.; // degree
    REAL8 fStep = fStepDegree * CST_PI / 180.;

    REAL8 tHalfStep = tStep / 2.;
    REAL8 rHalfStep = rStep / 2.;
    REAL8 tHalfStep2 = tHalfStep*tHalfStep;
    REAL8 rHalfStep2 = rHalfStep*rHalfStep;
    INT numr = 250;
    INT numhR = 2*numr-1;
    INT numt = 1000;
    INT numhT = 2*numt-1;
    INT numf = (360/fStepDegree)+1;
    /**
     * Here we use Lax-Wendroff iteration, which consists two steps
     * and we need half step.
     */
    REAL8Array *buffer = NULL;
    REAL8Array *arr = NULL;

#if 0
    // DEBUG: Test delta function
    arr = CreateREAL8Array(3, numt, numr, numf);
    for (k=0; k<numt; k++)
    {
        for (i=0; i<numr; i++)
        {
            for (j=0; j<numf; j++)
            {
                // arr->data[j + i*numf] = i*rStep;
                arr->data[j + i*numf + k * numf*numr] = Example_compute_source(k, i, j, tStep, rStep, fStep);
                // rCoord->data[i] = i*rStep;
                // rCoord->data[i + numr] = numerical_discrete_delta(0.0, i, 10.0, rStep);
            }
        }
    }
    tCoord = CreateREAL8Array(1, numt);
    rCoord = CreateREAL8Array(1, numr);
    fCoord = CreateREAL8Array(1, numf);
    for (i=0;i<numt;i++)
        tCoord->data[i] = i * tStep;
    for (i=0;i<numr;i++)
        rCoord->data[i] = i * rStep;
    for (i=0;i<numf;i++)
        fCoord->data[i] = i * fStep;
    sprintf(fsave, "%s/%s", params->prefix, "test_delta.h5");
    DumpREAL8ArrayTohdf5(fsave, "data", arr, 1);
    DumpREAL8ArrayTohdf5(fsave, "rCoord", rCoord, 0);
    DumpREAL8ArrayTohdf5(fsave, "fCoord", fCoord, 0);
    DumpREAL8ArrayTohdf5(fsave, "tCoord", tCoord, 0);
    goto QUIT;
#endif
    /**
     * 
     * Step 2. Set initial conditions
     * 
     */

    buffer = CreateREAL8Array(4, 2, numhT, numhR, numf);
    memset(buffer->data, 0, sizeof(REAL8) * buffer->size);
    /**
     *
     * Iteration
     * 
     */
    REAL8 u_tn_ri, u_tn_rip;
    INT idx_tn_ri;
    INT idx_r = numf;
    INT idx_t = numf * numhR, idx_t_itp = numf * numr;
    INT idx_u = idx_t * numhT, idx_u_itp = idx_t_itp * numt;
    REAL8 riph, ri;
    REAL8 force, d2df2U;
    REAL8 *u, *v;
    u = buffer->data;
    v = buffer->data + idx_u;
    INT kp1, kp2, km1, km2;
    REAL8 cc = 0.8;
    REAL8 c2 = cc*cc;
    REAL8 u0_mean = 0.0, v0_mean = 0.0;
    for (i=0; i<numt-1; i++)
    {
        // Step 1: first iteration
        for (j=0; j<numhR-1; j++)
        {
            for (k=0; k<numf-1; k++)
            {
                // t_i = 2i*hdt, r_j = j*hdr, f_k = k*df
                // t_i+1/2 = (2i+1)*hdt, r_j+1/2 = (j+1)*hdr
                riph = rHalfStep*(j+1);
                // Compute u0_ti+1/2_rj+1/2
                u[k+(j+1)*idx_r+(2*i+1)*idx_t] = 
                    (u[k+(j+2)*idx_r+(2*i)*idx_t] + u[k+j*idx_r+(2*i)*idx_t]) / 2. - 
                    tHalfStep * ( cc * (u[k+(j+2)*idx_r+(2*i)*idx_t] - u[k+j*idx_r+(2*i)*idx_t]) / rStep - 
                    v[k+(j+1)*idx_r+(2*i)*idx_t] );
                // Compure u1_ti+1/2_rj+1/2
                kp1 = k+1 > numf-2 ? k+2-numf : k+1;
                kp2 = k+2 > numf-2 ? k+3-numf : k+2;
                km1 = k-1 < 0 ? k-2+numf : k-1;
                km2 = k-2 < 0 ? k-3+numf : k-2;
                d2df2U = (-u[kp2+(j+1)*idx_r+(2*i)*idx_t]+
                    16*u[kp1+(j+1)*idx_r+(2*i)*idx_t]-
                    30*u[k+(j+1)*idx_r+(2*i)*idx_t]+
                    16*u[km1+(j+1)*idx_r+(2*i)*idx_t]-
                    u[km2+(j+1)*idx_r+(2*i)*idx_t])/12/fStep/fStep;
                force = Example_compute_source(2*i, (j+1), k, tHalfStep, rHalfStep, fStep);
                v[k+(j+1)*idx_r+(2*i+1)*idx_t] = 
                    (v[k+(j+2)*idx_r+(2*i)*idx_t] + v[k+(j)*idx_r+(2*i)*idx_t]) / 2. -
                    tHalfStep * ( (-c2* (u[k+(j+2)*idx_r+(2*i)*idx_t] - u[k+(j)*idx_r+(2*i)*idx_t])/riph -
                    cc*(v[k+(j+2)*idx_r+(2*i)*idx_t] - v[k+(j)*idx_r+(2*i)*idx_t]) ) / rStep - 
                    (force + c2*d2df2U / riph / riph ) );
            }
            u[numf-1+(j+1)*idx_r+(2*i+1)*idx_t] = u[(j+1)*idx_r+(2*i+1)*idx_t];
            v[numf-1+(j+1)*idx_r+(2*i+1)*idx_t] = v[(j+1)*idx_r+(2*i+1)*idx_t];
        }

        // Step 2: fix bondary condition
        u0_mean = 0.0;
        v0_mean = 0.0;
        for (k=0; k<numf; k++)
        {
            u0_mean += u[k+idx_r+(2*i+1)*idx_t]/numf;
            v0_mean += v[k+idx_r+(2*i+1)*idx_t]/numf;
            // u[k+(2*i+1)*idx_t] = (4*u[k+idx_r+(2*i+1)*idx_t] - u[k+2*idx_r+(2*i+1)*idx_t] ) / 3.;
            u[k+(numhR-1)*idx_r+(2*i+1)*idx_t] = (4*u[k+(numhR-2)*idx_r+(2*i+1)*idx_t] - u[k+(numhR-3)*idx_r+(2*i+1)*idx_t] ) / 3.;
            // v[k+(2*i+1)*idx_t] = (4*v[k+idx_r+(2*i+1)*idx_t] - v[k+2*idx_r+(2*i+1)*idx_t] ) / 3.;
            v[k+(numhR-1)*idx_r+(2*i+1)*idx_t] = (4*v[k+(numhR-2)*idx_r+(2*i+1)*idx_t] - v[k+(numhR-3)*idx_r+(2*i+1)*idx_t] ) / 3.;
            // u[k+(2*i+1)*idx_t] = 0.0;
            // u[k+(numhR-1)*idx_r+(2*i+1)*idx_t] = 0.0;
            // v[k+(2*i+1)*idx_t] = (4*v[k+idx_r+(2*i+1)*idx_t] - v[k+2*idx_r+(2*i+1)*idx_t] ) / 3.;
            // v[k+(numhR-1)*idx_r+(2*i+1)*idx_t] = (4*v[k+(numhR-2)*idx_r+(2*i+1)*idx_t] - v[k+(numhR-3)*idx_r+(2*i+1)*idx_t] ) / 3.;
        }
        for (k=0; k<numf; k++)
        {
            force = Example_compute_source(2*i+1, 0, k, tHalfStep, rHalfStep, fStep);
            if (i==0)
                u[k+(2*i+1)*idx_t] = (force*tHalfStep2*rHalfStep2 + (-2*u[k+2*i*idx_t])*rHalfStep2-4*c2*u0_mean*tHalfStep2 ) / 
                (rHalfStep2 - 4.*c2*tHalfStep2);
            else
                u[k+(2*i+1)*idx_t] = (force*tHalfStep2*rHalfStep2 + (u[k+(2*i-1)*idx_t]-2*u[k+2*i*idx_t])*rHalfStep2-4*c2*u0_mean*tHalfStep2 ) / 
                    (rHalfStep2 - 4.*c2*tHalfStep2);
            v[k+(2*i+1)*idx_t] = v0_mean;
        }

        // Step 3: second iteration
        for (j=1; j<numhR-1; j++)
        {
            for (k=0; k<numf-1; k++)
            {
                ri = j*rHalfStep;
                // Compute u0_ti+1_j
                u[k+(j)*idx_r+(2*i+2)*idx_t] = u[k+(j)*idx_r+(2*i)*idx_t] - 
                    tStep * ( cc*(u[k+(j+1)*idx_r+(2*i+1)*idx_t] - u[k+(j-1)*idx_r+(2*i+1)*idx_t]) / rStep - 
                    v[k+(j)*idx_r+(2*i+1)*idx_t]);
                // Compure u1_ti+1/2_rj+1/2
                kp1 = k+1 > numf-1 ? k+2-numf : k+1;
                kp2 = k+2 > numf-1 ? k+3-numf : k+2;
                km1 = k-1 < 0 ? k-2+numf : k-1;
                km2 = k-2 < 0 ? k-3+numf : k-2;
                d2df2U = (-u[kp2+(j)*idx_r+(2*i+1)*idx_t]+
                    16*u[kp1+(j)*idx_r+(2*i+1)*idx_t]-
                    30*u[k+(j)*idx_r+(2*i+1)*idx_t]+
                    16*u[km1+(j)*idx_r+(2*i+1)*idx_t]-
                    u[km2+(j)*idx_r+(2*i+1)*idx_t])/12/fStep/fStep;
                force = Example_compute_source((2*i+1), j, k, tHalfStep, rHalfStep, fStep);
                // if(force > 0.0)
                //     print_debug("f[%d, %d, %d] = %f\n", i, j, k, force);
                v[k+j*idx_r+(2*i+2)*idx_t] = v[k+j*idx_r+(2*i)*idx_t] - 
                    tStep * ( ( -c2*(u[k+(j+1)*idx_r+(2*i+1)*idx_t] - u[k+(j-1)*idx_r+(2*i+1)*idx_t])/ri - 
                        cc*(v[k+(j+1)*idx_r+(2*i+1)*idx_t] - v[k+(j-1)*idx_r+(2*i+1)*idx_t]) ) / rStep - 
                        (force + c2*d2df2U / ri / ri));
            }
            u[numf-1+(j)*idx_r+(2*i+2)*idx_t] = u[(j)*idx_r+(2*i+2)*idx_t];
            v[numf-1+(j)*idx_r+(2*i+2)*idx_t] = v[(j)*idx_r+(2*i+2)*idx_t];
        }
        u0_mean = 0.0;
        v0_mean = 0.0;
        // Step 4: fix bondary condition
        for (k=0; k<numf; k++)
        {
            u0_mean += u[k+idx_r+(2*i+2)*idx_t]/numf;
            v0_mean += v[k+idx_r+(2*i+2)*idx_t]/numf;
            // u[k+(2*i+2)*idx_t] = (4*u[k+idx_r+(2*i+2)*idx_t] - u[k+2*idx_r+(2*i+2)*idx_t] ) / 3.;
            u[k+(numhR-1)*idx_r+(2*i+2)*idx_t] = (4*u[k+(numhR-2)*idx_r+(2*i+2)*idx_t] - u[k+(numhR-3)*idx_r+(2*i+2)*idx_t] ) / 3.;
            // v[k+(2*i+2)*idx_t] = (4*v[k+idx_r+(2*i+2)*idx_t] - v[k+2*idx_r+(2*i+2)*idx_t] ) / 3.;
            v[k+(numhR-1)*idx_r+(2*i+2)*idx_t] = (4*v[k+(numhR-2)*idx_r+(2*i+2)*idx_t] - v[k+(numhR-3)*idx_r+(2*i+2)*idx_t] ) / 3.;
            // u[k+(2*i+2)*idx_t] = 0.0;
            // u[k+(numhR-1)*idx_r+(2*i+2)*idx_t] = 0.0;
            // v[k+(2*i+2)*idx_t] = (4*v[k+idx_r+(2*i+2)*idx_t] - v[k+2*idx_r+(2*i+2)*idx_t] ) / 3.;
            // v[k+(numhR-1)*idx_r+(2*i+2)*idx_t] = (4*v[k+(numhR-2)*idx_r+(2*i+2)*idx_t] - v[k+(numhR-3)*idx_r+(2*i+2)*idx_t] ) / 3.;
        }
        for (k=0; k<numf; k++)
        {
            force = Example_compute_source(2*i+2, 0, k, tHalfStep, rHalfStep, fStep);
            u[k+(2*i+2)*idx_t] = (force*tHalfStep2*rHalfStep2 + (u[k+(2*i)*idx_t]-2*u[k+(2*i+1)*idx_t])*rHalfStep2-4*c2*u0_mean*tHalfStep2 ) / 
                (rHalfStep2 - 4.*c2*tHalfStep2);
            v[k+(2*i+2)*idx_t] = v0_mean;
        }
    }

    /**
     *
     * Interpolation
     * 
     */
    arr = CreateREAL8Array(4, 2, numt, numr, numf);
    for (i=1; i<numt; i++)
    {
        for (j=0; j<numr; j++)
        {
            for (k=0; k<numf; k++)
            {
                arr->data[k + j*idx_r + i*idx_t_itp] = buffer->data[k + 2*j*idx_r + 2*i*idx_t];
                // print_debug("arr[%d, %d, %d] = %f\n", k, j, i, buffer->data[k + 2*j*idx_r + 2*i*idx_t]);
                arr->data[k + j*idx_r + i*idx_t_itp + idx_u_itp] = buffer->data[k + 2*j*idx_r + 2*i*idx_t + idx_u];
            }
        }
    }
    // dump data
    tCoord = CreateREAL8Array(1, numt);
    rCoord = CreateREAL8Array(1, numr);
    fCoord = CreateREAL8Array(1, numf);
    for (i=0;i<numt;i++)
        tCoord->data[i] = i * tStep;
    for (i=0;i<numr;i++)
        rCoord->data[i] = i * rStep;
    for (i=0;i<numf;i++)
        fCoord->data[i] = i * fStep;
    sprintf(fsave, "%s/%s", params->prefix, "plane_osci.h5");
    DumpREAL8ArrayTohdf5(fsave, "data", arr, 1);
    DumpREAL8ArrayTohdf5(fsave, "tCoord", tCoord, 0);
    DumpREAL8ArrayTohdf5(fsave, "rCoord", rCoord, 0);
    DumpREAL8ArrayTohdf5(fsave, "fCoord", fCoord, 0);

QUIT:
    STRUCTFREE(buffer, REAL8Array);
    STRUCTFREE(arr, REAL8Array);
    STRUCTFREE(tCoord, REAL8Array);
    STRUCTFREE(rCoord, REAL8Array);
    STRUCTFREE(fCoord, REAL8Array);
    if (failed)
        return CEV_FAILURE;
    return CEV_SUCCESS;
}
