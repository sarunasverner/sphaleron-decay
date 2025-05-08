/*******************************************************************************
 * sphaleron_solve.c
 * Written by Sarunas Verner
 *
 * Solves eq1-eq5 in radial coordinates, storing:
 *   fA(r,t), fB(r,t), fC(r,t), H(r,t), KK(r,t)
 *   plus their radial derivatives: fA'(r,t), fB'(r,t), fC'(r,t), H'(r,t), K'(r,t)
 *   plus their time derivatives: dot{fA}(r,t), dot{fB}(r,t), dot{fC}(r,t), dot{H}(r,t), dot{K}(r,t)
 *
 * Domain: r in [0.0001, 30], Nx=3001 points
 * Time:   t in [0,50], integrated with dt=0.0005 (RK4)
 * Output: every 0.1 seconds => 5001 snapshots
 *
 * We'll store and then print all arrays.
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* ------------------------------------------------------------------ */
/* Physical parameters (dimensionless units: v = 1)                   */
/* Change MW and MH below to explore different mass ratios.           */
/* ------------------------------------------------------------------ */
#define MW  1.0       /* W-boson mass  (current default = 1) */
#define MH  1.556     /* Higgs mass    (current default = 1.556) */
#define MW2 (MW*MW)
#define MH2 (MH*MH)


/* -- Radial domain -- */
#define RMIN  0.001
#define RMAX  30.0
#define Nx    3001

/* -- Time domain & stepping -- */
#define TSTART  0.0
#define TFINAL  50.0
#define DT      0.00001

/* We'll output every 0.1 seconds => 101 times */
#define TOUT    0.1
#define NT_OUT  ((int)((TFINAL - TSTART)/TOUT + 0.5) + 1)

/* Global arrays for the real-time solution (RK4) */
static double rVals[Nx];
static double dr;

/* The fields + velocities at the "current" time: */
static double fA[Nx], vA[Nx];
static double fB[Nx], vB[Nx];
static double fC[Nx], vC[Nx];
static double H [Nx], vH[Nx];
static double KK[Nx], vK[Nx];

/* 2D arrays to store solution at each output time */
static double tOutVals[NT_OUT];  /* times at which we store data */

/* Fields: fA, fB, fC, H, KK */
static double fAOut  [NT_OUT][Nx];
static double fBOut  [NT_OUT][Nx];
static double fCOut  [NT_OUT][Nx];
static double HOut   [NT_OUT][Nx];
static double KKOut  [NT_OUT][Nx];

/* Radial derivatives: fA', fB', fC', H', KK' */
static double fAprimeOut [NT_OUT][Nx];
static double fBprimeOut [NT_OUT][Nx];
static double fCprimeOut [NT_OUT][Nx];
static double HprimeOut  [NT_OUT][Nx];
static double KKprimeOut [NT_OUT][Nx];

/* Time derivatives: vA = dot{fA}, vB = dot{fB}, etc. */
static double vAOut [NT_OUT][Nx];
static double vBOut [NT_OUT][Nx];
static double vCOut [NT_OUT][Nx];
static double vHOut [NT_OUT][Nx];
static double vKOut [NT_OUT][Nx];


/* sech function */
static double sech(double x) {
    return 1.0 / cosh(x);
}

/*******************************************************************************
 * PDE boundary functions for fA, KK
 ******************************************************************************/
static double fA_boundary(double r) {
    return 1.0 - 2*(1 - sech(1.150155*r));
}

static double KK_boundary(double r) {
    return tanh(1.0463536*r);
}

/* New unstable mode profiles for the remaining fields */
static double fB_boundary(double r) {
    return sqrt(0.01) * 10.2199 * exp(-2.74187*r) * pow(r, 1.6506);
}

static double fC_boundary(double r) {
    return (sqrt(2)/r) * sqrt(0.01) * 7.48017 * exp(-2.86292*r) * pow(r, 1.67796);
}

static double H_boundary(double r) {
    return (1/sqrt(2)) * (1/r) * sqrt(0.01) * 4.00337 * exp(-2.62178*r) * pow(r, 1.5913);
}


/*******************************************************************************
 * Initial conditions
 ******************************************************************************/
static void setInitialConditions() {
    for(int i=0; i<Nx; i++){
        double rr = rVals[i];
        /* fA, KK from boundary expressions across domain: */
        fA[i] = fA_boundary(rr);
        KK[i] = KK_boundary(rr);

        /* fB, fC, H from exponent/r^power profiles */
        fB[i] = sqrt(0.01) *  10.243 * exp(-2.7439*rr) * pow(rr, 1.652);
        fC[i] = (sqrt(2)/rr) * sqrt(0.01) * 7.52242 * exp(-2.8682*rr) * pow(rr, 1.68214);
        H[i] = (1/sqrt(2)) * (1/rr) * sqrt(0.01) * 4.0121  * exp(-2.62378*rr) * pow(rr, 1.59307);


        
        /* Zero velocities at t=0 */
        vA[i] = 0.0;
        vB[i] = 0.0;
        vC[i] = 0.0;
        vH[i] = 0.0;
        vK[i] = 0.0;
    }
}

/*******************************************************************************
 * Finite difference helpers
 ******************************************************************************/
static inline double first_derivative(double *f, int i){
    /* Central difference for i in [1..Nx-2], handle boundary i=0 or i=Nx-1 carefully */
    return (f[i+1] - f[i-1])/(2.0*dr);
}
static inline double first_derivative_boundary(double *f, int i){
    /* If i=0, forward difference; if i=Nx-1, backward difference. */
    if(i==0){
        return (f[i+1] - f[i])/dr;
    } else if(i == Nx-1){
        return (f[i] - f[i-1])/dr;
    } else {
        return first_derivative(f,i);
    }
}
static inline double second_derivative(double *f, int i){
    return (f[i+1] - 2.0*f[i] + f[i-1])/(dr*dr);
}
static inline double second_deriv_r_times_field(double *f, int i){
    double rp_i   = rVals[i]*f[i];
    double rp_ip1 = rVals[i+1]*f[i+1];
    double rp_im1 = rVals[i-1]*f[i-1];
    return (rp_ip1 - 2.0*rp_i + rp_im1)/(dr*dr);
}

/*******************************************************************************
 * Compute PDE Right-Hand Side => partial_t vA = fA_tt, etc.
 ******************************************************************************/
static void computePDERHS(
    double *fA_, double *vA_,
    double *fB_, double *vB_,
    double *fC_, double *vC_,
    double *H_,  double *vH_,
    double *KK_, double *vK_,
    double *dfA, double *dvA,
    double *dfB, double *dvB,
    double *dfC, double *dvC,
    double *dH,  double *dvH,
    double *dKK, double *dvK
){
    /* 1) Dirichlet BC for fA, KK at r=RMIN, RMAX */
    fA_[0]  = fA_boundary(rVals[0]);
    KK_[0]  = KK_boundary(rVals[0]);
    vA_[0]  = 0.0;
    vK_[0]  = 0.0;

    fA_[Nx-1] = fA_boundary(rVals[Nx-1]);
    KK_[Nx-1] = KK_boundary(rVals[Nx-1]);
    vA_[Nx-1] = 0.0;
    vK_[Nx-1] = 0.0;

    /* For fB, fC, H => 0 at boundaries, also zero velocities. */
    fB_[0] = fB_boundary(rVals[0]);
    fC_[0] = fC_boundary(rVals[0]);
    H_[0] = H_boundary(rVals[0]); 
    
    vB_[0] = 0.0;
    vC_[0] = 0.0;
    vH_[0] = 0.0;
    
       

    fB_[Nx-1] = fB_boundary(rVals[Nx-1]);
    fC_[Nx-1] = fC_boundary(rVals[Nx-1]);
    H_[Nx-1]  = H_boundary(rVals[Nx-1]);
    
    
    vB_[Nx-1] = 0.0;
    vC_[Nx-1] = 0.0;
    vH_[Nx-1] = 0.0;

    /* 2) PDE interior: i=1..Nx-2 */
    for(int i=1; i<Nx-1; i++){
        double r = rVals[i];

        /* eq1 => dvA[i] = fA_tt */
        double fA_rr = second_derivative(fA_, i);
        double dfB_dr = first_derivative(fB_, i);
        double dfC_dr = first_derivative(fC_, i);

        double eq1RHS = fA_rr
            - ((fA_[i]*fA_[i] + fB_[i]*fB_[i] - 1.0)*fA_[i])/(r*r)
            - (MW2*((H_[i]*H_[i] + KK_[i]*KK_[i])*fA_[i] + KK_[i]*KK_[i] - H_[i]*H_[i]))
            - (fA_[i]*fC_[i]*fC_[i])
            + 2.0*dfB_dr*fC_[i]
            + fB_[i]*dfC_dr;

        dvA[i] = eq1RHS;

        /* eq2 => dvB[i] = fB_tt */
        double fB_rr = second_derivative(fB_, i);
        double dfA_dr = first_derivative(fA_, i);

        double eq2RHS = fB_rr
            - ((fA_[i]*fA_[i] + fB_[i]*fB_[i] - 1.0)*fB_[i])/(r*r)
            - (MW2*((H_[i]*H_[i] + KK_[i]*KK_[i])*fB_[i] - 2.0*KK_[i]*H_[i]))
            - (fB_[i]*fC_[i]*fC_[i])
            - 2.0*dfA_dr*fC_[i]
            - fA_[i]*dfC_dr;

        dvB[i] = eq2RHS;

        /* eq3 => dvC[i] = fC_tt = -[...] */
        double dH_dr = first_derivative(H_, i);
        double dK_dr = first_derivative(KK_, i);

        double eq3RHS = -(
            (2.0/(r*r))*(fA_[i]*fA_[i] + fB_[i]*fB_[i])*fC_[i]
            + (MW2*(H_[i]*H_[i] + KK_[i]*KK_[i])*fC_[i])
            + (2.0*MW2)*( dH_dr*KK_[i] - H_[i]*dK_dr )
            + (2.0/(r*r))*( dfA_dr*fB_[i] - fA_[i]*dfB_dr )
        );
        
        dvC[i] = eq3RHS;

        /* eq4 => dvH[i] = H_tt = (1.0/r)*d^2/dr^2(rH) - bracket_4 */
        double d2_rH = second_deriv_r_times_field(H_, i);
        double dfCdr = dfC_dr;
        double bracket_4 =
            (1.0/(2.0*r*r))*(fA_[i]*fA_[i] + fB_[i]*fB_[i] + 1.0)*H_[i]
            - (1.0/(r*r))*(H_[i]*fA_[i] + KK_[i]*fB_[i])
            + 0.5*MH2*(H_[i]*H_[i] + KK_[i]*KK_[i] - 1.0)*H_[i]
            - (1.0/r)*KK_[i]*fC_[i]
            + 0.25*H_[i]*fC_[i]*fC_[i]
            - fC_[i]*dK_dr
            - 0.5*KK_[i]*dfCdr;

        dvH[i] = (1.0/r)*d2_rH - bracket_4;

        /* eq5 => dvK[i] = KK_tt = (1.0/r)*d^2/dr^2(r*KK) - bracket_5 */
        double d2_rK = second_deriv_r_times_field(KK_, i);
        double bracket_5 =
            (1.0/(2.0*r*r))*(fA_[i]*fA_[i] + fB_[i]*fB_[i] + 1.0)*KK_[i]
            + (1.0/(r*r))*( KK_[i]*fA_[i] - H_[i]*fB_[i] )
            + 0.5*MH2*(H_[i]*H_[i] + KK_[i]*KK_[i] - 1.0)*KK_[i]
            + (1.0/r)*H_[i]*fC_[i]
            + 0.25*KK_[i]*fC_[i]*fC_[i]
            + fC_[i]*dH_dr
            + 0.5*H_[i]*dfCdr;

        dvK[i] = (1.0/r)*d2_rK - bracket_5;
    }

    /* partial_t fA[i] = vA_[i], etc. */
    for(int i=0; i<Nx; i++){
        dfA[i] = vA_[i];
        dfB[i] = vB_[i];
        dfC[i] = vC_[i];
        dH [i] = vH_[i];
        dKK[i] = vK_[i];
    }

    /* BC points => pinned => zero derivatives. */
    dfA[0] = 0.0;  dvA[0] = 0.0;
    dfA[Nx-1] = 0.0; dvA[Nx-1] = 0.0;

    dfB[0] = 0.0;  dvB[0] = 0.0;
    dfB[Nx-1] = 0.0; dvB[Nx-1] = 0.0;

    dfC[0] = 0.0;  dvC[0] = 0.0;
    dfC[Nx-1] = 0.0; dvC[Nx-1] = 0.0;

    dH[0] = 0.0;   dvH[0] = 0.0;
    dH[Nx-1] = 0.0; dvH[Nx-1] = 0.0;

    dKK[0] = 0.0;  dvK[0] = 0.0;
    dKK[Nx-1] = 0.0; dvK[Nx-1] = 0.0;
}

/*******************************************************************************
 * 6) Compute radial derivatives f'(r) for each field
 ******************************************************************************/
static void computeRadialDerivatives(
    double *fA_, double *fB_, double *fC_, double *H_, double *KK_,
    double *fAprime, double *fBprime, double *fCprime, double *Hprime, double *KKprime
){
    for(int i=0; i<Nx; i++){
        fAprime[i] = first_derivative_boundary(fA_, i);
        fBprime[i] = first_derivative_boundary(fB_, i);
        fCprime[i] = first_derivative_boundary(fC_, i);
        Hprime[i]  = first_derivative_boundary(H_,  i);
        KKprime[i] = first_derivative_boundary(KK_, i);
    }
}

/*******************************************************************************
 * 7) RK4 integrator step
 ******************************************************************************/
static double k1fA[Nx], k1vA[Nx], k1fB[Nx], k1vB[Nx], k1fC[Nx], k1vC[Nx], k1H[Nx], k1vH[Nx], k1K[Nx], k1vK[Nx];
static double k2fA[Nx], k2vA[Nx], k2fB[Nx], k2vB[Nx], k2fC[Nx], k2vC[Nx], k2H[Nx], k2vH[Nx], k2K[Nx], k2vK[Nx];
static double k3fA[Nx], k3vA[Nx], k3fB[Nx], k3vB[Nx], k3fC[Nx], k3vC[Nx], k3H[Nx], k3vH[Nx], k3K[Nx], k3vK[Nx];
static double k4fA[Nx], k4vA[Nx], k4fB[Nx], k4vB[Nx], k4fC[Nx], k4vC[Nx], k4H[Nx], k4vH[Nx], k4K[Nx], k4vK[Nx];

static double tmpfA[Nx], tmpvA[Nx], tmpfB[Nx], tmpvB[Nx], tmpfC[Nx], tmpvC[Nx], tmpH[Nx], tmpvH[Nx], tmpK[Nx], tmpvK[Nx];

static void rk4_step(double dt){
    /* 1) k1 */
    computePDERHS(fA, vA, fB, vB, fC, vC, H, vH, KK, vK,
                  k1fA, k1vA, k1fB, k1vB, k1fC, k1vC, k1H, k1vH, k1K, k1vK);

    /* 2) k2 */
    for(int i=0; i<Nx; i++){
        tmpfA[i] = fA[i] + 0.5*dt*k1fA[i];
        tmpvA[i] = vA[i] + 0.5*dt*k1vA[i];
        tmpfB[i] = fB[i] + 0.5*dt*k1fB[i];
        tmpvB[i] = vB[i] + 0.5*dt*k1vB[i];
        tmpfC[i] = fC[i] + 0.5*dt*k1fC[i];
        tmpvC[i] = vC[i] + 0.5*dt*k1vC[i];
        tmpH[i]  = H[i]  + 0.5*dt*k1H[i];
        tmpvH[i] = vH[i] + 0.5*dt*k1vH[i];
        tmpK[i]  = KK[i] + 0.5*dt*k1K[i];
        tmpvK[i] = vK[i] + 0.5*dt*k1vK[i];
    }
    computePDERHS(tmpfA, tmpvA, tmpfB, tmpvB, tmpfC, tmpvC, tmpH, tmpvH, tmpK, tmpvK,
                  k2fA, k2vA, k2fB, k2vB, k2fC, k2vC, k2H, k2vH, k2K, k2vK);

    /* 3) k3 */
    for(int i=0; i<Nx; i++){
        tmpfA[i] = fA[i] + 0.5*dt*k2fA[i];
        tmpvA[i] = vA[i] + 0.5*dt*k2vA[i];
        tmpfB[i] = fB[i] + 0.5*dt*k2fB[i];
        tmpvB[i] = vB[i] + 0.5*dt*k2vB[i];
        tmpfC[i] = fC[i] + 0.5*dt*k2fC[i];
        tmpvC[i] = vC[i] + 0.5*dt*k2vC[i];
        tmpH[i]  = H[i]  + 0.5*dt*k2H[i];
        tmpvH[i] = vH[i] + 0.5*dt*k2vH[i];
        tmpK[i]  = KK[i] + 0.5*dt*k2K[i];
        tmpvK[i] = vK[i] + 0.5*dt*k2vK[i];
    }
    computePDERHS(tmpfA, tmpvA, tmpfB, tmpvB, tmpfC, tmpvC, tmpH, tmpvH, tmpK, tmpvK,
                  k3fA, k3vA, k3fB, k3vB, k3fC, k3vC, k3H, k3vH, k3K, k3vK);

    /* 4) k4 */
    for(int i=0; i<Nx; i++){
        tmpfA[i] = fA[i] + dt*k3fA[i];
        tmpvA[i] = vA[i] + dt*k3vA[i];
        tmpfB[i] = fB[i] + dt*k3fB[i];
        tmpvB[i] = vB[i] + dt*k3vB[i];
        tmpfC[i] = fC[i] + dt*k3fC[i];
        tmpvC[i] = vC[i] + dt*k3vC[i];
        tmpH[i]  = H[i]  + dt*k3H[i];
        tmpvH[i] = vH[i] + dt*k3vH[i];
        tmpK[i]  = KK[i] + dt*k3K[i];
        tmpvK[i] = vK[i] + dt*k3vK[i];
    }
    computePDERHS(tmpfA, tmpvA, tmpfB, tmpvB, tmpfC, tmpvC, tmpH, tmpvH, tmpK, tmpvK,
                  k4fA, k4vA, k4fB, k4vB, k4fC, k4vC, k4H, k4vH, k4K, k4vK);

    /* Combine => y_{n+1} = y_n + dt/6*(k1 + 2k2 + 2k3 + k4) */
    for(int i=0; i<Nx; i++){
        fA[i] += (dt/6.0)*(k1fA[i] + 2.0*k2fA[i] + 2.0*k3fA[i] + k4fA[i]);
        vA[i] += (dt/6.0)*(k1vA[i] + 2.0*k2vA[i] + 2.0*k3vA[i] + k4vA[i]);

        fB[i] += (dt/6.0)*(k1fB[i] + 2.0*k2fB[i] + 2.0*k3fB[i] + k4fB[i]);
        vB[i] += (dt/6.0)*(k1vB[i] + 2.0*k2vB[i] + 2.0*k3vB[i] + k4vB[i]);

        fC[i] += (dt/6.0)*(k1fC[i] + 2.0*k2fC[i] + 2.0*k3fC[i] + k4fC[i]);
        vC[i] += (dt/6.0)*(k1vC[i] + 2.0*k2vC[i] + 2.0*k3vC[i] + k4vC[i]);

        H[i]  += (dt/6.0)*(k1H[i] + 2.0*k2H[i] + 2.0*k3H[i] + k4H[i]);
        vH[i] += (dt/6.0)*(k1vH[i] + 2.0*k2vH[i] + 2.0*k3vH[i] + k4vH[i]);

        KK[i] += (dt/6.0)*(k1K[i] + 2.0*k2K[i] + 2.0*k3K[i] + k4K[i]);
        vK[i] += (dt/6.0)*(k1vK[i] + 2.0*k2vK[i] + 2.0*k3vK[i] + k4vK[i]);
    }
}

/*******************************************************************************
 * 8) Main
 ******************************************************************************/
int main(void){
    /* 1) Build radial array, set dr */
    dr = (RMAX - RMIN)/(Nx - 1);
    for(int i=0; i<Nx; i++){
        rVals[i] = RMIN + i*dr;
    }

    /* 2) Set initial conditions */
    setInitialConditions();

    /* 3) Build array of times at which we'll store the solution */
    for(int n=0; n<NT_OUT; n++){
        tOutVals[n] = TSTART + n*TOUT;  /* e.g. 0,0.1,0.2,...,10.0 */
    }

    /* 4) We'll do the integration from t=0..10, store solutions at tOutVals */
    int outIndex = 0;
    double nextOutTime = tOutVals[outIndex];
    double t = TSTART;

    /* (A) Store the initial solution (t=0) right away */
    {
        /* Compute radial derivatives for the initial fields */
        static double fAprime[Nx], fBprime[Nx], fCprime[Nx], Hprime[Nx], KKprime[Nx];
        computeRadialDerivatives(fA, fB, fC, H, KK, fAprime, fBprime, fCprime, Hprime, KKprime);

        /* Copy everything into the 2D arrays at outIndex=0 */
        for(int i=0; i<Nx; i++){
            fAOut [outIndex][i] = fA[i];
            fBOut [outIndex][i] = fB[i];
            fCOut [outIndex][i] = fC[i];
            HOut  [outIndex][i] = H [i];
            KKOut [outIndex][i] = KK[i];

            fAprimeOut[outIndex][i] = fAprime[i];
            fBprimeOut[outIndex][i] = fBprime[i];
            fCprimeOut[outIndex][i] = fCprime[i];
            HprimeOut [outIndex][i] = Hprime[i];
            KKprimeOut[outIndex][i] = KKprime[i];

            vAOut[outIndex][i] = vA[i];
            vBOut[outIndex][i] = vB[i];
            vCOut[outIndex][i] = vC[i];
            vHOut[outIndex][i] = vH[i];
            vKOut[outIndex][i] = vK[i];
        }
    }
    outIndex++;
    if(outIndex < NT_OUT){
        nextOutTime = tOutVals[outIndex];
    } else {
        nextOutTime = 1.0e30;
    }

    /* (B) Main time loop: integrate until TFINAL, store solutions every 0.1s */
    while(t < TFINAL){
        /* Do one RK4 step: */
        rk4_step(DT);
        t += DT;

        /* Check if we've passed the next output time */
        if(t >= nextOutTime - 1.0e-12 && outIndex < NT_OUT){
            /* Compute radial derivatives for the new fields */
            static double fAprime[Nx], fBprime[Nx], fCprime[Nx], Hprime[Nx], KKprime[Nx];
            computeRadialDerivatives(fA, fB, fC, H, KK, fAprime, fBprime, fCprime, Hprime, KKprime);

            /* Store them in the 2D arrays at outIndex */
            for(int i=0; i<Nx; i++){
                fAOut [outIndex][i] = fA[i];
                fBOut [outIndex][i] = fB[i];
                fCOut [outIndex][i] = fC[i];
                HOut  [outIndex][i] = H [i];
                KKOut [outIndex][i] = KK[i];

                fAprimeOut[outIndex][i] = fAprime[i];
                fBprimeOut[outIndex][i] = fBprime[i];
                fCprimeOut[outIndex][i] = fCprime[i];
                HprimeOut [outIndex][i] = Hprime[i];
                KKprimeOut[outIndex][i] = KKprime[i];

                vAOut[outIndex][i] = vA[i];
                vBOut[outIndex][i] = vB[i];
                vCOut[outIndex][i] = vC[i];
                vHOut[outIndex][i] = vH[i];
                vKOut[outIndex][i] = vK[i];
            }

            outIndex++;
            if(outIndex < NT_OUT){
                nextOutTime = tOutVals[outIndex];
            } else {
                nextOutTime = 1.0e30; /* no more output times */
            }
        }
    }

    /* 5) Print out the stored data */
    /* For each output time => we print Nx lines with
     *  r, fA, fB, fC, H, KK, fA', fB', fC', H', KK', vA, vB, vC, vH, vK
     */
    for(int n=0; n<outIndex; n++){
        double tt = tOutVals[n];
        printf("# time = %.4f\n", tt);
        for(int i=0; i<Nx; i++){
            double rr = rVals[i];
            printf("%.10g  "   /* r */
                   "%.10g  "   /* fA */
                   "%.10g  "   /* fB */
                   "%.10g  "   /* fC */
                   "%.10g  "   /* H */
                   "%.10g  "   /* KK */
                   "%.10g  "   /* fA' */
                   "%.10g  "   /* fB' */
                   "%.10g  "   /* fC' */
                   "%.10g  "   /* H' */
                   "%.10g  "   /* KK' */
                   "%.10g  "   /* vA = dot{fA} */
                   "%.10g  "   /* vB */
                   "%.10g  "   /* vC */
                   "%.10g  "   /* vH */
                   "%.10g\n",  /* vK */
                   rr,
                   fAOut[n][i],
                   fBOut[n][i],
                   fCOut[n][i],
                   HOut[n][i],
                   KKOut[n][i],
                   fAprimeOut[n][i],
                   fBprimeOut[n][i],
                   fCprimeOut[n][i],
                   HprimeOut[n][i],
                   KKprimeOut[n][i],
                   vAOut[n][i],
                   vBOut[n][i],
                   vCOut[n][i],
                   vHOut[n][i],
                   vKOut[n][i]
            );
        }
    }

    return 0;
}
