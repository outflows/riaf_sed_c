#include "libs.h"

double get_rdomega(double cs, double Omega, double y[4])
{
    double ww, w1, w2, w3, rdomega;

    ww = -alpha/y[1]*cs*cs;
    w1 = -ww/y[1];
    w2 = ww/(y[2]/MUE + y[3]/MUI);
    w3 = ww/y[0] - 2.*Omega;

    rdomega = w1*dy[1] + w2*(dy[2]/MUE + dy[3]/MUI) + w3;
    if(rdomega > 0.) {
        fprintf(stderr, "ERROR: rdomega > 0.\n");
        exit(1);
    }

    return(rdomega);
}


void get_advec(double rho, double p0, double y[4], double gam, double k2,
    double sint, double advec, double advec2)
{

    double dlnok, a1, a2, a3, rhodr;
    double Thetae, Thetai, xe, xi, sumTh;
    double K1, K2, K3;
    double Qadv, ate, Qadv_e, Qie, Qvis;

    dlnok = -0.5/y[0] - 1./(y[0] - 2.);
    a1 = -1./y[1];
    a2 = -0.5/(y[2]/MUE + y[3]/MUI));
    //       a3=dlnok-1./y[0];
    a3 = dlnok + (p0 - 1.)/y[0];
    rhodr = a1*dy[1] + a2*(1./MUI*dy[3] + 1./MUE*dy[2]) + a3;

    Qadv = rho*y[1]*RU*(dy[3]/(MUI*(gam - 1.)) - y[3]/MUI*rhodr); // rate of energy advection
    //Qie = 37.85/MBH*rho*rho*(y[3] - y[2])*(sqrt(2./M_PI) + sqrt(sumTh))/pow(sumTh, 1.5); // where does this come from?
    Qie = coulomb_cooling(y[2], y[3]);
    Qvis = Qadv + Qie; // rate of viscous heating
    advec = Qadv/Qvis; // f(r) in Yuan+2005 (eq 3)

    Thetae = (KBOL/(ME*CL*CL))*y[2]; // dimensionless electron temperature
    Thetai = (KBOL/(MP*CL*CL))*y[3]; // dimensionless ion temperature
    xe = 1./Thetae;
    xi = 1./Thetai;
    sumTh = Thetae + Thetai;

    K1 = gsl_sf_bessel_Kn(1, xe);
    K2 = gsl_sf_bessel_Kn(2, xe);
    K3 = gsl_sf_bessel_Kn(3, xe);

    ate = xe*((3.*K3 + K1)/(4.*K2) - 1.);
    Qadv_e = rho*y[1]*RU*(ate/MUE*dy[2] - MUE*y[2]*rhodr); // electron advection
    advec2 = (Qadv + Qadve)/Qvis;

    // CASE: s is NOT equal to 0., i.e., outflow's case
    if (advec < 0.) {
        p0 = 0.;
    } else if(advec > 0.) {
        p0 = advec*pp0;
    }

    // the following calculates the integration for 's'
    // sint denotes the integration for 's'
    if (k2 == 0) {
        p0 = 0.;
    }
    sint = sint + p0/y[0]*H_disc;

    if (k2 == 0) {
        Theta = 2.3536744d0;
        ate = 2.73821d0 + 0.12773*(Thetae - Theta) -
              0.0901154*pow((Thetae - Theta), 2.) +
              0.0227991*pow((Thetae - Theta), 3.) -
              0.00219443*pow((Thetae - Theta), 4.) +
              7.07116d-5*pow((Thetae - Theta), 5.);
        //ate=1./Thetae*((3.*mbsl4(3,1.d0/Thetae)+mbsl4(1,1.d0/Thetae))/4./mbsl4(2,1.d0/Thetae)-1.d0)
        Qadv_e = rho*y[1]*RU*(ate/MUE*dy[2] - MUE*y[2]*rhodr);
    }

}
