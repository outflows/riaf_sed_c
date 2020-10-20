double get_jnu(double Thetae, double Ne, double nu)
{
    double jnutotal, jnus, jnub;

    double jnu_synch(double Thetae, double Ne, double nu);
    double jnu_brems(double Thetae, double Ne, double nu);

    jnus = jnu_synch(Thetae, Ne, nu);
    jnub = jnu_brems(Thetae, Ne, nu);

    jnutotal = jnus + jnub;

    return(jnutotal);
}

// BREMSSTRAHLUNG ROUTINES
double cooling_brems(double Thetae, double Ne)
{
    /*
        Returns the cooling rate per unit volume for bremsstrahlung.
    */

    double Fei, qee, qie, qbr;
    double fac1 = Ne*Ne*CL*RE*RE*ALPHA_F*CL*CL;
    double fac2 = (20./(9.*sqrt(M_PI)))*(44. - 3.*M_PI*M_PI);
    double fac3 = 1.25*Ne*Ne*SIGMA_THOMSON*CL*ALPHA_F*ME*CL*CL;

    if (Thetae < 1.) {
        Fei = 4.*sqrt(2.*Thetae/pow(M_PI,3.))*(1. + 1.781*pow(Thetae,1.34)); // eq. 29 of MANMOTO+1997
        //qee = 1.73*pow(Thetae, 1.5);
        qee = fac1*fac2*pow(Thetae, 1.5);
        qee *= (1. + 1.1*Thetae + Thetae*Thetae - 1.25*pow(Thetae, 2.5)); // eq. 31 of MANMOTO+1997
        qie = fac3*Fei; // eq. 28 of MANMOTO+1997
    } else if (Thetae > 1.) {
        Fei = 9.*Thetae/(2.*M_PI)*(log(1.123*Thetae + 0.48) + 1.5); // eq. 30 of MANMOTO+1997
        qee = fac1*24.*Thetae*(log(1.123*Thetae) + 1.28); // eq. 32 of MANMOTO+1997
        qei = fac3*Fei; // eq. 28 of MANMOTO+1997
    }
    qbr = qie + qee; // eq. 27 of MANMOTO+1997

    return(qbr);
}

double jnu_brems(double Thetae, double Ne, double nu)
{
    /*
        Returns bremsstrahlung emissivity, eq. 33 of MANMOTO+1997
    */

    double Te, x, invx, gaunt, jnu, qbr;
    double cooling_brems(double Thetae, double Ne);

    Te = ME*CL*CL*Thetae/KBOL;
    x = HPL*nu/(KBOL*Te);
    invx = KBOL*Te/(HPL*nu);

    if (invx < 1) {
        gaunt = HPL/(KBOL*Te) * sqrt(3.*invx/M_PI); // eq. 34 of MANMOTO+1997
    } else {
        gaunt = HPL/(KBOL*Te) * \sqrt(3)/M_PI * log(4.*invx/1.781); // eq. 35 of MANMOTO+1997
    }

    qbr = cooling_brems(Thetae, Ne);
    jnu = qbr*gaunt*exp(x);

    return(jnu);
}

// SYNCHROTRON ROUTINES
double jnu_synch(double Thetae, double Ne, double nu)
{
    /*
        Returns synchrotron emissivity jnu, eq. 36 of MANMOTO+1997
    */

    double Te, Bmag, x;
    double x12, x13, x14, x16, invThetae;
    double I, K2, jnu;

    invThetae = 1./Thetae;
    Te = ME*CL*CL*Thetae/KBOL;
    Bmag = sqrt(8.*M_PI*Ne*KBOL*Te/beta);
    x = 4.*M_PI*ME*CL*nu/(3.*ELECTRON*Bmag*Thetae*Thetae);
    x12 = pow(x, 1./2.);
    x13 = pow(x, 1./3.);
    x14 = pow(x, 1./4.);
    x16 = pow(x, 1./6.);
    I = (4.0505/x16); // eq. 37 of MANMOTO+1997
    I *= (1. + 0.4/x14 + 0.5316/x12); // eq. 37 of MANMOTO+1997
    I *= exp(-1.8899*x13); // eq. 37 of MANMOTO+1997

    K2 = gsl_sf_bessel_Kn(2, 1. / Thetae); // use function here!

    jnu = 4.43e-30*(4.*M_PI*Ne*nu/K2)*I; // eq. 36 of MANMOTO+1997

    return(jnu);
}


// COMPTON
double get_eta(double Thetae, double Ne, double nu, double h, double kappanu, double taunu)
{
    /*
        This returns the energy enhancement factor eta, defined
        in equation 38 of MANMOTO+1997
    */

    double s, a_temp, eta_max, jm, tau_eff, tau_nu, max, tau_es;
    double Qgamma1, Qgamma2, eta;

    a_temp = 1. + 4.*Thetae + 16.*Thetae*Thetae; // eq. 39 of MANMOTO+1997
    eta_max = 3.*KBOL*Thetae/(HPL*nu); // eq. 39 of MANMOTO+1997
    jm = log(eta_max)/log(a_temp); // eq. 39 of MANMOTO+1997

    tau_eff = taunu*sqrt(1. + Ne*SIGMA_THOMSON/kappanu);
    max = 1.;
    if (1./tau_eff > max) {
        max = 1./tau_eff;
    }
    tau_es = 2.*Ne*SIGMA_THOMSON*h*max; // eq. 40 of MANMOTO+1997
    s = tau_es + tau_es*tau_es;

    // USE GSL TO FIND NORMALIZED INCOMPLETE GAMMA FUNCTION gsl_sf_gamma_inc_Q_e
    Qgamma1 = gsl_sf_gamma_inc_Q_e(jm + 1., a_temp*s);
    Qgamma2 = gsl_sf_gamma_inc_Q_e(jm + 1., *s);

    eta = exp(s*(a_temp - 1.))*(1. - Qgamma1) + eta_max*Qgamma2;

    return(eta);
}

double get_kappanu(double jnu, double Bnu)
{
    return(jnu/(4.M_PI*Bnu));
}

double get_taunu(double kappanu, double h)
{
    double fac = sqrt(M_PI)/2.;

    return(fac*kappanu*h);
}

double get_Fnu(double Thetae, double Ne, double nu, double h, double jnu, double Bnu, double taunu)
{
    /*
        Returns energy flux Fnu, eq. 26 of MANMOTO+1997
    */

    double Fnu = (2.*M_PI/3.)*Bnu*(1. - exp(-2.*sqrt(3)*taunu)); // eq. 26 of MANMOTO+1997

    return(Fnu);
}

double func(double Thetae, double Ne, double nu, double h, double kappanu, double taunu, double Fnu)
{
    /*
        Integrand of eq. 41 of MANMOTO+1997
    */

    double energy_eta(double Thetae, double Ne, double nu, double h, double kappanu, double taunu);

    return(2.*Fnu*energy_eta(Thetae, Ne, nu, h, kappanu, taunu));
}


double get_Bnu(double Thetae, double nu)
{
	double x;

	x = HPL*nu/(ME*CL*CL*Thetae);
	if (x < 1.e-3) // Taylor expand
        return ((2.*HPL/(CL*CL))/
			(x/24.*(24. + x*(12. + x*(4. + x)))));
	else
		return ((2.*HPL/(CL*CL))/(exp(x) - 1.));
}


double coulomb_cooling(double Te, double Ti)
{
    double xe, xi, sumTh, productTh, ratio, Thetae, Thetai, lnL;
    double K2e, K2i, K1, K0;

    Thetae = KBOL*Te/(ME*CL*CL);
    Thetai = KBOL*Ti/(MP*CL*CL);
    xe = 1./Thetae;
    xi = 1./Thetai;
    sumTh = Thetae + Thetai;
    productTh = Thetae*Thetai;
    ratio = sumTh/productTh;
    lnL = 20.; // Coulomb logarithm

    K2e = gsl_sf_bessel_Kn(2, xe);
    K2i = gsl_sf_bessel_Kn(2, xi);
    K1 = gsl_sf_bessel_Kn(1, ratio);
    K0 = gsl_sf_bessel_Kn(0, ratio);

    if (xi > 7.e2) {
        // rate of Coulomb collision cooling
        Qie = 37.85/MBH*rho*rho*(Te - Ti)*(sqrt(2./M_PI) + sqrt(sum))/pow((sum), 1.5); // where does this come from?
    } else {
        //Qie = 2.366*20.*0.75/MBH*rho*rho*(y[3] - y[2])/
        //        mbsl4(2,1./Thetae)/mbsl4(2,1./Thetai)*
        //        ((2.*(pow(Thetae + Thetai), 2.) + 1.)/
        //        (Thetae + Thetai)*mbsl4(1,ratio)+2.*mbsl4(0,ratio)); // Equation 13 of Manmoto+1997

        Qie = 1.25*(3./2.)*(ME/MP)*Ne*Ne*SIGMA_THOMSON*CL;
        Qie *= KBOL*(y[3] - y[2])/(K2e*K2i)*lnL;
        Qie *= (K1*(2.*sum*sum + 1.)/sum + 2.*K0);
    }

    return(Qie);

}
