#include "libs.h"

void dynamics(double sl0)
{
    int i, k2, inform;
    double y[4];
    double cs, p0, sint, Omegar2, OmegaKr2, rho;
    double sk, Omega, OmegaK, Omegar2;
    double tau, y0, H_disc;
    double rdomega, advec, advec2;
    double dh, tao, compy, pgas, Bmag, norm, Mdout_out;
    double sumlum, eddlum;
    double r[1000], H_disc_arr[1000], sum[9000], sigma[1000]; //H_disc_arr is disc thickness
    FILE *fp1, *fp2;

    //M_Edd_dot = 3.44E-15*MBH; //M_edd_dot
    cs = sqrt(RU/beta*(T_i/MUI + T_e/MUE)); // speed of sound
    p0 = 0.;
    sint = 0.; // integral of s, eq. 1 of YUAN+2005

    y[0] = rout;
	y[1] = -cs*vcs;
	y[2] = T_e;
	y[3] = T_i;

    sk = 1./((sqrt(y[0])*(y[0] - 2.))); // ?
	Omegar2 = -alpha*y[0]*cs*cs/y[1] + sl0; // Eq. 3 of YUAN+2008, why is it called sout in the fortran code?
	Omega = Omegar2/(y[0]*y[0]); // From eq. 3 of YUAN+2008
	OmegaK = sk; // Keplerian velocity at the equator

	if (Omega/OmegaK > 1.) {
        fprintf(stderr, "Must increase vcs. Omega/OmegaK = %lf > 1\n", Omega/OmegaK);
        exit(1);
    }

    i = 2;
    k2 = 0;
    inform = i + 10;
    tau = 0.;
	y0 = rtran;

    do {
        if(k2 < 160) H_disc = -1.e1;
    	if(k2 >= 160) H_disc = -1.e0;
        if(k2 >= 500) H_disc = -1.e-1;

        H_disc = -y[0]/100.;

        mers(4, inform, H_disc, 1.e-4, y, dy, w, qn, qd, tau, sint); // AWFUL FUNCTION FULL OF GOTOS
        // this function will change inform, y, dy, qn, qd, tau, sint

        cs = sqrt(RU/beta*(y[3]/MUI + y[2]/MUE)); // speed of sound
        Omegar2 = -alpha*cs*cs*y[0]/y[1] + sl0; // Eq. 3 of YUAN+2008, why is it called ssll in fortran code?
    	Omega = Omegar2/(y[0]*y[0]); // angular velocity
    	OmegaK = 1./((sqrt(y[0])*(y[0] - 2.))); // Keplerian angular velocity
    	OmegaKr2 = OmegaK*y[0]*y[0]; // why is it called slk in fortran code? k is for keplerian, what about sl?
     	rho = -M0dot*exp(sint)*OmegaK/(4.*M_PI*y[0]*cs*y[1]); // gas density
    	tau = tau + rho*(y0 - y[0])*3.65e16/MBH; // understand this 3.65e16/MBH, what is tau?
        y0 = y[0];

        rdomega = get_rdomega(cs, Omega, y);

        advec = 0.;
        advec2 = 0.;
        get_advec(rho, p0, y, gam, k2, sint, advec, advec2);

        dh = cs/OmegaK;
        tao = 2.5*dh*rho*3.65e16/MBH; // understand this 3.65e16/MBH, also what is tao?!?! where does 2.5 come from?
        //dotmr = -4.*M_PI*dh*rho*y[0]*y[1]/M_Edd_dot;
        // now tao is in units of cgs
        //WTF?!?!?! this never appears anywhere!!!!!
        //poten = Omegar2*Omegar2/(2.*y[0]*y[0]) - 1./(y[0] - 2.);
        //Be = 0.5*y[1]*y[1] + 1.5/0.5*RU*y[2] + 2.5*RU*y[3] + Omegar2*Omegar2/2./(y[0] - 2.);
        compy = 4.*tao*(KBOL/(ME*CL*CL))*y[2]; // what is this?!?! 4 tao Thetae...
        pgas = rho*6.e5/(MBH*MBH) * KBOL/MP*(y[2] + y[3]); // where does this 6.e5/(MBH*MBH) come from?
        Bmag = sqrt((1. - beta)/beta*pgas*8.*M_PI); // magnetic field, p=pgas/beta, pmag=(1-beta)*p=B2/8pi

        k2++; // should this really be here and not in the end of the loop?
        if( (fabs(k2/10 - k2/10.) < 1.e-5) || (k2 == 1) ) {
            norm = sqrt((4. + 1.33*alpha*alpha)/(8./3.)); // what is this? where do the factors come from?
            Mdot_out = -rho*4.*M_PI*y[0]*cs*y[1]/(OmegaK*M_Edd_dot*0.05); // Mdot outflow, what is 0.05???

            // Modified the standard output
            // the two last columns are the density and electron number density in cgs units

            fprintf(stderr, "y[0] = %lf\n", y[0]);
            fprintf(stderr, "-y[1]/(cs*norm) = %lf\n", -y[1]/(cs*norm));
            fprintf(stderr, "log10(y[2]) = %lf\n", log10(y[2]));
            fprintf(stderr, "log10(y[3]) = %lf\n", log10(y[3]));
            fprintf(stderr, "advec = %lf\n", advec);
            fprintf(stderr, "advec2 = %lf\n", advec2);
            fprintf(stderr, "cs = %lf\n", cs);
            fprintf(stderr, "dh = %lf\n", dh);
            fprintf(stderr, "qn = %lf\n", qn);
            fprintf(stderr, "tau = %lf\n", tau);
            fprintf(stderr, "log10(OmegaKr2) = %lf\n", log10(OmegaKr2));
            fprintf(stderr, "log10(Omegar2) = %lf\n", log10(Omegar2));
            fprintf(stderr, "Bmag = %lf\n", Bmag);
            fprintf(stderr, "log10(rho*6.1749e5/(MBH*MBH)) = %lf\n", log10(rho*6.1749e5/(MBH*MBH)));
            fprintf(stderr, "log10(rho*6.1749e5/(MBH*MBH*1.67d-24*MUE)) = %lf\n", log10(rho*6.1749e5/(MBH*MBH*1.67d-24*MUE)));

            fp1 = fopen("x.dat", "w");
            fprintf(fp1, "%lf\t", log10(y[0]/2.));
            fprintf(fp1, "%lf\t", -y[1]);
            fprintf(fp1, "%lf\t", log10(y[2]));
            fprintf(fp1, "%lf\t", log10(y[3]));
            fprintf(fp1, "%lf\t", log10(OmegaKr2));
            fprintf(fp1, "%lf\t", log10(Omegar2));
            fprintf(fp1, "%lf\t", log10(tau));
            fprintf(fp1, "%lf\t", Bmag);
            fprintf(fp1, "%lf\t", -8.e8/(y[0]/y[1]/0.203*MBH)/Bmag/Bmag); // electron cooling Lorentz factor
            fprintf(fp1, "%lf\t", compy);
            fprintf(fp1, "%lf\t", advec);
            fprintf(fp1, "%lf\n", qn*y[0]*y[0]*dh*1.e18));
            fclose(fp1);

            fp2 = fopen("hottem.dat", "w");
            fprintf(fp2, "%lf\t%lf\t%lf\t%lf\n", y[0], dh, qn, tau);
            fclose(fp2);

        }
        if( (y[0] < 1.6) && (mmm == 0) ) {
            y[3] = (3./16.)*y[1]*y[1]/RU;
            y[2] = y[3];
            y[1] = y[1]/4.;
            inform = 12;
            mmm = 1;
            fprintf(stderr, "something went wrong, line 277\n");
        }

        // Luminosity emitted by the hot disc
        r[k2] = y[0];
        H_disc_arr[k2] = cs/OmegaK;
        rho = -dotm*OmegaK/(4.*M_PI*y[0]*cs*y[1]); // gas density
        sum[k2] = qn*8.87e-27*MBH*MBH*MBH;

    } while (y[0] >= 2.2e0);

    sumlum = 0.;
    for (i = 1; i < k2-2; i++) {
        sigma[i] = M_PI*( pow((r[i-1] + r[i]), 2.)/4. - (r[i] + r[i+1])**2./4.);
        sumlum = sumlum + sigma[i]*2.*H_disc_arr[i]*sum[i]]/2.*2.;
        // in the above formula, a factor '/2.d0' should be present. This is because
        //  what we observe is only half of the total emitted radiation: no!!!!

        eddlum = 1.3e38*1.e6*MBH*5.59e-61/MBH/2.03e-1*MBH;

        eddlum = (4.*M_PI*GNEWT*MP*CL/SIGMA_THOMSON)*MBH/MSUN; // erg/s
    }

    fprintf(stderr, "Total luminosity = %lf% lf\n", sumlum/eddlum, sumlum/(5.59e-61/MBH/2.03e-1*MBH);
}



/***********************************************************************/

void dery(int n, double y[4], double dy[4], double qn, double qd, double tau, double sint)
{
    double OmegaK, dlnok, cs, Omegar2;
    double Omega, rho, p, h;
    double Qie, ate;
    double jnu, Bnu, kappanu, taunu, Fnu, eta, result, Qrad;
    double Qssd, Qtot, qneg;
    double Thetae, Thetai, xe, xi, K1, K2, K3;
    double ww, w1, w2, w3;
    double a1, a2, a3;
    double b1, b2, b3, b4;
    double c1, c2, c3, c4;
    double d1, d2, d3, d4;
    double dd, dd1, dd2, dd3;

    if ( (y[2] < 0) || (y[3] < 0) ) {
        return;
    }

    OmegaK = 1./ ((sqrt(y[0]) * (y[0] - 2.)));
    dlnok = -0.5/y[0] - 1./(y[0] - 2.);
    cs = sqrt(RU/beta*(y[3]/MUI + y[2]/MUE)); // speed of sound
    Omegar2 = -alpha*cs*cs*y[0]/y[1] + sl0; // Eq. 3 of YUAN+2008

    if (Omegar2 < 0.) {
        fprintf(stderr, "ERROR: Omegar2 <0\n");
        exit(1);
    }

    Omega = Omegar2/(y[0]*y[0]);
    rho = -M0dot*exp(sint)*OmegaK/(4.*M_PI*y[0]*cs*y[1]); // gas density
    p = rho*cs*cs; // total pressure
    h = cs/(OmegaK*(6.77e-12*MBH)); // what is this factor?? transforming to cgs? how?

    ww = -alpha/y[1]*cs*cs;
    w1 = -ww/y[1];
    w2 = ww/(y[2]/MUE + y[3]/MUI);
    w3 = ww/y[0] - 2.*Omega;

    a1 = -1./y[1];
    a2 = -1./2./(y[2]/MUE + y[3]/MUI);
    //a3=dlnok-1./y[0];
    a3 = dlnok + (p0 - 1.)/y[0];

    // note that in the following the viscous dissipation
    // heat electron is considered
    b1 = rho*y[1]*RU*a1*y[3]/MUI - alpha*p*w1*(1. - delta);
    b2 = -alpha*p*w2*(1. - delta)/MUE + rho*y[1]*RU*a2*y[3] / MUI / MUE;
    b3 = b2*MUE/MUI - 1./(gam - 1.)*rho*y[1]*RU/MUI;
    Qie = coulomb_cooling(y[2], y[3]);
    b4 = Qie + alpha*p*w3*(1. - delta) - rho*y[1]*RU*a3*y[3]/MUI;

    Thetae = KBOL*Te/(ME*CL*CL);
    Thetai = KBOL*Ti/(MP*CL*CL);
    xe = 1./Thetae;
    xi = 1./Thetai;
    K1 = gsl_sf_bessel_Kn(1, xe);
    K2 = gsl_sf_bessel_Kn(2, xe);
    K3 = gsl_sf_bessel_Kn(3, xe);

    c1 = a1*y[2] - alpha*delta*p*w1/(rho*y[1]*RU/MUE);
    ate = xe*((3.*K3 + K1)/(4.*K2) - 1.);
    c2 = a2*y[2]/MUE - ate - alpha*delta*p*w2/(MUE*rho*y[1]*RU/MUE);
    c3 = a2*y[2]/MUI - alpha*delta*p*w2/(MUI*rho*y[1]*RU/MUE);

    jnu = get_jnu(Thetae, Ne, nu);
    Bnu = get_Bnu(Thetae, nu);
    kappanu = get_kappanu(jnu, Bnu);
    taunu = get_taunu(kappanu, h);
    Fnu = get_Fnu(Thetae, Ne, nu, h, jnu, Bnu, taunu);
    eta = get_eta(Thetae, Ne, nu, h, kappanu, taunu);
    //func = func(Thetae, Ne, nu, h, kappanu, taunu, Fnu);

    gsl_integration_workspace *ww = gsl_integration_romberg_alloc (1000);
    gsl_function F;
    struct my_params = {Thetae; Ne; nu; h; kappanu; taunu; Fnu};
    F.function = &func;
    F.params = &my_params;

    gsl_integration_romberg(F, 1.e9, 6.e10*y[2], 0, 1e-7, result, 1000, ww);

    Qrad = result;

    gsl_integration_romberg_free(ww);

    // Take into account Compton cooling rate Qssd due to the soft photons
    // from the outside SSD. First, convert radii to CGS
    //rr = y[0]/6.77e-12*MBH;
    //sup = 10.*rtran/6.77e-12*MBH;
    //slow = rtran/6.77e-12*MBH;
    //simp(slow,sup,10,Qssd,hh,rr,tau,rhorho,Thetae);
    Qssd = 0.;

    Qtot = Qrad + Qssd;

    // transfer qrad to the dissipation rate per volum
    qneg = Qtot/(2.*H_disc); // factor of 2 to account for both sides of disc
    qn = qneg;

    // Convert from cgs units to the c=G=M=1 units
    qneg = qn*8.87e-27*MBH*MBH*MBH;
    c4 = (qneg - Qie)/(rho*y[1]*RU/MUE) - a3*y[2] + alpha*delta*w3*p/(rho*y[1]*RU/MUE);

    q1 = cs*cs;
    d1 = y[1] + q1*a1;
    d2 = q1*a2/MUE + 1.5*RU/(beta*MUE);
    d3 = d2*(MUE/MUI);
    d4 = (Omega*Omega - OmegaK*OmegaK)*y[0] - q1*a3;

    dd  = c1*d2*b3 + d1*b2*c3 + c2*d3*b1 - b1*d2*c3 - c2*d1*b3 - b2*d3*c1;
    dd1 = c4*d2*b3 + d4*b2*c3 + c2*d3*b4 - b4*d2*c3 - c2*d4*b3 - b2*d3*c4;
    dd2 = c1*d4*b3 + d1*b4*c3 + c4*d3*b1 - b1*d4*c3 - d1*c4*b3 - b4*d3*c1;
    dd3 = c1*d2*b4 + d1*b2*c4 + c2*d4*b1 - b1*d2*c4 - d1*c2*b4 - b2*d4*c1;

    dy[0] = 1.;
    dy[1] = dd1/dd;
    dy[2] = dd2/dd;
    dy[3] = dd3/dd;
}



/***********************************************************************/

subroutine mers(n,inform,h,eps,y,dy,w,qn,qd,tau,sint)
    //dimension y(n),dy(n)

    double w[4][1];

    double a[4] = {0.3333333333,0.16666666667,0.375,2.0};
    double b[4] = {1.,0.444444444444,0.9375,8.};
    double c[4] = {0.,0.25,-0.2,-1.125};

    l = 0;
    k = w[2][0];

    if (inform > 10) { // goto 14
        inform -= 10;
        k = 1;
        l += 1; // goto 10...
        dery(n,y,dy,qn,qd,tau,sint);
        if (sig == -1.) { //goto 17
            inform = 2;
            return;
        }
        for (j = 1; j < n; j++) {
            w[0][j] = dy[j];
            w[1][j] = 0.0;
        }
        if (l != k) { // goto 1
            for (j = 0; j < n; j++) {
                w[2][j] = y[j];
                w[3][j] = dy[j];
            }
            for (i = 0; i < 4; i++) {
                for (j = 0; j < n; j++) {
                    y[j] = y[j] + ht*a[i]*(dy[j]-w[1][j]);
                    w[1][j] = b[i]*(dy[j]+c[i]*w[0][j]);
                    if(i == 3) w[0][j] = dy[j] - 0.22222222220*w[0][j];
                }
                dery(n,y,dy,qn,qd,tau,sint);
            }
            if (sig == -1.) { //goto 17
                inform = 2;
                return;
            }
        }
    }

    if (inform >= 3) { // goto 3
        for (i = 0; i < 4; i++) {
            for (j = 0; j < n; j++) {
                y[j] = y[j] + ht*a[i]*(dy[j]-w[1][j]);
                w[1][j] = b[i]*(dy[j]+c[i]*w[0][j]);
                if(i == 3) w[0][j] = dy[j] - 0.22222222220*w[0][j];
            }
            dery(n,y,dy,qn,qd,tau,sint);
        }
        if (sig == -1.) { //goto 17
            inform = 2;
            return;
        }
    }

    epsl = 0.0781250*eps;
    for (j = 0; j < n; j++) {
        w[2][j] = y[j];
        w[3][j] = dy[j];
    }
    for (i = 0; i < 4; i++) {
        for (j = 0; j < n; j++) {
            y[j] = y[j] + ht*a[i]*(dy[j]-w[1][j]);
            w[1][j] = b[i]*(dy[j]+c[i]*w[0][j]);
            if(i == 3) w[0][j] = dy[j] - 0.22222222220*w[0][j];
        }
        dery(n,y,dy,qn,qd,tau,sint);
    }
    if (sig == -1.) { //goto 17
        inform = 2;
        return;
    }

    i = 1;
    for (j = 0; j < n; j++) {
        r = 0.16666666670*ht*(dy[j] - w[1][j]);
        y[j] = y[j] + r;
    }



        l=0
        k=w(3,1)
        if(inform.ge.10) goto 14
        ht=h/w(3,1)
        if(inform.ge.3) goto 3
        epsl=0.0781250*eps
 1      do 2 j=1,n
        w(3,j)=y(j)
 2      w(4,j)=dy(j)
  3     do 5 i=1,4
        do 4 j=1,n
        y(j)=y(j)+ht*a(i)*(dy(j)-w(2,j))
        w(2,j)=b(i)*(dy(j)+c(i)*w(1,j))
 4      if(i.eq.3) w(1,j)=dy(j)-0.22222222220*w(1,j)
 5      call dery(n,y,dy,qn,qd,tau,sint)
	if(sig.eq.-1.) goto 17
        i=1
        do 9 j=1,n
        r=0.16666666670*ht*(dy(j)-w(2,j))
        y(j)=y(j)+r
        goto (7,6,9),inform
  6     if(y(j).ne.0.0) r=r/y(j)
 7      if(abs(r*0.2).gt.eps) if(y(1)-(0.125*ht+y(1))) 12,8,12
 8      if(abs(r).gt.epsl) i=3
 9      continue
        if(i.ne.l-l/2*2) goto 15
        k=k/2
        l=1+l/2
        ht=ht*2
 10     call dery(n,y,dy,qn,qd,tau,sint)
        if(sig.eq.-1.) goto 17
        do 11 j=1,n
        w(1,j)=dy(j)
 11     w(2,j)=0.0
        if(l.ne.k) goto 1
        w(3,1)=k
17	inform=2
         return
 12     ht=0.5*ht


	cshu=cshu+1
        if(cshu.gt.150.) then
         inform=3
         print*,'###########################'
        endif

        l=l+l
        k=k+k
        do 13 j=1,n
        y(j)=w(3,j)
        dy(j)=w(4,j)
        w(1,j)=dy(j)
 13     w(2,j)=0.
        goto 3
 14     inform=inform-10
        k=1
 15     l=l+1
        goto 10
        end
