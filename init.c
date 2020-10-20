#include "libs.h"

void read_params(char *fname)
{
    /*
        Reads parameters from parameter file and stores in global variables.
    */
    
    FILE *fp;

    fp = fopen(fname, "r");
    if (fp == NULL) {
		fprintf(stderr, "ERROR: Parameter file not found.\n");
		exit(1);
	} else {
        fprintf(stderr, "Parameter file found: '%s'.\n", fname);
		fprintf(stderr, "Successfully opened parameter file.\n");
	}
    fprintf(stderr, "Reading parameter file...\n");

    // this is readParam in dyn.pl
    fscanf(param_file, "%lf\n", &gam); // adiabatic index
    fscanf(param_file, "%lf", &MBH); // BH mass
    fscanf(param_file, "%lf", &beta); // plasma beta
    fscanf(param_file, "%lf", &alfa); // alpha viscosity
    fscanf(param_file, "%lf", &delta); // fraction of turbulent dissipation that directly heats electrons
    fscanf(param_file, "%lf", &epsilon); // accretion efficiency
    fscanf(param_file, "%lf", &M0_dot); // Mdot_out in Eddington units
    fscanf(param_file, "%lf", &rout); // Torus outside radius (in R_S)
    fscanf(param_file, "%lf", &pp0); // strength of wind
    fscanf(param_file, "%lf", &sl0i);
    fscanf(param_file, "%lf", &sl0f);
    fscanf(param_file, "%d", &nmodels);
    fscanf(param_file, "%lf", &ti); // ion temperature
    fscanf(param_file, "%lf", &te); // electron temperature
    fscanf(param_file, "%lf", &vcs); // Mach speed (v_R/v_cs)
    fscanf(param_file, "%lf", &diag);

    fprintf(stderr, "Finished reading parameter file.\n");
    fclose(fp);
}

void set_BCs()
{
    /*
        Self-consistently determine boundary conditions at the outer boundary
        of the accretion flow.
        WARNING: the formulas below are valid in the range R=100-10000 R_S
    */

    double T_vir;

    if ((rout < 100.) || (rout > 10000.)) {
        fprintf(stderr, "ERROR: ADAF outer radius must be 100 <= Rout <= 10000 R_S\n");
        exit(1);
    }

    fprintf(stderr, "Setting boundary conditions...\n");

    // Virial temperature at outer boundary
    //tvir = 0.5444091492e13 * (gam - 1.)/rout;
    T_vir = 3.6e12/rout;

    T_vir = (gam - 1.)*GNEWT*MBH*MP/(KBOL*rout);

    // Ion temperature
    T_i = T_i*T_vir;

    // Electron temperature
    T_e = T_e*T_vir;

    // Mach speed
    //vcs="0.5d0";
    //vcs=0.503 - 3.03e-5*rout;
    //vcs=vcs;

    M_Edd_dot = (4.*M_PI*GNEWT*MP*CL/SIGMA_THOMSON)*(LSUN/MSUN)*MBH/(epsilon*CL*CL);
    M0dot *= M_Edd_dot;
    rout *= 2.;
    rtran = rout*1.001; // transition radius
}
