#define SYNCHROTRON (1)
#define BREMSSTRAHLUNG (1)
#define COMPTON (1)
#define NONTHERMAL (1)

#define ME 9.1093826e-28   // electron mass (g)
#define MP 1.67262171e-24  // proton mass (g)
#define CL 2.99792458e10   // speed of light (cm/s)
#define GNEWT 6.6742e-8    // gravitational constant (cgs)
#define KBOL 1.3806505e-16 // Boltzmann constant (erg/K)
#define SIGMA_THOMSON 6.65245873e-25 // Thomson cross-section (cm^2)
#define ALPHA_F 0.007297351 // fine structure constant
#define RE 2.8179403227e-13 // Classical electron radius (cm)
#define MSUN 1.989e33      // solar mass (g)
#define LSUN 3.827e33      // solar luminosity (erg/s)

#define MUI = 1.24 // Mean molecular weight, ions
#define MUE = 1.13 // Mean moleculaer weight, electron
#define RU 9.2e-14 // k/M_mu -- what is this?!?!

extern double Ne; // electron number density
extern double M_Edd_dot // Mdot Eddington
extern double cs; // speed of sound
extern double rtran; // transition radius
extern double M0dot;

extern double gam; // adiabatic index
extern double MBH; // BH mass (in 10^6 solar masses)
extern double beta; // plasma beta
extern double alfa; // alpha viscosity
extern double delta; // fraction of turbulent dissipation that directly heats electrons
extern double dotm0; // Mdot_out in Eddington units
extern double rout; // Torus outside radius (in R_S)
extern double pp0; // strength of wind
extern double sl0i; // initial eigenvalue
extern double sl0f; // final eigenvalue
extern int nmodels; // number of models to be computed
extern double T_i; // ion temperature
extern double T_e; // electron temperature
extern double vcs; // Mach number=v_R/v_cs (radial velocity/sound speed)
extern char diag[100]; // name of log file
