



rhorho = rho*6.1749e5/(MBH*MBH);
hh = h/6.77e-12*MBH;

// the first, the bremsstrahlung radiation
if (Thetae < 1.) {
    F_ei = 4.*sqrt(2.*Thetae/M_PI)/M_PI*(1. + 1.781*pow(Thetae,1.34)) + // Eq. 29 of Manmoto+1997
    1.73*pow(Thetae, 1.5)*(1. + 1.1*Thetae + Thetae*Thetae - 1.25*pow(Thetae,2.5)); // Eq. 31 of Manmoto+1997
} else {
    F_ei = 9.*Thetae/2./M_PI*(log(1.123*Thetae + 0.48) + 1.5) + // Eq. 30 of Manmoto+1997
        2.3*Thetae*(log(1.123*Thetae) + 1.28); // Eq. 32 of Manmoto+1997
}

qbr = 5.3e25*rhorho*rhorho*F_ei/MUE/MUE;
