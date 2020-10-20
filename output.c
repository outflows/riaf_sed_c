void write_header(FILE *fname)
{
    FILE *fp;

    fp = fopen(fname, "w");

    fprintf(stderr, "Writing header into file '%s'...\n", fname);

    fprintf(fp, "gam,MBH,beta,alpha,delta,dotm0,rout,pp0,ti/tvir,te/tvir,vcs,sl0,tvir/1e9\n");
    fprintf(fp, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
            gam,MBH,beta,alpha,delta,dotm0,rout,pp0,ti/tvir,te/tvir,vcs,sl0,tvir/1.e9);

    fclose(fp);
}
