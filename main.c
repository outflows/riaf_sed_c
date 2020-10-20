#include "libs.h"

int main() {
    double sl0, d_sl0, lastok;
    int iterat, brackets, ops;

    read_params();
    set_BCs();

    iterat = 1;
    brackets = 0;

    // Calculates the increment in sl0 given the desired number of models
    sl0 = sl0i;
    d_sl0 = (sl0f - sl0i)/nmodels;

    // Main loop
    while (sl0 <= sl0f + d_sl0) {
        fprintf(stderr, "Eigenvalue = %lf Iteration = %d Brackets = %d\n\n", sl0, iterat, brackets);
        // Creates header of each log file
        // The output file keeps being rewritten until the final run. -- WHY??
        write_header();

        // Calls adaf Fortran code and computes dynamical solution
        dynamics(sl0);

        // Diagnose if the computed solution is physical or not
        diagnose();

        // Decides if the solution is OK or not, and what should be done next
        if ((increase==1) && (discont==0) && (weirdam==0) && (sonic==0) &&
            (nan==0) && (no_output==0)) {
            fprintf(stderr, "Solution is good. Iterations = %d Brackets = %d\n\n",
                    iterat, brackets);
            ops = 0; // signal that solution was found, used after end of loop
            break;

        } else if ((increase==0) && (sonic==1)) {
            // stores the last eigenvalue computed with no stops or jumps
            lastok = sl0;
            sl0 = sl0 + d_sl0;
            printf("Subsonic solution.\n");
            ops = 1;
            iterat++;

        /*
        } else if (no_output==1) {
            fprintf(stderr, "No global solution. Change the OBCs!\n");
            ops = 1;
            break;
        */

        } else if (no_output==1 ||
                 (increase==1 && sonic==1) ||
                 (increase==1 && discont==1) ||
                 (weirdam==1) ||
                 (discont==1 && sonic==0)) {
            //  || (nan==1) || (increase==0 && $sonic==0) {}
            if (no_output==1) {
                fprintf(stderr, "No global solution. Verify the OBCs!\n\n");
            }
            if (iterat==1) {
                fprintf(stderr, "Bad eigenvalue at 1st iteration! Decrease the lower limit.\n\n";)
                ops = 1;
                break;
            }

            // stores the "bad" eigenvalue which caused a jump or hang
            bad = sl0;
            // BRACKET the solution between the last OK eigenvalue and the current BAD one
            sl0f = sl0;
            d_sl0 = (bad - lastok)/nmodels;
            sl0 = lastok + d_sl0;
            fprintf(stderr, "Bracketing condition found.\n\n");
            iterat++;
            brackets++;

        } else {
            fprintf(stderr, "Something weird happened! Check conditions!\n\n");
            ops = 1;
            break;
        }
        // Prints results of diagnostics
        fprintf(stderr,"Diagnostics: \n");
        fprintf(stderr,"Increasing: %d\n", increase);
        fprintf(stderr,"Jumps: %d\n", discont);
        fprintf(stderr, "R sonic: %d\n", sonic);
        fprintf(stderr, "Weird AM: %d\n", weirdam);
        fprintf(stderr, "NaN: %d\n", nan);
        fprintf(stderr, "Failed: %d\n", failed);
        fprintf(stderr, "No output: %d\n", no_output);
        fprintf(stderr, "Mach max: %lf\n", largest);
        fprintf(stderr, "R max: %lf\n\n", largestR);
    }

    if (ops == 1) {
        //Print diagnostics for the bad solution
        fprintf(stderr,"Bad solution diagnostics: \n");
        fprintf(stderr,"Increasing: %d\n", increase);
        fprintf(stderr,"Jumps: %d\n", discont);
        fprintf(stderr, "R sonic: %d\n", sonic);
        fprintf(stderr, "Weird AM: %d\n", weirdam);
        fprintf(stderr, "NaN: %d\n", nan);
        fprintf(stderr, "Failed: %d\n", failed);
        fprintf(stderr, "No output: %d\n", no_output);
        fprintf(stderr, "Mach max: %lf\n", largest);
        fprintf(stderr, "R max: %lf\n", largestR);
        fprintf(stderr, "No solution found within the eigenvalue interval. Try changing the OBCs, or range/step of eigenvalues.\n");
    } else {
        //Print diagnostics for the good solution
        fprintf(stderr,"Good solution diagnostics: \n");
        fprintf(stderr,"Increasing: %d\n", increase);
        fprintf(stderr,"Jumps: %d\n", discont);
        fprintf(stderr, "R sonic: %d\n", sonic);
        fprintf(stderr, "Weird AM: %d\n", weirdam);
        fprintf(stderr, "NaN: %d\n", nan);
        fprintf(stderr, "Failed: %d\n", failed);
        fprintf(stderr, "No output: %d\n", no_output);
        fprintf(stderr, "Mach max: %lf\n", largest);
        fprintf(stderr, "R max: %lf\n", largestR);
        fprintf(stderr, "Number of lines in output file: %d\n\n", linesout);
    }

    // Benchmarking execution time
    //bench1 = new Benchmark;
    //dbench = timediff($bench1, $bench0);
    //fprintf(stderr, "Total running time: \n");

    return 0;
}
