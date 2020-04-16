#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

int main() {
    // Create low-resolution data
    const int ni = 40, nj=40;
    double *x = (double *) malloc(ni * sizeof(double));
    double *y = (double *) malloc(nj * sizeof(double));
    double *z = (double *) malloc(ni * nj * sizeof(double));
    for (int i = 0; i < ni; i++) {
        for (int j = 0; j < nj; j++) {
            x[i] = i;
            y[j] = j;
            z[j*ni + i] = sin(x[i]) * cos(y[j]); // note: column major
        }
    }

    // Prepare the interpolation library
    const gsl_interp2d_type *T = gsl_interp2d_bicubic;
    gsl_spline2d *mySpline = gsl_spline2d_alloc(T, ni, nj);
    gsl_interp_accel *xacc = gsl_interp_accel_alloc();
    gsl_interp_accel *yacc = gsl_interp_accel_alloc();
    gsl_spline2d_init(mySpline, x, y, z, ni, nj);

    // Call the library at an arbitrary set of (x,y) positions
    const int upSamp = 5;
    for (int i = 0; i < upSamp*ni; i++) {
        for (int j = 0; j < upSamp*nj; j++) {
            double xx = double(i) / upSamp;
            double yy = double(j) / upSamp;
            double zz = 0;
            if (xx >= x[0] && xx <= x[ni-1] && yy >= y[0] && yy <= y[nj-1])
                zz = gsl_spline2d_eval(mySpline, xx, yy, xacc, yacc);
            printf("%f %f %f\n", xx, yy, zz);
        }
    }

    // Clean up and exit
    gsl_spline2d_free(mySpline);
    gsl_interp_accel_free(xacc);
    gsl_interp_accel_free(yacc);
    free(x);
    free(y);
    free(z);

    return 0;
}

