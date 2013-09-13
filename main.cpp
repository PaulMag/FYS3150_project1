#include <iostream>
#include <cmath>
#include <fstream>
#include <armadillo>
#include "time.h"

using namespace std;
using namespace arma;

double F(double x) {
    return 100 * exp(-10*x);
}

double U(double x) {
    return 1 - (1 - exp(-10)) * x - exp(-10*x);
}

double difference(double v, double u) {
    return log10(abs( (v - u) / u ));
}

mat matProduct(mat b, mat c) {

    int n = b.n_cols; // assumes both b and c are quare and same size
    mat a = zeros<mat>(n, n);

    for (int i=0 ; i<n ; i++) {
        cout << i << " of " << n << endl;
        for (int j=0 ; j<n ; j++) {
            for (int k=0 ; k<n ; k++) {
                a(i,j) += b(i,k) * c(k,j);
            }
        }
    }
    return a;
}

// int main(int argc, char* argv[]) // command line arguments
int main() {

    /******************** (b) ********************/

    int n = 1e1; // CHANGE THIS TO CHANGE STEP LENGTH
    double a = -1;
    double b =  2;
    double c = -1;

    double v_0 = 0;
    double v_n_1 = 0;

    double *v = new double[n]; // numeric
    double *u = new double[n]; // analytic
    double *x = new double[n];
    double *fBar = new double[n]; // fBar instead of b in Av = b
    double *fHat = new double[n];
    double *bHat = new double[n];

    clock_t *start = new clock_t[10];
    clock_t *finish = new clock_t[10];
    int measurements = 0;

    v[0]   = v_0;
    v[n-1] = v_n_1;

    double h = 1. / (n-1);

    // Fill in arrays:
    for (int i=0; i < n; i++) {
        x[i] = i*h;
        u[i] = U(x[i]);
        fBar[i] = F(x[i]) * (h*h);
    }

    // The algorithm:

    // *****TiC*****
    start[measurements] = clock();
        bHat[0] = b;
        fHat[0] = fBar[0];
        double factor;

        for (int i=1; i<n; i++) {
            factor = a / bHat[i-1]; // avoids doing this twice
            bHat[i] = b - c * factor;
            fHat[i] = fBar[i] - fHat[i-1] * factor;
        }

        v[n-1] = (fBar[n-1]) / bHat[n-1];
        for (int i=n-2; i>=0; i--) {
            v[i] = (fHat[i] - v[i+1] * c) / bHat[i];
        }
    finish[measurements] = clock();
    measurements++;
    // *****TOC*****

    // Writing data to files:
    ofstream info;
    ofstream x_data;
    ofstream analytical_data;
    ofstream numerical_data;

    info.open("info.txt");
    x_data.open("x_data.txt");
    analytical_data.open("analytical_data.txt");
    numerical_data.open("numerical_data.txt");

    info << n << endl;
    for (int i=0; i < n; i++) {
        x_data << x[i] << endl;
        analytical_data << u[i] << endl;
        numerical_data << v[i] << endl;
    }

    info.close();
    analytical_data.close();
    numerical_data.close();


    /******************** (c) ********************/

    double *eps = new double[n];
    double maxError = difference(v[1], u[1]);

    for (int i=2; i < n-1; i++) { // avoid end because zero divison
        eps[i] = difference(v[i], u[i]);
        if (eps[i] > maxError) {
            maxError = eps[i];
        }
    }
    cout << "Largest error: " << maxError << endl;


    /******************** (d) ********************/
    // Compare with armadillo's solver:

    mat A_a(n,n);
    vec v_a(n);
    vec f_a(n);

    A_a.zeros();
    A_a(0, 0) = b;
    f_a(0)  = fBar[0];

    for (int i=1; i<n; i++) {
        A_a(i-1, i) = c;
        A_a(i,i)    = b;
        A_a(i, i-1) = a;
        f_a(i) = fBar[i];
    }

    // *****TIC*****
    start[measurements] = clock();
        v_a = solve(A_a, f_a);
    finish[measurements] = clock();
    measurements++;
    // *****TOC*****

    // End of main project.
    delete [] v;
    delete [] u;
    delete [] x;
    delete [] fHat;
    delete [] fBar;
    delete [] bHat;


    /******************** (e) ********************/

    n = 3e2; // CHANGE THIS TO CHANGE MATRIX SIZES
    mat B = randu<mat>(n, n);
    mat C = randu<mat>(n, n);
    mat A(n, n);

    // *****TIC*****
    start[measurements] = clock();
        A = matProduct(B, C); // my slow algorithm
    finish[measurements] = clock();
    measurements++;
    // *****TOC*****

    // *****TIC*****
    start[measurements] = clock();
        /* Armadillo's (hopefully) fast algorithm. I am not 100% certain
         * that this uses BLAS, but as far as I understood that is the default.
         */
        A = B * C;
    finish[measurements] = clock();
    measurements++;
    // *****TOC*****

    // Write algorithm times into file:
    ofstream timer_data;
    timer_data.open("timer_data.txt");

    for (int i=0; i<measurements; i++) {
        timer_data << ((finish[i] - start[i]) * 1.0 / CLOCKS_PER_SEC) << endl;
    }
    timer_data.close();

    delete [] start;
    delete [] finish;

    return 0;
}
