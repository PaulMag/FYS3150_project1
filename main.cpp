#include <iostream>
#include <cmath>
#include <fstream>
#include <armadillo>

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

int main(int argc, char* argv[])
{
    int n = 6;
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

    v[0]   = v_0;
    v[n-1] = v_n_1;

    double h = 1. / (n-1);

    // Fill in arrays:
    for (int i=0; i < n; i++) {
        x[i] = i*h;
        u[i] = U(x[i]);
        fBar[i] = F(x[i]) * (h*h);
    }

    // The wrong algorithm:
    /*for (int i=1; i < n+1; i++) {
        bHat[i] = a*v[i-1] + b*v[i] + c*v[i+1];
    }*/

    // The right algorithm:
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

    // Compare with armadillo's solver:
    mat A(n,n);
    vec v_a(n);
    vec f_a(n);

    A.zeros();
    A(0, 0) = b;
    f_a(0)  = fBar[0];

    for (int i=1; i<n; i++) {
        A(i-1, i) = c;
        A(i,i)    = b;
        A(i, i-1) = a;
        f_a(i) = fBar[i];
    }

    v_a = solve(A, f_a);


    // (c)

    double *eps = new double[n];
    double maxError = 0;
    for (int i=1; i < n-1; i++) { // avoid end because zero divison
        eps[i] = difference(v[i], u[i]);
        if (eps[i] > maxError) {
            maxError = eps[i];
        }
    }
    cout << maxError << endl;


    // End:
    delete [] v;
    delete [] u;
    delete [] x;
    delete [] fHat;
    delete [] fBar;
    delete [] bHat;

    // (e)

    return 0;
}
