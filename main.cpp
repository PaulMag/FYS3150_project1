#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

double F(double x) {
    return 100 * exp(-10*x);
}

double U(double x) {
    return 1 - (1 - exp(-10)) * x - exp(-10*x);
}

int main()
{
    int n = 10000;
    int a = -1;
    int b =  2;
    int c = -1;

    double *v = new double[n+2];
    double *x = new double[n+2];
    double *f = new double[n+2];

    double start = 0;
    double end = 1;
    double h = (end - start) / (n+1);
    x[0]   = start;
    x[n+1] = end;

    for (int i=0; i < n+2; i++) {
        x[i] = start + i*h;
        v[i] = U(x[i]);
        f[i] = F(x[i]);
    }

    double * bBar = new double[n+2];
    for (int i=1; i < n+1; i++) {
        bBar[i] = a*v[i-1] + b*v[i] + c*v[i+1];
    }

    ofstream x_data;
    ofstream analytical_data;
    ofstream numerical_data;
    x_data.open("x_data.txt");
    analytical_data.open("analytical_data.txt");
    numerical_data.open("numerical_data.txt");

    for (int i=0; i < n+2; i++) {
        x_data << x[i] << endl;
        analytical_data << f[i] << endl;
        numerical_data << bBar[i] / (h*h) << endl;
    }

    analytical_data.close();
    numerical_data.close();
    delete [] v;
    delete [] x;
    delete [] f;
    delete [] bBar;
    return 0;
}
