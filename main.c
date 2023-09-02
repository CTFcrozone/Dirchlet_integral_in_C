#include <stdio.h>
#include <complex.h>
#include <math.h>

double dirichlet_integral() {
    double a = 0.0;
    double b = 2.0 * M_PI;
    int N = 1000000;

    double complex result = 0.0;

    for (int n = -N; n <= N; n++) {
        double theta = (n * (b - a) / (2.0 * N)) + ((b + a) / 2.0);
        double complex z = cos(theta) + I * sin(theta);
        double complex residue = cexp(I * z) / (z * z);
        result += residue;
    }

    result *= (a-b) / (2.0 * N);
    return creal(result);
}

int main() {
    double integral_value = dirichlet_integral();
    printf("Dirichlet Integral Value: %.5f\n", integral_value);
    return 0;
}
