#include <stdio.h>
#include "first_method.h"
int first_method(int n, double a, double b, double *X, double *F, double *koeff, double derivative1, double derivativeN)
{
    double c1, c2, c3, c4, d1, d2, f_help;
    double h = (b - a) / (n - 1);
    d1 = derivative1;
    d2 = (F[2] - F[0]) / (2 * h);
    for (int i = 0; i < n; ++i) {
        f_help = (F[i + 1] - F[i]) / (X[i + 1] - X[i]);
        c1 = F[i];
        c2 = d1;
        c3 = (3 * f_help - 2 * d1 - d2) / (X[i + 1] - X[i]);
        c4 = (d1 + d2 - 2 * f_help) / (X[i + 1] - X[i]) / (X[i + 1] - X[i]);
        koeff[i * 4] = c1;
        koeff[1 + i * 4] = c2;
        koeff[2 + i * 4] = c3;
        koeff[3 + i * 4] = c4;
        d1 = d2;
        if (i < n - 2) {
            d2 = (F[i + 2] - F[i]) / (2 * h);
        }
        else d2 = derivativeN;
    }
    return 0;
}
