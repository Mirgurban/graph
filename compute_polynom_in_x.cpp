#include "compute_polynom_in_x.h"
double compute_pol(int n, double a, double b, double x, double *koeff)
{
    int buf = 0;
    double h;
    double x_buf;
    h = (b - a) / (double)(n - 1);
    for (int i = -1; i < n; i++)
    {
        if (x < a + (i + 1) * h && x > a + i * h)
        {
            buf = i * 4;
        }
    }
    x_buf = a + (buf / 4) * h;
    return koeff[buf] + koeff[buf + 1] * (x - x_buf) + koeff[buf + 2] * (x - x_buf) * (x - x_buf) + koeff[buf + 3] * (x - x_buf) * (x - x_buf) * (x - x_buf);
}
