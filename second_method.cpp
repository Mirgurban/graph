#include <stdio.h>
#include "second_method.h"

int second_method(int n, double *X, double *F, double *koeff, double derivative1, double derivativeN)
{
    double *d, *matrix, *column_matrix;
    double c1, c2, c3, c4, f_help, f_help1, f_help2, matrix_koeff;
    d = (double *)calloc(n, sizeof(double));
    column_matrix = (double *)calloc(n, sizeof(double));
    matrix = (double *)calloc(n * n, sizeof(double));

    // решение системы уравнений
    matrix[0] = 1;
    column_matrix[0] = derivative1;
    for (int j = 1; j < n; ++j) {
        matrix[j] = 0;
    }
    matrix[n * n - 1] = 1;
    column_matrix[n - 1] = derivativeN;
    for (int j = 0; j < n - 1; ++j) {
        matrix[j + n * (n - 1)] = 0;
    }

    for (int i = 1; i < n - 1; ++i) {
        for (int j = 0; j < n; ++j) {
            if (j == i - 1) {
                matrix[j + i * n] = X[i + 1] - X[i];
            } else if (j == i) {
                matrix[j + i * n] = 2 * (X[i + 1] - X[i - 1]);
            } else if (j == i + 1) {
                matrix[j + i * n] = X[i] - X[i - 1];
            }
            else matrix[j + i * n] = 0;
        }
    }

    for (int i = 1; i < n - 1; ++i) {
        f_help1 = (F[i] - F[i - 1]) * (X[i + 1] - X[i]) / (X[i] - X[i - 1]);
        f_help2 = (F[i + 1] - F[i]) * (X[i] - X[i - 1]) / (X[i + 1] - X[i]);
        column_matrix[i] = 3 * f_help1 + 3 * f_help2;
    }
    // алгоритм Гаусса
    for (int k = 0; k < n; ++k) {
        for (int j = k + 1; j < n; ++j) {
            matrix_koeff = matrix[j + k * n] / matrix[k + k * n];
            for (int i = k; i < n; ++i) {
                matrix[j + n * i] = matrix[j + n * i] - matrix_koeff * matrix[k + n * i];
            }
            column_matrix[j] = column_matrix[j] - matrix_koeff * column_matrix[k];
        }
    }
    // обратный ход
    double d_need, s;
    for (int k = n - 1; k >= 0; k--) {
        d_need = 0;
        for (int j = k + 1; j < n; ++j) {
            s = matrix[k + n * j] * d[j];
            d_need = d_need + s;
        }
        d[k] = (column_matrix[k] - d_need) / matrix[k + n * k];
    }


    for (int i = 0; i < n; ++i) {
        f_help = (F[i + 1] - F[i]) / (X[i + 1] - X[i]);
        c1 = F[i];
        c2 = d[i];
        c3 = (3 * f_help - 2 * d[i] - d[i + 1]) / (X[i + 1] - X[i]);
        c4 = (d[i] + d[i + 1] - 2 * f_help) / (X[i + 1] - X[i]) / (X[i + 1] - X[i]);
        koeff[i * 4] = c1;
        koeff[1 + i * 4] = c2;
        koeff[2 + i * 4] = c3;
        koeff[3 + i * 4] = c4;
    }
    free(d);
    free(column_matrix);
    free(matrix);
    return 0;
}
