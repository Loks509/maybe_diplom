#include "matrix.h"
#include "integrals.h"
#include <stdio.h> 
#include <time.h> 

const double pi = acos(-1);
const double  k0 = 1, k1 = 1.5 * k0;

int n = 2, N = n * n;

double R = 5;
double lambda = 1;

//отрезок для двумерных интегральных уравнений
double A = 0, B = 1;
double C = 0, D = 1;
double E = 0, F = 1;

//шаг для одномерных интегральных уравнений
double h = (B - A) / n;

//отрезок для двумерных интегральных уравнений
double h1 = (B - A) / n, h2 = (D - C) / n, h3 = (F - E) / n;


inline double u(double y1, double y2, double y3) {     
    return y1 + y2 + y3;
}


void mk(double**& var, size_t type, size_t dim_s = 1) {
    //метод Коллокаций
    //var - матрица, в которую будут записываться данные 
    //type - тип уравнения Фредгольма, первого или второго рода
    //dim_s - размер измерения в котором считается метод Галёркина
    size_t M = _msize(var) / sizeof(var[0]);
    int m = (int) sqrt(M);

    if(dim_s == 1)
        for (size_t i = 0; i < M; i++){
            double ksi = A + (i + 0.5) * h;

            for (size_t j = 0; j < M; j++){
                double a = A + j * h, b = A + (j + 1.0) * h;
                
                var[i][j] = base_func(i, j) * (type - 1.0) - lambda * I_k(100, a, b, ksi);
            }
            var[i][M] = func(ksi);
        }
    else
        for (size_t i = 0; i < M; i++){
            short i1 = i / m, i2 = i % m;
            double ksi1 = A + (i1 + 0.5) * h1, ksi2 = A + (i2 + 0.5) * h2;

            for (size_t j = 0; j < M; j++) {
                short j1 = j / m, j2 = j % m;
                double a = A + j1 * h1, b = A + (j1 + 1.0) * h1, c = C + j2 * h2, d = C + (j2 + 1.0) * h2;

                var[i][j] = base_func(i, j) * (type - 1.0) - lambda * I_k(100, a, b, c, d, ksi1, ksi2);     
            }
            var[i][M] = func(ksi1, ksi2);
        }
}

void mg(double**& var, size_t type, size_t dim_s = 1) {
    //метод Галёркина
    //var - матрица, в которую будут записываться данные 
    //type - тип уравнения Фредгольма, первого или второго рода
    //dim_s - размер измерения в котором считается метод Галёркина
    size_t M = _msize(var) / sizeof(var[0]);
    int m = (int)sqrt(M);

    if (dim_s == 1) {
        for (size_t i = 0; i < M; i++) {
            double a = A + i * h, b = A + (i + 1.0) * h;

            for (size_t j = 0; j < M; j++) {
                double c = A + j * h, d = A + (j + 1.0) * h;

                var[i][j] = h * base_func(i, j) * (type - 1.0) - lambda * I<double>(100, k, a, b, c, d);
            }
            var[i][M] = I<double>(100, a, b);
        }
    }
    else if(dim_s == 2) {
        for (size_t i = 0; i < M; i++) {
            short i1 = i / m, i2 = i % m;
            double a = A + i1 * h1, b = A + (i1 + 1.0) * h1;
            double c = C + i2 * h2, d = C + (i2 + 1.0) * h2;

            //printf("i=%d\n", i);
            for (size_t j = 0; j < M; j++) {
                short j1 = j / m, j2 = j % m;
                double e = A + j1 * h1, f = A + (j1 + 1.0) * h1;
                double g = C + j2 * h2, l = C + (j2 + 1.0) * h2;

                //printf("j=%d\n", j);
                var[j][i] = h1 * h2 * base_func(i, j) * (type - 1.0) - lambda * I<double>(10, k, a, b, c, d, e, f, g, l);
            }
            var[i][M] = I<double>(100, func, a, b, c, d);
        }
    }
    else if (dim_s == 3) {
        int m = (int)pow(M, 1.0 / 3.0);

        for (size_t i = 0; i < M; i++) {
            int i1 = i % m, i2 = i / m - m * (i / m / m), i3 = i / m / m;

            double a = A + i1 * h1, b = a + h1;
            double c = C + i2 * h2, d = c + h2;
            double e = E + i3 * h3, f = e + h3;

            if(i % 100 == 0) printf("i=%d\n", i);
            for (size_t j = 0; j < M; j++) {
                int j1 = j % m, j2 = j / m - m * (j / m / m), j3 = j / m / m;

                double g = A + j1 * h1, l = g + h1;
                double o = C + j2 * h2, p = o + h2;
                double q = E + j3 * h3, s = q + h3;

                var[j][i] = h1 * h2 * h3 * base_func(i, j) * (type - 1.0) - lambda * I<double>(3, k, a, b, c, d, e, f, g, l, o, p, q, s);
            }
            var[i][M] = I<double>(3, func, a, b, c, d, e, f);
        }
    }
}


int main1()
{
    return 0;
    double** a = createm<double>(n, n + 1.0);

    //mk(a, 2, 2);
    mg(a, 2);

    cout << gm(a) << "\n";

    //print(a);

    space(1);

    /*for (size_t i = 0; i < N; i++)
    {
        short i1 = i / n, i2 = i % n;
        double ksi1 = A + (i1 + 0.5) * h1, ksi2 = A + (i2 + 0.5) * h2;
        
        printf("%f %f %f\n", ksi1, ksi2, a[i][N]);
    }*/

    for (size_t i = 0; i < n; i++)
    {
        printf("%f\n", a[i][n]);
    }
    return 0;
}

int main2() {
    return 0;
    double **a = createm<double>(N, N + 1.0), **a1 = createm<double>(N, N + 1.0);

    mg(a, 2.0, 2);
    
    mk(a1, 2.0, 2);

    //printf("mg\n");
    //print(a);
    //space(1);

    //printf("mk\n");
    //print(a1);
    //space(1);


    gm(a);


    gm(a1);

    space(1);

    printf("mg\n");
    for (size_t i = 0; i < N; i++)
    {
        if (i % n == 0 && i != 0) printf("\n");
        printf("%f ", a[i][N]);
        
    }
    space(1);

    printf("mk\n");
    for (size_t i = 0; i < N; i++)
    {
        /*short i1 = i / n, i2 = i % n;
        double ksi1 = A + (i1 + 0.5) * h1, ksi2 = A + (i2 + 0.5) * h2;
        printf("%f %f %f\n", ksi1, ksi2, a1[i][N]);*/
        
        if (i % n == 0 && i != 0) printf("\n");
        printf("%f ", a1[i][N]);
    }

    //printf("\n\n");
    //for (size_t i = 0; i < N; i++)
    //{
    //    int i1 = i/ n, i2 = i % n;
    //    printf("%f %f %f\n", A + (i1 + 0.5) * h1, C + (i2 + 0.5) * h2, a[i][N]);
    //}

}

int main() {
    //return 0;
    N *= n;
    double** a = createm<double>(N, N + 1.0), **b = createm<double>(N, N + 1.0);

    mg(a, 2.0, 3); 
    mg(b, 2.0, 3);
    //print(a, "g");
    //space(1);
    printf("cond: %f\n", cond(a));

    gm(a);

    double* init_approx = createv<double>(N);
    for (size_t i = 0; i < N; i++)
    {
        int i1 = i % n, i2 = i / n - n * (i / n / n), i3 = i / n / n;
        double x1 = A + (i1 + 0.5) * h1, 
               x2 = C + (i2 + 0.5) * h2,
               x3 = E + (i3 + 0.5) * h3;
        init_approx[i] = u(x1, x2, x3);
        printf("%f %f %f %f u=%f fabs()=%f\n", x1, x2, x3, a[i][N], init_approx[i], fabs(init_approx[i] - a[i][N]));
    }
    cout << "\n\n";


    double *ores = simple_iteration<double>(b, init_approx, 100);

    for (size_t i = 0; i < N; i++)
    {
        printf("%f u=%f fabs()=%f\n", ores[i], init_approx[i], fabs(init_approx[i] - ores[i]));
    }
}


