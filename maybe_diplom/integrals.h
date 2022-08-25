#pragma once
#ifndef INTEGRALS_H
#define INTEGRALS_H
#include <math.h>
#include <stdarg.h>

double mach_eps = sqrt(2.22045e-16);

inline double k(double y1, ...) {
    double y2, x1, x2;
    va_list args;

    va_start(args, y1);
    y2 = va_arg(args, double);
    x1 = va_arg(args, double);
    x2 = va_arg(args, double);
    va_end(args);


    return x1 + x2 - y1 * y2;
    //return y1 + y2 - x1 * x2;
    //return y1 - y2;
}

inline double func(double x1, ...) {
    va_list args;
    va_start(args, x1);
    double x2 = va_arg(args, double);
    va_end(args);

    //return sin(10.0 * x1 + 10.0 * x2) + (1.0 / 10000.0) * (-200.0 * x1 - 200.0 * x2 - 2.0) * sin(10.0) + (1.0 / 10000.0) * (100.0 * x1 + 100.0 * x2 - 99.0) * sin(20.0) + (1.0 / 500.0) * cos(10.0) - (1.0 / 500.0) * cos(20.0);
    //return (1.0 / 6.0) * (-6.0 * x1 - 6.0 * x2 + 3.0) * cos(1.0) + (1.0 / 2.0) * x1 + (3.0 / 2.0) * x2 - (1.0 / 2.0) * sin(1.0) - sin(x1) + 1.0 / 6.0;
    return 1.0 / 3.0;
    //return pow(x1, 2) - 1.0 * x1 / 3.0 + 1.0 / 4.0;
}

inline double base_func(int i, int j) {
    return i == j ? 1.0 : 0.0;
}

inline void grad(double(*function)(double, ...), double*& res, double x1, double x2) {
    res[0] = ((*function)(x1 * (1 + mach_eps), x2) - (*function)(x1 * (1 - mach_eps), x2)) / 2.0 / x1 / mach_eps;
    res[1] = ((*function)(x1, x2 * (1 + mach_eps)) - (*function)(x1, x2 * (1 - mach_eps))) / 2.0 / x2 / mach_eps;
}

inline double div(double(*function)(double, ...), double x1, double x2) {
    double* res = createv<double>(2);

    grad((*function), res, x1, x2);
    return res[0] + res[1];
}


inline double I_k(int N_i, double a, double b, double ksi) {
    //одномерный интеграл для метода Коллокаций 
    double h_int = (b - a) / (N_i * 1.0);
    double Sum = 0.0;
    for (size_t i = 0; i < N_i; i++){
        double l = (i + 0.5) * h_int + a;
        Sum += k(ksi, l) * h_int;
    }
    return Sum;
}

inline double I_k(int N_i, double a, double b, double c, double d, double ksi1, double ksi2) {
    //двумерный интеграл для метода Коллокаций
    double h_i = (b - a) / N_i, h_j = (d - c) / N_i;
    double Sum = 0;
    for (size_t i = 0; i < N_i; i++) {
        for (size_t j = 0; j < N_i; j++) {
            double l1 = a + (i + 0.5) * h_i, l2 = c + (j + 0.5) * h_j;
            Sum += k(l1, l2, ksi1, ksi2);
        }
    }
    return h_i * h_j * Sum;
}

inline double I(int N_i, double a, double b) {
    //одномерный интеграл
    double h_int = (b - a) / (N_i * 1.0);
    double Sum = 0.0;
    for (size_t i = 0; i < N_i; i++) {
        double l = (i + 0.5) * h_int + a;
        Sum += func(l);
    }
    return Sum * h_int;
}

inline double I(int N_i, double(*function)(double, ...), double a, double b, double c, double d) {
    //двумерный интеграл
    double h_i = (b - a) / N_i, 
           h_j = (d - c) / N_i;
    double Sum = 0;
    for (size_t i = 0; i < N_i; i++) {
        for (size_t j = 0; j < N_i; j++) {
            double l1 = a + (i + 0.5) * h_i, 
                   l2 = c + (j + 0.5) * h_j;
            Sum += (*function)(l1, l2);
        }
    }
    return h_i * h_j * Sum;
}

inline double I(int N_i, double(*function)(double, ...), double a, double b, double c, double d, double e, double f) {
    //трёхмерный интеграл
    double h_i = (b - a) / N_i,
           h_j = (d - c) / N_i,
           h_k = (f - e) / N_i;

    double Sum = 0;
    for (size_t i = 0; i < N_i; i++) {
        for (size_t j = 0; j < N_i; j++) {
            for (size_t k = 0; k < N_i; k++) {
                double l1 = a + (i + 0.5) * h_i,
                       l2 = c + (j + 0.5) * h_j,
                       l3 = e + (k + 0.5) * h_k;

                Sum += (*function)(l1, l2, l3);
            }
        }
    }
    return h_i * h_j * h_k * Sum;
}

inline double I(int N_i, double(*function)(double, ...), double a, double b, double c, double d, double e, double f, double g, double h) {
    //четырёхмерный интеграл
    double h_i = (b - a) / N_i, 
           h_j = (d - c) / N_i,
           h_k = (f - e) / N_i,
           h_z = (h - g) / N_i;

    double Sum = 0;
    for (size_t i = 0; i < N_i; i++) {
        for (size_t j = 0; j < N_i; j++) {
            for (size_t k = 0; k < N_i; k++) {
                for (size_t z = 0; z < N_i; z++) {
                    
                    double l1 = a + (i + 0.5) * h_i,
                           l2 = c + (j + 0.5) * h_j,
                           l3 = e + (k + 0.5) * h_k,
                           l4 = g + (z + 0.5) * h_z;
                    Sum += (*function)(l1, l2, l3, l4);
                }
            }
        }
    }
    return h_i * h_j * h_k * h_z * Sum;
}

inline double I(int N_i, double(*function)(double, ...), double a, double b, double c, double d, double e, double f, double g, double h, double l, double o, double p, double q) {
    //шестимерый интеграл
    double h_i = (b - a) / N_i,
           h_j = (d - c) / N_i,
           h_k = (f - e) / N_i,
           h_z = (h - g) / N_i,
           h_v = (o - l) / N_i,
           h_w = (q - p) / N_i;

    double Sum = 0;
    for (size_t i = 0; i < N_i; i++) {
        for (size_t j = 0; j < N_i; j++) {
            for (size_t k = 0; k < N_i; k++) {
                for (size_t z = 0; z < N_i; z++) {
                    for (size_t v = 0; v < N_i; v++) {
                        for (size_t w = 0; w < N_i; w++) {

                            double l1 = a + (i + 0.5) * h_i,
                                   l2 = c + (j + 0.5) * h_j,
                                   l3 = e + (k + 0.5) * h_k,
                                   l4 = g + (z + 0.5) * h_z,
                                   l5 = l + (v + 0.5) * h_v,
                                   l6 = p + (w + 0.5) * h_w;
                            Sum += (*function)(l1, l2, l3, l4, l5, l6);
                        }
                    }
                }
            }
        }
    }
    return h_i * h_j * h_k * h_z * h_v * h_w * Sum;
}

#endif INTEGRALS_H
