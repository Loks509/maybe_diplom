#pragma once


#ifndef MATRIX_H
#define MATRIX_H
#include <Windows.h>
#include <iostream>
#include <malloc.h>
#include <fstream>
#include <complex>
#include <string>
#include <math.h>

using namespace std;


template<class type_matrix>
inline type_matrix** createm(size_t M, size_t N, bool mod = false) {
    //выделение памяти под указатель
    //mod - нужноли при создании создавать единичную матрицу

    type_matrix** var = (type_matrix**)malloc(M * sizeof(type_matrix*));
    for (size_t i = 0; i < M; i++)
        var[i] = (type_matrix*)malloc(N * sizeof(type_matrix));



    if (mod)
        for (size_t i = 0; i < M; i++)
            for (size_t j = 0; j < N; j++)
                if (i == j) var[i][j] = 1;
                else var[i][j] = 0;

    return var;
}

template<class type_vector>
inline type_vector* createv(size_t N, bool mod = false) {
    //выделение памяти под указатель
    //mod - нужноли при создании создавать единичный вектор
    type_vector* var = (type_vector*)malloc(N * sizeof(type_vector));

    if (mod)
        for (size_t i = 0; i < N; i++)
            var[i] = 1;

    return var;
}

template<class type_matrix_print>
inline void print(type_matrix_print** var, string c = "", int cM = -1, int cN = -1) {
    //вывод указателя var в консоль
    //c (color) - цвет главной диагонали при выводе на экран: G, g - зелёный, B, b - синий, R, r - красный, I,i - интенсивнее серого
    size_t M, N;
    if (cM == -1 && cN == -1) {
        M = _msize(var) / sizeof(var[0]);
        N = _msize(var[0]) / sizeof(var[0][0]);
    }
    else {
        M = cM;
        N = cN;
    }

    const char* type = "";
    if (sizeof(type_matrix_print) == sizeof(int))  type = "%d ";
    if (sizeof(type_matrix_print) == sizeof(double))  type = "%f ";
    if (sizeof(type_matrix_print) == sizeof(complex<double>))  type = "complex";

    HANDLE hConsoleHandle = GetStdHandle(STD_OUTPUT_HANDLE);
    for (size_t i = 0; i < M; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            if (i == j) {
                if (c == "g" || c == "G")SetConsoleTextAttribute(hConsoleHandle, FOREGROUND_GREEN);
                if (c == "b" || c == "B")SetConsoleTextAttribute(hConsoleHandle, FOREGROUND_BLUE);
                if (c == "r" || c == "R")SetConsoleTextAttribute(hConsoleHandle, FOREGROUND_RED);
                if (c == "i" || c == "I")SetConsoleTextAttribute(hConsoleHandle, FOREGROUND_INTENSITY);

                if (type != "complex") printf(type, var[i][j]);
                else cout << var[i][j] << " ";

                fflush(stdout);
                SetConsoleTextAttribute(hConsoleHandle, 15);
            }
            else {
                if (type != "complex") printf(type, var[i][j]);
                else cout << var[i][j] << " ";

                fflush(stdout);
            }

        }
        printf("\n");
        fflush(stdout);
    }
    printf("\n");
    fflush(stdout);
}

template<class type_vector_print>
inline void print(type_vector_print* var) {
    //вывод указателя var в консоль
    size_t M = _msize(var) / sizeof(var[0]);

    const char* type = "";
    if (sizeof(type_vector_print) == sizeof(int))  type = "%d\n";
    if (sizeof(type_vector_print) == sizeof(double))  type = "%f\n";
    if (sizeof(type_vector_print) == sizeof(complex<double>))  type = "complex";

    for (size_t i = 0; i < M; i++)
    {
        if (type != "complex") printf(type, var[i]);
        else cout << var[i] << "\n";

        fflush(stdout);
    }
    printf("\n");
    fflush(stdout);
}

void space(size_t k = 0) {
    //вывод пробелов в консоль
    for (size_t ind_k = 0; ind_k < k; ind_k++)
    {
        printf("\n");
        fflush(stdout);
    }
}

template<class type_matrix_size>
int size(type_matrix_size** var) {
    //вычисляет объём занимаемой памяти указателем var
    size_t M = _msize(var) / sizeof(var[0]);
    size_t sum = 0;
    for (size_t i = 0; i < M; i++)
        sum += _msize(var[i]);
    return sum;
}

template<class type_vector_size>
int size(type_vector_size* var) {
    //вычисляет объём занимаемой памяти указателем var
    return _msize(var);
}

template<class type_gm>
string gm(type_gm**& var) {
    //метод Гаусса
    size_t M = _msize(var) / sizeof(var[0]);
    size_t N = _msize(var[0]) / sizeof(var[0][0]);

    if (M + 1 != N) return "Error";

    for (size_t k = 0; k < M; k++) {
       /* if (k == 0) printf("\nk=%i\n", k);
        if (k != 0) printf("k=%i\n", k);
        fflush(stdout);*/
        type_gm ed = 1;
        if (var[k][k] != ed) {
            type_gm T = var[k][k];
            for (size_t j = k; j < N; j++) {
                var[k][j] = var[k][j] / T;
            }
        }
        for (size_t i = k; i < M; i++) {
            if ((var[i][k] != ed) && (i != k)) {
                type_gm T = var[i][k];
                var[i][k] = 0;
                for (size_t j = k + 1; j < N; j++) {
                    var[i][j] -= var[k][j] * T;
                }
            }
        }
    }
    for (int i = M - 1; i >= 0; i--) {
        type_gm Sum = var[i][M];
        for (size_t j = i + 1.0; j < M; j++) {
            Sum -= var[i][j] * var[j][M];
        }
        var[i][M] = Sum;
    }

    return "Successfully";
}

template<class type_fill>
type_fill* fill(int M, type_fill item = 0.0) {
    type_fill* var = createv<type_fill>(M);
    for (size_t i = 0; i < M; i++) {
        var[i] = item;
    }
    return var;
}

template<class type_simple_iteration>
type_simple_iteration* simple_iteration(type_simple_iteration** var = 0, type_simple_iteration* init_approx = 0, size_t k = 0)
{
    size_t M = _msize(var) / sizeof(var[0]);
    type_simple_iteration* res = createv<type_simple_iteration>(M);
    size_t k_ind = 0;
    while (true) {
        res = fill<type_simple_iteration>(M);
        k_ind++;
        for (size_t i = 0; i < M; i++)
        {
            if (var[i][i] > 0) {
                for (size_t j = 0; j < M; j++)
                {
                    if (i == j)res[i] -= (var[i][j] - 1) * init_approx[j];
                    if (i != j)res[i] -= var[i][j] * init_approx[j];
                }
                res[i] += var[i][M];
            }
            if (var[i][i] < 0) {
                for (size_t j = 0; j < M; j++)
                {
                    if (i == j) res[i] += (var[i][j] + 1) * init_approx[j];
                    if (i != j) res[i] += var[i][j] * init_approx[j];
                }
                res[i] -= var[i][M];
            }
        }
        if (k_ind <= k) {
            init_approx = res;
            res = fill<type_simple_iteration>(M);
        }
        else {
            return res;
        }
    }
}


template<typename type_matrix_del>
void del(type_matrix_del**& var) {
    //очистка памяти указателя var
    size_t M = _msize(var) / sizeof(var[0]);

    for (size_t i = 0; i < M; i++) free(var[i]);
    free(var);
}

template<typename type_vector_del>
void del(type_vector_del*& var) {
    //очистка памяти указателя var
    free(var);
}

//ниже возможно ошибки

template<class type_LU>
void LU(type_LU** var, type_LU**& L, type_LU**& U) {
    size_t M = _msize(var) / sizeof(var[0]);

    for (size_t i = 0; i < M; i++) {
        for (size_t j = 0; j < M; j++) {
            if (i == j) L[i][i] = 1;
            else L[i][j] = 0;
            U[i][j] = 0;
        }
    }
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            type_LU sumU = 0, sumL = 0;
            if (i <= j) {
                for (int z = 0; z <= i - 1; z++) {
                    sumU += L[i][z] * U[z][j];
                }
                U[i][j] = var[i][j] - sumU;
            }

            if (i > j) {
                for (int z = 0; z <= j - 1; z++) {
                    sumL += L[i][z] * U[z][j];
                }
                L[i][j] = (var[i][j] - sumL) / U[j][j];
            }
        }
    }
}

template<class type_mult>
inline void mult(type_mult** var1, type_mult** var2, type_mult**& res) {
    size_t M = _msize(var1) / sizeof(var1[0]);
    size_t N = _msize(var1[0]) / sizeof(var1[0][0]);

    for (size_t i = 0; i < M; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            type_mult sum = 0;
            for (size_t k = 0; k < M; k++)
            {
                sum += var1[i][k] * var2[k][j];
            }
            res[i][j] = sum;
        }
    }
}

template<class type_diag>
type_diag** diag(type_diag** var) {
    size_t M = _msize(var) / sizeof(var[0]);

    type_diag** var_D = createm<type_diag>(M, M), ** var_L = createm<type_diag>(M, M), ** var_U = createm<type_diag>(M, M);

    LU(var, var_L, var_U);
    mult(var_U, var_L,var_D);

    for (size_t k = 1; k <= 10; k++) {
        LU(var_D, var_L, var_U);
        mult(var_U, var_L, var_D);
    }

    del(var_L); del(var_U);
    for (size_t i = 0; i < M; i++)
    {
        for (size_t j = 0; j < M; j++)
        {
            if (i != j) var_D[i][j] = 0;
        }
    }
    return var_D;
}

template<class type_eigenvalues>
type_eigenvalues* eigenvalues(type_eigenvalues** var) {
    size_t M = _msize(var) / sizeof(var[0]);

    type_eigenvalues** var_D = diag<type_eigenvalues>(var), * res = createv<type_eigenvalues>(M);
    for (size_t i = 0; i < M; i++)
        res[i] = var_D[i][i];

    del(var_D);
    return res;
}

template<class type_cond>
type_cond cond(type_cond** var) {
    size_t M = _msize(var) / sizeof(var[0]);


    type_cond* var1 = eigenvalues(var);

    type_cond max = 0, min = 0;
    for (size_t i = 0; i < M; i++)
    {
        if (abs(var1[i]) >= max) max = abs(var1[i]);
        if (abs(var1[i]) <= max) min = abs(var1[i]);
    }

    del(var1);
    return max / min;
}
#endif MATRIX_H