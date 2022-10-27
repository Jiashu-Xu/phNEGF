//
// Created by xjs on 2022/10/23.
//

#include "surfGFph3.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <malloc.h>
#define eps pow(2,-52)
void getNewMatrix(double ** expand_mat, double ** new_mat, int len);
void initExpandMatrix(double ** mat, double ** expand_mat, int len);
void calculateExpandMatrix(double ** expand_mat, int len);
int adjustMatrix(double ** expand_mat, int len);
double** inverse(double ** mat, int len) {
    // 扩展矩阵定义 //
    double **expand_matrix = (double**)malloc(sizeof(double *) * len);
    for (int i = 0; i < len; i++) {
        expand_matrix[i] = (double *)malloc(sizeof(double ) * (len * 2));
    }

    // 逆矩阵定义 //
    double **new_matrix = (double**)malloc(sizeof(double *) * len);
    for (int i = 0; i < len; i++) {
        new_matrix[i] = (double *)malloc(sizeof(double ) * len);
    }

    // init
    initExpandMatrix(mat, expand_matrix, len);

    // adjust
    int canAdjust = adjustMatrix(expand_matrix, len);

    if (canAdjust == 0) {
        return NULL;
    }

    // calc expand
    calculateExpandMatrix(expand_matrix, len);

    // 取后面的N*N矩阵，即为所求 //
    getNewMatrix(expand_matrix, new_matrix, len);

    return new_matrix;
}


// init expand_matrix
void initExpandMatrix(double ** mat, double ** expand_mat, int len) {
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < len * 2; j++) {
            if (j < len) {
                expand_mat[i][j] = mat[i][j];
            } else {
                if (j == len + i) {
                    expand_mat[i][j] = 1;
                } else {
                    expand_mat[i][j] = 0;
                }
            }
        }
    }
}


// adjust expand matrix
int adjustMatrix(double ** expand_mat, int len) {

    for (int i = 0; i < len; i++) {
        if (expand_mat[i][i] == 0) {
            int j;
            for (j = 0; j < len; j++) {
                if (expand_mat[j][i] != 0) {
                    double* tmp = expand_mat[i];
                    expand_mat[i] = expand_mat[j];
                    expand_mat[j] = tmp;
                    break;
                }
            }
            if (j >= len) {
                printf("Inv Matrix does not exists\n");
                return false;
            }
        }
    }
    return true;
}


// calc
void calculateExpandMatrix(double ** expand_mat, int len) {
    for (int i = 0; i < len; i++) {
        double fir_ele = expand_mat[i][i];
        for (int j = 0; j < len * 2; j++) {
            expand_mat[i][j] /= fir_ele;  // 该行所有元素除以首元素 //
        }
        for (int m = 0; m < len; m++) {
            if (m == i) {
                continue;
            }
            // 倍数 //
            double times = expand_mat[m][i];
            for (int n = 0; n < len * 2; n++) {
                expand_mat[m][n] -= expand_mat[i][n] * times;
            }
        }
    }
}


// get res
void getNewMatrix(double ** expand_mat, double ** new_mat, int len) {
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < len * 2; j++) {
            if (j >= len) {
                new_mat[i][j - len] = expand_mat[i][j];
            }
        }
    }
}




// matrix multiplying
double** getProductMatrix(double** init_mat, double** new_mat, int len) {
    double** product_mat = (double**)malloc(sizeof(double*) * len);
    for (int i = 0; i < len; i++) {
        product_mat[i] = (double*)malloc(sizeof(double) * len);
    }
    // need initializing to zero
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < len; j++) {
            product_mat[i][j] = 0;
        }
    }
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < len; j++) {
            for (int k = 0; k < len; k++) {
                product_mat[i][j] += init_mat[i][k] * new_mat[k][j];
            }
        }
    }
    return product_mat;
}
double **add(double **c,double **d) {
    int i,j,n,m;
    double **a;
    n = sizeof(c[0]) / sizeof(c[0][0]);
    m = sizeof(c) / sizeof(c[0]);
    a=(double *) malloc(sizeof(double)*n);
    for(int i=0;i<n;i++)
        a[i]=(double *)malloc(sizeof(double)*m);
    for(i=0; i<n; i++) {
        for(j=0; j<m; j++) {
            a[i][j]=d[i][j]+c[i][j];
        }
    }
    return a;
}
double **minus(double **c,double **d) {
    int i,j,n,m;
    double **a;
    n = sizeof(c[0]) / sizeof(c[0][0]);
    m = sizeof(c) / sizeof(c[0]);
    a=(double *) malloc(sizeof(double)*n);
    for(int i=0;i<n;i++)
        a[i]=(double *)malloc(sizeof(double)*m);
    for(i=0; i<n; i++) {
        for(j=0; j<m; j++) {
            a[i][j]=c[i][j]-d[i][j];
        }
    }
    return a;
}
double **multiply(double **c,double **d) {
    int i,j,k,m,n,r;
    double **a;
    m = sizeof(c[0]) / sizeof(c[0][0]);
    r = sizeof(c) / sizeof(c[0]);
    n = sizeof(d) / sizeof(d[0]);
    a=(double *) malloc(sizeof(double)*m);
    for(int i=0;i<n;i++)
        a[i]=(double *)malloc(sizeof(double)*n);
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            a[i][j] = 0;
        }
    }
    for(i=0; i<m; i++) {
        for(j=0; j<n; j++) {
            for(k=0; k<r; k++) {
                a[i][j] += c[i][k]*d[k][j];
            }
        }
    }
    return a;
}
double maxabs(double **a){
    int i,j,count=0,m,n;
    m = sizeof(a[0]) / sizeof(a[0][0]);
    n = sizeof(a) / sizeof(a[0]);
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            if(abs(a[i][j])>count)
            {
                count = abs(a[i][j]);
            }
        }
    }
    return count;
}
double** surfGFph3(double ene, double eta, double** V10, double** V01, double** h0) {
    int i, j, k;
    int dim = sizeof(h0[0]) / sizeof(h0[0][0]);
    double** c, ** d, ** e, ** f, ** I;
    I =(double **)malloc(sizeof(double *) * dim);
    for (i=0; i < dim; i++)
        I[i]=(double *)malloc(sizeof(double) * dim);
    for (i = 0; i < dim; i++)
    {
        for (j = 0; j < dim; j++)
        {
            if (i == j)
            {
                I[i][j] = 1;
            }
            else I[i][j] = 0;
        }
    }
    double** l;
    l = (double**)malloc(sizeof(double*) * dim);
    for (i=0; i < dim; i++)
        l[i] = (double*)malloc(sizeof(double) * dim);
    for (i = 0; i < dim; i++)
    {
        for (j = 0; j < dim; j++)
        {
            l[i][j] = (ene +  eta) * I[i][j];
        }
    }
    c = minus(l, h0);
    double** cs, ** c0;
    cs = c;
    double** a, ** b;
    a = V10;
    b = V01;
    int n = 0;
    while (maxabs(b) > eps || maxabs(a) > eps)
    {
        c0 = c;
        d = multiply(multiply(b, inverse(c0, sizeof(c0[0]) / sizeof(c0[0][0]))), a);
        cs = minus(cs, d);
        c = minus(minus(c0, multiply(multiply(a, inverse(c0, sizeof(c0[0]) / sizeof(c0[0][0]))), b)), multiply(multiply(b, inverse(c0, sizeof(c0[0]) / sizeof(c0[0][0]))), a));
        a = multiply((a, inverse(c0, sizeof(c0[0]) / sizeof(c0[0][0]))), a);
        b = multiply((b, inverse(c0, sizeof(c0[0]) / sizeof(c0[0][0]))), b);
        n++;
    }
    /*if (rcond(a) < 1e-15)//抄的matlab
    {
        k = sizeof(ene[0]) / sizeof(ene[0][0]);
        for (i = 0; i < k; i++)
        {
            for (j = 0; j < k; j++)
            {
                if (i == j)
                {
                    e[i][j] = 1;
                }
                else e[i][j] = 0;
            }
        }
    }*/
    f = (double**)malloc(sizeof(double*) * dim);
    for (i=0; i < dim; i++)
        f[i] = (double*)malloc(sizeof(double) * dim);
    for (i = 0; i < dim; i++)
    {
        for (j = 0; j < dim; j++)
        {
            if (i == j)
            {
                f[i][j] = 1;
            }
            else f[i][j] = 0;
        }
    }
    return multiply(f, inverse(cs, sizeof(cs[0]) / sizeof(cs[0][0])));
}