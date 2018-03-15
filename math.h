

#ifndef __MATH__M__
#define __MATH__M__
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
/* rand numbers*/
double ran1(long *idum);
double gamdev(int ia, long *idum);
double gasdev(long *idum);
/* matrix and vector functions -----------------------------------------------*/
extern double *mat  (int n, int m);
extern int    *imat (int n, int m);
extern double *zeros(int n, int m);
extern double *eye  (int n);
extern double dot (const double *a, const double *b, int n);
extern double norm(const double *a, int n);
extern void cross3(const double *a, const double *b, double *c);
extern int  normv3(const double *a, double *b);
extern void matcpy(double *A, const double *B, int n, int m);
extern void matmul(const char *tr, int n, int k, int m, double alpha,
                   const double *A, const double *B, double beta, double *C);
extern int  matinv(double *A, int n);
extern int  solve (const char *tr, const double *A, const double *Y, int n,
                   int m, double *X);
extern int  lsq   (const double *A, const double *y, int n, int m, double *x,
                   double *Q);
extern int  lusolve(double*A, double*Y, int n, double*X);
extern int  filter(double *x, double *P, const double *H, const double *v,
                   const double *R, int n, int m);
extern int  smoother(const double *xf, const double *Qf, const double *xb,
                     const double *Qb, int n, double *xs, double *Qs);
extern void matprint (const double *A, int n, int m, int p, int q);
extern void matfprint(const double *A, int n, int m, int p, int q, FILE *fp);
extern double str2num(const char *s, int i, int n);
#endif