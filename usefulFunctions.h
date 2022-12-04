#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <string.h>
#include <lapacke.h>

#define eps 1e-15


typedef struct matrice {
    int nb_ligne;
    int nb_colonne;
    int taille;
    double* a;
}mat;


// === FONCTIONS POUR INITIALISER ET MANIPULER LES MATRICES === //
mat *initMat(const int nb_ligne, const int nb_colonne);
mat* LectureMatrice(const char *filename);
void AfficheMatrice(mat* m);
double* transpose(int row, int col, double* A);

// ======== OPÉRATIONS BLAS 1 ============//
void VecAdd(double* x, double* y, int n);
double DotProduct(double* x, double* y, int taille) ;
void VecScalProd(double* x, double alpha, int n);
double Norme(double* x, int n);


// ========= OPÉRATIONS BLAS 2 ===========//
void MatVecProd(mat* A, double* v, double* Av);


// ======== OPÉRATIONS BLAS 3 =============//
double* MatMatprod(double* A, double* B, int rowA, int colA, int colB, double* AB);