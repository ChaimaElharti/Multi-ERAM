#include "usefulFunctions.h"



// === FONCTIONS POUR INITIALISER ET MANIPULER LES MATRICES === //

mat *initMat(const int nb_ligne, const int nb_colonne){
    mat *m = (mat*)malloc(sizeof(mat));
    if (!m) {
        return NULL;
    }

    m->nb_ligne = nb_ligne;
    m->nb_colonne = nb_colonne;
    m->taille = nb_ligne * nb_colonne;

    m->a = (double*)(malloc(nb_ligne * nb_colonne * sizeof(double)));
	if (!m->a) {
    	return NULL;
	}

	return m;
}


mat* LectureMatrice(const char *filename){

    if (!filename) {
        printf("error: invalid function argument (filename is null)\n");
        exit(EXIT_FAILURE);
    }
    int nb_ligne = 0, nb_colonne = 0;
    double value = 0.0;

    // Open file
    FILE *fp = fopen(filename, "rb");
    if (!fp) {
    	printf("error: failed to open file %s\n", filename);
    	exit(EXIT_FAILURE);
    }

    // Read matrix dimensions from file
    fscanf(fp, "%u %u", &nb_ligne, &nb_colonne);
    if (!nb_ligne || !nb_colonne) {
        printf("error: failed to read from file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    mat* m = initMat(nb_ligne, nb_colonne);

    for (int i = 0; i < m->taille; i++) {
    	// read values from file
    	fscanf(fp, "%lf ", &value);
        m->a[i] = value;
	}
	fclose(fp);
    return m;
}


void AfficheMatrice(mat* m) {
    printf("La matrice est de taille %d * %d\n", m->nb_ligne, m->nb_colonne);
    
    for (int i = 0; i < m->nb_ligne; i++) {
    	for (int j = 0; j < m->nb_colonne; j++) {
        	printf("%f\t", m->a[i * m->nb_colonne + j]);
        }
        printf("\n");
    }
}

void FreeMatrice(mat* m) {
    if (!m) {
        printf("error: matrix is already null\n");
        return;
    }

    free(m->a);
    free(m);
}


// Transposition d'un matrice conforme au format requis pour l'opération LAPACK nécessaire à la résolution du problème
double* transpose(int row, int col, double* A){
    double* AT = calloc(row*col, sizeof(double));
    for(int i=0; i<row; i++) {
        for(int j=0; j<col; j++) {
            AT[j*row+i] = A[i*col+j];
        }
    }
    return AT;
}


// ======== OPÉRATIONS BLAS 1 ============//

void VecAdd(double* x, double* y, int n){
    for (int i = 0; i < n; i++){
        y[i] = x[i] + y[i];
    }
}

// Produit scalaire de deux vecteurs x, y de longueur n :
double DotProduct(double* x, double* y, int taille) {
    int sum = 0;
    for(int i = 0; i<taille; i++){
        sum += x[i] * y[i];
    }
    return sum;
}

// Produit d'un scalaire par un vecteur de longueur n
void VecScalProd(double* x, double alpha, int n) {
    for (int i = 0; i < n; i++){
        x[i] *=  alpha;
    }
}


// Norme d’un vecteur x de longueur n 
double Norme(double* x, int n){
    return sqrt(DotProduct(x, x, n));
}


// ========= OPÉRATIONS BLAS 2 ===========//

void MatVecProd(mat* A, double* v, double* Av){
    for(int i=0; i<A->nb_ligne; i++) {
        for(int j=0; j<A->nb_colonne; j++) {
            Av[i] += A->a[i*(A->nb_colonne)+j] * v[j];
        }
    }
}



// ======== OPÉRATIONS BLAS 3 =============//

double* MatMatprod(double* A, double* B, int rowA, int colA, int colB, double* AB){
    if(AB == NULL)
        AB = calloc(rowA*colB, sizeof(double));
    for(int i=0; i<rowA; i++) {
        for(int k=0; k<colB; k++) {
            for(int j=0; j<colA; j++) {
                AB[i*colB+k] += A[i*colA+j] * B[j*colB+k];
            }
        }
    }
    return AB;
}


