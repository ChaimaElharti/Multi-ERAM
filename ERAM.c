#include "usefulFunctions.h"


// Réduction de degré m d’Arnoldi d’une matrice de taille n 
void ArnoldiProjectionMGS(int m, double* V, mat* A, double* H){

    double* z = (double*)malloc(A->nb_ligne * sizeof(double));
    

    double tempNorm = Norme(V, A->nb_ligne);

    for (int i = 0; i < A->nb_ligne; i++){
        V[i] = V[i]/tempNorm;
    }

    // Boucle principale sur j
    for (int j = 0; j < m; j++) {

        MatVecProd(A, &V[j*(A->nb_ligne)], &V[(j+1)*(A->nb_ligne)]);

    	for (int i = 0; i < j+1; i++) {

            H[i*m+j] = DotProduct(&V[(j+1)*(A->nb_ligne)], &V[i*(A->nb_ligne)], A->nb_ligne);

            for(int k=0; k<m+1; k++){
                V[(j+1)*(A->nb_ligne)+k] -= V[i*(A->nb_ligne)+k] * H[i*m+j];
            }
        }


        H[(j+1)*m+j] = Norme(&V[(j+1)*(A->nb_ligne)], (A->nb_ligne));


        for(int i=0;i<m+1; i++){
            V[(j+1)*(A->nb_ligne)+i] /= H[(j+1)*m+j];
        }
    }

    /*for (int i = 0; i < A->nb_ligne; i++){
        printf("%d\t", V[i]);
    }*/
}


void ERAM(int n, int m, double* v0, mat* A){

    double *Vm = calloc(n*(m+1),sizeof(double));
    double *H = calloc((m+1)*m,sizeof(double));

    // Vecteur propres
    double *wr = malloc(m*sizeof(double));
    double *wi = malloc(m*sizeof(double));
    // Matrices associées aux vecteurs propres
    double *vr = malloc(m*m*sizeof(double));
    // Vecteurs propres calculés
    double *u = calloc(n*m,sizeof(double));
    // Erreur
    double *rho = malloc(m*sizeof(double));
    double ritz;

    double* VmT;
    double* uT;

    lapack_int info;

    int converged = 0, max_it = 1, nb_it = 0;

    while(converged == 0 && nb_it < max_it){
        memcpy(Vm, v0, n * sizeof(double));

        ArnoldiProjectionMGS(m, Vm, A, H);

        // Calcul Valeurs/Vecteurs propres du sous espace H(m,m)

        info = LAPACKE_dgeev(LAPACK_ROW_MAJOR,'N','V', m, H, m, wr, wi, NULL, m, vr, m);
        if(info != 0){
            if(info < 0) { printf("the %d-th argument had a wrong value.", info); }
            else { printf("the QR algorithm failed to compute all the eigenvalues; elements %d:%d of WR and WI contain eigenvalues which have converged.", info, n); }
        }


        // Retour vers l'espace 

        // u = Vm . y
        VmT = transpose(m, n, Vm);
        MatMatprod(VmT, vr, n, m, m, u);

        // Calcul des coefficients de Ritz
        for (size_t i = 0; i < m; i++)
        {
            rho[i] = vr[(m-1)*m+i] * H[m*m+m-1];
            ritz += abs(rho[i]);
        }

        // Test de convergence
        if(ritz <= eps){
            converged = 1;
        }else{
            uT =  transpose(n, m, u);
            for (size_t i = 0; i < m; i++)
            {
                VecScalProd(&uT[i*n], wr[i], n);
            }
            for (size_t i = 0; i < m-1 ; i++)
            {
                VecAdd(&uT[i*n], &uT[(i+1)*n], n);
            }

            memcpy(v0, &uT[(m-1)*n], n* sizeof(double));
            
            memset(Vm, 0,  n*(m+1) * sizeof(double));
            memset(H, 0,  m*(m+1) * sizeof(double));
            memset(u, 0,  m*n * sizeof(double));
            ritz = 0;
        }
        nb_it += 1;
    }

    free(H);
    free(Vm);
}


int main(int argc, char **argv) {
    /*
    if (argc < 2) {
        return printf("usage: %s MATRIX\n ", argv[0]), 1;
    }

	mat* m = LectureMatrice(argv[1]);
	if (!m) {
    	return printf("error: failed to read matrix from file\n"), 2;
	}

    AfficheMatrice(m);

    // Initialisation du vecteur unitaire V1
    double* Vi = (double*)malloc(m->nb_ligne * sizeof(double));

    Vi[0] = 1;
    for(int i =1; i<m->nb_ligne;i++){Vi[i] = 0;}*/

   

    int it_m = 100; 
    int n = 3, m = 2;
    mat *A = initMat(n,n);
    A->a[0] = 2;A->a[1] = 5;A->a[2] = 5;
    A->a[3] = 5;A->a[4] = 7;A->a[5] = 1;
    A->a[6] = 0;A->a[7] = 0;A->a[8] = 7;

     double *v0 = calloc(n,sizeof(double));
    // Initialiser V0 : vecteur initial
    for (int i = 0; i < n; i++){
        v0[i] = 1;
    }

    ERAM(n, m, v0, A);




	//FreeMatrice(m);
    return 0;
}

