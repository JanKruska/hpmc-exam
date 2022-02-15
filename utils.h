void RandomMatrix( int m, int n, double *ap, int lda );
void ZeroMatrix( int m, int n, double *ap, int lda );
void RandomMatrixSymmetric(int n, double *ap, int lda );
void printMatrix(int m, int n, double *ap, int lda);
int str2int(char * string);
void resymmetrize(int m, int n, double *ap, int lda, char uplo);
int issymetric(double *ap, int n, int m, int lda);