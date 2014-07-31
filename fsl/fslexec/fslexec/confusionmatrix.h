void createConfusionMatrix(int rank, int **matrix);
void zeraMatrix(int *matrix, int rank);
void copiaMatrix(int rank, int *destino, int *origem);
void signalResult(int desej, int obt, int rank, int *matrix);
float Kappa(int *matrix, int rank);
float chance(int rank);
float acertos(int *matrix, int rank);
int exemplos(int *matrix, int rank);
float erroRefclasse(int classe, int *matrix, int rank);
void saveMatrix(int *matrix, int rank, char **nomes, char *arquivo);

