typedef struct
{
	int inicio, fim;
	char nome[255];
} Intervalo;

void filtraclasses(int *classes, int tamclasses, int classeescolhida);
void deslocaarray(int *classes, int tamclasses, int desloc);
void salvanomes(char *arquivo, char **condicoes, int tamconds);
void salvanomes(char *arquivo, Intervalo *intervalos, int tamintervalos);
void preparaarquivocondicoes(char *arquivo, char *saida, float TR, float tamint, int numvolumes, int numcondicoes);
int PegaIndices(Intervalo *intervalos, int tamintervalos, int offset, int inicio, int fim, int **indices);
void performance(char *arquivo, int *classes, int numcondicoes, float **performances, int inicio=1);
void identificaclasses(char *arquivo, int *classes, int numcondicoes, int centros, int **maps, int *indices, int tamindices);
int retornaarrayclasses(Intervalo *intervalos, int tamintervalos, char **condicoes, int tamcondicoes, int **classes);
int pegalistacondicoes(Intervalo * intervalos, int tam, char ***lista);
void desalocacondicoes(char **condicoes, int tamconds);
int indicestring(char **lista, int tam, char* nome);
int pegadesenho(float * &desenho, Intervalo *intervalos, int tam, char *condicao, int inicio=0, int fim=0);
int retornaintervalo(int indice, Intervalo *intervalos, int tamintervalos);
int LeIntervalo(char *arquivo, Intervalo **intervalos);
int pegaindicesclasse(int *classes, int classeescolhida, int *indices, int tamindices, int **indicesclasse);
int juntaindices(int *indices1, int tamindices1, int *indices2, int tamindices2, int **indices);

void extractfilepath(char *file, char *saida);
void extractfilename(char *file, char *saida);
int fileexists(char *arquivo);