#include "intervalos.h"
#include <string>
#include <string.h>
#include <stdlib.h>
#include <stdio.h> 
#include <math.h> 
#include <fstream>
#include <list>
#include <algorithm>

using namespace std;
using std::list;

int fileexists(char *arquivo)
{
	FILE *f = fopen(arquivo, "rb");
	if (f != NULL)
	{
		fclose(f);
		return 1;
	}
	else return 0;
}

void extractfilepath(char *file, char *saida)
{
	strcpy(saida, file);
	for(int i=strlen(saida)-1;i--; i >= 0)
		if (saida[i] == '\\') 
		{
			saida[i]=0;
			break;
		}
}

void extractfilename(char *file, char *saida)
{
	for(int i=strlen(file)-1;i--; i >= 0)
		if (file[i] == '\\') 
		{
			strcpy(saida, (file + i+1));
			break;
		}
}

int LeIntervalo(char *arquivo, Intervalo **intervalos)
{
	ifstream in;
	string linha;
	in.open(arquivo);
	
    int tamanho=0;
	*intervalos=NULL;
	while (getline(in, linha))
	{
       *intervalos= (Intervalo *) realloc(*intervalos, ++tamanho*sizeof(Intervalo));
	   if (*intervalos != NULL)
	   {
           char num[20];
		   int pos = linha.find("-");
		   int pos2= linha.find(",");
		   linha.copy(num, pos);
		   num[pos]='\0';
           (*intervalos)[tamanho-1].inicio = atoi(num);
		   linha.copy(num, pos2-pos-1, pos+1);
		   num[pos2-pos-1]='\0';
           (*intervalos)[tamanho-1].fim = atoi(num);
		   linha.copy((*intervalos)[tamanho-1].nome, linha.length()-pos2-1, pos2+1);
		   (*intervalos)[tamanho-1].nome[linha.length()-pos2-1] = '\0';
	   }
	}
	in.close();
	return tamanho;
}

int retornaintervalo(int indice, Intervalo *intervalos, int tamintervalos)
{
	int ret = -1;
	for (int i=0;i<tamintervalos;i++)
	{
		if ((indice >= intervalos[i].inicio) && (indice <= intervalos[i].fim))
		{
			ret = i;
			break;
		}
	}
	return ret;
}

int pegadesenho(float * &desenho, Intervalo *intervalos, int tam, char *condicao, int inicio, int fim)
{
	if (inicio==0) inicio = 1;
	if (fim==0) fim = intervalos[tam-1].fim;
	desenho=(float *) malloc((fim-inicio+1) * sizeof(float));
	for (int i=inicio;i<=fim;i++)
	{
		int interv = retornaintervalo(i, intervalos, tam);
		if (interv >= 0)
		{
		   if (strcmp(condicao, intervalos[interv].nome)==0) desenho[i-inicio] = 1;
			else desenho[i-inicio]=0;
		}
	}
	return (fim-inicio+1);
}

int indicestring(char **lista, int tam, char* nome)
{
	int r = -1;
	for (int i = 0; i < tam; i++)
	{
		if (strcmp(lista[i], nome) == 0)
		{
			r = i;
			break;
		}
	}
	return r;
}

void desalocacondicoes(char **condicoes, int tamconds)
{
	for (int i=0; i< tamconds; i++)	free(condicoes[i]);
	free(condicoes);
}

int pegalistacondicoes(Intervalo * intervalos, int tam, char ***lista)
{
	int tamanho = 0;
	*lista = NULL;
	for (int i = 0; i< tam; i++)
	{
		if (tamanho == 0) 
		{
			*lista = (char **) realloc(*lista, ++tamanho * (sizeof(char *)));
			(*lista)[tamanho-1] = (char *) malloc((strlen(intervalos[i].nome)+1) * sizeof(char));
			strcpy((*lista)[tamanho-1], intervalos[i].nome);
		}
		else if (indicestring(*lista, tamanho, intervalos[i].nome) == -1)
		{
			*lista = (char **) realloc(*lista, ++tamanho * (sizeof(char *)));
			(*lista)[tamanho-1] = (char *) malloc((strlen(intervalos[i].nome)+1) * sizeof(char));
			strcpy((*lista)[tamanho-1], intervalos[i].nome);
		}
	}
	return tamanho;
}

int retornaarrayclasses(Intervalo *intervalos, int tamintervalos, char **condicoes, int tamcondicoes, int **classes)
{
	int tamanho = intervalos[tamintervalos-1].fim;
	*classes = (int *) malloc(tamanho * sizeof(int));
	for (int i=0;i<tamintervalos;i++)
	{
		int classe = indicestring(condicoes, tamcondicoes, intervalos[i].nome)+1;
		for (int j=intervalos[i].inicio;j<=intervalos[i].fim;j++) (*classes)[j-1] = classe;
	}
	return tamanho;
}

void identificaclasses(char *arquivo, int *classes, int numcondicoes, int centros, int **maps, int *indices, int tamindices)
{
	FILE *f;
	char linha[255];
	int predicao;
	int *predicoes;
	float valor;
	int **mapeamentos;
	int *map;
	int *totais;
	int i=0;

	map = (int *) malloc((centros) * sizeof(int));
	totais = (int *) malloc((centros) * sizeof(int));
	mapeamentos = (int **) malloc((centros) * sizeof(int *));
	predicoes = (int *) malloc(tamindices * sizeof(int));
	for (int t=0;t<centros;t++)
	{
		totais[t] = 0;
		mapeamentos[t] = (int *)malloc((numcondicoes) * sizeof(int));
	    for (int j=0;j<numcondicoes;j++) mapeamentos[t][j]=0;
	}

	f=fopen(arquivo, "rt");
	if (f != NULL)
	{
		while (fgets(linha, 255, f))
		{
		   sscanf(linha, "%d %f", &predicao, &valor);
		   mapeamentos[predicao][classes[indices[i]-1]-1]++;
		   totais[predicao]++;
		   predicoes[i]=predicao;
		   i++;
		}
		fclose(f);
		for(int t=0;t<centros;t++) 
		{
			int winner=-1;
			float maxval=0;
			for (int j=0;j<numcondicoes;j++)
			{
				float valor = 0;
				if (totais[t]!=0) valor = (100.0*mapeamentos[t][j]/(float) totais[t]);
				if ((winner==-1) || (maxval < valor))
				{
					maxval = valor;
				 	winner = j;
				}
			}
			if (maxval > 0) map[t] = winner;
			else map[t] = -1;
		}
	}
	*maps=(int *) malloc(tamindices*sizeof(int));
	for (int t=0;t<tamindices;t++)
	{
		if (map[predicoes[t]]==classes[t]-1) (*maps)[t] = classes[t];
		else (*maps)[t]=0-classes[t];
	}
	free(totais);
	free(map);
	for (int t=0;t<numcondicoes;t++) free(mapeamentos[t]);
	free(mapeamentos);
}

void performance(char *arquivo, int *classes, int numcondicoes, float **performances, int inicio)
{
	FILE *f;
	f=fopen(arquivo, "rt");
	char linha[255];
	float predicao, valor;
	int *acertos, *totais;
	int i=0;

	*performances = (float *) malloc((numcondicoes+1) * sizeof(float));

	acertos = (int *) malloc((numcondicoes+1) * sizeof(int));
	totais  = (int *) malloc((numcondicoes+1) * sizeof(int));
	for(int t=0;t<=numcondicoes;t++)
	{
		acertos[t]=0;
		totais[t]=0;
	}

	inicio--;
	if (f != NULL)
	{
		while (fgets(linha, 255, f))
		{
		   sscanf(linha, "%f %f", &predicao, &valor);
		   if (classes[inicio+i]==predicao)
		   {
			   acertos[0]++;
			   acertos[classes[inicio+i]]++;
		   }
		   totais[0]++;
		   totais[classes[inicio+i]]++;
		   i++;
		}
		for(int t=0;t<=numcondicoes;t++) (*performances)[t] = ((float)acertos[t]*100.0/(float)totais[t]);
		fclose(f);
	}
}

int PegaIndices(Intervalo *intervalos, int tamintervalos, int offset, int inicio, int fim, int **indices)
{
	int tamanho=0;
	tamanho = (fim-inicio+1);
	*indices = (int *) malloc(tamanho * sizeof(int));

	int h=0;
	tamanho=0;
	for (int i=0;i<tamintervalos;i++)
	{
		for (int j=(intervalos[i].inicio+offset);j<=intervalos[i].fim;j++) 
			if ((j >= inicio) && (j <= fim)) 
			{
				(*indices)[h++] = j;
				tamanho++;
			}
	}
	return tamanho;
}

void preparaarquivocondicoes(char *arquivo, char *saida, float TR, float tamint, int numvolumes, int numcondicoes)
{
   typedef struct
   {
	   float *onsets;
	   char nome[255];
	   int tamanho, indice;
   } condicao;

    condicao *condicoes;

	condicoes = (condicao *) malloc(numcondicoes * sizeof(condicao));

    FILE *f;
	f=fopen(arquivo, "rt");
	char linha[5000];
	if (f != NULL)
	{
		int indlinha=0;
		while (fgets(linha, 5000, f))
		{
		    int tamanho=0;
			char *auxlinha, *ultimo=NULL;
			sscanf(linha, "%80s", condicoes[indlinha].nome);
			condicoes[indlinha].onsets = NULL;

			auxlinha = strchr(linha, ' ');
			do
			{
				auxlinha++;
				if (auxlinha != NULL)
				{
				   ++tamanho;
				   condicoes[indlinha].onsets = (float *) realloc(condicoes[indlinha].onsets, tamanho * sizeof(float));
				   sscanf(auxlinha, "%f", &condicoes[indlinha].onsets[tamanho-1]);
				   ultimo=auxlinha;
				}
			}
			while ((auxlinha !=NULL) && (auxlinha=strchr(auxlinha, ' ')));
		/*
			if (strlen(ultimo))
			{
				++tamanho;
				condicoes[indlinha].onsets = (float *) realloc(condicoes[indlinha].onsets, tamanho * sizeof(float));
				sscanf(ultimo, "%f", &condicoes[indlinha].onsets[tamanho-1]);
			}
			*/
			condicoes[indlinha].tamanho = tamanho;
			condicoes[indlinha].indice = 0;
			indlinha++;
		}
		fclose(f);
	}

	f=fopen(saida, "wt+");
	int ultimofim=0;
	if (f != NULL)
	{
		while (1)
		{
			int minind=-1;
			float minval=0;
	        for (int t=0; t < numcondicoes; t++)
	        {
			   if (condicoes[t].indice < condicoes[t].tamanho)
			   {
			      if ((minind == -1) || (condicoes[t].onsets[condicoes[t].indice] < minval))
			      {
					  minind=t;
					  minval=condicoes[t].onsets[condicoes[t].indice];
			      }
			   }
			}
			if (minind ==-1) break;
			float inicio = condicoes[minind].onsets[condicoes[minind].indice++];
			int Inicio = floor(max((float)1.0, (float)(inicio+1.0) / TR) + 0.5);
			int Fim = floor(max((float)1.0, (float)(inicio+tamint) / TR));
			if ((Inicio-ultimofim) != 1) fprintf(f, "ARQUIVO INVALIDO !!!!!!!!!!!!!!!!!!!!!");
			fprintf(f, "%d-%d,%s\n", Inicio, min(numvolumes, Fim), condicoes[minind].nome);
			ultimofim=Fim;
		}
	    fclose(f);
	}

	for (int t=0; t<numcondicoes; t++) free(condicoes[t].onsets);
	free(condicoes);
}

void salvanomes(char *arquivo, Intervalo *intervalos, int tamintervalos)
{
	FILE *f;
	f=fopen(arquivo, "wt+");
	if (f!=NULL)
	{
		for (int t=0;t<tamintervalos;t++)
			fprintf(f, "%s\n", intervalos[t].nome);
		fclose(f);
	}
}

void salvanomes(char *arquivo, char **condicoes, int tamconds)
{
	FILE *f;
	f=fopen(arquivo, "wt+");
	if (f!=NULL)
	{
		for (int t=0;t<tamconds;t++)
			fprintf(f, "%s\n", condicoes[t]);
		fclose(f);
	}
}

void deslocaarray(int *classes, int tamclasses, int desloc)
{
	for (int t=tamclasses-1;t>=0;t--)
	{
		if (t-desloc>=0) classes[t]=classes[t-desloc];
		else classes[t]=classes[0];
	}
}

void filtraclasses(int *classes, int tamclasses, int classeescolhida)
{
	for (int t=0;t<tamclasses;t++)
		if (classes[t]==classeescolhida) classes[t]=1;
		else classes[t]=2;
}

int pegaindicesclasse(int *classes, int classeescolhida, int *indices, int tamindices, int **indicesclasse)
{
	int i;
	int numindices;
	numindices=0;
	int *temp;
	for (int t=0;t<tamindices;t++)
	{
		if (classes[indices[t]-1]==classeescolhida) numindices++;
	}
	temp = NULL;
	temp = (int *) malloc(numindices * sizeof(int));
	i = 0;
	for (int t=0;t<tamindices;t++)
	{
		if (classes[indices[t]-1]==classeescolhida) 
			temp[i++]=indices[t];
	}
	*indicesclasse=temp;
	return numindices;
}

int juntaindices(int *indices1, int tamindices1, int *indices2, int tamindices2, int **indices)
{
	int tamindices = tamindices1 + tamindices2;
	int *temp = (int *) malloc(tamindices * sizeof(int));
	for (int t=0;t<tamindices1;t++)	temp[t]=indices1[t];
	for (int t=tamindices1;t<tamindices;t++) temp[t]=indices2[t-tamindices1];
	*indices=temp;
	return tamindices;
}