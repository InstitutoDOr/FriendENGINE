#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void createConfusionMatrix(int rank, int **matrix)
{
	*matrix = (int *) malloc(rank * rank * sizeof(int));
}

void zeraMatrix(int *matrix, int rank)
{
	int maximo = rank*rank;
	for (int i=0; i<maximo;i++) matrix[i] = 0;
}

void copiaMatrix(int rank, int *destino, int *origem)
{
	int maximo = rank*rank;
	for (int i=0; i<maximo;i++) destino[i] = origem[i];
}

void signalResult(int desej, int obt, int rank, int *matrix)
{
	matrix[(desej-1) * rank + (obt-1)]++;
}

float chance(int rank)
{
	return (float) (100.0 / (float) rank);
}

float acertos(int *matrix, int rank)
{
   int total=0, acc=0;

   for (int i=0;i<rank;i++)
	   for (int j=0;j<rank;j++)
	   {
		   if (i==j) acc = acc + matrix[i*rank + j];
		   total = total + matrix[i*rank + j];
	   };
   return ((float)acc * 100) / (float) total;
}

float Kappa(int *matrix, int rank) // Porcentagem de acerto obtida - Porcentagem da chance / 1 - Porcentagem da chance
{
   return (acertos(matrix, rank) - chance(rank)) / (100-chance(rank));
}

int exemplos(int *matrix, int rank)
{
	int maximo = rank*rank;
	int total=0;
	for (int i=0; i<maximo;i++) total+=matrix[i];
	return total;
}

float erroclasse(int classe, int *matrix, int rank)
{
   int acc=0, total=0;
   for (int i=0;i<rank;i++)
   {
	   if (i == classe) acc = matrix[i*rank+classe-1];
	   total += matrix[i*rank+classe-1];
   }
   return (float)(acc * 100) / (float)total;
}

float erroRefclasse(int classe, int *matrix, int rank)
{
   int acc=0, total=0;
   for (int i=0;i<rank;i++)
   {
      if (i == classe) acc = matrix[(classe-1)*rank+i];
	  total += matrix[(classe-1)*rank+i];
   }
   return (float)(acc * 100) / (float)total;
}

void saveMatrix(int *matrix, int rank, char **nomes, char *arquivo)
{
	FILE *f=fopen(arquivo, "wt+");
	char buffer[255];
	if (f!=NULL)
	{
		fprintf(f, "-------------------------------------------------------------------------------------------------------------------------------\n");
        fprintf(f, "|%20s|", "Nomes");
		for (int i=0;i<rank;i++) fprintf(f, "%20s|", nomes[i]);
        fprintf(f, "%20s|\n", "Acuracia");
		fprintf(f, "-------------------------------------------------------------------------------------------------------------------------------\n");

		for (int i=0;i<rank;i++) 
		{
			fprintf(f, "|%20s|", nomes[i]);
			int acc=0, total=0;
			for (int j=0;j<rank;j++)
			{
				fprintf(f, "%20d|", matrix[i*rank+j]);
				if (i==j) acc = matrix[i*rank+j];
				total += matrix[i*rank+j];

			}
			if (total == 0) sprintf(buffer, "%.2f", 0);
			else sprintf(buffer, "%.2f", (float)acc*100/(float)total);
			strcat(buffer, "%");
			fprintf(f, "%20s|\n", buffer);
		}

		fprintf(f, "-------------------------------------------------------------------------------------------------------------------------------\n");
		fprintf(f, "|%20s|", "Especificidade");

		for (int j=0;j<rank;j++)
		{
			int acc=0, total=0;
			for (int i=0;i<rank;i++)
			{
				if (i==j) acc = matrix[i*rank+j];
				total += matrix[i*rank+j];
			}
			if (total == 0) sprintf(buffer, "%.2f", 0);
			else sprintf(buffer, "%.2f", (float)acc*100/(float)total);
			strcat(buffer, "%");
			fprintf(f, "%20s|", buffer);
		}
        fprintf(f, "%20s|\n", "");
		fprintf(f, "-------------------------------------------------------------------------------------------------------------------------------\n");
		fprintf(f, "\n");

		sprintf(buffer, "%.2f", acertos(matrix, rank));
		strcat(buffer, "%");

		fprintf(f, "Valor Kappa    : %.2f\n", Kappa(matrix, rank));
		fprintf(f, "Acurácia Geral : %s\n", buffer);
		fprintf(f, "Total Exemplos : %d\n", exemplos(matrix, rank));
		fclose(f);
	}
}

