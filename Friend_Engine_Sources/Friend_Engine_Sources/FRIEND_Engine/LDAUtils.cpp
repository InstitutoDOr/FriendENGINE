#include "LDAUtils.h"
#include <algorithm>

string generateLDALine(strings &vocabulary, strings &list)
{
	string line;
	line = to_string(list.size()) + " ";

	for (int t = 0; t < list.size(); t++)
		line += to_string(stringIndex(vocabulary, list[t])) + ":1 ";
	return line;
}


typedef struct
{
	int index;
	double value;
} doublePoint;

bool greaterPoint(doublePoint i, doublePoint j) { return (i.value > j.value); }

void saveBestTopics(char *outputDir, strings &vocabulary)
{
	char betaFile[1024], outputFile[1024];
	sprintf(betaFile, "%s%cfinal.beta", outputDir, PATHSEPCHAR);

	ifstream entrada;
	strings vetor;
	vector<doublePoint> betas;

	// reading data
	entrada.open(betaFile, fstream::in);
	int topic = 1;

	sprintf(outputFile, "%s.txt", outputDir);
	FILE *saida = fopen(outputFile, "wt+");
	while (entrada.good())
	{
		string line;
		getline(entrada, line);
		line = trimString(line);
		if (line.size() > 0)
		{
			vetor.clear();
			betas.clear();
			splitLine(line, vetor, ' ');
			for (int t = 0; t < vetor.size(); t++)
			{
				doublePoint p;
				p.index = t;
				p.value = stod(vetor[t]);
				betas.push_back(p);
			}
			std::sort(betas.begin(), betas.end(), greaterPoint);

			fprintf(saida, "Topico %d : ", topic);
			for (int t = 0; t < 5; t++)
				fprintf(saida, "%s - (%f) / ", vocabulary[betas[t].index].c_str(), betas[t].value);

			fprintf(saida, "\n");
			topic++;
		}
	}
	entrada.close();
	fclose(saida);
}
