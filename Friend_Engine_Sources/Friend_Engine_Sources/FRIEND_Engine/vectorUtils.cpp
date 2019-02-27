#include "vectorUtils.h"

void readVector(char *filename, float_vect &vetor)
{
	ifstream entrada;
	entrada.open(filename, std::ifstream::in);

	while (entrada.good())
	{
		double valor;
		entrada >> valor;
		vetor.push_back(valor);
	}
	entrada.close();
}

void saveVector(char *filename, float_vect &vetor)
{
	ofstream saida;
	saida.open(filename, std::ofstream::out | std::ofstream::trunc);

	for (int t = 0; t < vetor.size(); t++)
		saida << vetor[t] << "\n";
	saida.close();
}

float_vect detrendVector(float_vect &unfiltered, int w, int degree)
{
	float_vect filtered = sg_smooth(unfiltered, w, degree);
	for (int t = 0; t < filtered.size(); t++) filtered[t] = unfiltered[t] - filtered[t];
	return filtered;
}

