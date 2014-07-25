// fslexec.cpp : Defines the entry point for the console application.
//

#include "fslio.h"
#include "intervalos.h"
#include <string>
#include <string.h>
#include <stdlib.h>
#include <stdio.h> 
#include <math.h> 
#include "dirfuncs.h"
#include <fstream>
#include <list>
#include "confusionmatrix.h"

extern "C" __declspec(dllimport) int _stdcall fazsusan(char *vol, char *saida, char *mascara, float FWHM);
extern "C" __declspec(dllimport) int _stdcall Filtra4D(char *Vol4D, int *indices, int tamindices, char *Saida);
extern "C" __declspec(dllimport) int _stdcall NormalizaPorBaseLine(char *Vol4D, char *arqintervalos, char *condicaobasal, char *Saida);
extern "C" __declspec(dllimport) int _stdcall fsl_glm(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall film_gls(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall fslroi(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall flirt(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall mcflirt(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall fslmerge(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall fslsplit(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall applywarp(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall applyxfm4D(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall bet(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall cluster(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall connectedcomp(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall convert_xfm(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall convertwarp(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall createlut(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall fnirt(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall fsl_boxplot(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall fsl_histogram(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall fsl_tsplot(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall fslcc(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall fslchfiletype(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall fslcomplex(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall fslcorrecthd(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall fslcpgeom(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall fslfft(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall fslhd(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall fslmaths(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall fslmeants(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall fslorient(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall fslstats(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall fslswapdim(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall fslcreatehd(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall fugue(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall fast(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall melodic(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall overlay(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall slicer(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall slicetimer(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall smoothest(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall susan(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall feat_model(char *CmdLn);


extern "C" __declspec(dllimport) int _stdcall convolui(int tamanho, float Mult, float *valores, float *saida);
extern "C" __declspec(dllimport) int _stdcall retornaparams(char *matrixfile, char *volnam, float *rx, float *ry, float *rz, float *tx, float *ty, float *tz);

extern "C" __declspec(dllimport) int _stdcall CalculaAtivacao(int Inicio, int ini, int fim, int tamb, char *prefixo, char *Basal, char *Mascara, char *Saida);
extern "C" __declspec(dllimport) int _stdcall CalculaAtivacao4D(int Inicio, int ini, int fim, int tamb, char *vol4D, char *Basal, char *Mascara, char *Saida);
extern "C" __declspec(dllimport) int _stdcall axial(char *innam, char *outnam);
extern "C" __declspec(dllimport) int _stdcall SetaOrigem(char *Volume, char *Saida, short dx, short dy, short dz);
extern "C" __declspec(dllimport) int _stdcall PegaOrigem(char *Volume, short *dx, short *dy, short *dz);
extern "C" __declspec(dllimport) int _stdcall Realinha(char *Anatomico, char *sAnatomico);
extern "C" __declspec(dllimport) int _stdcall AnatomicoResampleado(char *Anatomico, char *sAnatomico, float dx, float dy, float dz, float TR, int nn);
extern "C" __declspec(dllimport) int _stdcall StandarizaVolume(char *base, char *mask, char *saida, int NN);
extern "C" __declspec(dllimport) void _stdcall SalvaArquivoSVM(char *Volume, char *mascara, char *svm, float Minimo, int inicio, int fim, int *classes);
extern "C" __declspec(dllimport) void _stdcall SalvaArquivoSVMIndices(char *Volume, char *mascara, char *svm, float Minimo, int *indices, int tamindices, int *classes);
extern "C" __declspec(dllimport) void _stdcall SalvaArquivoSVMIndices2(char *Volume, char *mascara, char *svm, float Minimo, int *indices, int tamindices, int *classes);

extern "C" __declspec(dllexport) int _stdcall fslswapdim_rt(char *CmdLn, FSLIO *src);
extern "C" __declspec(dllexport) FSLIO * _stdcall fslioopen(char *arquivo);
extern "C" __declspec(dllexport) void _stdcall fslioclose(FSLIO *src);


extern "C" __declspec(dllimport) int _stdcall SVMtrain(char *CmdLn);
extern "C" __declspec(dllimport) int _stdcall SVMpredict(char *CmdLn);

//extern "C" __declspec(dllimport) int _stdcall RodaKMeans(char *nomealg, char *Volume, char *mascara, float minimo, int *indices, int tamindices, int centros, char *saida);
extern "C" __declspec(dllimport) void _stdcall demean(char *arq, char *saida);
extern "C" __declspec(dllimport) void _stdcall stdone(char *arq, char *mascara, char *saida);

extern "C" __declspec(dllimport) void _stdcall NormalizacaoZ(char *arquivobase, char *mascara, char *saida);
extern "C" __declspec(dllimport) void _stdcall NormalizacaoZ2(char *arquivobase, char *mascara, char *saida);
extern "C" __declspec(dllimport) void _stdcall CalculaMedias(char *arquivobase, char*arqintervalos);
extern "C" __declspec(dllimport) int _stdcall PegaTamanhoMascara(char *Mascara, float Minimo);

#define TAMBASELINE 3
#define GAUSSFHWM 5
#define normsvm 1

void FiltroTemporal(char *Arquivo, char *ArquivoSaida, float TR, int LowP, int HighP, float HpCutOff = 100)
{
   float hp_sigma_vol, lp_sigma_vol;
   float hp_sigma_sec, lp_sigma_sec;
   char cmd[5000];

   if (LowP || HighP)
   {
      hp_sigma_vol = -1;
      if (HighP) 
	  {
         hp_sigma_sec = (float) HpCutOff / 2;
         hp_sigma_vol = (float) hp_sigma_sec / TR;
	  }

      lp_sigma_vol = -1;

      if (LowP)
	  {
         lp_sigma_sec = (float) 2.8;
         lp_sigma_vol = (float) lp_sigma_sec / TR;
	  }
	  sprintf(cmd, "fslmaths %s -bptf %f %f %s", Arquivo, hp_sigma_vol, lp_sigma_vol, ArquivoSaida);
      fslmaths(cmd);
   }
}

using namespace std;
using std::list;

// contrastes
// timeseries

void performanceindices(char *arquivo, int *classes, char **condicoes, int numcondicoes, float **performances, int *indices, char *arquivocm)
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

	if (f != NULL)
	{
		int *matrix;
	    createConfusionMatrix(numcondicoes, &matrix);
	    zeraMatrix(matrix, numcondicoes);
		while (fgets(linha, 255, f))
		{
		   sscanf(linha, "%f %f", &predicao, &valor);
		   signalResult(classes[indices[i]-1], predicao, numcondicoes, matrix);
		   if (classes[indices[i]-1]==predicao)
		   {
			   acertos[0]++;
			   acertos[classes[indices[i]-1]]++;
		   }
		   totais[0]++;
		   totais[classes[indices[i]-1]]++;
		   i++;
		}
		for(int t=0;t<=numcondicoes;t++) (*performances)[t] = ((float)acertos[t]*100.0/(float)totais[t]);
		fclose(f);
	    saveMatrix(matrix, numcondicoes, condicoes, arquivocm);
	}
}

void preprocessamento(char *arquivo, char *arqintervalo, float FHWM, int highpass, int ativacao, int kernelbox, char *condicaobasal, char **saida)
{
          char DirT[255];
		  char Arq[255];
		  char SaidaMcFlirt[255];
		  char SaidaMcFlirtAct[255];
		  char SaidaFiltroT[255];
		  char SaidaBet[255];
		  char MascaraBet[255];
		  char SaidaGauss[255];
		  char SaidaAtivacao[255];
		  char SaidaKernelBox[255];
		  char Media[255];
		  char prefixobasal[255];
		  int tamanhovolume;
		  
		  char cmd[2550];

		  char Arquivo[255];
		  strcpy(Arquivo, arquivo);
		  extractfilepath(Arquivo, DirT);
		  if (strlen(DirT)>0) strcat(DirT, "\\");
		  extractfilename(Arquivo, Arq);
          sprintf(SaidaMcFlirt, "%s%s%s", DirT, "MC_", Arq);
          sprintf(SaidaMcFlirtAct, "%s%s%s", DirT, "MC_Act_", Arq);
          sprintf(SaidaFiltroT, "%s%s%s", DirT, "Tempo_", Arq);
          sprintf(SaidaBet, "%s%s", DirT, "bet.nii");
          sprintf(MascaraBet, "%s%s", DirT, "mascara.nii");
          sprintf(SaidaGauss, "%s%s", DirT, "bet_gauss.nii");
          sprintf(SaidaAtivacao, "%s%s", DirT, "bet_gauss_act.nii");
          sprintf(SaidaKernelBox, "%s%s", DirT, "bet_gauss_box.nii");
          sprintf(Media, "%s%s", DirT, "media.nii");
          sprintf(prefixobasal, "%s%s", DirT, "media");

		  FSLIO *arq = fslioopen(arquivo);
		  tamanhovolume=arq->niftiptr->nt;
		  fslioclose(arq);

		  if (!fileexists(SaidaMcFlirt))
		  {
             sprintf(cmd, "mcflirt -in %s -refvol 1 -out %s -plots -mats -rmsabs", Arquivo, SaidaMcFlirt);
		     mcflirt(cmd);
		  }

		  CalculaMedias(SaidaMcFlirt, arqintervalo);

		  if (!fileexists(Media))
		  {
		     sprintf(cmd, "fslmaths %s -Tmean %s -odt float", SaidaMcFlirt, Media);
		     fslmaths(cmd);
		  }
   
		  if (!fileexists(MascaraBet))
		  {
		     sprintf(cmd, "bet %s %s -f 0.3 -n -m", Media, MascaraBet);
		     bet(cmd);
		  }

		  NormalizaPorBaseLine(SaidaMcFlirt, arqintervalo, condicaobasal, SaidaMcFlirtAct);

		  if (highpass > 0)
		  {
		     if (!fileexists(SaidaFiltroT))  
			 {
	            sprintf(cmd, "fslmaths %s -bptf %f %f %s -odt float", SaidaMcFlirtAct, 2, highpass, SaidaFiltroT);
                fslmaths(cmd);
			 }
		  }
		  else strcpy(SaidaFiltroT, SaidaMcFlirtAct);

		  if (!fileexists(SaidaBet))
		  {
		     sprintf(cmd, "fslmaths %s -mas %s %s -odt float", SaidaFiltroT, MascaraBet, SaidaBet);
		     fslmaths(cmd);
		  }

		  if (!fileexists(SaidaGauss))
		  {
		     sprintf(cmd, "fslmaths %s -kernel gauss %f -fmean %s -odt float", SaidaBet, (float) (FHWM/2.3548), SaidaGauss);
		     fslmaths(cmd);
		  }

		  strcpy(*saida, SaidaGauss);
	      if (ativacao == 1) 
		  {
		     CalculaAtivacao4D(1, 1, tamanhovolume, TAMBASELINE, *saida, NULL, NULL, SaidaAtivacao);
		     strcpy(*saida, SaidaAtivacao);
		  }

		  if (kernelbox == 1) 
		  {
			  char cmd[500];
		      sprintf(cmd, "fslmaths %s -kernel boxv 3 -fmeanu %s -odt float", *saida, SaidaKernelBox);
		      fslmaths(cmd);
			  strcpy(*saida, SaidaKernelBox);
		  }
}

void gravamat(char *arquivo, Intervalo *intervalos, int tamint, char **condicoes, int tamconds, int inicio=0, int fim=0)
{
	int tamfloat;
	float **desenhos, **convolucoes;
	desenhos = (float **) malloc(tamconds * sizeof(float *));
	convolucoes = (float **) malloc(tamconds * sizeof(float *));
	if (inicio == 0) inicio = 1;
	if (fim == 0) fim = intervalos[tamint-1].fim;
	tamfloat = fim-inicio+1;
	for (int i = 0; i<tamconds; i++)
	{
		desenhos[i] = NULL;
		convolucoes[i] = (float *) malloc(tamfloat * sizeof(float));
		pegadesenho(desenhos[i], intervalos, tamint, condicoes[i], inicio, fim);
		convolui(tamfloat, 0.72, desenhos[i], convolucoes[i]);
	}
	FILE *f;
    f = fopen(arquivo, "wt+");
	char linha[255];

	for (int i = 0; i < tamfloat; i++)
	{
		for (int j = 0; j < tamconds; j++)
		{
		   if (j==0) sprintf(linha, "%1.2f", convolucoes[j][i]);
		   else sprintf(linha, "%s\t%1.2f", linha, convolucoes[j][i]);
		}
		fprintf(f, "%s\n", linha);
	}
	fclose(f);

	for (int i = 0; i<tamconds; i++)
	{
		free(desenhos[i]);
		free(convolucoes[i]);
	}
	free(desenhos);
	free(convolucoes);

}

void inserecondicao(char *linha, int condicao)
{
	if (strlen(linha) == 0) sprintf(linha, "%d", condicao);
	else sprintf(linha, "%s\t%d", linha, condicao);
}

void fazarquivocontraste(char *arquivo, char **condicoes, int tamconds, int baseline)
{
	FILE *f;
	f = fopen(arquivo, "wt+");
	char linha[255];
	for (int i = 0; i<tamconds; i++)
	{
		linha[0] = '\0';
		if (i != baseline)
		{
		   for (int j = 0; j < tamconds; j++)
		   {
			   if (j==i) inserecondicao(linha, 1);
			   else if (j==baseline) inserecondicao(linha, -1);
			   else inserecondicao(linha, 0);
		   }
		   fprintf(f, "%s\n", linha);
		}
	}
	fclose(f);
}

/*
# Slice timing correction
# 0 : None
# 1 : Regular up (0, 1, 2, 3, ...)
# 2 : Regular down
# 3 : Use slice order file
# 4 : Use slice timings file
# 5 : Interleaved (0, 2, 4 ... 1, 3, 5 ... )
*/

void SliceTimeCorrection(char *arquivo, char *arquivosaida, int tipo, float TR, int direcao = 0, char *ordem = NULL)
{
	char opt[255];
	if (tipo == 0) return;

	switch (tipo)
	{ 
	case 1 :
		strcpy(opt, "");
		break;
	case 2 :
		strcpy(opt, "--down");
		break;
	case 3 :
		if (ordem != NULL) sprintf(opt, "--ocustom=%s", ordem);
		break;
	case 4 :
		if (ordem != NULL) sprintf(opt, "--tcustom=%s", ordem);
		break;
	case 5 :
		strcpy(opt, "--odd");
		break;
	}

	char cmd[500];

    if (direcao > 0) sprintf(cmd, "slicetimer -i %s --out=%s -r %f %s -d %d", arquivo, arquivosaida, TR, opt, direcao);
	else sprintf(cmd, "slicetimer -i %s --out=%s -r %f %s", arquivo, arquivosaida, TR, opt);
	slicetimer(cmd);
}

int main(int argc, char* argv[])
{

	char cmd[255];
//	strcpy(cmd, "fsl_tsplot -i C:\\projetos\\GITO\\FIG\\prefiltered_func_data_mcf.par -t \"MCFLIRT estimated rotations (radians)\" -u 1 --start=1 --finish=3 -a x,y,z -w 640 -h 144 -o C:\\projetos\\GITO\\FIG\\rot.png");
//	fsl_tsplot(cmd);
/*	strcpy(cmd, "fsl_tsplot -i C:\\projetos\\GITO\\FIG\\prefiltered_func_data_mcf.par -t \"MCFLIRT estimated translations (mm)\"   -u 1 --start=4 --finish=6 -a x,y,z -w 640 -h 144 -o C:\\projetos\\GITO\\FIG\\trans.png");
	fsl_tsplot(cmd);
    strcpy(cmd, "fsl_tsplot -i C:\\projetos\\GITO\\FIG\\prefiltered_func_data_mcf_abs.rms,C:\\projetos\\GITO\\FIG\\prefiltered_func_data_mcf_rel.rms -t \"MCFLIRT estimated mean displacement (mm)\" -u 1 -w 640 -h 144 -a absolute,relative -o C:\\projetos\\GITO\\FIG\\disp.png");
	fsl_tsplot(cmd);
*/
	if (argc > 1)
	{
       char cmdln[50000];
	   char Cmd[500];
	   strcpy(cmdln, argv[1]);
	   for (int i = 2; i<argc;i++) 
	   {
		   strcat(cmdln, " \"");
		   strcat(cmdln, argv[i]);
		   strcat(cmdln, "\"");
	   }

	   char *cmd = _strupr(argv[1]);
	   if (strcmp(cmd, "FSL_GLM") == 0) fsl_glm(cmdln);
	   else if (strcmp(cmd, "FILM_GLS") == 0) 
	   {
		   film_gls(cmdln);
	   }
	   else if (strcmp(cmd, "FLIRT") == 0) flirt(cmdln);
	   else if (strcmp(cmd, "MCFLIRT") == 0) mcflirt(cmdln);
	   else if (strcmp(cmd, "FSLMERGE") == 0)
	   {
          if (argc<5) fslmerge(cmdln);
		  else if (strchr(argv[4], '*')!=NULL)
		  {
             flist  list = { 0, 0, NULL };
		     char nome[255];
		     strcpy(nome, argv[4]);

	         size_t origsize = strlen(nome) + 1;
	         size_t convertedChars = 0;
	         wchar_t root[255];
	         char rootc[255];
	         extractfilepath(nome, rootc);
	         if (strlen(rootc) > 0) strcat(rootc, "\\");

	         mbstowcs_s(&convertedChars, root, origsize, nome, _TRUNCATE);
	         search(list, root);
	         sort(list);
	         char testef[50000];
	         strcpy(testef, argv[1]);
			 strcat(testef, " ");
			 strcat(testef, argv[2]);
			 strcat(testef, " ");
			 strcat(testef, argv[3]);
             for (int i = 0; i < list.num_entries; i++)
	         {
                char tempfile[255];     
   	            char filename[255];
		        origsize=wcslen(list.files[i].cFileName) + 1;
		        wcstombs_s(&convertedChars, filename, origsize, list.files[i].cFileName, _TRUNCATE);
		        sprintf(tempfile, " %s%s", rootc, filename);
		        strcat(testef, tempfile);
	         }
		     fslmerge(testef);
		  }
		  else fslmerge(cmdln);
	   }
	   else if (strcmp(cmd, "FSLSPLIT") == 0) fslsplit(cmdln);
	   else if (strcmp(cmd, "FSLROI") == 0) fslroi(cmdln);
	   else if (strcmp(cmd, "APPLYWARP") == 0) applywarp(cmdln);
	   else if (strcmp(cmd, "APPLYXFM4D") == 0) applyxfm4D(cmdln);
	   else if (strcmp(cmd, "BET") == 0) bet(cmdln);
	   else if (strcmp(cmd, "CLUSTER") == 0) cluster(cmdln);
	   else if (strcmp(cmd, "CONNECTEDCOMP") == 0) connectedcomp(cmdln);
	   else if (strcmp(cmd, "CONVERT_XFM") == 0) convert_xfm(cmdln);
	   else if (strcmp(cmd, "CONVERTWARP") == 0) convertwarp(cmdln);
	   else if (strcmp(cmd, "CREATELUT") == 0) createlut(cmdln);
	   else if (strcmp(cmd, "FNIRT") == 0) fnirt(cmdln);
	   else if (strcmp(cmd, "FSL_BOXPLOT") == 0) fsl_boxplot(cmdln);
	   else if (strcmp(cmd, "FSL_HISTOGRAM") == 0) fsl_histogram(cmdln);
	   else if (strcmp(cmd, "FSL_TSPLOT") == 0) fsl_tsplot(cmdln);
	   else if (strcmp(cmd, "FSLCC") == 0) fslcc(cmdln);
	   else if (strcmp(cmd, "FSLCHFILETYPE") == 0) fslchfiletype(cmdln);
	   else if (strcmp(cmd, "FSLCOMPLEX") == 0) fslcomplex(cmdln);
	   else if (strcmp(cmd, "FSLCORRECTHD") == 0) fslcorrecthd(cmdln);
	   else if (strcmp(cmd, "FSLCPGEOM") == 0) fslcpgeom(cmdln);
	   else if (strcmp(cmd, "FSLFFT") == 0) fslfft(cmdln);
	   else if (strcmp(cmd, "FSLHD") == 0) fslhd(cmdln);
	   else if (strcmp(cmd, "FEAT_MODEL") == 0) 
	   {
		   feat_model(cmdln);
	   }
	   else if (strcmp(cmd, "FSLMATHS") == 0) fslmaths(cmdln);
	   else if (strcmp(cmd, "FSLMEANTS") == 0) fslmeants(cmdln);
	   else if (strcmp(cmd, "FSLORIENT") == 0) fslorient(cmdln);
	   else if (strcmp(cmd, "FSLSTATS") == 0) fslstats(cmdln);
	   else if (strcmp(cmd, "FSLSWAPDIM") == 0) fslswapdim(cmdln);
	   else if (strcmp(cmd, "FSLCREATEHD") == 0) fslcreatehd(cmdln);
	   else if (strcmp(cmd, "FUGUE") == 0) fugue(cmdln);
	   else if (strcmp(cmd, "FAST") == 0) fast(cmdln);
	   else if (strcmp(cmd, "MELODIC") == 0) melodic(cmdln);
	   else if (strcmp(cmd, "OVERLAY") == 0) overlay(cmdln);
	   else if (strcmp(cmd, "SLICER") == 0) slicer(cmdln);
	   else if (strcmp(cmd, "SLICETIMER") == 0) slicetimer(cmdln);
	   else if (strcmp(cmd, "SMOOTHEST") == 0) smoothest(cmdln);
	   else if (strcmp(cmd, "SVMTRAIN") == 0) SVMtrain(cmdln);
	   else if (strcmp(cmd, "SVMPREDICT") == 0) SVMpredict(cmdln);
	   else if (strcmp(cmd, "SUSAN") == 0) 
	   {
			   susan(cmdln);
	   }

	   else if ((strcmp(cmd, "AXIAL") == 0) && (argc > 3)) axial(argv[2], argv[3]);
	   else if ((strcmp(cmd, "TAMANHOMASCARA") == 0) && (argc > 3)) 
	   {
		   int tamanho = PegaTamanhoMascara(argv[2], (float) atof(argv[3]));
		   short x, y, z, v;
//           FSLIO *hdr = FslOpen(argv[2], "r");
//           FslGetDim(hdr, &x, &y, &z, &v);
//		   printf("Size = %d Percentage=%3.2f", tamanho, 100*(tamanho/(x*y*z)));
//           FslClose(hdr);
	   }
	   else if ((strcmp(cmd, "SETAORIGEM") == 0) && (argc > 6)) SetaOrigem(argv[2], argv[3], atoi(argv[4]), atoi(argv[5]), atoi(argv[6]));
	   else if ((strcmp(cmd, "REALINHA") == 0) && (argc > 3)) Realinha(argv[2], argv[3]);
	   else if ((strcmp(cmd, "ANATOMICORESAMPLEADO") == 0) && (argc > 7)) AnatomicoResampleado(argv[2], argv[3], (float) atof(argv[4]), (float) atof(argv[5]), (float) atof(argv[6]), (float) atof(argv[7]), 0);
	   else if ((strcmp(cmd, "STANDARIZAVOLUME") == 0) && (argc > 4)) 
	   {
		   StandarizaVolume(argv[2], argv[3], argv[4], 1);
		   return 0;
	   }
	   else if ((strcmp(cmd, "RUNMELODIC") == 0) && (argc > 0)) 
	   {
		  char DirSaida[255];
		  extractfilepath(argv[2], DirSaida);

		  if (strlen(DirSaida) > 0) strcat(DirSaida, "\\");
		  sprintf(Cmd, "melodic -i %s -o %s -v --nobet --bgthreshold=1 --tr=2 -d 10 --mmthresh=0.5 --nomask --report", argv[2], DirSaida);
		  melodic(Cmd);
	   }
	   else if ((strcmp(cmd, "RUNGLM") == 0) && (argc > 3)) 
	   {
          char SaidaGauss[255];
		  char DirSaida[255];
		  char ArqMatriz[255];
		  char ArqContraste[255];
		  char MascaraBet[255];
		  char MascaraBin[255];
		  char Saida[255];
		  char Saidat[255];
		  char Saidaz[255];
		  char cope[255];
		  char treino[255];

		  extractfilepath(argv[2], DirSaida);
		  if (strlen(DirSaida) > 0) strcat(DirSaida, "\\");
		  strcpy(SaidaGauss, argv[2]);

		  sprintf(ArqMatriz, "%s%s", DirSaida, "matriz.txt");
		  sprintf(ArqContraste, "%s%s", DirSaida, "contraste.con");
		  sprintf(Saida, "%s%s", DirSaida, "glmsaida");
		  sprintf(Saidat, "%s%s", DirSaida, "glmsaidat");
		  sprintf(Saidaz, "%s%s", DirSaida, "glmsaidaz");
          sprintf(MascaraBet, "%s%s", DirSaida, "mascara.nii");
          sprintf(MascaraBin, "%s%s", DirSaida, "mascarabin.nii");
          sprintf(treino, "%s%s", DirSaida, "4Dtreino.nii");
          sprintf(cope, "%s%s", DirSaida, "cope.nii");

		  //printf("passo1");
		  Intervalo *intervalos;
	      char **condicoes;
		  int inicio=0, fim;

	      int tamintervalos, tamcondicoes;
	      tamintervalos=LeIntervalo(argv[3], &intervalos);
	      tamcondicoes=pegalistacondicoes(intervalos, tamintervalos, &condicoes);
		  fim=intervalos[tamintervalos-1].fim / 2;
		  //printf("passo2");
		  if (argc > 4) inicio=atoi(argv[4]);
		  if (argc > 5) fim=atoi(argv[5]);
		  if (inicio==0) inicio=1;
		  if (fim==0) fim=intervalos[tamintervalos-1].fim;
		  gravamat(ArqMatriz, intervalos, tamintervalos, condicoes, tamcondicoes, inicio, fim);
	      if (!fileexists(ArqContraste)) fazarquivocontraste(ArqContraste, condicoes, tamcondicoes, 0);

		  //printf("passo3");
		  if (!fileexists(MascaraBet)) 
		  {
             sprintf(Cmd, "fslmaths %s -Tmean %s", argv[2], MascaraBet);
		     fslmaths(Cmd);
		  }

		  if (!fileexists(MascaraBin))
		  {
             sprintf(Cmd, "fslmaths %s -bin %s -odt char", MascaraBet, MascaraBin);
		     fslmaths(Cmd);
		  }

		  if (!fileexists(treino))
		  {
		     sprintf(cmd, "fslroi %s %s %d %d", SaidaGauss, treino, inicio-1, fim-inicio+1);
		     fslroi(cmd);
		  }

		  //printf("passo4");
		  sprintf(Cmd, "fsl_glm -i %s -d %s -c %s -m %s -o %s --out_t=%s --out_z=%s --out_cope=%s", treino, ArqMatriz, ArqContraste, MascaraBin, Saida, Saidat, Saidaz, cope);
		  fsl_glm(Cmd);
		  
		  //printf("passo5");
		  free(intervalos);
		  //printf("passo6");
		  desalocacondicoes(condicoes, tamcondicoes);
		  //printf("passo7");
	   }
	   else if ((strcmp(cmd, "UNEMASCARAS") == 0) && (argc > 0)) 
	   {
		  sprintf(Cmd, "fslmaths %s -Tmax %s", argv[2], argv[3]);
		  fslmaths(Cmd);

		  sprintf(Cmd, "fslmaths %s -thr %s %s", argv[3], argv[4], argv[3]);
		  fslmaths(Cmd);
	   }
	   else if ((strcmp(cmd, "SVM") == 0) && (argc >= 7)) 
	   {
		   int *classes, *classes2;

		   Intervalo *intervalos;
	       char **condicoes;
		   float *performances;

		   char DirSaida[255];
		   char arqtreino[255];
		   char arqativacao[255];
		   char arqteste[255];
		   char modelo[255];
		   char mascara[255];
		   char arqresultado[255];
		   char arqresultadotreino[255];
		   char arqperformance[255];
		   char treino[255];
		   char ztreino[255];
		   char teste[255];
		   char zteste[255];
		   int meio, total;

		   extractfilepath(argv[2], DirSaida);
		   if (strlen(DirSaida) > 0) strcat(DirSaida, "\\");

		   sprintf(arqtreino, "%s%s", DirSaida, "treino.txt");
		   sprintf(arqteste, "%s%s", DirSaida, "teste.txt");
		   sprintf(modelo, "%s%s", DirSaida, "modelo.svm");
		   sprintf(arqresultadotreino, "%s%s", DirSaida, "resultadotreino.txt");
		   sprintf(arqresultado, "%s%s", DirSaida, "resultado.txt");
		   sprintf(arqperformance, "%s%s", DirSaida, "performance.txt");
		   sprintf(mascara, "%s%s", DirSaida, "mascara.nii");
		   sprintf(treino, "%s%s", DirSaida, "4Dtreinosvm.nii");
		   sprintf(ztreino, "%s%s", DirSaida, "4Dtreinosvmz.nii");
		   sprintf(teste, "%s%s", DirSaida, "4Dtestesvm.nii");
		   sprintf(zteste, "%s%s", DirSaida, "4Dtestesvmz.nii");
		   strcpy(arqativacao, argv[2]);

		   int tamintervalos, tamcondicoes;
	       tamintervalos=LeIntervalo(argv[3], &intervalos);
	       tamcondicoes=pegalistacondicoes(intervalos, tamintervalos, &condicoes);
		   retornaarrayclasses(intervalos, tamintervalos, condicoes, tamcondicoes, &classes);

		   total = intervalos[tamintervalos-1].fim;
		   meio = total / 2;
//		   deslocaarray(classes, total, 1);
		   classes2 = (int *) malloc(total * sizeof(int));

		   //for (int t=0;t<tamcondicoes;t++)
		   {
               int *indices;
		       int tamindices = PegaIndices(intervalos, tamintervalos, atoi(argv[5]), 1, meio, &indices);
			   for(int j=0;j< total; j++) classes2[j]=classes[j];
			   if (1)
			   {
			       int *indices1=NULL;
			       int *indices2=NULL;
				   int classeatual=indicestring(condicoes, tamcondicoes, "Positivo")+1;
				   int tamindices1;
				   tamindices1=pegaindicesclasse(classes2, classeatual, indices, tamindices, &indices1);

				   classeatual=indicestring(condicoes, tamcondicoes, "Negativo")+1;
				   int tamindices2=pegaindicesclasse(classes2, classeatual, indices, tamindices, &indices2);
				   free(indices);
				   indices=NULL;
				   tamindices=juntaindices(indices1, tamindices1, indices2, tamindices2, &indices);
				   free(indices1);
				   free(indices2);
			   }
			   //filtraclasses(classes2, total, t+1);

			   char arquivo[255];
			   sprintf(arquivo, "%s%s.txt", DirSaida, "treinoi");
			   FILE *f=fopen(arquivo, "wt+");
			   for (int h=0;h<tamindices;h++) fprintf(f, "%d\t%d\n", indices[h], classes2[indices[h]-1]);
			   fclose(f);
			   Filtra4D(arqativacao, indices, tamindices, treino);
			   if (normsvm) NormalizacaoZ(treino, mascara, ztreino);
			   else strcpy(ztreino, treino);

			   SalvaArquivoSVMIndices2(ztreino, argv[4], arqtreino, 0, indices, tamindices, classes2);
			   free(indices);
 		       sprintf(Cmd, "svmtrain -t 0 %s %s", arqtreino, modelo);
			   SVMtrain(Cmd);

			   sprintf(Cmd, "svmpredict %s %s %s", arqtreino, modelo, arqresultadotreino);
		       SVMpredict(Cmd);
		   }

           int *indices;
           int tamindices = PegaIndices(intervalos, tamintervalos, atoi(argv[5]), meio+1, intervalos[tamintervalos-1].fim, &indices);

		   //for (int t=0;t<tamcondicoes;t++)
		   {
			   for(int j=0;j< total; j++) classes2[j]=classes[j];
			   //filtraclasses(classes2, total, t+1);
			   if (1)
			   {
				   int classeatual=indicestring(condicoes, tamcondicoes, "Positivo")+1;
				   int *indices1;
				   int tamindices1=pegaindicesclasse(classes2, classeatual, indices, tamindices, &indices1);

				   classeatual=indicestring(condicoes, tamcondicoes, "Negativo")+1;
				   int *indices2;
				   int tamindices2=pegaindicesclasse(classes2, classeatual, indices, tamindices, &indices2);
				   free(indices);
				   tamindices=juntaindices(indices1, tamindices1, indices2, tamindices2, &indices);
				   free(indices1);
				   free(indices2);
			   }

			   char arquivo[255];
			   sprintf(arquivo, "%s%s.txt", DirSaida, "testei");
			   FILE *f=fopen(arquivo, "wt+");
			   for (int h=0;h<tamindices;h++) fprintf(f, "%d\t%d\n", indices[h], classes2[indices[h]-1]);
			   fclose(f);

			   Filtra4D(arqativacao, indices, tamindices, teste);
			   if (normsvm) NormalizacaoZ2(teste, mascara, zteste);
			   else strcpy(zteste, teste);

			   SalvaArquivoSVMIndices2(zteste, argv[4], arqteste, 0, indices, tamindices, classes2);
 		       sprintf(Cmd, "svmpredict %s %s %s", arqteste, modelo, arqresultado);
		       SVMpredict(Cmd);
		   }

		   performanceindices(arqresultado, classes, condicoes, tamcondicoes, &performances, indices, argv[6]);
    	   free(indices);
		   FILE *f;
		   f=fopen(arqperformance, "wt+");
		   if (f != NULL)
		   {
			   fprintf(f, "Performance Global\t\t: %.2f\n", performances[0]);
			   for (int t=0;t<tamcondicoes; t++) fprintf(f, "Performance da classe %s\t: %.2f\n", condicoes[t], performances[t+1]);
			   fclose(f);
		   }
		   free(intervalos);
		   free(classes);
		   free(classes2);
		   free(performances);
		   desalocacondicoes(condicoes, tamcondicoes);
	   }
	   else if ((strcmp(cmd, "ANA2NII") == 0) && (argc >= 6)) 
	   {
		   FSLIO * base=fslioopen(argv[2]);
		   int volumes=atoi(argv[3]);
		   
		   char cmd[255];
		   char infile[255];
		   char outfile[255];
		   char prefixoin[255];
		   char prefixoout[255];

		   strcpy(prefixoin, argv[4]);
		   strcpy(prefixoout, argv[5]);
		   for (int t=1; t<=volumes; t++)
		   {
			   sprintf(infile, "%s%.5d.img", prefixoin, t);
			   sprintf(outfile, "%s%.5d.nii", prefixoout, t);
			   sprintf(cmd, "fslswapdim %s x -y z %s", infile, outfile);
			   fslswapdim_rt(cmd, base);
		   }
		   fslioclose(base);
	   }
	   else if ((strcmp(cmd, "MOTIONCORRECT") == 0) && (argc >= 6)) 
	   {
		   int volumes=atoi(argv[3]);
		   
		   char cmd[500];
		   char infile[255];
		   char outfile[255];
		   char prefixoin[255];
		   char prefixoout[255];

		   strcpy(prefixoin, argv[4]);
		   strcpy(prefixoout, argv[5]);
		   for (int t=1; t<=volumes; t++)
		   {
			   sprintf(infile, "%s%.5d.nii", prefixoin, t);
			   sprintf(outfile, "%s%.5d.nii", prefixoout, t);
			   sprintf(cmd, "mcflirt -in %s -reffile %s -out %s", infile, argv[2], outfile);
			   mcflirt(cmd);
		   }
	   }
	   else if ((strcmp(cmd, "GAUSS") == 0) && (argc >= 5)) 
	   {
		   int volumes=atoi(argv[2]);
		   
		   char cmd[500];
		   char infile[255];
		   char outfile[255];
		   char prefixoin[255];
		   char prefixoout[255];

		   strcpy(prefixoin, argv[3]);
		   strcpy(prefixoout, argv[4]);
		   for (int t=1; t<=volumes; t++)
		   {
			   sprintf(infile, "%s%.5d.nii", prefixoin, t);
			   sprintf(outfile, "%s%.5d.nii", prefixoout, t);
		       sprintf(cmd, "fslmaths %s -kernel gauss %f -fmean %s", infile, (float) (GAUSSFHWM/2.3548), outfile);
		       fslmaths(cmd);
		   }
	   }
	   else if ((strcmp(cmd, "GAUSSF") == 0) && (argc >= 5)) 
	   {
          char cmd[500];
		  sprintf(cmd, "fslmaths %s -kernel gauss %f -fmean %s", argv[2], (float) (atof(argv[4])/2.3548), argv[3]);
		  fslmaths(cmd);
	   }
	   else if ((strcmp(cmd, "ATIVACAO4D") == 0) && (argc >= 7)) 
	   {
		   char *basal;
		   char *mascara;

		   if (argc >= 8) basal = argv[7];
		   else basal = NULL;
		   if (argc >= 9) mascara = argv[8];
		   else mascara = NULL;

		   CalculaAtivacao4D(1, atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), argv[2], basal, mascara, argv[6]);
           
		   char cmd[500];
		   sprintf(cmd, "fslmaths %s -kernel boxv 3 -fmeanu %s", argv[6], argv[6]);
		   fslmaths(cmd);
	   }
	   else if ((strcmp(cmd, "KMEANS") == 0) && (argc >= 5)) 
	   {
		  int *classes;
		  int *map;
		  Intervalo *intervalos;
	      char **condicoes;
		  float *performances;
	      int tamintervalos, tamcondicoes;
		  char DirSaida[255];
		  char assignment[255];
		  char arqperformance[255];

		  extractfilepath(argv[2], DirSaida);
		  if (strlen(DirSaida) > 0) strcat(DirSaida, "\\");

		  sprintf(assignment, "%s%s", DirSaida, "assignment.txt");
		  sprintf(arqperformance, "%s%s", DirSaida, "perf.txt");

          tamintervalos=LeIntervalo(argv[3], &intervalos);
          tamcondicoes=pegalistacondicoes(intervalos, tamintervalos, &condicoes);
   	      retornaarrayclasses(intervalos, tamintervalos, condicoes, tamcondicoes, &classes);

		  int tam = intervalos[tamintervalos-1].fim;
		  int  centros = tam / 30;
	      int *indices=(int *) malloc(tam * sizeof(int));
	      for (int t=0;t<tam;t++) indices[t] = t+1;

//	      RodaKMeans("hybrid", argv[2], argv[4], 0, indices, tam, (int) centros, assignment);
		  identificaclasses(assignment, classes, tamcondicoes, centros, &map, indices, tam);
		  FILE *f;
		  f=fopen(arqperformance, "wt+");
		  if (f != NULL)
		  {
			  fprintf(f, "Mapeamentos KMeans :\n");
			  char status[20];
			  for (int t=0;t<tam;t++)
			  {
				  if (map[t] < 0) strcpy(status, "transicao");
				  else strcpy(status, "bloco");
				  fprintf(f, "Indice : %d Classe : %s Status : %s\n", indices[t], condicoes[abs(map[t])-1], status);
			  }
			  fclose(f);
		  }

	      free(indices);
		  free(intervalos);
		  free(classes);
		  free(map);
		  desalocacondicoes(condicoes, tamcondicoes);
	   }
	   else if ((strcmp(cmd, "PREPROC") == 0) && (argc >= 9)) 
	   {
		   char *saida;

		   saida = (char *) malloc(500*sizeof(char));
		   preprocessamento(argv[2], argv[3], atof(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), argv[8], &saida);
		   free(saida);
	   }
	   else if ((strcmp(cmd, "DEMEAN") == 0) && (argc >= 4)) 
	   {
		   demean(argv[2], argv[3]);
	   }
	   else if ((strcmp(cmd, "STDONE") == 0) && (argc >= 5)) 
	   {
		   stdone(argv[2], argv[3], argv[4]);
	   }
	   else if ((strcmp(cmd, "OUPROC") == 0) && (argc >= 5)) 
	   {
		   FSLIO *arq;
		   int total, meio;

		   char parte1[500];
		   char parte2[500];
		   char DirSaida[500];
		   char cmd[1000];

		   extractfilepath(argv[2], DirSaida);
		   if (strlen(DirSaida) > 0) strcat(DirSaida, "\\");
		   arq = fslioopen(argv[2]);
		   total=arq->niftiptr->nt;
		   meio = total / 2;
		   fslioclose(arq);

		   sprintf(parte1, "%s%s", DirSaida, "parte1.nii");
		   sprintf(parte2, "%s%s", DirSaida, "parte2.nii");

		   sprintf(cmd, "fslroi %s %s %d %d", argv[2], parte1, 0, meio);
		   fslroi(cmd);
		   sprintf(cmd, "fslroi %s %s %d %d", argv[2], parte2, meio, total-meio);
		   fslroi(cmd);

		   demean(parte1, parte1);
		   stdone(parte1, argv[3], parte1);
		   CalculaAtivacao4D(1, 1, meio, TAMBASELINE, parte1, NULL, NULL, parte1);

		   demean(parte2, parte2);
		   stdone(parte2, argv[3], parte2);
		   CalculaAtivacao4D(1, 1, total-meio, TAMBASELINE, parte2, NULL, NULL, parte2);

		   sprintf(cmd, "fslmerge -t %s %s %s", argv[4], parte1, parte2);
		   fslmerge(cmd);
	   }
	   else if ((strcmp(cmd, "OOPROC") == 0) && (argc >= 5)) 
	   {
		   FSLIO *arq;
		   int total, meio;

		   char parte1[500];
		   char parte2[500];
		   char DirSaida[500];
		   char cmd[1000];

		   extractfilepath(argv[2], DirSaida);
		   if (strlen(DirSaida) > 0) strcat(DirSaida, "\\");
		   arq = fslioopen(argv[2]);
		   total=arq->niftiptr->nt;
		   meio = total / 2;
		   fslioclose(arq);

		   sprintf(parte1, "%s%s", DirSaida, "parte1.nii");
		   sprintf(parte2, "%s%s", DirSaida, "parte2.nii");

		   sprintf(cmd, "fslroi %s %s %d %d", argv[2], parte1, 0, meio);
		   fslroi(cmd);
		   sprintf(cmd, "fslroi %s %s %d %d", argv[2], parte2, meio, total-meio);
		   fslroi(cmd);

		   demean(parte1, parte1);
		   stdone(parte1, argv[3], parte1);
		   CalculaAtivacao4D(1, 1, meio, TAMBASELINE, parte1, NULL, NULL, parte1);

		   demean(parte2, parte2);
		   stdone(parte2, argv[3], parte2);
		   CalculaAtivacao4D(1, 1, total-meio, TAMBASELINE, parte2, NULL, NULL, parte2);

		   sprintf(cmd, "fslmerge -t %s %s %s", argv[4], parte1, parte2);
		   fslmerge(cmd);
	   }
	   else if ((strcmp(cmd, "CALCULAMEDIAS") == 0) && (argc >= 4)) 
	   {
		   CalculaMedias(argv[2], argv[3]);
	   }
	   else if ((strcmp(cmd, "OPROC") == 0) && (argc >= 5)) 
	   {
   		  //printf("passo7.5");
		   NormalizacaoZ(argv[2], argv[4], argv[3]);
   		  //printf("passo7.6");
	   }
	   else if ((strcmp(cmd, "SEGMENTAGM") == 0) && (argc >= 3)) 
	   {
	      char cmd[500];
		  char dirprog[255];
		  char DirSaida[255];
		  char betfile[255];

		  extractfilepath(argv[0], dirprog);
		  if (strlen(dirprog) > 0) strcat(dirprog, "\\");

		  extractfilepath(argv[2], DirSaida);
		  if (strlen(DirSaida) > 0) strcat(DirSaida, "\\");

		  sprintf(betfile, "%s%s", DirSaida, "bethires.nii");
		  if (!fileexists(betfile))
		  {
		     sprintf(cmd, "bet %s %s -f 0.3 -n -m", argv[2], betfile);
		     bet(cmd);
		  }

          sprintf(cmd,"fast -A %stissuepriors\\avg152T1_csf %stissuepriors\\avg152T1_gray %stissuepriors\\avg152T1_white %s", dirprog, dirprog, dirprog, betfile);
		  fast(cmd);
	   }
	   else if ((strcmp(cmd, "FAZSUSAN") == 0) && (argc >= 5)) 
	   {
		   fazsusan(argv[2], argv[3], argv[4], atof(argv[5]));
	   }

	}
  //printf("passo8");
	return 0;
}
