#include "fMRIProcessing.h"
#include "correlationUtils.h"
#include <iomanip>

int correlation(char *file1, char *file2, int run1, int run2, vector<correlationConsolidation> &globalList, char *saida, float tthresh = 4.0, int demean = 1, double thresh = 0.3, int precision = 2, int absvalue = 0)
{
	vector<vector<correlationShadow>> list;
	volume4D<float> input_volume1, input_volume2;
	string input_name1(file1);
	string input_name2(file2);
	read_volume4D(input_volume1, input_name1);
	read_volume4D(input_volume2, input_name2);

	if (tthresh > 0)
	{
		input_volume1.threshold(tthresh, input_volume1.max() + 1, inclusive);
		input_volume2.threshold(tthresh, input_volume2.max() + 1, inclusive);
	}
	if (input_volume1.maxx() != input_volume2.maxx() || input_volume1.maxy() != input_volume2.maxy() || input_volume1.maxz() != input_volume2.maxz())
	{
		cerr << "Error: Mismatch in image dimensions" << endl;
		return 1;
	}

	int semmascara = 1;
	volume<float> mask;
	{
		mask = input_volume1[0];
		mask = 1;
	}


	if (demean) {
		for (int t1 = 0; t1 <= input_volume1.maxt(); t1++)
			input_volume1[t1] -= input_volume1[t1].mean(mask);
		for (int t2 = 0; t2 <= input_volume2.maxt(); t2++)
			input_volume2[t2] -= input_volume2[t2].mean(mask);
	}


	list.resize(input_volume1.maxt() + 1);
	for (int t1 = 0; t1 <= input_volume1.maxt(); t1++)
	{
		double ss1 = sqrt(input_volume1[t1].sumsquares(mask));
		for (int t2 = 0; t2 <= input_volume2.maxt(); t2++)
		{
			double ss2 = sqrt(input_volume2[t2].sumsquares(mask));
			double score = 0;
			for (int k = 0; k <= input_volume1.maxz(); k++)
				for (int j = 0; j <= input_volume1.maxy(); j++)
					for (int i = 0; i <= input_volume1.maxx(); i++)
						if (semmascara || mask(i, j, k)>0)
							score += (double)input_volume1(i, j, k, t1)*(double)input_volume2(i, j, k, t2);
			if (absvalue)
				score = fabs(score);
			score /= (ss1*ss2);
			if (score > thresh)
			{
				correlationShadow c;
				c.run = run2;
				c.volume = t2;
				c.correlation = score;
				list[t1].push_back(c);
			}
		}
	}

	// write actual results
	remove(saida);
	ofstream ofs;
	ofs.open(saida, std::ofstream::out | std::ofstream::app);
	for (int t1 = 0; t1 < list.size(); t1++)
	{
		if (list[t1].size() > 0)
		{
			sort(list[t1].begin(), list[t1].end(), maiorQue());
			ofs << setw(3) << t1 + 1 << " : ";
			for (int t2 = 0; t2 < list[t1].size(); t2++)
				ofs << std::setw(3) << list[t1][t2].volume + 1 << " " << setiosflags(ios::fixed) << setprecision(precision) << list[t1][t2].correlation << "; ";
			ofs << endl;
		}
	}
	ofs.close();

	// consolidate results so far
	if (run1 == 1)
	{
		if (run2 == 2)
		{
			for (int t1 = 0; t1 < list.size(); t1++)
			{
				if (list[t1].size() > 0)
				{
					correlationConsolidation c;
					c.volume = t1;
					c.run = run1;
					for (int t2 = 0; t2 < list[t1].size(); t2++)
						c.list.push_back(list[t1][t2]);
					globalList.push_back(c);
				}
			}
		}
		else
		{
			for (int t1 = 0; t1 < globalList.size(); t1++)
			{
				int volume = globalList[t1].volume;
				if (list[volume].size() > 0)
				{
					for (int t2 = 0; t2 < list[t1].size(); t2++)
						globalList[t1].list.push_back(list[t1][t2]);
				}
				else globalList.erase(globalList.begin() + t1);
			}
		}

	}
	return 0;
}

int correlation2(char *file1, char *sinal, char *saida)
{
	Matrix signal = read_ascii_matrix(string(sinal));
	volume4D<float> input_volume;
	volume<float> output_volume;
	read_volume4D(input_volume, string(file1));

	output_volume.reinitialize(input_volume.xsize(), input_volume.ysize(), input_volume.zsize());
	output_volume.copyproperties(input_volume[0]);

	volume<float> mask;
	mask = input_volume[0];

	for (int t1 = 0; t1 <= input_volume.maxt(); t1++)
		input_volume[t1] -= input_volume[t1].mean(mask);

	double mean = 0;

	for (int t = 1; t <= signal.Nrows(); t++)
		mean += signal(t, 1);

	mean /= signal.Nrows();

	double ss2 = 0;
	for (int t = 1; t <= signal.Nrows(); t++)
	{
		signal(t, 1) -= mean;
		ss2 += MISCMATHS::pow(signal(t, 1), 2);
	}
	ss2 = sqrt(ss2);

	double score = 0;
	float demeaned[500];
	for (int k = 0; k <= input_volume.maxz(); k++)
		for (int j = 0; j <= input_volume.maxy(); j++)
			for (int i = 0; i <= input_volume.maxx(); i++)
			{
				if (mask(i, j, k) != 0)
				{
					double ss1 = 0, mean1 = 0;
					for (int t1 = 0; t1 <= input_volume.maxt(); t1++)
					{
						mean1 += input_volume(i, j, k, t1) / input_volume.maxt();
						demeaned[t1] = input_volume(i, j, k, t1);
					}

					score = 0;
					for (int t1 = 0; t1 <= input_volume.maxt(); t1++)
					{
						demeaned[t1] -= mean1;
						ss1 += MISCMATHS::pow(demeaned[t1], 2);
						score += demeaned[t1] * signal(t1 + 1, 1);
					}
					ss1 = sqrt(ss1);
					output_volume(i, j, k) = score / (ss1*ss2);
				}
				else output_volume(i, j, k) = 0;
			}

	save_volume(output_volume, string(saida));
	return 0;
}

void generateActivationMaps(char *dirStudy, char *subject, char *feedback, float threshold = 4.0, int precision = 2)
{
	char fname[1024];
	char dirRun1[1024], dirRun2[1024];
	char dirSuj[1024];
	char dirMelodic[1024];
	char file1[1024], file2[1024];
	vector<correlationConsolidation> globalList;
	float globalThreshold = 0.51;

	sprintf(dirSuj, "%s\\MNI\\%s", dirStudy, subject);
	sprintf(dirMelodic, "%s\\Melodic", dirSuj);
	for (int i = 1; i < 5; i++)
	{
		sprintf(dirRun1, "RUN%.2dNFB0%c", i, feedback[i]);
		sprintf(file1, "%s\\sl%s.ica\\melodic_IC.nii.gz", dirMelodic, dirRun1);
		for (int j = i + 1; j < 5; j++)
		{
			char saida[1024];

			sprintf(dirRun2, "RUN%.2dNFB0%c", j, feedback[j]);
			sprintf(saida, "%s\\correlation_%s_%s.txt", dirMelodic, dirRun1, dirRun2);
			sprintf(file2, "%s\\sl%s.ica\\melodic_IC.nii.gz", dirMelodic, dirRun2);
			correlation(file1, file2, i, j, globalList, saida, threshold);
		}
	}

	vector<correlationConsolidation> filteredList;
	for (int t1 = 0; t1 < globalList.size(); t1++)
	{
		char runCheck[5] = "1000";
		int count = 0;
		for (int t2 = 0; t2 < globalList[t1].list.size(); t2++)
		{
			if (globalList[t1].list[t2].correlation > globalThreshold)
			{
				count++;
				if (globalList[t1].list[t2].run == 2) runCheck[1] = '2';
				else if (globalList[t1].list[t2].run == 3) runCheck[2] = '3';
				else if (globalList[t1].list[t2].run == 4) runCheck[3] = '4';
			}
		}

		if ((count > 0) && (strcmp("1234", runCheck) == 0))
		{
			correlationConsolidation c;
			c.volume = globalList[t1].volume;
			c.run = globalList[t1].run;
			for (int t2 = 0; t2 < globalList[t1].list.size(); t2++)
				if (globalList[t1].list[t2].correlation > globalThreshold)
					c.list.push_back(globalList[t1].list[t2]);
			filteredList.push_back(c);
		}
	}

	// write the results
	char saida[BUFF_SIZE];
	sprintf(saida, "%s\\consolidation.txt", dirMelodic);
	remove(saida);
	ofstream ofs;
	ofs.open(saida, std::ofstream::out | std::ofstream::app);
	for (int t1 = 0; t1 < filteredList.size(); t1++)
	{
		ofs << setw(3) << filteredList[t1].volume + 1 << " : ";
		for (int t2 = 0; t2 < filteredList[t1].list.size(); t2++)
			ofs << std::setw(3) << filteredList[t1].list[t2].run << " - " << std::setw(3) << filteredList[t1].list[t2].volume + 1 << " " << setiosflags(ios::fixed) << setprecision(precision) << filteredList[t1].list[t2].correlation << "; ";
		ofs << endl;
	}
	ofs.close();

	// generate the maps
	volume4D<float>runs[4];
	for (int i = 1; i < 5; i++)
	{
		sprintf(dirRun1, "RUN%.2dNFB0%c", i, feedback[i]);
		sprintf(file1, "%s\\sl%s.ica\\melodic_IC.nii.gz", dirMelodic, dirRun1);
		read_volume4D(runs[i - 1], string(file1));
	}

	volume4D<float>maps;
	maps.reinitialize(runs[0].xsize(), runs[0].xsize(), runs[0].xsize(), filteredList.size());
	maps.copyproperties(runs[0][0]);

	sprintf(saida, "%s\\consolidation.nii", dirMelodic);
	remove(saida);
	for (int t1 = 0; t1 < filteredList.size(); t1++)
	{
		int count = filteredList[t1].list.size();
		volume<float> vol;
		vol = runs[0][filteredList[t1].volume];
		for (int t2 = 0; t2 < count; t2++)
			vol += runs[filteredList[t1].list[t2].run - 1][filteredList[t1].list[t2].volume];
		vol /= count;

		//vol.threshold(threshold);
		maps[t1] = vol;
	}
	save_volume4D(maps, string(saida));
}
