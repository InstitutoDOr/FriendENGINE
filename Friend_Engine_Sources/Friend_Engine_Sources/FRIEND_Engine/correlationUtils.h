#include <vector>

using namespace std;

typedef struct
{
	int volume, run;
	double correlation;

} correlationShadow;

typedef struct
{
	int volume, run;
	vector <correlationShadow>list;
} correlationConsolidation;

struct menorQue
{
	inline bool operator() (const correlationShadow& struct1, const correlationShadow& struct2)
	{
		return (struct1.correlation < struct2.correlation);
	}
};

struct maiorQue
{
	inline bool operator() (const correlationShadow& struct1, const correlationShadow& struct2)
	{
		return (struct1.correlation > struct2.correlation);
	}
};
