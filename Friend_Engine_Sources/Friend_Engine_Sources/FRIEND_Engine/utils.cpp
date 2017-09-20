#include "utils.h"
#include <time.h>

#ifdef _WIN32
#include <direct.h>
#include <sys/timeb.h>
#elif
#include <sys/time.h>
#include <sys/stat.h>
#endif

void replaceAll(string& str, const string& from, const string& to) 
{
	if (from.empty())
		return;
	string wsRet;
	wsRet.reserve(str.length());
	size_t start_pos = 0, pos;
	while ((pos = str.find(from, start_pos)) != string::npos) {
		wsRet += str.substr(start_pos, pos - start_pos);
		wsRet += to;
		pos += from.length();
		start_pos = pos;
	}
	wsRet += str.substr(start_pos);
	str.swap(wsRet); // faster than str = wsRet;
}

#ifdef WINDOWS
int getMilliCount()
{
	timeb tb;
	ftime(&tb);
	int nCount = tb.millitm + (tb.time & 0xfffff) * 1000;
	return nCount;
}

int getMilliSpan(int nTimeStart)
{
	int nSpan = getMilliCount() - nTimeStart;
	if (nSpan < 0)
		nSpan += 0x100000 * 1000;
	return nSpan;
}
#endif

double GetWallTime()
{

#ifdef WINDOWS
	return (double)getMilliCount() / 1000.0;
#else
	struct timeval time;
	if (gettimeofday(&time, NULL))
	{
		//  Handle error
		return 0;
	}
	return (double)time.tv_sec + (double)time.tv_usec * .000001;
#endif
}
