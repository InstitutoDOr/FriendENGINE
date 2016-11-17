#include <Windows.h>

class Timer
{
private:
	double start, end;

public:
	Timer() { startCounter(); }

	double elapsed() {
		end = GetTickCount();
		return (end-start) / 1000.;
	}

	void startCounter() { start = GetTickCount(); }
};
