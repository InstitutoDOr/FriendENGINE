#include "logObject.h"
#include "stdarg.h"
#include <ctime>

// create log file
void LogObject::initializeLogFile(char *filename)
{
	fileName = filename;
	outputStream.open(filename, fstream::out | fstream::trunc);
	outputStream << buffer.str().c_str();
}

// closes the log file
void LogObject::closeLogFile()
{
	if (fileName != "")
		outputStream.close();
}

// writes the message in the log file
void LogObject::writeLog(int inScreen, const char * format, ...)
{
	char auxBuffer[3000];
	char timestamp[100];

	time_t time_now;
	time(&time_now);
	struct tm *timeinfo = localtime(&time_now);
	strftime(timestamp, 100, "%Y-%m-%d %H:%M:%S : ", timeinfo);
	va_list args, screen_args;

	va_start(args, format);
	vsprintf(auxBuffer, format, args);
	va_end(args);

	if (fileName != "")
		outputStream << timestamp << auxBuffer;
	else
		buffer << timestamp << auxBuffer;

	va_start(screen_args, format);
	if (inScreen)
		//fprintf(stderr, "%s\n", timestamp);
		vfprintf(stderr, format, screen_args);
	va_end(screen_args);
}
