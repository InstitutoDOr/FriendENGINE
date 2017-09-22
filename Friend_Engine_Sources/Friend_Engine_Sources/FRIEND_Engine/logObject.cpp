#include "logObject.h"
#include "stdarg.h"

// create temporary log file
void LogObject::initializeTempFile(char *filename)
{
	tempFileName = filename;
	outputStream.open(filename, fstream::out | fstream::trunc);
	outputStream << buffer.str().c_str();
}

// create log file
void LogObject::initializeLogFile(char *filename)
{
	fileName = filename;
	if (outputStream.is_open()) outputStream.close();
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
	va_list args;
	va_start(args, format);
	vsprintf(auxBuffer, format, args);
	if (fileName != "")
		outputStream << auxBuffer;
	else
	{
		buffer << auxBuffer;
		if (tempFileName != "")
			outputStream << auxBuffer;
	}

	if (inScreen)
		fprintf(stderr, "%s", auxBuffer);
	va_end(args);
}
