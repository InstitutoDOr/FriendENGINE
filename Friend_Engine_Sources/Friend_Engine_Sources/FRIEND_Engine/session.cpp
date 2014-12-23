//
//  session.cpp
//  
//
//  Created by IDOR on 16/04/14.
//
//

#include "session.h"
#include "parser.h"
#include <string>
#include <string.h>
#include <stdlib.h>
#include <stdio.h> 

void Session::processGraphMessage(const char *msg)
{
   int count;
   char **parts;
   GraphParams params;
   
   parser(msg, count, parts, ';');
   if (count == 9)
   {
      int index = atoi(parts[1]);
      
      params.rotX = atof(parts[2]);
      params.rotY = atof(parts[3]);
      params.rotZ = atof(parts[4]);
      
      params.transX = atof(parts[5]);
      params.transY = atof(parts[6]);
      params.transZ = atof(parts[7]);
      
      params.rmsValue = atof(parts[8]);
      
	  int newItem = 0;

	  if (index > graphParamsList.size()) newItem = 1;

	  graphParamsList.resize(index);
      graphParamsList[index-1] = params;

	  if (newItem)
	     graphParamsList[index - 1].additionalParams[0] = 0;
   };
   
   freeparser(count, parts);
}

void Session::processFeedback(int index, float returnedClass, float returnedPercentage)
{
	int newItem = 0;

	if (index > feedbacksList.size()) newItem = 1;
    feedbacksList.resize(index);
    feedbacksList[index-1].classNumber = (int) returnedClass;
    feedbacksList[index-1].percentage = returnedPercentage;

	if (newItem)
		feedbacksList[index - 1].additionalResponse[0] = 0;
}

void Session::processAdditionalGraphInfo(int index, const char *info)
{
	if (strlen(info) > 0)
	{
		graphParamsList.resize(index);
		strcpy(graphParamsList[index - 1].additionalParams, info);
		strcat(graphParamsList[index - 1].additionalParams, "\n");
	}
}

void Session::processAdditionalFeedBackInfo(int index, const char *info)
{
	if (strlen(info) > 0)
	{
		feedbacksList.resize(index);
		strcpy(feedbacksList[index - 1].additionalResponse, info);
		strcat(feedbacksList[index - 1].additionalResponse, "\n");
	}
}

void Session::getGraphResponse(int index, char *msg)
{
   if (graphParamsList.size() < index) sprintf(msg, "GRAPHPARS\n");
   else
   {
	   sprintf(msg, "GRAPHPARS;%d;%f;%f;%f;%f;%f;%f;%f\n%s", index, graphParamsList[index - 1].rotX, graphParamsList[index - 1].rotY, graphParamsList[index - 1].rotZ, graphParamsList[index - 1].transX, graphParamsList[index - 1].transY, graphParamsList[index - 1].transZ, graphParamsList[index - 1].rmsValue, graphParamsList[index - 1].additionalParams);
   }
}

void Session::getFeedbackResponse(int index, char *msg)
{
   if (feedbacksList.size() < index) sprintf(msg, "0\n0\n");
   else sprintf(msg, "%d\n%f\n%s", feedbacksList[index - 1].classNumber, feedbacksList[index - 1].percentage, feedbacksList[index - 1].additionalResponse);
}

void Session::getCommandResponse(string command, char *response)
{
   std::map<string, int>::iterator it;
   it = responseList.find(command);
   if (it != responseList.end())
   {
      if (responseList[command] == 0) sprintf(response, "NYT\n");
      else sprintf(response, "OK\n");
   }
   else sprintf(response, "NOK\n");
}

void Session::setVDBPointer(void *vdb)
{
   vdbPtr = vdb;
   
}

void *Session::getVDBPointer()
{
   return vdbPtr;
}

void Session::setCommandResponse(string command, int status)
{
   responseList[command] = status;
}
