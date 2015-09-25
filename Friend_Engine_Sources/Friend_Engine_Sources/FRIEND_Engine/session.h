//
//  session.h
//  
//
//  Created by IDOR on 16/04/14.
//
//

#ifndef ____session__
#define ____session__

#include <iostream>
#include <vector>
#include <string>
#include <map>

using namespace std;
#define stringInfoSize 500

typedef struct
{
   
   double rotX, rotY, rotZ;
   double transX, transY, transZ;
   double rmsValue;
   char additionalParams[stringInfoSize];
   
} GraphParams;

typedef struct
{
   
   int classNumber;
   float percentage;
   char additionalResponse[stringInfoSize];
   
} FeedBackResponse;



class Session
{
   char ip[200];
   int port;
   int runSize;
   int terminationStatus;
   
   vector<GraphParams> graphParamsList;
   vector<FeedBackResponse> feedbacksList;
   map<string, int> responseList;
   void *vdbPtr;
   
public:
   int getFeedbackResponses;
   void processGraphMessage(const char *msg);
   void processFeedback(int index, float returnedClass, float returnedPercentage);

   void processAdditionalGraphInfo(int index, const char *info);
   void processAdditionalFeedBackInfo(int index, const char *info);
   
   void getGraphResponse(int index, char *msg);
   void getFeedbackResponse(int index, char *msg);
   void getCommandResponse(string command, char *response);
   void setCommandResponse(string command, int status);
   void setVDBPointer(void *vdb);
   void *getVDBPointer();
   void terminate();
   int getTerminateState();
   
   Session() { getFeedbackResponses = 0; terminationStatus = 0; };
};

#endif /* defined(____session__) */
