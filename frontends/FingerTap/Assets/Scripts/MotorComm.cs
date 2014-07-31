using System;
using System.Collections;
namespace AssemblyCSharp
{

public class MotorComm : FriendEngineComm {
		public double firstRoiMean = 0, secondRoiMean = 0;

		protected override void writePluginInfo()
		{
			mainThread.writeSocket ("PLUGIN");
			mainThread.writeSocket ("libMotor"); //libMotor    libconnectivity  libBrainDecoding
			mainThread.writeSocket ("no");
			mainThread.writeSocket ("processMotorROI");
			mainThread.writeSocket ("initializeMotorProcessing");
			mainThread.writeSocket ("finalizeMotorProcessing");
			mainThread.writeSocket ("no");
			mainThread.writeSocket ("no");
		}

		protected override void handleGetFeedBack() 
		{
			String response; 
			if (actualState == 1) {
				responseThread.setupSocket(HostData, Port);
				responseThread.writeSocket("SESSION");
				responseThread.writeSocket(sessionID);
				FeedbackFailed = 0;
				actualState = 2;
			}
			
			if (actualState == 2) {
				response = responseThread.readSocket();
				if (response == "OK") 
				{
					responseThread.writeSocket ("TEST");
					responseThread.writeSocket (actualVolume.ToString());
					actualState = 3;
				} 
				else if (response != "") 
				{
					responseThread.closeSocket();
					actualState = 0;
					FeedbackFailed = 1;
					lastGraphResponse = response;
				}
			}
			
			if (actualState == 3) {
				lastFeedBackClass = responseThread.readSocket();
				if (lastFeedBackClass != "") 
				{
					actualState = 4;
				}
				else if ((lastFeedBackClass == "") && (!responseThread.Connected()))
				{
					responseThread.closeSocket();
					FeedbackFailed = 1;
					actualState = 0;
				}
			}
			
			if (actualState == 4) {
				lastFeedBack = responseThread.readSocket();
				if (lastFeedBack != "") 
				{
					actualState = 5;
					firstRoiMean = double.Parse(lastFeedBack);
				}
				else if ((lastFeedBack == "") && (!responseThread.Connected()))
				{
					responseThread.closeSocket();
					actualState = 0;
					FeedbackFailed = 1;
					firstRoiMean = 0;
				}
			}
			
			if (actualState == 5) {
				lastFeedBack = responseThread.readSocket();
				if (lastFeedBack != "") 
				{
					actualState = 6;
					if (!double.TryParse(lastFeedBack, out secondRoiMean))
					{
						// no second roi mean. Closing
						responseThread.closeSocket();
						actualState = 0;
						secondRoiMean = 0;
					}
				}
				else if ((lastFeedBack == "") && (!responseThread.Connected()))
				{
					responseThread.closeSocket();
					actualState = 0;
					FeedbackFailed = 1;
					secondRoiMean = 0;
				}
			}

			if (actualState == 6) {
				response = responseThread.readSocket();
				if ((response != "") || (!responseThread.Connected()))
				{
					responseThread.closeSocket();
					actualState = 0;
				}
			}
			
			if (actualState == 0) operation = 0;
		}

}
}
