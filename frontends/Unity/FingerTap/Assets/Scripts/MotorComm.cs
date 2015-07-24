using System;
using System.Collections;
using System.IO;

namespace AssemblyCSharp
{

public class MotorComm : FriendEngineComm {
		public double firstRoiMean = 0, secondRoiMean = 0;
		public int blockLength = 22;

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

		public override void doSetup()
		{
			feedbackRun = 1;
			timeBaseVolumeIndex = true;
			setRunsize (296);

			addConfigurationPair ("MNIMask", "studydirhmat_spm_final.nii");
			addConfigurationPair ("MNITemplate", "studydirMNI152_T1_1mm_brain.nii.gz");
			addConfigurationPair ("Prefix", "outputdirRUN01" + Path.DirectorySeparatorChar + "DRIN-");
			addConfigurationPair ("ActivationLevel", "0.01");
			addConfigurationPair ("CurrentRunSuffix", "RUN01");
			startBlockIndexes = new int[] {01, 16, 38, 53, 75, 90, 112, 127, 149, 164, 186, 201, 223, 238, 260, 275};
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
					responseThread.writeSocket (actualQueryVolume.ToString());
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
