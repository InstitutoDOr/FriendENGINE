using UnityEngine;
using System.Collections;
using System;

public class movimentaPrimeiroBloco : MonoBehaviour {
	public GUIText texto; 
	private AssemblyCSharp.FriendEngineComm comm;
	private String HostData = "10.36.1.140";
	private Int32 Port = 5678;
	private float TR = 2;
	private float tryInterval;
	private float timeCounter;
	private Int32 actualState = 0, actualStepState = 0;
	private int actualVolume;
	private int baselineBlockLength = 10;
	private GameObject actualRock, toDisappear;
	private String actualPath;
	private double meanActivation=0;
	private int[] initialsIndex = new int[] {11, 36, 61, 86};
	private int[] triggerIndex = new int[] {24, 49, 74, 99};
	private int actualPhase=0;
	private String[] paths = new String[] {"firstPath", "secondPath", "thirdPath", "fourthPath"};
	private String[] rocks = new String[] {"firstRock", "secondRock", "thirdRock", "fourthRock"};


	void calculateMeanActivation ()
	{
	   meanActivation = 0;
	   int firstIndex = actualVolume - 8;
	   if (firstIndex < 0) firstIndex = 0;
	   for (int i = 0; i< 8; i++) meanActivation += comm.feedbackValues[firstIndex + i];
	   meanActivation /= 8;
    }

	IEnumerator stateManager()
	{
		// connecting and creating a new session
		if (actualState == 0) 
		{
			texto.text = "REST " + actualPhase.ToString(); 

			comm.connect();
			comm.createSession();
			actualState = 1;
			yield return null;
		}

		// waiting the session creation
		if (actualState == 1)
		{
			if (comm.stateManager() == 0)
				actualState = 2;
			yield return null;
		}

		// setting the plug-in information
		if (actualState == 2)
		{
			comm.setupPlugIn();
			actualState = 3;
			yield return null;
		}

		// waiting the plugin configuration
		if (actualState == 3)
		{
			if (comm.stateManager() == 0)
				actualState = 4;
			yield return null;
		}

		// Issuing the Preproc Command  
		if (actualState == 4)
		{
			comm.issueCommand("NBPREPROC");
			actualState = 5;
			yield return null;
		};

		// waiting the command acknowledge
		if (actualState == 5) 
		{
			if (comm.stateManager() == 0)
				actualState = 6;
			yield return null;
		};

		// initiate preprocessing query
		if (actualState == 6) 
		{
			// sending FEEDBACK command
			comm.getSessionCommandState("PREPROC");
			actualState = 7;
			yield return null;
		}

		// waiting the end of preprocessing 
		if (actualState == 7) 
		{
			if (comm.stateManager() == 0)
				actualState = 8;
			yield return null;
		}

		// initiate feedback processing
		if (actualState == 8) 
		{
			// sending FEEDBACK command
			comm.issueCommand("NBFEEDBACK");
			actualState = 9;
			yield return null;
		}

		// waiting acknowledge
		if (actualState == 9)
		{
			if (comm.stateManager() == 0)
				actualState = 10;
			yield return null;
		}

		// getting the first graph parameter
		if (actualState == 10)
		{
			actualVolume = 1;
			comm.getGraphParams(actualVolume);
			actualState = 11;
			yield return null;
		}

		// waiting the end of graph call
		if (actualState == 11)
		{
			if (comm.stateManager() == 0)
				actualState = 12;
			yield return null;
		}

		// starting the first path
		if (actualState == 12)
		{
			if (comm.lastGraphResponse == "GRAPHPARS") actualState = 10;
			else
			{
				texto.text = "";

				actualPhase = 0;
				actualPath = paths[actualPhase];
				actualRock = GameObject.Find (rocks[actualPhase]);
				actualStepState = 0;
				startPath();
				actualState = 15;
			}
			yield return null;
		}

		if (actualState == 35)
		{
			if (actualVolume >= triggerIndex[actualPhase]) actualState = 36;
		}

		if (actualState == 36) 
		{
			actualPhase = 1;
			texto.text = "REST " + actualPhase.ToString(); 
			toDisappear = actualRock;
			iTween.MoveTo (actualRock, iTween.Hash("y", 500, "time", 2));

			actualPath = paths[actualPhase];
			actualRock = GameObject.Find (rocks[actualPhase]);
			startPath();
			actualState = 37;
			yield return null;
		}


		if (actualState == 45)
		{
			if (actualVolume >= triggerIndex[actualPhase]) actualState = 46;
		}

		if (actualState == 46) 
		{
			actualPhase = 2;
			toDisappear.renderer.enabled = false;
			texto.text = "REST " + actualPhase.ToString(); 
			toDisappear = actualRock;
			iTween.MoveTo (actualRock, iTween.Hash("y", 500, "time", 2));

			actualPath = paths[actualPhase];
			actualRock = GameObject.Find (rocks[actualPhase]);
			startPath();
			actualState = 47;
			yield return null;
		}

		if (actualState == 55)
		{
			if (actualVolume >= triggerIndex[actualPhase]) actualState = 56;
		}

		if (actualState == 56) 
		{
			actualPhase = 3;
			toDisappear.renderer.enabled = false;
			texto.text = "REST " + actualPhase.ToString(); 
			toDisappear = actualRock;
			iTween.MoveTo (actualRock, iTween.Hash("y", 500, "time", 2));

			actualPath = paths[actualPhase];
			actualRock = GameObject.Find (rocks[actualPhase]);
			startPath();
			actualState = 57;
			yield return null;
		}

		if (actualState == 65)
		{
			if (actualVolume >= triggerIndex[actualPhase]) actualState = 66;
		}

		if (actualState == 66) 
		{
			toDisappear.renderer.enabled = false;
			yield return null;
		};

		if ((actualState > 11) && (actualState < 100))
		{
			if (actualStepState == 0)
			{
				comm.getGraphParams(actualVolume);
				actualStepState = 1;
				yield return null;
			}

			if (actualStepState == 1)
			{
				if (comm.stateManager() == 0)
				   actualStepState = 2;
				yield return null;
			}

			if (actualStepState == 2)
			{
				comm.getFeedBack(actualVolume);
				actualStepState = 3;
				yield return null;
			}

			if (actualStepState == 3)
			{
				if (comm.stateManager() == 0)
					actualStepState = 4;
				yield return null;
			}

			if (actualStepState == 4)
			{
			   if ((comm.lastGraphResponse != "GRAPHPARS") && (comm.lastGraphResponse != "END"))
			   {
				  double Feedback = comm.feedbackValues[actualVolume-1];
				  actualVolume++;
			   }
			   actualStepState = 0;  
				yield return null;
			}
			if (comm.lastGraphResponse == "END") 
			{
				actualState = 100;
				comm.endSession();
				StopCoroutine("stateManager");
			}
		};
	}

	void startPath()
	{
		iTween.MoveTo(gameObject, iTween.Hash("path", iTweenPath.GetPath(actualPath), "orientToPath", true, "time", (baselineBlockLength-1) * TR, 
		                                      "easeType", iTween.EaseType.easeInOutSine, "oncomplete", "BaselineBlockEnd"));
	}
	// Use this for initialization
	void Start () {
		tryInterval = TR / 2;
		actualState = 0;
		timeCounter = tryInterval;

		comm = new AssemblyCSharp.FriendEngineComm();
		comm.setupConnection(HostData, Port);
		comm.setRunsize (200);
	}


	// Update is called once per frame
	void Update () 
	{
		timeCounter -= Time.deltaTime;
		if (timeCounter <= 0) 
		{
			StartCoroutine("stateManager");
			timeCounter = tryInterval;
		};
	}

	void BaselineBlockEnd() {
		texto.text = "FINGER TAP"; 
		if (actualState < 30) actualState = 35;
		else if (actualState < 40) actualState = 45;
		else if (actualState < 50) actualState = 55;
		else if (actualState < 60) actualState = 65;
	}
}
