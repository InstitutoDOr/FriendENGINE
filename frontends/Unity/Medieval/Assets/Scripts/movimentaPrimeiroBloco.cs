using UnityEngine;
using System.Collections;
using System;

public class movimentaPrimeiroBloco : MonoBehaviour {
	public GUIText texto; 
	private AssemblyCSharp.FriendEngineComm comm;
	private String HostData = "127.0.0.1";
	private Int32 Port = 5678;
	private float TR = 2;
	private float tryInterval;
	private float timeCounter;
	private int baselineBlockLength = 10;
	private GameObject actualRock, toDisappear;
	private String actualPath;
	private double meanActivation=0;
	private int[] initialsIndex = new int[] {11, 36, 61, 86};
	private int[] triggerIndex = new int[] {24, 49, 74, 99};
	private String[] paths = new String[] {"firstPath", "secondPath", "thirdPath", "fourthPath"};
	private String[] rocks = new String[] {"firstRock", "secondRock", "thirdRock", "fourthRock"};


	void calculateMeanActivation ()
	{
	   meanActivation = 0;
	   int firstIndex = comm.volume() - 8;
	   if (firstIndex < 0) firstIndex = 0;
	   for (int i = 0; i< 8; i++) meanActivation += comm.feedbackValues[firstIndex + i];
	   meanActivation /= 8;
    }

	void handleNextPath()
	{
		texto.text = "REST "; 
		toDisappear = actualRock;
		iTween.MoveTo (actualRock, iTween.Hash("y", 500, "time", 2));
		
		actualPath = paths[comm.phase()];
		actualRock = GameObject.Find (rocks[comm.phase()]);
		startPath();
	}

	IEnumerator stateManager()
	{
		if (comm.state() == 0)
			texto.text = "REST "; 

		comm.coreCommunication();

		// starting the first path
		if (comm.state() == 15)
		{
			texto.text = "";
				
			comm.setPhase(0);
			actualPath = paths[comm.phase()];
			actualRock = GameObject.Find (rocks[comm.phase()]);
			startPath();
			comm.setState(16);
		}

		if (comm.state() == 35)
		{
			if (comm.volume() >= triggerIndex[comm.phase()]) comm.setState(36);
		}

		if (comm.state() == 36) 
		{
			comm.setPhase(1);
			texto.text = "REST "; 
			handleNextPath();
			comm.setState(37);
		}


		if (comm.state() == 45)
		{
			if (comm.volume() >= triggerIndex[comm.phase()]) comm.setState(46);
		}

		if (comm.state() == 46) 
		{
			comm.setPhase(2);
			toDisappear.renderer.enabled = false;
			texto.text = "REST " ; 
			handleNextPath();
			comm.setState(47);
		}

		if (comm.state() == 55)
		{
			if (comm.volume() >= triggerIndex[comm.phase()]) comm.setState(56);
		}

		if (comm.state() == 56) 
		{
			comm.setPhase(3);
			toDisappear.renderer.enabled = false;
			handleNextPath();
			comm.setState(57);
		}

		if (comm.state() == 65)
		{
			if (comm.volume() >= triggerIndex[comm.phase()]) comm.setState(66);
		}

		if (comm.state() == 66) 
		{
			toDisappear.renderer.enabled = false;
		};

		comm.handleGraphInformation();
		if (comm.state() == 100) 
		{
			StopCoroutine("stateManager");
			texto.text = "END";
		};

		yield return null;
	}

	void startPath()
	{
		iTween.MoveTo(gameObject, iTween.Hash("path", iTweenPath.GetPath(actualPath), "orientToPath", true, "time", (baselineBlockLength-1) * TR, 
		                                      "easeType", iTween.EaseType.easeInOutSine, "oncomplete", "BaselineBlockEnd"));
	}
	// Use this for initialization
	void Start () {
		tryInterval = TR / 2;
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
		if (comm.state() < 30) comm.setState(35);
		else if (comm.state() < 40) comm.setState(45);
		else if (comm.state() < 50) comm.setState(55);
		else if (comm.state() < 60) comm.setState(65);
	}
}
