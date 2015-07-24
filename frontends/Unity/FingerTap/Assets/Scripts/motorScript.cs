using UnityEngine;
using System.Collections;
using System.IO;
using System;

public class motorScript : MonoBehaviour {
	public GUIText texto; 

	private Vector3 leftStart, leftEnd;
	private Vector3 rightStart, rightEnd;

	private GameObject leftHand, rightHand;
	private GameObject leftTarget, rightTarget;

	private float leftPercentage, rightPercentage;
	private bool leftMove, rightMove;
	private double timeCounter, tryInterval;
	private double armTimeCounter, armTryInterval;
	private double handTimeCounter, handTryInterval;
	private bool leftHandClose, rightHandClose;
	private bool leftHandAnimate, rightHandAnimate;
	private float incr = 5;

	private AssemblyCSharp.MotorComm comm;


	void resetState()
	{
		leftHandAnimate = false;
		rightHandAnimate = false;
		
		leftMove = false;
		rightMove = false;

		leftPercentage = 0.2f;
		rightPercentage = 0.2f;

		moveHand(leftTarget, ref leftPercentage, leftStart, leftEnd, incr);
		moveHand(rightTarget, ref rightPercentage, rightStart, rightEnd, incr);
	}

	bool isOdd(int value)
	{
		return value % 2 != 0;
	}

	IEnumerator stateManager()
	{
		if (comm.state() == 0)
			texto.text = "RELAX"; 

		comm.coreCommunication();

		if ((isOdd(comm.actualBlock)) && (comm.isBlockStart()))
		{
			texto.text = "RELAX";
			resetState();
		}
		else if ((!isOdd(comm.actualBlock)))
		{
			if (comm.state() < 100)
			{
				if (comm.isBlockStart()) texto.text = "FINGERTAP ONE HAND ONLY";
				if ((comm.firstRoiMean != 0) || (comm.secondRoiMean != 0))
				{
					if (comm.firstRoiMean > comm.secondRoiMean) 
					{
						leftHandAnimate = true;
						rightHandAnimate = false;
						leftPercentage = (float) comm.firstRoiMean;
						
						leftMove = true;
						rightMove = false;
					}
					else
					{
						leftHandAnimate = false;
						rightHandAnimate = true;
						rightPercentage = (float) comm.secondRoiMean;
						
						leftMove = false;
						rightMove = true;
					}
					if (!rightMove) 
					{
						rightPercentage = 0.2f;
						moveHand(rightTarget, ref rightPercentage, rightStart, rightEnd, incr);
					}
					
					if (!leftMove) 
					{
						leftPercentage = 0.2f;
						moveHand(leftTarget, ref leftPercentage, leftStart, leftEnd, incr);
					}
				}
			}
		}
		if (comm.state () >= 100) 
		{
			StopCoroutine("stateManager");
			texto.text = "END";
			
			resetState ();
		}
		yield return null;
	}

	// Use this for initialization
	void Start () {
		leftStart = new Vector3 (999.752f, 0.96305f, 1000.09f );
		leftEnd = new Vector3 (999.767f, 1.29232f, 1000.349f );

		rightStart = new Vector3 (1000.303f, 0.96305f, 1000.09f );
		rightEnd = new Vector3 (1000.318f, 1.29232f, 1000.349f );

		leftTarget = GameObject.Find("InputTarget_L");
		rightTarget = GameObject.Find("InputTarget_R");

		leftHand = GameObject.Find("Hand_L");
		rightHand = GameObject.Find("Hand_R");

		leftHandClose = false;
		rightHandClose = false;

		leftHandAnimate = false;
		rightHandAnimate = false;

		leftMove = false;
		rightMove = false;

		leftPercentage = 20f;
		rightPercentage = 20f;

		leftTarget.transform.position = leftStart;
		rightTarget.transform.position = rightStart;

		comm = new AssemblyCSharp.MotorComm();
		comm.setupExperiment();

		comm.TR = 1.0f;
		tryInterval = comm.TR  / 5f;
		timeCounter = 0;
		
		armTryInterval = comm.TR / 5f;
		armTimeCounter = 0;
		
		handTryInterval = comm.TR / 10f;
		handTimeCounter = 0;
	}

	void toggleHand(GameObject handObj, ref bool HandClose)
	{
		HandClose = !HandClose;
		
		VR_Hand hand = handObj.GetComponent<VR_Hand>();
		if (HandClose)
		{
			hand.thumbInput = 1;
			hand.indexInput = 1;
			hand.middleInput = 1;
			hand.ringInput = 1;
			hand.littleInput = 1;
		}
		else
		{
			hand.thumbInput = 0;
			hand.indexInput = 0;
			hand.middleInput = 0;
			hand.ringInput = 0;
			hand.littleInput = 0;
		}
	}

	void moveHand(GameObject target, ref float percentage, Vector3 start, Vector3 end, float increment)
	{
		float tempValue = percentage;
		if (tempValue > 1) tempValue = 1;
		else if (tempValue < 0) tempValue = 0;

		target.transform.position = start + (end-start) * tempValue;
	}

	// Update is called once per frame
	void Update () {
		if ( Input.GetKey(KeyCode.Space ) ) 
			comm.initRun ();
		
		if (comm.experimentStarted())
			comm.updateVolume();
		else 
			if (comm.checkParallelPort())
				comm.initRun();

		timeCounter -= Time.deltaTime;
		if (timeCounter <= 0) 
		{
			if (comm.experimentStarted())
			   StartCoroutine("stateManager");
			timeCounter = tryInterval;
		};

		armTimeCounter -= Time.deltaTime;
		if (armTimeCounter <= 0) 
		{
			if (leftMove)
				moveHand(leftTarget, ref leftPercentage, leftStart, leftEnd, incr);

			if (rightMove)
				moveHand(rightTarget, ref rightPercentage, rightStart, rightEnd, incr);

			armTimeCounter = armTryInterval;
		};

		handTimeCounter -= Time.deltaTime;
		if (handTimeCounter <= 0) 
		{
			if (leftHandAnimate)
				toggleHand(leftHand, ref leftHandClose);

			if (rightHandAnimate)
				toggleHand(rightHand, ref rightHandClose);

			handTimeCounter = handTryInterval;
		};
	}
}
