using UnityEngine;
using System.Collections;
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
	private int[] triggerIndex = new int[] {16, 53, 90, 127, 164, 201, 238, 275};
	private int blockLength = 22;

	private String HostData = "127.0.0.1";
	private Int32 Port = 5678;
	private float TR = 1.2f;


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

	IEnumerator stateManager()
	{
		if (comm.state() == 0)
			texto.text = "RELAX"; 

		comm.coreCommunication();
		if (comm.state() == 15)
		{
			texto.text = "RELAX";
			resetState();
			if (comm.volume() >= triggerIndex[comm.phase()]) comm.setState(16);
		}
		
		if (comm.state() == 16) 
		{
			if (comm.phase() > 7) comm.setState(100);
			else
			{
				texto.text = "FINGERTAP ONE HAND ONLY";
				comm.setState(17);
			}
		}
		
		if (comm.state() == 17)
		{
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
			
			if (comm.volume () >= triggerIndex[comm.phase()]+blockLength) 
			{
				comm.setPhase(comm.phase() + 1);
				comm.setState(15);
			}
		}

		if (comm.state () == 100) 
		{
			StopCoroutine("stateManager");
			texto.text = "END";
			
			resetState ();
		}
		comm.handleGraphInformation();
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

		tryInterval = TR  / 5f;
		timeCounter = 0;

		armTryInterval = TR / 5f;
		armTimeCounter = 0;

		handTryInterval = TR / 10f;
		handTimeCounter = 0;

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
		comm.setupConnection(HostData, Port);

		comm.addConfigurationPair ("MNIMask", "studydirhmat_spm_final.nii");
		comm.addConfigurationPair ("MNITemplate", "studydirMNI152_T1_1mm_brain.nii.gz");
		comm.addConfigurationPair ("Prefix", "outputdirRUN01\\DRIN-");
		comm.addConfigurationPair ("ActivationLevel", "0.01");
		comm.addConfigurationPair ("CurrentRunSuffix", "RUN01");

		comm.setRunsize (296);
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
		timeCounter -= Time.deltaTime;
		if (timeCounter <= 0) 
		{
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
