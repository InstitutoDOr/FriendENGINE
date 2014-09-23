/* VR Body Movements Free
 * author: Pascal Serrarnes
 * email: unity@serrarens.nl
 * version: 2.1.2
 * date: May 17, 2014
 */

using UnityEngine;
using System.Collections;

public class VR_BodyMovementsFree: MonoBehaviour {
	
	public enum BodySide {
		Left,
		Right
	};

	public Transform rightHandTarget = null;
	public Transform leftHandTarget = null;
	public Transform rightFootTarget, leftFootTarget;
	public BodySide side;

	private GameObject rightPalmTarget, leftPalmTarget;
	private Transform rightHandOTarget, leftHandOTarget;


	public bool enableTorso = true;
	public bool enableLegs = true;

	public Transform realMe;
	private Transform characterTransform;
	
	private Arm rightArm, leftArm;
	private Torso torso;
	private Leg rightLeg;
	private Leg leftLeg;

	private Transform leftUpperArm, rightUpperArm;
	
	private static float maxHipAngle = 70;
	
	private float palmLength = 0f;
	
	void Start() {
		characterTransform = this.transform;
		Animator animator = characterTransform.GetComponentInChildren<Animator>();
		
		Transform upperArm, forearm, hand;
		
		upperArm = animator.GetBoneTransform(HumanBodyBones.RightUpperArm);
		forearm = animator.GetBoneTransform(HumanBodyBones.RightLowerArm);
		hand = animator.GetBoneTransform(HumanBodyBones.RightHand);
		rightArm = new Arm(upperArm, forearm, hand, characterTransform, BodySide.Right);

		upperArm = animator.GetBoneTransform(HumanBodyBones.LeftUpperArm);
		forearm = animator.GetBoneTransform(HumanBodyBones.LeftLowerArm);
		hand = animator.GetBoneTransform(HumanBodyBones.LeftHand);
		leftArm = new Arm(upperArm, forearm, hand, characterTransform, BodySide.Left);

		leftUpperArm = animator.GetBoneTransform(HumanBodyBones.LeftUpperArm);
		rightUpperArm = animator.GetBoneTransform(HumanBodyBones.RightUpperArm);

		Transform neck, spine, hips;
		
		neck = animator.GetBoneTransform(HumanBodyBones.Neck);
		spine = animator.GetBoneTransform(HumanBodyBones.Spine);
		hips = animator.GetBoneTransform(HumanBodyBones.Hips);
		torso = new Torso(neck, spine, hips, characterTransform, rightArm);
		
		
		Transform upperLeg, lowerLeg, foot;
		
		upperLeg = animator.GetBoneTransform(HumanBodyBones.RightUpperLeg);
		lowerLeg = animator.GetBoneTransform(HumanBodyBones.RightLowerLeg);
		foot = animator.GetBoneTransform(HumanBodyBones.RightFoot);
		rightLeg = new Leg(upperLeg, lowerLeg, foot, characterTransform);
		
		upperLeg = animator.GetBoneTransform(HumanBodyBones.LeftUpperLeg);
		lowerLeg = animator.GetBoneTransform(HumanBodyBones.LeftLowerLeg);
		foot = animator.GetBoneTransform(HumanBodyBones.LeftFoot);
		leftLeg = new Leg(upperLeg, lowerLeg, foot, characterTransform);
		
		
		if (rightHandTarget != null) {
			float targetScaleZ = rightHandTarget.transform.localScale.z;
			rightPalmTarget = new GameObject();
			rightPalmTarget.name = "Palm_R_Target";
			rightPalmTarget.transform.parent = rightHandTarget.transform;
			rightPalmTarget.transform.localPosition = new Vector3(0, 0, -palmLength/targetScaleZ);
			rightPalmTarget.transform.localEulerAngles = Vector3.zero;

			rightHandOTarget = rightHandTarget;
			rightHandTarget = rightPalmTarget.transform;
		}

		if (leftHandTarget != null) {
			float targetScaleZ = leftHandTarget.transform.localScale.z;
			leftPalmTarget = new GameObject();
			leftPalmTarget.name = "Palm_L_Target";
			leftPalmTarget.transform.parent = leftHandTarget.transform;
			leftPalmTarget.transform.localPosition = new Vector3(0, 0, palmLength/targetScaleZ);
			leftPalmTarget.transform.localEulerAngles = Vector3.zero;
			
			leftHandOTarget = leftHandTarget;
			leftHandTarget  = leftPalmTarget.transform;
		}

		if (enableLegs) {
			if (rightFootTarget == null) {
				GameObject rightFootTargetGO = new GameObject("Foot_R_Target");
				rightFootTarget = rightFootTargetGO.transform;
				rightFootTarget.parent = characterTransform;
				rightFootTarget.position = rightLeg.foot.position;
				rightFootTarget.rotation = characterTransform.rotation;
			}
			
			if (leftFootTarget == null) {
				GameObject leftFootTargetGO = new GameObject("Foot_L_Target");
				leftFootTarget = leftFootTargetGO.transform;
				leftFootTarget.parent = characterTransform;
				leftFootTarget.position = leftLeg.foot.position;
				leftFootTarget.rotation = characterTransform.rotation;
			}
		}
	}
	
	bool crouching = false;
	float bendAngle = 0;
	
	void LateUpdate () {
		
		CalculateTargetPositions();
		
		if (torso.userNeckTarget) 
		{
			if (enableTorso)
				torso.CalculateHorizontal();
			
			if (enableLegs)
				torso.CalculateVertical();

			rightArm.Calculate(rightHandTarget);
			leftArm.Calculate(leftHandTarget);
		} 
		else {
			if (bendAngle <= 0 && !crouching) {
				rightArm.Calculate(rightHandTarget);
				leftArm.Calculate(leftHandTarget);
			}
			
			if (enableTorso && !crouching) 
			{
				bendAngle = torso.AutoHorizontal(rightHandOTarget, rightHandTarget, rightArm, leftUpperArm);
				bendAngle = torso.AutoHorizontal(leftHandOTarget, leftHandTarget, leftArm, rightUpperArm);
			}
			
			if (enableLegs && bendAngle >= maxHipAngle)
			{
				crouching = torso.AutoVertical(rightHandOTarget, rightHandTarget, rightArm);
				crouching = torso.AutoVertical(leftHandOTarget, leftHandTarget, leftArm);
			}
		}
		
		if (enableLegs) {
			rightLeg.Calculate(rightFootTarget.transform);
			leftLeg.Calculate(leftFootTarget.transform);
		}
		
		realMe.rotation = characterTransform.rotation;
	}
	
	void CalculateTargetPositions() {
		realMe.transform.position = characterTransform.position;
	}
	
	public class Torso {
		private Transform characterTransform;
		
		public Transform neck;
		public Transform spine;
		public Transform hips;
		
		public Quaternion fromNorm;
		
		public float length;
		
		public Transform neckStart;
		public Transform hipsStart;
		private Vector3 spineStartOrientation;
		public Quaternion spineStartRotation;
		
		public bool userNeckTarget;
		private Vector3 spineAxis;
		
		public Torso(Transform neck_in, Transform spine_in, Transform hips_in, Transform characterTransform_in, Arm arm) {
			characterTransform = characterTransform_in;
			
			neck = neck_in;
			
			GameObject neckStartGO = new GameObject("Neck_Start");
			neckStart = neckStartGO.transform;
			neckStart.position = neck.position;
			neckStart.parent = characterTransform;
			
			spine = spine_in;
			spineStartRotation = spine.rotation;
			spineStartOrientation = neck.position - spine.position;
			
			hips = hips_in;
			GameObject hipStartGO = new GameObject("Hip_Start");
			hipsStart = hipStartGO.transform;
			hipsStart.position = hips.position;
			hipsStart.parent = characterTransform;
			
			Vector3 neckAtUpperArm = new Vector3(neck.position.x, arm.upperArm.position.y, neck.position.z);
			length = Vector3.Distance(spine.position, neckAtUpperArm);
			
			fromNorm = Quaternion.Inverse(Quaternion.LookRotation(neck.position - spine.position)) * spine.rotation;
			
			spineAxis = spine.InverseTransformDirection(characterTransform.right);
		}
		
		public void CalculateHorizontal() {
			spine.LookAt(neck, characterTransform.forward);
			spine.rotation *= fromNorm;			
		}
		
		public float AutoHorizontal(Transform HandOTarget, Transform HandTarget, Arm oneArm, Transform opositeUpperArm) {
			float bendAngle = 0;
			Vector3 torsoTarget = Vector3.zero;
			Vector3 dShoulderNeck = (opositeUpperArm.position - oneArm.upperArm.position) / 2;
			
			Vector3 ToTarget = HandOTarget.position - oneArm.upperArmStart.position;
			float Over = ToTarget.magnitude - oneArm.length;
			
			if (Over > 0) {
				torsoTarget = HandOTarget.position + dShoulderNeck;
				bendAngle = BendAngle(torsoTarget, oneArm); // arm should be left or right
				spine.rotation = spineStartRotation * Quaternion.AngleAxis(bendAngle, spineAxis);
			} else
				spine.rotation = spineStartRotation;

			oneArm.Calculate(HandTarget);
			
			return bendAngle;
		}
		
		public float BendAngle(Vector3 torsoTarget, Arm arm) {
			float baseAngle = Vector3.Angle (spineStartOrientation, torsoTarget - spine.position);
			
			float dSpine2Target = Vector3.Distance (spine.position, torsoTarget);
			float spineAngle = Mathf.Acos((dSpine2Target * dSpine2Target + length*length - arm.length * arm.length)/ (2 * dSpine2Target * length)) * Mathf.Rad2Deg;
			if (float.IsNaN(spineAngle)) spineAngle = 0;
			
			return Mathf.Min(baseAngle - spineAngle, maxHipAngle);
		}
		
		public void CalculateVertical() {
			float dY = neck.transform.position.y - neck.position.y;
			if (hips.position.y + dY < hipsStart.position.y)  {
				hips.Translate(0, dY, 0, Space.World);
			} else if (hips.position.y + dY > hipsStart.position.y) {
				hips.position = hipsStart.position;
			}
		}
		
		public bool AutoVertical(Transform HandOTarget, Transform HandTarget, Arm oneArm) {
			Vector3 neckDelta = Vector3.zero;
			
			Vector3 ToTarget = HandOTarget.position - oneArm.upperArmStart.position;
			float Over = ToTarget.magnitude - oneArm.length;
			if (Over > 0)
				neckDelta = ToTarget.normalized * Over;

			float dY = neckStart.position.y + neckDelta.y - neck.position.y;

			if (hips.position.y + dY < hipsStart.position.y)  {
				hips.Translate(0, dY, 0, Space.World);
			} else if (hips.position.y + dY > hipsStart.position.y) {
				hips.position = hipsStart.position;
			}

			oneArm.Calculate(HandTarget);

			return (hips.position.y < hipsStart.position.y);
		}
		
	}
	
	public class Arm {

		public Transform upperArm;
		public Transform forearm;
		public Transform hand;
		
		private Quaternion fromNormUpperArm;
		private Quaternion fromNormForearm;
		private Quaternion fromNormHand;
		
		private float upperArmLength, forearmLength;
		public float length;
		
		public Transform upperArmStart;
		
		private float upperArmLength2, forearmLength2;
		
		public Arm(Transform upperArm_in, Transform forearm_in, Transform hand_in, Transform characterTransform, BodySide side) {

			upperArm = upperArm_in;
			forearm = forearm_in;
			hand = hand_in;
			Quaternion rotation;
			
			fromNormUpperArm = Quaternion.Inverse(Quaternion.LookRotation(forearm.position - upperArm.position)) * upperArm.rotation;
			fromNormForearm = Quaternion.Inverse(Quaternion.LookRotation(hand.position - forearm.position)) * forearm.rotation;
			fromNormHand = Quaternion.Inverse(Quaternion.LookRotation(hand.position - forearm.position)) * hand.rotation;

			if (side == BodySide.Left)
			{
			   fromNormUpperArm *= Quaternion.Euler(180, 0, 0);
			   fromNormForearm *= Quaternion.Euler(180, 0, 0);
			   fromNormHand *= Quaternion.Euler(180, 0, 0);
			}

			GameObject upperArmStartGO = new GameObject("Upper_Start");
			upperArmStart = upperArmStartGO.transform;
			upperArmStart.position = upperArm.position;
			upperArmStart.parent = characterTransform;
			
			upperArmLength = Vector3.Distance(upperArm.position, forearm.position);
			forearmLength = Vector3.Distance(forearm.position, hand.position);
			length = Vector3.Distance(upperArm.position, hand.position);
			
			upperArmLength2 = upperArmLength * upperArmLength;
			forearmLength2 = forearmLength * forearmLength;
		}
		
		public void Calculate(Transform handTarget) {
			if (handTarget != null) {
				float dShoulderTarget = Vector3.Distance (upperArm.position, handTarget.position);
				float shoulderAngle = Mathf.Acos((dShoulderTarget * dShoulderTarget + upperArmLength2 - forearmLength2)/ (2 * dShoulderTarget * upperArmLength)) * Mathf.Rad2Deg;
				if (float.IsNaN(shoulderAngle)) shoulderAngle = 0;

				Vector3 handTargetUp;

				handTargetUp = new Vector3(handTarget.up.x, Mathf.Max (handTarget.up.y, 0.1f), handTarget.up.z);
				upperArm.LookAt(handTarget.position, handTargetUp);
				upperArm.rotation = Quaternion.AngleAxis(shoulderAngle, upperArm.up) * upperArm.rotation;
				upperArm.rotation *= fromNormUpperArm;
				
				forearm.LookAt(handTarget.position, handTarget.up);
				forearm.rotation *= fromNormForearm;

				hand.rotation = handTarget.rotation * fromNormHand;
			}
		}
		
		public void StretchToTarget(Transform handTarget) {
			upperArm.LookAt(handTarget.position, handTarget.up);
			upperArm.rotation *= fromNormUpperArm;
			forearm.LookAt(handTarget.position, handTarget.up);
			forearm.rotation *= fromNormForearm;
			
			hand.rotation = handTarget.rotation * fromNormHand;
		}
	}
	
	public class Leg {
		private Transform characterTransform;
		public Transform upperLeg;
		public Transform lowerLeg;
		public Transform foot;
		
		private Quaternion fromNormUpperLeg;
		private Quaternion fromNormLowerLeg;
		private Quaternion fromNormFoot;
		
		private float upperLegLength, lowerLegLength;
		private float upperLegLength2, lowerLegLength2;
		
		public Leg(Transform upperLeg_in, Transform lowerLeg_in, Transform foot_in, Transform characterTransform_in) {
			characterTransform = characterTransform_in;
			
			upperLeg = upperLeg_in;
			lowerLeg = lowerLeg_in;
			foot = foot_in;
			
			fromNormUpperLeg = Quaternion.Inverse(Quaternion.LookRotation (foot.position - upperLeg.position, characterTransform.forward)) * upperLeg.rotation;
			fromNormLowerLeg = Quaternion.Inverse(Quaternion.LookRotation(foot.position - lowerLeg.position, characterTransform.forward)) * lowerLeg.rotation;
			fromNormFoot = Quaternion.Inverse(Quaternion.LookRotation (characterTransform.forward)) * foot.rotation;
			
			upperLegLength = Vector3.Distance(upperLeg.position, lowerLeg.position);
			lowerLegLength = Vector3.Distance(lowerLeg.position, foot.position);
			
			upperLegLength2 = upperLegLength * upperLegLength;
			lowerLegLength2 = lowerLegLength * lowerLegLength;
			
		}
		
		public void Calculate(Transform footTarget) {
			float dHipTarget = Vector3.Distance (upperLeg.position, footTarget.position);
			float hipAngle = Mathf.Acos((dHipTarget * dHipTarget + upperLegLength2 - lowerLegLength2)/ (2 * upperLegLength * dHipTarget)) * Mathf.Rad2Deg;
			if (float.IsNaN(hipAngle)) hipAngle = 0;
			
			upperLeg.LookAt (footTarget.position, characterTransform.forward);
			upperLeg.rotation = Quaternion.AngleAxis(-hipAngle, upperLeg.right) * upperLeg.rotation;
			upperLeg.rotation *= fromNormUpperLeg;
			
			lowerLeg.LookAt (footTarget.position, -characterTransform.up + characterTransform.forward);
			lowerLeg.rotation *= fromNormLowerLeg;
			
			foot.rotation = footTarget.rotation * fromNormFoot;
		}
	}	
}