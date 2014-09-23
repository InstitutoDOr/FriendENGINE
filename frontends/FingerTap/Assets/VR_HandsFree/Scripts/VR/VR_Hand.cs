/* VR Gabbing Hand
 * author: Pascal Serrarnes
 * email: unity@serrarens.nl
 * version: 0.3.1
 * date:  June 7, 2014
 * 
 * changes:
 * - minor performance improvements
 */

using UnityEngine;
using System.Collections;
using System.Collections.Generic;

public class VR_Hand : MonoBehaviour {
	public Animator animator;

	public enum WhichHands { Left, Right };
	public WhichHands whichHand;

	private GameObject hand;
	private Transform shoulder;
	public Transform handTarget;
	private Transform cHandTarget;
	private Vector3 dPalm;

	public float thumbInput;
	public float indexInput;
	public float middleInput;
	public float ringInput;
	public float littleInput;

	private Thumb thumb;
	private Finger indexFinger;
	private Finger middleFinger;
	private Finger ringFinger;
	private Finger littleFinger;

	private Digit[] digits;

	private GameObject grabbedObject = null;
	private Rigidbody grabbedRigidbody = null;
	private GrabbedRigidbody grabbedRBfunc;
	
	private Collider collidedObject;

	private Vector3 palmOffset;
	private GameObject go;

	private enum LeftHandBones {
		ThumbProximal = 24,
		ThumbIntermediate = 25,
		ThumbDistal = 26,
		IndexProximal = 27,
		IndexIntermediate = 28,
		IndexDistal = 29,
		MiddleProximal = 30,
		MiddleIntermediate = 31,
		MiddleDistal = 32,
		RingProximal = 33,
		RingIntermediate = 34,
		RingDistal = 35, 
		LittleProximal = 36,
		LittleIntermediate = 37,
		LittleDistal = 38
	};
	private enum RightHandBones {
		ThumbProximal = 39,
		ThumbIntermediate = 40,
		ThumbDistal = 41,
		IndexProximal = 42,
		IndexIntermediate = 43,
		IndexDistal = 44,
		MiddleProximal = 45,
		MiddleIntermediate = 46,
		MiddleDistal = 47,
		RingProximal = 48,
		RingIntermediate = 49,
		RingDistal = 50, 
		LittleProximal = 51,
		LittleIntermediate = 52,
		LittleDistal = 53
	};

	void Start () {
		hand = this.gameObject;
		thumb = new Thumb();
		indexFinger = new Finger();
		middleFinger = new Finger();
		ringFinger = new Finger();
		littleFinger = new Finger();
		
		if (whichHand == WhichHands.Left) {
			thumb.transform = animator.GetBoneTransform((HumanBodyBones) LeftHandBones.ThumbIntermediate);
			indexFinger.transform = animator.GetBoneTransform((HumanBodyBones) LeftHandBones.IndexProximal);
			middleFinger.transform = animator.GetBoneTransform((HumanBodyBones) LeftHandBones.MiddleProximal);
			ringFinger.transform = animator.GetBoneTransform((HumanBodyBones) LeftHandBones.RingProximal);
			littleFinger.transform = animator.GetBoneTransform((HumanBodyBones) LeftHandBones.LittleProximal);

			shoulder = animator.GetBoneTransform(HumanBodyBones.LeftUpperArm);

		} else {
			thumb.transform = animator.GetBoneTransform((HumanBodyBones) RightHandBones.ThumbIntermediate);
			indexFinger.transform = animator.GetBoneTransform((HumanBodyBones) RightHandBones.IndexProximal);
			middleFinger.transform = animator.GetBoneTransform((HumanBodyBones) RightHandBones.MiddleProximal);
			ringFinger.transform = animator.GetBoneTransform((HumanBodyBones) RightHandBones.RingProximal);
			littleFinger.transform = animator.GetBoneTransform((HumanBodyBones) RightHandBones.LittleProximal);

			shoulder = animator.GetBoneTransform(HumanBodyBones.RightUpperArm);
		}

		GameObject x = new GameObject("Hand Palm");
		x.transform.position = hand.transform.position;
		x.transform.LookAt(indexFinger.transform.transform, Vector3.up);
		Vector3 handRightAxis = hand.transform.InverseTransformDirection(x.transform.right);

		palmOffset = (indexFinger.transform.transform.position - hand.transform.position) * 0.85f;
		palmOffset = hand.transform.InverseTransformDirection(palmOffset);
		palmOffset += new Vector3(0, 0, -0.03f); // to get it into the hand palm HACK

		x.transform.parent = hand.transform;

		if (handTarget.childCount > 0) {
			cHandTarget = handTarget.GetChild(0); // retrieves the palm target
		} else {
			GameObject cHandTargetGO = new GameObject();
			cHandTargetGO.name = "Corrected Target";
			cHandTarget = cHandTargetGO.transform;
			cHandTarget.transform.parent = handTarget;
			cHandTarget.transform.localPosition = Vector3.zero;
			cHandTarget.transform.localEulerAngles = Vector3.zero;
			cHandTarget.transform.localScale = new Vector3(1, 1, 1);
		}
		
		digits = new Digit[5];
		digits[0] = thumb;
		digits[1] = indexFinger;
		digits[2] = middleFinger;
		digits[3] = ringFinger;
		digits[4] = littleFinger;
		
		thumb.Init (hand, handRightAxis, whichHand);
		for (int i = 1; i < digits.Length; i++)
			digits[i].Init(hand, handRightAxis, whichHand);

		MakeTrigger(handTarget.gameObject);
	}

	void MakeTrigger(GameObject targetGO) {
		targetGO.layer = 10;

		SphereCollider sc = targetGO.AddComponent<SphereCollider>();
		sc.isTrigger = true;
		sc.radius = 0.1f / targetGO.transform.localScale.z; //should take localscale into account
		sc.center = new Vector3(0, 0, 0.1f / targetGO.transform.localScale.z);

		Rigidbody rb = targetGO.AddComponent<Rigidbody>();
		rb.isKinematic = true;

		VR_HandTarget ht = targetGO.AddComponent<VR_HandTarget>();
		ht.handScript = this;
	}

	void Update() {
		thumb.Update (thumbInput);
		digits[1].Update(indexInput);
		digits[2].Update(middleInput);
		digits[3].Update(ringInput);
		digits[4].Update(littleInput);

		if (collidedObject != null) 
			ProjectHand(collidedObject);

		Vector3 palmOffsetWorld = hand.transform.TransformDirection(palmOffset);
		if (grabbedRigidbody != null && grabbedRBfunc) {
			grabbedRBfunc.GrabbedUpdate(hand.transform.position + palmOffsetWorld, go);
		}
	}

	bool ProjectHand(Collider collidedObject) {
		Vector3 direction = handTarget.position - shoulder.position;
		RaycastHit[] hits;
		RaycastHit hit;
		bool hitFound = false;
		float minDistance = 100;

		Vector3 orgPosition = transform.position;
		transform.position = shoulder.position;

		cHandTarget.position = handTarget.position;

		hits = rigidbody.SweepTestAll(direction, direction.magnitude);
		for (int i = 0; i < hits.Length; i++) {
			if (hits[i].rigidbody == null) { //only collide with colliders without rigidbodies (often static)
				hitFound = true;
				if (hits[i].distance < minDistance)
					minDistance = hits[i].distance;
			}
		}

		if (hitFound) {
			Vector3 contactPoint = transform.position + (minDistance * direction.normalized);
			cHandTarget.position = contactPoint;
			transform.position = orgPosition;
			return true;
		} else {
			cHandTarget.position = handTarget.position;
			transform.position = orgPosition;
			return false;
		}
	}

	public void OnTargetTriggerEnter(Collider other) {
		if (other.gameObject.isStatic == true)
			collidedObject = other;
	}


	public void OnTargetTriggerExit(Collider other) {
		if (other = collidedObject) {
			if (ProjectHand(other) == false) {
				cHandTarget.position = handTarget.position;
				collidedObject = null;
			}
		}
	}

	void FixedUpdate() {
		if (grabbedObject != null) {
			bool fingersGrabbing = false;
			for (int i = 1; i < digits.Length; i++) {
				Finger finger = (Finger) digits[i];
				if (finger.hasGrabbed())
					fingersGrabbing = true;
			}
			
			if (!fingersGrabbing) {
				GrabbedRigidbody grb = grabbedObject.GetComponent<GrabbedRigidbody>();
				grabbedRBfunc = grb;
				grb.Unkinematize();
				Destroy(grb);
				Destroy(go);

				grabbedObject = null;
				grabbedRigidbody = null;
			}
		}

	}

	void OnCollisionStay (Collision otherCollider) {
		Transform thisTransform, otherTransform;
		
		if (grabbedObject == null) {
			bool fingersCollided = false;
			bool thumbCollided = false;
			
			int ncontacts = otherCollider.contacts.Length;
			for (int i = 0; i < ncontacts; i++ ) {
				thisTransform = otherCollider.contacts[i].thisCollider.transform;
				otherTransform = otherCollider.contacts[i].otherCollider.transform;
				if (thisTransform == thumb.transform || otherTransform == thumb.transform)
					thumbCollided = true;
				for (int j = 1; j < digits.Length; j++) {
					Finger finger = (Finger) digits[j];
					if (thisTransform == finger.transform || otherTransform == finger.transform) {
						fingersCollided = true;
					}
				}
			}
			
			bool grabbed = false;
			// We are touching it both sides
			if (fingersCollided && thumbCollided) {
				for (int j = 1; j < digits.Length; j++) {
					Finger finger = (Finger) digits[j];
					if (finger.input > 0 && !finger.hasGrabbed()) {
						finger.Grab();
				
						if (!grabbed) {
							grabbedObject = otherCollider.gameObject;
							grabbedRigidbody = otherCollider.collider.attachedRigidbody;
							Grab (grabbedObject, grabbedRigidbody);
							grabbed = true;
						}
					}
				}
			}
		}
	}

	void Grab(GameObject obj, Rigidbody rb) {
		Vector3 palmOffsetWorld = hand.transform.TransformDirection(palmOffset);

		GrabbedRigidbody grb = obj.AddComponent<GrabbedRigidbody>();
		grabbedRBfunc = grb;

		go = new GameObject();
		go.transform.parent = hand.transform;
		go.transform.localPosition = Enscale(palmOffset, go.transform.localScale);
		if (whichHand == WhichHands.Right)
			go.transform.localRotation = Quaternion.Euler(20,270,0); // hack, depends on rigs hand axes, need to use norm axes
		else
			go.transform.localRotation = Quaternion.Euler(-20,270,180); // hack, depends on rigs hand axes, need to use norm axes

		grb.Kinematize(hand.transform, whichHand == WhichHands.Right, hand.transform.position + palmOffsetWorld, go);

		//GrabWithHinge(obj);
	}
	
	private Vector3 Enscale(Vector3 position, Vector3 scale) {
		return new Vector3(
			position.x * scale.x,
			position.y * scale.y,
			position.z * scale.z);
	}

	public abstract class Digit {
		public abstract void Init(GameObject hand, Vector3 handRightAxis, WhichHands whichHand);
		
		public abstract void Update(float inputValue);
	}
	
	[System.Serializable]
	public class Finger : Digit {
		public Transform transform;
		public float input;
		private float val;
		private Vector3 axis;
		private float grabAmount = 0;
		private int nPhalanges;
		private Phalanx[] phalanges;
		
		public override void Init(GameObject hand, Vector3 axis_in, WhichHands whichHand) {
			if (transform != null) {
				axis = axis_in;

				phalanges = new Phalanx[3];
				phalanges[0] = new Phalanx(transform);
				if (phalanges[0].transform.childCount == 1) {
					phalanges[1] = new Phalanx(phalanges[0].transform.GetChild(0).transform);
					if (phalanges[1].transform.childCount == 1) {
						phalanges[2] = new Phalanx(phalanges[1].transform.GetChild(0).transform);
						nPhalanges = 3;
					} else {
						phalanges[2] = null;
						nPhalanges = 2;
					}
				} else {
					phalanges[1] = null;
					nPhalanges = 1;
				}
			}
		}

		private float grabSpeed = 0.1f;
		public override void Update(float inputValue) {
			if (transform != null) {
				input = inputValue;
				if (grabAmount > 0)
					val = grabAmount;
				else {
					float d = input - val;
					if (d > grabSpeed)
						val += grabSpeed;
					else if (d < -grabSpeed)
						val -= grabSpeed;
					else
						val += d;
				}
				Bend();
				
				if (grabAmount > 0) {
					if (input < grabAmount) {
						// drop the object
						grabAmount = 0f;
					}
				}
			}
		}
		
		private void Bend() {
			switch (nPhalanges) {
			case 3:
				phalanges[0].Bend(val * 45, axis);
				phalanges[1].Bend(val * 90, axis);
				phalanges[2].Bend(val * 90, axis);
				break;
			case 2:
				phalanges[0].Bend(val * 90, axis);
				phalanges[1].Bend(val * 90, axis);
				break;
			case 1:
				phalanges[0].Bend(val * 135, axis);
				break;
			}
		}
		
		public bool hasGrabbed() {
			return (grabAmount > 0);
		}
		
		public void Grab() {
			grabAmount = val;
		}
	}
	
	[System.Serializable]
	public class Thumb : Digit {
		public Transform transform;
		public float input;
		private float val, grabAmount = 0;
		private Vector3 axis1, axis2;
		private Quaternion startRotation;
		
		public override void Init(GameObject hand, Vector3 handRightAxis, WhichHands whichHand) {
			if (transform != null) {
				axis1 = transform.InverseTransformDirection(hand.transform.up);
				if (whichHand == WhichHands.Left)
					axis2 = transform.InverseTransformDirection(-hand.transform.forward);
				else
					axis2 = transform.InverseTransformDirection(hand.transform.forward);
				startRotation = transform.localRotation;
				startRotation *= Quaternion.AngleAxis(45, -axis1);
			}
		}

		private float grabSpeed = 0.1f;
		public override void Update(float inputValue) {
			if (transform != null) {
				input = inputValue;
				if (grabAmount > 0)
					val = grabAmount;
				else {
					float d = input - val;
					if (d > grabSpeed)
						val += grabSpeed;
					else if (d < -grabSpeed)
						val -= grabSpeed;
					else
						val += d;
				}
				transform.localRotation = startRotation * Quaternion.AngleAxis(val * 45, axis1);
				transform.localRotation *= Quaternion.AngleAxis(val * 30, axis2);
			}
		}
	}
	
	[System.Serializable]
	public class Phalanx {
		public Transform transform;
		public Quaternion startRotation;
		
		public Phalanx(Transform newTransform) {
			transform = newTransform;
			startRotation = transform.localRotation;
		}
		
		public void Bend(float bendFactor, Vector3 axis) {
			transform.localRotation = startRotation * Quaternion.AngleAxis(bendFactor, axis);
		}
	}
}
