using UnityEngine;
using System.Collections;

public class GrabbedRigidbody : MonoBehaviour {

	private Transform originalParent = null;
	private int originalLayer = 0;
	static int handLayer = 11;

	private bool grabbedHandle;
	private GameObject go2;

	public void Kinematize(Transform handTransform, bool isRightHand, Vector3 targetPos, GameObject go) {
		if (originalParent == null) { // We havent grabbed it yet
			// Make it kinematic as we are holding it now
			rigidbody.isKinematic = true;
			
			// Put it in the MyHands layer (== layer 11) such that it does not collide with my body
			originalLayer = gameObject.layer;
			gameObject.layer = handLayer;
			// However, we need to make everything connected (children, hinges) also to layer 11...
			
			// Store the origGLnal parent in order to restore it when letting loose
			originalParent = transform.parent;

			Vector3 storePos = transform.position;
			transform.position = handTransform.position;
			transform.position = storePos;

			// check for handle
			VRhand_handle handle = gameObject.GetComponent<VRhand_handle>();
			if (handle != null) {

				Vector3 handlePos = transform.TransformPoint(handle.position);
				float d = Vector3.Distance(targetPos, handlePos);
				if (d < handle.range) {
					grabbedHandle = true;
					go2 = new GameObject();
					transform.parent = go2.transform;
					transform.localRotation = Quaternion.Inverse(Quaternion.Euler(handle.rotation));
					transform.localPosition = transform.localRotation * -handle.position ;

				}else {
					grabbedHandle = false;
				}
			}
			if (handle == null | grabbedHandle == false) {
				go.transform.position = transform.position;
				go.transform.rotation = transform.rotation;
			}
		}
	}

	public Vector3 GrabbedUpdate(Vector3 targetPos, GameObject go) {

		VRhand_handle handle = gameObject.GetComponent<VRhand_handle>();
		if (handle != null && grabbedHandle) {
			go2.transform.position = go.transform.position;
			go2.transform.rotation = go.transform.rotation;
		} else {
			transform.position = go.transform.position;
			transform.rotation = go.transform.rotation;
		}
		return Vector3.zero;
	}

	public void Unkinematize() {
		//Destroy(go);
		Destroy(go2);
		// make it non-kinematic again
		rigidbody.isKinematic = false;
		// unparent it from the thumb
		transform.parent = originalParent;
		originalParent = null;
		// restore the layer it was in
		gameObject.layer = originalLayer;
		// and clear the grabbed object
	}
}
