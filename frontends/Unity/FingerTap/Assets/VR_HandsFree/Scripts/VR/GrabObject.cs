using UnityEngine;
using System.Collections;

public class GrabObject : MonoBehaviour {

	public GameObject hand;

	private GameObject grabbedObject = null;
	private Transform objectOrgParent = null;

	void OnCollisionStay(Collision collision) {
		if (Input.GetMouseButton(0) == true) { // Hand is closed
			// the object we grab
			grabbedObject = collision.gameObject;
			// Make it kinematic as we are holding it now
			grabbedObject.rigidbody.isKinematic = true;
			// Store the original parent to restore it when letting loose
			objectOrgParent = grabbedObject.transform.parent;
			// And parent it to the hand
			grabbedObject.transform.parent = hand.transform;
		}
	}

	void FixedUpdate() {
		if (grabbedObject != null) { // are we holding an object?
			if (Input.GetMouseButton(0) == false) { // Hand is open
				// unparent it from the thumb
				grabbedObject.transform.parent = objectOrgParent;
				// make it non-kinematic again
				grabbedObject.rigidbody.isKinematic = false;
				// and clear the grabbed object
				grabbedObject = null;
			}
		}
	}
}
