using UnityEngine;
using System.Collections;

public class VR_HandTarget_1 : MonoBehaviour {
	public Transform inputTarget;
	public Transform shoulder;

	private bool collided;
	private Collider otherCollider;

	public void OnTriggerEnter(Collider other) {
		if (!collided && other.gameObject.isStatic == true) {
			collided = true;
			otherCollider = other;
		}
	}
	
	public void OnTriggerExit(Collider other) {
		if (collided && other == otherCollider) {
			collided = false;
			otherCollider = null;
		}
	}

	void Update() {
		if (collided) {
			Ray ray = new Ray(shoulder.position, inputTarget.position - shoulder.position);
			RaycastHit hit;
			
			Vector3 dir = inputTarget.position - shoulder.position;
			int rcLayers = Physics.DefaultRaycastLayers & ~(1<<9);
			if (Physics.Raycast(ray, out hit, dir.magnitude + 0.05f, rcLayers)) {
				Vector3 contactPoint = hit.point - dir.normalized * 0.05f;
				transform.position = contactPoint;
			} else
				transform.position = inputTarget.position;
		} else
			transform.position = inputTarget.position;
		transform.rotation = inputTarget.rotation;
	}
	
}
