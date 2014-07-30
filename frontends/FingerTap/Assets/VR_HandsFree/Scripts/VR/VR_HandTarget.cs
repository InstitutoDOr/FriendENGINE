using UnityEngine;
using System.Collections;

public class VR_HandTarget : MonoBehaviour {

	public VR_Hand handScript;

	void OnTriggerEnter(Collider other) {
		handScript.OnTargetTriggerEnter(other);
	}

	void OnTriggerExit(Collider other) {
		handScript.OnTargetTriggerExit(other);
	}

}
