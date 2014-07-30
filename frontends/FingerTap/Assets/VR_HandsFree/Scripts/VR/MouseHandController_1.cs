using UnityEngine;
using System.Collections;

public class MouseHandController_1 : MonoBehaviour {
	
	public float sensitivityX, sensitivityY, sensitivityZ;
	

	private Vector3 startPosition;
	
	// Use this for initialization
	void Start () {
	}
	
	// Update is called once per frame
	void Update () {
		float x = transform.localPosition.x + Input.GetAxis("Mouse X") * sensitivityX;
		float y = transform.localPosition.y + Input.GetAxis("Mouse ScrollWheel") * sensitivityY;
		float z = transform.localPosition.z + Input.GetAxis("Mouse Y") * sensitivityZ;

		transform.localPosition = new Vector3(x, y, z);
	}
}
