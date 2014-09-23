using UnityEngine;
using System.Collections;

public class MouseHandController : MonoBehaviour {
	
	public float sensitivityX, sensitivityY, sensitivityZ;
	public VR_Hand vrHand;

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

		float ax = transform.localEulerAngles.x;
		float ay = transform.localEulerAngles.y + Input.GetAxis("Mouse X") * sensitivityY * 100;
		float az = transform.localEulerAngles.z;

		transform.localEulerAngles = new Vector3(ax, ay, az);

		if (Input.GetButton("Fire1")) {
			vrHand.thumbInput = 1;
		} else {
			vrHand.thumbInput = 0;
		}
		if (Input.GetButton("Fire2")) {
			vrHand.indexInput = 1;
			vrHand.middleInput = 1;
			vrHand.ringInput = 1;
			vrHand.littleInput = 1;
		} else {
			vrHand.indexInput = 0;
			vrHand.middleInput = 0;
			vrHand.ringInput = 0;
			vrHand.littleInput = 0;
		}
	}
}
