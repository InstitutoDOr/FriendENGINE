using UnityEngine;
using System.Collections;

public class VRhand_handle : MonoBehaviour {
	public Vector3 position = Vector3.zero;
	public Vector3 rotation = Vector3.zero;
	public float range = 0.1f;

	void OnDrawGizmos() {
		Matrix4x4 m = Matrix4x4.identity;
		Vector3 p = transform.TransformPoint(position);
		Quaternion q = Quaternion.Euler(rotation);
		m.SetTRS(p, transform.rotation * q, Vector3.one); //transform.localScale);
		Gizmos.color = Color.yellow;
		Gizmos.matrix = m;

		Gizmos.DrawCube(Vector3.zero, new Vector3(0.03f, 0.10f,0.04f));
		Gizmos.DrawWireSphere(Vector3.zero, range);
	}

	void Update() {
		Quaternion q = Quaternion.Euler(rotation);
		Debug.DrawLine(transform.position, transform.position + transform.TransformDirection(q * position));
	}
}
