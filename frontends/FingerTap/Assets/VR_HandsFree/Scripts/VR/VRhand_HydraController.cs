/* VR Hands Hydra Controller
 * author: Pascal Serrarnes
 * email: unity@serrarens.nl
 * version: 0.1
 * date: May 10, 2014
 */

using UnityEngine;
using System.Collections;

public class VRhand_HydraController : SixenseObjectController
{
	public VR_Hand vrGrabHand;

	protected override void Start() 
	{
		base.Start();
	}
	
	protected override void UpdateObject( SixenseInput.Controller controller )
	{
		if (controller.Enabled)
			UpdateHand(controller);

		if (!m_enabled && controller.Trigger > 0.01f) {
			m_enabled = true;
			
			// delta controller position is relative to this point
			m_baseControllerPosition = new Vector3( controller.Position.x * Sensitivity.x,
			                                       controller.Position.y * Sensitivity.y,
			                                       controller.Position.z * Sensitivity.z );
			
			// this is the new start position
			m_initialPosition = this.gameObject.transform.localPosition;
		}
		if (m_enabled && controller.GetButton(SixenseButtons.START)) {
			m_enabled = false;
		}

		if ( m_enabled ) {
			UpdatePosition( controller );
			UpdateRotation( controller );
		}
	}
	
	
	void OnGUI()
	{
		if (Hand == SixenseHands.UNKNOWN )
			return;

		if (!m_enabled) {
			int labelWidth = 250;
			int labelPadding = 120;
			int horizOffset = Hand == SixenseHands.LEFT ? -labelWidth - labelPadding  : labelPadding;
			
			string handStr = Hand == SixenseHands.LEFT ? "left" : "right";
			GUI.Box( new Rect( Screen.width / 2 + horizOffset, Screen.height - 40, labelWidth, 30 ),  "Pull " + handStr + " TRIGGER to control hand");		
		}		
	}

	// Updates the animated object from controller input.
	private void UpdateHand(SixenseInput.Controller controller)
	{

		if (Hand == SixenseHands.RIGHT && vrGrabHand != null) {
			vrGrabHand.thumbInput = 0f;
			if (controller.GetButton(SixenseButtons.ONE))
				vrGrabHand.thumbInput = 0.5f;
			if (controller.GetButton(SixenseButtons.THREE))
				vrGrabHand.thumbInput = 1f;

			vrGrabHand.indexInput = controller.GetButton(SixenseButtons.BUMPER) ? 1 : 0;
			vrGrabHand.middleInput = controller.Trigger;
			vrGrabHand.ringInput = controller.Trigger;
			vrGrabHand.littleInput = controller.Trigger;
		}

		if (Hand == SixenseHands.LEFT && vrGrabHand != null) {
			vrGrabHand.thumbInput = 0f;
			if (controller.GetButton(SixenseButtons.TWO))
				vrGrabHand.thumbInput = 0.5f;
			if (controller.GetButton(SixenseButtons.FOUR))
				vrGrabHand.thumbInput = 1f;

			vrGrabHand.indexInput = controller.GetButton(SixenseButtons.BUMPER) ? 1 : 0;
			vrGrabHand.middleInput = controller.Trigger;
			vrGrabHand.ringInput = controller.Trigger;
			vrGrabHand.littleInput = controller.Trigger;
		}
	}

	new private void UpdatePosition( SixenseInput.Controller controller )
	{
		Vector3 controllerPosition = new Vector3( controller.Position.x * Sensitivity.x,
		                                         controller.Position.y * Sensitivity.y,
		                                         controller.Position.z * Sensitivity.z );
		
		// distance controller has moved since enabling positional control
		Vector3 vDeltaControllerPos = controllerPosition - m_baseControllerPosition;
		
		// update the localposition of the object
		this.gameObject.transform.localPosition = m_initialPosition + vDeltaControllerPos;
	}
	
	
	new private void UpdateRotation( SixenseInput.Controller controller )
	{
		this.gameObject.transform.localRotation = controller.Rotation * m_initialRotation;
	}
}

