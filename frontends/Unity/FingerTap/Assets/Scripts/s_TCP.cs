using UnityEngine;
using System.Collections;
using System;
using System.IO;
using System.Net.Sockets;

public class s_TCP : MonoBehaviour {
	internal Boolean socketReady = false;
	
	TcpClient mySocket;
	NetworkStream theStream;
	StreamWriter theWriter;
	StreamReader theReader;
	
	void Start () {
	}
	void Update () {
	}
	// **********************************************
	public void setupSocket(String Host, Int32 Port) {
		try {
			mySocket = new TcpClient(Host, Port);
			theStream = mySocket.GetStream();
			theWriter = new StreamWriter(theStream);
			theReader = new StreamReader(theStream);

			socketReady = true;
		}
		catch (Exception e) {
			Debug.Log("Socket error: " + e);
		}
	}
	public bool Connected() 
	{
		return mySocket.Connected;
	}
	public void writeSocket(string theLine) {
		if (!socketReady)
			return;
		String foo = theLine + "\r\n";
		theWriter.Write(foo);
		theWriter.Flush();
	}
	public String readSocket() {
		if (!socketReady)
			return "";
		try {
			return theReader.ReadLine ();
		} catch  {
			return "";
		}
	}

	public void closeSocket() {
		if (!socketReady)
			return;
		theWriter.Close();
		theReader.Close();
		mySocket.Close();
		socketReady = false;
	}
} // end class s_TCP