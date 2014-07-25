//
//  socket2.h
//  
//
//  Created by IDOR on 10/03/14.
//
//

#ifndef ____socket2__
#define ____socket2__

#ifndef WINDOWS
#include <netdb.h>
#include <sys/ioctl.h>
#endif

#include <sstream>


// class responsible for handling the TCP/IP communication
// for FRIEND Engine, based on a previously created socket
class Socket2
{
private:
   int buffersize;
   
   // verifies if there is more data to read in socket
   int nextReadSize();
   
   // core function to read `size` bytes from communication buffer
   int read(char *buffer, int size);
   
   // core function to read all socket data avaiable to a stream object
   void readStream();
   
   // return the size of last read operation
   int lastRead();
   
public:
   std::stringstream ostr;
   
   int Socketfd, TCPIPTimeOut;
   
   // assigns a previously created socket to this class
   void setSocketfd(int Sock);
   
   // Set timeout for the communications
   void setTimeOut(int TimeOut);
   
   //sets the receiving and sending buffer
   void setBufferSize(int size);
   
   // Reads a string terminated with a new line from the communication buffer. The new line caracter is excluded
   int readLine(char *buffer, int size);
   
   // send a string through the socket
   int writeString(const char *buffer);
   
   // read `size` bytes from communication buffer and verifies if everything is ok
   int readBuffer(char *buffer, int size);
   
   // function to check the health of the tcp connection
   int	clientAlive(void);
   
	Socket2() {buffersize=0;};
	virtual	~Socket2() {};
};

#endif /* defined(____socket2__) */
