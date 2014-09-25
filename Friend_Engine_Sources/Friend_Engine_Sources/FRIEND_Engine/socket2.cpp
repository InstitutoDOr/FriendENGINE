//
//  socket2.cpp
//  
//
//  Created by IDOR on 10/03/14.
//
//

#include "socket2.h"
#ifdef WINDOWS
#include <WinSock2.h>
#else
#include <errno.h>
#endif
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

/// searches for new line and carriage return message delimiters in a string
bool hasNewline(const char *text)
{
   if ((strchr(text, 10) != NULL) || (strchr(text, 13) != NULL)) return 1;
   else return 0;
}

// verifies if the client is alive. Not ready right now
int Socket2::clientAlive(void)
{
   return 1;
}

// verifies if there is more data to read in socket
int Socket2::nextReadSize()
{
#ifdef WIN32
   u_long result = -1;
#else
   long int result = -1;
#endif
   
   int status;
   
#ifdef WIN32
   status = ioctlsocket(Socketfd, FIONREAD, &result);
#else
   status = ioctl(Socketfd, FIONREAD, &result);
#endif
   if (status != 0) connectionProblem = 1;
   return result;
}

// core function to read `size` bytes from communication buffer
int Socket2::read(char *buffer, int size)
{
   if (buffersize<size) readStream();
   ostr.read(buffer, size);
   return 0;
}

int Socket2::readToBuffer(char *buffer, int &dataSize)
{
	int moreToRead = 1, count = 0;
	if (dataSize > sizeof(buffer)) dataSize = sizeof(buffer); // reads for now the minimum between the local buffer and data avaiable
	if (dataSize > 0) // something more ?
	{
		count = recv(Socketfd, buffer, dataSize, 0);
		if (count > 0)
		{
			ostr.write(buffer, count); // writes the data in the stream object
			buffersize += count; // updates the stream buffer size variable
		}
		else
		{
			if (count = 0)
			{
				fprintf(stderr, "connection closed.\n");
				connectionProblem = 1;
			}
			else
			{
/*			
#ifdef WINDOWS
				int error = WSAGetLastError();
#else
				int error = errno();
#endif
*/
				connectionProblem = 1;
			}
			moreToRead = 0;
		}
	}
	return moreToRead;
}

// core function to read all socket data avaiable to a stream object
void Socket2::readStream()
{
   char buffer[1024];
   int dataSize;
   
#ifndef	WINDOWS
   fd_set	fds;
#else
   struct	fd_set	fds;
#endif
   struct	timeval	tv	= {TCPIPTimeOut, 0};	// poll
   
#ifdef	MAC
   memset((void*)&fds, 0, sizeof(struct fd_set));
#else
   FD_ZERO(&fds);
#endif
   FD_SET(Socketfd, &fds);

   int descriptors = select(Socketfd + 1, &fds, NULL, NULL, &tv);
   connectionProblem = 0;
   if(descriptors==1)
   {
      // verifying if exists something to read
	  dataSize = nextReadSize();

	  if (dataSize == 0) // nothing to read. Trying to read something anyway till timeout.
	  {
		  dataSize = 1;
  	      readToBuffer(buffer, dataSize);
	  }
      while ((dataSize=nextReadSize())) // reading all data available
      {
		  if (!readToBuffer(buffer, dataSize)) 
			  break;
      }
   }
   else if (descriptors < 0) // socket potentially invalid
	   connectionProblem = 1;
}

// return the size of last read operation
int Socket2::lastRead()
{
   return ostr.gcount();
}

// assigns a previously created socket to this class
void Socket2::setSocketfd(int Sock)
{
	Socketfd = Sock;
}

// set timeout for the communications
void Socket2::setTimeOut(int TimeOut)
{
	TCPIPTimeOut = TimeOut;
}

//sets the receiving and sending buffer
void Socket2::setBufferSize(int size)
{
   int ret;
   
   ret = setsockopt(Socketfd, SOL_SOCKET, SO_SNDBUF, (char *)&size, sizeof(size));
   
   ret = setsockopt(Socketfd, SOL_SOCKET, SO_RCVBUF, (char *)&size, sizeof(size));
}

// reads a string terminated with a new line from the communication buffer. The new line caracter is excluded.
int Socket2::readLine(char *buffer, int size)
{
   buffer[0] = 0;
   ostr.getline(buffer, size); //reading a line
   buffersize-=ostr.gcount();
   if ((strlen(buffer) >= size-1) || (strlen(buffer)< 1)) // If no newline in buffer or buffer size too short, then trying to read more data from the socket
   {
      // verifying status of the stream object
      if (!ostr.good())
         ostr.clear();
      
      readStream(); //reading socket data
      char *temp = (char *) malloc((size+1) * sizeof(char));
      temp[0] = 0;
      ostr.getline(temp, size); //reading the rest of the line
      buffersize-=ostr.gcount();
      
      strcat(buffer, temp); // concatenating with the previously readed data
      free(temp);
   };
   return 0;
}

// send a string through the socket
int Socket2::writeString(const char *buffer)
{
#ifdef WINDOWS
   char value = 1;
   setsockopt(Socketfd, IPPROTO_TCP, TCP_NODELAY, &value, sizeof(value));
#endif
   int result=send(Socketfd, buffer, strlen(buffer), 0);
   if (result <0) return 0;
   else return 1;
}

// read `size` bytes from communication buffer and verifies if everything is ok
int Socket2::readBuffer(char *buffer, int size)
{
   read(buffer, size);
   if ((lastRead() == size) && (size > 0)) return 1;
   else return 0;
}
