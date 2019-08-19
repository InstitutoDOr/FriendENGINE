#include <map>
#include <string>
#include "session.h"
#include "cctypes.h"

using namespace std;

// Class to responsible for receiving commands and thread management of the FRIEND processing
class	friendEngine
	{
	public:
                map<string, Session *> sessionList;
      // here the action actually happens
                virtual	bool	serverChild ( int );
      
      // just a function to start a thread
                bool		serverChildThread ( int );
      
      // This is the entry point function. It starts a server that listen in `port` port and calls a function to handle a FRONTEND request (in a threading way or not, depending on variable noThread)
                bool		server ( BYTE * );
		int		childSocketFd;
                char     workingDir[500];
		volatile int	lock;
#ifndef __GNUC__
		friendEngine () { childSocketFd = 0; lock = 0; };
#else //The GCC way
		friendEngine ():childSocketFd(0), lock(0) {};
#endif
		virtual ~friendEngine() {};
	};


