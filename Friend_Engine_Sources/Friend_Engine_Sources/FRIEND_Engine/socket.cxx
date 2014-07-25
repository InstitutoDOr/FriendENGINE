/*
20010426	ljz	Added 'setsockopt' in 'Listen' and 'Open'
			Call 'shutdown' before closing a socket
20010720	ljz     Changed 'print' to 'fprint(stderr,...)'
20050116	mvh     Adapted for LINUX compile
20050119        mvh     Blocked out SOCKET END ERROR message
20070406	mvh	Use reentrant gethostbyname_r when available
20070415	mvh	Small fix there for windows
20090620	jf	Added unistd.h for UNIX
20090824	mvh	Debug output of adress connecting to socket
20091231	bcb	Changed char* to const char* for gcc4.2 warnings
20100111	mvh	Merged
20100122	mvh     Gethostbyname: fixed bug reported by Arash 2station would not be found
20100122	mvh     Added Poll() method: returns TRUE if data available
20100125	mvh     Linux warning
20100309	bcb     Commented out or inited to fix unused variable warning (gcc4.2 Warnings)
20100619	bcb     Fix gcc4 warnings, improve speed and made TimeOut local.
20100707	mvh     Merged
*/
/****************************************************************************
          Copyright (C) 1995, University of California, Davis

          THIS SOFTWARE IS MADE AVAILABLE, AS IS, AND THE UNIVERSITY
          OF CALIFORNIA DOES NOT MAKE ANY WARRANTY ABOUT THE SOFTWARE, ITS
          PERFORMANCE, ITS MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR
          USE, FREEDOM FROM ANY COMPUTER DISEASES OR ITS CONFORMITY TO ANY
          SPECIFICATION. THE ENTIRE RISK AS TO QUALITY AND PERFORMANCE OF
          THE SOFTWARE IS WITH THE USER.

          Copyright of the software and supporting documentation is
          owned by the University of California, and free access
          is hereby granted as a license to use this software, copy this
          software and prepare derivative works based upon this software.
          However, any distribution of this software source code or
          supporting documentation or derivative works (source code and
          supporting documentation) must include this copyright notice.
****************************************************************************/

/***************************************************************************
 *
 * University of California, Davis
 * UCDMC DICOM Network Transport Libraries
 * Version 0.1 Beta
 *
 * Technical Contact: mhoskin@ucdavis.edu
 *
 ***************************************************************************/

/*******************************************************************************
 *
 * socket.C
 *
 *******************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <cstring>
#include <stdio.h>
#include <time.h>

#ifdef UNIX
#       include <unistd.h> //DUCKHEAD92
#endif

#ifndef WINDOWS
#include <arpa/inet.h>
#include <sys/ioctl.h>
#else
#include <ws2tcpip.h>
#endif

#include "socket.hpp"
#include <cerrno>
#include "process.h"
#include <iostream>

int TCPIPTimeOut=60;

#ifdef	WINDOWS
#define SD_BOTH         0x02		/* From winsock2.h */
static	WORD	WinSockUsed	=	0;
static	WSADATA	wsaData;

Socket	::	Socket ()
	{
	WORD		VersionNeeded;
	FARPROC		FarProc;

	TimeOut = 60 * 5;	// 5 minutes	
	UDP = 0;
	if(! WinSockUsed )
		{
		/* WinSock.dll not WSAStartup(), so call it and init the library */
		VersionNeeded = 0x0101;
		
		Error = WSAStartup( VersionNeeded, &wsaData);
		if( Error )
			{
			/* no socket library */
			return;
			}
		/*FarProc = (FARPROC) SocketBlockingHook;
		if(!WSASetBlockingHook(FarProc))
			{
			OutputDebugString("Failed to Install Blocking Hook\n");
			}*/
		}
	++WinSockUsed;
	Socketfd = 0;
	ListenSocketfd = 0;
	Connected = 0;
	Listened = 0;
	}

Socket	::	~Socket ()
	{

	if( ! WinSockUsed )
		return;	// must have failed installation
	
	if ( Socketfd )
		{
		shutdown(Socketfd, SD_BOTH); 
		closesocket(Socketfd);
		}
	if ( ListenSocketfd )
		{
		shutdown(ListenSocketfd, SD_BOTH); 
		closesocket(ListenSocketfd);
		}
	
	--WinSockUsed;
	if(! WinSockUsed )
		{
		//WSAUnhookBlockingHook();
		WSACleanup();
		}
	
	}
#else

#define SYSV_SOCKET

#	define	closesocket(xxx)	close(xxx)

/*
int	closesocket(int	sock)
	{
	int	flags;
	struct	linger	sl;

	sl.l_onoff=1;
	sl.l_linger=30;
	errno=0;
	setsockopt(sock, SOL_SOCKET, SO_LINGER, &sl, sizeof(linger));
	fprintf(stderr, "%d\n", errno);

	errno = 0;
	fprintf(stderr, "CloseSocket(%d) : %d ", sock, close(sock));
	fprintf(stderr, "%d\n", errno);
	return(0);
	}
*/

Socket	::	Socket()
#ifdef __GNUC__ //Faster with member initialization.
:Error(0),
tulong(0),
tulongptr(NULL),
Socketfd(0),
ListenSocketfd(0),
UDP(0),
Connected(0),
Listened(0),
sa(),
hes(),
servs(),
TimeOut(300) {}	// 5 minutes
#else
	{
	Socketfd = 0;
	Connected = 0;
	Listened = 0;
	ListenSocketfd = 0;
	TimeOut = 60 * 5;	// 5 minutes
	UDP = 0;
	}
#endif

Socket	::	~Socket()
	{
	if ( Socketfd )
		{
		closesocket(Socketfd);
		}
	if ( ListenSocketfd )
		{
		closesocket(ListenSocketfd);
		}
	}

#endif
	
#ifndef gethostbyname_r
#define gethostbyname_r(name, a, b, c, d, e) *d=gethostbyname(name)
#endif


struct	hostent	*	Socket	::	Gethostbyname(char	*name1)
	{
	struct	hostent		*he;
	struct  hostent he2;
	unsigned	long	ip;
//	unsigned	long	i1,i2,i3,i4;
//	char		*dot2, *dot3, *dot4;
	char		*name;
	char		ndat[128], buf[64];
	int		dum;
	
	he2.h_length = 0;//Added for compiler unused warning. bcb
	buf[63] = 0;//Added for compiler unused warning. bcb
	dum = 0;//Added for compiler unused warning. bcb

	strcpy(ndat, name1);
	name = &ndat[0];
	
	if(!atoi(name1) || strchr(name1, '.')==NULL)		  // fix for name like 2station
		gethostbyname_r(name, &he2, buf, 64, &he, &dum);
	else
		he = NULL;					  // has . and starts with number
	if(he)
		{
		memcpy ((void *) &hes, (void *) he, sizeof(struct hostent));
		return(&hes);
		}
#ifdef	MAC
	struct	in_addr 	ipa	= inet_addr(name);
	ip = ipa.s_addr;
#else
	ip = inet_addr(name);
#endif
/******************************
	dot2 = strchr(name, '.');
	if(!dot2)
		return ( NULL );
	(*dot2) = '\0';
	++dot2;
	dot3 = strchr(dot2, '.');
	if(!dot3)
		return ( NULL );
	(*dot3) = '\0';
	++dot3;
	dot4 = strchr(dot3, '.');
	if(!dot4)
		return ( NULL );
	(*dot4) = '\0';
	++dot4;
	i1 = atoi(name);
	i2 = atoi(dot2);
	i3 = atoi(dot3);
	i4 = atoi(dot4);
	
#if	NATIVE_ENDIAN == BIG_ENDIAN
	i1 = i4 << 24;
	i2 = i3 << 16;
	i3 = i2 << 8;
#else	
	i4 = i4 << 24;
	i3 = i3 << 16;
	i2 = i2 << 8;
#endif
	
	ip = i1 | i2 | i3 | i4;
********************************/
	hes.h_addr_list = (char **) &tulongptr;
	tulongptr = &tulong;
	tulong = ip;
	return ( &hes );	
	}
	
struct	servent	*	Socket	::	Getservbyname(const char	*name, char	*prot)
	{
	unsigned	short	s;
//	unsigned	short	s1;
	struct	servent		*se;
	
	s = atoi(name);
	if(!s)
		{
		se = getservbyname(name, prot);
		if(!se)
			return(se);
		memcpy ((void *) &servs, (void *) se, sizeof(struct servent));
		return(&servs);
		}
#ifdef	MAC
	servs.s_port = s;
#else
	servs.s_port = htons(s);
#endif
/****************************
	s1 = s & 0x00ff;
	s1 = s1 << 8;
	s = s >> 8;
	s1 = s1 | s;
	servs.s_port = s1;
*******************************/	
//	fprintf(stderr, "returning: %d\n", servs.s_port);
	return ( &servs );
	}

BOOL	Socket	::	Listen ( char 	*port )
	{
	UINT32				tuint32;
//	struct	hostent		*he;
//	struct	hostent		hecopy;
	struct	servent		*se;
	struct	servent		secopy;
	struct	linger		Linger;

	Connected = 0;
	Listened = 0;

	if ( Socketfd )
		closesocket(Socketfd);

//	fprintf(stderr, "[Listen: %s]\n", port); fflush(stderr);
	/*he = Gethostbyname("hampson");
	
	if ( !he )
		{
		return ( FALSE );	// could not resolve host name
		}
	*/   
	//memcpy ((void *) &hecopy, (void *) he, sizeof(struct hostent));
	
	se = Getservbyname(port, NULL);
	
	if ( !se )
		return ( FALSE );	// could not resolv port
	memcpy ((void *) &secopy, (void *) se, sizeof(struct servent));

	Socketfd = socket(PF_INET, SOCK_STREAM, IPPROTO_TCP);
	if(!Socketfd)
		return ( FALSE );	// could not create socket
	Linger.l_onoff  = 1;
	Linger.l_linger = 0;
	setsockopt(Socketfd, SOL_SOCKET, SO_LINGER, (char*)&Linger, sizeof(struct linger));
	//memcpy((void *) &sa.sin_addr.S_un.S_addr , (void *) hecopy.h_addr_list[0], 4);	// adddress
	// should bind to the default address
	tuint32 = 0;
#ifdef	SYSV_SOCKET
	sa.sin_addr.s_addr = tuint32;
#else
	memcpy((void *) &sa.sin_addr.S_un.S_addr, (void *) &tuint32, 4);
#endif
	memcpy((void *) &sa.sin_port , (void *) &secopy.s_port, 2);			// port
	sa.sin_port = secopy.s_port;
	sa.sin_family = AF_INET;

	if(bind(Socketfd, (struct sockaddr *) &sa, sizeof(struct sockaddr_in)))
		{
//		fprintf(stderr, "Erro (bind) : %d\n", GetLinkError());
#ifdef	FORCEDIO
/*		while(errno==125)
			{
			errno=0;
			bind(Socketfd, (struct sockaddr *) &sa, sizeof(struct
				sockaddr_in));
			}*/
#endif
		if(GetLinkError())
			{
			closesocket(Socketfd);
			Error = GetLinkError();
			Socketfd = 0;
			return ( FALSE );
			}
		}
	if(listen(Socketfd, 5))
		{
		//fprintf(stderr, "Erro (bind): %d\n", GetLinkError());fflush(stderr);
		closesocket(Socketfd);
		Socketfd = 0;
		return ( FALSE );
		}

	Listened = 1;
	ListenSocketfd = Socketfd;
	Socketfd = 0;
	Connected = 0;
	return ( TRUE );
	}

BOOL	Socket	::	Accept()
	{
#ifdef	SYSV_SOCKET
	socklen_t	Size;
#else
	int		Size;
#endif
        char str[INET_ADDRSTRLEN];
        time_t Now;

	//fprintf(stderr, "Accept()\n");fflush(stderr);
	if(!Listened)
		return ( FALSE );
	if(!ListenSocketfd)
		return ( FALSE );

	Size = sizeof(struct sockaddr_in);
	Socketfd = accept(ListenSocketfd, (struct sockaddr *) &sa, &Size);
	if(Socketfd>0)
		{
         	    Connected = TRUE;
                Now = time(NULL);
#ifdef WINDOWS
				//disable nagle on the client's socket 
				if (!UDP)
				{
					char value = 1;
					setsockopt(Socketfd, IPPROTO_TCP, TCP_NODELAY, &value, sizeof(value));
				}

                printf("Connected in %s", ctime(&Now));
#else
                inet_ntop(AF_INET, &(sa.sin_addr), str, INET_ADDRSTRLEN);
                printf("Connected by address: %s in %s", str, ctime(&Now));
#endif
		return ( TRUE );
		}
	printf("Error (accept) : %d\n", errno);
	closesocket(ListenSocketfd);
	Listened = 0;
	ListenSocketfd = 0;
	Socketfd = 0;
	Connected = FALSE;
	return ( FALSE );
	}
	

BOOL	Socket	::	Open ( char	*ip, char	*port)
	{
	struct	hostent		*he;
	struct	hostent		hecopy;
	struct	servent		*se;
	struct	servent		secopy;
	struct	linger		Linger;

	if ( Socketfd )
		closesocket(Socketfd);

	Connected = 0;

	he = Gethostbyname(ip);
	
	if ( !he )
		{
	//	fprintf(stderr, "Could not resolve host\n");
		return ( FALSE );	// could not resolve host name
		}
	   
	memcpy ((void *) &hecopy, (void *) he, sizeof(struct hostent));
	
	se = Getservbyname(port, NULL);
	
	if ( !se )
		{
//		fprintf(stderr, "Could not resolve port\n");
		return ( FALSE );	// could not resolv port
		}
	
	memcpy ((void *) &secopy, (void *) se, sizeof(struct servent));
	if(!hecopy.h_addr_list)
		return ( FALSE );
	
	Socketfd = socket(PF_INET, SOCK_STREAM, IPPROTO_TCP);
	if(Socketfd<0)
		return ( FALSE );	// could not create socket
	Linger.l_onoff  = 1;
	Linger.l_linger = 0;
	setsockopt(Socketfd, SOL_SOCKET, SO_LINGER, (char*)&Linger, sizeof(struct linger));
		
#ifdef	SYSV_SOCKET
	memcpy((void *) &sa.sin_addr.s_addr , (void *) hecopy.h_addr_list[0], 4);	// adddress
#else
	memcpy((void *) &sa.sin_addr.S_un.S_addr , (void *) hecopy.h_addr_list[0], 4);	// adddress
#endif
	memcpy((void *) &sa.sin_port , (void *) &secopy.s_port, 2);			// port
	sa.sin_family = AF_INET;
	sa.sin_port = secopy.s_port;

	Error = connect(Socketfd, (struct sockaddr *) &sa, sizeof(struct sockaddr_in));
	if( ! Error )
		{
		Connected = TRUE;
		return ( TRUE );
		}

	closesocket(Socketfd);
	Socketfd = 0;
	return ( FALSE );
	}

BOOL	Socket	::	Close()
	{


	if(ListenSocketfd)
		{
		closesocket(ListenSocketfd);
		}
	if(Socketfd)
		{
		closesocket(Socketfd);
		}
	ListenSocketfd = 0;
	Socketfd = 0;
	Connected = 0;
	Listened = 0;
	UDP = 0;
	return ( TRUE );
	}

BOOL	Socket	::	Poll(void)
	{
#ifdef	SYSV_SOCKET
	fd_set	fds;
#else
	struct	fd_set	fds;
#endif
	struct	timeval	tv	= {0, 0};	// return immediately

	if ( ! Connected )
		return ( FALSE );
#ifdef	MAC
	memset((void*)&fds, 0, sizeof(struct fd_set));
#else
	FD_ZERO(&fds);
#endif
	FD_SET(Socketfd, &fds);
	return (select(Socketfd + 1, &fds, NULL, NULL, &tv)==1);
	}

BOOL		Socket	::	SendBinary(BYTE *s, UINT	count)
	{
	if ( !Connected)
		return ( 0 );
	
	if(send ( Socketfd, (char *) s, count, 0)<0)
		{
//		fprintf(stderr, "SOCKET END ERROR: %d\n", GetLastError());
		}
	return ( TRUE );
	}
	
INT	Socket	::	ReadBinary(BYTE *s, UINT	count)
	{
#ifdef	SYSV_SOCKET
	fd_set	fds;
#else
	struct	fd_set	fds;
#endif
	struct	timeval	tv	= {TimeOut, 0};	// poll

	if ( ! Connected )
		return ( -1 );
#ifdef	MAC
	memset((void*)&fds, 0, sizeof(struct fd_set));
#else
	FD_ZERO(&fds);
#endif
	FD_SET(Socketfd, &fds);
	if(select(Socketfd + 1, &fds, NULL, NULL, &tv)==1)
		{
		if((count = recv(Socketfd, (char *) s, count, 0))==0)
			{
/*#ifdef	WINDOWS
			if(WSAGetLastError()==WSAEDISCON)
				return(-1);
			else
				return(0);
#else
			return ( -1 );
#endif
*/
			Close ();			
			return ( -1 );
			}
		return ( count );
		}
	if ( TimeOut )
		return ( -1 );
	return ( 0 );
	}
	
int		Socket	::	GetLinkError()
	{
#ifdef	WINDOWS
	return ( WSAGetLastError() );
#else
	return ( errno );
#endif
	}


BOOL	Socket	::	SetTimeOut(int	lTimeOut)
	{
	this->TimeOut = lTimeOut;
	return ( TRUE );
	}	

BOOL
Socket	::	BindUDPServer (
	char	*port)
	{
	Close();
	struct	servent		*se;
	struct	servent		secopy;
	se = Getservbyname(port, NULL);
	
	if ( !se )
		{
		return ( FALSE );	// could not resolv port
		}
	
	memcpy ((void *) &secopy, (void *) se, sizeof(struct servent));
	Socketfd = socket(PF_INET, SOCK_DGRAM, IPPROTO_UDP);
	if(Socketfd<0)
		return ( FALSE );	// could not create socket
		
#ifdef	SYSV_SOCKET
	memset((void *) &sa.sin_addr.s_addr , 0, 4);
#else
	memset((void *) &sa.sin_addr.S_un.S_addr , 0, 4);
#endif
	memcpy((void *) &sa.sin_port , (void *) &secopy.s_port, 2);			// port
	sa.sin_family = AF_INET;
	sa.sin_port = secopy.s_port;
	
	if(bind(Socketfd, (struct sockaddr *) &sa, sizeof(struct sockaddr_in)))
		{
		closesocket(Socketfd);
		Error = GetLinkError();
		Socketfd = 0;
		return ( FALSE );
		}
	Connected = 1;
	UDP = 1;
	return ( TRUE );
	}

BOOL	
Socket	::	 BindUDPClient (
	char	*ip,
	const char	*port)
	{
	struct	hostent		*he;
	struct	hostent		hecopy;
	struct	servent		*se;
	struct	servent		secopy;
	

	if ( Socketfd )
		closesocket(Socketfd);
	if ( ListenSocketfd )
		closesocket(ListenSocketfd);
	
	Listened = 0;
	Connected = 0;
	Socketfd = 0;
	ListenSocketfd = 0;
	
	he = Gethostbyname(ip);
	
	if ( !he )
		{
		fprintf(stderr, "Could not resolve host\n");
		return ( FALSE );	// could not resolve host name
		}
	   
	memcpy ((void *) &hecopy, (void *) he, sizeof(struct hostent));
	
	se = Getservbyname(port, NULL);
	
	if ( !se )
		{
		fprintf(stderr, "Could not resolve port\n");
		return ( FALSE );	// could not resolv port
		}
	
	memcpy ((void *) &secopy, (void *) se, sizeof(struct servent));
	if(!hecopy.h_addr_list)
		return ( FALSE );
	
	Socketfd = socket(PF_INET, SOCK_DGRAM, IPPROTO_UDP);
	if(Socketfd<0)
		return ( FALSE );	// could not create socket
		
#ifdef	SYSV_SOCKET
	memcpy((void *) &sa.sin_addr.s_addr , (void *) hecopy.h_addr_list[0], 4);	// adddress
#else
	memcpy((void *) &sa.sin_addr.S_un.S_addr , (void *) hecopy.h_addr_list[0], 4);	// adddress
#endif
	memcpy((void *) &sa.sin_port , (void *) &secopy.s_port, 2);			// port
	sa.sin_family = AF_INET;
	sa.sin_port = secopy.s_port;
	Error = connect(Socketfd, (struct sockaddr *) &sa, sizeof(struct sockaddr_in));
	if( ! Error )
		{
		Connected = TRUE;
		UDP = 1;
		return ( TRUE );
		}
	closesocket(Socketfd);
	Socketfd = 0;
	return ( FALSE );
	}

/*#endif*/

