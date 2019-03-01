#ifndef SOCKET_H
#define SOCKET_H
// 20091231 bcb     Changed char* to const char* for gcc4.2 warnings
// 20100111 mvh     Merged
// 20100122 mvh     Added Poll function
// 20100619 bcb     Added no-copy to the Socket class.
// 20100707 mvh     Merged

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

****************************************************************************
 *
 * University of California, Davis
 * UCDMC DICOM Network Transport Libraries
 * Version 0.1 Beta
 *
 * Technical Contact: mhoskin@ucdavis.edu
 *
 ***************************************************************************/

/*********************************************************************
 *
 * Connection Class
 *
 *
 *********************************************************************/

// The library used to link with the DICOM Lib has a naming conflict with
// socket.  Therefor we simply remap our socket class into MacSocket.
 
#ifdef	MAC
#	define	Socket	MacSocket
#endif

#ifndef WINDOWS
#include <netinet/in.h>
#include <netdb.h>
#endif

#ifdef WINDOWS

#ifndef _WINDOWS_
#ifndef VOID
#define VOID void
typedef char CHAR;
typedef short SHORT;
typedef long LONG;
#if !defined(MIDL_PASS)
typedef int INT;
#endif
#endif

#ifndef BASETYPES
#define BASETYPES
typedef unsigned long ULONG;
typedef ULONG *PULONG;
typedef unsigned short USHORT;
typedef USHORT *PUSHORT;
typedef unsigned char UCHAR;
typedef UCHAR *PUCHAR;
typedef char *PSZ;
#endif  /* !BASETYPES */

#include <inaddr.h>

typedef struct sockaddr_in {

#if(_WIN32_WINNT < 0x0600)
    short   sin_family;    
#else //(_WIN32_WINNT < 0x0600)
    ADDRESS_FAMILY sin_family;
#endif //(_WIN32_WINNT < 0x0600)

    USHORT sin_port;
    IN_ADDR sin_addr;
    CHAR sin_zero[8];
} SOCKADDR_IN, *PSOCKADDR_IN;

struct  hostent {
        char    FAR * h_name;           /* official name of host */
        char    FAR * FAR * h_aliases;  /* alias list */
        short   h_addrtype;             /* host address type */
        short   h_length;               /* length of address */
        char    FAR * FAR * h_addr_list; /* list of addresses */
#define h_addr  h_addr_list[0]          /* address, for backward compat */
};

struct  servent {
        char    FAR * s_name;           /* official service name */
        char    FAR * FAR * s_aliases;  /* alias list */
#ifdef _WIN64
        char    FAR * s_proto;          /* protocol to use */
        short   s_port;                 /* port # */
#else
        short   s_port;                 /* port # */
        char    FAR * s_proto;          /* protocol to use */
#endif
};

#endif

#endif

#include <sstream>
#include "cctypes.h"

class	Socket
	{
	int					Error;
	unsigned	long	tulong,*tulongptr;
	public:
		int					Socketfd;
		int					ListenSocketfd;
		BOOL				UDP;
		BOOL				Connected;
		BOOL				Listened;
	private:
	struct	sockaddr_in	sa;
	struct	hostent		hes;
	struct	servent		servs;

		int					TimeOut;					
	public:
		virtual	BOOL	SendBinary(BYTE *, UINT);
		virtual	int	ReadBinary(BYTE *, UINT);
		virtual	BOOL	Poll(void);
		struct	hostent	*Gethostbyname(char	*);
		struct	servent	*Getservbyname(const char	*, char *);
		int		GetLinkError();
				Socket();
		virtual	~Socket();		
		virtual	BOOL	Open ( char *ip, char *port );
		BOOL	Close ();
		int		GetLastError()	{ return ( Error ); };
		BOOL	Listen (char *port);
		BOOL	Accept ();
		BOOL	SetTimeOut(int);
		
		// UDP Extensions
		BOOL	BindUDPServer(char *port);
		BOOL	BindUDPClient(char *host, const char *port);
#ifdef __GNUC__
	private:// This will prevent it from being copied (it has a pointer)
		Socket(const Socket&);
		const	Socket & operator = (const Socket&);
#endif
	};

#endif
