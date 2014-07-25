// various definitions
#ifndef	UNIX
#	define	PATHSEPCHAR	'\\'
#else
#	define	PATHSEPCHAR	'/'
#endif

#define	MAXQUERYLENGTH 310000

#ifdef UNIX

#define stricmp(s1, s2) strcasecmp(s1, s2)
#define memicmp(s1, s2, n) strncasecmp(s1, s2, n)

#define O_BINARY 0              /* LINE UNIX ONLY? */
#define _SH_DENYRW      0x10    /* deny read/write mode */
#define _SH_DENYWR      0x20    /* deny write mode */
#define _SH_DENYRD      0x30    /* deny read mode */
#define _SH_DENYNO      0x40    /* deny none mode */
#define SH_DENYNO       _SH_DENYNO
#define sopen(a, b, c, d) open(a, b) /* LINE UNIX ONLY? */
#define WINAPI /* LINE UNIX ONLY */
#define DeleteCriticalSection(a) pthread_mutex_destroy(a);  /* LINE UNIX ONLY */
#define	closesocket(xxx)	close(xxx) /* LINE UNIX ONLY? */
#define strnicmp(s1, s2, n) strncasecmp(s1, s2, n)
#endif

#define BUFF_SIZE 500

// uppercase a string
void strToUpper(char* str);

// lowercase a string
void strToLower(char* str);

// just removing spaces from an array of char
void removeSpace(char *str);

// just removing `\n` from an array of char
void stripReturns(char *str);
