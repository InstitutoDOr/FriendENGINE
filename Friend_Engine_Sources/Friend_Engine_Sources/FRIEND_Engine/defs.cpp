#include <cctype>
#include <string.h>
#include <stdio.h>

// uppercase a string
void strToUpper(char* str)
{
    while (*str)
    {
        *str = toupper(*str);
        str++;
    }
}

// lowercase a string
void strToLower(char* str)
{
    while (*str)
    {
        *str = tolower(*str);
        str++;
    }
}

// just removing spaces from an array of char
void removeSpace(char *str)
{
   char *p1 = str, *p2 = str;
   bool noExit;
   do
   {
      while (*p2 == ' ') p2++;
      noExit = (*p1++ = *p2++);
   }
   while (noExit);
}

// just removing `\n` from an array of char
void stripReturns(char *str)
{
   for(unsigned int i=0; i<strlen(str);i++)
   if ((str[i] == 13u) || (str[i] == 10u))
   {
      str[i]=0;
      break;
   }
}

