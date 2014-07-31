#include <stdio.h> 
#include <stdlib.h> 
#include "dirfuncs.h"

void errormessage(void)
{
  LPVOID  lpMsgBuf;
  FormatMessage
  (
    FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
    NULL,
    GetLastError(),
    MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),  // Default language
    (LPTSTR) & lpMsgBuf,
    0,
    NULL
  );
  fprintf(stderr, "%s\n", (char*)lpMsgBuf);
  LocalFree(lpMsgBuf);
}

void format_time(FILETIME *t, char *buff)
{
  FILETIME    localtime;
  SYSTEMTIME  sysloc_time;
  FileTimeToLocalFileTime(t, &localtime);
  FileTimeToSystemTime(&localtime, &sysloc_time);
  sprintf
  (
    buff,
    "%4d/%02d/%02d %02d:%02d:%02d",
    sysloc_time.wYear,
    sysloc_time.wMonth,
    sysloc_time.wDay,
    sysloc_time.wHour,
    sysloc_time.wMinute,
    sysloc_time.wSecond
  );
}

void format_attr(DWORD attr, char *buff)
{
  *buff++ = attr & FILE_ATTRIBUTE_ARCHIVE ? 'A' : '-';
  *buff++ = attr & FILE_ATTRIBUTE_SYSTEM ? 'S' : '-';
  *buff++ = attr & FILE_ATTRIBUTE_HIDDEN ? 'H' : '-';
  *buff++ = attr & FILE_ATTRIBUTE_READONLY ? 'R' : '-';
  *buff++ = attr & FILE_ATTRIBUTE_DIRECTORY ? 'D' : '-';
  *buff++ = attr & FILE_ATTRIBUTE_ENCRYPTED ? 'E' : '-';
  *buff++ = attr & FILE_ATTRIBUTE_COMPRESSED ? 'C' : '-';
  *buff = '\0';
}

void addfile(flist *list, WIN32_FIND_DATA data)
{
  if (list->num_entries == list->max_entries)
  {
    int             newsize = list->max_entries == 0 ? 16 : list->max_entries * 2;
    WIN32_FIND_DATA *temp = (WIN32_FIND_DATA *) realloc(list->files, newsize * sizeof(WIN32_FIND_DATA));
    if (temp == NULL)
    {
      fprintf(stderr, "Out of memory\n");
      exit(1);
    }
    else
    {
      list->max_entries = newsize;
      list->files = temp;
    }
  }

  list->files[list->num_entries++] = data;
}

int sortfiles(const void *a, const void *b)
{
  const WIN32_FIND_DATA *pa = (WIN32_FIND_DATA *) a;
  const WIN32_FIND_DATA *pb = (WIN32_FIND_DATA *) b;
  int                   a_is_dir = pa->dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY;
  int                   b_is_dir = pb->dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY;
  if (a_is_dir ^ b_is_dir)
  {
    // one of each, prefer the directories first
    if (a_is_dir) return(-1);
    if (b_is_dir) return(+1);
    return(0);
  }
  else
  {
    // both files, or both directories - return the strcmp result
    return(wcscmp(pa->cFileName, pb->cFileName));
  }
}

void search(flist &list, wchar_t *root)
{
  HANDLE          h;
  WIN32_FIND_DATA info;

  // build a list of files
  h = FindFirstFile(root, &info);
  if (h != INVALID_HANDLE_VALUE)
  {
    do
    {
      if (!(wcscmp(info.cFileName, L".") == 0 || wcscmp(info.cFileName, L"..") == 0))
      {
        addfile(&list, info);
      }
    } while (FindNextFile(h, &info));
    if (GetLastError() != ERROR_NO_MORE_FILES) errormessage();
    FindClose(h);
  }
  else
  {
    errormessage();
  }
}

void sort(flist &list)
{
  // sort them
  qsort(list.files, list.num_entries, sizeof(list.files[0]), sortfiles);
}

void print(flist *list, wchar_t *root)
{
  int             i;
  // print out in sorted order
  int numdirs = 0;
  for (i = 0; i < list->num_entries; i++)
  {
    char  t1[50], t2[50], t3[50], a[10];
    format_time(&list->files[i].ftCreationTime, t1);
    format_time(&list->files[i].ftLastAccessTime, t2);
    format_time(&list->files[i].ftLastWriteTime, t3);
    format_attr(list->files[i].dwFileAttributes, a);
    if (list->files[i].dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)
    {
      // 'null' date for directory access times, which change each time
      // we run this tool
      sprintf(t2, "%4d/%02d/%02d %02d:%02d:%02d", 2000, 1, 1, 0, 0, 0);
    }

    wprintf(L"%s%s\n", root, list->files[i].cFileName);
    if (list->files[i].dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) numdirs++;
  }
  return;

  // now process all the sub-dirs
  // free all the files first, to save a bit of space
  // the sort function will have put them all at the end.
  list->files = (WIN32_FIND_DATA *) realloc(list->files, numdirs * sizeof(WIN32_FIND_DATA));
  for (i = 0; i < numdirs; i++)
  {
    char  newroot[MAX_PATH];
    sprintf(newroot, "%s\\%s", root, list->files[i].cFileName);
    SetCurrentDirectory(list->files[i].cFileName);
//    doit(newroot);
    SetCurrentDirectory(L"..");
  }

  // free the remainder
  free(list->files);
}

void banner(void)
{
  //       12345678901234567890
  printf("Attribs ");
  printf("      Size ");
  printf("       CreationTime ");
  printf("     LastAccessTime ");
  printf("      LastWriteTime ");
  printf("Filename\n");
}

