#include <windows.h> 

struct flist
{
  int             num_entries;
  int             max_entries;
  WIN32_FIND_DATA *files;
};

void search(flist &list, wchar_t *root);
void sort(flist &list);
void print(flist *list, wchar_t *root);
