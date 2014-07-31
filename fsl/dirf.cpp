/*
#include "boost/filesystem.hpp"     // includes all needed Boost.Filesystem declarations
#include <iostream>                 // for std::cout
using namespace boost::filesystem;  // for ease of tutorial presentation;
using namespace std;                //  a namespace alias is preferred practice in real code

/*
void show_files( const path & directory, bool recurse_into_subdirs = true )
{
  if( exists( directory ) )
  {
    directory_iterator end ;
    for( directory_iterator iter(directory) ; iter != end ; ++iter )
      if ( is_directory( *iter ) )
      {
        cout << iter->native_directory_string() << " (directory)\n" ;
        if( recurse_into_subdirs ) show_files(*iter) ;
      }
      else 
        cout << iter->native_file_string() << " (file)\n" ;
  }
}


void onefile(path wildcard, string &file)
{
  string mask = wildcard.leaf();
  path directory = wildcard.branch_path();

  if( exists( directory ) )
  {
    directory_iterator end ;
    for( directory_iterator iter(directory) ; iter != end ; ++iter )
      if (is_directory(*iter)) onefile(*iter, file);
      else if (1) file = (*iter).filename();
  }
}
*/

#include <stdio.h>
#include <windows.h>
#include <string>
using namespace std;                

void onefile(string wildcard, string &file)
{
   string path = wildcard;
   WIN32_FIND_DATA f;
   HANDLE h = FindFirstFile(wildcard.c_str(), &f);
   path.erase(path.find_last_of('\\')+1, path.size());
   if(h != INVALID_HANDLE_VALUE)
   {
	  file = path + f.cFileName;
	  FindClose(&f);
   }
   else file = "";
}