#define WANT_STREAM
#include "include.h"
#include "newran.h"

#ifdef use_namespace
using namespace NEWRAN;
#endif

void test4(int)
{
   cout << endl << endl;

   {
      cout << "Doing permutations" << endl;
      RandomPermutation RP;
      int i, j;
      int p[10];

      cout << "... select 10 items from 100...119 without replacement" << endl;
      for (i = 1; i <= 10; i++)
      {
         RP.Next(20,10,p,100);
         for (j = 0; j < 10; j++) cout << p[j] << " ";
         cout << "\n";
      }
      cout << "\n";

      cout << "... select 10 items from 100...109 without replacement" << endl;
      for (i = 1; i <= 10; i++)
      {
         RP.Next(10,10,p,100);
         for (j = 0; j < 10; j++) cout << p[j] << " ";
         cout << "\n";
      }
      cout << "\n";
   }

   cout << endl << endl;

   {
      cout << "Doing combinations" << endl;
      RandomCombination RC;
      int i, j;
      int p[10];

      cout << "... select 10 items from 100...119 without replacement" << endl;
      for (i = 1; i <= 10; i++)
      {
         RC.Next(20,10,p,100);
         for (j = 0; j < 10; j++) cout << p[j] << " ";
         cout << "\n";
      }
      cout << "\n";

      cout << "... select 10 items from 100...109 without replacement" << endl;
      for (i = 1; i <= 10; i++)
      {
         RC.Next(10,10,p,100);
         for (j = 0; j < 10; j++) cout << p[j] << " ";
         cout << "\n";
      }
      cout << "\n";
   }
}


