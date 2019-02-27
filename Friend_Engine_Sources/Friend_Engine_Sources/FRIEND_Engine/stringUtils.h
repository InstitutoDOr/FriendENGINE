#ifndef _____STRINGUTILS_H
#define _____STRINGUTILS_H
#include <string>
#include <vector>
#include <sstream>

using namespace std;

#define strings vector<string>
#define astrings vector<strings>

string trimString(const std::string &s);
void saveStrings(char *file, strings &list);
void splitLine(string &line, strings &list, char separator = ',');
int stringIndex(strings &list, string token);
void insertToken(strings &vocabulary, string &token);
void deleteOccurrences(string &name, string &toDel);
void mergeStrings(strings &merged, strings &s1, strings &s2);
void clearstrings(strings &list);
void clearastrings(astrings &list);
#endif