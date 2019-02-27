#include "stringUtils.h"

void saveStrings(char *file, strings &list)
{
	FILE *f;
	f = fopen(file, "wt+");
	if (f != NULL)
	{
		for (int t = 0; t < list.size(); t++)
			fprintf(f, "%s\n", list[t].c_str());
		fclose(f);
	}
}

void deleteOccurrences(string &name, string &toDel)
{
	for (basic_string<int>::size_type t = name.find(toDel); t != basic_string<int>::npos; t = name.find(toDel)) name.erase(t, 1);
}

void mergeStrings(strings &merged, strings &s1, strings &s2)
{
	merged.clear();
	for (int t = 0; t < s1.size(); t++) merged.push_back(s1[t]);
	for (int t = 0; t < s2.size(); t++) merged.push_back(s2[t]);
}
void splitLine(string &line, strings &list, char separator)
{
	stringstream ss(line);

	list.clear();
	while (ss.good())
	{
		string token;
		getline(ss, token, separator);
		list.push_back(token);
	}
}

// returns the index of name in vector list
int stringIndex(strings &list, string token)
{
	int r = -1;
	for (int i = 0; i < list.size(); i++)
	{
		if (list[i] == token)
		{
			r = i;
			break;
		}
	}
	return r;
}

void insertToken(strings &vocabulary, string &token)
{
	std::size_t found = token.find("?");
	if (found != std::string::npos)
		token.erase(found, 1);
	if (stringIndex(vocabulary, token) == -1) vocabulary.push_back(token);
}

// remove spaces from string
string trimString(const std::string &s)
{
	string::const_iterator it = s.begin();
	while (it != s.end() && isspace(*it))
		it++;

	string::const_reverse_iterator rit = s.rbegin();
	while (rit.base() != it && isspace(*rit))
		rit++;

	return string(it, rit.base());
}

void clearstrings(strings &list)
{
	list.resize(0);
}

void clearastrings(astrings &list)
{
	for (int t = 0; t < list.size(); t++) list[t].resize(0);
	list.resize(0);
}
