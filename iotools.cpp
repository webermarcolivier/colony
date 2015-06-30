/***************************************************************************//**
 * Project: Colony
 *
 * \file    iotools.cpp
 * \author  Marc Weber\n
 *          The SiMBioSys group (CosmoLab)\n
 *          Parc Científic de Barcelona\n
 *          Barcelona, Spain.\n
 *          http://www.thesimbiosys.com
 * \version 1.0
 * \date    11/2009
 *
 *          Copyright 2009 by Marc Weber
 ******************************************************************************/

#include "iotools.h"

//------------------------------------------------------------------------------

template <class T>
void readOneVariable(ifstream& inputfile, T& x)
{
  string line;
  stringstream stream;
  getline(inputfile,line);
  stream << line;
  stream >> x;
}

/**
 * To avoid linker errors, we have to define the template functions we will use.
 * For more information, see the C++ FAQ:\n
 * http://www.parashift.com/c++-faq-lite/templates.html#faq-35.13
 */
template void readOneVariable<int>   (ifstream& inputfile, int& x);
template void readOneVariable<double>(ifstream& inputfile, double& x);
template void readOneVariable<string>(ifstream& inputfile, string& x);
template void readOneVariable<bool>  (ifstream& inputfile, bool& x);

//------------------------------------------------------------------------------

void jumpNLines(ifstream& inputStream, int n)
{
  int i;
  for (i=0; i<n; i++)
  {
    inputStream.ignore(numeric_limits<streamsize>::max(),'\n');
  }
}

//------------------------------------------------------------------------------

void readWordsInOneLine(ifstream& inputStream, vector<string>& stringVector)
{
  vector<string>::iterator it;

  // Get the line from the inputStream and write it in a stringstream.
  string line;
  getline(inputStream,line);
  stringstream lineStream(line);

  // Erase the content of the string vector.
  stringVector.clear();

  // Write the string elements in the vector.
  string str;
  while (lineStream.good())
  {
    str.clear();
    lineStream >> str;
    stringVector.push_back(str);
  }
}

//------------------------------------------------------------------------------

string insertNumberSuffixToFileName(const string& fileName, int i)
{
  stringstream sStream;
  sStream << setfill('0') << setw(3) << i;
  sStream.seekg (0, ios::beg);
  string suffix;
  sStream >> suffix;

  string fileNameSuffix (fileName);
  fileNameSuffix.insert(fileNameSuffix.size()-4,suffix);

  return fileNameSuffix;
}

//------------------------------------------------------------------------------

ostream& printTimeDifference(ostream& out, double timeDifferenceSeconds)
{
  out << setw(2) << setfill('0') <<  int(timeDifferenceSeconds)/3600 << "h"
      << setw(2) << setfill('0') << (int(timeDifferenceSeconds)%3600)/60 << "m"
      << setw(2) << setfill('0') <<  int(timeDifferenceSeconds)%60 << "s";
  out << setfill(' ');
  return out;
}

//------------------------------------------------------------------------------



