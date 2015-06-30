/***************************************************************************//**
 * Project: Colony
 *
 * \file    iotools.h
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

#ifndef IOTOOLS_H
#define IOTOOLS_H

#include "debug.h"

// standard C++ header files
#include <iostream>   // input/output interface
#include <iomanip>    // input/output formatting
#include <fstream>    // manipulate files using streams
#include <string>     // manipulate strings of characters
#include <sstream>    // manipulate strings as stream
#include <vector>
#include <limits>

// libraries header files
#include <blitz/array.h>          // Blitz++ library:
  /* Blitz++ is a C++ class library for scientific computing which provides
  performance on par with Fortran 77/90. It uses template techniques to achieve
  high performance. The current versions provide dense arrays and vectors,
  random number generators, and small vectors and matrices. */

// namespaces
using blitz::Array;
using namespace std;



template <class T>
void readOneVariable(ifstream& inputfile, T& x);

void jumpNLines(ifstream& inputStream, int n);

void readWordsInOneLine(ifstream& inputStream, vector<string>& stringVector);

string insertNumberSuffixToFileName(const string& fileName, int i);

/**
 * Prints a time difference in seconds in the the format 00h00m00s.
 */
ostream& printTimeDifference(ostream& out, double timeDifferenceSeconds);


#endif // IOTOOLS_H
