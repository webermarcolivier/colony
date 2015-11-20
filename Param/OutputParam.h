/***************************************************************************//**
 * Project: Colony
 *
 * \file    OutputParam.h
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

#ifndef OUTPUTPARAM_H
#define OUTPUTPARAM_H

#include "../compilation_options.h"

// standard C++ header files
#include<string>

// libraries header files

// user header files

using std::string;



class OutputParam
{
public:

  string outputPath_;

  int nbTimeStepsIntervalConsoleOutput_;

};

#endif // OUTPUTPARAM_H
