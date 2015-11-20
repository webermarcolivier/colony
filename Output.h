/***************************************************************************//**
 * Project: Colony
 *
 * \file    Output.h
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

#ifndef OUTPUT_H
#define OUTPUT_H

#include "compilation_options.h"

// standard C++ header files
#include <iostream>   // input/output interface
#include <iomanip>    // input/output formatting
#include <fstream>    // manipulate files using streams
#include <string>     // manipulate strings of characters
#include <sstream>    // manipulate strings as stream
#include <utility>
#include <list>

// libraries header files
#include <blitz/array.h>
#include <boost/filesystem.hpp>   // Create, read, delete files, etc.

// user header files
#include "Input.h"
#include "TimeSlice.h"
#include "TimeSliceCellsAvg.h"
#include "CellLineageGeneration.h"
#include "Param/OutputParam.h"

// namespaces
using blitz::Array;
using blitz::Range;
using std::cout;
using std::ostream;
using std::endl;
using std::setw;
using std::ifstream;
using std::ofstream;
using std::string;
using std::stringstream;
using std::ios_base;
using std::pair;
using std::list;
using std::ios;
using std::scientific;
using std::setprecision;
using std::fixed;
using std::setfill;
namespace bf = boost::filesystem;  // Create an alias for the namespace.



/**
 * %Output class with methods to write output files.
 */
class Output
{
public:

//---- LIFECYCLE

  /**
   * Default constructor.
   */
  Output();

  ~Output();

  void initialize(const Input& p);


//---- OPERATORS


//---- ACCESS


//---- OPERATIONS

  void writeTimeMeshTrajectory
  (
    list<TimeSlice>& timeMeshTrajectory,
    int iTraj,
    bool isEnabledComputeSpatialDynamics
  );

  void writeTimeMeshCellsAvg
  (
    Array<TimeSliceCellsAvg,1>& timeMeshCellsAvg,
    int iTraj
  );

  void writeCellLineage
  (
    list<CellLineageGeneration> cellLineage,
    int iTraj
  );

  void writeTimeMeshCellsTrajectoriesAvg
  (
    Array<TimeSliceCellsAvg,1>& timeMeshCellsTrajectoriesAvg
  );

  void writeCellCyclePhases(Array<double,1> cellCyclePhases);

  void writeFirstPassageTime(Array<double,1>& firstPassageTimeValues);

  void close();


private:

//---- DATA

  string path_;

  int nExistingTrajFiles_;

  string simRunningFileNameString_;

  string xFileName_;

  string xConcFileName_;

  string volumeFileName_;

  string cellsAvgFileName_;

  string cellLineageFileName_;

  string positionFileName_;

  string angleFileName_;

  string cellsTrajectoriesAvgFileName_;

  string cellCyclePhasesFileName_;

  string firstPassageTimeFileName_;


};

#endif // OUTPUT_H
