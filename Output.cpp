/***************************************************************************//**
 * Project: Colony
 *
 * \file    Output.cpp
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

#include "Output.h"

//------------------------------------------------------------------------------

Output::Output()
{}

//------------------------------------------------------------------------------

Output::~Output()
{
  // Deleting dummy file indicating running simulation in output directory.
  bf::path file1(simRunningFileNameString_);
  bf::remove(file1);
}

//------------------------------------------------------------------------------

void Output::initialize(const Input& p)
{
  // Write the new directory path into a string.
  path_ = p.out.outputPath_;

  // Output directory with name corresponding to the global parameters value.
  stringstream pathSS;
  pathSS << path_ << "/";
  vector<string>::iterator it;
  int i;
  for (it=Input::listGlobalParametersName.begin(),i=0;
       it!=Input::listGlobalParametersName.end();
       it++,i++)
  {
    if (i!=0)
    {
      pathSS << "_";
    }
    // Write the name of the directory with the initial value of the global parameters.
    //pathSS << *it << scientific << setprecision(3) << Input::globalParameter(i);
    // Write the name of the directory with the final value of the global parameters.
    int lastIndex = Input::globalParameterTimePoints.ubound(blitz::secondDim);
    pathSS << *it << scientific << setprecision(3) << Input::globalParameterTimePoints( 0, lastIndex, i );
  }
  pathSS >> path_;

  // Test if directory exists, create directory and check for existing files.

  nExistingTrajFiles_ = 0;

  bf::path dir1(path_);   // Boost Filesystem path object.

  if ( bf::exists(dir1) )
  {
    if ( bf::is_directory(dir1) )
    {
      // Find the number of trajectories already written in the directory.
      bf::directory_iterator end_itr; // default construction yields past-the-end
      // Iterate on the file list in directory
      for ( bf::directory_iterator itr( dir1 ); itr != end_itr; ++itr )
      {
        if ( !is_directory(itr->status()) )
        {
          // Search for "xConc_trajXXX.dat" file and extract trajectory number
          if ( itr->path().filename().string().substr(0,10).compare("xConc_traj") == 0)
          {
            int iTraj;
            stringstream ssTraj(itr->path().filename().string().substr(10,3));
            ssTraj >> iTraj;
            nExistingTrajFiles_ = max(nExistingTrajFiles_,iTraj+1);
          }
        }
      }

    } else {
      // dir1 is not a directory, delete it.
      bf::remove(dir1);
      // Create directory.
      bf::create_directory(dir1);
    }
  } else {
    // Create directory.
    bf::create_directory(dir1);
  }

  if (nExistingTrajFiles_ > 0)
  {
    cout << "Found existing simulation output files in directory \"" << dir1.string() << "\"\n";
    cout << "Number of existing trajectories: " << nExistingTrajFiles_ << endl;
    cout << "Simulation output trajectories data files will start with number (counter starts at 0): " <<
            nExistingTrajFiles_ << "\n\n\n" << endl;
  }


  // Writing a dummy file to indicate that simulation will write output in this directory.
  simRunningFileNameString_ = path_ + "/waiting_for_simulation_output";
  const char* simRunningFileName (simRunningFileNameString_.c_str());
  ofstream simRunningOutputFile (simRunningFileName);
  simRunningOutputFile.close();

  xFileName_                    = path_ + "/x_traj.dat";
  xConcFileName_                = path_ + "/xConc_traj.dat";
  volumeFileName_               = path_ + "/volumes_traj.dat";
  cellsAvgFileName_             = path_ + "/cells_avg_traj.dat";
  cellLineageFileName_          = path_ + "/cell_lineage_traj.dat";
  positionFileName_             = path_ + "/position_traj.dat";
  angleFileName_                = path_ + "/angle_traj.dat";
  cellsTrajectoriesAvgFileName_ = path_ + "/cells_trajectories_avg.dat";
  cellCyclePhasesFileName_      = path_ + "/cell_cycle_phases.dat";
  firstPassageTimeFileName_     = path_ + "/firstPassageTime.dat";
}

//------------------------------------------------------------------------------

void Output::close()
{
  // Deleting dummy file indicating running simulation in output directory.
  bf::path file1(simRunningFileNameString_);
  bf::remove(file1);
}

//------------------------------------------------------------------------------

void Output::writeTimeMeshTrajectory(list<TimeSlice>& timeMeshTrajectory,
                                     int iTraj,
                                     bool isEnabledComputeSpatialDynamics)
{
  // Note: I should write a general method to write the data that has the
  // same structure.

  string filenameString2 (
    insertNumberSuffixToFileName(xFileName_,iTraj + nExistingTrajFiles_) );
  const char* filename2 (filenameString2.c_str());
  ofstream outputfile2 (filename2);

  string filenameString5 (
    insertNumberSuffixToFileName(xConcFileName_,iTraj + nExistingTrajFiles_) );
  const char* filename5 (filenameString5.c_str());
  ofstream outputfile5 (filename5);

  string filenameString4 (
    insertNumberSuffixToFileName(positionFileName_,iTraj + nExistingTrajFiles_) );
  const char* filename4 (filenameString4.c_str());
  ofstream outputfile4 (filename4);

  string filenameString6 (
    insertNumberSuffixToFileName(angleFileName_,iTraj + nExistingTrajFiles_) );
  const char* filename6 (filenameString6.c_str());
  ofstream outputfile6 (filename6);


  if ( outputfile2.is_open() && outputfile4.is_open() && outputfile5.is_open() && outputfile6.is_open())
  {

    list<TimeSlice>::iterator it;
    for (it = timeMeshTrajectory.begin();
         it!= timeMeshTrajectory.end();
         it++)
    {
      outputfile2 << "   0, " << setw(7) << fixed << (*it).nCells_ << ", "
      << scientific << setprecision(10) <<  (*it).time_;
      outputfile5 << "   0, " << setw(7) << fixed << (*it).nCells_ << ", "
      << scientific << setprecision(10) <<  (*it).time_;

      if (isEnabledComputeSpatialDynamics)
      {
        outputfile4 << "   0, " << setw(7) << fixed << (*it).nCells_ << ", "
        << scientific << setprecision(10) <<  (*it).time_;
        outputfile6 << "   0, " << setw(7) << fixed << (*it).nCells_ << ", "
        << scientific << setprecision(10) <<  (*it).time_;
      }

      int i;

      for ( i=0; i < (*it).t0State_.size(); ++i)
      {
        outputfile2 << fixed <<  ", " << setw(16) << (*it).t0State_(i);
      }
      for ( i=0; i < (*it).xConc_.size(); ++i)
      {
        outputfile5 << fixed <<  ", " << setw(16) << (*it).xConc_(i);
      }

      if (isEnabledComputeSpatialDynamics)
      {
        for ( i=0; i < (*it).positionArray_.size(); ++i)
        {
          outputfile4 << fixed <<  ", " << setw(16) << (*it).positionArray_(i)(0);
          outputfile4 << fixed <<  ", " << setw(16) << (*it).positionArray_(i)(1);
        }
        for ( i=0; i < (*it).angleArray_.size(); ++i)
        {
          outputfile6 << fixed <<  ", " << setw(16) << (*it).angleArray_(i);
        }
      }

      outputfile2 << '\n';
      outputfile5 << '\n';
      if (isEnabledComputeSpatialDynamics)
      {
        outputfile4 << '\n';
        outputfile6 << '\n';
      }
    }

  outputfile2 << endl;
  outputfile5 << endl;
  if (isEnabledComputeSpatialDynamics)
  {
    outputfile4 << endl;
    outputfile6 << '\n';
  }

  } else {
    cout << "ERROR: Unable to open file \"" << filename2 << '\"'  << " and \"" << filename4 << '\"'
         << " and \"" << filename5 << '\"' << " and \"" << filename6 << '\"' << endl;
  }

  if (!isEnabledComputeSpatialDynamics)
  {
    // Delete the files positions and angles, which are empty
    outputfile4.close();
    outputfile6.close();
    bf::path file4(filenameString4);
    bf::remove(file4);
    bf::path file6(filenameString6);
    bf::remove(file6);
  }


  #ifdef TIME_DEPENDENT_PROPENSITIES

  string filenameString3 (
    insertNumberSuffixToFileName(volumeFileName_,iTraj + nExistingTrajFiles_) );
  const char* filename3 (filenameString3.c_str());
  ofstream outputfile3 (filename3);

  if ( outputfile3.is_open() )
  {

    list<TimeSlice>::iterator it;
    for (it = timeMeshTrajectory.begin();
         it!= timeMeshTrajectory.end();
         it++)
    {
      outputfile3 << "   0, " << setw(7) << fixed << (*it).nCells_ << ", "
      << scientific << setprecision(10) <<  (*it).time_;

      int i;

      for ( i=0; i < (*it).volumeArray_.size(); i++)
      {
        outputfile3 << scientific << setprecision(10) << ", " << (*it).volumeArray_(i);
      }

      outputfile3 << '\n';
    }

  outputfile3 << endl;

  } else {
    cout << "ERROR: Unable to open file \"" << filename3 << '\"' << endl;
  }

  #endif
}

//------------------------------------------------------------------------------

void Output::writeTimeMeshCellsAvg
(
  Array<TimeSliceCellsAvg,1>& timeMeshCellsAvg,
  int iTraj
)
{
  string filenameString1 (
    insertNumberSuffixToFileName(cellsAvgFileName_,iTraj + nExistingTrajFiles_) );

  const char* filename (filenameString1.c_str());
  ofstream outputfile (filename);

  if ( outputfile.is_open() )
  {
    int i;
    for (i=0; i<timeMeshCellsAvg.size(); i++)
    {
      outputfile << setw(4) << fixed
      << timeMeshCellsAvg(i).nSamples_
      << ", " << setw(7) << setprecision(2)
      << timeMeshCellsAvg(i).nCellsAvg_ /
           double(timeMeshCellsAvg(i).nSamples_)
      << ", " << scientific << setprecision(10)
      << timeMeshCellsAvg(i).timeStart_;

      int j;
      for ( j=0; j < timeMeshCellsAvg(i).stateAvg_.size(); j++)
      {
        outputfile << scientific << setprecision(10) << ", "
        << timeMeshCellsAvg(i).stateAvg_(j) /
           double(timeMeshCellsAvg(i).nSamples_);
      }
      outputfile << '\n';
    }

  outputfile << endl;

  } else {
    cout << "ERROR: Unable to open file \"" << filename << endl;
  }

}

//------------------------------------------------------------------------------

void Output::writeCellLineage
(
  list<CellLineageGeneration> cellLineage,
  int iTraj
)
{
  string filenameString1 (
    insertNumberSuffixToFileName(cellLineageFileName_,iTraj + nExistingTrajFiles_) );

  const char* filename1 (filenameString1.c_str());
  ofstream outputfile1 (filename1);

  if ( outputfile1.is_open() )
  {

    list<CellLineageGeneration>::iterator it;
    for (it = cellLineage.begin();
         it!= cellLineage.end();
         it++)
    {
      outputfile1
      << scientific << setprecision(10) <<  (*it).time_
      << ", " << setw(7) << fixed << (*it).nCells_
      << ", " << setw(7) << fixed << (*it).nCellsDelta_;

      int i;
      for ( i=0; i < (*it).motherCellsIndices_.size(); i++)
      {
        outputfile1 << ", " << setw(7) << fixed << (*it).motherCellsIndices_(i);
      }

      outputfile1 << '\n';
    }

  outputfile1 << endl;

  } else {
    cout << "ERROR: Unable to open file \"" << filename1 << '\"' << endl;
  }


}

//------------------------------------------------------------------------------

void Output::writeTimeMeshCellsTrajectoriesAvg
(
  Array<TimeSliceCellsAvg,1>& timeMeshCellsTrajectoriesAvg
)
{
  const char* filename (cellsTrajectoriesAvgFileName_.c_str());
  ofstream outputfile (filename);

  if ( outputfile.is_open() )
  {
    int i;
    for (i=0; i<timeMeshCellsTrajectoriesAvg.size(); i++)
    {
      outputfile << setw(4) << fixed
      << timeMeshCellsTrajectoriesAvg(i).nSamples_
      << ", " << setw(7) << setprecision(2)
      << timeMeshCellsTrajectoriesAvg(i).nCellsAvg_ /
           double(timeMeshCellsTrajectoriesAvg(i).nSamples_)
      << ", " << scientific << setprecision(10)
      << timeMeshCellsTrajectoriesAvg(i).timeStart_;

      int j;
      for ( j=0; j < timeMeshCellsTrajectoriesAvg(i).stateAvg_.size(); j++)
      {
        outputfile << scientific << setprecision(10) << ", "
        << timeMeshCellsTrajectoriesAvg(i).stateAvg_(j) /
           double(timeMeshCellsTrajectoriesAvg(i).nSamples_);
      }
      outputfile << '\n';
    }

  outputfile << endl;

  } else {
    cout << "ERROR: Unable to open file \"" << filename << "\"." << endl;
  }
}

//------------------------------------------------------------------------------

void Output::writeCellCyclePhases(Array<double,1> cellCyclePhases)
{
  const char* filename (cellCyclePhasesFileName_.c_str());
  ofstream outputfile (filename);

  if ( outputfile.is_open() )
  {
    int i;
    for (i=0; i<cellCyclePhases.size(); i++)
    {
      outputfile
      << scientific << setprecision(10)
      << cellCyclePhases(i);
      outputfile << '\n';
    }

  outputfile << endl;

  } else {
    cout << "ERROR: Unable to open file \"" << filename << "\"." << endl;
  }
}

//------------------------------------------------------------------------------

void Output::writeFirstPassageTime
(
  Array<double,1>& firstPassageTimeValues
)
{
  const char* filename (firstPassageTimeFileName_.c_str());
  ofstream outputfile (filename);

  if ( outputfile.is_open() )
  {
    int i;
    for (i=0; i<firstPassageTimeValues.size(); i++)
    {
      outputfile << scientific << setprecision(10) << firstPassageTimeValues(i) << '\n';
    }

  } else {
    cout << "ERROR: Unable to open file \"" << filename << endl;
  }

}

//------------------------------------------------------------------------------








