/***************************************************************************//**
 * Project: Colony
 *
 * \file    Simulator.h
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

#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "compilation_options.h"

// standard C++ header files
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <time.h>
#include <limits.h>
#include <vector>
#include <cmath>

// libraries header files
#include <blitz/array.h>
#include <gsl/gsl_math.h>
#ifdef GUI
  #include <QApplication>
  #include <QMainWindow>
#endif // GUI

// user header files
#include "Output.h"
#include "Input.h"
#include "Param/SimulatorParam.h"
#include "Param/GraphicsCellCompositeParam.h"
#include "TimeSlice.h"
#include "TimeSliceCellsAvg.h"
#include "CellLineageGeneration.h"
#include "CellCollection.h"
#include "IntegratorContext.h"
#include "IntegratorGillespie.h"
#include "IntegratorGillespieModified.h"
#include "IntegratorChemicalLangevin.h"
#include "SpatialIntegratorContext.h"
#include "SpatialIntegratorODE.h"

// namespaces
using blitz::Array;
using blitz::Range;
using std::cout;
using std::ostream;
using std::endl;
using std::setw;
using std::min;




/**
  * Core of the simulation: interaction between the
  * integrator and the cell collection.
  *
  * This class controls the simulation loop
  * over time. It has a cell collection and an interface to the integrator used
  * for computing the dynamics. It outputs the trajectories and other data to
  * the Output object.
  */
class Simulator
{
public:

//---- LIFECYCLE

  void init();

  /**
   * Default constructor initialize the cell collection with input object.
   * Initializes the cell collection, the integrator context,
   * the output object and the time variable.
   */
  Simulator(string inputFile);
  Simulator(string inputFile, GraphicsCellCompositeParam& graphicsCellParam);

  ~Simulator();

  void initMainWindowParam(GraphicsCellCompositeParam& graphicsCellParam);

  void initGraphicsCellODEParam();


//---- OPERATORS


//---- ACCESS

  double getTime() const;
  double getSimulationTime() const;
  Array<double,1> getTimeMesh() const;
  int getITimeSlice() const;
  int getNTimeSlice() const;
  int getNTrajectory() const;
  int getITrajectory() const;
  int getNParameterSet() const;
  int getIParameterSet() const;
  int getNGlobalParameter() const;
  Array<double,1> getGlobalParameterArray() const;
  vector<string> getlistGlobalParametersName() const;
  Array<double,1> getCellsAngle() const;
  Array<TinyVector<double,3>,1> getCellsPosition() const;
  double getSpatialDynamicsEquilibrationTime() const;



//---- OPERATIONS

  void initializeCellCollection();

  void increaseTime(double t1);

  void setTime(double t1);

  void setIsEquilibrationSteps(bool isEquilibrationSteps);

  void computeSpatialIntegration(double timeStep);

  /**
   * Compute the evolution over a time slice of length #timeMeshSliceLength_. Output the state to the output object.
   * @param iSlice [in] Index of the time slice (defined in #timeMesh_),
   * 0 <= iSlice <= nTimeSlices-1.
   */
  void computeTrajectoryTimeSlice(int iSlice);

  /**
   * Compute a step of simulation, in general one time slice, in order to return the control to the GUI.
   */
  bool computeSimulationStep();

  /**
    * Write the position array into the last time slice. This method is used to update the position
    * of the cells in the time slice when the positions are updated from outside the simulator,
    * for example if the dynamics of the cells are computed in the GUI module.
    */
  void updatePositionLastTimeSlice();

  /**
    * Write the angle array into the last time slice. This method is used to update the angles
    * of the cells in the time slice when the angles are updated from outside the simulator,
    * for example if the dynamics of the cells are computed in the GUI module.
    */
  void updateAngleLastTimeSlice();

private:
  void initializeTimeMeshCellsTrajectoriesAvg();
  void initializeTimeMeshCellsAvg();
  // The following methods are designed in order to run correctly the #computeSimulationStep() method.
  void computeTrajectoryInit(int iTraj);
  void computeTrajectoryEnd(int iTraj);
  void computeNTrajectoriesInit();
  void computeNTrajectoriesEnd();
  void computeInit();
  void computeEnd();
  void checkTrajectoryStopCondition();
public:


//---- DATA

public:
  CellCollection cellCollection_;

  double proportionNegativeSteps_;


private:
  /**
   * The integrator context used for the simulation. The IntegratorContext class
   * is a template class with template type defining the type of integrator
   * (Gillespie, modified Gillespie, etc).
   */
  #ifndef USE_CHEMICAL_LANGEVIN
    #ifndef TIME_DEPENDENT_PROPENSITIES
      IntegratorContext<IntegratorGillespie> integratorContext_;
    #else
      IntegratorContext<IntegratorGillespieModified> integratorContext_;
    #endif
  #else
    IntegratorContext<IntegratorChemicalLangevin> integratorContext_;
  #endif

  /**
   * Integrate the spatial dynamics (cell movement).
   */
  SpatialIntegratorContext<SpatialIntegratorODE> spatialIntegratorContext_;
  //void integrateSpatialDynamics(double timeStep);

  void computeSpatialDynamicsEquilibration();

  bool isEnabledComputeSpatialDynamics_;

  double time_;

  Output output_;

  Input input_;

  /**
   * Length of time slices used to compute time mesh averages.
   */
  double timeMeshSliceLength_;

  /**
   * Time step used in the integration of the chemical langevin algorithm
   */
  double chemicalLangevinTimeStep_;

  /**
   * Number of time slices in the time mesh.
   */
  int nTimeSlice_;
  /**
   * Iterator over time slices.
   */
  int iTimeSlice_;

  /**
   * Total time of simulation. All the trajectories will be computed for time
   * 0 < t < simulationTime_ (with some inprecision due to the time mesh
   * definition).
   */
  double simulationTime_;

  /**
   * Number of trajectories to compute.
   */
  int nTrajectory_;
  /**
   * Iterator over trajectories.
   */
  int iTrajectory_;

  /**
   * Array containing the time mesh points.
   */
  blitz::Array<double,1> timeMesh_;

  /**
   * Object that contains the data of one time slice.
   */
  TimeSlice timeSlice_;

  /**
   * Object that contains the data of one time slice:
   * average of the state over all the cells.
   */
  TimeSliceCellsAvg timeSliceCellsAvg_;

  /**
   * List of the time slices data for one trajectory.
   */
  list<TimeSlice> timeMeshTrajectory_;

  /**
   * List of the time slices data for one trajectory:
   * average of the state over all the cells.
   */
  Array<TimeSliceCellsAvg,1> timeMeshCellsAvg_;

  /**
   * List of the time slices data averaged over all the trajectories:
   * average of the state over all the cells.
   */
  Array<TimeSliceCellsAvg,1> timeMeshCellsTrajectoriesAvg_;

  /**
   * Number of parameter sets (defined in input object).
   */
  int nParameterSet_;
  /**
   * Iterator over parameter sets.
   */
  int iParameterSet_;

  /**
   * Boolean indicating the end of the simulation.
   */
  bool simulationEnd_;

  /**
    * Boolean value used to stop the trajectory simulation loop.
    */
  bool stopTrajectoryFlag_;
  /**
    * First passage simulation mode: Index of the GFP species used in the condition to stop trajectory simulation.
    */
  int iGFPSpeciesIndexThreshold_;
  /**
    * First passage simulation mode: Threshold of GFP concentration used in the condition to stop trajectory simulation.
    */
  double GFPconcentrationThreshold_;
  /**
    * First passage simulation mode: direction of the first passage condition, passing GFP threshold from ON to OFF (true) or from OFF to ON (false)
    */
  bool isFirstPassageTimeDirectionOnToOff_;

  /**
    * For simulation of unlimited growth, stop simulation when the external volume gets to zero, i.e. when the cells fill all the total volume.
    */
  bool stopSimulationWhenReachingMaximumVolume_;
  /**
    * Array of size nTrajectories containing the list of first-passage times
    * (when the cell with index 0 passes the GFP concentration threshold).
    */
  Array<double,1> firstPassageTimeValues_;

  // Time variables for calculating the computational time of the simulation.
  // We use the ctime standard library requiring #include <time.h>
  time_t trajTimeStart_;
  time_t trajTimeEnd_;
  double trajTimeDiff_;
  time_t nTrajTimeStart_;
  time_t nTrajTimeEnd_;
  double nTrajTimeDiff_;
  time_t simTimeStart_;
  time_t simTimeEnd_;
  double simTimeDiff_;

};


inline double Simulator::getTime() const
  {return time_;}

inline double Simulator::getSimulationTime() const
  {return simulationTime_;}

inline void Simulator::increaseTime(double t1)
  {time_ += t1;}

inline void Simulator::setTime(double t1)
  {time_ = t1;}

inline int Simulator::getITimeSlice() const
  {return iTimeSlice_;}
inline int Simulator::getNTimeSlice() const
  {return nTimeSlice_;}

inline int Simulator::getNTrajectory() const
  {return nTrajectory_;}
inline int Simulator::getITrajectory() const
  {return iTrajectory_;}

inline int Simulator::getNParameterSet() const
  {return nParameterSet_;}
inline int Simulator::getIParameterSet() const
  {return iParameterSet_;}



#endif // SIMULATOR_H
