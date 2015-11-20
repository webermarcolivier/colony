/***************************************************************************//**
 * Project: Colony
 *
 * \file    Simulator.cpp
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

#include "Simulator.h"

//------------------------------------------------------------------------------

void Simulator::init()
{
  /// - Initialize parameters of simulation.
  timeMeshSliceLength_ = input_.sim.timeMeshSliceLength_;
  chemicalLangevinTimeStep_ = input_.sim.chemicalLangevinTimeStep_;
  simulationTime_      = input_.sim.simulationTime_;
  nTrajectory_         = input_.sim.nTrajectory_;
  iTrajectory_         = 0;
  nParameterSet_       = input_.getNParameterSet();
  iParameterSet_       = 0;
  simulationEnd_ = false;
  stopTrajectoryFlag_ = false;
  stopSimulationWhenReachingMaximumVolume_ = input_.sim.stopSimulationWhenReachingMaximumVolume_;
  isEnabledComputeSpatialDynamics_ = input_.sim.isEnabledComputeSpatialDynamics_;
  time_ = 0.0;

  output_.initialize(input_);

  nTimeSlice_ = int( simulationTime_ / timeMeshSliceLength_ );
  iTimeSlice_ = 0;

  /// - Initialize the time mesh.
  #ifndef FIRST_PASSAGE_TIME_COMPUTATION
    // In the computation of the first passage time, the total time of a trajectory
    // can be very large, and the timeMesh uses too much memory. We don't need
    // this variable since we only record the time of transition.

    timeMesh_.resize( nTimeSlice_+1 );

    int i;
    for (i=0; i<nTimeSlice_+1; i++)
    {

      timeMesh_(i) = i * timeMeshSliceLength_;
    }
  #endif

  // Find the index of the species GFP
  vector<string>::iterator it;
  it = find (input_.p.cellInitParam_.cellBaseInitParam_.stateInitParam_.listSpeciesName_.begin(),
             input_.p.cellInitParam_.cellBaseInitParam_.stateInitParam_.listSpeciesName_.end(),
             "Y");
  iGFPSpeciesIndexThreshold_ = 0;
  for ( ; it!=input_.p.cellInitParam_.cellBaseInitParam_.stateInitParam_.listSpeciesName_.begin(); --it)
  {
    iGFPSpeciesIndexThreshold_++;
  }

  GFPconcentrationThreshold_ = input_.sim.GFPconcentrationThreshold_;
  isFirstPassageTimeDirectionOnToOff_ = input_.sim.isFirstPassageTimeDirectionOnToOff_;
  firstPassageTimeValues_.resize(nTrajectory_);

  /// - Initialize the timeMeshCellsAvg_ and the timeMeshCellsTrajectoriesAvg_.
  //initializeTimeMeshCellsTrajectoriesAvg();
  //initializeTimeMeshCellsAvg();

  computeInit();
}

//------------------------------------------------------------------------------

void Simulator::initMainWindowParam(GraphicsCellCompositeParam& graphicsCellParam)
{
  /// - Pass the initialization parameters for GraphicsCell classes (pointer to the graphics scene, etc.
  input_.p.cellInitParam_.graphicsCellParam_.copy(graphicsCellParam);
}

//------------------------------------------------------------------------------

void Simulator::initGraphicsCellODEParam()
{
  /// - Pass the ODE world and space references as parameters to initialize correctly
  ///   the GraphicsCellODE class, as well as a pointer to the spatialIntegratorODE.
  input_.p.cellInitParam_.graphicsCellParam_.space_ = spatialIntegratorContext_.integrator_.getSpace();
  input_.p.cellInitParam_.graphicsCellParam_.world_ = spatialIntegratorContext_.integrator_.getWorld();
  input_.p.cellInitParam_.graphicsCellParam_.spatialIntegratorODE_ = &(spatialIntegratorContext_.integrator_);
}

//------------------------------------------------------------------------------

Simulator::Simulator(string inputFile)
  : integratorContext_(this),
    spatialIntegratorContext_(this),
    input_(inputFile)
{
  initGraphicsCellODEParam();
  init();
  srand ( time(NULL) );
}

//------------------------------------------------------------------------------

Simulator::Simulator(string inputFile, GraphicsCellCompositeParam& graphicsCellParam)
  : integratorContext_(this),
    spatialIntegratorContext_(this),
    input_(inputFile)
{
  initMainWindowParam(graphicsCellParam);
  initGraphicsCellODEParam();
  init();
  srand ( time(NULL) );
}

//------------------------------------------------------------------------------

Simulator::~Simulator()
{}

//------------------------------------------------------------------------------

void Simulator::initializeCellCollection()
{
  /// - Initialize the cell collection.
  cellCollection_.initializeCellCollection(input_);
}

//------------------------------------------------------------------------------

void Simulator::computeInit()
{
  cout << "\n###### SIMULATION START ###### " << endl;

  simTimeDiff_ = 0.0;
  time (&simTimeStart_);
  struct tm * timePrint;
  timePrint = localtime (&simTimeStart_);
  cout << "###### Current Time: " << asctime(timePrint);
  cout << "\n" << endl;

  computeNTrajectoriesInit();
  computeTrajectoryInit(0);
}

//------------------------------------------------------------------------------

void Simulator::computeNTrajectoriesInit()
{
  /// Parameters values are stored in a public static array of class Input.
  /// The computePropensities function has direct access to these values.

  input_.setIndexParameterSet(iParameterSet_);
  input_.setGlobalParameter(0.0);
  input_.setX0(iParameterSet_);

  // Initialize the output directory.
  output_.initialize(input_);

  // Initialize the cell collection.
  initializeCellCollection();

  /// - Initialize the timeMeshCellsTrajectoriesAvg_.
  #ifndef FIRST_PASSAGE_TIME_COMPUTATION
    // In the computation of the first passage time, the total time of a trajectory
    // can be very large, and the time mesh uses too much memory. We don't need
    // this variable since we only record the time of transition.
    initializeTimeMeshCellsTrajectoriesAvg();
  #endif

  // Initialize the first passage time array
  firstPassageTimeValues_ = -1.0;

  cout << "#### Computing parameter set #" << fixed << setw(4) << setfill(' ') << iParameterSet_+1
       << "/"  << setw(4) << setfill(' ') << nParameterSet_ << " ####" << endl;

  int i;
  vector<string>::iterator it;
  for ( i=0, it=Input::listGlobalParametersName.begin();
                 it!=Input::listGlobalParametersName.end(); it++,i++ )
  {
    cout << *it << " " << Input::globalParameter(i) << endl;
  }
  cout << endl;

  nTrajTimeDiff_ = 0.0;
  time (&nTrajTimeStart_);
  struct tm * timePrint;
  timePrint = localtime (&nTrajTimeStart_);
  cout << "#### Current Time: " << asctime(timePrint);
  cout << endl;
}

//------------------------------------------------------------------------------

void Simulator::computeTrajectoryInit(int iTraj)
{
  cout << "## Computing trajectory #" << setw(4) << setfill(' ') << iTraj+1
       << "/" << setw(4) << setfill(' ') << nTrajectory_ << " ##" << endl;
  trajTimeDiff_ = 0.0;
  time (&trajTimeStart_);
  struct tm * timePrint;
  timePrint = localtime (&trajTimeStart_);
  //cout << "## Current Time: " << asctime(timePrint);

  time_ = 0.0;

  /// - Initialize the cell collection.
  initializeCellCollection();

  /// - Clear timeMeshTrajectory_.
  timeMeshTrajectory_.clear();

  /// - Initialize timeMeshCellsAvg_.
  #ifndef FIRST_PASSAGE_TIME_COMPUTATION
    // In the computation of the first passage time, the total time of a trajectory
    // can be very large, and the time mesh uses too much memory. We don't need
    // this variable since we only record the time of transition.
    initializeTimeMeshCellsAvg();
  #endif

  #ifdef USE_CHEMICAL_LANGEVIN
    integratorContext_.integrator_.setFixedNCells(cellCollection_.getNCells());
  #endif

  // Reset stop trajectory flag
  stopTrajectoryFlag_ = false;

  /// - At the beginning of each trajectory, we integrate
  ///   the spatial dynamics for equilibration of the cells.
  if(isEnabledComputeSpatialDynamics_)
  {
    computeSpatialDynamicsEquilibration();
  }
  cout << cellCollection_[0].getTimeDependentVolume(0) << endl;
}

//------------------------------------------------------------------------------

bool Simulator::computeSimulationStep()
{
  if (iTimeSlice_ < nTimeSlice_)
  {
    if ( iTimeSlice_ % input_.out.nbTimeStepsIntervalConsoleOutput_ == 0)
    {
      #ifndef FIRST_PASSAGE_TIME_COMPUTATION
        cout << "time " <<  setw(10) << setfill(' ') << fixed << time_
             << " iTimeSlice " << setw(6) << setfill(' ') << fixed << iTimeSlice_+1
             << " / " << setw(6) << setfill(' ') << fixed << nTimeSlice_+1
             //<< "\n" << cellCollection_.getCellsPosition()
             //<< "external volume = " << cellCollection_.getTimeDependentVolumeMilieu(time_) << " stopSimulationWhenReachingMaximumVolume_ = " << stopSimulationWhenReachingMaximumVolume_
             << endl;
      #endif
    }
    // Update the global parameters.
    input_.setGlobalParameter(time_);
    // Compute one time slice.
    double time0 = time_;
    computeTrajectoryTimeSlice(iTimeSlice_);
    double timeStep = time_ - time0;
    // Integrate the spatial dynamics
    if(isEnabledComputeSpatialDynamics_)
    {
      if (timeStep > 0)
      {
        computeSpatialIntegration(timeStep);
        updatePositionLastTimeSlice();
      }
    }
    ++iTimeSlice_;
  }
  else
  {
    #ifndef FIRST_PASSAGE_TIME_COMPUTATION
        cout << "time " <<  setw(10) << setfill(' ') << fixed << time_
             << " iTimeSlice " << setw(6) << setfill(' ') << fixed << iTimeSlice_+1
             << " / " << setw(6) << setfill(' ') << fixed << nTimeSlice_+1
             << endl;
    #endif

    computeTrajectoryEnd(iTrajectory_);
    iTimeSlice_ = 0;

    if (iTrajectory_ < nTrajectory_ - 1)
    {
      ++iTrajectory_;
      computeTrajectoryInit(iTrajectory_);
    }
    else
    {
      computeNTrajectoriesEnd();
      iTrajectory_ = 0;

      if (iParameterSet_ < nParameterSet_-1)
      {
        ++iParameterSet_;
        computeNTrajectoriesInit();
        computeTrajectoryInit(iTrajectory_);
      }
      else
      {
        // END OF SIMULATION
        simulationEnd_ = true;
        iTrajectory_ = nTrajectory_ - 1;
        iTimeSlice_ = nTimeSlice_ - 1;
        time_ = simulationTime_;
        computeEnd();
      }

    }
  }

  return simulationEnd_;
}

//------------------------------------------------------------------------------

void Simulator::computeTrajectoryTimeSlice(int iSlice)
{

  /// - Integration loop over the time slice.
  #ifndef FIRST_PASSAGE_TIME_COMPUTATION
    // In the computation of the first passage time, the total time of a trajectory
    // can be very large, and the time mesh uses too much memory. We don't need
    // this variable since we only record the time of transition.

    //double timeStart = timeMesh_(iSlice);
    //double timeEnd   = timeMesh_(iSlice+1);
  #endif
  double timeStart = iTimeSlice_*timeMeshSliceLength_;
  double timeEnd   = (iTimeSlice_+1)*timeMeshSliceLength_;


  ///   - Initialize the timeSliceCellsAvg_ object.
  int nSpecies       = cellCollection_[0].getNSpecies();
  int nSpeciesMilieu = cellCollection_.getMilieu().getNSpecies();
  int nSpeciesTot    = nSpecies + nSpeciesMilieu;
  timeSliceCellsAvg_.initialize( timeStart, nSpeciesTot );

  ///   The division events define subslices (t0,t1) dividing the time slice.
  ///   We initialize t0 = timeMesh_(iSlice) and t1 = 0.
  double t0 = timeStart;
  double t1 = 0.0;


  /// - Loop over the subslices.
  while ( t1 < timeEnd )
  {
    ///   - t1 = min( timeEnd, timeNextDivision ). Note: if there is no cell growth/division,
    ///     the timeNextDivision is equal to infinite.
    t1 = min( timeEnd, cellCollection_.getTimeNextDivision() );

    ///   - Initialize the timeSlice_ object.
    int nCells = cellCollection_.getNCells();
    // Creating an array with first element = volume of the milieu, and
    // following elements are volume table of all cells
    Array<double,1> volumeArrayTemp(cellCollection_.getTimeDependentVolumes(time_).copy());
    Array<double,1> volumeArray( volumeArrayTemp.extent(blitz::firstDim) + 1 );
    volumeArray( blitz::Range(1,blitz::toEnd) ) = volumeArrayTemp;
    volumeArray(0) = cellCollection_.getTimeDependentVolumeMilieu(time_);

    Array<TinyVector<double,3>,1> positionArray(cellCollection_.getCellsPosition().copy());

    #ifndef USE_CHEMICAL_LANGEVIN
      timeSlice_.initialize(nCells, t0, volumeArray, cellCollection_.getGlobalX(),
                            cellCollection_.getXConc(t0), positionArray,
                            cellCollection_.getCellsAngle());
    #else
      timeSlice_.initialize(nCells, t0, volumeArray, cellCollection_.getGlobalX(),
                            cellCollection_.getGlobalXConc(), positionArray,
                            cellCollection_.getCellsAngle());
    #endif


    double previousTime = t0;
    double timeStep = 0.0;


    ///   - Loop over time in the subslice.
    while ( time_ < t1 )
    {      
      /*
        Remark: if there is a cell division event, input parameter t1 has
        been set to the exact division event time. In this case, condition of
        while-loop in next step will evaluate to false (because time == t1).
      */

      #ifndef USE_CHEMICAL_LANGEVIN
        //timeStep = time_ - previousTime;
        //previousTime = time_;
        integratorContext_.integrateOneStep();
      #else
        timeStep = min(t1 - time_, chemicalLangevinTimeStep_);
        integratorContext_.integrateOneStep(timeStep);
      #endif
    }

    cellCollection_.updateTimeDependentVolumeMilieu(time_);

    #ifndef FIRST_PASSAGE_TIME_COMPUTATION
      // Accumulate the state of the cells only for the normal computation mode.

      ///   - Insert the time slice data to the trajectory list.
      // Remark: we have to insert a copy of the time slice average because of the
      // stateAvg_ array. Otherwise it just insert a reference to this array.
      TimeSlice timeSliceCopy( timeSlice_ );
      timeMeshTrajectory_.push_back( timeSliceCopy );

      /// - Compute the average over the cells of the time slice.
      int k;
      double weightFactor = (t1 - t0) / (timeEnd - timeStart) / double(nCells);
      double weightFactorNCellsAvg = (t1 - t0) / (timeEnd - timeStart);

      /// - Accumulate the average number of cells.
      timeSliceCellsAvg_.nCellsAvg_ += weightFactorNCellsAvg * double(nCells);

      /// - Accumulate the state of the milieu.
      timeSliceCellsAvg_.stateAvg_( Range(0, nSpeciesMilieu-1) ) +=
        //weightFactor * nCells * timeSlice_.t0State_( Range(0, nSpeciesMilieu-1) );
        weightFactor * nCells * timeSlice_.xConc_( Range(0, nSpeciesMilieu-1) );

      for (k=0; k<nCells; k++)
      {
        ///   - Accumulate the average over the time slice for each cell in the
        ///     timeSliceCellsAvg_ object and increment the samples count by nCells.
        timeSliceCellsAvg_.stateAvg_( Range( nSpeciesMilieu,
                                             nSpeciesMilieu + nSpecies - 1) )
          += weightFactor*
                /*timeSlice_.t0State_( Range( nSpeciesMilieu +  k   *nSpecies,
                                            nSpeciesMilieu + (k+1)*nSpecies - 1) );*/
                timeSlice_.xConc_( Range( nSpeciesMilieu +  k   *nSpecies,
                                          nSpeciesMilieu + (k+1)*nSpecies - 1) );
      }

      checkTrajectoryStopCondition();

    #else
      // In the first-passage computation mode, do not accumulate the state of the cells,
      // only check for the trajectory stop condition and record time value.

      // - Check for trajectory stop condition
      checkTrajectoryStopCondition();

    #endif

    /// - t0 = t1.
    t0 = t1;

  }


  #ifndef FIRST_PASSAGE_TIME_COMPUTATION
    // Accumulate the state of the cells only for the normal computation mode.

    ///   - Insert the time slice data to the cells average array.
    ///     Uses the + operator of class TimeSliceCellsAvg.
    timeMeshCellsAvg_(iSlice).initialize( timeMesh_(iSlice), nSpeciesTot );
    timeMeshCellsAvg_(iSlice) + timeSliceCellsAvg_;

    ///   - Insert the time slice data to the cells-trajectories average array.
    ///     Uses the + operator of class TimeSliceCellsAvg.
    timeMeshCellsTrajectoriesAvg_(iSlice) + timeSliceCellsAvg_;

  #endif

  #ifndef USE_CHEMICAL_LANGEVIN
    /// Compute all the propensities in order to reduce the errors generated by
    /// the optimization method.
    #ifdef OPTIMIZE_COMPUTE_PROPENSITIES_UPDATE
      // Initialize the temporary sum of propensities.
      #ifdef TIME_DEPENDENT_PROPENSITIES
        cellCollection_.initializeSumPropensities(time_);
      #else
        cellCollection_.initializeSumPropensities();
      #endif
    #endif
  #endif // USE_CHEMICAL_LANGEVIN
}

//------------------------------------------------------------------------------

void Simulator::computeTrajectoryEnd(int iTraj)
{
  /// - Output the averaged trajectory.
  #ifdef WRITE_ALL_TRAJECTORIES
    int nTrajWrite = nTrajectory_;
  #else
    int nTrajWrite = 5;
  #endif
  #ifndef FIRST_PASSAGE_TIME_COMPUTATION
    if (iTraj < nTrajWrite)
    {
      output_.writeTimeMeshTrajectory(timeMeshTrajectory_,iTraj,isEnabledComputeSpatialDynamics_);
      #ifdef WRITE_CELLS_AVG
        output_.writeTimeMeshCellsAvg(timeMeshCellsAvg_,iTraj);
      #endif
      #ifdef WRITE_CELL_LINEAGE
        output_.writeCellLineage(cellCollection_.getCellLineage(),iTraj);
      #endif
    }
  #endif

  #ifndef FIRST_PASSAGE_TIME_COMPUTATION
    cout << "## End of trajectory ##" << endl;
  #endif
  time(&trajTimeEnd_);
  struct tm * timePrint;
  timePrint = localtime (&trajTimeEnd_);
  //cout << "## Current Time: " << asctime(timePrint);
  trajTimeDiff_ = difftime(trajTimeEnd_,trajTimeStart_);
  cout << "## Computing time: ";
  printTimeDifference(cout, trajTimeDiff_);
  cout << endl << endl;
}

//------------------------------------------------------------------------------

void Simulator::computeNTrajectoriesEnd()
{
  /// - Output the cells-averaged trajectory.
  #ifndef FIRST_PASSAGE_TIME_COMPUTATION
    #ifdef WRITE_CELLS_TRAJ_AVG
      output_.writeTimeMeshCellsTrajectoriesAvg(timeMeshCellsTrajectoriesAvg_);
    #endif
  #endif

  /// - Output the phases of the cells (to check the distribution of cell cycles.
  #ifdef WRITE_CELL_CYCLE_PHASES_TIME_END
    output_.writeCellCyclePhases( cellCollection_.getCellCyclePhases(time_) );
  #endif

  // Output the first passage time array
  #ifdef FIRST_PASSAGE_TIME_COMPUTATION
    output_.writeFirstPassageTime( firstPassageTimeValues_ );
  #endif

  output_.close();

  cout << "#### End of parameter set ####" << endl;
  time(&nTrajTimeEnd_);
  struct tm * timePrint;
  timePrint = localtime (&nTrajTimeEnd_);
  cout << "#### Current Time: " << asctime(timePrint);
  nTrajTimeDiff_ = difftime(nTrajTimeEnd_,nTrajTimeStart_);
  cout << "#### Computing time: ";
  printTimeDifference(cout, nTrajTimeDiff_);
  cout << endl << endl << endl;
}
//------------------------------------------------------------------------------

void Simulator::computeEnd()
{
  cout << "###### SIMULATION END ######" << endl;

  time(&simTimeEnd_);
  struct tm * timePrint;
  timePrint = localtime (&simTimeEnd_);
  cout << "###### Current Time: " << asctime(timePrint);
  simTimeDiff_ = difftime(simTimeEnd_,simTimeStart_);
  cout << "###### Computing time: ";
  printTimeDifference(cout, simTimeDiff_);
  cout << "\n\n" << endl;
  #ifdef USE_CHEMICAL_LANGEVIN
    cout << "Proportion of correction steps = " << proportionNegativeSteps_ << endl;
  #endif
}

//------------------------------------------------------------------------------

void Simulator::initializeTimeMeshCellsTrajectoriesAvg()
{
  timeMeshCellsTrajectoriesAvg_.resize(nTimeSlice_);
  int nSpecies       = cellCollection_[0].getX().size();
  int nSpeciesMilieu = cellCollection_.getMilieu().getNSpecies();
  int nSpeciesTot    = nSpecies + nSpeciesMilieu;
  int i;
  for (i=0; i<nTimeSlice_; i++)
  {
    timeMeshCellsTrajectoriesAvg_(i).initialize(timeMesh_(i), nSpeciesTot);
  }
}

//------------------------------------------------------------------------------

void Simulator::initializeTimeMeshCellsAvg()
{
  timeMeshCellsAvg_.resize(nTimeSlice_);
  int nSpecies       = cellCollection_[0].getX().size();
  int nSpeciesMilieu = cellCollection_.getMilieu().getNSpecies();
  int nSpeciesTot    = nSpecies + nSpeciesMilieu;
  int i;
  for (i=0; i<nTimeSlice_; i++)
  {
    timeMeshCellsAvg_(i).initialize(timeMesh_(i), nSpeciesTot);
  }
}

//------------------------------------------------------------------------------

void Simulator::checkTrajectoryStopCondition()
{
  // Do not retest the condition if the flag has been activated.
  // This is to ensure that we do not rewrite the time value.
  // The flag has to be reset at the beginning of each trajectory.
  if (!stopTrajectoryFlag_)
  {
    #ifdef FIRST_PASSAGE_TIME_COMPUTATION
      // Test if the concentration of GFP in cell #0 is above or under a fixed threshold (depend on the direction).
      if ( (cellCollection_[0].getXconc(time_)(iGFPSpeciesIndexThreshold_) > GFPconcentrationThreshold_) ^ isFirstPassageTimeDirectionOnToOff_)
      {
        // Set the flag to stop the trajectory
        stopTrajectoryFlag_ = true;
        // Set the time slice to the end of the trajectory
        iTimeSlice_ = nTimeSlice_;
        // Record the time value
        firstPassageTimeValues_( iTrajectory_ ) = time_;

        cout << "Cell has passed the threshold. Time = " << time_ << endl;
      }
    #else

      if ( ( cellCollection_.getTimeDependentVolumeMilieu(time_) <= 2.*cellCollection_[0].getVolume0() ) && stopSimulationWhenReachingMaximumVolume_)
      {
        cout << "Cells have reached maximum volume. Time = " << time_ << endl;
        // Set the flag to stop the trajectory
        stopTrajectoryFlag_ = true;

        // Set the time slice to the end of the trajectory
        iTimeSlice_ = nTimeSlice_;
      }
    #endif
  }

}

//------------------------------------------------------------------------------

Array<double,1> Simulator::getTimeMesh() const
{
  Array<double,1> timeMesh(timeMesh_.copy()  );
  return timeMesh;
}

//------------------------------------------------------------------------------

int Simulator::getNGlobalParameter() const
{
  return Input::globalParameter.extent(firstDim);
}

//------------------------------------------------------------------------------

Array<double,1> Simulator::getGlobalParameterArray() const
{
  return Input::globalParameter;
}

//------------------------------------------------------------------------------

vector<string> Simulator::getlistGlobalParametersName() const
{
  return Input::listGlobalParametersName;
}

//------------------------------------------------------------------------------

void Simulator::setIsEquilibrationSteps(bool isEquilibrationSteps)
{
  spatialIntegratorContext_.setIsEquilibrationSteps(isEquilibrationSteps);
}

//------------------------------------------------------------------------------

void Simulator::computeSpatialDynamicsEquilibration()
{
  cout << "ODE equilibration of colony... start" << endl;

  float refreshTime = 0.1;
  float spatialEquilibrationTime = getSpatialDynamicsEquilibrationTime();
  // When there is only 1 cell, the ODE integration somehow explodes, therefore
  // we just skip the equilibration integration. Anyway, with 1 cell there is no
  // need to move the cells prior to the simulation start.
  if (cellCollection_.getNCells() > 1) {
      for (float time = 0.0; time < getSpatialDynamicsEquilibrationTime();
           time+=refreshTime)
      {
        setIsEquilibrationSteps(true);
        computeSpatialIntegration(refreshTime);
        if ( int(time / refreshTime) % 20 == 0)
        {
          cout << "ODE equilibration of colony... " << time << "/"
               << getSpatialDynamicsEquilibrationTime() << endl;
        }
      }
      cout << "ODE equilibration of colony... end" << endl;
  }

  setIsEquilibrationSteps(false);
}

//------------------------------------------------------------------------------

void Simulator::computeSpatialIntegration(double timeStep)
{
  spatialIntegratorContext_.integrate(timeStep);
}


//------------------------------------------------------------------------------

Array<double,1> Simulator::getCellsAngle() const
{
  Array<double,1> cellsAngleArray( cellCollection_.getCellsAngle().copy() );
  return cellsAngleArray;
}

//------------------------------------------------------------------------------

Array<TinyVector<double,3>,1> Simulator::getCellsPosition() const
{
  Array<TinyVector<double,3>,1> cellsPositionArray( cellCollection_.getCellsPosition().copy() );
  return cellsPositionArray;
}

//------------------------------------------------------------------------------

double Simulator::getSpatialDynamicsEquilibrationTime() const
{
  return input_.ODEEquilibrationTime;
}

//------------------------------------------------------------------------------

void Simulator::updatePositionLastTimeSlice()
{
  if ( !timeMeshTrajectory_.empty() )
  {
    Array<TinyVector<double,3>,1> cellsPositionArray( getCellsPosition().copy() );
    timeMeshTrajectory_.back().setPosition( cellsPositionArray );
  }
}

//------------------------------------------------------------------------------

void Simulator::updateAngleLastTimeSlice()
{
  if ( !timeMeshTrajectory_.empty() )
  {
    Array<double,1> cellsAngleArray( getCellsAngle().copy() );
    timeMeshTrajectory_.back().setAngle( cellsAngleArray );
  }
}

//------------------------------------------------------------------------------
