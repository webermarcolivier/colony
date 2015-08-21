/***************************************************************************//**
 * Project: Colony
 *
 * \file    Input.cpp
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

#include "Input.h"

//------------------------------------------------------------------------------

Array<double,1> Input::globalParameter;

Array<double,3> Input::globalParameterTimePoints;

vector<string> Input::listGlobalParametersName;

double Input::ODEEquilibrationTime;

//------------------------------------------------------------------------------

Input::Input(string inputFile)
{
  // Set path of the input files.
  string path = "./";
  string filepath;
  const char* filename;

  // Read the input file.
  filepath = path + inputFile;
  filename = filepath.c_str();
  ifstream inputfile (filename);

  if (inputfile.is_open())
  {
    jumpNLines(inputfile,1);
    readOneVariable<int>   (inputfile,seed);
    RandomNumberGenerator::initialize(seed);

    jumpNLines(inputfile,1);
    readOneVariable<string>(inputfile,out.outputPath_);

    jumpNLines(inputfile,1);
    readOneVariable<double>(inputfile,sim.GFPconcentrationThreshold_);

    jumpNLines(inputfile,1);
    readOneVariable<bool>(inputfile,sim.isFirstPassageTimeDirectionOnToOff_);

    jumpNLines(inputfile,1);
    readOneVariable<bool>(inputfile,sim.stopSimulationWhenReachingMaximumVolume_);

    jumpNLines(inputfile,1);
    readOneVariable<double>(inputfile,sim.timeMeshSliceLength_);

    jumpNLines(inputfile,1);
    readOneVariable<double>(inputfile,sim.chemicalLangevinTimeStep_);

    jumpNLines(inputfile,1);
    readOneVariable<bool>(inputfile,sim.isEnabledComputeSpatialDynamics_);

    jumpNLines(inputfile,1);
    readOneVariable<double>(inputfile,ODEEquilibrationTime);

    jumpNLines(inputfile,1);
    readOneVariable<int>   (inputfile,sim.nTrajectory_);

    jumpNLines(inputfile,1);
    readOneVariable<int>   (inputfile,out.nbTimeStepsIntervalConsoleOutput_);

    jumpNLines(inputfile,1);
    inputfile >> simTimePoints_;

    jumpNLines(inputfile,1);

    // Set total simulation time to the last simulation time point.
    sim.simulationTime_ = simTimePoints_(simTimePoints_.ubound());

    // Global parameters
    jumpNLines(inputfile,1);
    readWordsInOneLine(inputfile,listGlobalParametersName);

    jumpNLines(inputfile,1);
    inputfile >> globalParameterTimePoints;
    jumpNLines(inputfile,1);

    // The array globalParameter contains the present value of the global
    // parameters and is used in the computePropensities function.
    // When reaching the simulation time points, the values in this array
    // are updated.
    globalParameter.resize( globalParameterTimePoints.extent(thirdDim) );
    globalParameter = globalParameterTimePoints( 0, 0, Range::all() );
    indexParameterSet_ = 0;
    nParameterSet_ = globalParameterTimePoints.extent(firstDim);

    // definition of the cells
    jumpNLines(inputfile,1);
    readOneVariable<int>   (inputfile,p.nCells_);

    jumpNLines(inputfile,1);
    readOneVariable<bool>  (inputfile,p.constantCellDensity_);

    jumpNLines(inputfile,1);
    readOneVariable<double>(inputfile,p.totalVolume_);

    jumpNLines(inputfile,1);
    readOneVariable<double>(inputfile,
      p.cellInitParam_.cellBaseInitParam_.stateInitParam_.cellVolume0_);

    jumpNLines(inputfile,1);
    readOneVariable<double>(inputfile,
      p.cellInitParam_.cellBaseInitParam_.stateInitParam_.cellLength0_);

    jumpNLines(inputfile,1);
    readOneVariable<double>(inputfile,
      p.cellInitParam_.cellBaseInitParam_.stateInitParam_.cellHeight0_);

    jumpNLines(inputfile,1);
    readOneVariable<double>(inputfile,
      p.cellInitParam_.cellCycleDeterministic_);

    jumpNLines(inputfile,1);
    readOneVariable<double>(inputfile,
      p.cellInitParam_.cellCycleGamma_);

    jumpNLines(inputfile,1);
    readWordsInOneLine(inputfile,
      p.cellInitParam_.cellBaseInitParam_.stateInitParam_.listSpeciesName_);

    jumpNLines(inputfile,1);
    readOneVariable<int>   (inputfile,
      p.cellInitParam_.cellBaseInitParam_.chemicalSystemInitParam_.nSpecies_);

    jumpNLines(inputfile,1);
    readOneVariable<int>   (inputfile,
      p.cellInitParam_.cellBaseInitParam_.chemicalSystemInitParam_.nChannels_);

    jumpNLines(inputfile,1);
    inputfile >> x0CellChangingWithGlobalParameter_;
    jumpNLines(inputfile,1);

    jumpNLines(inputfile,1);
    inputfile >> x0CellHeterogenous_;
    jumpNLines(inputfile,1);

    jumpNLines(inputfile,1);
    if (!x0CellHeterogenous_)
    {

      x0CellTable_.resize( nParameterSet_, p.cellInitParam_.cellBaseInitParam_.chemicalSystemInitParam_.nSpecies_ );

      if (x0CellChangingWithGlobalParameter_)
      {
        // Read all the initial conditions vectors x0Cell.
        int iPar;
        for (iPar=0; iPar<nParameterSet_; iPar++)
        {
          Array<int,1> dummyTable;
          inputfile >> dummyTable;
          jumpNLines(inputfile,1);
          x0CellTable_(iPar,blitz::Range::all()) = dummyTable;
        }
      } else {

        // Copy the unique x0 into the x0CellTable for all iPar.
        Array<int,1> dummyTable;
        inputfile >> dummyTable;
        int iPar;
        for (iPar=0; iPar<nParameterSet_; iPar++)
        {
          x0CellTable_(iPar,blitz::Range::all()) = dummyTable;
        }
        jumpNLines(inputfile,1);
      }

    } else {
      x0CellHeterogenousTable_.resize( nParameterSet_, p.nCells_,
                                       p.cellInitParam_.cellBaseInitParam_.chemicalSystemInitParam_.nSpecies_ );
      initialPhaseTable_.resize( nParameterSet_, p.nCells_ );

      if (x0CellChangingWithGlobalParameter_)
      {
        // Read the initial conditions arrays x0Cell
        int iPar;
        for (iPar=0; iPar<nParameterSet_; iPar++)
        {
          Array<int,2> dummyTable;
          inputfile >> dummyTable;
          jumpNLines(inputfile,1);
          x0CellHeterogenousTable_(iPar, blitz::Range::all(), blitz::Range::all()) = dummyTable;
        }
        // Read the initial conditions vectors initialPhase
        for (iPar=0; iPar<nParameterSet_; iPar++)
        {
          Array<double,1> dummyTable;
          inputfile >> dummyTable;
          jumpNLines(inputfile,1);
          initialPhaseTable_(iPar, blitz::Range::all()) = dummyTable;
        }

      } else {

        // Read the initial conditions arrays x0Cell
        {
          Array<int,2> dummyTable;
          inputfile >> dummyTable;
          int iPar;
          for (iPar=0; iPar<nParameterSet_; iPar++)
          {
            x0CellHeterogenousTable_(iPar, blitz::Range::all(), blitz::Range::all()) = dummyTable;
          }
          jumpNLines(inputfile,1);
        }

        // Read the initial conditions vectors initialPhase
        {
          Array<double,1> dummyTable;
          inputfile >> dummyTable;
          int iPar;
          for (iPar=0; iPar<nParameterSet_; iPar++)
          {
            initialPhaseTable_(iPar, blitz::Range::all()) = dummyTable;
          }
          jumpNLines(inputfile,1);
        }

      }
    }

    // definition of the milieu
    jumpNLines(inputfile,1);
    readWordsInOneLine(inputfile,
        p.milieuInitParam_.stateInitParam_.listSpeciesName_);

    jumpNLines(inputfile,1);
    readOneVariable<int>   (inputfile,
      p.milieuInitParam_.chemicalSystemInitParam_.nSpecies_);

    jumpNLines(inputfile,1);
    readOneVariable<int>   (inputfile,
      p.milieuInitParam_.chemicalSystemInitParam_.nChannels_);

    jumpNLines(inputfile,1);
    inputfile >> x0MilieuChangingWithGlobalParameter_;
    jumpNLines(inputfile,1);

    jumpNLines(inputfile,1);
    x0MilieuTable_.resize( nParameterSet_, p.milieuInitParam_.chemicalSystemInitParam_.nSpecies_ );
    if (x0MilieuChangingWithGlobalParameter_)
    {
      // Read all the initial conditions vectors x0Milieu.
      int iPar;
      for (iPar=0; iPar<nParameterSet_; iPar++)
      {
        Array<int,1> dummyTable;
        inputfile >> dummyTable;
        jumpNLines(inputfile,1);
        x0MilieuTable_(iPar,blitz::Range::all()) = dummyTable;
      }
    } else {
      // Copy the unique x0Milieu into the x0MilieuTable for all iPar.
      Array<int,1> dummyTable;
      inputfile >> dummyTable;
      int iPar;
      for (iPar=0; iPar<nParameterSet_; iPar++)
      {
        x0MilieuTable_(iPar,blitz::Range::all()) = dummyTable;
      }
      jumpNLines(inputfile,1);
    }

    setX0(0);

    jumpNLines(inputfile,1);
    inputfile >>
      p.cellInitParam_.dnaSpecies_;
    jumpNLines(inputfile,1);

    jumpNLines(inputfile,1);
    inputfile >>
      p.cellInitParam_.cellBaseInitParam_.chemicalSystemInitParam_.stoichMatrix_;
    jumpNLines(inputfile,1);

    #ifndef USE_CHEMICAL_LANGEVIN
      #ifndef TIME_DEPENDENT_PROPENSITIES
        p.cellInitParam_.cellBaseInitParam_.chemicalSystemInitParam_.computePropensitiesPtr_
          = computePropensitiesCell;
      #else
        p.cellInitParam_.cellBaseInitParam_.chemicalSystemInitParam_.computePropensitiesTimeDependentPtr_
          = computePropensitiesCellTimeDependent;
      #endif // TIME_DEPENDENT_PROPENSITIES
    #endif // USE_CHEMICAL_LANGEVIN


    // definition of the cell milieu reactions
    jumpNLines(inputfile,1);
    readOneVariable<int>   (inputfile,
    p.cellInitParam_.cellMilieuChemicalSystemInitParam_.nChannelsCellMilieu_);

    jumpNLines(inputfile,1);
    inputfile >>
      p.cellInitParam_.cellMilieuChemicalSystemInitParam_.stoichMatrixCellMilieu_;
    jumpNLines(inputfile,1);

    #ifndef USE_CHEMICAL_LANGEVIN
      #ifndef TIME_DEPENDENT_PROPENSITIES
        p.cellInitParam_.cellMilieuChemicalSystemInitParam_.computePropensitiesInterCellularPtr_
          = computePropensitiesCellMilieu;
      #else
        p.cellInitParam_.cellMilieuChemicalSystemInitParam_.computeTimeDependentPropensitiesInterCellularPtr_
          = computePropensitiesTimeDependentDiffusion;
      #endif // TIME_DEPENDENT_PROPENSITIES
    #endif // USE_CHEMICAL_LANGEVIN

    jumpNLines(inputfile,1);
    inputfile >>
      p.milieuInitParam_.chemicalSystemInitParam_.stoichMatrix_;
    jumpNLines(inputfile,1);

    #ifndef USE_CHEMICAL_LANGEVIN
      #ifndef TIME_DEPENDENT_PROPENSITIES
        p.milieuInitParam_.chemicalSystemInitParam_.computePropensitiesPtr_
          = computePropensitiesMilieu;
      #else
        p.milieuInitParam_.chemicalSystemInitParam_.computePropensitiesTimeDependentPtr_
          = computePropensitiesMilieuTimeDependent;
      #endif // TIME_DEPENDENT_PROPENSITIES
    #endif // USE_CHEMICAL_LANGEVIN
cout << "finish reading input" << endl;
    inputfile.close();

  } else {
    cout << "ERROR: Unable to open file \"" << filename << '\"' << endl;
  }

}

//------------------------------------------------------------------------------

Input::~Input()
{}

//------------------------------------------------------------------------------

void Input::setGlobalParameter(const double time)
{
  int it = 0;
  while ( simTimePoints_(it) < time )
  {
    it++;
  }

  globalParameter = globalParameterTimePoints( indexParameterSet_, it, Range::all() );
}

//------------------------------------------------------------------------------

void Input::setX0(const int iParameterSet)
{
  p.milieuInitParam_.stateInitParam_.x_.resize(x0MilieuTable_.extent(secondDim));
  p.milieuInitParam_.stateInitParam_.x_ =
                                 x0MilieuTable_( iParameterSet, blitz::Range::all() );

  if (!x0CellHeterogenous_)
  {
    p.x0CellHeterogenous_ = false;
    p.cellInitParam_.cellBaseInitParam_.stateInitParam_.x_.resize(x0CellTable_.extent(secondDim));
    p.cellInitParam_.cellBaseInitParam_.stateInitParam_.x_ =
                                   x0CellTable_( iParameterSet, blitz::Range::all() );
  } else {
    p.x0CellHeterogenous_ = true;

    p.x0CellHeterogenousTable_.resize(p.nCells_,
                                      p.cellInitParam_.cellBaseInitParam_.chemicalSystemInitParam_.nSpecies_ );
    p.x0CellHeterogenousTable_ = x0CellHeterogenousTable_( iParameterSet, blitz::Range::all(), blitz::Range::all() );
    p.cellInitParam_.cellBaseInitParam_.stateInitParam_.x_.resize
        (p.cellInitParam_.cellBaseInitParam_.chemicalSystemInitParam_.nSpecies_);

    p.initialPhaseTable_.resize(p.nCells_);
    p.initialPhaseTable_ = initialPhaseTable_( iParameterSet, blitz::Range::all() );
  }
}

//------------------------------------------------------------------------------

void Input::setIndexParameterSet(const int i)
{
  indexParameterSet_ = i;
}

//------------------------------------------------------------------------------

int Input::getNParameterSet() const
{
  return nParameterSet_;
}

//------------------------------------------------------------------------------



