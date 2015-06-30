/***************************************************************************//**
 * Project: Colony
 *
 * \file    apptest.cpp
 * \author  Marc Weber\n
 *          The Si.M.Bio.Sys. Group (CosmoLab)\n
 *          Parc Científic de Barcelona\n
 *          Barcelona, Spain.\n
 *          http://thesimbiosys.isgreat.org
 * \version 0.1
 * \date    11/2009
 *
 *          Copyright 2009 by Marc Weber
 ******************************************************************************/

#include "debug.h"

// standard C++ header files
#include <iostream>
#include <iomanip>
#include <utility>

// libraries header files
#include <blitz/array.h>

// user header files
#include "Simulator.h"

// namespaces
using blitz::Array;
using blitz::Range;
using std::cout;
using std::ostream;
using std::endl;
using std::setw;



int main()
{

//------------------------------------------------------------------------------

//#define DEBUG_CELLDIVISION

#ifdef DEBUG_CELLDIVISION
  cout << "Creating cell collection" << endl;

  CellCollection cellCollection;

  Input input;

  cellCollection.initializeHomogeneousCellCollection(input);

  cout << cellCollection;

  cout << "### DUPLICATING CELL#0" << endl;
  cellCollection.duplicateCell(0);
  cout << cellCollection;

  cout << "### DELETING CELL#0" << endl;
  cellCollection.deleteCell(0);
  cout << cellCollection;

#endif //DEBUG_CELLDIVISION


//------------------------------------------------------------------------------

//#define DEBUG_PROPENSITIES

#ifdef DEBUG_PROPENSITIES
  cout << "Creating cell collection" << endl;

  CellCollection cellCollection;

  Input input;

  cellCollection.initializeHomogeneousCellCollection(input);

  cout << cellCollection;

  cout << "#### MODIFYING ELEMENT IN THE PROPENSITIES ARRAY OF"
          " CELL-MILIEU IN CELL#0 (VALUE=777):" << endl;
  cellCollection[0].DEBUG_MODIFY1(777);

  cout << "#### PRINTING PROPENSITIES ARRAYS:" << endl;
  cout << "cell 0 P: \n" << cellCollection[0].getPropensities() << endl;
  cout << "cell 0 CellBase P: \n" << cellCollection[0].CellBase::getPropensities() << endl;
  cout << "cell 0 Cell Milieu chemical system P: \n" << cellCollection[0].getP_CM() << endl;
  cout << "milieu P: \n" << cellCollection.getMilieu().getPropensities() << endl;
  cout << "cell 1 P: \n" << cellCollection[1].getPropensities() << endl;
  cout << "global array P: \n" << cellCollection.getGlobalPropensities() << endl;

  cout << "### DUPLICATING CELL#0" << endl;
  cellCollection.duplicateCell(0);
  cout << "#### PRINTING PROPENSITIES ARRAYS:" << endl;
  cout << "global array P: \n" << cellCollection.getGlobalPropensities() << endl;

  cout << "#### MODIFYING ELEMENT IN THE PROPENSITIES ARRAY OF"
          " CELL-MILIEU IN CELL#0 (VALUE=888):" << endl;
  cellCollection[0].DEBUG_MODIFY1(888);

  cout << "#### PRINTING PROPENSITIES ARRAYS:" << endl;
  cout << "cell 0 P: \n" << cellCollection[0].getPropensities() << endl;
  cout << "cell 0 CellBase P: \n" << cellCollection[0].CellBase::getPropensities() << endl;
  cout << "cell 0 Cell Milieu chemical system P: \n" << cellCollection[0].getP_CM() << endl;
  cout << "milieu P: \n" << cellCollection.getMilieu().getPropensities() << endl;
  cout << "cell 1 P: \n" << cellCollection[1].getPropensities() << endl;
  cout << "global array P: \n" << cellCollection.getGlobalPropensities() << endl;

#endif //DEBUG_PROPENSITIES

//------------------------------------------------------------------------------

//#define DEBUG_APPLY_REACTION

#ifdef DEBUG_APPLY_REACTION

  CellCollection cellCollection;

  Input input;

  cellCollection.initializeHomogeneousCellCollection(input);

  cout << cellCollection;

  cout << "#### APPLY DIFFUSION REACTION" << endl;
  cellCollection.applyReaction(2+52+0);
  cout << cellCollection;
  cout << "#### APPLY REVERSE DIFFUSION REACTION" << endl;
  cellCollection.applyReaction(2+52+1); //apply reverse reaction
  cout << cellCollection;

  cout << "#### APPLY MILIEU REACTION: DEGRADATION OF AI" << endl;
  cellCollection.applyReaction(0);
  cout << cellCollection;
  cellCollection.applyReaction(1); //apply reverse reaction

  cout << "#### APPLY CELL REACTION: DIMERIZATION OF AL" << endl;
  cellCollection.applyReaction(2+1);
  cout << cellCollection;
  cellCollection.applyReaction(2+2); //apply reverse reaction

#endif // DEBUG_APPLY_REACTION

//------------------------------------------------------------------------------

//#define DEBUG_COMPUTEPROPENSITIES

#ifdef DEBUG_COMPUTEPROPENSITIES

  /* Remark: there is no test condition for the positivity of the number of
   * molecules or the propensities. If we have only second order reaction, it
   * should be ok. When x = 0, we have x*(x-1) = 0*(-1) = 0. But in the case
   * that we have 3rd-order reactions, then x*(x-1)*(x-2) = 0*-1*-2 = 2. !!!!!
   * even if propensity should be exactly 0.
   * Also, CONDITION ABOUT POSITIVITY OF NUMBER OF MOLECULES SHOULD BE CHECKED.
   */

  CellCollection cellCollection;

  Input input;

  cellCollection.initializeHomogeneousCellCollection(input);

  cout << cellCollection;

  cout << "#### COMPUTE GLOBAL PROPENSITIES" << endl;
  cout << "Global propensities = " << endl;
  cout << cellCollection.computeAndGetPropensities() << endl;
  cout << cellCollection.getMilieu().getPropensities() << endl;
  cout << cellCollection[0].getPropensities() << endl;

#endif // DEBUG_COMPUTEPROPENSITIES

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//#define DEBUG_SIMULATOR

#ifdef DEBUG_SIMULATOR

  Simulator simulator;

  cout << simulator.cellCollection_;

  int i;
  for (i=0; i<20; i++)
  {
    simulator.integratorContext_.integrateOneStep();
    cout << "simulator time = " << simulator.time_ << endl;
    cout << simulator.cellCollection_;
  }

#endif // DEBUG_SIMULATOR

//------------------------------------------------------------------------------

//#define DEBUG_SIMULATOR2

#ifdef DEBUG_SIMULATOR2

  Simulator simulator;

  cout << simulator.cellCollection_;

  while (simulator.time_ < 300.0)
  {
    simulator.integratorContext_.integrateOneStep();
  }

  cout << simulator.cellCollection_ << endl;
  cout << "Propensities = \n";
  cout << simulator.cellCollection_.computeAndGetPropensities() << endl;


#endif // DEBUG_SIMULATOR2

//------------------------------------------------------------------------------

//#define DEBUG_SIMULATOR4

#ifdef DEBUG_SIMULATOR4

  Simulator simulator;

  cout << simulator.cellCollection_;

  simulator.computeTrajectory(1);

  cout << simulator.cellCollection_ << endl;

#endif // DEBUG_SIMULATOR4

//------------------------------------------------------------------------------

#define DEBUG_SIMULATOR_NTRAJ

#ifdef DEBUG_SIMULATOR_NTRAJ

  Simulator simulator;

  simulator.computeNTrajectories();

#endif // DEBUG_SIMULATOR_NTRAJ

//------------------------------------------------------------------------------

  return 0;
}
