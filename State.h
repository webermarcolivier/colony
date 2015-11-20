/***************************************************************************//**
 * Project: Colony
 *
 * \file    State.h
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

#ifndef STATE_H
#define STATE_H

#include "compilation_options.h"

// standard C++ header files
#include <limits>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <algorithm>

// libraries header files
#include <blitz/array.h>

#include "Param/StateInitParam.h"

// namespaces
using blitz::Array;
using blitz::TinyVector;
using std::numeric_limits;
using std::vector;
using std::string;
using std::cout;
using std::ostream;
using std::endl;
using std::setw;
using std::max;


/**
 * %State of a cell.
 *
 * This class defines all the variables for the state of a cell.
 *
 * Remark 1: These variables will be accessed and updated frequently by the
 * integrator algorithm.
 *
 * Remark 2: All the variable are public. Since this class is intended to be used
 * as a data structure for another class, it is easier to make the attributes
 * public.
*/

class State
{
public:

//---- LIFECYCLE

  State ();

  ~State ();

  /**
   * Remark: The timePreviousDivision_ and cellCycleDuration_ variables are
   * not initialized by this function.
   * There initial values are set to -infinity and +infinity respectively,
   * and can be changed when initializing
   * an object of the Cell class (they cannot be changed from CellBase).
   * Note: Then these variables should be put in the Cell class, that would be
   * the best way to do it...
   */
  void initialize(const StateInitParam& p);


//---- OPERATORS

  /**
   * Assignment operator.
   */
  State& operator = (const State& s2);

  /**
   * Operator that outputs the content of a Cell state. Remark: it is declared
   * as a friend operator, so that it has access to the private part of the
   * State class.
   */
  friend ostream& operator<< (ostream& out, const State& s);


//---- ACCESS

  double getTimeDependentVolume(const double time) const;

  TinyVector<double,3> getTimeDependentCellDimensions(const double time);


public:

//---- DATA

// Rem: All data are PUBLIC.

  /**
   * Vector containing the number of molecules of each chemical species.
   */
  Array<int,1> x_;

  /**
   * Vector containing the concentration of each chemical species.
   */
  Array<double,1> xConc_;

  /**
   * List of the species names.
   */
  vector<string> listSpeciesName_;

  /**
   * Vector containing the propensities of each reaction channels, depends on
   * the chemical system.
   */
  Array<double,1> propensities_;

  TinyVector<double,3> position_;
  double angle_;
  /**
   * Volume of the cell at the beginning of the cell cycle.
   */
  double cellVolume0_;
  double cellLength_;
  double cellLength0_;
  double cellHeight_;
  double cellHeight0_;

  /**
   * Absolute time of the previous division event (birth of the cell).
   */
  double timePreviousDivision_;

  /**
   * Duration of the current cell cycle. The cell cycle can vary from on cycle
   * to another, thus this variable is updated after the cell divison.
   */
  double cellCycleDuration_;


private:

  /**
   * Preventing copy constructor use.
   */
   State(const State& state2);

  /**
   * Copy method used in the assigment operator and the copy constructor.
   */
  void copy(const State& s2);

};

#endif // STATE_H
