/***************************************************************************//**
 * Project: Colony
 *
 * \file    CellBase.h
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

#ifndef CELLBASE_H
#define CELLBASE_H

#include "debug.h"

// standard C++ header files
#include <iostream>   // input/output interface
#include <vector>
#include <string>

// libraries header files
#include <blitz/array.h>

// user header files
#include "State.h"
#include "ChemicalSystem.h"
#include "RandomNumberGenerator.h"
#include "Param/CellBaseInitParam.h"

// namespaces
using blitz::Array;
using blitz::Range;
using blitz::firstDim;
using blitz::secondDim;
using blitz::TinyVector;
using std::cout;
using std::ostream;
using std::endl;
using std::vector;
using std::string;


/**
  * Base class for cell.
  *
  * This class defines the basic Cell objects. They contains an internal state
  * and a chemical system, as well as methods to compute the propensities and
  * apply a reaction.
  */

class CellBase
{
public:

//---- LIFECYCLE

  /**
   * Default constructor is protected.
   * CellBase objects cannot be instantiated directly by the user. The protected
   * access assures that only the subclasses can instanciate.
   */
  protected:  CellBase ( );
  public:

  virtual ~CellBase ( );

  /**
   * Initialize the cell and define all its members.
   */
  void initialize (CellBaseInitParam& p);


//---- OPERATORS

  /**
   * Operator that outputs the content of a Cell state. Remark: it is declared
   * as a friend operator, so that it has access to the private part of the
   * CellBase class.
   */
  friend ostream& operator<< (ostream& out, const CellBase& c);


//---- ACCESS

  int getNChannels() const;

  int getNSpecies() const;

  double getVolume0() const;

  double getCellLength0() const;

  double getCellHeight0() const;

  virtual double getTimeDependentVolume(const double time) const = 0;

  TinyVector<double,3> getTimeDependentCellDimensions(const double time);

  Array<int,1>& getX();

  Array<double,1>& getXConc();

  Array<double,1> getXconc(const double time) const;

  Array<double,1>& getPropensities();

  double getTimeNextDivision() const;

  double getTimePreviousDivision() const;

  double getCellCycleDuration() const;

  double getCellCyclePhase(double time) const;

  double getAngle() const;

  TinyVector<double,3> getPosition() const;

  vector<string> getListSpeciesName() const;


//---- INQUIRY


//---- OPERATIONS

  void changeVolume0By(double volume0Change);

  /**
   * Apply the reaction given by the channel number. Change the state of the cell
   * based on the stoichiometric matrix.
   */
  void applyReaction(int mu);

  /**
   * Apply a chemical reaction that is defined outside the Cell class. This
   * function takes as input the stoichVector of the reaction and apply it to
   * the state of the cell.
   */
  void applyInterCellularReaction(Array<int,1> stoichVector);

  /**
   * Compute the propensities of all reaction channels.
   */
  void computePropensities();

  /**
   * Compute the propensities of all reaction channels, some of them having
   * time-dependent reaction rates.
   */
  void computePropensities(const double time);

  /**
   * Compute the propensities ONLY of the time-dependent reaction channels.
   */
  void computePropensitiesTimeDependentReactions(const double time);

  /**
   * Set the state array #x_ to be a reference to another array.
   */
  void setReferenceX(Array<int,1>& refArray);

  /**
   * Set the state array #xConc_ to be a reference to another array.
   */
  void setReferenceXConc(Array<double,1>& refArray);

  /**
   * Set the propensities array #propensities_ to be a reference to another array.
   */
  void setReferencePropensities(Array<double,1> refArray);

  /**
   * Check that the number of molecules in state vector x is positive.
   */
  void checkPositivity() const;

  /**
   * Change the angle of the cell.
   */
  void setAngle(double angle);

  /**
   * Change the position of the cell.
   */
  virtual void setPosition(TinyVector<double,3> position);



protected:

//---- DATA

  /**
   * Internal state of the cell.
   * Remark: derived classes have access to this attribute.
   */
  State state_;

  /**
   * Chemical system of the cell.
   * Remark: derived classes have access to this attribute.
   */
  ChemicalSystem chemicalSystem_;

  /**
   * Copy method used in the assigment operator and the copy constructor.
   * This copy() method is not used directly in the CellBase class but it is
   * called in its derived class Cell.
   */
  void copy(const CellBase& c2);


private:

  /**
   * Preventing copy constructor use.
   */
  CellBase(const CellBase& c2);


};

//------------------------------------------------------------------------------

// INLINE ACCESS METHODS

inline int CellBase::getNChannels() const
  {return chemicalSystem_.nChannels_;}

inline int CellBase::getNSpecies() const
  {return chemicalSystem_.nSpecies_;}

inline double CellBase::getVolume0() const
  {return state_.cellVolume0_;}

inline Array<int,1>& CellBase::getX()
  {return state_.x_;}

inline Array<double,1>& CellBase::getXConc()
  {return state_.xConc_;}

inline Array<double,1>& CellBase::getPropensities()
  {return state_.propensities_;}

inline double CellBase::getTimeNextDivision() const
  {return state_.timePreviousDivision_ + state_.cellCycleDuration_;}

inline double CellBase::getTimePreviousDivision() const
  {return state_.timePreviousDivision_;}

inline double CellBase::getCellCycleDuration() const
  {return state_.cellCycleDuration_;}

inline double CellBase::getCellCyclePhase(double time) const
  {return (time - state_.timePreviousDivision_) / state_.cellCycleDuration_;}

inline double CellBase::getAngle() const
  {return state_.angle_;}

inline TinyVector<double,3> CellBase::getPosition() const
  {return state_.position_;}

inline double CellBase::getCellLength0() const
  {return state_.cellLength0_;}

inline double CellBase::getCellHeight0() const
  {return state_.cellHeight0_;}

#endif // CELLBASE_H
