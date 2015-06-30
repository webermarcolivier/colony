/***************************************************************************//**
 * Project: Colony
 *
 * \file    CellMilieuChemicalSystem.h
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

#ifndef CELLMILIEUCHEMICALSYSTEM_H
#define CELLMILIEUCHEMICALSYSTEM_H

#include "debug.h"

// libraries header files
#include <blitz/array.h>
#include "Param/CellMilieuChemicalSystemInitParam.h"

// namespaces
using blitz::Array;
using blitz::Range;
using std::cout;
using std::ostream;
using std::endl;



/**
 * Use forward declaration of class Cell. The header file "Cell.h" is included
 * after the declaration of the CellMilieuChemicalSystem class.
 * This is necessary since the CellMilieuChemicalSystem class contains a pointer to
 * the Cell class.
 */
class Cell;
/**
 * Use forward declaration of class Milieu. The header file "Milieu.h" is included
 * after the declaration of the CellMilieuChemicalSystem class.
 * This is necessary since the CellMilieuChemicalSystem class contains a pointer to
 * the Milieu class.
 */
class Milieu;






/**
  * This Class defines the cell-milieu reactions like diffusion, signalling, etc. It
  * is contained in a Cell object and has a pointer to the milieu Cell object of the
  * colony. The definition of the reactions are given by the chemicalSystemCellMilieu_
  * attribute.
  */
class CellMilieuChemicalSystem
{
public:

  friend class Cell;

//---- LIFECYCLE

  CellMilieuChemicalSystem();

  virtual ~CellMilieuChemicalSystem();

  void initialize(const CellMilieuChemicalSystemInitParam& p);


//---- OPERATORS

  /**
   * Assignment operator.
   */
  CellMilieuChemicalSystem& operator = (const CellMilieuChemicalSystem& b);


//---- ACCESS

  int getNChannels() const;

  Array<double,1>& getPropensities();

//---- INQUIRY


//---- OPERATIONS

  /**
   * Set the propensities array propensities_ to be a reference to another array.
   */
  void setReferencePropensities(Array<double,1> refArray);

  /**
   * Compute the propensities of cell-milieu reaction channels.
   */
  void computePropensities();

  /**
   * Compute the propensities of cell-milieu reaction channels, some of them having
   * time-dependent reaction rates.
   */
  void computePropensities(const double time);

  /**
   * Apply the reaction given by the channel number. Change the state of the
   * cell and the milieu based on the stoichiometric matrix.
   */
  void applyReaction(int mu);

  /**
   * Check that the number of molecules in state vector x is positive.
   */
  void checkPositivity() const;


private:

//---- DATA

  /**
   * Pointer to the milieu Cell object.
   */
  Milieu* milieuPtr_;

  /**
   * Pointer to the cell Cell object.
   */
  Cell* cellPtr_;

  /**
   * nSpeciesCell + nSpeciesMilieu.
   */
  int nTotalSpecies_;

  int nSpeciesCell_;

  int nSpeciesMilieu_;

  /**
   * Number of number of cell-milieu reaction channels.
   */
  int nChannels_;

  /**
   * Stoichiometric matrix of the cell-milieu reactions. Its dimension should be
   * nChannels x (nSpeciesCell+nSpeciesMilieu)
   * The first 'nSpeciesCell' columns refer to the species in the cell, the
   * following 'nSpeciesMilieu' columns refer to the species in the milieu.
   */
  Array<int,2> stoichMatrix_;

  /**
   * Pointer to the function that computes the propensities of the different
   * reaction channels. The function is meant to be written by an external
   * program. Otherwise, it can be directly included in the ChemicalSystem
   * class.
   * @param x1 [in] vector of size #nSpeciesCell_ that contains the number of
   *        molecules of the species in the cell.
   * @param x2 [in] vector of size #nSpeciesMilieu_ that contains the number of
   *        molecules of the species in the milieu.
   * @param volume1 [in] volume of the cell. This is needed for computing the
   *        diffusion rate.
   * @param volume2 [in] volume of the milieu. This is needed for computing the
   *        diffusion rate.
   * @param propensities [in,out] vector of size #nChannels_ that contains the
   *        computed propensities.
   */
  void (*computePropensitiesInterCellularPtr_)
  (
    const Array<int,1>& x1,
    const Array<int,1>& x2,
    const double volume1,
    const double volume2,
    Array<double,1>& propensities
  );

  /**
   * Pointer to the function that computes the propensities of the different
   * reaction channels. The function is meant to be written by an external
   * program. Otherwise, it can be directly included in the ChemicalSystem
   * class.
   *
   * @see computePropensitiesTimeDependentDiffusion()
   *
   * @param tc [in] The duration of the cell cycle.
   * @param x1 [in] %State vector of the cell.
   * @param x2 [in] %State vector of the milieu.
   * @param V0 [in] Volume of the cell at the beginning of the cell cycle (constant).
   * @param Vext [in] Volume of the external milieu at time t.
   * @param t [in] Time t of the evaluation of the propensities.
   * @param t0 [in] Time of the previous cell division (time of cell's birth).
   * @param tc [in] Duration of the cell cycle.
   * @param a [out] Array of the computed propensities.
   */
  void (*computeTimeDependentPropensitiesInterCellularPtr_)
  (
    const Array<int,1>& x1,
    const Array<int,1>& x2,
    const double volume,
    const double volume0,
    const double volumeExt,
    Array<double,1>& a
  );

  /**
   * This propensities array is just a reference to a part of the total array
   * of the Cell object:
   * totalPropensities_( nChannels , nTotalChannels-1 )
   */
  Array<double,1> propensities_;

  /**
   * Preventing copy constructor use.
   */
  CellMilieuChemicalSystem(const CellMilieuChemicalSystem& c2);

  /**
   * Copy method used in the assigment operator and the copy constructor.
   */
  void copy(const CellMilieuChemicalSystem& b);


};


inline int CellMilieuChemicalSystem::getNChannels() const
  {return nChannels_;}

inline Array<double,1>& CellMilieuChemicalSystem::getPropensities()
  {return propensities_;}



#include "Cell.h"
#include "Milieu.h"

#endif // CELLMILIEUCHEMICALSYSTEM_H
