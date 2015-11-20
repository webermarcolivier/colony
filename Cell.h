/***************************************************************************//**
 * Project: Colony
 *
 * \file    Cell.h
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

#ifndef CELL_H
#define CELL_H

#include "compilation_options.h"

// standard C++ header files
#include <vector>
#include <string>

// libraries header files
#include <blitz/array.h>

// user header files
#include "CellBase.h"
#include "GraphicsCellComposite.h"
#include "CellMilieuChemicalSystem.h"
#include "RandomNumberGenerator.h"
#include "Param/CellInitParam.h"

// namespaces
using blitz::Array;
using blitz::Range;
using std::cout;
using std::ostream;
using std::endl;
using std::vector;
using std::string;



/**
  * %Cell able to divide and communicate with the milieu.
  *
  * This class is derived from CellBase and defines the cells able to live in a
  * milieu and divide.
  *
  * Remark: it has a propensities array of size nTotalChannels that contains the
  * propensities of intracellular reactions and cell-milieu reactions.
  */

class Cell : public CellBase, public GraphicsCellComposite
{
public:

//---- LIFECYCLE

  Cell();

  virtual ~Cell();

  /**
   * Copy constructor. We need the copy constructor in order to use a Blitz++
   * array of Cell objects (unkown details of the implementation of the
   * Blitz++ library). It is also used in the #duplicate() method to return
   * by value a copy of a Cell object.
   */
   Cell(const Cell& c2);

  /**
   * Initialize and define a Cell object, including dnaSpecies array
   * and cell-milieu reactions.
   *
   * The propensities array will have the size nChannels + nChannelsMilieu.
   * The propensities array for internal reactions (from CellBase class) is
   * referenced to the subarray #totalPropensities_(0,nChannels-1) and the prop.
   * array for cell-milieu reactions is referenced to the subarray
   * #totalPropensities_(nChannels,nTotalChannels).
   */
  void initialize(CellInitParam& p);



//---- OPERATORS

  /**
   * Assignment operator.
   */
  Cell& operator = (const Cell& c2);

  /**
   * Operator that outputs the content of a Cell state. Remark: it is declared
   * as a friend operator, so that it has access to the private part of the
   * Cell class. Overloading the << operator of class CellBase.
   */
  friend ostream& operator<< (ostream& out, const Cell& c);


//---- ACCESS

  double getTimeDependentVolume(const double time) const;

  /**
   * Get the total number of Cell objects.
   */
  static int getNCellObjects() {return nCellObjects_;}

  /**
   * Returns the total number of channels (internal + cell-milieu).
   */
  int getNChannels() const;

  /**
   * Returns the total array of propensities (internal + cell-milieu).
   */
  Array<double,1>& getPropensities();

  /**
   * Returns the array of propensities of the cell reactions.
   */
  Array<double,1>& getPropensitiesCell();

  /**
   * Returns the array of propensities of the cell-milieu reactions.
   */
  Array<double,1>& getPropensitiesCellMilieu();

  double getAngle() const;


//---- INQUIRY

  bool isReactionInternal(int mu) const;


//---- OPERATIONS

  /**
   * Set the propensities array to be a reference to another array.
   * This method set the #totalPropensities_ array to be a reference to
   * another array and it also set the reference of the sub arrays in the
   * CellBase class CellBase#state_.propensities_  and in the CellMilieuChemicalSystem
   * class #cellMilieuChemicalSystem_.propensities_ .
   */
  void setReferencePropensities(Array<double,1> refArray);


  /**
   * Apply the reaction given by the channel number. Change the state of the cell
   * based on the stoichiometric matrix.
   * There are two cases depending on the channel number: it can concern an internal
   * reaction (CellBase::applyReaction()) or a cell-milieu reaction
   * (cellMilieuChemicalSystem_.applyReaction()).
   * @param mu [in] channel number.
   */
  void applyReaction(int mu);

  /**
   * Compute the propensities of all reaction channels, including the internal and
   * cell-milieu reactions.
   */
  void computePropensities();

  /**
   * Compute the propensities of all reaction channels, including the internal and
   * cell-milieu reactions, some of them having time-dependent reaction rates.
   */
  void computePropensities(double time);

  /**
   * Compute the propensities <b>only</b> of the intracellular reaction channels.
   */
  void computePropensitiesCell();
  void computePropensitiesCell(const double time);

  /**
   * Compute the propensities <b>only</b> of the cell-milieu reaction channels.
   */
  void computePropensitiesCellMilieu();
  void computePropensitiesCellMilieu(const double time);

  /**
   * Duplicate the cell and returns by value the daughter cell.
   * - Detach proteins from DNA.
   * - Make a copy of the cell.
   * - Distribute proteins between the two cells.
   * - Set the time of next division for both cells.
   * - Return by value the new daughter cell (needs Cell copy constructor).
   */
  Cell duplicate();

  /**
    * Set the cell index in the GraphicsCellComposite class. The cell index should be
    * directly used in the Cell class, it is just used by the graphical interface to know
    * which cell it is when interacting with it.
    */
  void setCellIndex(int i);

  void setAngle(double angle);

  /**
    * Set the position of the cell, both in CellBase and GraphicsCellComposite classes.
    */
  virtual void setPosition(TinyVector<double,3> position);

  /**
    * Copy the cell dimensions and position from the CellBase::state_ to the GraphicsCellComposite class.
    */
  void updateGraphicsCell(const double time);


private:

//---- OPERATIONS

  /**
   * Detach the proteins from DNA. The species that contains genetic material are
   * released in solution. Species containing genetic material are defined in
   * the #dnaSpecies_ array.
   */
  void detachProteinsFromDna();


//---- DATA

  /**
   * Total number of Cell objects.
   */
  static int nCellObjects_;

  /**
   * Deterministic value for the cell cycle distribution.
   * @see RandomNumberGenerator::getCellCycle()
   */
  double cellCycleDeterministic_;
  /**
   * Stochastic to deterministic weight parameter for the cell cycle distribution.
   * @see RandomNumberGenerator::getCellCycle()
   */
  double cellCycleGamma_;

  /**
   * Set of cell-milieu reactions (for example, diffusion processes,
   * signalling, etc).
   */
  CellMilieuChemicalSystem cellMilieuChemicalSystem_;

  /**
   * Total number of reactions channels, including internal reactions and
   * cell-milieu reactions.
   */
  int nTotalChannels_;

  /**
   * Array containing the propensities of all reactions channels, including
   * internal reactions and cell-milieu reactions. The corresponding
   * propensities arrays for internal and cell-milieu reactions are referenced
   * to a part of totalPropensities_.
   */
  Array<double,1> totalPropensities_;

  /**
   *  Array defining genetic species.
   *  int array of size nSpecies x nSpecies.\n
   *  This array defines which species are genetic material and which are
   *  proteins. It is used to define which species are halved and which species
   *  are copied during a division event. It is defined the following way:
   *
   *  On the diagonal, for \c i == \c j :
      <table>
      <tr><td> \code dnaSpecies_(i,i) = 0 \endcode
      <td>  if species \c i is a protein => will be halved during division event. </tr>

      <tr><td> \code dnaSpecies_(i,i) = 1 \endcode
      <td> if species \c i is a protein-DNA complex => will be
           unbound before division event. </tr>

      <tr><td> \code dnaSpecies_(i,i) = -1 \endcode
      <td> if species \c i is a DNA molecule => will be
           copied during division event. </tr>
      </table>

      In the bulk matrix, for \c i != \c j :
      <table>
      <tr><td> \code dnaSpecies_(i,j) = n_ij \endcode
      <td> If species \c i is a protein or a DNA molecule, then \code n_ij = 0
          \endcode for all \c j.
          If species \c i is a protein-DNA complex, then \c n_ij
          define the stoichiometric vector of the
          unbounding reaction. In other words, it means
          that the protein-DNA complex is formed of:
          \code n_i0*species0 + n_i1*species1 + ... \endcode </tr>
      </table>
   */
  Array<int,2> dnaSpecies_;

  /**
   * Copy method used in the assigment operator and the copy constructor.
   */
  void copy(const Cell& c2);


};

//------------------------------------------------------------------------------

// INLINE ACCESS METHODS

inline int Cell::getNChannels() const
  {return nTotalChannels_;}

inline Array<double,1>& Cell::getPropensities()
  {return totalPropensities_;}

inline Array<double,1>& Cell::getPropensitiesCell()
  {return CellBase::getPropensities();}

inline Array<double,1>& Cell::getPropensitiesCellMilieu()
  {return cellMilieuChemicalSystem_.getPropensities();}

inline bool Cell::isReactionInternal(int mu) const
{
  if (mu <= chemicalSystem_.nChannels_)
  {
    return true;
  } else {
    return false;
  }
}

#endif // CELL_H
