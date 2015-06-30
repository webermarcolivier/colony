/***************************************************************************//**
 * Project: Colony
 *
 * \file    GlobalArray.h
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

#ifndef GLOBALARRAY_H
#define GLOBALARRAY_H

#include "debug.h"

// standard C++ header files
#include <iostream>
#include <set>
#include <list>
#include <utility>

// libraries header files
#include <blitz/array.h>

// user header files

// namespaces
using blitz::Array;
using blitz::Range;
using std::cout;
using std::ostream;
using std::endl;
using std::set;
using std::pair;
using std::list;



/**
 * Template Class that contains a pointer to a local array <T,1> and its size.
 */
template<class T>
class LocalArrayTarget {
public:
  Array<T,1>* ptr;
  int size;
};



/**
  * This class defines a global array of type T which stores consecutively a series
  * of local arrays, called nodes. Each node is itself a smaller array of type T,
  * and can have an arbitrary size.
  */

// Don't look into this class definition if you don't want to get a headache...

template<class T>
class GlobalArray
{
public:

//---- LIFECYCLE

  GlobalArray();

//  GlobalArray(const GlobalArray& g2);

  ~GlobalArray();

  /**
   * Construct/Reconstruct the global array based on the local arrays information.
   * The method goes like this:
   * - Set nNodes.
   * - Loop over local array list:
   *   - Calculate the node table and the total size of the global array.
   * - Resize global array.
   * - Loop over local array list:
   *   - Make a copy of local array.
   *   - Reference local array to corresponding part of the global array.
   *   - Copy back values into local array.
   *
   * @param  [in,out] localArraysList list of the local arrays. Each element of the
   * list is a LocalArrayTarget<T> containing a <B>pointer</B> to the array and its size.
   */
  void build (list< LocalArrayTarget<T> >& localArrayList);

// ... I warned you...

  /**
   * Construct/Reconstruct the global array based on the local arrays information,
   * let the responsability of local array referencing to the interface.
   * The method goes like this:
   * - Set nNodes.
   * - Loop over local array list:
   *   - Calculate the node table and the total size of the global array.
   * - Resize global array.
   * - Return a list of global array parts corresponding to the local arrays.
   *   This list can be used as a reference targets list for referencing the
   *   local arrays.
   *
   * Remark 1: The advantage of letting the responsability of the referencing
   * to the interface class is that the referencing method called in the concerned
   * objects (e.g. Cell objects) can control the referencing. If the local array
   * is depending on other members, if it is itself sharing reference to other arrays,
   * then the method in the object can resolve the local references and dependencies.
   * In the #build() method, we change the reference directly by dereferencing the
   * pointer to the local array. This is dangerous since we are changing directly
   * the local array in an object.
   *
   * Remark 2:
   * This method only needs to know the size of the local arrays, because the
   * referencing method is let to the GlobalArrayInterface class. However, it
   * takes as input a LocalArrayTarget list containing the pointers to the local
   * arrays. This is just for sake of generality: this way, we only need to
   * define one method to get the list of local arrays in the interface class.
   *
   * @param  [in,out] localArraysList list of the local arrays. Each element of the
   * list is a LocalArrayTarget<T> containing a <B>pointer</B> to the array and its size.
   */
  list<Array<T,1> > buildGlobalArrayRangesList
  (
    list< LocalArrayTarget<T> >& localArrayList
  );


//---- OPERATORS

//---- ACCESS

  /**
   * Get the global array (by reference).
   */
  Array<T,1>& getGlobalArray()
    {return globalArray_;}

  /**
   * Returns the array corresponding to node \c k.
   * Remark: the size of the returned array can vary.\n
   * Remark: \code 0 <= k <= nNodes_ - 1 \endcode
   */
   Array<T,1> getNode(int k);

  /**
   * Get the local indices (node number \c k, position inside the node \c j)
   * from the global index \c i.
   */
   pair<int,int> getLocalIndicesFromGlobal(int i) const;


//---- INQUIRY

//---- OPERATIONS


private:

//---- DATA

  /**
   * Total size (number of elements of type T) of the global array.
   */
  int size_;

  /**
   * Number of nodes in the global array.
   */
  int nNodes_;

  /**
    * Stores consecutively all the local arrays (which are referenced to a part
    * of the global array).
    */
  Array<T,1> globalArray_;

  /**
   * A table of the indices of the first element of each node. Nodes can have
   * a different size and this table allows for a direct access to a node in the
   * global array.
   */
  set<int> nodesTable_;


};



#endif // GLOBALARRAY_H
