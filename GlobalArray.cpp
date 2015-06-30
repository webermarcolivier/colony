/***************************************************************************//**
 * Project: Colony
 *
 * \file    GlobalArray.cpp
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

#include "GlobalArray.h"

//------------------------------------------------------------------------------

template <class T>
GlobalArray<T>::GlobalArray()
  : size_(0),
    nNodes_(0)
{}

//------------------------------------------------------------------------------

template <class T>
GlobalArray<T>::~GlobalArray()
{
  globalArray_.free();
}

//------------------------------------------------------------------------------

template <class T>
void GlobalArray<T>::build(list< LocalArrayTarget<T> >& localArrayList)
{
  // Set the number of nodes.
  nNodes_ = localArrayList.size();

  int nodeIndex = 0;

  set<int>::iterator itset;
  /*
   * The 'typename' keyword is needed for the compiler to understand that
   *  this statement is a declaration of the variable 'itlist' with type
   *  'list< LocalArrayTarget<T> >::iterator'.
   */
  typename list< LocalArrayTarget<T> >::iterator itlist;
  for (itlist =localArrayList.begin(), itset=nodesTable_.begin();
       itlist!=localArrayList.end();
       itlist++, itset++)
  {
    // Insert the position of the first element of the local array in the
    // nodes table.
    nodesTable_.insert(itset, nodeIndex); // Rem: order of these 2 instructions
    nodeIndex += (*itlist).size;          //      is important
  }
  size_ = nodeIndex; // Total size of the global array.

  // Resize global array (all data are garbage).
  globalArray_.resize(size_);

  int a = 0, b = 0;
  for (itlist =localArrayList.begin(), itset=nodesTable_.begin();
       itlist!=localArrayList.end();
       itlist++)
  {
    // Make a copy of local array.
    Array<T,1> localArrayCopy( (*itlist).ptr->copy() );

    // Set local array to be a reference to a part of the global array.
    a = (*itset);
    itset++;
    b = (*itset) - 1;
    (*itlist).ptr->reference( globalArray_( Range(a, b) ) );

    // Copy back values into local array.
    (*(*itlist).ptr) = localArrayCopy;
  }
}

//------------------------------------------------------------------------------

template <class T>
list<Array<T,1> > GlobalArray<T>::buildGlobalArrayRangesList
  (list< LocalArrayTarget<T> >& localArrayList)
{
  // Set the number of nodes.
  nNodes_ = localArrayList.size();

  int nodeIndex = 0;

  set<int>::iterator itset;
  /*
   * The 'typename' keyword is needed for the compiler to understand that the
   *  this statement is a declaration of the variable 'itlist' with type
   *  'list< LocalArrayTarget<T> >::iterator'.
   */
  typename list< LocalArrayTarget<T> >::iterator itlist;
  for (itlist =localArrayList.begin(), itset=nodesTable_.begin();
       itlist!=localArrayList.end();
       itlist++)
  {
    /* Insert the position of the first element of the local array in the
       nodes table. */

    // Insert the nodeIndex value in the nodes table.
    itset = nodesTable_.insert(itset, nodeIndex);

    // Update the node index with the size of the local array.
    nodeIndex += (*itlist).size;
  }
  size_ = nodeIndex; // Total size of the global array.

  // Resize global array (all data are garbage).
  globalArray_.resize(size_);

  // Build the list of global array parts.

  list<Array<T,1> > globalArrayPartsList;

  int a = 0, b = 0;
  for (itlist =localArrayList.begin(), itset=nodesTable_.begin();
       itlist!=localArrayList.end();
       itlist++)
  {
    // Build the list of the parts of the global array, to be sent to the local arrays to be a referenced.
    a = (*itset);
    itset++;
    if (itset != nodesTable_.end())
    {
      b = (*itset) - 1;
    } else{
      b = size_-1;
    }

    globalArrayPartsList.push_back( globalArray_( Range(a, b) ) );
  }

  return globalArrayPartsList;
}

//------------------------------------------------------------------------------

template <class T>
Array<T,1> GlobalArray<T>::getNode(int k)
{
  if (k > nNodes_ - 1)
  {
    cout << "ERROR: class GlobalArray, getNode(), index k is out of bounds"
         << " k = " << k << " nNodes_ - 1 = " << nNodes_ - 1 << endl;
    exit(1);
  } else
  {
    set<int>::iterator it;
    it = nodesTable_.begin();

    // Jump k elements in the list.
    int j;
    for (j=0; j<k; j++) {it++;}

    int a = *it;
    it++;
    int b = *it - 1;

//    Array<T,1> nodeArray( globalArray_( Range(a, b) ) );
    return globalArray_( Range(a, b) );
  }
}

//------------------------------------------------------------------------------

template <class T>
pair<int,int> GlobalArray<T>::getLocalIndicesFromGlobal(int i) const
{
  // Set an iterator pointing to the node to which belongs element i.
  set<int>::iterator it;
  it = nodesTable_.upper_bound(i);
  it--;

  // Index inside the node
  int j = i - *it;

  // Count the number of nodes from the beginning of the list (reverse order)
  int k;
  for (k = 0; it != nodesTable_.begin(); it--, k++) {} // Node number starts at 0.

  pair<int,int> kj(k, j);
  return kj;
}

//------------------------------------------------------------------------------

/**
 * To avoid linker errors, we have to define the template class we will use.
 * For more information, see the C++ FAQ:\n
 * http://www.parashift.com/c++-faq-lite/templates.html#faq-35.15
 */
template class GlobalArray<double>;

template class GlobalArray<int>;

//------------------------------------------------------------------------------











