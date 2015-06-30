/***************************************************************************//**
 * Project: Colony
 *
 * \file    StringTableModel.cpp
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

#include "StringTableModel.h"

#ifdef GUI

//------------------------------------------------------------------------------

StringTableModel::StringTableModel(const Array<QString,2> &stringArray, QObject *parent)
  : QAbstractTableModel(parent), stringArray_( stringArray.copy() )
{
  headerHorizontal_.resize( stringArray_.extent(secondDim) );
  headerVertical_.resize( stringArray_.extent(firstDim) );

  for (int i=0; i<headerHorizontal_.size(); ++i)
  {
    headerHorizontal_(i).setNum(i);
  }
  emit headerDataChanged(Qt::Horizontal, 0, headerHorizontal_.size()-1);
  for (int i=0; i<headerVertical_.size(); ++i)
  {
    headerVertical_(i).setNum(i);
  }
  emit headerDataChanged(Qt::Vertical,0,headerVertical_.size()-1);
}

//------------------------------------------------------------------------------

int StringTableModel::rowCount(const QModelIndex &parent) const
{
  if (!parent.isValid())
  {
    return stringArray_.extent(firstDim);
  }
  else
  {
    return 0;
  }
}

//------------------------------------------------------------------------------


int StringTableModel::columnCount(const QModelIndex &parent) const
{
  if (!parent.isValid())
  {
    return stringArray_.extent(secondDim);
  }
  else
  {
    return 0;
  }
}

//------------------------------------------------------------------------------

QVariant StringTableModel::data(const QModelIndex &index, int role) const
{
  if (!index.isValid())
    return QVariant();

  if (index.row() >= stringArray_.extent(firstDim))
    return QVariant();

  if (index.column() >= stringArray_.extent(secondDim))
    return QVariant();

  if (role == Qt::DisplayRole)
    return stringArray_(index.row(),index.column());
  else
    return QVariant();
}

//------------------------------------------------------------------------------

Qt::ItemFlags StringTableModel::flags(const QModelIndex &index) const
{
  return Qt::ItemIsSelectable;
}

//------------------------------------------------------------------------------

bool StringTableModel::setData(const QModelIndex & index, const QVariant & value, int role)
{
  if (index.isValid() && role == Qt::EditRole)
  {
    stringArray_(index.row(),index.column()) = value.toString();
    emit dataChanged(index,index);
    return true;
  }
  return false;
}

//------------------------------------------------------------------------------

bool StringTableModel::setHeaderData ( int section, Qt::Orientation orientation, const QVariant & value, int role )
{
  if ( orientation == Qt::Horizontal )
  {
    if ( section >= headerHorizontal_.size() )
      return false;

    if (role == Qt::EditRole)
    {
      headerHorizontal_(section) = value.toString();
      emit headerDataChanged(Qt::Horizontal, section, section);
      return true;
    }
    else
      return false;
  }
  else if ( orientation == Qt::Vertical )
  {
    if ( section >= headerVertical_.size() )
      return false;

    if (role == Qt::EditRole)
    {
      headerVertical_(section) = value.toString();
      emit headerDataChanged(Qt::Vertical, section, section);
      return true;
    }
    else
      return false;
  }
  else
  {
    return false;
  }
}

//------------------------------------------------------------------------------

QVariant StringTableModel::headerData ( int section, Qt::Orientation orientation, int role ) const
{
  if ( orientation == Qt::Horizontal )
  {
    if ( section >= stringArray_.extent(secondDim) )
      return QVariant();

    if (role == Qt::DisplayRole)
      return headerHorizontal_(section);
    else
      return QVariant();
  }
  else if ( orientation == Qt::Vertical )
  {
    if ( section >= stringArray_.extent(firstDim) )
      return QVariant();

    if (role == Qt::DisplayRole)
      return headerVertical_(section);
    else
      return QVariant();
  }
  else
  {
    return QVariant();
  }
}

//------------------------------------------------------------------------------

#endif // GUI
