/***************************************************************************//**
 * Project: Colony
 *
 * \file    StringTableModel.h
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

#ifndef STRINGTABLEMODEL_H
#define STRINGTABLEMODEL_H

#include "debug.h"

#ifdef GUI

// standard C++ header files
#include <iostream>

// libraries header files
#include <blitz/array.h>
#include <QString>
#include <QObject>
#include <QModelIndex>
#include <QVariant>
#include <QAbstractTableModel>

// user header files
#include "Simulator.h"

// namespaces
using blitz::Array;
using std::cout;
using std::endl;


/**
  * String table model in the Qt Model/View framework.
  */
class StringTableModel : public QAbstractTableModel
{
  Q_OBJECT

public:

  StringTableModel(const Array<QString,2> &stringArray, QObject *parent=0);

  int rowCount(const QModelIndex &parent = QModelIndex()) const;
  int columnCount(const QModelIndex &parent = QModelIndex()) const;

  QVariant data(const QModelIndex &index, int role) const;

  Qt::ItemFlags flags(const QModelIndex &index) const;
  bool setData(const QModelIndex & index, const QVariant & value, int role = Qt::EditRole);
  QVariant headerData ( int section, Qt::Orientation orientation, int role = Qt::DisplayRole ) const;
  bool setHeaderData ( int section, Qt::Orientation orientation, const QVariant & value, int role = Qt::EditRole );

private:

  /**
    * Array that contain the strings of the table.
    */
  Array<QString,2> stringArray_;
  /**
    * Array that contains the strings of the horizontal header.
    */
  Array<QString,1> headerHorizontal_;
  /**
    * Array that contains the strings of the vertical header.
    */
  Array<QString,1> headerVertical_;
};

#endif // GUI

#endif // STRINGTABLEMODEL_H
