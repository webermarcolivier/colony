/***************************************************************************//**
 * Project: Colony
 *
 * \file    GraphicsCellGroup.h
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

#ifndef GRAPHICSCELLGROUP_H
#define GRAPHICSCELLGROUP_H

#include "debug.h"

#ifdef GUI

// standard C++ header files

// libraries header files
#include <QGraphicsItem>
#include <QGraphicsScene>

// user header files

// namespaces



class GraphicsCellGroup : public QGraphicsItem
{

public:

  GraphicsCellGroup(QGraphicsItem *parent = 0);

  QRectF boundingRect() const;
  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);


private:


};

#endif // GUI

#endif // GRAPHICSCELLGROUP_H
