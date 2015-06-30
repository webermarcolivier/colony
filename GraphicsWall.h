/***************************************************************************//**
 * Project: Colony
 *
 * \file    GraphicsWall.h
 * \author  Marc Weber\n
 *          The SiMBioSys group (CosmoLab)\n
 *          Parc Científic de Barcelona\n
 *          Barcelona, Spain.\n
 *          http://www.thesimbiosys.com
 * \version 1.0
 * \date    02/2011
 *
 *          Copyright 2009 by Marc Weber
 ******************************************************************************/

#ifndef GRAPHICSWALL_H
#define GRAPHICSWALL_H

#ifdef GUI

#define PI 3.14159265

// standard C++ header files
#include <math.h>

// libraries header files
#include <QGraphicsItem>
#include <QGraphicsRectItem>
#include <QGraphicsRotation>
#include <QObject>
#include <QEvent>
#include <QPainter>
#include <QStyleOption>
#include <QPen>
#include <QFocusEvent>

// user header files

// namespaces
class QGraphicsItem;
class QGraphicsScene;


/**
  * Wall in the Qt Graphics framework.
  * Draw a wall item in the graphics view. Can change color, orientation (angle),
  * length.
  */
class GraphicsWall : public QObject, public QGraphicsRectItem
{
  Q_OBJECT

public:

  GraphicsWall(float centerX, float centerY,
               float normalX, float normalY,
               float thickness, float length,
               QGraphicsItem *parent);
  ~GraphicsWall();
};

#endif // GUI

#endif // GRAPHICSWALL_H
