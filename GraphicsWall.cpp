/***************************************************************************//**
 * Project: Colony
 *
 * \file    GraphicsWall.cpp
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

#include "GraphicsWall.h"

#ifdef GUI

//------------------------------------------------------------------------------

GraphicsWall::~GraphicsWall()
{}

//------------------------------------------------------------------------------

GraphicsWall::GraphicsWall(float centerX, float centerY,
                           float normalX, float normalY,
                           float thickness, float length,
                           QGraphicsItem *parent)
    : QGraphicsRectItem(parent)
{
  setRect(0., 0., length, thickness);
  setTransformOriginPoint(length/2.0, 0.0);
  double alpha = atan( -normalY / normalX ) * 180.0 / PI;
  setRotation(180.0-(90.0-alpha)); // the 180.0-... is for the inverted y-axis of Qt,
  // the 90.0-... is because the rectangle is first created horizontal with normal (0,1).
  setPos(centerX-length/2.0, -centerY);
  QPen linePen;
  linePen.setWidthF(1.0);
  linePen.setColor(Qt::black);
  setPen(linePen);
  QBrush fillBrush(Qt::gray);
  setBrush(fillBrush);
}

//------------------------------------------------------------------------------

#endif // GUI
