/***************************************************************************//**
 * Project: Colony
 *
 * \file    GraphicsCellGroup.cpp
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

#include "GraphicsCellGroup.h"

#ifdef GUI

//------------------------------------------------------------------------------

GraphicsCellGroup::GraphicsCellGroup(QGraphicsItem *parent)
    : QGraphicsItem(parent)
{
}

//------------------------------------------------------------------------------

QRectF GraphicsCellGroup::boundingRect() const
{
  return childrenBoundingRect();
}

//------------------------------------------------------------------------------

void GraphicsCellGroup::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget)
{}

//------------------------------------------------------------------------------

#endif // GUI
