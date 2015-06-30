/***************************************************************************//**
 * Project: Colony
 *
 * \file    GraphicsCellBase.cpp
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

#include "GraphicsCellBase.h"

//------------------------------------------------------------------------------

GraphicsCellBase::~GraphicsCellBase() {}

//------------------------------------------------------------------------------

GraphicsCellBase::GraphicsCellBase()
  : cellLength_(0.0),
    cellHeight_(0.0),
    angle_(0.0),
    cellIndex_(0),
    posX_(0.0),
    posY_(0.0),
    posZ_(0.0),
    hasFocus_(false)
{}

//------------------------------------------------------------------------------

GraphicsCellBase::GraphicsCellBase(const GraphicsCellBase &cell)
{
  copy(cell);
}

//------------------------------------------------------------------------------

void GraphicsCellBase::copy(const GraphicsCellBase& c)
{
  cellIndex_ = c.cellIndex_;
  cellLength_ = c.cellLength_;
  cellHeight_ = c.cellHeight_;
  cellVolume_ = c.cellVolume_;
  angle_ = c.angle_;
  posX_ = c.posX_;
  posY_ = c.posY_;
  posZ_ = c.posZ_;
}

//------------------------------------------------------------------------------

void GraphicsCellBase::setCellLength(float length)
{
  cellLength_ = length;
}

//------------------------------------------------------------------------------

void GraphicsCellBase::setCellHeight(float height)
{
  cellHeight_ = height;
}

//------------------------------------------------------------------------------

void GraphicsCellBase::setCellIndex(int cellIndex)
{
  cellIndex_ = cellIndex;
}

//------------------------------------------------------------------------------

void GraphicsCellBase::setPos(float x, float y, float z)
{
  posX_ = x;
  posY_ = y;
  posZ_ = z;
}

//------------------------------------------------------------------------------

void GraphicsCellBase::setAngle(float angle)
{
  angle_ = angle;
}

//------------------------------------------------------------------------------

int GraphicsCellBase::getCellIndex() const
{
  return cellIndex_;
}

//------------------------------------------------------------------------------

float GraphicsCellBase::getCellLength() const
{
  return cellLength_;
}

//------------------------------------------------------------------------------

float GraphicsCellBase::getCellHeight() const
{
  return cellHeight_;
}

//------------------------------------------------------------------------------

float GraphicsCellBase::getAngle() const
{
  return angle_;
}

//------------------------------------------------------------------------------
