/***************************************************************************//**
 * Project: Colony
 *
 * \file    GraphicsCellComposite.cpp
 * \author  Marc Weber\n
 *          The SiMBioSys group (CosmoLab)\n
 *          Parc Científic de Barcelona\n
 *          Barcelona, Spain.\n
 *          http://www.thesimbiosys.com
 * \version 1.0
 * \date    05/2011
 *
 *          Copyright 2009 by Marc Weber
 ******************************************************************************/

#include "GraphicsCellComposite.h"

//------------------------------------------------------------------------------

GraphicsCellComposite::~GraphicsCellComposite() {}

//------------------------------------------------------------------------------

GraphicsCellComposite::GraphicsCellComposite()
{}

//------------------------------------------------------------------------------

#ifdef GUI
GraphicsCellComposite::GraphicsCellComposite(const GraphicsCellComposite &cell)
  : GraphicsCellBase(cell),
    GraphicsCellQt(cell),
    GraphicsCellODE(cell)
{}
#else
  GraphicsCellComposite::GraphicsCellComposite(const GraphicsCellComposite &cell)
    : GraphicsCellBase(cell),
      GraphicsCellODE(cell)
  {}
#endif

//------------------------------------------------------------------------------

#ifdef GUI
      void GraphicsCellComposite::copy(const GraphicsCellComposite& c)
      {
        GraphicsCellBase::copy(c);
        GraphicsCellQt::copy(c);
        GraphicsCellODE::copy(c);
      }
#else
  void GraphicsCellComposite::copy(const GraphicsCellComposite& c)
  {
    GraphicsCellBase::copy(c);
    GraphicsCellODE::copy(c);
  }
#endif

//------------------------------------------------------------------------------

#ifdef GUI  
  void GraphicsCellComposite::initialize(GraphicsCellCompositeParam& p)
  {
    GraphicsCellQt::initialize(p.cellGroup_, p.cellScene_);
    GraphicsCellODE::initialize(p.world_, p.space_, p.spatialIntegratorODE_);
  }
#else
  void GraphicsCellComposite::initialize(GraphicsCellCompositeParam& p)
  {
    GraphicsCellODE::initialize(p.world_, p.space_, p.spatialIntegratorODE_);
  }
#endif

//------------------------------------------------------------------------------

#ifdef GUI
  void GraphicsCellComposite::setPosition(float x, float y, float z)
  {
    GraphicsCellBase::setPos(x, y, z);
    GraphicsCellQt::setPos(x, y, z);
    GraphicsCellODE::setPos(x, y, z);
  }
#else
  void GraphicsCellComposite::setPosition(float x, float y, float z)
  {
    GraphicsCellBase::setPos(x, y, z);
    GraphicsCellODE::setPos(x, y, z);
  }
#endif

//------------------------------------------------------------------------------

#ifdef GUI
  void GraphicsCellComposite::setAngle(float angle)
  {
    float angleBase;
    if (angle < 0.0)
    {
      // Correct negative angles to the range [0,360].
     angleBase = 360.0 + angle;
    } else if (angle > 360.0) {
      angleBase = angle - 360.0;
    } else {
      angleBase = angle;
    }
    GraphicsCellBase::setAngle(angleBase);
    GraphicsCellQt::setAngle(angleBase);
    GraphicsCellODE::setAngle(angleBase);
  }
#else
  void GraphicsCellComposite::setAngle(float angle)
  {
    float angleBase;
    if (angle < 0.0)
    {
      // Correct negative angles to the range [0,360].
     angleBase = 360.0 + angle;
    } else if (angle > 360.0) {
      angleBase = angle - 360.0;
    } else {
      angleBase = angle;
    }
    GraphicsCellBase::setAngle(angleBase);
    GraphicsCellODE::setAngle(angleBase);
  }
#endif


//------------------------------------------------------------------------------

#ifdef GUI
  void GraphicsCellComposite::setCellLength(float length)
  {
    GraphicsCellBase::setCellLength(length);
    GraphicsCellQt::setCellLength(length);
    GraphicsCellODE::setCellLength(length);
  }
#else
  void GraphicsCellComposite::setCellLength(float length)
  {
    GraphicsCellBase::setCellLength(length);
    GraphicsCellODE::setCellLength(length);
  }
#endif

//------------------------------------------------------------------------------

#ifdef GUI
  void GraphicsCellComposite::setCellHeight(float height)
  {
    GraphicsCellBase::setCellHeight(height);
    GraphicsCellQt::setCellHeight(height);
    GraphicsCellODE::setCellHeight(height);
  }
#else
  void GraphicsCellComposite::setCellHeight(float height)
  {
    GraphicsCellBase::setCellHeight(height);
    GraphicsCellODE::setCellHeight(height);
  }
#endif

//------------------------------------------------------------------------------

#ifdef GUI
  void GraphicsCellComposite::setColor(QColor color)
  {
    GraphicsCellQt::setColor(color);
  }
#endif

//------------------------------------------------------------------------------
