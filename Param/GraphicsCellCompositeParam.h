/***************************************************************************//**
 * Project: Colony
 *
 * \file    GraphicsCellCompositeParam.h
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

#ifndef GRAPHICSCELLCOMPOSITEPARAM_H
#define GRAPHICSCELLCOMPOSITEPARAM_H

#include "../debug.h"

// standard C++ header files

// libraries header files
#include <ode/ode.h>

// user header files
#ifdef GUI
  #include "../GraphicsCellScene.h"
  #include "../GraphicsCellGroup.h"
#endif
class SpatialIntegratorODE;

// namespaces



class GraphicsCellCompositeParam
{
public:

  void copy(GraphicsCellCompositeParam& p)
  {
    #ifdef GUI
      cellScene_ = p.cellScene_;
      cellGroup_ = p.cellGroup_;
    #endif

    world_ = p.world_;
    space_ = p.space_;
  }

  #ifdef GUI
    void init(QGraphicsScene* cellScene, GraphicsCellGroup* cellGroup)
    {
      cellScene_ = cellScene;
      cellGroup_ = cellGroup;
    }
  #endif

  #ifdef GUI
    QGraphicsScene* cellScene_; ///< Pointer to the Qt colony cell scene
    GraphicsCellGroup* cellGroup_; ///< Pointer to the Qt graphic object group
  #endif

  dWorldID world_;
  dSpaceID space_;
  SpatialIntegratorODE* spatialIntegratorODE_;
  
};

#endif // GRAPHICSCELLCOMPOSITEPARAM_H
