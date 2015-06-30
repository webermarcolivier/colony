/***************************************************************************//**
 * Project: Colony
 *
 * \file    GraphicsCellSceneODE.h
 * \author  Marc Weber\n
 *          The Si.M.Bio.Sys. Group (CosmoLab)\n
 *          Parc Científic de Barcelona\n
 *          Barcelona, Spain.\n
 *          http://thesimbiosys.isgreat.org
 * \version 0.1
 * \date    11/2009
 *
 *          Copyright 2009 by Marc Weber
 ******************************************************************************/

#ifndef GRAPHICSCELLSCENEODE_H
#define GRAPHICSCELLSCENEODE_H

#include "debug.h"
#include "GraphicsCellScene.h"
#include "ODESimulation.h"

// namespaces


class MainWindow;


/**
  * Graphics scene of the cell colony in the Qt Graphics framework with Open Dynamics Engine (ODE) simulation.
  * Derives from GraphicsCellScene, ODESimulation.
  */
class GraphicsCellSceneODE : public GraphicsCellScene, public ODESimulation
{

public:

  GraphicsCellSceneODE(QObject *parent = 0);

  /**
    * Set the position of cell with index cellIndex.
    * Input coordinates are in the simulation coordinates system with range [0,1]
    * and conversion to Graphics coordinate is automatically done.
    */
  void setCellPosition(int cellIndex, double x, double y);

  /**
    * Add a new cell to the colony scene.
    */
  virtual GraphicsCell *addOneCell(int cellIndex);

  /**
    * Duplicate cell.
    */
  virtual void duplicateCell(int motherCellIndex);

  virtual void deleteOneCell(int cellIndex);

  /**
    * Compute timestep of the ODE simulation.
    */
  void computeODEStep(float timeStep, bool equilibration);

  GraphicsWall* addWall(float centerX, float centerY,
                            float normalX, float normalY,
                            float thickness, float length);

private:

  QList<dGeomID> wallsODEList_;

};

#include "mainwindow.h"

#endif // GRAPHICSCELLSCENEODE_H
