/***************************************************************************//**
 * Project: Colony
 *
 * \file    GraphicsCellSceneODE.cpp
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

#include "GraphicsCellSceneODE.h"

//------------------------------------------------------------------------------

GraphicsCellSceneODE::GraphicsCellSceneODE(QObject *parent)
    : GraphicsCellScene(parent), ODESimulation()
{
}

//------------------------------------------------------------------------------

GraphicsWall* GraphicsCellSceneODE::addWall(float centerX, float centerY,
                          float normalX, float normalY,
                          float thickness, float length)
{
 GraphicsWall *newWall = new GraphicsWall (centerX, centerY, normalX, normalY,
                                           thickness, length, 0);
 wallsList_.append( newWall );
 addItem(newWall);

  // Add ODE object
  dGeomID newWallODE = dCreateBox (space_, length, thickness, 100.0);
  // - Rotate object.
  dMatrix3 R;
  double alpha = atan( normalY / normalX ) * 180.0 / PI;
  dRFromAxisAndAngle(R, 0.0, 0.0, 1.0, PI*(90.0-alpha)/180.0);
  dGeomSetRotation(newWallODE, R);
  dGeomSetPosition (newWallODE, centerX - cos(alpha)*thickness/2.0, centerY - sin(alpha)*thickness/2.0, 0.0);
  wallsODEList_.append(newWallODE);

  return newWall;
}

//------------------------------------------------------------------------------

void GraphicsCellSceneODE::setCellPosition(int cellIndex, double x, double y)
{
  TinyVector<double,2> normalPos;
  normalPos = fromSimToNormal(x, y);
  qreal xNormal, yNormal;
  xNormal = normalPos(0);
  yNormal = normalPos(1);

  static_cast<GraphicsCellODE*>((*this)[cellIndex])->GraphicsCellODE::setPos(xNormal, yNormal);
}

//------------------------------------------------------------------------------

GraphicsCell* GraphicsCellSceneODE::addOneCell(int cellIndex)
{
  GraphicsCell *newCell = new GraphicsCellODE( cellIndex, cellGroup_, this, world_, space_);
  cellsList_.append( newCell );
  ++nCells_;
  return newCell;
}

//------------------------------------------------------------------------------

void GraphicsCellSceneODE::duplicateCell(int motherCellIndex)
{
  //cout << "GraphicsCellSceneODE::duplicateCell start" << endl;
  int nCellsBeforeDuplication = nCells_;

  // Get the ODE world and space of the mother cell.
  dWorldID world = static_cast<GraphicsCellODE*>(cellsList_[motherCellIndex])->ODE_getWorld();
  dSpaceID space = static_cast<GraphicsCellODE*>(cellsList_[motherCellIndex])->ODE_getSpace();
  // Allocate new ODE cell.
  GraphicsCell *newCell = new GraphicsCellODE(nCellsBeforeDuplication, cellGroup_, this, world, space);
  // Duplicate mother cell and write the daughter cell in the newCell pointer.
  static_cast<GraphicsCellODE*>(cellsList_[motherCellIndex])->GraphicsCellODE::duplicateCell(newCell, cellLength0_);
  cellsList_.append( newCell );
  ++nCells_;
  //cout << "GraphicsCellSceneODE::duplicateCell end. nCells = " << nCells_ << endl;
}

//------------------------------------------------------------------------------

void GraphicsCellSceneODE::deleteOneCell(int cellIndex)
{
  int i;
  for (i=0; i<cellsList_.size(); ++i)
  {
    if (cellsList_[i]->getCellIndex() == cellIndex )
    {
      cellsList_[i]->clearFocus();
      removeItem(cellsList_[i]);
      static_cast<GraphicsCellODE*>(cellsList_[i])->~GraphicsCellODE();
      cellsList_.removeAt(i);
    }
  }
  --nCells_;

  for (i=0; i<cellsList_.size(); ++i)
  {
    cellsList_[i]->setCellIndex(i);
  }
}


//------------------------------------------------------------------------------

void GraphicsCellSceneODE::computeODEStep(float timeStep, bool equilibration)
{
  /// - Compute the ODE step.
  ODESimulation::computeStep(timeStep, &cellsList_, colonySize_, equilibration,
                             mainWindow_->getCellLength0()*1.442695041,
                             mainWindow_->getCellHeight0());

  /// - Send the cell positions and angle to the simulation/CellCollection object.
  TinyVector<double,2> position;
  double angle;
  int cellIndex;
  for (int i=0; i<cellsList_.size(); ++i)
  {
    // Get the cell index.
    cellIndex = cellsList_[i]->getCellIndex();

    // Get the position of the cell from the ODE object.
    GraphicsCellODE* cell = static_cast<GraphicsCellODE*>(cellsList_[i]);
    position = cell->ODE_getPosition();
    // Set the position of the cell to the new position (both in GraphicsCell and ODE objects).
    qreal x = position(0);
    qreal y = position(1);
    cell->GraphicsCell::setPos(x, y);

    // Get the angle of the cell from the ODE object.
    angle = cell->ODE_getAngle();
    // Set the angle of the cell to the new angle (both in GraphicsCell and ODE objects).
    cell->GraphicsCell::setAngle(angle);

    // Adapt position from normal to simulation coordinate systems.
    position = fromNormalToSim(x, y);

    // Send the new position and angle to the simulator.
    mainWindow_->setSimulatorCellPosition(cellIndex, position);
    mainWindow_->setSimulatorCellAngle(cellIndex, angle);
  }
  mainWindow_->updateSimulatorPositionAndAngleTimeSlice();
}

//------------------------------------------------------------------------------
