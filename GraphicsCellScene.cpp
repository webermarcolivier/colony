/***************************************************************************//**
 * Project: Colony
 *
 * \file    GraphicsCellScene.cpp
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


#include "GraphicsCellScene.h"

#ifdef GUI

//------------------------------------------------------------------------------

GraphicsCellScene::GraphicsCellScene(QObject *parent)
    : QGraphicsScene(parent)
{
  wallsList_.clear();
  cellGroup_ = new GraphicsCellGroup();
  addItem(cellGroup_);
}

//------------------------------------------------------------------------------

GraphicsCellScene::~GraphicsCellScene()
{
  removeItem( cellGroup_ );
  delete cellGroup_;
}

//------------------------------------------------------------------------------

void GraphicsCellScene::clearCellGroup()
{
  removeItem( cellGroup_ );
  delete cellGroup_;
  cellGroup_ = new GraphicsCellGroup();
  addItem(cellGroup_);
}

//------------------------------------------------------------------------------

void GraphicsCellScene::init(MainWindow *mainWindow)
{
  mainWindow_ = mainWindow;
  clearCellGroup();
}

//------------------------------------------------------------------------------

GraphicsCellGroup* GraphicsCellScene::getCellGroup()
{
  return cellGroup_;
}

//------------------------------------------------------------------------------

// GraphicsWall* GraphicsCellScene::addWall(float centerX, float centerY,
//                           float normalX, float normalY,
//                           float thickness, float length)
//{
//  GraphicsWall *newWall = new GraphicsWall (centerX, centerY, normalX, normalY,
//                                            thickness, length, 0);
//  wallsList_.append( newWall );
//  addItem(newWall);
//  return newWall;
//}

//------------------------------------------------------------------------------

void GraphicsCellScene::drawPoint(float x, float y)
{
  QGraphicsEllipseItem* point = new QGraphicsEllipseItem ( x-1.0, -y-1.0, 2.0, 2.0);
  point->setBrush(Qt::green);
  point->setZValue(10.0);
  addItem(point);
}

//------------------------------------------------------------------------------

void GraphicsCellScene::on_cellGainedFocus(GraphicsCellQt *cell)
{
  //cell->toggleHasBeenFocused();
  emit cellGainedFocus(cell);
}

//------------------------------------------------------------------------------

void GraphicsCellScene::on_cellLostFocus(GraphicsCellQt *cell)
{
  //cell->toggleHasBeenFocused();
  emit cellLostFocus(cell);
}

//------------------------------------------------------------------------------

#endif // GUI
