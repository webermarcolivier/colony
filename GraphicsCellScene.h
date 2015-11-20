/***************************************************************************//**
 * Project: Colony
 *
 * \file    GraphicsCellScene.h
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

#ifndef GRAPHICSCELLSCENE_H
#define GRAPHICSCELLSCENE_H

#include "compilation_options.h"

#ifdef GUI

// standard C++ header files
#include <list>

// libraries header files
#include <blitz/array.h>
#include <QGraphicsItem>
#include <QGraphicsScene>
#include <QGraphicsRectItem>
#include <QGraphicsRotation>
#include <QList>
#include <math.h>

// user header files
#include "GraphicsWall.h"

// namespaces
using std::cout;
using std::endl;
using std::list;
using blitz::Array;
using blitz::TinyVector;

class GraphicsCellGroup;
class GraphicsCellQt;
class MainWindow;



/**
  * Graphics scene of the cell colony in the Qt Graphics framework.
  * Contains the list of pointers to the GraphicsCell objects,
  * defines method to add and remove cells,
  * defines access to a cell through the [] operator.
  */
class GraphicsCellScene : public QGraphicsScene
{
  Q_OBJECT

public:

  GraphicsCellScene(QObject *parent = 0);
  ~GraphicsCellScene();

  void init(MainWindow *mainWindow);

  GraphicsCellGroup* getCellGroup();

  void clearCellGroup();

//  GraphicsWall* addWall(float centerX, float centerY,
//                        float normalX, float normalY,
//                        float thickness, float length);

  void drawPoint(float x, float y);

signals:

  void cellGainedFocus(GraphicsCellQt *cell);
  void cellLostFocus(GraphicsCellQt *cell);

public slots:

  void on_cellGainedFocus(GraphicsCellQt *cell);
  void on_cellLostFocus(GraphicsCellQt *cell);

protected:

  /**
    * GraphicsItem group of all the cells.
    */
  GraphicsCellGroup *cellGroup_;
  /**
    * List of the walls graphics items in the scene.
    */
  QList<GraphicsWall*> wallsList_;

  MainWindow *mainWindow_;


};

#include "GraphicsCellGroup.h"
#include "GraphicsCellQt.h"

#endif // GUI

#endif // GRAPHICSCELLSCENE_H
