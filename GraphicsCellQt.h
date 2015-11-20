/***************************************************************************//**
 * Project: Colony
 *
 * \file    GraphicsCellQt.h
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

#ifndef GRAPHICSCELLQT_H
#define GRAPHICSCELLQT_H

#include "compilation_options.h"

#ifdef GUI

#define PI 3.14159265

// standard C++ header files
#include <math.h>

// libraries header files
#include <blitz/array.h>
#include <QGraphicsItem>
#include <QObject>
#include <QEvent>
#include <QPainter>
#include <QGraphicsScene>
#include <QStyleOption>
#include <QPen>
#include <QFocusEvent>

// user header files
#include "GraphicsCellBase.h"

// namespaces
using blitz::Array;
using blitz::TinyVector;
class QFocusEvent;
class QGraphicsItem;
class QGraphicsScene;
class GraphicsCellScene;


/**
  * %Cell in the Qt Graphics framework.
  * Draw a cell item in the graphics view. Can change color, orientation (angle),
  * length and can be selected (by emitting the focus gained/lost events).
  */
class GraphicsCellQt : public QObject, public QGraphicsItem, virtual public GraphicsCellBase
{
  Q_OBJECT

public:

  GraphicsCellQt();
  GraphicsCellQt(const GraphicsCellQt &cell);
  ~GraphicsCellQt();

  void initialize(QGraphicsItem *parent, QGraphicsScene *colonyScene);
  void copy(const GraphicsCellQt& cell);

  QRectF boundingRect() const;
  QPainterPath shape() const;
  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
  void setColor(QColor color);

  /**
    * Converts Qt Graphics Scene coordinates to normal coordinates (y axis inverted).
    */
  TinyVector<double,2> fromSceneToNormal(double x, double y) const;

  /**
    * Converts normal coordinates to Qt Graphics Scene coordinates (y axis inverted).
    */
  TinyVector<double,2> fromNormalToScene(double x, double y) const;

signals:

  void gainedFocus(GraphicsCellQt *cell);
  void lostFocus(GraphicsCellQt *cell);

protected:

  void focusInEvent(QFocusEvent *event);
  void focusOutEvent(QFocusEvent *event);

  /**
    * Overloaded version: set the position of GraphicsCellQt object,
    * warning: coordinates in argument should be in normal system.
    */
  void setPos(float x, float y, float z);
  /**
    * Overloaded version: set the angle of GraphicsCellQt object.
    */
  void setAngle(float angle);
  /**
    * Overloaded version: set the length of GraphicsCellQt object.
    */
  void setCellLength(float length);
  /**
    * Overloaded version: set the height of GraphicsCellQt object.
    */
  void setCellHeight(float height);


private:

  /**
    * Compute the cell shape from the parameters cellLength_ and cellHeight_.
    * The cell shape is a cylinder with spherical caps.
    */
  void computeCellShape();

  /**
    * Angle in the Qt environment.
    */
  QColor color_;
  QColor contourColor_;
  QColor contourColorNotFocused_;
  QColor contourColorFocused_;
  QPainterPath cellShape_;
  qreal penWidth_;

};

#include "GraphicsCellScene.h"

#endif // GUI

#endif // GRAPHICSCELLQT_H
