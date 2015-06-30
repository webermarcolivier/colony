/***************************************************************************//**
 * Project: Colony
 *
 * \file    GraphicsCellQt.cpp
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

#include "GraphicsCellQt.h"

#ifdef GUI

//------------------------------------------------------------------------------

GraphicsCellQt::~GraphicsCellQt() {}

//------------------------------------------------------------------------------

GraphicsCellQt::GraphicsCellQt()
{}

//------------------------------------------------------------------------------

void GraphicsCellQt::initialize(QGraphicsItem *parent, QGraphicsScene *colonyScene)
{
  QGraphicsItem::setParentItem(parent);

  color_ = Qt::red;
  penWidth_ = GraphicsCellBase::getCellHeight() / 15.0;

  contourColorNotFocused_.setHsv(0,0,65); // light gray
  contourColorFocused_ = Qt::red;
  contourColor_ = contourColorNotFocused_;

  computeCellShape();

  setFlag(QGraphicsItem::ItemIsSelectable, true);
  setFlag(QGraphicsItem::ItemIsFocusable, true);
  setAcceptedMouseButtons(Qt::LeftButton);

  connect (this, SIGNAL(gainedFocus(GraphicsCellQt*)), colonyScene, SIGNAL(cellGainedFocus(GraphicsCellQt*)));
  connect (this, SIGNAL(lostFocus(GraphicsCellQt*)), colonyScene, SIGNAL(cellLostFocus(GraphicsCellQt*)));
}

//------------------------------------------------------------------------------

GraphicsCellQt::GraphicsCellQt(const GraphicsCellQt &cell)
    : QGraphicsItem(cell.parentItem()),
      color_(cell.color_),
      contourColor_(cell.contourColor_),
      contourColorNotFocused_(cell.contourColorNotFocused_),
      contourColorFocused_(cell.contourColorFocused_),
      penWidth_(cell.penWidth_),
      cellShape_(cell.cellShape_)
{
  copy(cell);
}

//------------------------------------------------------------------------------

void GraphicsCellQt::copy(const GraphicsCellQt& cell)
{
  GraphicsCellQt::initialize(cell.parentItem(), cell.scene());
  setTransform(cell.transform());
  QGraphicsItem::setPos(cell.pos());
  color_ = cell.color_;
  contourColor_ = cell.contourColor_;
  contourColorNotFocused_ = cell.contourColorNotFocused_;
  contourColorFocused_ = cell.contourColorFocused_;
  penWidth_ = cell.penWidth_;
  cellShape_ = cell.cellShape_;
}

//------------------------------------------------------------------------------

void GraphicsCellQt::computeCellShape()
{
  prepareGeometryChange();
  cellShape_ = QPainterPath();
  cellShape_.moveTo( -(cellLength_ - cellHeight_)/2.0,   (cellHeight_ - penWidth_)/2.0);
  cellShape_.lineTo(  (cellLength_ - cellHeight_)/2.0,   (cellHeight_ - penWidth_)/2.0);
  cellShape_.cubicTo( (cellLength_ - cellHeight_)/2.0 + (8./6.)*(cellHeight_ - penWidth_)/2.0,
                       (cellHeight_ - penWidth_)/2.0, // first bezier curve point
                       (cellLength_ - cellHeight_)/2.0 + (8./6.)*(cellHeight_ - penWidth_)/2.0,
                      -(cellHeight_ - penWidth_)/2.0, // second bezier curve point
                       (cellLength_ - cellHeight_)/2.0,
                      -(cellHeight_ - penWidth_)/2.0); // final point
  cellShape_.lineTo( -(cellLength_ - cellHeight_)/2.0, - (cellHeight_ - penWidth_)/2.0);
  cellShape_.cubicTo(-(cellLength_ - cellHeight_)/2.0 - (8./6.)*(cellHeight_ - penWidth_)/2.0,
                      -(cellHeight_ - penWidth_)/2.0, // first bezier curve point
                      -(cellLength_ - cellHeight_)/2.0 - (8./6.)*(cellHeight_ - penWidth_)/2.0,
                       (cellHeight_ - penWidth_)/2.0, // second bezier curve point
                      -(cellLength_ - cellHeight_)/2.0,   (cellHeight_ - penWidth_)/2.0); // final point
}

//------------------------------------------------------------------------------

QRectF GraphicsCellQt::boundingRect() const
{
  return cellShape_.boundingRect().adjusted(penWidth_/2.0,penWidth_/2.0,penWidth_/2.0,penWidth_/2.0);
}

//------------------------------------------------------------------------------

void GraphicsCellQt::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget)
{
  computeCellShape();
  // Body
  QPen bodyPen;
  bodyPen.setWidthF(penWidth_);
  bodyPen.setColor(contourColor_);
  painter->setPen(bodyPen);
  painter->setBrush(color_);
  painter->drawPath(cellShape_);
}

//------------------------------------------------------------------------------

QPainterPath GraphicsCellQt::shape() const
{
  return cellShape_;
}

//------------------------------------------------------------------------------

TinyVector<double,2> GraphicsCellQt::fromSceneToNormal(double x, double y) const
{
  TinyVector<double,2> normalPos;
  normalPos(0) = x;
  normalPos(1) = - y;
  return normalPos;
}

//------------------------------------------------------------------------------

TinyVector<double,2> GraphicsCellQt::fromNormalToScene(double x, double y) const
{
  TinyVector<double,2> scenePos;
  scenePos(0) = x;
  scenePos(1) = - y;
  return scenePos;
}

//------------------------------------------------------------------------------

void GraphicsCellQt::setPos(float x, float y, float z)
{
  TinyVector<float,2> scenePos;
  scenePos = fromNormalToScene(x, y);
  QGraphicsItem::setPos( scenePos(0), scenePos(1) );
}


//------------------------------------------------------------------------------

void GraphicsCellQt::setAngle(float angle)
{
  // The y axis is inverted, so the angle is the opposite.
  float angleQt = - angle;
  QTransform transformMatrix;
  transformMatrix.rotate(angleQt);
  setTransform(transformMatrix, false);
}

//------------------------------------------------------------------------------

void GraphicsCellQt::setColor(QColor color)
{
  color_ = color;
}

//------------------------------------------------------------------------------

void GraphicsCellQt::setCellLength(float length)
{
  computeCellShape();
}

//------------------------------------------------------------------------------

void GraphicsCellQt::setCellHeight(float height)
{
  computeCellShape();
}

//------------------------------------------------------------------------------

void GraphicsCellQt::focusInEvent(QFocusEvent *event)
{
  GraphicsCellBase::hasFocus_ = true;
  contourColor_ = contourColorFocused_;
  emit gainedFocus(this);
  QGraphicsItem::focusInEvent(event);
}

//------------------------------------------------------------------------------

void GraphicsCellQt::focusOutEvent(QFocusEvent *event)
{
  // We allow the cell to loose focus only when the reason is a mouse click.
  // In this way, we keep the focus on the cell when selecting antoher window, etc.
  if ( event->reason() == Qt::MouseFocusReason )
  {
    GraphicsCellBase::hasFocus_ = false;
    contourColor_ = contourColorNotFocused_;
    emit lostFocus(this);
    QGraphicsItem::focusOutEvent(event);
  }
}

//------------------------------------------------------------------------------

#endif // GUI
