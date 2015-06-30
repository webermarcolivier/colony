/***************************************************************************//**
 * Project: Colony
 *
 * \file    GraphicsCellGL.cpp
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

#include "GraphicsCellGL.h"

#ifdef GUI

#include "OpenGLViewer.h"

//------------------------------------------------------------------------------

GraphicsCellGL::GraphicsCellGL()
{}

//------------------------------------------------------------------------------

GraphicsCellGL::~GraphicsCellGL()
{
  /// Delete reference in the list of cell in OpenGLViewer class.
  if (isAddedInOpenGLViewerCellsList_)
  {
    openGLViewer_->getCellsList()->removeAll(this);
    isAddedInOpenGLViewerCellsList_ = false;
  }

  delete cylinder_;
  delete arrowCylinder_;
}

//------------------------------------------------------------------------------

void GraphicsCellGL::copy(const GraphicsCellGL& cell)
{
  initialize(cell.openGLViewer_);
}

//------------------------------------------------------------------------------

GraphicsCellGL::GraphicsCellGL(const GraphicsCellGL &cell)
{
  copy(cell);
}

//------------------------------------------------------------------------------

void GraphicsCellGL::initialize(OpenGLViewer* openGLViewer)
{
  openGLViewer_ = openGLViewer;

  if (!isAddedInOpenGLViewerCellsList_)
  {
    openGLViewer_->getCellsList()->append(this);
    isAddedInOpenGLViewerCellsList_ = true;
  }

  /// - Initiate the openGL cylinder (quadric pointer).
  cylinder_ = gluNewQuadric();
  arrowCylinder_ = gluNewQuadric();
  gluQuadricNormals(cylinder_, GLU_SMOOTH);
  gluQuadricTexture(cylinder_, GL_TRUE);
  gluQuadricNormals(arrowCylinder_, GLU_SMOOTH);
  gluQuadricTexture(arrowCylinder_, GL_TRUE);

  /// - Set a random color.
  color_[0] = rand()%255;
  color_[1] = rand()%255;
  color_[2] = rand()%255;
}

//------------------------------------------------------------------------------

void GraphicsCellGL::setPos(float x, float y, float z)
{
  render();
}

//------------------------------------------------------------------------------

void GraphicsCellGL::setAngle(float angle)
{
  render();
}

//------------------------------------------------------------------------------

void GraphicsCellGL::setCellLength(float length)
{
  render();
}

//------------------------------------------------------------------------------

void GraphicsCellGL::setCellHeight(float height)
{
  render();
}

//------------------------------------------------------------------------------

void GraphicsCellGL::setColor(int r, int g, int b)
{
  color_[0] = r;
  color_[1] = g;
  color_[2] = b;
}

//------------------------------------------------------------------------------

void GraphicsCellGL::render() const
{
  // Render the cell
  glPushMatrix();
  glTranslatef(posX_, posY_, posZ_);
  glRotatef(-90.0, 1.0, 0.0, 0.0);
  glRotatef(90.0 - angle_, 0.0, 1.0, 0.0);
  renderInLocalFrame();
  glPopMatrix();

  // If cell is focused, draw an arrow above it.
  if ( hasFocus_ )
  {
    glPushMatrix();
    glTranslatef(posX_, posY_, posZ_);
    float radius = GraphicsCellBase::getCellHeight();
    float arrowTotalLength = 5.0*radius;
    float arrowRadius = arrowTotalLength/20.0;
    float arrowConeHeight = arrowTotalLength/5.0;
    float arrowConeRadius = 2.0*arrowRadius;
    float arrowCylinderLength = arrowTotalLength - arrowConeHeight;
    float arrowDistanceFromCell = 1.2*radius;
    //glRotatef(-90.0, cos(PI*(angle_+90.0)/180.0), sin(PI*(angle_+90.0)/180.0), 0.0);
    glTranslatef(0.0, 0.0, arrowDistanceFromCell + arrowConeHeight);
    gluCylinder(arrowCylinder_, arrowRadius, arrowRadius, arrowCylinderLength, 32, 32);
    glRotatef(180.0, 1.0, 0.0, 0.0);
    glutSolidCone(arrowConeRadius, arrowConeHeight, 32, 32);
    glPopMatrix();
  }

}

//------------------------------------------------------------------------------

void GraphicsCellGL::renderInLocalFrame() const
{

  float radius = GraphicsCellBase::getCellHeight()/2.0;
  float length = GraphicsCellBase::getCellLength() - 2.0*radius;

  // render the capped cylinder.
  glPushMatrix();
  glColor3ubv(color_);
  glTranslatef(0.0f, 0.0f, -length/2.0);
  gluCylinder(cylinder_, radius, radius, length, 32, 32);
  glutSolidSphere(radius, 32, 32);
  glTranslatef(0.0f, 0.0f, length);
  glutSolidSphere(radius, 32, 32);
  glPopMatrix();
}

//------------------------------------------------------------------------------

#endif // GUI
