/***************************************************************************//**
 * Project: Colony
 *
 * \file    OpenGLViewer.cpp
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

#include "OpenGLViewer.h"

#ifdef GUI

#include "GraphicsCellGL.h"

using namespace std;
using namespace qglviewer;

//------------------------------------------------------------------------------

OpenGLViewer::OpenGLViewer(QWidget *parent)
  : QGLViewer(QGLFormat(QGL::SampleBuffers |
                        QGL::DoubleBuffer |
                        QGL::DepthBuffer |
                        QGL::Rgba |
                        QGL::AlphaChannel |
                        QGL::StencilBuffer),parent)
{
  cellsList_.clear();
  resize(800,600);
  setWindowTitle("Quorum Sensing Simulation: OpenGLViewer");
  init();
}

//------------------------------------------------------------------------------

OpenGLViewer::~OpenGLViewer()
{}

//------------------------------------------------------------------------------

void OpenGLViewer::init()
{
  setSceneRadius(100.0);
  Vec center(0.0,0.0,0.0);
  setSceneCenter(center);

  // Try to restore previous state or focus the whole scene
  //if (!restoreStateFromFile())
  {
    showEntireScene();
  }
  setAxisIsDrawn(false);

  // Setup OpenGL
  GLfloat ambient[] = { 0.5f, 0.5f, 0.5f };
  GLfloat diffuse[] = { 0.5f, 0.5f, 0.5f , 1.0f};
  GLfloat specular[] = { 0.0f, 0.0f, 0.0f , 1.0f};
  glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
  glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
}

//------------------------------------------------------------------------------

void OpenGLViewer::draw()
{
  // Draw the objects
  foreach (GraphicsCellGL *cell, cellsList_)
  {
    cell->render();
  }

  /*
  // Draw the xy plane (with an offset of cellHeight/2.0)
  glEnable(GL_DEPTH_TEST);
  glShadeModel(GL_FLAT);
  glPushMatrix();
  glBegin(GL_QUADS);
  glColor3f(1.0, 1.0, 1.0);
  float l = colonySize_;
  float z = - cellHeight0_ / 2.0;
  glVertex3f(-l/2.,-l/2., z);
  glVertex3f(-l/2., l/2., z);
  glVertex3f( l/2., l/2., z);
  glVertex3f( l/2.,-l/2., z);
  glEnd();
  glPopMatrix();
  */
}

//------------------------------------------------------------------------------

QList<GraphicsCellGL*>* OpenGLViewer::getCellsList()
{
  return &cellsList_;
}

//------------------------------------------------------------------------------

#endif // GUI
