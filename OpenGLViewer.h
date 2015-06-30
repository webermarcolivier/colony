/***************************************************************************//**
 * Project: Colony
 *
 * \file    OpenGLViewer.h
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

#ifndef OPENGLVIEWER_H
#define OPENGLVIEWER_H

#include "debug.h"

#ifdef GUI

#include <QGLViewer/qglviewer.h>
#include <GL/glut.h>
#include <QList>

class GraphicsCellGL;


/**
  * A small openGL viewer of the GraphicsCell objects.
  * @see QGLViewer
  */
class OpenGLViewer : public QGLViewer
{

  Q_OBJECT

public:
  OpenGLViewer(QWidget *parent=NULL);
  ~OpenGLViewer();
  QList<GraphicsCellGL*>* getCellsList();

protected:

  virtual void init();
  virtual void draw();

private:

  float _aabb[6];
  QList<GraphicsCellGL*> cellsList_;

};

#endif // GUI

#endif // OPENGLVIEWER_H
