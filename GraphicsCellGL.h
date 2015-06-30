/***************************************************************************//**
 * Project: Colony
 *
 * \file    GraphicsCellGL.h
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


#ifndef GRAPHICSCELLGL_H
#define GRAPHICSCELLGL_H

#include "debug.h"

#ifdef GUI

// standard C++ header files
#include <cmath>

// libraries header files
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <stdlib.h>
#include <time.h>

// user header files
#include "GraphicsCellBase.h"

class OpenGLViewer;

// namespaces



class GraphicsCellGL : virtual public GraphicsCellBase
{

public:

  GraphicsCellGL();
  GraphicsCellGL(const GraphicsCellGL &cell);
  ~GraphicsCellGL();

  void initialize(OpenGLViewer* openGLViewer);
  void copy(const GraphicsCellGL& cell);

  void render() const;
  void renderInLocalFrame() const;


protected:

  void setPos(float x, float y, float z);
  void setAngle(float angle);
  void setCellLength(float length);
  void setCellHeight(float height);
  void setColor(int r, int g, int b);


private:

  unsigned char color_[3];
  GLUquadric *cylinder_;
  GLUquadric *arrowCylinder_;
  bool isAddedInOpenGLViewerCellsList_;
  OpenGLViewer* openGLViewer_;

};

#endif // GUI

#endif // GRAPHICSCELLGL_H
