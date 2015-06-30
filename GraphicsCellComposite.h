/***************************************************************************//**
 * Project: Colony
 *
 * \file    GraphicsCellComposite.h
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

#ifndef GRAPHICSCELLCOMPOSITE_H
#define GRAPHICSCELLCOMPOSITE_H

#define PI 3.14159265

#include "debug.h"

// user header files
#include "GraphicsCellBase.h"
#include "GraphicsCellODE.h"
#ifdef GUI
  #include "GraphicsCellQt.h"
  #include <QColor>
#endif
class GraphicsCellCompositeParam;


// namespaces

#ifdef GUI
/**
  * %Cell class packing together all the graphical representation of the cell (Qt, OpenGl, etc).
  */
  class GraphicsCellComposite : public GraphicsCellODE, public GraphicsCellQt
#else
  class GraphicsCellComposite : public GraphicsCellODE
#endif
{

public:

  GraphicsCellComposite();
  GraphicsCellComposite(const GraphicsCellComposite &cell);
  ~GraphicsCellComposite();

  void initialize(GraphicsCellCompositeParam& p);
  void copy(const GraphicsCellComposite& c);

  #ifdef GUI
    void setColor(QColor color);
  #endif


protected:

  void setPosition(float x, float y, float z);
  void setAngle(float angle);
  void setCellLength(float length);
  void setCellHeight(float height);

};

#include "Param/GraphicsCellCompositeParam.h"

#endif // GRAPHICSCELLCOMPOSITE_H
