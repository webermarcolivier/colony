/***************************************************************************//**
 * Project: Colony
 *
 * \file    GraphicsCellBase.h
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

#ifndef GRAPHICSCELLBASE_H
#define GRAPHICSCELLBASE_H

#define PI 3.14159265

#include "compilation_options.h"

// standard C++ header files

// libraries header files

// user header files

// namespaces

/**
  * Base class for the graphical representations of the cell.
  * Define virtual functions as set position, update position, etc.
  */
class GraphicsCellBase
{

protected:

  GraphicsCellBase();
  GraphicsCellBase(const GraphicsCellBase &cell);
  ~GraphicsCellBase();
  void copy(const GraphicsCellBase& c);

public:

  int getCellIndex() const;
  float getCellLength() const;
  float getCellHeight() const;
  float getAngle() const;

protected:

  void setCellLength(float length);
  void setCellHeight(float height);
  void setCellIndex(int cellIndex);
  void setPos(float x, float y, float z);
  void setAngle(float angle);

  int cellIndex_;
  float angle_;
  /**
    * Cell length in microns.
    */
  float cellLength_;
  /**
    * Cell height in microns.
    */
  float cellHeight_;
  float cellVolume_;
  float posX_;
  float posY_;
  float posZ_;
  bool hasFocus_;

};

#endif // GRAPHICSCELLBASE_H
