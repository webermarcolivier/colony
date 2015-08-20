/***************************************************************************//**
 * Project: Colony
 *
 * \file    GraphicsCellODE.h
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

#ifndef GRAPHICSCELLODE_H
#define GRAPHICSCELLODE_H

#include "debug.h"

// standard C++ header files

// libraries header files
#include <ode/ode.h>
#include <blitz/tinyvec2.h>

// user header files
#include "GraphicsCellBase.h"
class SpatialIntegratorODE;

// namespaces
using blitz::TinyVector;



/**
  * %Cell as ODE object definition.
  */
class GraphicsCellODE : virtual public GraphicsCellBase
{

public:

  GraphicsCellODE();
  GraphicsCellODE(const GraphicsCellODE &cell);
  ~GraphicsCellODE();

  void initialize(dWorldID world, dSpaceID space, SpatialIntegratorODE* spatialIntegratorODE);
  void copy(const GraphicsCellODE& cell);

  dWorldID getWorld() const;
  dSpaceID getSpace() const;
  dBodyID getBody();
  double getAngle() const;
  TinyVector<double,3> getPositionODE() const;
  /**
    * Correct the angular velocity and rotation matrix to verify the xy plane constrain.
    */
  void alignToZAxis();
  void limitVelocity(float maxVelocitySquared, float maxRotVelocitySquared);


protected:

  void setPos(float x, float y, float z);
  void setAngle(float angle);
  void setCellLength(float length);
  void setCellHeight(float height);
  void updateGeometry();


private:

  //TinyVector<double,2> getPosition() const;

  dWorldID world_;
  dSpaceID space_;
  SpatialIntegratorODE* spatialIntegratorODE_;
  bool isAddedInSpatialIntegratorODECellsList_;
  /**
    * ODE body (mass, inertia).
    */
  dBodyID body_;
  /**
    * ODE geometry (shape).
    */
  dGeomID geom_;
  /**
    * A plane joint for restricting the cells movement to the xy plane.
    */
  dJointID planeJointID_;

};

#endif // GRAPHICSCELLODE_H
