/***************************************************************************//**
 * Project: Colony
 *
 * \file    SpatialIntegratorODE.h
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

#ifndef SPATIALINTEGRATORODE_H
#define SPATIALINTEGRATORODE_H

#include "compilation_options.h"

// standard C++ header files
#include <iostream>
#include <limits> ///< numerical limits for each of the fundamental types

// libraries header files
#include "QList"
#include "QVector"
#include "blitz/array.h"
#include "blitz/tinyvec2.h"
//#include "blitz/tinyvec-et.h" // file does not exist in last version of Blitz++
#include <ode/ode.h>

// user header files
class Simulator;
class GraphicsCellODE;

// namespaces
using std::cout;
using std::ostream;
using std::endl;
using std::numeric_limits;
using blitz::Array;
using blitz::Range;
using blitz::TinyVector;



/**
  * Open Dynamics Engin (ODE) for simulation of rigid body dynamics.
  */
class SpatialIntegratorODE
{
public:

  /**
   * Initiate the IntegratorContext with the pointer to the Simulator object
   * as argument.
   */
  SpatialIntegratorODE(Simulator* simulatorPtr);

  ~SpatialIntegratorODE();

  dWorldID getWorld() const;
  dSpaceID getSpace() const;
  QList<GraphicsCellODE*>* getCellsList();
  void init();
  void handleCollisionBetween(dGeomID o0, dGeomID o1);
  void integrate(double timeStep);
  void setIsEquilibrationSteps(bool isEquilibrationSteps);

protected:

  dWorldID world_;
  dSpaceID space_;
  void computeForces(bool lowFriction);
  TinyVector<double,4> computeColonyCenterAndRadius();

private:

  Simulator* simulatorPtr_;
  bool isEquilibrationSteps_;

  QList<dGeomID> wallsODEList_;
  QList<GraphicsCellODE*> cellsList_;

  dJointGroupID contactgroup_;
  QVector<dGeomID> planes_;
  //QTime time_;
  float nbSecondsByStep_;
  float interactionForce_;
  float centralForce0_;
  float centralForce_;
  float viscousCoefficient_;
  float rotationalViscousCoefficient_;
  float aligningAngle_;
  float aligningTorqueCoeff_;
  float aligningRotationalViscousCoeff_;
  float brownianForce_;
  float maxVelocitySquared_;
  float maxRotVelocitySquared_;
  float erp_;
  float cfm_;
  float mu_;
  int contactMode_;

};

#include "Simulator.h"

#endif // SPATIALINTEGRATORODE_H
