/***************************************************************************//**
 * Project: Colony
 *
 * \file    ODESimulation.h
 * \author  Marc Weber\n
 *          The Si.M.Bio.Sys. Group (CosmoLab)\n
 *          Parc Científic de Barcelona\n
 *          Barcelona, Spain.\n
 *          http://thesimbiosys.isgreat.org
 * \version 0.1
 * \date    11/2009
 *
 *          Copyright 2009 by Marc Weber
 ******************************************************************************/

#ifndef ODESIMULATION_H
#define ODESIMULATION_H

#include <iostream>
#include <QTime>
#include <QVector>
#include <ode/ode.h>
#include <math.h>
#include "RandomNumberGenerator.h"

using std::cout;
using std::endl;


class GraphicsCell;
class GraphicsCellODE;


/**
  * Open Dynamics Engine (ODE) simulation.
  * Contains and defines the ODE world and the basics and the ODE simulation:
  * simulation step method, computes forces between objects, handle collision method.
  */
class ODESimulation
{

public:

  ODESimulation();
  ~ODESimulation();

  dWorldID getWorld() const;
  void init();
  void handleCollisionBetween(dGeomID o0, dGeomID o1);

protected:

  dWorldID world_;
  dSpaceID space_;
  virtual void computeStep(float timeStep, QList<GraphicsCell*> *cellsList,
                           double colonySize, bool equilibration,
                           double averageCellLength, double cellHeight);
  void computeStepWithAligningForce(float timeStep, QList<GraphicsCell*> *cellsList, double colonySize);
  void computeForces(QList<GraphicsCell*> *cellsList, double colonySize, bool lowFriction);
  void computeAligningForces(QList<GraphicsCell*> *cellsList, double colonySize);
  double computeColonyRadius(QList<GraphicsCell*> *cellsList, double colonySize);

private:

  dJointGroupID contactgroup_;
  QVector<dGeomID> planes_;
  QTime time_;
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

#include "GraphicsCell.h"
#include "GraphicsCellODE.h"

#endif // ODESIMULATION_H
