/***************************************************************************//**
 * Project: Colony
 *
 * \file    SpatialIntegratorODE.cpp
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

#include "SpatialIntegratorODE.h"
#include "GraphicsCellODE.h"


namespace
{
  void nearCallback(void *data, dGeomID o0, dGeomID o1)
  {
    reinterpret_cast<SpatialIntegratorODE*>(data)->handleCollisionBetween(o0,o1);
  }
}

//------------------------------------------------------------------------------

SpatialIntegratorODE::SpatialIntegratorODE(Simulator* simulatorPtr)
  : simulatorPtr_(simulatorPtr)
{
  cellsList_.clear();
  init();
}

//------------------------------------------------------------------------------

SpatialIntegratorODE::~SpatialIntegratorODE()
{
  dSpaceDestroy(space_);
  dWorldDestroy(world_);
  dCloseODE();
}

//------------------------------------------------------------------------------

void SpatialIntegratorODE::init()
{
  isEquilibrationSteps_ = false;

  dInitODE();
  world_ = dWorldCreate();
  space_ = dHashSpaceCreate(0);
  dWorldSetGravity (world_, 0.0f, 0.0f, 0.0f);
  contactgroup_ = dJointGroupCreate(0);

  //nbSecondsByStep_ = 0.001;
  nbSecondsByStep_ = 0.002;
  //contactMode_ = dContactSoftCFM | dContactSoftERP;
  contactMode_ = dContactSoftERP;
  erp_ = 0.1;//0.2
  cfm_ = 0.02;//0.0
  mu_ = 0;
  interactionForce_ = 0.0;
  //centralForce0_ = 400.0f;
  centralForce0_ = 10.0;
  //viscousCoefficient_ = 200.0;
  viscousCoefficient_ = 50.0;
  //rotationalViscousCoefficient_ = 600.0;
  rotationalViscousCoefficient_ = 50.0;
  aligningAngle_ = 35.0;
  aligningTorqueCoeff_ = 5.0;
  aligningRotationalViscousCoeff_ = 50.0;
  brownianForce_ = 0.0;
  //maxVelocitySquared_ = 0.001*0.001;
  //maxRotVelocitySquared_ = 0.05*0.05;

  // Set the damping coefficients of the ODE library in the world
  //dWorldSetDamping ( world_, 0.1, 0.2);
  //dWorldSetLinearDampingThreshold ( world_, 0.001);
  //dWorldSetAngularDampingThreshold ( world_, 0.001);
  dWorldSetMaxAngularSpeed	(	world_, 1.0 );
  dWorldSetContactMaxCorrectingVel	(	world_, 0.8);
}

//------------------------------------------------------------------------------

dWorldID SpatialIntegratorODE::getWorld() const
{
  return world_;
}

//------------------------------------------------------------------------------

dSpaceID SpatialIntegratorODE::getSpace() const
{
  return space_;
}

//------------------------------------------------------------------------------

QList<GraphicsCellODE*>* SpatialIntegratorODE::getCellsList()
{
  return &cellsList_;
}

//------------------------------------------------------------------------------

void SpatialIntegratorODE::setIsEquilibrationSteps(bool isEquilibrationSteps)
{
  isEquilibrationSteps_ = isEquilibrationSteps;
}

//------------------------------------------------------------------------------

void SpatialIntegratorODE::integrate(double timeStep)
{
  // We tweak the central force intensity in function of the radius of the colony.
  {
    TinyVector<double,4> centerAndRadius = computeColonyCenterAndRadius();

    int n = cellsList_.size();

    // This is more or less the radius of a totally packed colony
    double cellLength0 = simulatorPtr_->cellCollection_[0].getCellLength0();
    double cellHeight0 = simulatorPtr_->cellCollection_[0].getCellHeight0();
    double denseColonyRadius = sqrt((1.44*cellLength0*cellHeight0*n)/3.141592654) / 2.0;

    // Correct central force in function of colony radius
    centralForce_ = (centerAndRadius(3)/denseColonyRadius)*
                    (centerAndRadius(3)/denseColonyRadius)*centralForce0_;
  }

  int nbStepsToPerform = static_cast<int>(timeStep/nbSecondsByStep_);

  // Make these steps to advance world time
  for (int i=0;i<nbStepsToPerform;++i)
  {
    // Add forces
    computeForces(isEquilibrationSteps_);
    // Detect collision
    dSpaceCollide(space_, this, &nearCallback);
    // Step world
    dWorldQuickStep(world_, nbSecondsByStep_);
    //dWorldStep(world_, nbSecondsByStep_);

    foreach (GraphicsCellODE *cell, cellsList_)
    {
      // Apply correction from 2D constraint
      cell->alignToZAxis();
    }

    if (!isEquilibrationSteps_)
    {
      foreach (GraphicsCellODE *cell, cellsList_)
      {
        // Apply maximum velocity limit
        //cell->limitVelocity(maxVelocitySquared_,maxRotVelocitySquared_);
      }
    }

    // Remove all temporary collision joints now that the world has been stepped
    dJointGroupEmpty(contactgroup_);
  }

  // Send the new position and angle of the cells back to the simulator
  for (int i=0; i < cellsList_.size(); ++i)
  {
    TinyVector<double,3> pos = cellsList_.at(i)->getPositionODE();
    simulatorPtr_->cellCollection_[i].Cell::setPosition(pos);
    double angle = cellsList_.at(i)->getAngle();
    simulatorPtr_->cellCollection_[i].Cell::setAngle(angle);
  }

}

//------------------------------------------------------------------------------

void SpatialIntegratorODE::handleCollisionBetween(dGeomID o0, dGeomID o1)
{
    // Create an array of dContact objects to hold the contact joints
    static const int MAX_CONTACTS = 10;
    dContact contact[MAX_CONTACTS];

    /*
    By adjusting the values of ERP and CFM, you can achieve various effects.
    For example you can simulate springy constraints, where the two bodies
    oscillate as though connected by springs. Or you can simulate more spongy
    constraints, without the oscillation. In fact, ERP and CFM can be selected
    to have the same effect as any desired spring and damper constants. If you
    have a spring constant kp and damping constant kd, then the corresponding
    ODE constants are:
    ERP = h kp / (h kp + kd)
    CFM = 1 / (h kp + kd)
    where h is the time step
    or
    (spring) kp = ECM / CFM h
    (damping) kd = (ECM - 1) / CFM
    */
    for (int i = 0; i < MAX_CONTACTS; i++)
    {
      contact[i].surface.mode = contactMode_;
      contact[i].surface.soft_erp = erp_;
      contact[i].surface.soft_cfm = cfm_;
      contact[i].surface.mu = mu_;
    }
    if (int numc = dCollide(o0, o1, MAX_CONTACTS, &contact[0].geom, sizeof(dContact)))
    {
      // Get the dynamics body for each geom
      dBodyID b1 = dGeomGetBody(o0);
      dBodyID b2 = dGeomGetBody(o1);
      // To add each contact point found to our joint group we call dJointCreateContact which is just one of the many
      // different joint types available.
      for (int i = 0; i < numc; i++)
      {
        // dJointCreateContact needs to know which world and joint group to work with as well as the dContact
        // object itself. It returns a new dJointID which we then use with dJointAttach to finally create the
        // temporary contact joint between the two geom bodies.
        dJointID c = dJointCreateContact(world_, contactgroup_, contact + i);
        dJointAttach(c, b1, b2);
      }
    }
}

//------------------------------------------------------------------------------

void SpatialIntegratorODE::computeForces(bool lowFriction)
{
  int n = cellsList_.size();
  double cellLength0 = simulatorPtr_->cellCollection_[0].getCellLength0();
  float cutoff = 5.0*cellLength0;
  float forceMax = interactionForce_;
  float center[3] = {0.0, 0.0, 0.0};
  float viscousCoefficient;
  float rotationalViscousCoefficient;
  if (!lowFriction)
  {
    viscousCoefficient = viscousCoefficient_;
    rotationalViscousCoefficient = rotationalViscousCoefficient_;
  } else {
    viscousCoefficient = 0.05*viscousCoefficient_;
    rotationalViscousCoefficient = 0.05*rotationalViscousCoefficient_;
  }

  for (int i=0; i<n; ++i)
  {
    TinyVector<double,3> pos_i = cellsList_.at(i)->getPositionODE();

    // Interaction force.
    /*
    for (int j=0; j<i; ++j)
    {
      TinyVector<double,3> pos_j = cellsList_.at(j)->getPosition();

      float force[3], x[3], r, r2;
      r2 = 0.0;
      for (int k=0; k<3; ++k)
      {
        x[k] = pos_j[k] - pos_i[k];
        r2 += x[k]*x[k];
      }
      r = sqrt(r2);
      for (int k=0; k<3; ++k)
      {
        x[k] /= r;
        if (r > cutoff)
        {
          force[k] = 0.0;
        } else {
          force[k] = (forceMax/cutoff)*(cutoff - r) * x[k];
        }
      }
      dBodyAddForce(cellsList_.at(i)->getBody(),  force[0],  force[1],  force[2]);
      dBodyAddForce(cellsList_.at(j)->getBody(), -force[0], -force[1], -force[2]);
    }
    */

    // Central force.
    {
      float force[3], x[3], r, r2;
      r2 = 0.0;
      for (int k=0; k<3; ++k)
      {
        x[k] = center[k] - pos_i(k);
        r2 += x[k]*x[k];
      }
      r = sqrt(r2);
      for (int k=0; k<3; ++k)
      {
        if (r != 0)
        {
          x[k] /= r;
        } else {
          x[k] = 0.0;
        }
        //force[k] = max(1.0,r/(colonySize/4.0)) * centralForce_ * x[k];
        force[k] = centralForce_ * x[k];
        //force[k] = max(0.0, cos((3.141592654*r)/(2.0*20.0*cellLength0))) * centralForce_ * x[k];
      }
      dBodyAddForce(cellsList_[i]->getBody(),  force[0],  force[1],  0.0);
    }

    // Viscous force.
    {
      float force[3];
      const dReal *velocity = dBodyGetLinearVel(cellsList_.at(i)->getBody());
      for (int k=0; k<3; ++k)
      {
        force[k] = - viscousCoefficient * velocity[k];
      }
      dBodyAddForce(cellsList_[i]->getBody(),  force[0],  force[1],  0.0);
    }

    // Rotational viscous force.
    {
      float torque[3];
      const dReal *angularVelocity = dBodyGetAngularVel(cellsList_.at(i)->getBody());
      torque[0] = 0.0;
      torque[1] = 0.0;
      torque[2] = - rotationalViscousCoefficient* angularVelocity[2];
      dBodyAddTorque(cellsList_[i]->getBody(),  torque[0],  torque[1],  torque[2]);
    }
  }
}

//------------------------------------------------------------------------------

TinyVector<double,4> SpatialIntegratorODE::computeColonyCenterAndRadius()
{
  TinyVector<double,4> centerAndRadius;
  int n = cellsList_.size();
  double r2Sum = 0.0;
  double colonyRadius;

  // Compute the center of mass
  TinyVector<double,3> center = (0.0, 0.0, 0.0);
  QList<GraphicsCellODE*>::const_iterator it;
  for (it = cellsList_.constBegin(); it != cellsList_.constEnd(); ++it)
  {
    TinyVector<double,3> pos = (*it)->getPositionODE();
    center += pos;
  }
  center /= n;
  centerAndRadius(0) = center(0);
  centerAndRadius(1) = center(1);
  centerAndRadius(2) = center(2);

  // Compute the radius of gyration
  double radius = 0.0;
  for (it = cellsList_.constBegin(); it != cellsList_.constEnd(); ++it)
  {
    TinyVector<double,3> pos = (*it)->getPositionODE();
    radius += sum( (pos - center)*(pos - center) );
  }
  radius = sqrt(radius/n);
  centerAndRadius(3) = radius;

  return centerAndRadius;
}

//------------------------------------------------------------------------------
