/***************************************************************************//**
 * Project: Colony
 *
 * \file    ODESimulation.cpp
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

#include "ODESimulation.h"

//------------------------------------------------------------------------------

namespace
{
  void nearCallback(void *data, dGeomID o0, dGeomID o1)
  {
    reinterpret_cast<ODESimulation*>(data)->handleCollisionBetween(o0,o1);
  }
}

//------------------------------------------------------------------------------

ODESimulation::ODESimulation()
{
  dInitODE();
  world_ = dWorldCreate();
  space_ = dHashSpaceCreate(0);
  dWorldSetGravity (world_, 0.0f, 0.0f, 0.0f);
  contactgroup_ = dJointGroupCreate(0);

  //nbSecondsByStep_ = 0.001;
  nbSecondsByStep_ = 0.001;
  //contactMode_ = dContactSoftCFM | dContactSoftERP;
  contactMode_ = dContactSoftERP;
  erp_ = 0.2;
  cfm_ = 0.0005;
  mu_ = 0.2;
  interactionForce_ = 20.0f;
  //centralForce0_ = 400.0f;
  centralForce0_ = 100.0f;
  centralForce_ = 400.0f;
  //viscousCoefficient_ = 200.0;
  viscousCoefficient_ = 50.0;
  //rotationalViscousCoefficient_ = 600.0;
  rotationalViscousCoefficient_ = 200.0;
  aligningAngle_ = 35.0;
  aligningTorqueCoeff_ = 5.0;
  aligningRotationalViscousCoeff_ = 50.0;
  brownianForce_ = 0.0;
  //maxVelocitySquared_ = 0.001*0.001;
  maxVelocitySquared_ = 0.01*0.01;
  //maxRotVelocitySquared_ = 0.05*0.05;
  maxRotVelocitySquared_ = 0.1*0.1;

  // Set the damping coefficients of the ODE library in the world
  dWorldSetDamping ( world_, 0.1, 0.2);
  dWorldSetLinearDampingThreshold ( world_, 0.001);
  dWorldSetAngularDampingThreshold ( world_, 0.001);
  dWorldSetMaxAngularSpeed	(	world_, 1.0 );
  dWorldSetContactMaxCorrectingVel	(	world_, 4.0);

  init();
}

//------------------------------------------------------------------------------

void ODESimulation::init()
{
}

//------------------------------------------------------------------------------

ODESimulation::~ODESimulation()
{
  dSpaceDestroy(space_);
  dWorldDestroy(world_);
  dCloseODE();
}

//------------------------------------------------------------------------------

dWorldID ODESimulation::getWorld() const
{
  return world_;
}

//------------------------------------------------------------------------------

void ODESimulation::computeStep(float timeStep, QList<GraphicsCell*> *cellsList,
                                double colonySize, bool equilibration,
                                double averageCellLength, double cellHeight)
{
  // We tweak the central force intensity in function of the radius of the colony.
  {
    double colonyRadius = computeColonyRadius(cellsList, colonySize);
    //cout << "colony radius = " << colonyRadius << endl;
    int n = cellsList->size();
    // This is more or less the radius of a totally packed colony
    double denseColonyRadius = sqrt(1.1*averageCellLength*cellHeight*n/3.141592654);
    //cout << "dense colony radius = " << denseColonyRadius << endl;

    // correct central force in function of colony radius
    if ( colonyRadius <  denseColonyRadius/1.414213562 )
    {
      centralForce_ = 0.5 * centralForce0_;
    }
    else
    {
      centralForce_ = (colonyRadius/denseColonyRadius)*(colonyRadius/denseColonyRadius)*
                      (colonyRadius/denseColonyRadius)*centralForce0_;
    }
    //cout << "centralForce_ = " << centralForce_ << endl;
  }

  int nbStepsToPerform = static_cast<int>(timeStep/nbSecondsByStep_);

  // Make these steps to advance world time
  for (int i=0;i<nbStepsToPerform;++i)
  {
    // Add forces
    bool lowFriction = equilibration;
    computeForces(cellsList, colonySize, lowFriction);
    // Detect collision
    dSpaceCollide(space_, this, &nearCallback);
    // Step world
    dWorldQuickStep(world_, nbSecondsByStep_);
    //dWorldStep(world_, nbSecondsByStep_);
    foreach (GraphicsCell *o, *cellsList)
    {
      // Apply correction from 2D constraint
      static_cast<GraphicsCellODE*>(o)->ODE_AlignToZAxis();
    }

    if (!equilibration)
    {
      foreach (GraphicsCell *o, *cellsList)
      {
        // Apply maximum velocity limit
        //static_cast<GraphicsCellODE*>(o)->ODE_limitVelocity(maxVelocitySquared_,maxRotVelocitySquared_);
      }
    }

    // Remove all temporary collision joints now that the world has been stepped
    dJointGroupEmpty(contactgroup_);
  }
}

//------------------------------------------------------------------------------

void ODESimulation::computeStepWithAligningForce(float timeStep, QList<GraphicsCell*> *cellsList, double colonySize)
{
  /* For an obscure reason the aligning force produce an error in ODE:
     "ODE INTERNAL ERROR 1: assertion "bNormalizationResult" failed in _dNormalize4() [../../include/ode/odemath.h]" */

  static float nbSecondsByStep = 0.01;
  int nbStepsToPerform = static_cast<int>(timeStep/nbSecondsByStep);

  // Make these steps to advance world time
  for (int i=0;i<nbStepsToPerform;++i)
  {
    // Add forces
    computeAligningForces(cellsList, colonySize);
    computeForces(cellsList, colonySize, false);
    // Detect collision
    dSpaceCollide(space_, this, &nearCallback);
    // Step world
    dWorldQuickStep(world_, nbSecondsByStep);
    //dWorldStep(world_, nbSecondsByStep);
    foreach (GraphicsCell *o, *cellsList)
    {
      // Apply correction from 2D constraint
      static_cast<GraphicsCellODE*>(o)->ODE_AlignToZAxis();
    }

    // Remove all temporary collision joints now that the world has been stepped
    dJointGroupEmpty(contactgroup_);
  }
}

//------------------------------------------------------------------------------

void ODESimulation::handleCollisionBetween(dGeomID o0, dGeomID o1)
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

void ODESimulation::computeForces(QList<GraphicsCell*> *cellsList, double colonySize, bool lowFriction)
{
  int n = cellsList->size();
  float cutoff = colonySize/10.0;
  float forceMax = interactionForce_;
  float center[3] = {0.0*colonySize,0.0*colonySize,0.0};
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

    GraphicsCellODE *cell_i = static_cast<GraphicsCellODE*>((*cellsList)[i]);
    const dReal *pos_i = dBodyGetPosition(cell_i->ODE_getBody());

    // Interaction force.
    /*
    for (int j=0; j<i; ++j)
    {
      GraphicsCellODE *cell_j = static_cast<GraphicsCellODE*>((*cellsList)[j]);
      const dReal *pos_j = dBodyGetPosition(cell_j->ODE_getBody());

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
      dBodyAddForce(cell_i->ODE_getBody(),  force[0],  force[1],  force[2]);
      dBodyAddForce(cell_j->ODE_getBody(), -force[0], -force[1], -force[2]);
    }
    */

    // Central force.
    {
      float force[3], x[3], r, r2;
      r2 = 0.0;
      for (int k=0; k<3; ++k)
      {
        x[k] = center[k] - pos_i[k];
        r2 += x[k]*x[k];
      }
      r = sqrt(r2);
      if (r > 0)
      {
        for (int k=0; k<3; ++k)
        {
          x[k] /= r;
          //force[k] = max(1.0,r/(colonySize/4.0)) * centralForce_ * x[k];
          //force[k] = centralForce_ * x[k];
          force[k] = max(0.0, cos((3.141592654*r)/(2.*colonySize))) * centralForce_ * x[k];
        }
        dBodyAddForce(cell_i->ODE_getBody(),  force[0],  force[1],  force[2]);
      }
    }

    // Viscous force.
    {
      float force[3];
      const dReal *velocity = dBodyGetLinearVel(cell_i->ODE_getBody());
      for (int k=0; k<3; ++k)
      {
        force[k] = - viscousCoefficient * velocity[k];
      }
      dBodyAddForce(cell_i->ODE_getBody(),  force[0],  force[1],  force[2]);
    }

    // Rotational viscous force.
    {
      float torque[3];
      const dReal *angularVelocity = dBodyGetAngularVel(cell_i->ODE_getBody());
      torque[0] = 0.0;
      torque[1] = 0.0;
      torque[2] = - rotationalViscousCoefficient* angularVelocity[2];
      dBodyAddTorque(cell_i->ODE_getBody(),  torque[0],  torque[1],  torque[2]);
    }
  }
}

//------------------------------------------------------------------------------

void ODESimulation::computeAligningForces(QList<GraphicsCell*> *cellsList, double colonySize)
{
  int n = cellsList->size();
  for (int i=0; i<n; ++i)
  {

    GraphicsCellODE *cell_i = static_cast<GraphicsCellODE*>((*cellsList)[i]);
    const dReal *pos_i = dBodyGetPosition(cell_i->ODE_getBody());

    // Aligning Force
    // alpha0 is in [0,180).
    {
      float angle = cell_i->getAngle();
      float alpha0 = aligningAngle_;
      float alpha1;
      if (angle < 0)
      {
        angle = 360.0 + angle;
      }
      if ( angle - alpha0 > 90.0 )
      {
        alpha1 = angle - 180.0 - alpha0;
      }
      else if ( angle - alpha0 < -90.0 )
      {
        alpha1 = angle + 180.0 - alpha0;
      }

      float torque[3];
      torque[0] = 0.0;
      torque[1] = 0.0;
      torque[2] = - aligningTorqueCoeff_* alpha1;
      dBodyAddTorque(cell_i->ODE_getBody(),  torque[0],  torque[1],  torque[2]);
    }

    // Rotational viscous force.
    {
      float torque[3];
      const dReal *angularVelocity = dBodyGetAngularVel(cell_i->ODE_getBody());
      torque[0] = 0.0;
      torque[1] = 0.0;
      torque[2] = - aligningRotationalViscousCoeff_* angularVelocity[2];
      dBodyAddTorque(cell_i->ODE_getBody(),  torque[0],  torque[1],  torque[2]);
    }

    // Brownian force (to facilitate the equilibration).
    {
      float force[3], x[3];
      // Random direction.
      float angle = 2.0*PI*RandomNumberGenerator::getUniform();
      x[0] = cos(angle);
      x[1] = sin(angle);
      x[2] = 0.0;
      for (int k=0; k<3; ++k)
      {
        force[k] = brownianForce_ * x[k];
      }
      dBodyAddForce(cell_i->ODE_getBody(),  force[0],  force[1],  force[2]);
    }
  }
}

//------------------------------------------------------------------------------

double ODESimulation::computeColonyRadius(QList<GraphicsCell*> *cellsList, double colonySize)
{
  int n = cellsList->size();
  float center[3] = {0.0*colonySize,0.0*colonySize,0.0};
  double r2Sum = 0.0;
  double colonyRadius;

  for (int i=0; i<n; ++i)
  {
    GraphicsCellODE *cell_i = static_cast<GraphicsCellODE*>((*cellsList)[i]);
    const dReal *pos_i = dBodyGetPosition(cell_i->ODE_getBody());

    float x[3], r2;
    r2 = 0.0;
    for (int k=0; k<3; ++k)
    {
      x[k] = center[k] - pos_i[k];
      r2 += x[k]*x[k];
    }
    r2Sum = r2Sum + r2;
  }

  colonyRadius = sqrt(2.0*r2Sum/n);
  return colonyRadius;
}

//------------------------------------------------------------------------------
