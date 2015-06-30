/***************************************************************************//**
 * Project: Colony
 *
 * \file    GraphicsCellODE.cpp
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

#include "GraphicsCellODE.h"
#include "SpatialIntegratorODE.h"

//------------------------------------------------------------------------------

GraphicsCellODE::GraphicsCellODE()
  : isAddedInSpatialIntegratorODECellsList_(false)
{}

//------------------------------------------------------------------------------

GraphicsCellODE::~GraphicsCellODE()
{
  /// Delete reference in the list of cell in SpatialIntegratorODE class.
  if (isAddedInSpatialIntegratorODECellsList_)
  {
    spatialIntegratorODE_->getCellsList()->removeAll(this);
    isAddedInSpatialIntegratorODECellsList_ = false;
  }

  dJointDestroy(planeJointID_);
  dBodyDestroy (body_);
  dGeomDestroy (geom_);
}

//------------------------------------------------------------------------------

void GraphicsCellODE::initialize(dWorldID world, dSpaceID space, SpatialIntegratorODE* spatialIntegratorODE)
{
  /// - Set the ODE world and space.
  world_ = world;
  space_ = space;
  spatialIntegratorODE_ = spatialIntegratorODE;
  planeJointID_ = dJointCreatePlane2D( world_, 0);

  /// - Get the cell height from the GraphicsCell object.
  float radius = cellHeight_/2.0;
  if (radius <= 0) radius = 0.5;
  float cylinderLength = cellLength_ - 2.0*radius;
  if (cylinderLength <= 0) cylinderLength = 1.0;

  geom_ = dCreateCCylinder(space_, radius, cylinderLength);

  /// - Initiate body (mass, inertia).
  body_ = dBodyCreate(world_);
  dMass m;
  dMassSetCappedCylinderTotal(&m, 1.0f, 3, radius, cylinderLength);
  dBodySetMass(body_,&m);

  /// - Associate body and geometry.
  dGeomSetBody(geom_,body_);

  dGeomSetPosition(geom_, 0.0, 0.0, 0.0);
  dBodySetLinearVel(body_, 0.0, 0.0, 0.0);
  dBodySetAngularVel(body_, 0.0, 0.0, 0.0);

  /// - Restrict movement to the plane.
  dJointAttach( planeJointID_, body_, 0 );

  /// - Add the cell to the list in the spatialIntegratorODE
  if (!isAddedInSpatialIntegratorODECellsList_)
  {
    spatialIntegratorODE_->getCellsList()->append(this);
    isAddedInSpatialIntegratorODECellsList_ = true;
  }
}

//------------------------------------------------------------------------------

GraphicsCellODE::GraphicsCellODE(const GraphicsCellODE& cell)
{
  copy(cell);
}

//------------------------------------------------------------------------------

void GraphicsCellODE::copy(const GraphicsCellODE& cell)
{
  isAddedInSpatialIntegratorODECellsList_ = false;

  /// - Set the ODE world and space.
  world_ = cell.world_;
  space_ = cell.space_;
  spatialIntegratorODE_ = cell.spatialIntegratorODE_;
  planeJointID_ = dJointCreatePlane2D( world_, 0);

  /// - Get body position and orientation.
  const dReal *posConstPointer = dGeomGetPosition(cell.geom_);
  dReal pos[3];
  for(int i=0; i<3; ++i)
  {
    pos[i] = posConstPointer[i];
  }
  dQuaternion quat;
  dGeomGetQuaternion(cell.geom_, quat);

  /// - Get the body linear velocity and angular velocity.
  const dReal *linearVelocity = dBodyGetLinearVel(cell.body_);
  const dReal *angularVelocity = dBodyGetAngularVel(cell.body_);

  /// - Get the cell height from the GraphicsCell object.
  float radius = cell.cellHeight_/2.0;
  if (radius <= 0) radius = 0.5;
  float cylinderLength = cell.cellLength_ - 2.0*radius;
  if (cylinderLength <= 0) cylinderLength = 1.0;
  float angle = cell.angle_;

  /// - Initiate geometry (object shape).
  geom_ = dCreateCCylinder(space_, radius, cylinderLength);

  /// - Initiate body (mass, inertia).
  body_ = dBodyCreate(world_);
  dMass m;
  dMassSetCappedCylinderTotal(&m, 1.0f, 3, radius, cylinderLength);
  dBodySetMass(body_,&m);

  /// - Associate body and geometry.
  dGeomSetBody(geom_,body_);

  /// - Set geometry and body position to the same as before.
  dGeomSetPosition(geom_, pos[0], pos[1], pos[2]);
  /// - Set geometry and body rotation to the same as before.
  dGeomSetQuaternion(geom_, quat);
  /// - Set geometry and body linear velocity to the same as before.
  dBodySetLinearVel(body_, linearVelocity[0], linearVelocity[1], linearVelocity[2]);
  /// - Set geometry and body angular velocity to the same as before.
  dBodySetAngularVel(body_, angularVelocity[0], angularVelocity[1], angularVelocity[2]);

  /// - Restrict movement to the plane.
  dJointAttach( planeJointID_, body_, 0 );

  /// - Add the cell to the list in the spatialIntegratorODE
  if (!isAddedInSpatialIntegratorODECellsList_)
  {
    spatialIntegratorODE_->getCellsList()->append(this);
    isAddedInSpatialIntegratorODECellsList_ = true;
  }
}

//------------------------------------------------------------------------------

void GraphicsCellODE::setPos(float x, float y, float z)
{
  dBodySetPosition(body_, x, y, z);
}

//------------------------------------------------------------------------------

void GraphicsCellODE::setAngle(float angle)
{
  /// - Set the rotation matrix R corresponding to the orientation defined by the angle parameter.
  dMatrix3 R;
  /// - The axis of rotation is in the xy plane and is perpendicular to the cell orientation.
  ///   The angle of rotation is 90 degrees.
  dRFromAxisAndAngle(R, cos(PI*(angle_+90.0)/180.0), sin(PI*(angle_+90.0)/180.0), 0.0, PI/2.0);
  dBodySetRotation(body_, R);
}

//------------------------------------------------------------------------------

void GraphicsCellODE::setCellLength(float length)
{
  GraphicsCellBase::setCellLength(length);
  updateGeometry();
}

//------------------------------------------------------------------------------

void GraphicsCellODE::setCellHeight(float height)
{
  GraphicsCellBase::setCellHeight(height);
  updateGeometry();
}

//------------------------------------------------------------------------------

void GraphicsCellODE::updateGeometry()
{
  /// - Get body position and orientation.
  const dReal *posConstPointer = dGeomGetPosition(geom_);
  dReal pos[3];
  for(int i=0; i<3; ++i)
  {
    pos[i] = posConstPointer[i];
  }
  dQuaternion quat;
  dGeomGetQuaternion(geom_, quat);

  /// - Get the body linear velocity and angular velocity.
  const dReal *linearVelocity = dBodyGetLinearVel(body_);
  const dReal *angularVelocity = dBodyGetAngularVel(body_);

  /// - Get the cell height from the GraphicsCell object.
  float radius = cellHeight_/2.0;
  if (radius <= 0) radius = 0.5;
  float cylinderLength = cellLength_ - 2.0*radius;
  if (cylinderLength <= 0) cylinderLength = 1.0;
  float angle = angle_;

  /// - Initiate geometry (object shape).
  dGeomDestroy(geom_);
  geom_ = dCreateCCylinder(space_, radius, cylinderLength);

  /// - Initiate body (mass, inertia).
  dMass m;
  dMassSetCappedCylinderTotal(&m, 1.0f, 3, radius, cylinderLength);
  dBodySetMass(body_,&m);

  /// - Associate body and geometry.
  dGeomSetBody(geom_,body_);

  /// - Set geometry and body position to the same as before.
  dGeomSetPosition(geom_, pos[0], pos[1], pos[2]);
  /// - Set geometry and body rotation to the same as before.
  dGeomSetQuaternion(geom_, quat);
  /// - Set geometry and body linear velocity to the same as before.
  dBodySetLinearVel(body_, linearVelocity[0], linearVelocity[1], linearVelocity[2]);
  /// - Set geometry and body angular velocity to the same as before.
  dBodySetAngularVel(body_, angularVelocity[0], angularVelocity[1], angularVelocity[2]);

  /// - Restrict movement to the plane.
  dJointAttach( planeJointID_, body_, 0 );
}

//------------------------------------------------------------------------------

dWorldID GraphicsCellODE::getWorld() const
{
  return world_;
}

//------------------------------------------------------------------------------

dSpaceID GraphicsCellODE::getSpace() const
{
  return space_;
}

//------------------------------------------------------------------------------

dBodyID GraphicsCellODE::getBody()
{
  return body_;
}

//------------------------------------------------------------------------------

TinyVector<double,3> GraphicsCellODE::getPositionODE() const
{
  TinyVector<double,3> position;
  const dReal *positionVector = dGeomGetPosition (geom_);
  position(0) = positionVector[0];
  position(1) = positionVector[1];
  position(2) = positionVector[2];
  return position;
}

//------------------------------------------------------------------------------

double GraphicsCellODE::getAngle() const
{
  /// - Get the quaternion matrix of the ODE geometry.
  const dReal *quat = dBodyGetQuaternion( body_ );

  /// - Remark: a quaternion is equal to (cos(alpha/2),sin(alpha/2)*u), where u is a unit vector of the rotation axis.
  /// - Remark: the rotation axis should already be in the xy plane and the rotation angle should be 90 degrees.

  /// - Calculate the orientation of the rotation axis in the xy plane, Euler psi angle.
  float psi;
  psi = atan2(2*(quat[0]*quat[3]+quat[1]*quat[2]), 1-2*(quat[2]*quat[2]+quat[3]*quat[3]));

  /// - The orientation of the cell is perpendicular to the rotation axis.
  double angle;
  angle = psi - (PI/2.0);

  /// - Transform from radians to degrees.
  angle = angle*180 / PI;
  return angle;
}

//------------------------------------------------------------------------------

void GraphicsCellODE::alignToZAxis()
{
  const dReal *rot = dBodyGetAngularVel( body_ );
  const dReal *quat_ptr;
  dReal quat[4], quat_len;
  quat_ptr = dBodyGetQuaternion( body_ );
  quat[0] = 1.0 / sqrt(2.0);
  quat[1] = quat_ptr[1];
  quat[2] = quat_ptr[2];
  quat[3] = 0;
  quat_len = quat[1] * quat[1] + quat[2] * quat[2] ;
  quat[1] = 0.5 * quat[1] / quat_len;
  quat[2] = 0.5 * quat[2] / quat_len;
  dBodySetQuaternion( body_, quat );
  dBodySetAngularVel( body_, 0, 0, rot[2] );
}

//------------------------------------------------------------------------------

void GraphicsCellODE::limitVelocity(float maxVelocitySquared, float maxRotVelocitySquared)
{
  const dReal *vel = dBodyGetLinearVel( body_ );
  float velNormSquared = vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2];
  if ( velNormSquared > maxVelocitySquared)
  {
    dReal newVel[3];
    newVel[0] = sqrt(maxVelocitySquared/velNormSquared)*vel[0];
    newVel[1] = sqrt(maxVelocitySquared/velNormSquared)*vel[1];
    newVel[2] = sqrt(maxVelocitySquared/velNormSquared)*vel[2];
    dBodySetLinearVel  (body_, newVel[0], newVel[1], newVel[2]);
  }

  const dReal *rotVel = dBodyGetAngularVel( body_ );
  float rotVelNormSquared = rotVel[0]*rotVel[0] + rotVel[1]*rotVel[1] + rotVel[2]*rotVel[2];
  if ( rotVelNormSquared > maxRotVelocitySquared)
  {
    dReal newRotVel[3];
    newRotVel[0] = sqrt(maxRotVelocitySquared/rotVelNormSquared)*rotVel[0];
    newRotVel[1] = sqrt(maxRotVelocitySquared/rotVelNormSquared)*rotVel[1];
    newRotVel[2] = sqrt(maxRotVelocitySquared/rotVelNormSquared)*rotVel[2];
    dBodySetAngularVel  (body_, newRotVel[0], newRotVel[1], newRotVel[2]);
  }
}

//------------------------------------------------------------------------------
