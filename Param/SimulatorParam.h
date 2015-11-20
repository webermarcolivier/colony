/***************************************************************************//**
 * Project: Colony
 *
 * \file    SimulatorParam.h
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

#ifndef SIMULATORPARAM_H
#define SIMULATORPARAM_H

#include "../compilation_options.h"

// standard C++ header files

// libraries header files

// user header files

// namespaces



class SimulatorParam
{
public:

  double simulationTime_;

  double timeMeshSliceLength_;

  double chemicalLangevinTimeStep_;

  int nTrajectory_;

  double GFPconcentrationThreshold_;

  bool isFirstPassageTimeDirectionOnToOff_;

  bool isEnabledComputeSpatialDynamics_;

  bool stopSimulationWhenReachingMaximumVolume_;

};

#endif // SIMULATORPARAM_H
