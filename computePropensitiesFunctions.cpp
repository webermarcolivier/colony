/***************************************************************************//**
 * Project: Colony
 *
 * \file    computePropensitiesFunctions.cpp
 * \author  Marc Weber\n
 *          The Si.M.Bio.Sys. Group (CosmoLab)\n
 *          Parc Cienti­fic de Barcelona\n
 *          Barcelona, Spain.\n
 *          http://thesimbiosys.isgreat.org
 * \version 0.1
 * \date    03/2009
 *
 *          Output file from Mathematica.\n
 *          Copyright 2009 by Marc Weber
 ******************************************************************************/

#include "compilation_options.h"

#ifndef USE_CHEMICAL_LANGEVIN

#include "computePropensitiesFunctions.h"

void computePropensitiesCell
(
  const Array<int,1>& x,
  Array<double,1>& a
)
{
  #ifdef COMPUTEPROPENSITIES_ARRAY_SIZE_CHECK
  if ( 4 != a.size()) {
  cout << "ERROR: function computePropensitiesCell(x, propensities),"
          " array \"propensities\" does not have "
          "the same size as the number of reactions";
  exit(1);
  }

  if ( 2 != x.size()) {
  cout << "ERROR: function computePropensitiesCell(x, propensities),"
          " array \"x\" does not have "
          "the same size as the number of species";
  exit(1);
  }
  #endif //COMPUTEPROPENSITIES_ARRAY_SIZE_CHECK

  a(0) = 0.9033209999999998*(0.7499999999999998 + Input::globalParameter(0) + 1.8750000000000002/(0.125 + 1.3566683085495015*pow(x(1),3)));
  a(1) = 0.9033209999999998*(0.7499999999999998 + 1.8750000000000002/(0.125 + 1.3566683085495015*pow(x(0),3)));

  a(2) = 1.*x(0);
  a(3) = 1.*x(1);


  #ifdef COMPUTEPROPENSITIES_POSITIVITY_CHECK
  int i;
  for (i=0; i<4; i++)
  {
    if ( a(i) < 0.0 )
    {
      cout << "WARNING: function computePropensitiesCell(x, propensities),"
              " negative propensity detected: a(" << i << ") = " << a(i)
           << ". Resetting a to 0.0" << endl;
      a(i) = 0.0;
    }
  }
  #endif //COMPUTEPROPENSITIES_POSITIVITY_CHECK
}


void computePropensitiesMilieu
(
  const Array<int,1>& x,
  Array<double,1>& a
)
{
  #ifdef COMPUTEPROPENSITIES_ARRAY_SIZE_CHECK
  if ( 4 != a.size()) {
  cout << "ERROR: function computePropensitiesMilieu(x, propensities),"
          " array \"propensities\" does not have "
          "the same size as the number of reactions";
  exit(1);
  }

  if ( 2 != x.size()) {
  cout << "ERROR: function computePropensitiesMilieu(x, propensities),"
          " array \"x\" does not have "
          "the same size as the number of species";
  exit(1);
  }
  #endif //COMPUTEPROPENSITIES_ARRAY_SIZE_CHECK

  a(0) = 1.*x(0);
  a(1) = 1.*x(1);

  a(2) = 0.;
  a(3) = 0.;


  #ifdef COMPUTEPROPENSITIES_POSITIVITY_CHECK
  int i;
  for (i=0; i<4; i++)
  {
    if ( a(i) < 0.0 )
    {
      cout << "WARNING: function computePropensitiesMilieu(x, propensities),"
              " negative propensity detected: a(" << i << ") = " << a(i)
           << ". Resetting a to 0.0" << endl;
      a(i) = 0.0;
    }
  }
  #endif //COMPUTEPROPENSITIES_POSITIVITY_CHECK
}


void computePropensitiesCellTimeDependent
(
  const Array<int,1>& x,
  const double volume,
  const double volume0,
  Array<double,1>& a
)
{


  #ifdef COMPUTEPROPENSITIES_TIME_CYCLE_CHECK
  if ( volume < 0.0 || volume > 2.0*volume0 )
  {
  cout << "ERROR: function computePropensitiesCellTimeDependent(x, propensities),"
          " volume value is negative or bigger than the 2*V0.";
  exit(1);
  }
  #endif //COMPUTEPROPENSITIES_TIME_CYCLE_CHECK

  #ifdef COMPUTEPROPENSITIES_ARRAY_SIZE_CHECK
  if ( 4 != a.size()) {
  cout << "ERROR: function computePropensitiesCellTimeDependent(x, propensities),"
          " array \"propensities\" does not have "
          "the same size as the number of reactions";
  exit(1);
  }

  if ( 2 != x.size()) {
  cout << "ERROR: function computePropensitiesCellTimeDependent(x, propensities),"
          " array \"x\" does not have "
          "the same size as the number of species";
  exit(1);
  }
  #endif //COMPUTEPROPENSITIES_ARRAY_SIZE_CHECK

  a(0) = (0.9033209999999998*volume*(0.7499999999999998 + Input::globalParameter(0) + 1.8750000000000002/(0.125 + 1.3566683085495015*pow(x(1),3))))/volume0;
  a(1) = (0.9033209999999998*volume*(0.7499999999999998 + 1.8750000000000002/(0.125 + 1.3566683085495015*pow(x(0),3))))/volume0;

  a(2) = x(0);
  a(3) = x(1);


  #ifdef COMPUTEPROPENSITIES_POSITIVITY_CHECK
  int i;
  for (i=0; i<4; i++)
  {
    if ( a(i) < 0.0 )
    {
      cout << "WARNING: function computePropensitiesCellTimeDependent(x, propensities),"
              " negative propensity detected: a(" << i << ") = " << a(i)
           << ". Resetting a to 0.0" << endl;
      a(i) = 0.0;
    }
  }
  #endif //COMPUTEPROPENSITIES_POSITIVITY_CHECK
}


void computePropensitiesMilieuTimeDependent
(
  const Array<int,1>& x,
  const double volume,
  const double volume0,
  Array<double,1>& a
)
{


  #ifdef COMPUTEPROPENSITIES_TIME_CYCLE_CHECK
  if ( volume < 0.0 || volume > 2.0*volume0 )
  {
  cout << "ERROR: function computePropensitiesMilieuTimeDependent(x, propensities),"
          " volume value is negative or bigger than the 2*V0.";
  exit(1);
  }
  #endif //COMPUTEPROPENSITIES_TIME_CYCLE_CHECK

  #ifdef COMPUTEPROPENSITIES_ARRAY_SIZE_CHECK
  if ( 4 != a.size()) {
  cout << "ERROR: function computePropensitiesMilieuTimeDependent(x, propensities),"
          " array \"propensities\" does not have "
          "the same size as the number of reactions";
  exit(1);
  }

  if ( 2 != x.size()) {
  cout << "ERROR: function computePropensitiesMilieuTimeDependent(x, propensities),"
          " array \"x\" does not have "
          "the same size as the number of species";
  exit(1);
  }
  #endif //COMPUTEPROPENSITIES_ARRAY_SIZE_CHECK

  a(0) = x(0);
  a(1) = x(1);

  a(2) = 0.;
  a(3) = 0.;


  #ifdef COMPUTEPROPENSITIES_POSITIVITY_CHECK
  int i;
  for (i=0; i<4; i++)
  {
    if ( a(i) < 0.0 )
    {
      cout << "WARNING: function computePropensitiesMilieuTimeDependent(x, propensities),"
              " negative propensity detected: a(" << i << ") = " << a(i)
           << ". Resetting a to 0.0" << endl;
      a(i) = 0.0;
    }
  }
  #endif //COMPUTEPROPENSITIES_POSITIVITY_CHECK
}


void computePropensitiesCellMilieu
(
  const Array<int,1>& x1,
  const Array<int,1>& x2,
  const double volume1,
  const double volume2,
  Array<double,1>& a
)
{
  #ifdef COMPUTEPROPENSITIES_ARRAY_SIZE_CHECK
  if ( 4 != a.size()) {
  cout << "ERROR: function computePropensitiesCellMilieu(x1, x2, propensities),"
          " array \"propensities\" does not have "
          "the same size as the number of reactions";
  exit(1);
  }

  if ( 2 != x1.size()) {
  cout << "ERROR: function computePropensitiesCellMilieu(x1, x2, propensities),"
          " array \"x1\" does not have "
          "the same size as the number of species";
  exit(1);
  }

  if ( 2 != x2.size()) {
  cout << "ERROR: function computePropensitiesCellMilieu(x1, x2, propensities),"
          " array \"x2\" does not have "
          "the same size as the number of species";
  exit(1);
  }
  #endif //COMPUTEPROPENSITIES_ARRAY_SIZE_CHECK

  a(0) = 1.*x1(0);
  a(1) = 1.*x1(1);

  a(2) = (volume1/volume2)*1.*x2(0);
  a(3) = (volume1/volume2)*1.*x2(1);


  #ifdef COMPUTEPROPENSITIES_POSITIVITY_CHECK
  int i;
  for (i=0; i<4; i++)
  {
    if ( a(i) < 0.0 )
    {
      cout << "WARNING: function computePropensitiesCellMilieu(x, propensities),"
              " negative propensity detected: a(" << i << ") = " << a(i)
           << ". Resetting a to 0.0" << endl;
      a(i) = 0.0;
    }
  }
  #endif //COMPUTEPROPENSITIES_POSITIVITY_CHECK
}


void computePropensitiesTimeDependentDiffusion
(
  const Array<int,1>& x1,
  const Array<int,1>& x2,
  const double volume,
  const double volume0,
  const double volumeExt,
  Array<double,1>& a
)
{
  /**
  * Time dependent propensities of diffusion reaction processes.
  *
  * We assume the cell to be a cylinder whose volume is growing exponentially,
  * with its radius constant.
  * The volume at the end of the cell cycle, \f$ V(t_0+t_c) = 2V_0 \f$, where
  * \f$ V_0 \f$ is the initial volume at the beginning of the cell cycle and
  * \f$ t_0 \f$ the time of the previous division event. The time evolution
  * of the volume writes
  * \f[ V(t) = V_0 2^{\frac{t-t_0}{tc}} \f]
  *
  * We can calculate the evolution of the cell's surface area \f$ \Sigma(t) \f$,
  * \f[ \Sigma(t) = \Sigma_0 \theta(t) \f]
  * with
  * \f[ \theta(t) = 1 + \frac{L_0}{L_0+r_0}\left( 2^{\frac{t-t_0}{tc}}-1 \right) \f]
  * where \f$ L_0 \f$ is the initial length of the cylinder and \f$ r_0 \f$ the
  * radius of the cylinder.
  *
  * However, for the sake of simplicity, we use the following approximation:
  * \f[ \theta(t) = V(t)/V_1 \f]
  * where \f$V_1\f$ is the average of the cell volume during a cell cycle, i.e. for exponential growth is
  * \f$ V1 = V_0/ln(2) \f$
  *
  * The expression of the rates of reaction for the diffusion processes in and out
  * the cell, are given by
  * \f[ k_{out}(t) = D_x \theta(t) \f]
  * \f[ k_{in}(t) = D_x \theta(t) \frac{V(t)}{V_{ext}} \f]
  * where we included \f$ \Sigma_0 \f$ in the diffusion coefficient \f$ D_x \f$.
  *
  * @param x1 [in] %State vector of the cell.
  * @param x2 [in] %State vector of the milieu.
  * @param volume [in] Volume of the cell at time t.
  * @param volume0 [in] Volume of the cell at the beginning of the cell cycle (constant).
  * @param volumeExt [in] Volume of the external milieu at time t.
  * @param a [out] Array of the computed propensities.
  */


  #ifdef COMPUTEPROPENSITIES_ARRAY_SIZE_CHECK
  if ( 4 != a.size()) {
  cout << "ERROR: function computePropensitiesTimeDependentDiffusion(x1, x2, propensities),"
          " array \"propensities\" does not have "
          "the same size as the number of reactions";
  exit(1);
  }

  if ( 2 != x1.size()) {
  cout << "ERROR: function computePropensitiesTimeDependentDiffusion(x1, x2, propensities),"
          " array \"x1\" does not have "
          "the same size as the number of species";
  exit(1);
  }

  if ( 2 != x2.size()) {
  cout << "ERROR: function computePropensitiesTimeDependentDiffusion(x1, x2, propensities),"
          " array \"x2\" does not have "
          "the same size as the number of species";
  exit(1);
  }
  #endif //COMPUTEPROPENSITIES_ARRAY_SIZE_CHECK

  a(0) = (volume/(volume0/0.693147180559945))*1.*x1(0);
  a(1) = (volume/(volume0/0.693147180559945))*1.*x1(1);

  a(2) = (volume/volumeExt)*1.*x2(0);
  a(3) = (volume/volumeExt)*1.*x2(1);


  #ifdef COMPUTEPROPENSITIES_POSITIVITY_CHECK
  int i;
  for (i=0; i<4; i++)
  {
    if ( a(i) < 0.0 )
    {
      cout << "WARNING: function computePropensitiesTimeDependentDiffusion(x, propensities),"
              " negative propensity detected: a(" << i << ") = " << a(i)
           << ". Resetting a to 0.0" << endl;
      a(i) = 0.0;
    }
  }
  #endif //COMPUTEPROPENSITIES_POSITIVITY_CHECK
}


#endif // USE_CHEMICAL_LANGEVIN
