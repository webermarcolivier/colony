/***************************************************************************//**
 * Project: Colony
 *
 * \file    CellLineageGeneration.cpp
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

#include "CellLineageGeneration.h"

//------------------------------------------------------------------------------

CellLineageGeneration::CellLineageGeneration()
{}

//------------------------------------------------------------------------------

CellLineageGeneration::~CellLineageGeneration()
{}

//------------------------------------------------------------------------------

CellLineageGeneration::CellLineageGeneration(const CellLineageGeneration& t)
  : time_(t.time_),
    nCells_(t.nCells_),
    nCellsDelta_(t.nCellsDelta_),
    motherCellsIndices_(t.motherCellsIndices_.copy())
{}

//------------------------------------------------------------------------------

void CellLineageGeneration::initialize
(
  double time,
  int nCells,
  int nCellsDelta,
  Array<int,1> motherCellsIndices
)
{
  time_ = time;
  nCells_ = nCells;
  nCellsDelta_ = nCellsDelta;
  motherCellsIndices_.resize(motherCellsIndices.size());
  motherCellsIndices_ = motherCellsIndices;
}

//------------------------------------------------------------------------------
