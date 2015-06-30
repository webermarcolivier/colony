/***************************************************************************//**
 * Project: Colony
 *
 * \file    mainwindow.cpp
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

#include "mainwindow.h"

#ifdef GUI

#include "ui_mainwindow.h"

//------------------------------------------------------------------------------

MainWindow::MainWindow(string inputFile, QWidget *parent)
    : QMainWindow(parent),  ui(new Ui::MainWindow)
{
  colonyScene_.init(this);

  ui->setupUi(this);
  ui->speciesColorComboBox->setToolTip("Species to display in cell color code");
  ui->speciesBackgroundColorComboBox->setToolTip("Species to display in background color code");
  ui->speciesColorMaxValueSpinBox->setToolTip("Maximum value of species concentration (nM) for cell color code");
  ui->speciesBackgroundColorMaxValueSpinBox->setToolTip("Maximum value of species concentration (nM) for background color code");
  ui->colorCodeGraphicsView->setToolTip("Cell color code gradient");
  ui->colorCodeBackgroundGraphicsView->setToolTip("Background color code gradient");
  ui->playButton->setToolTip("Play simulation");
  ui->pauseButton->setToolTip("Pause simulation");
  // Spanning the layout with the play/pause controls in the first two cells of the grid layout.
  ui->gridLayout_2->removeItem(ui->horizontalLayout);
  ui->horizontalLayout->setParent(0);
  ui->gridLayout_2->addLayout(ui->horizontalLayout, 0, 0, 1, 2);
  ui->gridLayout_2->setAlignment(ui->horizontalLayout, Qt::AlignLeft);

  graphicsCellParam_.init((QGraphicsScene*) &colonyScene_, colonyScene_.getCellGroup());

  simulator_ = new Simulator(inputFile, graphicsCellParam_);

  simulationEnd_ =  false;
  simulationPaused_ = true;

  speciesCellColorIndex_ = 0;
  speciesCellColorMaxValue_ = 100.0;
  speciesBackgroundColorIndex_ = 0;
  speciesBackgroundColorMaxValue_ = 100.0;
  gradColor1_.setHsv(0,0,0);
  gradColor2_.setHsv(120,255,255);
  gradColorBackground1_ = Qt::white;
  gradColorBackground2_ = Qt::blue;
  initiateGlobalParameterTable();
  cellStateIndex_ = -1;
  initiateCellStateTable(cellStateIndex_);
  initiateSpeciesColorComboBox();
  initiateSpeciesColorCodeView();
  initiateSpeciesColorBackgroundComboBox();
  initiateSpeciesColorBackgroundCodeView();

  connect( &colonyScene_, SIGNAL(cellGainedFocus(GraphicsCellQt*)), this, SLOT(cellGainedFocus(GraphicsCellQt*)));
  connect( &colonyScene_, SIGNAL(cellLostFocus(GraphicsCellQt*)), this, SLOT(cellLostFocus(GraphicsCellQt*)));

  initiateColonyView();
}

//------------------------------------------------------------------------------

MainWindow::~MainWindow()
{
  delete ui;
  delete simulator_;
}

//------------------------------------------------------------------------------

void MainWindow::runSimulation()
{
  while (!simulationEnd_ && !simulationPaused_)
  {
    /// - At the beginning of each trajectory, we initialize the colony view and integrate
    ///   the ODE dynamics for equilibration of the cells.
    if ( simulator_->getITimeSlice() == 0 )
    {
      initiateColonyView();
      cout << "ODE equilibration of colony... start" << endl;
      float refreshTime = 0.1;
      for (float time = 0.0; time < simulator_->getSpatialDynamicsEquilibrationTime();
           time+=refreshTime)
      {
        simulator_->setIsEquilibrationSteps(true);
        simulator_->computeSpatialIntegration(refreshTime);
        if ( int(time / refreshTime) % 20 == 0)
        {
          cout << "ODE equilibration of colony... " << time << "/"
               << simulator_->getSpatialDynamicsEquilibrationTime() << endl;
        }
        ui->colonyView->update();
        qApp->processEvents();
      }
      cout << "ODE equilibration of colony... end" << endl;
    }

    /// - Simulation step of the Simulator.
    int j;
    for (j=0;j<1;++j)
    {
      simulationEnd_ = simulator_->computeSimulationStep();
    }

    /// - Update MainWindow (all widgets update is done here).
    updateWindow();
  }

  if (simulationEnd_)
  {
    endSimulation();
  }
}

//------------------------------------------------------------------------------

void MainWindow::endSimulation()
{
}

//------------------------------------------------------------------------------

void MainWindow::playPauseSimulation()
{
  if (simulationPaused_ == true)
  {
    cout << "simulation resume " << endl;
    simulationPaused_ = false;
    runSimulation();
  }
  else
  {
    cout << "simulation pause " << endl;
    simulationPaused_ = true;
  }
}

//------------------------------------------------------------------------------

void MainWindow::playSimulation()
{
  if (simulationPaused_ == true)
  {
    cout << "simulation resume " << endl;
    simulationPaused_ = false;
    runSimulation();
  }
}

//------------------------------------------------------------------------------

void MainWindow::pauseSimulation()
{
  if (simulationPaused_ == false)
  {
    cout << "simulation pause " << endl;
    simulationPaused_ = true;
  }
}

//------------------------------------------------------------------------------

void MainWindow::fitInView()
{
  ui->colonyView->fitInView(colonyScene_.getCellGroup(),Qt::KeepAspectRatio);
  zoomValue_ = ui->colonyView->matrix().m11();
}

//------------------------------------------------------------------------------

void MainWindow::updateWindow()
{
  // update all information in window

  emit updateNCells( simulator_->cellCollection_.getNCells() );

  QString stringTimeSimCurrentAndTotal = QString("%1").arg(simulator_->getTime(), 6, 'f', 2) + '/' +
                                         QString("%1").arg(simulator_->getSimulationTime(), 6, 'f', 2);
  emit updateTimeSimulation( stringTimeSimCurrentAndTotal );

  QString stringTrajectoryCounter = QString("%1").arg(simulator_->getITrajectory()+1, 4, 10) + '/' +
                                    QString("%1").arg(simulator_->getNTrajectory(), 4, 10);
  emit updateTrajectoryCounter( stringTrajectoryCounter );

  QString stringParameterSetCounter = QString("%1").arg(simulator_->getIParameterSet()+1, 4, 10) + '/' +
                                      QString("%1").arg(simulator_->getNParameterSet(), 4, 10);
  emit updateParameterSetCounter( stringParameterSetCounter );

  emit updateSimulationProgressBar( calculateSimulationTotalProgress() );

  updateGlobalParameterTable();

  updateColonyView();

  updateCellStateTable();

  qApp->processEvents();

  repaint();
}

//------------------------------------------------------------------------------

int MainWindow::calculateSimulationTotalProgress()
{
  double wPar = 1.0 / double(simulator_->getNParameterSet());
  double iPar = double(simulator_->getIParameterSet());
  double wTraj = 1.0 / double(simulator_->getNTrajectory());
  double iTraj = double(simulator_->getITrajectory());
  double wTime = 1.0 / simulator_->getSimulationTime();
  double iTime = simulator_->getTime();

  //cout <<  wPar * ( iPar + wTraj*( iTraj + wTime*iTime ) ) << endl;
  return int(100 * wPar * ( iPar + wTraj*( iTraj + wTime*iTime ) ) );
}

//------------------------------------------------------------------------------

Array<QString,2> MainWindow::getGlobalParameterArray()
{
  // Declare data for the model: the 2 dimensional QString array.
  Array<QString,2> globalParameterStringArray;
  int nGlobalParameter = simulator_->getNGlobalParameter();
  globalParameterStringArray.resize( nGlobalParameter, 2 );

  // Get global parameter names and values from the simulator.
  vector<string> listGlobalParametersName = simulator_->getlistGlobalParametersName();
  Array<double,1> globalParameterArray = simulator_->getGlobalParameterArray();

  // Fill in the values of the QString array.
  int i;
  for(i=0;i<nGlobalParameter;++i)
  {
    globalParameterStringArray(i,0) = QString(listGlobalParametersName[i].c_str());
    globalParameterStringArray(i,1) = QString("%1").arg(globalParameterArray(i), 8, 'f', 3);
  }
  return globalParameterStringArray;
}

//------------------------------------------------------------------------------

void MainWindow::initiateGlobalParameterTable()
{

  Array<QString,2> globalParameterStringArray( getGlobalParameterArray().copy() );

  /// - Create the model with the global parameter array from the Simulator.
  globalParameterModel_ = new StringTableModel( globalParameterStringArray );

  /// - Set the columns headers.
  globalParameterModel_->setHeaderData(0, Qt::Horizontal, "Name");
  globalParameterModel_->setHeaderData(1, Qt::Horizontal, "Value");

  /// - Connect the TableView of the interface with the created model.
  ui->globalParameterTableView->setModel(globalParameterModel_);
  ui->globalParameterTableView->show();
}

//------------------------------------------------------------------------------

void MainWindow::updateGlobalParameterTable()
{
  Array<QString,2> globalParameterStringArray( getGlobalParameterArray().copy() );
  int nGlobalParameter = simulator_->getNGlobalParameter();

  int i;
  for(i=0;i<nGlobalParameter;++i)
  {
    globalParameterModel_->setData(globalParameterModel_->index(i,0), globalParameterStringArray(i,0));
    globalParameterModel_->setData(globalParameterModel_->index(i,1), globalParameterStringArray(i,1));
  }
}

//------------------------------------------------------------------------------

QColor MainWindow::gradientColor(double r)
{
  // Set the color gradient in HSV values.
  QColor color;
  int h1 = gradColor1_.hue();
  int s1 = gradColor1_.saturation();
  int v1 = gradColor1_.value();
  int h2 = gradColor2_.hue();
  int s2 = gradColor2_.saturation();
  int v2 = gradColor2_.value();

  double x = max(min(r,1.0),0.0);
  color.setHsv(int((1.-x)*h1+x*h2), int((1.-x)*s1+x*s2), int((1.-x)*v1+x*v2));

  return color;
}

//------------------------------------------------------------------------------

QColor MainWindow::gradientColorBackground(double r)
{
  // Set the color gradient in HSV values.
  QColor color;
  int h1 = gradColorBackground1_.hue();
  int s1 = gradColorBackground1_.saturation();
  int v1 = gradColorBackground1_.value();
  int h2 = gradColorBackground2_.hue();
  int s2 = gradColorBackground2_.saturation();
  int v2 = gradColorBackground2_.value();

  s1 = 0;
  s2 = 255;
  v1 = 255;
  v2 = 220;

  double x = max(min(r,1.0),0.0);
  color.setHsv(h2, int(x*x*s2), int((1.-x)*v1+x*v2));

  return color;
}

//------------------------------------------------------------------------------

void MainWindow::initiateColonyView()
{
  /// - Set #lastUpdateTime_ to 0.
  lastUpdateTime_ = 0.0;

  /// - Initiate the scene.
  colonyScene_.setFocusItem(0);
  // The scene rectangle is used only for indexing the graphical items. It should be
  // big enough to contain all the items.
  double sceneMaxSize = 400;
  colonyScene_.setSceneRect(-sceneMaxSize/2.0,-sceneMaxSize/2.0,sceneMaxSize,sceneMaxSize);
  colonyScene_.setItemIndexMethod(QGraphicsScene::NoIndex);
  /*
  /// - Add walls to the scene.
  double trapLength = 8.0; // in microns
  double trapWidth = 10.0; // in microns
  colonyScene_.addWall(-(trapLength/2.)*micronPixels_, 0.0*micronPixels_, 1.0, 0,
                       2.0*micronPixels_, trapWidth*micronPixels_);
  colonyScene_.addWall(+(trapLength/2.)*micronPixels_, 0.0*micronPixels_,-1.0, 0,
                       2.0*micronPixels_, trapWidth*micronPixels_);
  */
  updateColonyView();

  /// - Link the view to the scene and show it.
  ui->colonyView->setScene(&colonyScene_);
  ui->colonyView->show();
  QTimer::singleShot(100, this, SLOT(fitInView()));
}

//------------------------------------------------------------------------------

void MainWindow::updateColonyView()
{
  /// - Update cells length and color based on the volume and species concentration of simulator's cells.
  double time = simulator_->getTime();
  int nCells = simulator_->cellCollection_.getNCells();

  double r;
  for(int i=0;i<nCells;++i)
  {
    simulator_->cellCollection_[i].updateGraphicsCell(time);

    // Set the color value depending on the state of the cell.
    r = simulator_->cellCollection_[i].getXconc(time)(speciesCellColorIndex_) / speciesCellColorMaxValue_;
    simulator_->cellCollection_[i].setColor( gradientColor(r) );
  }

  /// - Update scene background color base on the species concentration of the simulator's milieu.
  r = simulator_->cellCollection_.getMilieu().getXconc(time)(speciesBackgroundColorIndex_) / speciesBackgroundColorMaxValue_;
  colonyScene_.setBackgroundBrush( gradientColorBackground(r) );

  /*
  #ifdef TIME_DEPENDENT_PROPENSITIES
    /// - Compute the ODE simulation step (cells spatial dynamics).
    double timeStep = simulator_->getTime() - lastUpdateTime_;
    colonyScene_.computeODEStep(timeStep, false);
  #endif
  */

  /// - Update #lastUpdateTime_.
  lastUpdateTime_ =  simulator_->getTime();
}

//------------------------------------------------------------------------------

void MainWindow::on_zoomSlider_valueChanged(int zoomValue)
{
  qreal zoomValueReal = pow(10.0,3.0*((zoomValue/50.0)-1.0)) * zoomValue_;
  QTransform transform;
  transform.scale(zoomValueReal, zoomValueReal);
  ui->colonyView->setTransform(transform);
}

//------------------------------------------------------------------------------

void MainWindow::initiateSpeciesColorComboBox()
{
  QStringList speciesNameList;
  vector<string> vectorList = simulator_->cellCollection_.getListSpeciesName();
  vector<string>::iterator it;
  for(it=vectorList.begin();it!=vectorList.end(); ++it)
  {
    speciesNameList << QString((*it).c_str());
  }
  speciesNameListModel_.setStringList(speciesNameList);

  ui->speciesColorComboBox->setModel(&speciesNameListModel_);
  ui->speciesColorComboBox->setCurrentIndex(speciesCellColorIndex_);
}

//------------------------------------------------------------------------------

void MainWindow::initiateSpeciesColorCodeView()
{
  int sceneWidth = 196;
  int sceneHeight = 26;
  speciesColorCodeScene_.items().clear();
  speciesColorCodeScene_.setSceneRect(0,0,sceneWidth,sceneHeight);

  int j, nGradientPoints = 50;
  for(j=0;j<nGradientPoints; ++j)
  {
    QGraphicsRectItem *rect = new QGraphicsRectItem( double(j)*sceneWidth/double(nGradientPoints), 0.0, sceneWidth/double(nGradientPoints), sceneHeight);
    QBrush brush( gradientColor((double(j)+0.5)/double(nGradientPoints)) );
    rect->setBrush(brush);
    QColor transparentColor(0,0,0,0);
    rect->setPen(QPen(transparentColor));
    speciesColorCodeScene_.addItem(rect);
  }
  ui->colorCodeGraphicsView->setScene(&speciesColorCodeScene_);
  ui->colorCodeGraphicsView->show();
}

//------------------------------------------------------------------------------

// This method is not working
void MainWindow::updateSpeciesColorCodeView()
{
  // Set the scene size to the same as the view's size.
  int sceneWidth = speciesColorCodeScene_.width();
  int viewWidth = ui->colorCodeGraphicsView->sizeHint().width();
cout << "sceneWidth " << sceneWidth << endl;
cout << "viewWidth " << viewWidth << endl;
  int sceneHeight = speciesColorCodeScene_.height();
  int viewHeight = ui->colorCodeGraphicsView->height();

  QTransform matrix;
  matrix.scale(viewWidth/sceneWidth,viewHeight/sceneHeight);
  ui->colorCodeGraphicsView->setTransform(matrix, false);
}

//------------------------------------------------------------------------------

// This method is not working
void MainWindow::initiateSpeciesColorBackgroundComboBox()
{
  QStringList speciesNameList;
  vector<string> vectorList = simulator_->cellCollection_.getMilieu().getListSpeciesName();
  vector<string>::iterator it;
  for(it=vectorList.begin();it!=vectorList.end(); ++it)
  {
    speciesNameList << QString((*it).c_str());
  }
  speciesMilieuNameListModel_.setStringList(speciesNameList);

  ui->speciesBackgroundColorComboBox->setModel(&speciesMilieuNameListModel_);
  ui->speciesBackgroundColorComboBox->setCurrentIndex(speciesBackgroundColorIndex_);
}

//------------------------------------------------------------------------------

void MainWindow::initiateSpeciesColorBackgroundCodeView()
{
  int sceneWidth = 196;
  int sceneHeight = 26;
  speciesColorCodeBackgroundScene_.items().clear();
  speciesColorCodeBackgroundScene_.setSceneRect(0,0,sceneWidth,sceneHeight);

  int j, nGradientPoints = 50;
  for(j=0;j<nGradientPoints; ++j)
  {
    QGraphicsRectItem *rect = new QGraphicsRectItem( double(j)*sceneWidth/double(nGradientPoints), 0.0, sceneWidth/double(nGradientPoints), sceneHeight);
    QBrush brush( gradientColorBackground((double(j)+0.5)/double(nGradientPoints)) );
    rect->setBrush(brush);
    QColor transparentColor(0,0,0,0);
    rect->setPen(QPen(transparentColor));
    speciesColorCodeBackgroundScene_.addItem(rect);
  }

  ui->colorCodeBackgroundGraphicsView->setScene(&speciesColorCodeBackgroundScene_);
  ui->colorCodeBackgroundGraphicsView->show();
}

//------------------------------------------------------------------------------

void MainWindow::updateSpeciesColorBackgroundCodeView()
{
  // Set the scene size to the same as the view's size.
  int sceneWidth = speciesColorCodeBackgroundScene_.width();
  int viewWidth = ui->colorCodeBackgroundGraphicsView->width();
cout << "viewWidth " << viewWidth << endl;
cout << "sceneWidth " << sceneWidth << endl;
  int sceneHeight = speciesColorCodeBackgroundScene_.height();
  int viewHeight = ui->colorCodeBackgroundGraphicsView->height();

  QTransform matrix;
  matrix.scale(viewWidth/sceneWidth,viewHeight/sceneHeight);
  ui->colorCodeBackgroundGraphicsView->setTransform(matrix);
}

//------------------------------------------------------------------------------

void MainWindow::setSpeciesCellColorIndex(int index)
{
  speciesCellColorIndex_ = index;
  updateColonyView();
}

//------------------------------------------------------------------------------

void MainWindow::setSpeciesCellColorMaxValue(double max)
{
  speciesCellColorMaxValue_ = max;
  updateColonyView();
}

//------------------------------------------------------------------------------

void MainWindow::setSpeciesBackgroundColorIndex(int index)
{
  speciesBackgroundColorIndex_ = index;
  updateColonyView();
}

//------------------------------------------------------------------------------

void MainWindow::setSpeciesBackgroundColorMaxValue(double max)
{
  speciesBackgroundColorMaxValue_ = max;
  updateColonyView();
}

//------------------------------------------------------------------------------

void MainWindow::initiateCellStateTable(int cellIndex)
{
  cellStateIndex_ = cellIndex;

  Array<QString,2> cellStateStringArray( getCellStateArray(cellIndex).copy() );

  /// - Create the model with the cell state array from the Simulator.
  cellStateModel_ = new StringTableModel( cellStateStringArray );

  /// - Set the columns headers.
  cellStateModel_->setHeaderData(0, Qt::Horizontal, "Species");
  cellStateModel_->setHeaderData(1, Qt::Horizontal, "Nb molecules");
  cellStateModel_->setHeaderData(2, Qt::Horizontal, "Concentration (nM)");

  /// - Connect the TableView of the interface with the created model.
  ui->cellStateTableView->setModel (cellStateModel_ );
  ui->cellStateTableView->setColumnWidth(0,80);
  ui->cellStateTableView->setColumnWidth(1,100);
  ui->cellStateTableView->setColumnWidth(2,100);
  ui->cellStateTableView->show();
}

//------------------------------------------------------------------------------

void MainWindow::updateCellStateTable()
{
  if (cellStateModel_ != 0)
  { 
    Array<QString,2> cellStateStringArray( getCellStateArray(cellStateIndex_).copy() );
    int n = cellStateStringArray.extent(firstDim);

    int i;
    for(i=0;i<n;++i)
    {
      cellStateModel_->setData(cellStateModel_->index(i,0), cellStateStringArray(i,0));
      cellStateModel_->setData(cellStateModel_->index(i,1), cellStateStringArray(i,1));
      cellStateModel_->setData(cellStateModel_->index(i,2), cellStateStringArray(i,2));
    }
  }
}

//------------------------------------------------------------------------------

Array<QString,2> MainWindow::getCellStateArray(int cellIndex)
{
  /// - Declare data for the model: a QString table.
  Array<QString,2> cellStateStringArray;
  int nSpecies;
  if (cellIndex >= 0)
  {
    nSpecies = simulator_->cellCollection_[cellIndex].getNSpecies();
  } else {
    nSpecies = simulator_->cellCollection_.getMilieu().getNSpecies();
  }
  cellStateStringArray.resize( nSpecies + 1, 3 );

  /// - Get species names, number of molecules and concentrations from the simulator.
  vector<string> listSpeciesName;
  Array<int,1> cellStateValueArray;
  Array<double,1> cellStateConcArray;
  if (cellIndex >= 0)
  {
    listSpeciesName = simulator_->cellCollection_[cellIndex].getListSpeciesName();
    cellStateValueArray.resize( simulator_->cellCollection_[cellIndex].getX().size()  );
    cellStateValueArray = simulator_->cellCollection_[cellIndex].getX();
    cellStateConcArray.resize( simulator_->cellCollection_[cellIndex].getXconc( simulator_->getTime() ).size() );
    cellStateConcArray = simulator_->cellCollection_[cellIndex].getXconc( simulator_->getTime() );
  } else {
    listSpeciesName = simulator_->cellCollection_.getMilieu().getListSpeciesName();
    cellStateValueArray.resize( simulator_->cellCollection_.getMilieu().getListSpeciesName().size()  );
    cellStateValueArray = simulator_->cellCollection_.getMilieu().getX();
    cellStateConcArray.resize( simulator_->cellCollection_.getMilieu().getXconc( simulator_->getTime() ).size() );
    cellStateConcArray = simulator_->cellCollection_.getMilieu().getXconc( simulator_->getTime() );
  }

  /// - The first line of the table contains the cell number.
  if (cellIndex >= 0)
  {
    cellStateStringArray(0,0) = QString("cell #");
    cellStateStringArray(0,1) = QString("%1").arg( cellIndex, 8 );
    cellStateStringArray(0,2) = QString("");
  } else {
    cellStateStringArray(0,0) = QString("Milieu");
    cellStateStringArray(0,1) = QString("");
    cellStateStringArray(0,2) = QString("");
  }

  /// - Fill in the values of the QString array with the data from the simulator.
  int i;
  for(i=0;i<nSpecies;++i)
  {
    cellStateStringArray(i+1,0) = QString(listSpeciesName[i].c_str());
    cellStateStringArray(i+1,1) = QString("%1").arg( cellStateValueArray(i), 8 );
                                //arg(double(cellStateValueArray(i)), 8, 'f', 3);
    cellStateStringArray(i+1,2) = QString("%1").arg(cellStateConcArray(i), 8, 'f', 3);
  }

  return cellStateStringArray;
}

//------------------------------------------------------------------------------

void MainWindow::cellLostFocus(GraphicsCellQt *cell)
{
  delete cellStateModel_;
  // Cell index -1 reads for the milieu.
  cellStateIndex_ = -1;
  int cellIndex = -1;
  initiateCellStateTable(cellIndex);
}

//------------------------------------------------------------------------------

void MainWindow::cellGainedFocus(GraphicsCellQt *cell)
{
  delete cellStateModel_;
  int cellIndex = cell->getCellIndex();
  initiateCellStateTable(cellIndex);
}

//------------------------------------------------------------------------------

int MainWindow::getSimulatorITimeSlice() const
{
  return simulator_->getITimeSlice();
}

//------------------------------------------------------------------------------

double MainWindow::getCellLength0() const
{
  return cellLength0_;
}

//------------------------------------------------------------------------------

double MainWindow::getCellHeight0() const
{
  return cellHeight0_;
}

//------------------------------------------------------------------------------

#endif // GUI
