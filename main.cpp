/***************************************************************************//**
 * Project: Colony
 *
 * \file    main.cpp
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

#include "compilation_options.h"

#ifdef GUI
  #include <QtGui>
  #include "mainwindow.h"
#else
  #include "Simulator.h"
#endif // GUI

#include <iostream>
#include <string>     // manipulate strings of characters


int main(int argc, char *argv[])
{

  // The program requires one command line argument, the name of the input file, which should be in the current directory
  if (argc != 2)
  {
    cout << "ERROR: main program, number of command line arguments should be 1 (name of the input file in the current directory)." << endl;
    exit(0);
  }
  string inputFile(argv[1]);

  #ifdef GUI

    // GUI version, the simulator class is created by the mainWindow.

    QApplication a(argc, argv);
    MainWindow w(inputFile);
    w.show();
    return a.exec();

  #else

    // No GUI version, we create here the simulator object and launch the simulation.

    Simulator *simulator_ = new Simulator(inputFile);

    bool simulationEnd_ = false;

    while (!simulationEnd_)
    {
      /// - Simulation step of the Simulator.
      simulationEnd_ = simulator_->computeSimulationStep();
    }

  #endif //GUI
}
