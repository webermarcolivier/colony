# -------------------------------------------------
# Project created by QtCreator 2010-03-22T18:58:59
# -------------------------------------------------
TARGET = COLONY
TEMPLATE = app
SOURCES += main.cpp \
    mainwindow.cpp \
    TimeSliceCellsAvg.cpp \
    TimeSlice.cpp \
    State.cpp \
    Simulator.cpp \
    RandomNumberGenerator.cpp \
    Output.cpp \
    Milieu.cpp \
    iotools.cpp \
    IntegratorGillespieModified.cpp \
    IntegratorGillespie.cpp \
    IntegratorContext.cpp \
    Input.cpp \
    GlobalArrayInterface.cpp \
    GlobalArray.cpp \
    ChemicalSystem.cpp \
    CellMilieuChemicalSystem.cpp \
    CellLineageGeneration.cpp \
    CellCollection.cpp \
    CellBase.cpp \
    Cell.cpp \
    computePropensitiesFunctions.cpp \
    StringTableModel.cpp \
    GraphicsCellScene.cpp \
    GraphicsCellGroup.cpp \
    GraphicsWall.cpp \
    GraphicsCellBase.cpp \
    GraphicsCellQt.cpp \
    GraphicsCellComposite.cpp \
    SpatialIntegratorContext.cpp \
    SpatialIntegratorODE.cpp \
    GraphicsCellODE.cpp \
    IntegratorChemicalLangevin.cpp \
    computeChemicalLangevinIncrement.cpp \
    chemicalLangevinComputeIncrement.cpp
HEADERS += mainwindow.h \
    Param/SimulatorParam.h \
    Param/OutputParam.h \
    Param/CellCollectionParam.h \
    TimeSliceCellsAvg.h \
    TimeSlice.h \
    State.h \
    Simulator.h \
    RandomNumberGenerator.h \
    Output.h \
    Milieu.h \
    iotools.h \
    IntegratorGillespieModified.h \
    IntegratorGillespie.h \
    IntegratorContext.h \
    Input.h \
    GlobalArrayInterface.h \
    GlobalArray.h \
    doxygen_documentation.h \
    debug.h \
    ChemicalSystem.h \
    CellMilieuChemicalSystem.h \
    CellLineageGeneration.h \
    CellCollection.h \
    CellBase.h \
    Cell.h \
    computePropensitiesFunctions.h \
    StringTableModel.h \
    GraphicsCellScene.h \
    GraphicsCellGroup.h \
    GraphicsWall.h \
    GraphicsCellBase.h \
    GraphicsCellQt.h \
    GraphicsCellComposite.h \
    Param/GraphicsCellCompositeParam.h \
    SpatialIntegratorContext.h \
    SpatialIntegratorODE.h \
    GraphicsCellODE.h \
    Param/CellInitParam.h \
    Param/CellBaseInitParam.h \
    Param/StateInitParam.h \
    Param/ChemicalSystemInitParam.h \
    Param/CellMilieuChemicalSystemInitParam.h \
    IntegratorChemicalLangevin.h
FORMS += mainwindow.ui
OTHER_FILES += icons/PlayGreenButton.png
RESOURCES += colonyResourceFile.qrc
QT += xml widgets
LIBS += -L/usr/lib -lblitz -lboost_system -lboost_filesystem -lgsl -lode -lopenblas -lstdc++ #-L/usr/local/atlas/lib -lgfortran -latlas -lf77blas -llapack #-lGLU -lglut
#INCLUDEPATH +=
#QMAKE_CFLAGS+=-pg
#QMAKE_CXXFLAGS+=-pg
#QMAKE_LFLAGS+=-pg
QMAKE_CXXFLAGS+=-static-libstdc++ -static-libgcc -static
QMAKE_CFLAGS+=-static-libgcc -static
DISTFILES +=
