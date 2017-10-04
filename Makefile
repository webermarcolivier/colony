#############################################################################
# Makefile for building: COLONY
# Generated by qmake (2.01a) (Qt 4.8.7) on: lun oct 2 15:31:19 2017
# Project:  COLONY.pro
# Template: app
# Command: /usr/lib/x86_64-linux-gnu/qt4/bin/qmake -o Makefile COLONY.pro
#############################################################################

####### Compiler, tools and options

CC            = gcc
CXX           = g++
DEFINES       = -DQT_NO_DEBUG -DQT_GUI_LIB -DQT_CORE_LIB -DQT_SHARED
CFLAGS        = -m64 -pipe -O2 -Wall -W -D_REENTRANT $(DEFINES)
CXXFLAGS      = -m64 -pipe -O2 -Wall -W -D_REENTRANT $(DEFINES)
INCPATH       = -I/usr/share/qt4/mkspecs/linux-g++-64 -I. -I/usr/include/qt4/QtCore -I/usr/include/qt4/QtGui -I/usr/include/qt4 -I. -I.
LINK          = g++
LFLAGS        = -m64 -Wl,-O1
LIBS          = $(SUBLIBS)  -L/usr/lib/x86_64-linux-gnu -L/usr/lib -lblitz -lboost_system -lboost_filesystem -lgsl -lode -lopenblas -lstdc++ -lQtGui -lQtCore -lpthread 
AR            = ar cqs
RANLIB        = 
QMAKE         = /usr/lib/x86_64-linux-gnu/qt4/bin/qmake
TAR           = tar -cf
COMPRESS      = gzip -9f
COPY          = cp -f
SED           = sed
COPY_FILE     = $(COPY)
COPY_DIR      = $(COPY) -r
STRIP         = strip
INSTALL_FILE  = install -m 644 -p
INSTALL_DIR   = $(COPY_DIR)
INSTALL_PROGRAM = install -m 755 -p
DEL_FILE      = rm -f
SYMLINK       = ln -f -s
DEL_DIR       = rmdir
MOVE          = mv -f
CHK_DIR_EXISTS= test -d
MKDIR         = mkdir -p

####### Output directory

OBJECTS_DIR   = ./

####### Files

SOURCES       = main.cpp \
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
		chemicalLangevinComputeIncrement.cpp moc_mainwindow.cpp \
		moc_StringTableModel.cpp \
		moc_GraphicsCellScene.cpp \
		moc_GraphicsWall.cpp \
		moc_GraphicsCellQt.cpp \
		qrc_colonyResourceFile.cpp
OBJECTS       = main.o \
		mainwindow.o \
		TimeSliceCellsAvg.o \
		TimeSlice.o \
		State.o \
		Simulator.o \
		RandomNumberGenerator.o \
		Output.o \
		Milieu.o \
		iotools.o \
		IntegratorGillespieModified.o \
		IntegratorGillespie.o \
		IntegratorContext.o \
		Input.o \
		GlobalArrayInterface.o \
		GlobalArray.o \
		ChemicalSystem.o \
		CellMilieuChemicalSystem.o \
		CellLineageGeneration.o \
		CellCollection.o \
		CellBase.o \
		Cell.o \
		computePropensitiesFunctions.o \
		StringTableModel.o \
		GraphicsCellScene.o \
		GraphicsCellGroup.o \
		GraphicsWall.o \
		GraphicsCellBase.o \
		GraphicsCellQt.o \
		GraphicsCellComposite.o \
		SpatialIntegratorContext.o \
		SpatialIntegratorODE.o \
		GraphicsCellODE.o \
		IntegratorChemicalLangevin.o \
		chemicalLangevinComputeIncrement.o \
		moc_mainwindow.o \
		moc_StringTableModel.o \
		moc_GraphicsCellScene.o \
		moc_GraphicsWall.o \
		moc_GraphicsCellQt.o \
		qrc_colonyResourceFile.o
DIST          = /usr/share/qt4/mkspecs/common/unix.conf \
		/usr/share/qt4/mkspecs/common/linux.conf \
		/usr/share/qt4/mkspecs/common/gcc-base.conf \
		/usr/share/qt4/mkspecs/common/gcc-base-unix.conf \
		/usr/share/qt4/mkspecs/common/g++-base.conf \
		/usr/share/qt4/mkspecs/common/g++-unix.conf \
		/usr/share/qt4/mkspecs/qconfig.pri \
		/usr/share/qt4/mkspecs/features/qt_functions.prf \
		/usr/share/qt4/mkspecs/features/qt_config.prf \
		/usr/share/qt4/mkspecs/features/exclusive_builds.prf \
		/usr/share/qt4/mkspecs/features/default_pre.prf \
		/usr/share/qt4/mkspecs/features/release.prf \
		/usr/share/qt4/mkspecs/features/default_post.prf \
		/usr/share/qt4/mkspecs/features/shared.prf \
		/usr/share/qt4/mkspecs/features/unix/gdb_dwarf_index.prf \
		/usr/share/qt4/mkspecs/features/warn_on.prf \
		/usr/share/qt4/mkspecs/features/qt.prf \
		/usr/share/qt4/mkspecs/features/unix/thread.prf \
		/usr/share/qt4/mkspecs/features/moc.prf \
		/usr/share/qt4/mkspecs/features/resources.prf \
		/usr/share/qt4/mkspecs/features/uic.prf \
		/usr/share/qt4/mkspecs/features/yacc.prf \
		/usr/share/qt4/mkspecs/features/lex.prf \
		/usr/share/qt4/mkspecs/features/include_source_dir.prf \
		COLONY.pro
QMAKE_TARGET  = COLONY
DESTDIR       = 
TARGET        = COLONY

first: all
####### Implicit rules

.SUFFIXES: .o .c .cpp .cc .cxx .C

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.cc.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.cxx.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.C.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.c.o:
	$(CC) -c $(CFLAGS) $(INCPATH) -o "$@" "$<"

####### Build rules

all: Makefile $(TARGET)

$(TARGET): ui_mainwindow.h $(OBJECTS)  
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS) $(OBJCOMP) $(LIBS)

Makefile: COLONY.pro  /usr/share/qt4/mkspecs/linux-g++-64/qmake.conf /usr/share/qt4/mkspecs/common/unix.conf \
		/usr/share/qt4/mkspecs/common/linux.conf \
		/usr/share/qt4/mkspecs/common/gcc-base.conf \
		/usr/share/qt4/mkspecs/common/gcc-base-unix.conf \
		/usr/share/qt4/mkspecs/common/g++-base.conf \
		/usr/share/qt4/mkspecs/common/g++-unix.conf \
		/usr/share/qt4/mkspecs/qconfig.pri \
		/usr/share/qt4/mkspecs/features/qt_functions.prf \
		/usr/share/qt4/mkspecs/features/qt_config.prf \
		/usr/share/qt4/mkspecs/features/exclusive_builds.prf \
		/usr/share/qt4/mkspecs/features/default_pre.prf \
		/usr/share/qt4/mkspecs/features/release.prf \
		/usr/share/qt4/mkspecs/features/default_post.prf \
		/usr/share/qt4/mkspecs/features/shared.prf \
		/usr/share/qt4/mkspecs/features/unix/gdb_dwarf_index.prf \
		/usr/share/qt4/mkspecs/features/warn_on.prf \
		/usr/share/qt4/mkspecs/features/qt.prf \
		/usr/share/qt4/mkspecs/features/unix/thread.prf \
		/usr/share/qt4/mkspecs/features/moc.prf \
		/usr/share/qt4/mkspecs/features/resources.prf \
		/usr/share/qt4/mkspecs/features/uic.prf \
		/usr/share/qt4/mkspecs/features/yacc.prf \
		/usr/share/qt4/mkspecs/features/lex.prf \
		/usr/share/qt4/mkspecs/features/include_source_dir.prf \
		/usr/lib/x86_64-linux-gnu/libQtGui.prl \
		/usr/lib/x86_64-linux-gnu/libQtCore.prl
	$(QMAKE) -o Makefile COLONY.pro
/usr/share/qt4/mkspecs/common/unix.conf:
/usr/share/qt4/mkspecs/common/linux.conf:
/usr/share/qt4/mkspecs/common/gcc-base.conf:
/usr/share/qt4/mkspecs/common/gcc-base-unix.conf:
/usr/share/qt4/mkspecs/common/g++-base.conf:
/usr/share/qt4/mkspecs/common/g++-unix.conf:
/usr/share/qt4/mkspecs/qconfig.pri:
/usr/share/qt4/mkspecs/features/qt_functions.prf:
/usr/share/qt4/mkspecs/features/qt_config.prf:
/usr/share/qt4/mkspecs/features/exclusive_builds.prf:
/usr/share/qt4/mkspecs/features/default_pre.prf:
/usr/share/qt4/mkspecs/features/release.prf:
/usr/share/qt4/mkspecs/features/default_post.prf:
/usr/share/qt4/mkspecs/features/shared.prf:
/usr/share/qt4/mkspecs/features/unix/gdb_dwarf_index.prf:
/usr/share/qt4/mkspecs/features/warn_on.prf:
/usr/share/qt4/mkspecs/features/qt.prf:
/usr/share/qt4/mkspecs/features/unix/thread.prf:
/usr/share/qt4/mkspecs/features/moc.prf:
/usr/share/qt4/mkspecs/features/resources.prf:
/usr/share/qt4/mkspecs/features/uic.prf:
/usr/share/qt4/mkspecs/features/yacc.prf:
/usr/share/qt4/mkspecs/features/lex.prf:
/usr/share/qt4/mkspecs/features/include_source_dir.prf:
/usr/lib/x86_64-linux-gnu/libQtGui.prl:
/usr/lib/x86_64-linux-gnu/libQtCore.prl:
qmake:  FORCE
	@$(QMAKE) -o Makefile COLONY.pro

dist: 
	@$(CHK_DIR_EXISTS) .tmp/COLONY1.0.0 || $(MKDIR) .tmp/COLONY1.0.0 
	$(COPY_FILE) --parents $(SOURCES) $(DIST) .tmp/COLONY1.0.0/ && $(COPY_FILE) --parents mainwindow.h Param/SimulatorParam.h Param/OutputParam.h Param/CellCollectionParam.h TimeSliceCellsAvg.h TimeSlice.h State.h Simulator.h RandomNumberGenerator.h Output.h Milieu.h iotools.h IntegratorGillespieModified.h IntegratorGillespie.h IntegratorContext.h Input.h GlobalArrayInterface.h GlobalArray.h doxygen_documentation.h ChemicalSystem.h CellMilieuChemicalSystem.h CellLineageGeneration.h CellCollection.h CellBase.h Cell.h computePropensitiesFunctions.h StringTableModel.h GraphicsCellScene.h GraphicsCellGroup.h GraphicsWall.h GraphicsCellBase.h GraphicsCellQt.h GraphicsCellComposite.h Param/GraphicsCellCompositeParam.h SpatialIntegratorContext.h SpatialIntegratorODE.h GraphicsCellODE.h Param/CellInitParam.h Param/CellBaseInitParam.h Param/StateInitParam.h Param/ChemicalSystemInitParam.h Param/CellMilieuChemicalSystemInitParam.h IntegratorChemicalLangevin.h compilation_options.h .tmp/COLONY1.0.0/ && $(COPY_FILE) --parents colonyResourceFile.qrc .tmp/COLONY1.0.0/ && $(COPY_FILE) --parents main.cpp mainwindow.cpp TimeSliceCellsAvg.cpp TimeSlice.cpp State.cpp Simulator.cpp RandomNumberGenerator.cpp Output.cpp Milieu.cpp iotools.cpp IntegratorGillespieModified.cpp IntegratorGillespie.cpp IntegratorContext.cpp Input.cpp GlobalArrayInterface.cpp GlobalArray.cpp ChemicalSystem.cpp CellMilieuChemicalSystem.cpp CellLineageGeneration.cpp CellCollection.cpp CellBase.cpp Cell.cpp computePropensitiesFunctions.cpp StringTableModel.cpp GraphicsCellScene.cpp GraphicsCellGroup.cpp GraphicsWall.cpp GraphicsCellBase.cpp GraphicsCellQt.cpp GraphicsCellComposite.cpp SpatialIntegratorContext.cpp SpatialIntegratorODE.cpp GraphicsCellODE.cpp IntegratorChemicalLangevin.cpp chemicalLangevinComputeIncrement.cpp .tmp/COLONY1.0.0/ && $(COPY_FILE) --parents mainwindow.ui .tmp/COLONY1.0.0/ && (cd `dirname .tmp/COLONY1.0.0` && $(TAR) COLONY1.0.0.tar COLONY1.0.0 && $(COMPRESS) COLONY1.0.0.tar) && $(MOVE) `dirname .tmp/COLONY1.0.0`/COLONY1.0.0.tar.gz . && $(DEL_FILE) -r .tmp/COLONY1.0.0


clean:compiler_clean 
	-$(DEL_FILE) $(OBJECTS)
	-$(DEL_FILE) *~ core *.core


####### Sub-libraries

distclean: clean
	-$(DEL_FILE) $(TARGET) 
	-$(DEL_FILE) Makefile


check: first

mocclean: compiler_moc_header_clean compiler_moc_source_clean

mocables: compiler_moc_header_make_all compiler_moc_source_make_all

compiler_moc_header_make_all: moc_mainwindow.cpp moc_StringTableModel.cpp moc_GraphicsCellScene.cpp moc_GraphicsWall.cpp moc_GraphicsCellQt.cpp
compiler_moc_header_clean:
	-$(DEL_FILE) moc_mainwindow.cpp moc_StringTableModel.cpp moc_GraphicsCellScene.cpp moc_GraphicsWall.cpp moc_GraphicsCellQt.cpp
moc_mainwindow.cpp: compilation_options.h \
		Simulator.h \
		Output.h \
		Input.h \
		iotools.h \
		RandomNumberGenerator.h \
		Param/CellCollectionParam.h \
		Param/CellBaseInitParam.h \
		Param/StateInitParam.h \
		Param/ChemicalSystemInitParam.h \
		Param/CellInitParam.h \
		Param/GraphicsCellCompositeParam.h \
		GraphicsCellScene.h \
		GraphicsWall.h \
		GraphicsCellGroup.h \
		GraphicsCellQt.h \
		GraphicsCellBase.h \
		Param/CellMilieuChemicalSystemInitParam.h \
		Param/OutputParam.h \
		Param/SimulatorParam.h \
		computePropensitiesFunctions.h \
		TimeSlice.h \
		TimeSliceCellsAvg.h \
		CellLineageGeneration.h \
		CellCollection.h \
		Milieu.h \
		CellBase.h \
		State.h \
		ChemicalSystem.h \
		Cell.h \
		GraphicsCellComposite.h \
		GraphicsCellODE.h \
		CellMilieuChemicalSystem.h \
		GlobalArrayInterface.h \
		GlobalArray.h \
		IntegratorContext.h \
		IntegratorGillespie.h \
		IntegratorGillespieModified.h \
		IntegratorChemicalLangevin.h \
		chemicalLangevinComputeIncrement.h \
		SpatialIntegratorContext.h \
		SpatialIntegratorODE.h \
		StringTableModel.h \
		mainwindow.h
	/usr/lib/x86_64-linux-gnu/qt4/bin/moc $(DEFINES) $(INCPATH) mainwindow.h -o moc_mainwindow.cpp

moc_StringTableModel.cpp: compilation_options.h \
		Simulator.h \
		Output.h \
		Input.h \
		iotools.h \
		RandomNumberGenerator.h \
		Param/CellCollectionParam.h \
		Param/CellBaseInitParam.h \
		Param/StateInitParam.h \
		Param/ChemicalSystemInitParam.h \
		Param/CellInitParam.h \
		Param/GraphicsCellCompositeParam.h \
		GraphicsCellScene.h \
		GraphicsWall.h \
		GraphicsCellGroup.h \
		GraphicsCellQt.h \
		GraphicsCellBase.h \
		Param/CellMilieuChemicalSystemInitParam.h \
		Param/OutputParam.h \
		Param/SimulatorParam.h \
		computePropensitiesFunctions.h \
		TimeSlice.h \
		TimeSliceCellsAvg.h \
		CellLineageGeneration.h \
		CellCollection.h \
		Milieu.h \
		CellBase.h \
		State.h \
		ChemicalSystem.h \
		Cell.h \
		GraphicsCellComposite.h \
		GraphicsCellODE.h \
		CellMilieuChemicalSystem.h \
		GlobalArrayInterface.h \
		GlobalArray.h \
		IntegratorContext.h \
		IntegratorGillespie.h \
		IntegratorGillespieModified.h \
		IntegratorChemicalLangevin.h \
		chemicalLangevinComputeIncrement.h \
		SpatialIntegratorContext.h \
		SpatialIntegratorODE.h \
		StringTableModel.h
	/usr/lib/x86_64-linux-gnu/qt4/bin/moc $(DEFINES) $(INCPATH) StringTableModel.h -o moc_StringTableModel.cpp

moc_GraphicsCellScene.cpp: compilation_options.h \
		GraphicsWall.h \
		GraphicsCellGroup.h \
		GraphicsCellQt.h \
		GraphicsCellBase.h \
		GraphicsCellScene.h \
		GraphicsCellScene.h
	/usr/lib/x86_64-linux-gnu/qt4/bin/moc $(DEFINES) $(INCPATH) GraphicsCellScene.h -o moc_GraphicsCellScene.cpp

moc_GraphicsWall.cpp: GraphicsWall.h
	/usr/lib/x86_64-linux-gnu/qt4/bin/moc $(DEFINES) $(INCPATH) GraphicsWall.h -o moc_GraphicsWall.cpp

moc_GraphicsCellQt.cpp: compilation_options.h \
		GraphicsCellBase.h \
		GraphicsCellScene.h \
		GraphicsWall.h \
		GraphicsCellGroup.h \
		GraphicsCellQt.h \
		GraphicsCellQt.h
	/usr/lib/x86_64-linux-gnu/qt4/bin/moc $(DEFINES) $(INCPATH) GraphicsCellQt.h -o moc_GraphicsCellQt.cpp

compiler_rcc_make_all: qrc_colonyResourceFile.cpp
compiler_rcc_clean:
	-$(DEL_FILE) qrc_colonyResourceFile.cpp
qrc_colonyResourceFile.cpp: colonyResourceFile.qrc \
		icons/realistik_player_play_64.png \
		icons/realistik_player_pause_48.png \
		icons/realistik_player_pause_16.png \
		icons/realistik_player_play_22.png \
		icons/realistik_player_pause_128.png \
		icons/realistik_player_pause_22.png \
		icons/realistik_player_play_32.png \
		icons/PlayGreenButton.png \
		icons/realistik_player_pause_32.png \
		icons/realistik_player_play_16.png \
		icons/realistik_player_play_48.png \
		icons/StopRedButton.png \
		icons/realistik_player_pause_64.png \
		icons/realistik_player_play_128.png
	/usr/lib/x86_64-linux-gnu/qt4/bin/rcc -name colonyResourceFile colonyResourceFile.qrc -o qrc_colonyResourceFile.cpp

compiler_image_collection_make_all: qmake_image_collection.cpp
compiler_image_collection_clean:
	-$(DEL_FILE) qmake_image_collection.cpp
compiler_moc_source_make_all:
compiler_moc_source_clean:
compiler_uic_make_all: ui_mainwindow.h
compiler_uic_clean:
	-$(DEL_FILE) ui_mainwindow.h
ui_mainwindow.h: mainwindow.ui
	/usr/lib/x86_64-linux-gnu/qt4/bin/uic mainwindow.ui -o ui_mainwindow.h

compiler_yacc_decl_make_all:
compiler_yacc_decl_clean:
compiler_yacc_impl_make_all:
compiler_yacc_impl_clean:
compiler_lex_make_all:
compiler_lex_clean:
compiler_clean: compiler_moc_header_clean compiler_rcc_clean compiler_uic_clean 

####### Compile

main.o: main.cpp compilation_options.h \
		mainwindow.h \
		Simulator.h \
		Output.h \
		Input.h \
		iotools.h \
		RandomNumberGenerator.h \
		Param/CellCollectionParam.h \
		Param/CellBaseInitParam.h \
		Param/StateInitParam.h \
		Param/ChemicalSystemInitParam.h \
		Param/CellInitParam.h \
		Param/GraphicsCellCompositeParam.h \
		GraphicsCellScene.h \
		GraphicsWall.h \
		GraphicsCellGroup.h \
		GraphicsCellQt.h \
		GraphicsCellBase.h \
		Param/CellMilieuChemicalSystemInitParam.h \
		Param/OutputParam.h \
		Param/SimulatorParam.h \
		computePropensitiesFunctions.h \
		TimeSlice.h \
		TimeSliceCellsAvg.h \
		CellLineageGeneration.h \
		CellCollection.h \
		Milieu.h \
		CellBase.h \
		State.h \
		ChemicalSystem.h \
		Cell.h \
		GraphicsCellComposite.h \
		GraphicsCellODE.h \
		CellMilieuChemicalSystem.h \
		GlobalArrayInterface.h \
		GlobalArray.h \
		IntegratorContext.h \
		IntegratorGillespie.h \
		IntegratorGillespieModified.h \
		IntegratorChemicalLangevin.h \
		chemicalLangevinComputeIncrement.h \
		SpatialIntegratorContext.h \
		SpatialIntegratorODE.h \
		StringTableModel.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o main.o main.cpp

mainwindow.o: mainwindow.cpp mainwindow.h \
		compilation_options.h \
		Simulator.h \
		Output.h \
		Input.h \
		iotools.h \
		RandomNumberGenerator.h \
		Param/CellCollectionParam.h \
		Param/CellBaseInitParam.h \
		Param/StateInitParam.h \
		Param/ChemicalSystemInitParam.h \
		Param/CellInitParam.h \
		Param/GraphicsCellCompositeParam.h \
		GraphicsCellScene.h \
		GraphicsWall.h \
		GraphicsCellGroup.h \
		GraphicsCellQt.h \
		GraphicsCellBase.h \
		Param/CellMilieuChemicalSystemInitParam.h \
		Param/OutputParam.h \
		Param/SimulatorParam.h \
		computePropensitiesFunctions.h \
		TimeSlice.h \
		TimeSliceCellsAvg.h \
		CellLineageGeneration.h \
		CellCollection.h \
		Milieu.h \
		CellBase.h \
		State.h \
		ChemicalSystem.h \
		Cell.h \
		GraphicsCellComposite.h \
		GraphicsCellODE.h \
		CellMilieuChemicalSystem.h \
		GlobalArrayInterface.h \
		GlobalArray.h \
		IntegratorContext.h \
		IntegratorGillespie.h \
		IntegratorGillespieModified.h \
		IntegratorChemicalLangevin.h \
		chemicalLangevinComputeIncrement.h \
		SpatialIntegratorContext.h \
		SpatialIntegratorODE.h \
		StringTableModel.h \
		ui_mainwindow.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o mainwindow.o mainwindow.cpp

TimeSliceCellsAvg.o: TimeSliceCellsAvg.cpp TimeSliceCellsAvg.h \
		compilation_options.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o TimeSliceCellsAvg.o TimeSliceCellsAvg.cpp

TimeSlice.o: TimeSlice.cpp TimeSlice.h \
		compilation_options.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o TimeSlice.o TimeSlice.cpp

State.o: State.cpp State.h \
		compilation_options.h \
		Param/StateInitParam.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o State.o State.cpp

Simulator.o: Simulator.cpp Simulator.h \
		compilation_options.h \
		Output.h \
		Input.h \
		iotools.h \
		RandomNumberGenerator.h \
		Param/CellCollectionParam.h \
		Param/CellBaseInitParam.h \
		Param/StateInitParam.h \
		Param/ChemicalSystemInitParam.h \
		Param/CellInitParam.h \
		Param/GraphicsCellCompositeParam.h \
		GraphicsCellScene.h \
		GraphicsWall.h \
		GraphicsCellGroup.h \
		GraphicsCellQt.h \
		GraphicsCellBase.h \
		Param/CellMilieuChemicalSystemInitParam.h \
		Param/OutputParam.h \
		Param/SimulatorParam.h \
		computePropensitiesFunctions.h \
		TimeSlice.h \
		TimeSliceCellsAvg.h \
		CellLineageGeneration.h \
		CellCollection.h \
		Milieu.h \
		CellBase.h \
		State.h \
		ChemicalSystem.h \
		Cell.h \
		GraphicsCellComposite.h \
		GraphicsCellODE.h \
		CellMilieuChemicalSystem.h \
		GlobalArrayInterface.h \
		GlobalArray.h \
		IntegratorContext.h \
		IntegratorGillespie.h \
		IntegratorGillespieModified.h \
		IntegratorChemicalLangevin.h \
		chemicalLangevinComputeIncrement.h \
		SpatialIntegratorContext.h \
		SpatialIntegratorODE.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o Simulator.o Simulator.cpp

RandomNumberGenerator.o: RandomNumberGenerator.cpp RandomNumberGenerator.h \
		compilation_options.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o RandomNumberGenerator.o RandomNumberGenerator.cpp

Output.o: Output.cpp Output.h \
		compilation_options.h \
		Input.h \
		iotools.h \
		RandomNumberGenerator.h \
		Param/CellCollectionParam.h \
		Param/CellBaseInitParam.h \
		Param/StateInitParam.h \
		Param/ChemicalSystemInitParam.h \
		Param/CellInitParam.h \
		Param/GraphicsCellCompositeParam.h \
		GraphicsCellScene.h \
		GraphicsWall.h \
		GraphicsCellGroup.h \
		GraphicsCellQt.h \
		GraphicsCellBase.h \
		Param/CellMilieuChemicalSystemInitParam.h \
		Param/OutputParam.h \
		Param/SimulatorParam.h \
		computePropensitiesFunctions.h \
		TimeSlice.h \
		TimeSliceCellsAvg.h \
		CellLineageGeneration.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o Output.o Output.cpp

Milieu.o: Milieu.cpp Milieu.h \
		compilation_options.h \
		CellBase.h \
		State.h \
		Param/StateInitParam.h \
		ChemicalSystem.h \
		Param/ChemicalSystemInitParam.h \
		RandomNumberGenerator.h \
		Param/CellBaseInitParam.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o Milieu.o Milieu.cpp

iotools.o: iotools.cpp iotools.h \
		compilation_options.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o iotools.o iotools.cpp

IntegratorGillespieModified.o: IntegratorGillespieModified.cpp IntegratorGillespieModified.h \
		compilation_options.h \
		RandomNumberGenerator.h \
		Simulator.h \
		Output.h \
		Input.h \
		iotools.h \
		Param/CellCollectionParam.h \
		Param/CellBaseInitParam.h \
		Param/StateInitParam.h \
		Param/ChemicalSystemInitParam.h \
		Param/CellInitParam.h \
		Param/GraphicsCellCompositeParam.h \
		GraphicsCellScene.h \
		GraphicsWall.h \
		GraphicsCellGroup.h \
		GraphicsCellQt.h \
		GraphicsCellBase.h \
		Param/CellMilieuChemicalSystemInitParam.h \
		Param/OutputParam.h \
		Param/SimulatorParam.h \
		computePropensitiesFunctions.h \
		TimeSlice.h \
		TimeSliceCellsAvg.h \
		CellLineageGeneration.h \
		CellCollection.h \
		Milieu.h \
		CellBase.h \
		State.h \
		ChemicalSystem.h \
		Cell.h \
		GraphicsCellComposite.h \
		GraphicsCellODE.h \
		CellMilieuChemicalSystem.h \
		GlobalArrayInterface.h \
		GlobalArray.h \
		IntegratorContext.h \
		IntegratorGillespie.h \
		IntegratorChemicalLangevin.h \
		chemicalLangevinComputeIncrement.h \
		SpatialIntegratorContext.h \
		SpatialIntegratorODE.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o IntegratorGillespieModified.o IntegratorGillespieModified.cpp

IntegratorGillespie.o: IntegratorGillespie.cpp IntegratorGillespie.h \
		compilation_options.h \
		RandomNumberGenerator.h \
		Simulator.h \
		Output.h \
		Input.h \
		iotools.h \
		Param/CellCollectionParam.h \
		Param/CellBaseInitParam.h \
		Param/StateInitParam.h \
		Param/ChemicalSystemInitParam.h \
		Param/CellInitParam.h \
		Param/GraphicsCellCompositeParam.h \
		GraphicsCellScene.h \
		GraphicsWall.h \
		GraphicsCellGroup.h \
		GraphicsCellQt.h \
		GraphicsCellBase.h \
		Param/CellMilieuChemicalSystemInitParam.h \
		Param/OutputParam.h \
		Param/SimulatorParam.h \
		computePropensitiesFunctions.h \
		TimeSlice.h \
		TimeSliceCellsAvg.h \
		CellLineageGeneration.h \
		CellCollection.h \
		Milieu.h \
		CellBase.h \
		State.h \
		ChemicalSystem.h \
		Cell.h \
		GraphicsCellComposite.h \
		GraphicsCellODE.h \
		CellMilieuChemicalSystem.h \
		GlobalArrayInterface.h \
		GlobalArray.h \
		IntegratorContext.h \
		IntegratorGillespieModified.h \
		IntegratorChemicalLangevin.h \
		chemicalLangevinComputeIncrement.h \
		SpatialIntegratorContext.h \
		SpatialIntegratorODE.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o IntegratorGillespie.o IntegratorGillespie.cpp

IntegratorContext.o: IntegratorContext.cpp IntegratorContext.h \
		compilation_options.h \
		Simulator.h \
		Output.h \
		Input.h \
		iotools.h \
		RandomNumberGenerator.h \
		Param/CellCollectionParam.h \
		Param/CellBaseInitParam.h \
		Param/StateInitParam.h \
		Param/ChemicalSystemInitParam.h \
		Param/CellInitParam.h \
		Param/GraphicsCellCompositeParam.h \
		GraphicsCellScene.h \
		GraphicsWall.h \
		GraphicsCellGroup.h \
		GraphicsCellQt.h \
		GraphicsCellBase.h \
		Param/CellMilieuChemicalSystemInitParam.h \
		Param/OutputParam.h \
		Param/SimulatorParam.h \
		computePropensitiesFunctions.h \
		TimeSlice.h \
		TimeSliceCellsAvg.h \
		CellLineageGeneration.h \
		CellCollection.h \
		Milieu.h \
		CellBase.h \
		State.h \
		ChemicalSystem.h \
		Cell.h \
		GraphicsCellComposite.h \
		GraphicsCellODE.h \
		CellMilieuChemicalSystem.h \
		GlobalArrayInterface.h \
		GlobalArray.h \
		IntegratorGillespie.h \
		IntegratorGillespieModified.h \
		IntegratorChemicalLangevin.h \
		chemicalLangevinComputeIncrement.h \
		SpatialIntegratorContext.h \
		SpatialIntegratorODE.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o IntegratorContext.o IntegratorContext.cpp

Input.o: Input.cpp Input.h \
		compilation_options.h \
		iotools.h \
		RandomNumberGenerator.h \
		Param/CellCollectionParam.h \
		Param/CellBaseInitParam.h \
		Param/StateInitParam.h \
		Param/ChemicalSystemInitParam.h \
		Param/CellInitParam.h \
		Param/GraphicsCellCompositeParam.h \
		GraphicsCellScene.h \
		GraphicsWall.h \
		GraphicsCellGroup.h \
		GraphicsCellQt.h \
		GraphicsCellBase.h \
		Param/CellMilieuChemicalSystemInitParam.h \
		Param/OutputParam.h \
		Param/SimulatorParam.h \
		computePropensitiesFunctions.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o Input.o Input.cpp

GlobalArrayInterface.o: GlobalArrayInterface.cpp GlobalArrayInterface.h \
		compilation_options.h \
		GlobalArray.h \
		CellCollection.h \
		Milieu.h \
		CellBase.h \
		State.h \
		Param/StateInitParam.h \
		ChemicalSystem.h \
		Param/ChemicalSystemInitParam.h \
		RandomNumberGenerator.h \
		Param/CellBaseInitParam.h \
		CellLineageGeneration.h \
		Cell.h \
		GraphicsCellComposite.h \
		GraphicsCellBase.h \
		GraphicsCellODE.h \
		GraphicsCellQt.h \
		GraphicsCellScene.h \
		GraphicsWall.h \
		GraphicsCellGroup.h \
		Param/GraphicsCellCompositeParam.h \
		CellMilieuChemicalSystem.h \
		Param/CellMilieuChemicalSystemInitParam.h \
		Param/CellInitParam.h \
		Param/CellCollectionParam.h \
		Input.h \
		iotools.h \
		Param/OutputParam.h \
		Param/SimulatorParam.h \
		computePropensitiesFunctions.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o GlobalArrayInterface.o GlobalArrayInterface.cpp

GlobalArray.o: GlobalArray.cpp GlobalArray.h \
		compilation_options.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o GlobalArray.o GlobalArray.cpp

ChemicalSystem.o: ChemicalSystem.cpp ChemicalSystem.h \
		compilation_options.h \
		Param/ChemicalSystemInitParam.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o ChemicalSystem.o ChemicalSystem.cpp

CellMilieuChemicalSystem.o: CellMilieuChemicalSystem.cpp CellMilieuChemicalSystem.h \
		compilation_options.h \
		Param/CellMilieuChemicalSystemInitParam.h \
		Cell.h \
		CellBase.h \
		State.h \
		Param/StateInitParam.h \
		ChemicalSystem.h \
		Param/ChemicalSystemInitParam.h \
		RandomNumberGenerator.h \
		Param/CellBaseInitParam.h \
		GraphicsCellComposite.h \
		GraphicsCellBase.h \
		GraphicsCellODE.h \
		GraphicsCellQt.h \
		GraphicsCellScene.h \
		GraphicsWall.h \
		GraphicsCellGroup.h \
		Param/GraphicsCellCompositeParam.h \
		Param/CellInitParam.h \
		Milieu.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o CellMilieuChemicalSystem.o CellMilieuChemicalSystem.cpp

CellLineageGeneration.o: CellLineageGeneration.cpp CellLineageGeneration.h \
		compilation_options.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o CellLineageGeneration.o CellLineageGeneration.cpp

CellCollection.o: CellCollection.cpp CellCollection.h \
		compilation_options.h \
		Milieu.h \
		CellBase.h \
		State.h \
		Param/StateInitParam.h \
		ChemicalSystem.h \
		Param/ChemicalSystemInitParam.h \
		RandomNumberGenerator.h \
		Param/CellBaseInitParam.h \
		CellLineageGeneration.h \
		Cell.h \
		GraphicsCellComposite.h \
		GraphicsCellBase.h \
		GraphicsCellODE.h \
		GraphicsCellQt.h \
		GraphicsCellScene.h \
		GraphicsWall.h \
		GraphicsCellGroup.h \
		Param/GraphicsCellCompositeParam.h \
		CellMilieuChemicalSystem.h \
		Param/CellMilieuChemicalSystemInitParam.h \
		Param/CellInitParam.h \
		GlobalArrayInterface.h \
		GlobalArray.h \
		Param/CellCollectionParam.h \
		Input.h \
		iotools.h \
		Param/OutputParam.h \
		Param/SimulatorParam.h \
		computePropensitiesFunctions.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o CellCollection.o CellCollection.cpp

CellBase.o: CellBase.cpp CellBase.h \
		compilation_options.h \
		State.h \
		Param/StateInitParam.h \
		ChemicalSystem.h \
		Param/ChemicalSystemInitParam.h \
		RandomNumberGenerator.h \
		Param/CellBaseInitParam.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o CellBase.o CellBase.cpp

Cell.o: Cell.cpp Cell.h \
		compilation_options.h \
		CellBase.h \
		State.h \
		Param/StateInitParam.h \
		ChemicalSystem.h \
		Param/ChemicalSystemInitParam.h \
		RandomNumberGenerator.h \
		Param/CellBaseInitParam.h \
		GraphicsCellComposite.h \
		GraphicsCellBase.h \
		GraphicsCellODE.h \
		GraphicsCellQt.h \
		GraphicsCellScene.h \
		GraphicsWall.h \
		GraphicsCellGroup.h \
		Param/GraphicsCellCompositeParam.h \
		CellMilieuChemicalSystem.h \
		Param/CellMilieuChemicalSystemInitParam.h \
		Milieu.h \
		Param/CellInitParam.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o Cell.o Cell.cpp

computePropensitiesFunctions.o: computePropensitiesFunctions.cpp compilation_options.h \
		computePropensitiesFunctions.h \
		Input.h \
		iotools.h \
		RandomNumberGenerator.h \
		Param/CellCollectionParam.h \
		Param/CellBaseInitParam.h \
		Param/StateInitParam.h \
		Param/ChemicalSystemInitParam.h \
		Param/CellInitParam.h \
		Param/GraphicsCellCompositeParam.h \
		GraphicsCellScene.h \
		GraphicsWall.h \
		GraphicsCellGroup.h \
		GraphicsCellQt.h \
		GraphicsCellBase.h \
		Param/CellMilieuChemicalSystemInitParam.h \
		Param/OutputParam.h \
		Param/SimulatorParam.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o computePropensitiesFunctions.o computePropensitiesFunctions.cpp

StringTableModel.o: StringTableModel.cpp StringTableModel.h \
		compilation_options.h \
		Simulator.h \
		Output.h \
		Input.h \
		iotools.h \
		RandomNumberGenerator.h \
		Param/CellCollectionParam.h \
		Param/CellBaseInitParam.h \
		Param/StateInitParam.h \
		Param/ChemicalSystemInitParam.h \
		Param/CellInitParam.h \
		Param/GraphicsCellCompositeParam.h \
		GraphicsCellScene.h \
		GraphicsWall.h \
		GraphicsCellGroup.h \
		GraphicsCellQt.h \
		GraphicsCellBase.h \
		Param/CellMilieuChemicalSystemInitParam.h \
		Param/OutputParam.h \
		Param/SimulatorParam.h \
		computePropensitiesFunctions.h \
		TimeSlice.h \
		TimeSliceCellsAvg.h \
		CellLineageGeneration.h \
		CellCollection.h \
		Milieu.h \
		CellBase.h \
		State.h \
		ChemicalSystem.h \
		Cell.h \
		GraphicsCellComposite.h \
		GraphicsCellODE.h \
		CellMilieuChemicalSystem.h \
		GlobalArrayInterface.h \
		GlobalArray.h \
		IntegratorContext.h \
		IntegratorGillespie.h \
		IntegratorGillespieModified.h \
		IntegratorChemicalLangevin.h \
		chemicalLangevinComputeIncrement.h \
		SpatialIntegratorContext.h \
		SpatialIntegratorODE.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o StringTableModel.o StringTableModel.cpp

GraphicsCellScene.o: GraphicsCellScene.cpp GraphicsCellScene.h \
		compilation_options.h \
		GraphicsWall.h \
		GraphicsCellGroup.h \
		GraphicsCellQt.h \
		GraphicsCellBase.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o GraphicsCellScene.o GraphicsCellScene.cpp

GraphicsCellGroup.o: GraphicsCellGroup.cpp GraphicsCellGroup.h \
		compilation_options.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o GraphicsCellGroup.o GraphicsCellGroup.cpp

GraphicsWall.o: GraphicsWall.cpp GraphicsWall.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o GraphicsWall.o GraphicsWall.cpp

GraphicsCellBase.o: GraphicsCellBase.cpp GraphicsCellBase.h \
		compilation_options.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o GraphicsCellBase.o GraphicsCellBase.cpp

GraphicsCellQt.o: GraphicsCellQt.cpp GraphicsCellQt.h \
		compilation_options.h \
		GraphicsCellBase.h \
		GraphicsCellScene.h \
		GraphicsWall.h \
		GraphicsCellGroup.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o GraphicsCellQt.o GraphicsCellQt.cpp

GraphicsCellComposite.o: GraphicsCellComposite.cpp GraphicsCellComposite.h \
		compilation_options.h \
		GraphicsCellBase.h \
		GraphicsCellODE.h \
		GraphicsCellQt.h \
		GraphicsCellScene.h \
		GraphicsWall.h \
		GraphicsCellGroup.h \
		Param/GraphicsCellCompositeParam.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o GraphicsCellComposite.o GraphicsCellComposite.cpp

SpatialIntegratorContext.o: SpatialIntegratorContext.cpp SpatialIntegratorContext.h \
		compilation_options.h \
		Simulator.h \
		Output.h \
		Input.h \
		iotools.h \
		RandomNumberGenerator.h \
		Param/CellCollectionParam.h \
		Param/CellBaseInitParam.h \
		Param/StateInitParam.h \
		Param/ChemicalSystemInitParam.h \
		Param/CellInitParam.h \
		Param/GraphicsCellCompositeParam.h \
		GraphicsCellScene.h \
		GraphicsWall.h \
		GraphicsCellGroup.h \
		GraphicsCellQt.h \
		GraphicsCellBase.h \
		Param/CellMilieuChemicalSystemInitParam.h \
		Param/OutputParam.h \
		Param/SimulatorParam.h \
		computePropensitiesFunctions.h \
		TimeSlice.h \
		TimeSliceCellsAvg.h \
		CellLineageGeneration.h \
		CellCollection.h \
		Milieu.h \
		CellBase.h \
		State.h \
		ChemicalSystem.h \
		Cell.h \
		GraphicsCellComposite.h \
		GraphicsCellODE.h \
		CellMilieuChemicalSystem.h \
		GlobalArrayInterface.h \
		GlobalArray.h \
		IntegratorContext.h \
		IntegratorGillespie.h \
		IntegratorGillespieModified.h \
		IntegratorChemicalLangevin.h \
		chemicalLangevinComputeIncrement.h \
		SpatialIntegratorODE.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o SpatialIntegratorContext.o SpatialIntegratorContext.cpp

SpatialIntegratorODE.o: SpatialIntegratorODE.cpp SpatialIntegratorODE.h \
		compilation_options.h \
		Simulator.h \
		Output.h \
		Input.h \
		iotools.h \
		RandomNumberGenerator.h \
		Param/CellCollectionParam.h \
		Param/CellBaseInitParam.h \
		Param/StateInitParam.h \
		Param/ChemicalSystemInitParam.h \
		Param/CellInitParam.h \
		Param/GraphicsCellCompositeParam.h \
		GraphicsCellScene.h \
		GraphicsWall.h \
		GraphicsCellGroup.h \
		GraphicsCellQt.h \
		GraphicsCellBase.h \
		Param/CellMilieuChemicalSystemInitParam.h \
		Param/OutputParam.h \
		Param/SimulatorParam.h \
		computePropensitiesFunctions.h \
		TimeSlice.h \
		TimeSliceCellsAvg.h \
		CellLineageGeneration.h \
		CellCollection.h \
		Milieu.h \
		CellBase.h \
		State.h \
		ChemicalSystem.h \
		Cell.h \
		GraphicsCellComposite.h \
		GraphicsCellODE.h \
		CellMilieuChemicalSystem.h \
		GlobalArrayInterface.h \
		GlobalArray.h \
		IntegratorContext.h \
		IntegratorGillespie.h \
		IntegratorGillespieModified.h \
		IntegratorChemicalLangevin.h \
		chemicalLangevinComputeIncrement.h \
		SpatialIntegratorContext.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o SpatialIntegratorODE.o SpatialIntegratorODE.cpp

GraphicsCellODE.o: GraphicsCellODE.cpp GraphicsCellODE.h \
		compilation_options.h \
		GraphicsCellBase.h \
		SpatialIntegratorODE.h \
		Simulator.h \
		Output.h \
		Input.h \
		iotools.h \
		RandomNumberGenerator.h \
		Param/CellCollectionParam.h \
		Param/CellBaseInitParam.h \
		Param/StateInitParam.h \
		Param/ChemicalSystemInitParam.h \
		Param/CellInitParam.h \
		Param/GraphicsCellCompositeParam.h \
		GraphicsCellScene.h \
		GraphicsWall.h \
		GraphicsCellGroup.h \
		GraphicsCellQt.h \
		Param/CellMilieuChemicalSystemInitParam.h \
		Param/OutputParam.h \
		Param/SimulatorParam.h \
		computePropensitiesFunctions.h \
		TimeSlice.h \
		TimeSliceCellsAvg.h \
		CellLineageGeneration.h \
		CellCollection.h \
		Milieu.h \
		CellBase.h \
		State.h \
		ChemicalSystem.h \
		Cell.h \
		GraphicsCellComposite.h \
		CellMilieuChemicalSystem.h \
		GlobalArrayInterface.h \
		GlobalArray.h \
		IntegratorContext.h \
		IntegratorGillespie.h \
		IntegratorGillespieModified.h \
		IntegratorChemicalLangevin.h \
		chemicalLangevinComputeIncrement.h \
		SpatialIntegratorContext.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o GraphicsCellODE.o GraphicsCellODE.cpp

IntegratorChemicalLangevin.o: IntegratorChemicalLangevin.cpp compilation_options.h \
		IntegratorChemicalLangevin.h \
		RandomNumberGenerator.h \
		chemicalLangevinComputeIncrement.h \
		Input.h \
		iotools.h \
		Param/CellCollectionParam.h \
		Param/CellBaseInitParam.h \
		Param/StateInitParam.h \
		Param/ChemicalSystemInitParam.h \
		Param/CellInitParam.h \
		Param/GraphicsCellCompositeParam.h \
		GraphicsCellScene.h \
		GraphicsWall.h \
		GraphicsCellGroup.h \
		GraphicsCellQt.h \
		GraphicsCellBase.h \
		Param/CellMilieuChemicalSystemInitParam.h \
		Param/OutputParam.h \
		Param/SimulatorParam.h \
		computePropensitiesFunctions.h \
		Simulator.h \
		Output.h \
		TimeSlice.h \
		TimeSliceCellsAvg.h \
		CellLineageGeneration.h \
		CellCollection.h \
		Milieu.h \
		CellBase.h \
		State.h \
		ChemicalSystem.h \
		Cell.h \
		GraphicsCellComposite.h \
		GraphicsCellODE.h \
		CellMilieuChemicalSystem.h \
		GlobalArrayInterface.h \
		GlobalArray.h \
		IntegratorContext.h \
		IntegratorGillespie.h \
		IntegratorGillespieModified.h \
		SpatialIntegratorContext.h \
		SpatialIntegratorODE.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o IntegratorChemicalLangevin.o IntegratorChemicalLangevin.cpp

chemicalLangevinComputeIncrement.o: chemicalLangevinComputeIncrement.cpp compilation_options.h \
		chemicalLangevinComputeIncrement.h \
		Input.h \
		iotools.h \
		RandomNumberGenerator.h \
		Param/CellCollectionParam.h \
		Param/CellBaseInitParam.h \
		Param/StateInitParam.h \
		Param/ChemicalSystemInitParam.h \
		Param/CellInitParam.h \
		Param/GraphicsCellCompositeParam.h \
		GraphicsCellScene.h \
		GraphicsWall.h \
		GraphicsCellGroup.h \
		GraphicsCellQt.h \
		GraphicsCellBase.h \
		Param/CellMilieuChemicalSystemInitParam.h \
		Param/OutputParam.h \
		Param/SimulatorParam.h \
		computePropensitiesFunctions.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o chemicalLangevinComputeIncrement.o chemicalLangevinComputeIncrement.cpp

moc_mainwindow.o: moc_mainwindow.cpp 
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o moc_mainwindow.o moc_mainwindow.cpp

moc_StringTableModel.o: moc_StringTableModel.cpp 
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o moc_StringTableModel.o moc_StringTableModel.cpp

moc_GraphicsCellScene.o: moc_GraphicsCellScene.cpp 
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o moc_GraphicsCellScene.o moc_GraphicsCellScene.cpp

moc_GraphicsWall.o: moc_GraphicsWall.cpp 
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o moc_GraphicsWall.o moc_GraphicsWall.cpp

moc_GraphicsCellQt.o: moc_GraphicsCellQt.cpp 
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o moc_GraphicsCellQt.o moc_GraphicsCellQt.cpp

qrc_colonyResourceFile.o: qrc_colonyResourceFile.cpp 
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o qrc_colonyResourceFile.o qrc_colonyResourceFile.cpp

####### Install

install:   FORCE

uninstall:   FORCE

FORCE:
