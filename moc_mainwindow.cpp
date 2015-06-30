/****************************************************************************
** Meta object code from reading C++ file 'mainwindow.h'
**
** Created: Tue May 24 15:05:17 2011
**      by: The Qt Meta Object Compiler version 62 (Qt 4.6.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "mainwindow.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'mainwindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.6.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_MainWindow[] = {

 // content:
       4,       // revision
       0,       // classname
       0,    0, // classinfo
      18,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       6,       // signalCount

 // signals: signature, parameters, type, tag, flags
      19,   12,   11,   11, 0x05,
      37,   11,   11,   11, 0x05,
      67,   11,   11,   11, 0x05,
     102,   11,   11,   11, 0x05,
     135,   11,   11,   11, 0x05,
     170,   11,   11,   11, 0x05,

 // slots: signature, parameters, type, tag, flags
     203,   11,   11,   11, 0x0a,
     219,   11,   11,   11, 0x0a,
     241,   11,   11,   11, 0x0a,
     258,   11,   11,   11, 0x0a,
     276,   11,   11,   11, 0x0a,
     305,  295,   11,   11, 0x0a,
     343,  337,   11,   11, 0x0a,
     377,  373,   11,   11, 0x0a,
     413,  337,   11,   11, 0x0a,
     449,  373,   11,   11, 0x0a,
     496,  491,   11,   11, 0x0a,
     527,  491,   11,   11, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_MainWindow[] = {
    "MainWindow\0\0nCells\0updateNCells(int)\0"
    "updateTimeSimulation(QString)\0"
    "updateTotalTimeSimulation(QString)\0"
    "updateTrajectoryCounter(QString)\0"
    "updateParameterSetCounter(QString)\0"
    "updateSimulationProgressBar(int)\0"
    "runSimulation()\0playPauseSimulation()\0"
    "playSimulation()\0pauseSimulation()\0"
    "openOpenGLViewer()\0zoomValue\0"
    "on_zoomSlider_valueChanged(int)\0index\0"
    "setSpeciesCellColorIndex(int)\0max\0"
    "setSpeciesCellColorMaxValue(double)\0"
    "setSpeciesBackgroundColorIndex(int)\0"
    "setSpeciesBackgroundColorMaxValue(double)\0"
    "cell\0cellLostFocus(GraphicsCellQt*)\0"
    "cellGainedFocus(GraphicsCellQt*)\0"
};

const QMetaObject MainWindow::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_MainWindow,
      qt_meta_data_MainWindow, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &MainWindow::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *MainWindow::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *MainWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_MainWindow))
        return static_cast<void*>(const_cast< MainWindow*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int MainWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: updateNCells((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: updateTimeSimulation((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 2: updateTotalTimeSimulation((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 3: updateTrajectoryCounter((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 4: updateParameterSetCounter((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 5: updateSimulationProgressBar((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 6: runSimulation(); break;
        case 7: playPauseSimulation(); break;
        case 8: playSimulation(); break;
        case 9: pauseSimulation(); break;
        case 10: openOpenGLViewer(); break;
        case 11: on_zoomSlider_valueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 12: setSpeciesCellColorIndex((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 13: setSpeciesCellColorMaxValue((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 14: setSpeciesBackgroundColorIndex((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 15: setSpeciesBackgroundColorMaxValue((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 16: cellLostFocus((*reinterpret_cast< GraphicsCellQt*(*)>(_a[1]))); break;
        case 17: cellGainedFocus((*reinterpret_cast< GraphicsCellQt*(*)>(_a[1]))); break;
        default: ;
        }
        _id -= 18;
    }
    return _id;
}

// SIGNAL 0
void MainWindow::updateNCells(int _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void MainWindow::updateTimeSimulation(QString _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}

// SIGNAL 2
void MainWindow::updateTotalTimeSimulation(QString _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 2, _a);
}

// SIGNAL 3
void MainWindow::updateTrajectoryCounter(QString _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 3, _a);
}

// SIGNAL 4
void MainWindow::updateParameterSetCounter(QString _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 4, _a);
}

// SIGNAL 5
void MainWindow::updateSimulationProgressBar(int _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 5, _a);
}
QT_END_MOC_NAMESPACE
