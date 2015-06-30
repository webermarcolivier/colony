/****************************************************************************
** Meta object code from reading C++ file 'GraphicsCellQt.h'
**
** Created: Tue May 24 15:05:25 2011
**      by: The Qt Meta Object Compiler version 62 (Qt 4.6.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "GraphicsCellQt.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'GraphicsCellQt.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.6.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_GraphicsCellQt[] = {

 // content:
       4,       // revision
       0,       // classname
       0,    0, // classinfo
       2,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       2,       // signalCount

 // signals: signature, parameters, type, tag, flags
      21,   16,   15,   15, 0x05,
      50,   16,   15,   15, 0x05,

       0        // eod
};

static const char qt_meta_stringdata_GraphicsCellQt[] = {
    "GraphicsCellQt\0\0cell\0gainedFocus(GraphicsCellQt*)\0"
    "lostFocus(GraphicsCellQt*)\0"
};

const QMetaObject GraphicsCellQt::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_GraphicsCellQt,
      qt_meta_data_GraphicsCellQt, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &GraphicsCellQt::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *GraphicsCellQt::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *GraphicsCellQt::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_GraphicsCellQt))
        return static_cast<void*>(const_cast< GraphicsCellQt*>(this));
    if (!strcmp(_clname, "QGraphicsItem"))
        return static_cast< QGraphicsItem*>(const_cast< GraphicsCellQt*>(this));
    if (!strcmp(_clname, "GraphicsCellBase"))
        return static_cast< GraphicsCellBase*>(const_cast< GraphicsCellQt*>(this));
    return QObject::qt_metacast(_clname);
}

int GraphicsCellQt::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: gainedFocus((*reinterpret_cast< GraphicsCellQt*(*)>(_a[1]))); break;
        case 1: lostFocus((*reinterpret_cast< GraphicsCellQt*(*)>(_a[1]))); break;
        default: ;
        }
        _id -= 2;
    }
    return _id;
}

// SIGNAL 0
void GraphicsCellQt::gainedFocus(GraphicsCellQt * _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void GraphicsCellQt::lostFocus(GraphicsCellQt * _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}
QT_END_MOC_NAMESPACE
