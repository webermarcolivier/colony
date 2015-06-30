/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created: Tue May 24 15:04:29 2011
**      by: Qt User Interface Compiler version 4.6.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QComboBox>
#include <QtGui/QDoubleSpinBox>
#include <QtGui/QGraphicsView>
#include <QtGui/QGridLayout>
#include <QtGui/QGroupBox>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QMainWindow>
#include <QtGui/QProgressBar>
#include <QtGui/QPushButton>
#include <QtGui/QSlider>
#include <QtGui/QStatusBar>
#include <QtGui/QTableView>
#include <QtGui/QToolBar>
#include <QtGui/QToolButton>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QWidget *centralWidget;
    QGridLayout *gridLayout_3;
    QGridLayout *gridLayout_2;
    QHBoxLayout *horizontalLayout;
    QToolButton *playButton;
    QToolButton *pauseButton;
    QProgressBar *simulationProgressBar;
    QPushButton *pushButton;
    QGridLayout *gridLayout;
    QLabel *NCells;
    QLabel *NCellsNum;
    QLabel *time;
    QLabel *timeNum;
    QLabel *iTraj;
    QLabel *iTrajNum;
    QLabel *iParameterSet;
    QLabel *iParameterSetNum;
    QVBoxLayout *verticalLayout;
    QGroupBox *groupBox;
    QHBoxLayout *horizontalLayout_5;
    QHBoxLayout *horizontalLayout_2;
    QGraphicsView *colorCodeGraphicsView;
    QDoubleSpinBox *speciesColorMaxValueSpinBox;
    QComboBox *speciesColorComboBox;
    QGroupBox *groupBox_2;
    QHBoxLayout *horizontalLayout_4;
    QHBoxLayout *horizontalLayout_3;
    QGraphicsView *colorCodeBackgroundGraphicsView;
    QDoubleSpinBox *speciesBackgroundColorMaxValueSpinBox;
    QComboBox *speciesBackgroundColorComboBox;
    QLabel *label_2;
    QLabel *label_3;
    QTableView *globalParameterTableView;
    QGraphicsView *colonyView;
    QTableView *cellStateTableView;
    QHBoxLayout *horizontalLayout_6;
    QLabel *label;
    QSlider *zoomSlider;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(1037, 767);
        centralWidget = new QWidget(MainWindow);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        gridLayout_3 = new QGridLayout(centralWidget);
        gridLayout_3->setSpacing(6);
        gridLayout_3->setContentsMargins(11, 11, 11, 11);
        gridLayout_3->setObjectName(QString::fromUtf8("gridLayout_3"));
        gridLayout_2 = new QGridLayout();
        gridLayout_2->setSpacing(6);
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setSpacing(6);
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        playButton = new QToolButton(centralWidget);
        playButton->setObjectName(QString::fromUtf8("playButton"));
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/icons/realistik_player_play_48.png"), QSize(), QIcon::Normal, QIcon::Off);
        playButton->setIcon(icon);
        playButton->setIconSize(QSize(48, 48));
        playButton->setAutoRaise(true);

        horizontalLayout->addWidget(playButton);

        pauseButton = new QToolButton(centralWidget);
        pauseButton->setObjectName(QString::fromUtf8("pauseButton"));
        QSizePolicy sizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(pauseButton->sizePolicy().hasHeightForWidth());
        pauseButton->setSizePolicy(sizePolicy);
        QIcon icon1;
        icon1.addFile(QString::fromUtf8(":/icons/realistik_player_pause_48.png"), QSize(), QIcon::Normal, QIcon::Off);
        pauseButton->setIcon(icon1);
        pauseButton->setIconSize(QSize(48, 48));
        pauseButton->setAutoRaise(true);

        horizontalLayout->addWidget(pauseButton);

        simulationProgressBar = new QProgressBar(centralWidget);
        simulationProgressBar->setObjectName(QString::fromUtf8("simulationProgressBar"));
        simulationProgressBar->setMaximumSize(QSize(300, 16777215));
        simulationProgressBar->setValue(0);

        horizontalLayout->addWidget(simulationProgressBar);


        gridLayout_2->addLayout(horizontalLayout, 0, 0, 1, 1);

        pushButton = new QPushButton(centralWidget);
        pushButton->setObjectName(QString::fromUtf8("pushButton"));
        sizePolicy.setHeightForWidth(pushButton->sizePolicy().hasHeightForWidth());
        pushButton->setSizePolicy(sizePolicy);
        pushButton->setMinimumSize(QSize(50, 15));
        pushButton->setMaximumSize(QSize(120, 30));
        pushButton->setBaseSize(QSize(120, 31));

        gridLayout_2->addWidget(pushButton, 0, 2, 1, 1);

        gridLayout = new QGridLayout();
        gridLayout->setSpacing(6);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        gridLayout->setContentsMargins(5, -1, 10, -1);
        NCells = new QLabel(centralWidget);
        NCells->setObjectName(QString::fromUtf8("NCells"));
        QFont font;
        font.setFamily(QString::fromUtf8("Sans Serif"));
        font.setPointSize(10);
        NCells->setFont(font);

        gridLayout->addWidget(NCells, 1, 0, 1, 1);

        NCellsNum = new QLabel(centralWidget);
        NCellsNum->setObjectName(QString::fromUtf8("NCellsNum"));
        QFont font1;
        font1.setFamily(QString::fromUtf8("Monospace"));
        font1.setPointSize(10);
        NCellsNum->setFont(font1);
        NCellsNum->setLayoutDirection(Qt::LeftToRight);
        NCellsNum->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(NCellsNum, 1, 1, 1, 1);

        time = new QLabel(centralWidget);
        time->setObjectName(QString::fromUtf8("time"));
        time->setFont(font);

        gridLayout->addWidget(time, 2, 0, 1, 1);

        timeNum = new QLabel(centralWidget);
        timeNum->setObjectName(QString::fromUtf8("timeNum"));
        timeNum->setFont(font1);
        timeNum->setLayoutDirection(Qt::LeftToRight);
        timeNum->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(timeNum, 2, 1, 1, 1);

        iTraj = new QLabel(centralWidget);
        iTraj->setObjectName(QString::fromUtf8("iTraj"));

        gridLayout->addWidget(iTraj, 3, 0, 1, 1);

        iTrajNum = new QLabel(centralWidget);
        iTrajNum->setObjectName(QString::fromUtf8("iTrajNum"));
        QFont font2;
        font2.setFamily(QString::fromUtf8("Monospace"));
        iTrajNum->setFont(font2);
        iTrajNum->setLayoutDirection(Qt::LeftToRight);
        iTrajNum->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(iTrajNum, 3, 1, 1, 1);

        iParameterSet = new QLabel(centralWidget);
        iParameterSet->setObjectName(QString::fromUtf8("iParameterSet"));

        gridLayout->addWidget(iParameterSet, 4, 0, 1, 1);

        iParameterSetNum = new QLabel(centralWidget);
        iParameterSetNum->setObjectName(QString::fromUtf8("iParameterSetNum"));
        iParameterSetNum->setFont(font2);
        iParameterSetNum->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(iParameterSetNum, 4, 1, 1, 1);


        gridLayout_2->addLayout(gridLayout, 1, 0, 1, 1);

        verticalLayout = new QVBoxLayout();
        verticalLayout->setSpacing(6);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        groupBox = new QGroupBox(centralWidget);
        groupBox->setObjectName(QString::fromUtf8("groupBox"));
        QSizePolicy sizePolicy1(QSizePolicy::Preferred, QSizePolicy::Fixed);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(groupBox->sizePolicy().hasHeightForWidth());
        groupBox->setSizePolicy(sizePolicy1);
        groupBox->setMinimumSize(QSize(0, 0));
        groupBox->setMaximumSize(QSize(700, 16777215));
        groupBox->setAlignment(Qt::AlignLeading|Qt::AlignLeft|Qt::AlignTop);
        horizontalLayout_5 = new QHBoxLayout(groupBox);
        horizontalLayout_5->setSpacing(6);
        horizontalLayout_5->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_5->setObjectName(QString::fromUtf8("horizontalLayout_5"));
        horizontalLayout_5->setContentsMargins(-1, 5, -1, 5);
        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setSpacing(6);
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        colorCodeGraphicsView = new QGraphicsView(groupBox);
        colorCodeGraphicsView->setObjectName(QString::fromUtf8("colorCodeGraphicsView"));
        sizePolicy1.setHeightForWidth(colorCodeGraphicsView->sizePolicy().hasHeightForWidth());
        colorCodeGraphicsView->setSizePolicy(sizePolicy1);
        colorCodeGraphicsView->setMinimumSize(QSize(30, 30));
        colorCodeGraphicsView->setMaximumSize(QSize(200, 16777215));
        colorCodeGraphicsView->setBaseSize(QSize(150, 30));
        colorCodeGraphicsView->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        colorCodeGraphicsView->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        colorCodeGraphicsView->setInteractive(false);

        horizontalLayout_2->addWidget(colorCodeGraphicsView);

        speciesColorMaxValueSpinBox = new QDoubleSpinBox(groupBox);
        speciesColorMaxValueSpinBox->setObjectName(QString::fromUtf8("speciesColorMaxValueSpinBox"));
        sizePolicy1.setHeightForWidth(speciesColorMaxValueSpinBox->sizePolicy().hasHeightForWidth());
        speciesColorMaxValueSpinBox->setSizePolicy(sizePolicy1);
        speciesColorMaxValueSpinBox->setMinimumSize(QSize(60, 30));
        speciesColorMaxValueSpinBox->setBaseSize(QSize(60, 30));
        speciesColorMaxValueSpinBox->setAccelerated(true);
        speciesColorMaxValueSpinBox->setDecimals(0);
        speciesColorMaxValueSpinBox->setMaximum(10000);
        speciesColorMaxValueSpinBox->setSingleStep(10);
        speciesColorMaxValueSpinBox->setValue(100);

        horizontalLayout_2->addWidget(speciesColorMaxValueSpinBox);

        speciesColorComboBox = new QComboBox(groupBox);
        speciesColorComboBox->setObjectName(QString::fromUtf8("speciesColorComboBox"));
        speciesColorComboBox->setMinimumSize(QSize(60, 30));
        speciesColorComboBox->setBaseSize(QSize(60, 30));

        horizontalLayout_2->addWidget(speciesColorComboBox);

        horizontalLayout_2->setStretch(0, 2);
        horizontalLayout_2->setStretch(1, 1);
        horizontalLayout_2->setStretch(2, 1);

        horizontalLayout_5->addLayout(horizontalLayout_2);


        verticalLayout->addWidget(groupBox);

        groupBox_2 = new QGroupBox(centralWidget);
        groupBox_2->setObjectName(QString::fromUtf8("groupBox_2"));
        sizePolicy1.setHeightForWidth(groupBox_2->sizePolicy().hasHeightForWidth());
        groupBox_2->setSizePolicy(sizePolicy1);
        groupBox_2->setMaximumSize(QSize(700, 16777215));
        horizontalLayout_4 = new QHBoxLayout(groupBox_2);
        horizontalLayout_4->setSpacing(6);
        horizontalLayout_4->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_4->setObjectName(QString::fromUtf8("horizontalLayout_4"));
        horizontalLayout_4->setContentsMargins(-1, 5, -1, 5);
        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setSpacing(6);
        horizontalLayout_3->setObjectName(QString::fromUtf8("horizontalLayout_3"));
        colorCodeBackgroundGraphicsView = new QGraphicsView(groupBox_2);
        colorCodeBackgroundGraphicsView->setObjectName(QString::fromUtf8("colorCodeBackgroundGraphicsView"));
        sizePolicy1.setHeightForWidth(colorCodeBackgroundGraphicsView->sizePolicy().hasHeightForWidth());
        colorCodeBackgroundGraphicsView->setSizePolicy(sizePolicy1);
        colorCodeBackgroundGraphicsView->setMinimumSize(QSize(30, 30));
        colorCodeBackgroundGraphicsView->setMaximumSize(QSize(200, 16777215));
        colorCodeBackgroundGraphicsView->setBaseSize(QSize(150, 30));
        colorCodeBackgroundGraphicsView->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        colorCodeBackgroundGraphicsView->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        colorCodeBackgroundGraphicsView->setInteractive(false);

        horizontalLayout_3->addWidget(colorCodeBackgroundGraphicsView);

        speciesBackgroundColorMaxValueSpinBox = new QDoubleSpinBox(groupBox_2);
        speciesBackgroundColorMaxValueSpinBox->setObjectName(QString::fromUtf8("speciesBackgroundColorMaxValueSpinBox"));
        sizePolicy1.setHeightForWidth(speciesBackgroundColorMaxValueSpinBox->sizePolicy().hasHeightForWidth());
        speciesBackgroundColorMaxValueSpinBox->setSizePolicy(sizePolicy1);
        speciesBackgroundColorMaxValueSpinBox->setMinimumSize(QSize(60, 30));
        speciesBackgroundColorMaxValueSpinBox->setDecimals(0);
        speciesBackgroundColorMaxValueSpinBox->setMaximum(10000);
        speciesBackgroundColorMaxValueSpinBox->setSingleStep(10);
        speciesBackgroundColorMaxValueSpinBox->setValue(100);

        horizontalLayout_3->addWidget(speciesBackgroundColorMaxValueSpinBox);

        speciesBackgroundColorComboBox = new QComboBox(groupBox_2);
        speciesBackgroundColorComboBox->setObjectName(QString::fromUtf8("speciesBackgroundColorComboBox"));
        speciesBackgroundColorComboBox->setMinimumSize(QSize(60, 30));
        speciesBackgroundColorComboBox->setBaseSize(QSize(60, 30));

        horizontalLayout_3->addWidget(speciesBackgroundColorComboBox);

        horizontalLayout_3->setStretch(0, 2);
        horizontalLayout_3->setStretch(1, 1);
        horizontalLayout_3->setStretch(2, 1);

        horizontalLayout_4->addLayout(horizontalLayout_3);


        verticalLayout->addWidget(groupBox_2);


        gridLayout_2->addLayout(verticalLayout, 1, 1, 1, 1);

        label_2 = new QLabel(centralWidget);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        QFont font3;
        font3.setFamily(QString::fromUtf8("Sans Serif"));
        font3.setBold(true);
        font3.setWeight(75);
        label_2->setFont(font3);
        label_2->setIndent(10);

        gridLayout_2->addWidget(label_2, 2, 0, 1, 1);

        label_3 = new QLabel(centralWidget);
        label_3->setObjectName(QString::fromUtf8("label_3"));
        label_3->setFont(font3);
        label_3->setIndent(10);

        gridLayout_2->addWidget(label_3, 2, 2, 1, 1);

        globalParameterTableView = new QTableView(centralWidget);
        globalParameterTableView->setObjectName(QString::fromUtf8("globalParameterTableView"));
        QSizePolicy sizePolicy2(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(globalParameterTableView->sizePolicy().hasHeightForWidth());
        globalParameterTableView->setSizePolicy(sizePolicy2);
        globalParameterTableView->setBaseSize(QSize(200, 500));
        globalParameterTableView->setEditTriggers(QAbstractItemView::NoEditTriggers);
        globalParameterTableView->horizontalHeader()->setVisible(false);
        globalParameterTableView->verticalHeader()->setVisible(false);

        gridLayout_2->addWidget(globalParameterTableView, 3, 0, 1, 1);

        colonyView = new QGraphicsView(centralWidget);
        colonyView->setObjectName(QString::fromUtf8("colonyView"));
        QSizePolicy sizePolicy3(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy3.setHorizontalStretch(0);
        sizePolicy3.setVerticalStretch(0);
        sizePolicy3.setHeightForWidth(colonyView->sizePolicy().hasHeightForWidth());
        colonyView->setSizePolicy(sizePolicy3);
        colonyView->setMaximumSize(QSize(700, 700));
        colonyView->setBaseSize(QSize(500, 500));
        colonyView->setFrameShape(QFrame::StyledPanel);
        colonyView->setFrameShadow(QFrame::Sunken);
        colonyView->setLineWidth(1);
        colonyView->setRenderHints(QPainter::Antialiasing|QPainter::TextAntialiasing);
        colonyView->setCacheMode(QGraphicsView::CacheBackground);
        colonyView->setViewportUpdateMode(QGraphicsView::SmartViewportUpdate);

        gridLayout_2->addWidget(colonyView, 3, 1, 1, 1);

        cellStateTableView = new QTableView(centralWidget);
        cellStateTableView->setObjectName(QString::fromUtf8("cellStateTableView"));
        sizePolicy2.setHeightForWidth(cellStateTableView->sizePolicy().hasHeightForWidth());
        cellStateTableView->setSizePolicy(sizePolicy2);
        cellStateTableView->setBaseSize(QSize(200, 500));

        gridLayout_2->addWidget(cellStateTableView, 3, 2, 1, 1);

        horizontalLayout_6 = new QHBoxLayout();
        horizontalLayout_6->setSpacing(6);
        horizontalLayout_6->setObjectName(QString::fromUtf8("horizontalLayout_6"));
        label = new QLabel(centralWidget);
        label->setObjectName(QString::fromUtf8("label"));
        QSizePolicy sizePolicy4(QSizePolicy::Minimum, QSizePolicy::Preferred);
        sizePolicy4.setHorizontalStretch(0);
        sizePolicy4.setVerticalStretch(0);
        sizePolicy4.setHeightForWidth(label->sizePolicy().hasHeightForWidth());
        label->setSizePolicy(sizePolicy4);
        label->setMinimumSize(QSize(50, 0));
        label->setMaximumSize(QSize(50, 16777215));
        label->setFont(font3);

        horizontalLayout_6->addWidget(label);

        zoomSlider = new QSlider(centralWidget);
        zoomSlider->setObjectName(QString::fromUtf8("zoomSlider"));
        QSizePolicy sizePolicy5(QSizePolicy::Expanding, QSizePolicy::Fixed);
        sizePolicy5.setHorizontalStretch(0);
        sizePolicy5.setVerticalStretch(0);
        sizePolicy5.setHeightForWidth(zoomSlider->sizePolicy().hasHeightForWidth());
        zoomSlider->setSizePolicy(sizePolicy5);
        zoomSlider->setMaximumSize(QSize(440, 16777215));
        zoomSlider->setValue(50);
        zoomSlider->setOrientation(Qt::Horizontal);

        horizontalLayout_6->addWidget(zoomSlider);


        gridLayout_2->addLayout(horizontalLayout_6, 4, 1, 1, 1);

        gridLayout_2->setColumnStretch(0, 1);
        gridLayout_2->setColumnStretch(1, 2);
        gridLayout_2->setColumnStretch(2, 1);

        gridLayout_3->addLayout(gridLayout_2, 0, 0, 1, 1);

        MainWindow->setCentralWidget(centralWidget);
        mainToolBar = new QToolBar(MainWindow);
        mainToolBar->setObjectName(QString::fromUtf8("mainToolBar"));
        MainWindow->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(MainWindow);
        statusBar->setObjectName(QString::fromUtf8("statusBar"));
        MainWindow->setStatusBar(statusBar);

        retranslateUi(MainWindow);
        QObject::connect(MainWindow, SIGNAL(updateSimulationProgressBar(int)), simulationProgressBar, SLOT(setValue(int)));
        QObject::connect(playButton, SIGNAL(clicked()), MainWindow, SLOT(playSimulation()));
        QObject::connect(MainWindow, SIGNAL(updateTrajectoryCounter(QString)), iTrajNum, SLOT(setText(QString)));
        QObject::connect(MainWindow, SIGNAL(updateNCells(int)), NCellsNum, SLOT(setNum(int)));
        QObject::connect(MainWindow, SIGNAL(updateTimeSimulation(QString)), timeNum, SLOT(setText(QString)));
        QObject::connect(MainWindow, SIGNAL(updateParameterSetCounter(QString)), iParameterSetNum, SLOT(setText(QString)));
        QObject::connect(zoomSlider, SIGNAL(valueChanged(int)), MainWindow, SLOT(on_zoomSlider_valueChanged(int)));
        QObject::connect(speciesColorComboBox, SIGNAL(activated(int)), MainWindow, SLOT(setSpeciesCellColorIndex(int)));
        QObject::connect(speciesColorMaxValueSpinBox, SIGNAL(valueChanged(double)), MainWindow, SLOT(setSpeciesCellColorMaxValue(double)));
        QObject::connect(speciesBackgroundColorMaxValueSpinBox, SIGNAL(valueChanged(double)), MainWindow, SLOT(setSpeciesBackgroundColorMaxValue(double)));
        QObject::connect(speciesBackgroundColorComboBox, SIGNAL(activated(int)), MainWindow, SLOT(setSpeciesBackgroundColorIndex(int)));
        QObject::connect(pauseButton, SIGNAL(clicked()), MainWindow, SLOT(pauseSimulation()));
        QObject::connect(pushButton, SIGNAL(clicked()), MainWindow, SLOT(openOpenGLViewer()));

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "Quorum Sensing Simulation", 0, QApplication::UnicodeUTF8));
        playButton->setText(QString());
        pauseButton->setText(QString());
        pushButton->setText(QApplication::translate("MainWindow", "OpenGL Viewer", 0, QApplication::UnicodeUTF8));
        NCells->setText(QApplication::translate("MainWindow", "NCells", 0, QApplication::UnicodeUTF8));
        NCellsNum->setText(QApplication::translate("MainWindow", "0", 0, QApplication::UnicodeUTF8));
        time->setText(QApplication::translate("MainWindow", "time", 0, QApplication::UnicodeUTF8));
        timeNum->setText(QApplication::translate("MainWindow", "0", 0, QApplication::UnicodeUTF8));
        iTraj->setText(QApplication::translate("MainWindow", "trajectory #", 0, QApplication::UnicodeUTF8));
        iTrajNum->setText(QApplication::translate("MainWindow", "0", 0, QApplication::UnicodeUTF8));
        iParameterSet->setText(QApplication::translate("MainWindow", "Parameter set #", 0, QApplication::UnicodeUTF8));
        iParameterSetNum->setText(QApplication::translate("MainWindow", "0", 0, QApplication::UnicodeUTF8));
        groupBox->setTitle(QApplication::translate("MainWindow", "Cell Color Code", 0, QApplication::UnicodeUTF8));
        groupBox_2->setTitle(QApplication::translate("MainWindow", "Background Color Code", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("MainWindow", "Global Parameters", 0, QApplication::UnicodeUTF8));
        label_3->setText(QApplication::translate("MainWindow", "Cell State", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        colonyView->setToolTip(QString());
#endif // QT_NO_TOOLTIP
        label->setText(QApplication::translate("MainWindow", "Zoom", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
