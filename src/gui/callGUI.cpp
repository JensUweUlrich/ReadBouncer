#include "mainwindow.h"
#include<QApplication>
#include "ibf_mainwindow.h"
#include "classify_mainwindow.h"
#include "connection_test_mainwindow.h"
#include "live_deplete_mainwindow.h"
#include "QDebugStream.h"


void caller (int argc, char *argv[]){


    QApplication a(argc, argv);
    MainWindow w;
    //w.QWidget::showMaximized();
    //w.resize(800, 120);
    w.show();
    a.exec();
}
