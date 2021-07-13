#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QFile>
#include <QFileDialog>
#include <QTextStream>
#include <QMessageBox>
#include <QString>
#include<QtPrintSupport/QPrinter>
#include<QtPrintSupport/QPrintDialog>
#include<QtWidgets>
#include <QPalette>
#include <QColor>


// include multiWindows
#include "ibf_mainwindow.h"
#include "classify_mainwindow.h"
#include "live_deplete_mainwindow.h"
#include "connection_test_mainwindow.h"


namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_pushButton_clicked();

    void on_action_triggered();

    void on_actionOpen_triggered();

    void on_actionExit_triggered();

    void on_actionNanoLive_triggered();

    void on_actionClassify_triggered();

    void on_actionLive_deplete_triggered();

    void on_actionConnection_test_triggered();

    void on_pushButton_2_clicked();

    void on_pushButton_3_clicked();

    void on_pushButton_4_clicked();

    void on_pushButton_5_clicked();


private:
    Ui::MainWindow *ui;
     QString currentFile = " ";
    //std::string currentFile = "";
     IBF_mainwindow *ibf_multi_window;
     Classify_mainwindow *classify_multi_window;
     live_deplete_mainwindow *live_deplete_multi_window;
     connection_test_mainwindow *connection_test_multi_window;

};

#endif // MAINWINDOW_H
