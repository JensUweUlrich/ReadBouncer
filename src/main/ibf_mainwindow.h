#ifndef IBF_MAINWINDOW_H
#define IBF_MAINWINDOW_H

#include <QDialog>
#include <QInputDialog>
#include <QMessageBox>
#include <QDebug>
#include <QFile>
#include <QFileDialog>
#include "QDebugStream.h"
#include <QTextEdit>
#include <QPlainTextEdit>
#include <QProgressBar>
#include <QObject>
#include <QtCore>
#include <QProgressDialog>

#include <string>
#include <iostream>
#include <fstream>
#include <chrono>

#include<QMoveEvent>





namespace Ui {
class IBF_mainwindow;
}

class IBF_mainwindow : public QDialog
{
    Q_OBJECT

public:
    explicit IBF_mainwindow(QWidget *parent = nullptr);

    ~IBF_mainwindow();

private slots:
    void on_pushButton_2_clicked();

    void on_pushButton_3_clicked();

    void on_spinBox_valueChanged(int arg1);//k

    void on_spinBox_2_valueChanged(int arg1);//t

    void on_spinBox_3_valueChanged(int arg1);//f

    void on_pushButton_clicked();

    void on_pushButton_4_clicked();

    void on_pushButton_6_clicked();

    void third_party();

    void clearResults();

    void proBar(void);


private:
    Ui::IBF_mainwindow *ui;
    int k {13};
    int t {1};
    int f {10000};

    QString input_reference = "";
    std::string refFile = "";// for test only!
    std::string ref_file_Name = "";

    QString output_file = "";
    std::string output_file_Name = "";


};

#endif // IBF_MAINWINDOW_H
