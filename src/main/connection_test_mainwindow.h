#ifndef CONNECTION_TEST_MAINWINDOW_H
#define CONNECTION_TEST_MAINWINDOW_H

#include <QDialog>
#include <QInputDialog>
#include <QMessageBox>
#include <QDebug>
#include <QFile>          // for opening the files
#include <QFileDialog>   // for pop up
#include <QTextEdit>
#include <QPointer>



// C++
#include <string>
#include <iostream>
#include <fstream>
#include <QtDebug>
#include <QObject>
#include <QPlainTextEdit>
#include <QMutex>
#include "QDebugStream.h"

namespace Ui {
class connection_test_mainwindow;
}

class connection_test_mainwindow : public QDialog
{
    Q_OBJECT

public:
    explicit connection_test_mainwindow(QWidget *parent = nullptr);
    ~connection_test_mainwindow();

private slots:

    void third_party();

    void on_pushButton_3_clicked();



    void on_spinBox_3_valueChanged(int port_);

    void on_lineEdit_3_textChanged(const QString &arg1);

    void on_lineEdit_2_textChanged(const QString &arg1);

    void on_pushButton_7_clicked();

    //void on_pushButton_2_clicked();

    void on_pushButton_4_clicked();

    void clearResults();

    void on_checkBox_toggled(bool checked);

    void on_pushButton_8_clicked();

private:
    Ui::connection_test_mainwindow *ui;
    bool unblock_all = false;

    int port {9501};

    QString host = "127.0.0.1";
    std::string host_name = "127.0.0.1";

    QString device{};
    std::string device_name{};

};

#endif // CONNECTION_TEST_MAINWINDOW_H
