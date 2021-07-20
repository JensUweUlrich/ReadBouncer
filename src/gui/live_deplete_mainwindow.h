#ifndef LIVE_DEPLETE_MAINWINDOW_H
#define LIVE_DEPLETE_MAINWINDOW_H

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
class live_deplete_mainwindow;
}

class live_deplete_mainwindow : public QDialog
{
    Q_OBJECT

public:
    explicit live_deplete_mainwindow(QWidget *parent = nullptr);
    ~live_deplete_mainwindow();

private slots:
    void on_pushButton_clicked();
    void third_party();

    void on_pushButton_3_clicked();

    void on_doubleSpinBox_valueChanged(double arg1);

    void on_doubleSpinBox_2_valueChanged(double arg1);

    void on_spinBox_3_valueChanged(int arg1);

    void on_lineEdit_3_textChanged(const QString &host_);

    void on_lineEdit_textChanged(const QString &weights_);

    void on_lineEdit_2_textChanged(const QString &device_);

    void on_pushButton_5_clicked();

    void on_pushButton_8_clicked();

    void on_pushButton_7_clicked();

    void on_pushButton_4_clicked();

    void clearResults();

    void on_pushButton_6_clicked();

    void on_spinBox_4_valueChanged(int arg1);

    void on_spinBox_5_valueChanged(int arg1);

    void on_pushButton_9_clicked();

private:
    Ui::live_deplete_mainwindow *ui;

    QString host = "127.0.0.1";
    std::string host_name = "127.0.0.1";

    QString device{};
    std::string device_name{};

    QString ibf_deplete_file{""};
    std::string ibf_deplete_file_name{""};

    QString ibf_target_file{""};
    std::string ibf_target_file_name{""};

    QString weights = "48";
    std::string weights_name = "48";
    int basecall_threads {1};
    int classify_threads {1};

    int port {9501};
    double kmer_significance{0.95};
    double error_rate {0.1};



};

#endif // LIVE_DEPLETE_MAINWINDOW_H
