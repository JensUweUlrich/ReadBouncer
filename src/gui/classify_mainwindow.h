#ifndef CLASSIFY_MAINWINDOW_H
#define CLASSIFY_MAINWINDOW_H


#include <QDialog>
#include <QInputDialog>
#include <QMessageBox>
#include <QDebug>
#include <QFile>
#include <QFileDialog>
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
class Classify_mainwindow;
}

class Classify_mainwindow : public QDialog
{
    Q_OBJECT

public:
    explicit Classify_mainwindow(QWidget *parent = nullptr);



    ~Classify_mainwindow();

private slots:

    void on_pushButton_3_clicked();

    void on_spinBox_4_valueChanged(int arg1);

    void on_doubleSpinBox_valueChanged(double arg1);

    void on_doubleSpinBox_2_valueChanged(double arg1);

    void on_spinBox_3_valueChanged(int arg1);

    void on_pushButton_5_clicked();

    void on_pushButton_8_clicked();

    //void on_pushButton_9_clicked();

    //void on_pushButton_10_clicked();

    void on_pushButton_6_clicked();

    void on_pushButton_7_clicked();

    void on_pushButton_4_clicked();

    //void on_pushButton_11_clicked();// to test third party

    void third_party();

    void on_lineEdit_3_textEdited(const QString &arg1);

    void on_lineEdit_4_textEdited(const QString &arg1);

    void clearResults();


    void on_pushButton_9_clicked();

private:

    Ui::Classify_mainwindow *ui;
    // verbose
    QString read_file{""};
    std::string read_file_name{""};

    QString ibf_deplete_file{""};
    std::string ibf_deplete_file_name{""};

    QString ibf_target_file{""};
    std::string ibf_target_file_name{""};

    QString classified_file{""};
    std::string classified_file_name{""};

    QString unclassified_file{""};
    std::string unclassified_file_name{""};

    std::string classified {""};
    std::string unclassified {""};

    double kmer_significance{0.95};
    double error_rate {0.1};
    int threads {1};
    int preLen {360};//prefix-length


};

#endif // CLASSIFY_MAINWINDOW_H


