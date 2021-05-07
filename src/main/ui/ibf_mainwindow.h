#ifndef IBF_MAINWINDOW_H
#define IBF_MAINWINDOW_H

#include <QDialog>
#include <QInputDialog>
#include <QMessageBox>
#include <QDebug>
#include <QFile>          // for opening the files
#include <QFileDialog>   // for pop up

// C++
#include <string>
#include <iostream>
#include <fstream>





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

    //void on_spinBox_4_valueChanged(int arg1);//s

    void on_pushButton_clicked();// open ref

    void on_pushButton_4_clicked();

    void on_pushButton_6_clicked();


private:
    Ui::IBF_mainwindow *ui;
    // to build:
    int k {13};
    int t {1};
    int f {10000};
    //int s {};
    QString input_reference = "";
    std::string refFile = ""; // make as std::string
    std::string ref_file_Name = ""; // for test .. not used!
    QString output_file = "";
    std::string output_file_Name = "";// for test


};

#endif // IBF_MAINWINDOW_H
