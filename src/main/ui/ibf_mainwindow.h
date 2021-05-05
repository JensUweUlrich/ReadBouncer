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

private:
    Ui::IBF_mainwindow *ui;
    QString input_reference = "";
    QString output_file = "";
};

#endif // IBF_MAINWINDOW_H
