#ifndef CLASSIFY_MAINWINDOW_H
#define CLASSIFY_MAINWINDOW_H

#include <QDialog>
#include <QInputDialog>
#include <QFile>
#include <QFileDialog>
#include <QMessageBox>
#include <QTextEdit>
#include <fstream>
#include <QDebug>
#include <QDebugStream.h>
//#include <experimental/filesystem>

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <filesystem>

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

    void check_params ();
    void on_pushButton_3_clicked();

    void on_pushButton_2_clicked();

    void on_pushButton_5_clicked();

    void on_doubleSpinBox_2_valueChanged(double arg1);

    void on_spinBox_9_valueChanged(int arg1);

    void on_spinBox_10_valueChanged(int arg1);

    void on_spinBox_3_valueChanged(int arg1);

    void on_spinBox_4_valueChanged(int arg1);

    void on_spinBox_6_valueChanged(int arg1);

    void on_classifyButton_clicked();

    void on_pushButton_4_clicked();

    void on_pushButton_8_clicked();

    void on_pushButton_7_clicked();

    void slot_control_std();


    void on_pushButton_6_clicked();

private:
    Ui::Classify_mainwindow *ui;
    // classify related params
    std::vector<std::filesystem::path>
    target_files{};
    std::vector<std::filesystem::path>
    deplete_files{};
    std::vector<std::filesystem::path>
    read_files{};

    std::filesystem::path output_dir;

    double error_rate = 0.1;
    int chunk_length = 360;
    int max_chunks = 1;

    // Building IBF related params
    int k = 13;
    int threads = 1;
    int fragment_size = 100000;

};

#endif // CLASSIFY_MAINWINDOW_H
