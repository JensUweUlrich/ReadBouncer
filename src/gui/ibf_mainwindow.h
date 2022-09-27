#ifndef IBF_MAINWINDOW_H
#define IBF_MAINWINDOW_H

#include <QDialog>
#include <QInputDialog>
#include <QFile>
#include <QFileDialog>
#include <QMessageBox>
#include <QTextEdit>
#include <fstream>
#include <QDebug>
#include <QDebugStream.h>
#include <filesystem>
#include <iostream>


namespace Ui {
class IBF_mainwindow;
}

class IBF_mainwindow : public QDialog
{
    Q_OBJECT

public:
    explicit IBF_mainwindow(QWidget *parent = nullptr);
    ~IBF_mainwindow();
    Ui::IBF_mainwindow *ui;


private slots:
    //void on_pushButton_clicked();

    void on_pushButton_2_clicked();

    void on_pushButton_3_clicked();

    void on_pushButton_4_clicked();

    void on_pushButton_5_clicked();

    void on_spinBox_3_valueChanged(int arg1);

    void on_spinBox_4_valueChanged(int arg1);

    void on_spinBox_5_valueChanged(int arg1);

    void on_spinBox_6_valueChanged(int arg1);

    void slot_control_std();

    void clean_results();

    void plainTextEditChange(QString&);
    void on_pushButton_6_clicked();

    void on_buildIBFbutton_clicked();

private:

    //IBF params
    std::vector<std::filesystem::path>
    reference_files{};

    std::filesystem::path output_dir;
    //std::string output_path, reference_file {};
    int k = 13;
    int threads = 1;
    int fragment_size = 100000;
    int filter_size = 0;

};

#endif // IBF_MAINWINDOW_H
