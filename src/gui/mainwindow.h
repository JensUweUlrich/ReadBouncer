#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QMessageBox>
#include <QPixmap>

#include "ibf_mainwindow.h"
//#include "ibfmainwindow.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
    IBF_mainwindow *ibf_multi_window;

private slots:
    void on_pushButton_clicked();

    void on_pushButton_10_clicked();

private:
    Ui::MainWindow *ui;
    //IBFMainWindow *ibf_main_window;
};

#endif // MAINWINDOW_H
