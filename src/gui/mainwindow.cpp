#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_clicked()
{
    QMessageBox msgBox;
    msgBox.setText("ReadBouncer Version 1!");
    msgBox.exec();
}


void MainWindow::on_pushButton_10_clicked()
{
    hide();
    ibf_multi_window = new IBF_mainwindow(this); // this: is mainwindow class
    ibf_multi_window -> show();
}


void MainWindow::on_pushButton_11_clicked()
{
    hide();
    classify_multi_window = new Classify_mainwindow(this); // this: is mainwindow class
    classify_multi_window -> show();
}

