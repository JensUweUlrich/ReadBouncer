#include "ibf_mainwindow.h"
#include "ui_ibf_mainwindow.h"

IBF_mainwindow::IBF_mainwindow(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::IBF_mainwindow)
{
    ui->setupUi(this);
}

IBF_mainwindow::~IBF_mainwindow()
{
    delete ui;
}

void IBF_mainwindow::on_pushButton_2_clicked()
{
    QApplication::quit();
}
