#include "live_deplete_mainwindow.h"
#include "ui_live_deplete_mainwindow.h"

live_deplete_mainwindow::live_deplete_mainwindow(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::live_deplete_mainwindow)
{
    ui->setupUi(this);
}

live_deplete_mainwindow::~live_deplete_mainwindow()
{
    delete ui;
}

void live_deplete_mainwindow::on_pushButton_clicked()
{
    QApplication::quit();
}
