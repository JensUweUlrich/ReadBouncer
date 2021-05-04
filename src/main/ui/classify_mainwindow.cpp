#include "classify_mainwindow.h"
#include "ui_classify_mainwindow.h"

Classify_mainwindow::Classify_mainwindow(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Classify_mainwindow)
{
    ui->setupUi(this);
}

Classify_mainwindow::~Classify_mainwindow()
{
    delete ui;
}

void Classify_mainwindow::on_pushButton_clicked()
{
    QApplication::quit();
}
