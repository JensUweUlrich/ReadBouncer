#include "connection_test_mainwindow.h"
#include "ui_connection_test_mainwindow.h"

connection_test_mainwindow::connection_test_mainwindow(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::connection_test_mainwindow)
{
    ui->setupUi(this);
}

connection_test_mainwindow::~connection_test_mainwindow()
{
    delete ui;
}
