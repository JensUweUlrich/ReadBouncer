#include "classify_mainwindow.h"
#include "ui_classify_mainwindow.h"

#include <iostream>
#include <experimental/filesystem>

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

void Classify_mainwindow::on_pushButton_3_clicked()
{
    close();
}


// back to mainWindow
void Classify_mainwindow::on_pushButton_2_clicked()
{
    this->hide();
    QWidget *parent = this->parentWidget();
    parent->show();
}


void Classify_mainwindow::on_pushButton_5_clicked()
{
    QStringList inputRef = QFileDialog::getOpenFileNames(this, "Reference sequence files in FASTA format");
    std::vector<std::experimental::filesystem::path> files;


}

