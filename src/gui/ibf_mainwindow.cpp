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

// Build IBF
void IBF_mainwindow::on_pushButton_clicked()
{
    if (IBF_mainwindow::k < 10){

        QMessageBox::warning(this , trUtf8("Warning") ,
                  "<P><FONT COLOR='#ffffff'>The selcted k-mer size is smaller than 10, we will use the default value 13</FONT>");
        IBF_mainwindow::k = 13;
    }

    std::cout<< IBF_mainwindow::output_path << '\n';
    std::cout<< IBF_mainwindow::reference_file << '\n';
    std::cout<< IBF_mainwindow::k << '\n';
    std::cout<< IBF_mainwindow::threads << '\n';
    std::cout<< IBF_mainwindow::fragment_size << '\n';
    std::cout<< IBF_mainwindow::filter_size << '\n';

    qDebug() << "C++ Style Debug Message";

}


void IBF_mainwindow::on_pushButton_2_clicked()
{
    this->hide();
    QWidget *parent = this->parentWidget();
    parent->show();
}


void IBF_mainwindow::on_pushButton_3_clicked()
{
    close();

}

void IBF_mainwindow::on_pushButton_4_clicked()
{
    QString out = QFileDialog::getSaveFileName(this, "Save the IBF as");
    IBF_mainwindow::output_path = out.toLocal8Bit().constData();

    // Create file (like ofstream)
    QFile out_file(out);

    if(!out_file.open(QFile::WriteOnly | QFile::Text)){

        QMessageBox::warning(this , trUtf8("Warning") ,
                             "<P><FONT COLOR='#ffffff'>Cannot save the selected file</FONT>");

    }
}


void IBF_mainwindow::on_pushButton_5_clicked()
{
    QString inputRef = QFileDialog::getOpenFileName(this, "Reference sequence file in FASTA format");
    IBF_mainwindow::reference_file =inputRef.toLocal8Bit().constData();
    QFile file(inputRef);
     if(!file.open(QIODevice::ReadOnly | QFile::Text)){
         QMessageBox::warning(this , trUtf8("Warning") ,
                   "<P><FONT COLOR='#ffffff'>Cannot open the selected file</FONT>");
         return;
      }
        setWindowTitle(inputRef);
        file.close();
}

// K-mer size
void IBF_mainwindow::on_spinBox_3_valueChanged(int arg1)
{
        IBF_mainwindow::k = arg1;

}


void IBF_mainwindow::on_spinBox_4_valueChanged(int arg1)
{
    IBF_mainwindow::threads = arg1;
}


void IBF_mainwindow::on_spinBox_5_valueChanged(int arg1)
{
    IBF_mainwindow::filter_size = arg1;
}


void IBF_mainwindow::on_spinBox_6_valueChanged(int arg1)
{
    IBF_mainwindow::fragment_size = arg1;
}

