#include "ibf_mainwindow.h"
#include "ui_ibf_mainwindow.h"
//#include "ibfbuild.hpp"


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

/*
 * void IBF_mainwindow::on_pushButton_3_clicked()
{

    //k = QInputDialog::getInt(this, "kmer-size", "Kmer size used for building the Interleaved Bloom Filter (default: 13)");
    //t = QInputDialog::getInt(this, "Threads", "Number of building threads");
    //f = QInputDialog::getInt(this, "Fragment-size", "Length of fragments from the reference that are put in one bin of the IBF (default: 100000)");
   // s = QInputDialog::getInt(this, "Filter-size", "IBF size in MB");

    QMessageBox::StandardButton ask;
     ask = QMessageBox::question(this, "Building IBF", "The IBF will be built, are you sure?",
                                   QMessageBox::Yes|QMessageBox::No);

     if (ask == QMessageBox::Yes) {
      // qDebug() << "Yes was clicked";
       //QApplication::quit();
         QMessageBox::StandardButton reply;
          reply = QMessageBox::question(this, "Verbose", "Show additional output as to what we are doing?",
                                        QMessageBox::Yes|QMessageBox::No);
          if (reply == QMessageBox::Yes) {
            qDebug() << "Yes was clicked";
            //QApplication::quit();

            /*  ibf_build_parser ibfbuild_parser{cli};
                buildIBF(ibfbuild_parser);
            */

            //struct ibf_build_parser ibf_build_test {output_file_Name, ref_file_Name, false, false, k, t, f, s, false};
            //buildIBF(ibf_build_test);
/*
          } else {
            qDebug() << "Yes was *not* clicked";
          }
     } else {
       qDebug() << "Yes was *not* clicked";
     }



}*/


void IBF_mainwindow::on_spinBox_valueChanged(int kmer_size)
{
    IBF_mainwindow::k = kmer_size;
}

void IBF_mainwindow::on_spinBox_2_valueChanged(int threads)
{
    IBF_mainwindow::t = threads;
}

void IBF_mainwindow::on_spinBox_3_valueChanged(int fragment_size)
{
    IBF_mainwindow::f = fragment_size;
}

/*void IBF_mainwindow::on_spinBox_4_valueChanged(int filer_size)
{
    IBF_mainwindow::s = filer_size;
}*/



void IBF_mainwindow::on_pushButton_clicked()
{
    QString inputRef = QFileDialog::getOpenFileName(this, "Reference sequence file (fasta format)");
    IBF_mainwindow::ref_file_Name =inputRef.toLocal8Bit().constData();
    QFile file(inputRef);
    input_reference = inputRef; // store filName
    if(!file.open(QIODevice::ReadOnly | QFile::Text)){
        // open a massage box for warning
        QMessageBox::warning(this, "Warning", "Cannot open file: " + file.errorString());
        return;
    }

    else {
        QTextStream in(&file);
        QString pointer = in.readAll();//copy the text
        this->refFile =pointer.toLocal8Bit().constData();
        //std::cout << refFile<<"\n";
    }
    setWindowTitle(inputRef);
    file.close();
}

void IBF_mainwindow::on_pushButton_4_clicked()
{
    //output_file = QInputDialog::getText(this, "Output file", "Output file of Interleaved Bloom Filter (required)");
    //std::string outputFile =output_file.toLocal8Bit().constData();
    //std::cout<<"The converted String is:"<<outputFile<<std::endl;

    // This section is for storing IBF in a file
    IBF_mainwindow::output_file = QFileDialog::getSaveFileName(this, "save the IBF as");
    output_file_Name =IBF_mainwindow::output_file.toLocal8Bit().constData();
    QFile out(IBF_mainwindow::output_file);

    if(!out.open(QFile::WriteOnly | QFile::Text)){
        // open a massage box for warning
        QMessageBox::warning(this, "Warning", "Cannot save file: " + out.errorString());
        return;
    }
    setWindowTitle(IBF_mainwindow::output_file);
    QTextStream out_1(&out);
    //QString text = ui->textEdit->toPlainText();
    // out << text;// wrte to this file.
    //file.close();
}

void IBF_mainwindow::on_pushButton_6_clicked()
{
    //std::cout<< "K-mer size: " << k << "\n";
    //std::cout<< "Number of threads: " << t << "\n";
    //std::cout<< "Fragment size: " << f << "\n";
    //std::cout<< "Filter size: " << s << "\n";
    //std::cout << "Name of ref: " <<  ref_file_Name << "\n";
    //std::cout << "And the sequence is : " <<  refFile << "\n";
    //std::cout<< "The Name of outpt File : " << output_file_Name << "\n";

    QString k_mer {"K-mer size: " + QString::number(k) + ", "};
    QString threadsN {"Number of Threads: " + QString::number(t)+ ", "};
    QString fragmentSize {"Fragment size: " + QString::number(f)+ ", "};
    //QString filterSize {"Filter size: " + QString::number(s)+ ", "};
    QString refF {"Name and location of reference: " + input_reference+ ", "};
    QString outF {"Name and location of output file: " + output_file+ "."};

    QString check {""};
    check.append(k_mer+ threadsN+ fragmentSize+ refF+ outF);
    QMessageBox::information(this, "The arguments for IBF", check);

}
