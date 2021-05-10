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

void Classify_mainwindow::on_pushButton_3_clicked()
{
    close();
    //QApplication::quit();
}



void Classify_mainwindow::on_spinBox_4_valueChanged(int threadsNr)
{
    Classify_mainwindow::threads = threadsNr;
}

void Classify_mainwindow::on_doubleSpinBox_valueChanged(double kmerSig)//kmer_significance
{
    Classify_mainwindow::kmer_significance = kmerSig;
}

void Classify_mainwindow::on_doubleSpinBox_2_valueChanged(double errorRa)
{
    Classify_mainwindow::error_rate = errorRa;
}

void Classify_mainwindow::on_spinBox_3_valueChanged(int prefixLen)
{
    Classify_mainwindow::preLen = prefixLen;
}


void Classify_mainwindow::on_pushButton_5_clicked()//Interleaved Bloom Filter file with depletion references
{

    ibf_deplete_file= QFileDialog::getOpenFileName(this, "Interleaved Bloom Filter file with depletion reference");
    Classify_mainwindow::ibf_deplete_file_name = ibf_deplete_file.toLocal8Bit().constData();
    QFile file(ibf_deplete_file);
    if(!file.open(QIODevice::ReadOnly | QFile::Text)){
        QMessageBox::warning(this, "Warning", "Cannot open file: " + file.errorString());
        return;
    }
    setWindowTitle(ibf_deplete_file);
    file.close();
}


void Classify_mainwindow::on_pushButton_8_clicked()//Interleaved Bloom Filter file with target references
{
    ibf_target_file= QFileDialog::getOpenFileName(this, "Interleaved Bloom Filter file with target references");
    Classify_mainwindow::ibf_target_file_name = ibf_target_file.toLocal8Bit().constData();
    QFile file(ibf_target_file);
    if(!file.open(QIODevice::ReadOnly | QFile::Text)){
        QMessageBox::warning(this, "Warning", "Cannot open file: " + file.errorString());
        return;
    }
    setWindowTitle(ibf_target_file);
    file.close();
}

void Classify_mainwindow::on_pushButton_9_clicked()
{
    classified_file= QFileDialog::getOpenFileName(this, "File with classified reads in FASTA format");
    Classify_mainwindow::classified_file_name = classified_file.toLocal8Bit().constData();
    QFile file(classified_file);
    if(!file.open(QIODevice::ReadOnly | QFile::Text)){
        QMessageBox::warning(this, "Warning", "Cannot open file: " + file.errorString());
        return;
    }
    setWindowTitle(classified_file);
    file.close();
}

void Classify_mainwindow::on_pushButton_10_clicked()
{
    unclassified_file= QFileDialog::getOpenFileName(this, "File with unclassified reads in FASTA format");
    Classify_mainwindow::unclassified_file_name = unclassified_file.toLocal8Bit().constData();
    QFile file(unclassified_file);
    if(!file.open(QIODevice::ReadOnly | QFile::Text)){
        QMessageBox::warning(this, "Warning", "Cannot open file: " + file.errorString());
        return;
    }
    setWindowTitle(unclassified_file);
    file.close();
}

void Classify_mainwindow::on_pushButton_6_clicked()
{
    read_file= QFileDialog::getOpenFileName(this, "File with reads to classify in FASTA or FASTQ format");
    Classify_mainwindow::read_file_name = read_file.toLocal8Bit().constData();
    QFile file(read_file);
    if(!file.open(QIODevice::ReadOnly | QFile::Text)){
        QMessageBox::warning(this, "Warning", "Cannot open file: " + file.errorString());
        return;
    }
    setWindowTitle(read_file);
    file.close();
}

void Classify_mainwindow::on_pushButton_7_clicked()// Verbose
{

     QString readF {"Input read file: " + read_file+"\n"};
     QString depF  {"Depletion IBF file: " +ibf_deplete_file+"\n"};
     QString targetF {"Target IBF file: " + ibf_target_file+"\n"};
     QString classF {"Classified reads file: " + classified_file+"\n"};
     QString unclassF {"Unclassified reads file: " + unclassified_file+"\n"};
     QString sigLe {"Significance level for confidence interval: " + QString::number(kmer_significance)+ "\n"};
     QString erRa {"Expected sequencing error rate: " + QString::number(error_rate)+ "\n"};
     QString rePrefixLe {"Length of read prefix used for classification: " + QString::number(preLen)+ "\n"};
     QString threadsN {"Building threads: " + QString::number(threads)+ "\n"};

    QString check {""};
    check.append(readF + depF + targetF + classF + unclassF + sigLe + erRa + rePrefixLe + threadsN);
    QMessageBox::information(this, "Classify Reads Arguments", check);

}

void Classify_mainwindow::third_party(){
    QDebugStream* test_1 = new QDebugStream(std::cerr, ui->text1);
    QDebugStream* test = new QDebugStream(std::cout, ui->text1);
}

// this part is only for private use and tests
/*void test1(){
std::cout<<"text1"<<std::endl;
std::cout<<"text2"<<std::endl;

}
void Classify_mainwindow::on_pushButton_11_clicked()
{
    third_party();
    test1();
}*/


