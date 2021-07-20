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
    close();
}

void live_deplete_mainwindow::third_party(){
    //QDebugStream* test_1 = new QDebugStream(std::cerr, ui->text1);
    //QDebugStream* test = new QDebugStream(std::cout, ui->text1);
}

void live_deplete_mainwindow::on_pushButton_3_clicked()
{
    close();
}

void live_deplete_mainwindow::on_doubleSpinBox_valueChanged(double kmerSig)
{
    live_deplete_mainwindow::kmer_significance = kmerSig;
}

void live_deplete_mainwindow::on_doubleSpinBox_2_valueChanged(double errorRa)
{
    live_deplete_mainwindow::error_rate = errorRa;
}

void live_deplete_mainwindow::on_spinBox_3_valueChanged(int port_)
{
    live_deplete_mainwindow::port = port_;
}

void live_deplete_mainwindow::on_lineEdit_3_textChanged(const QString &host_)
{
   live_deplete_mainwindow::host = host_;
    host_name = live_deplete_mainwindow::host.toLocal8Bit().constData();

}

void live_deplete_mainwindow::on_lineEdit_textChanged(const QString &weights_)
{
    live_deplete_mainwindow::weights = weights_;
    weights_name = live_deplete_mainwindow::weights.toLocal8Bit().constData();
}

void live_deplete_mainwindow::on_lineEdit_2_textChanged(const QString &device_)
{
     live_deplete_mainwindow::device = device_;
     device_name = live_deplete_mainwindow::device.toLocal8Bit().constData();
}


void live_deplete_mainwindow::on_pushButton_5_clicked()
{
    ibf_deplete_file= QFileDialog::getOpenFileName(this, "Interleaved Bloom Filter file with depletion reference");
    live_deplete_mainwindow::ibf_deplete_file_name = ibf_deplete_file.toLocal8Bit().constData();
    QFile file(ibf_deplete_file);
    if(!file.open(QIODevice::ReadOnly | QFile::Text)){
        QMessageBox::warning(this, "Warning", "Cannot open file: " + file.errorString());
        return;
    }
    setWindowTitle(ibf_deplete_file);
    file.close();
}

void live_deplete_mainwindow::on_pushButton_8_clicked()
{
    ibf_target_file= QFileDialog::getOpenFileName(this, "Interleaved Bloom Filter file with target references");
    live_deplete_mainwindow::ibf_target_file_name = ibf_target_file.toLocal8Bit().constData();
    QFile file(ibf_target_file);
    if(!file.open(QIODevice::ReadOnly | QFile::Text)){
        QMessageBox::warning(this, "Warning", "Cannot open file: " + file.errorString());
        return;
    }
    setWindowTitle(ibf_target_file);
    file.close();
}


void live_deplete_mainwindow::on_pushButton_7_clicked()
{
    QString host_v {"Host IP address: " +host+"\n"};
    QString port_v {"MinKNOW communication port: " + QString::number(port) +"\n"};
    QString device_v {"Device or Flowcell name: " + device+"\n"};
    QString depF  {"Depletion IBF file: " +ibf_deplete_file+"\n"};
    QString targetF {"Target IBF file: " + ibf_target_file+"\n"};
    QString sigLe {"Significance level for confidence interval: " + QString::number(kmer_significance)+ "\n"};
    QString erRa {"Expected sequencing error rate: " + QString::number(error_rate)+ "\n"};
    QString weights_v {"Deep Nano Weights: " +weights+ "\n"};
    QString classThreads {"Number of threads (base calling): " +QString::number(basecall_threads)+ "\n"};
    QString baseCallThreads {"Number of threads (classification):  " +QString::number(classify_threads)+ "\n"};

   QString check {""};
   check.append(host_v + port_v + device_v + depF + targetF + sigLe + erRa + weights_v +classThreads + baseCallThreads);
   QMessageBox::information(this, "Live Classification Arguments", check);
}



void live_deplete_mainwindow::clearResults(){
    //ui->text1->QTextEdit::clear();
}

void live_deplete_mainwindow::on_pushButton_6_clicked()
{
    this->hide();
    QWidget *parent = this->parentWidget();
    parent->show();
}

void live_deplete_mainwindow::on_spinBox_4_valueChanged(int classificationThreads)
{
    live_deplete_mainwindow::classify_threads = classificationThreads;
}

void live_deplete_mainwindow::on_spinBox_5_valueChanged(int basscallingThreads)
{
    live_deplete_mainwindow::basecall_threads = basscallingThreads;
}

