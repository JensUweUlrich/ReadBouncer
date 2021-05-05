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

void IBF_mainwindow::on_pushButton_3_clicked()
{
    int k = QInputDialog::getInt(this, "kmer-size", "Kmer size used for building the Interleaved Bloom Filter (default: 13)");
    int t = QInputDialog::getInt(this, "Threads", "Number of building threads");
    int f = QInputDialog::getInt(this, "Fragment-size", "Length of fragments from the reference that are put in one bin of the IBF (default: 100000)");
    int s = QInputDialog::getInt(this, "Filter-size", "IBF size in MB");
    std::cout<<k+t+f+s<<std::endl;//4 fun


    QMessageBox::StandardButton reply;
     reply = QMessageBox::question(this, "Verbose", "Show additional output as to what we are doing?",
                                   QMessageBox::Yes|QMessageBox::No);
     if (reply == QMessageBox::Yes) {
       qDebug() << "Yes was clicked";
       QApplication::quit();
     } else {
       qDebug() << "Yes was *not* clicked";
     }


     QString inputRef = QFileDialog::getOpenFileName(this, "Reference sequence file (fasta format)");
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
         std::string refFile =pointer.toLocal8Bit().constData();
         std::cout << refFile<<"\n";
     }
     setWindowTitle(inputRef);
     file.close();


     //this section is for saving the IBF in a file
     output_file = QInputDialog::getText(this, "Output file", "Output file of Interleaved Bloom Filter (required)");
     std::string outputFile =output_file.toLocal8Bit().constData();
     std::cout<<"The converted String is:"<<outputFile<<std::endl;

     // This section is for storing IBF in a file
     QString output_file = QFileDialog::getSaveFileName(this, "save the IBF as");
     QFile out(output_file);

     if(!out.open(QFile::WriteOnly | QFile::Text)){
         // open a massage box for warning
         QMessageBox::warning(this, "Warning", "Cannot save file: " + file.errorString());
         return;
     }
     setWindowTitle(output_file);
     QTextStream out_1(&out);
     //QString text = ui->textEdit->toPlainText();
    // out << text;// wrte to this file.
     file.close();

}
