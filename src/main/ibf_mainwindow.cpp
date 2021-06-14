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
    close();
    //QApplication::quit();
}


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

    QString k_mer {"K-mer size: " + QString::number(k) + "\n"};
    QString threadsN {"Building threads: " + QString::number(t)+ "\n"};
    QString fragmentSize {"Size of reference fragments per bin: " + QString::number(f)+  "\n"};
    QString filterSize {"IBF file size in MegaBytes: " + QString::number(0)+ "\n"};
    QString refF {"Input reference file: " + input_reference+ "\n"};
    QString outF {"Output IBF file: " + output_file+"\n"};

    QString check {""};
    check.append(k_mer+ threadsN+ fragmentSize + filterSize+ refF+ outF);
    QMessageBox::information(this, "Build Interleaved Bloom Filter Arguments", check);

}
void IBF_mainwindow::third_party(){
    QDebugStream* test_1 = new QDebugStream(std::cerr, ui->text1);
    QDebugStream* test = new QDebugStream(std::cout, ui->text1);

}

void IBF_mainwindow::clearResults(){
    ui->text1->QTextEdit::clear();
}

void IBF_mainwindow::proBar(void){

  int numTasks = 100000;

        QProgressDialog progress("Task in progress...", "Cancel", 0, numTasks, this);
        progress.setWindowModality(Qt::WindowModal);

        for (int i = 0; i < numTasks; i++) {
            progress.setValue(i);

            if (progress.wasCanceled())
                break;
        }
        progress.setValue(numTasks);
}





/*void IBF_mainwindow::on_pushButton_5_clicked()
{
    this->hide();
    QWidget *parent = this->parentWidget();
    parent->show();
}*/

void IBF_mainwindow::on_pushButton_7_clicked()
{
    this->hide();
    QWidget *parent = this->parentWidget();
    parent->show();
}