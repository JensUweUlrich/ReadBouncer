#include "mainwindow.h"
#include "ui_mainwindow.h"


// not anymore here --> in header
/*
#include "ibf_mainwindow.h"
#include "classify_mainwindow.h"
#include "live_deplete_mainwindow.h"
#include "connection_test_mainwindow.h"
*/

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{

    ui->setupUi(this);

}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_clicked()
{
    QApplication::quit();
}

void MainWindow::on_action_triggered()
{
    QString windowTitel = "Argument: Ibfbuild: ";
    QString boxInformations = "Build Interleaved Bloom Filter with given references sequences, for this option you "
                              "need an input-reference sequence, an output-file for storing the IBF, kmer-size is default 13, number "
                              "of threads and fragment-size. ";
    QMessageBox::information(this, windowTitel, boxInformations);
}

void MainWindow::on_actionOpen_triggered()
{
    QString fileName = QFileDialog::getOpenFileName(this, "Open the file");
    QFile file(fileName);
    currentFile = fileName; // store filName
    if(!file.open(QIODevice::ReadOnly | QFile::Text)){
        // open a massage box for warning
        QMessageBox::warning(this, "Warning", "Cannot open file: " + file.errorString());
        return;
    }
    setWindowTitle(fileName);
    QTextStream in(&file);
    QString text = in.readAll();//copy the text
    // put the text in text widget
    //ui->textEdit->setText(text);
    file.close();

}

/*void MainWindow::on_actionExit_triggered()
{
    QApplication::quit();
}*/

void MainWindow::on_actionExit_triggered()
{
    QApplication::quit();
}


void MainWindow::on_actionNanoLive_triggered()
{
    QString windowTitel = "About ReadBouncer Tool";

    QString boxInformations = "C++ based for live classification of Nanopore reads (aka adaptive sequencing) on Windows without the need for GPUs."
                              " The Toolkit uses Oxford Nanopore's Read Until functionality to unblock reads that match to a given reference sequence database."
                              " The database is indexed as Interleaved Bloom Filter for fast classification.";
    QMessageBox::information(this, windowTitel, boxInformations);
    //QMessageBox::setStyleSheet("background-color:red");
}



void MainWindow::on_actionClassify_triggered()
{
    QString windowTitel = "Argument: Classify ";
    QString boxInformations = "Classify nanopore reads based on a given IBF file, with:  "
                              "read-file, depletion-file, target-file, classified-file, unclassified-file"
                              "significance, error-rate, prefix-length and number of threads.";
    QMessageBox::information(this, windowTitel, boxInformations);
}

void MainWindow::on_actionLive_deplete_triggered()
{
    QString windowTitel = "Argument: Live-deplete ";
    QString boxInformations = "Live classification and rejection of nanopore reads, with:"
                              "device-name, host, port, depletion-file, target-file, significance"
                              " error-rate, and weights";
    QMessageBox::information(this, windowTitel, boxInformations);
}

void MainWindow::on_actionConnection_test_triggered()
{
    QString windowTitel = "Argument: Connection-test ";
    QString boxInformations = "Test connection to a working MinKNOW instance ";
    QMessageBox::information(this, windowTitel, boxInformations);
}


void MainWindow::on_pushButton_2_clicked()
{

    hide();
    ibf_multi_window = new IBF_mainwindow(this); // this: is mainwindow class
    ibf_multi_window -> show();
}

void MainWindow::on_pushButton_3_clicked()
{
    hide();
    classify_multi_window = new Classify_mainwindow(this); // this: is mainwindow class
    classify_multi_window -> show();
}

void MainWindow::on_pushButton_4_clicked()
{
    hide();
    live_deplete_multi_window = new live_deplete_mainwindow(this);
    live_deplete_multi_window -> show();
}

void MainWindow::on_pushButton_5_clicked()
{
    hide();
    connection_test_multi_window = new connection_test_mainwindow(this);
    connection_test_multi_window -> show();
}


