#include "mainwindow.h"
#include "ui_mainwindow.h"

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
    QString windowTitel = "Argument: Ibfbuild ";
    QString boxInformations = "Build Interleaved Bloom Filter with given references sequences: "
                              " -v, --verbose, "
                              "-o, --output-file <output-file>, "
                              "-i, --input-reference <input-reference>, "
                              " -k, --kmer-size <kmer-size>, "
                              "-t, --threads <threads>, "
                              "-f, --fragment-size <fragment-size>, "
                              "-s, --filter-size <filter-size>, ";
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
    ui->textEdit->setText(text);
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

/*void MainWindow::on_actionPrint_triggered()
{
    QPrinter printer; // def a printer
    printer.setPrinterName("Printer Name");
    QPrintDialog pDialog(&printer, this);

    if(pDialog.exec() == QDialog::Rejected){
        QMessageBox::warning(this, "Warning", "Cannot Access Printer");
        return;
    }
    ui->textEdit->print(&printer);
}*/

void MainWindow::on_actionNanoLive_triggered()
{
    QString windowTitel = "About NanoLive Tool";
    QString boxInformations = "C++ based for live classification of Nanopore reads (aka adaptive sequencing) on Windows without the need for GPUs."
                              " The Toolkit uses Oxford Nanopore's Read Until functionality to unblock reads that match to a given reference sequence database."
                              " The database is indexed as Interleaved Bloom Filter for fast classification.";

    QMessageBox::information(this, windowTitel, boxInformations);
}



void MainWindow::on_actionClassify_triggered()
{
    QString windowTitel = "Argument: Classify ";
    QString boxInformations = "Classify nanopore reads based on a given IBF file: "
                              "-v, --verbose, "
                              "-r, --read-file <read-file>, "
                              "-d, --depletion-file <ibf-file>, "
                              " -t, --target-file <ibf-file>, "
                              "-c, --classified-file <file>, "
                              "-u, --unclassified-file <file>, "
                              "-s, --significance <probability>, "
                              "-e, --error-rate <err>, "
                              " -p, --prefix-length <length>, "
                              "-n, --num-threads <threads>";
    QMessageBox::information(this, windowTitel, boxInformations);
}

void MainWindow::on_actionLive_deplete_triggered()
{
    QString windowTitel = "Argument: Live-deplete ";
    QString boxInformations = "Live classification and rejection of nanopore reads: "
                              "-v, --verbose, "
                              "-d, --device <device>, "
                              "-c, --host <host>, "
                              " -p, --port <port>, "
                              "-d, --depletion-file <ibf-file>, "
                              "-t, --target-file <ibf-file>, "
                              "-s, --significance <probability>, "
                              "-e, --error-rate <err>, "
                              " -w, --weights <weights>, ";
    QMessageBox::information(this, windowTitel, boxInformations);
}

void MainWindow::on_actionConnection_test_triggered()
{
    QString windowTitel = "Argument: Connection-test ";
    QString boxInformations = "Test connection to a working MinKNOW instance: "
                              "-v, --verbose, "
                              "-d, --device <device>, "
                              "-c, --host <host>, "
                              " -p, --port <port>, "
                              "-u, --unblock-all, ";
    QMessageBox::information(this, windowTitel, boxInformations);
}
