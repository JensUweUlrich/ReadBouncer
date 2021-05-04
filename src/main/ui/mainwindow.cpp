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
    // add pic
    QPixmap pix("C:/NanoLive_Qt-test/src/main/550285a-i2.jpg");
    ui -> label_pic ->setPixmap(pix);
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


void MainWindow::on_pushButton_2_clicked()
{
    // see the header file of mainwindow
    // here we are working with heap not stack
    /*
    IBF_mainwindow ibf_multi_window;
    ibf_multi_window.setModal(true);// but now we cannot access the mainWindow!
    ibf_multi_window.exec();
    */

    //hide(); // this will hide the first mainWindow
    ibf_multi_window = new IBF_mainwindow(this); // this: is mainwindow class
    ibf_multi_window -> show();
}

void MainWindow::on_pushButton_3_clicked()
{
    classify_multi_window = new Classify_mainwindow(this); // this: is mainwindow class
    classify_multi_window -> show();
}

void MainWindow::on_pushButton_4_clicked()
{
    live_deplete_multi_window = new live_deplete_mainwindow(this);
    live_deplete_multi_window -> show();
}

void MainWindow::on_pushButton_5_clicked()
{
    connection_test_multi_window = new connection_test_mainwindow(this);
    connection_test_multi_window -> show();
}
