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
}


// back to mainWindow
void Classify_mainwindow::on_pushButton_2_clicked()
{
    this->hide();
    QWidget *parent = this->parentWidget();
    parent->show();
}


// target_files
void Classify_mainwindow::on_pushButton_5_clicked()
{
    QStringList targetFilesQ = QFileDialog::getOpenFileNames(this, "Select target files");

    for (QString target : targetFilesQ){

        QFile file(target);

        if(!std::filesystem::exists(std::filesystem::path(target.toStdString()))){

             QString msg = "The following target file does not exist: ";
             msg.push_back(target);
             QMessageBox::warning(this , "Warning" , msg);

          }else if(!file.open(QIODevice::ReadOnly)){

              QString msg = "Cannot open the selected file: ";
              msg.push_back(target);
              QMessageBox::warning(this , "Warning" , msg);
              return;

           }else{

            Classify_mainwindow::target_files.emplace_back((std::filesystem::path(target.toStdString())).make_preferred());

        }
    }

}



void Classify_mainwindow::on_doubleSpinBox_2_valueChanged(double arg1)
{
    Classify_mainwindow::error_rate = arg1;
}


void Classify_mainwindow::on_spinBox_9_valueChanged(int arg1)
{
    Classify_mainwindow::chunk_length = arg1;

}


void Classify_mainwindow::on_spinBox_10_valueChanged(int arg1)
{
    Classify_mainwindow::max_chunks = arg1;
}


void Classify_mainwindow::on_spinBox_3_valueChanged(int arg1)
{
    Classify_mainwindow::k = arg1;

}


void Classify_mainwindow::on_spinBox_4_valueChanged(int arg1)
{
    Classify_mainwindow::threads = arg1;
}


void Classify_mainwindow::on_spinBox_6_valueChanged(int arg1)
{
    Classify_mainwindow::fragment_size = arg1;
}


// deplete_files
void Classify_mainwindow::on_pushButton_4_clicked()
{
    QStringList depleteFilesQ = QFileDialog::getOpenFileNames(this, "Select deplete files");

    for (QString deplete : depleteFilesQ){

        QFile file(deplete);

        if(!std::filesystem::exists(std::filesystem::path(deplete.toStdString()))){

             QString msg = "The following deplete file does not exist: ";
             msg.push_back(deplete);
             QMessageBox::warning(this , "Warning" , msg);

          }else if(!file.open(QIODevice::ReadOnly)){

              QString msg = "Cannot open the selected file: ";
              msg.push_back(deplete);
              QMessageBox::warning(this , "Warning" , msg);
              return;

           }else{

            Classify_mainwindow::deplete_files.emplace_back((std::filesystem::path(deplete.toStdString())).make_preferred());

        }
    }
}

// reads
void Classify_mainwindow::on_pushButton_8_clicked()
{
    QStringList readFilesQ = QFileDialog::getOpenFileNames(this, "Select read files");

    for (QString read : readFilesQ){

        QFile file(read);

        if(!std::filesystem::exists(std::filesystem::path(read.toStdString()))){

             QString msg = "The following read file does not exist: ";
             msg.push_back(read);
             QMessageBox::warning(this , "Warning" , msg);

          }else if(!file.open(QIODevice::ReadOnly)){

              QString msg = "Cannot open the selected file: ";
              msg.push_back(read);
              QMessageBox::warning(this , "Warning" , msg);
              return;

           }else{

            Classify_mainwindow::read_files.emplace_back((std::filesystem::path(read.toStdString())).make_preferred());

        }
    }
}

// output and log directory (same directory)
void Classify_mainwindow::on_pushButton_7_clicked()
{
    QString output_log_file =QFileDialog::getExistingDirectory(this, "Select log and output directory");
    Classify_mainwindow::output_dir = output_log_file.toStdString();
    Classify_mainwindow::output_dir.make_preferred();

    if(!std::filesystem::is_directory(Classify_mainwindow::output_dir) || !std::filesystem::exists(Classify_mainwindow::output_dir)){

        std::filesystem::create_directory(Classify_mainwindow::output_dir);
    }
}

void Classify_mainwindow::slot_control_std()
{
    QDebugStream* out_cerr = new QDebugStream(std::cerr, ui->output_window);
    QDebugStream* out_cout = new QDebugStream(std::cout, ui->output_window);

}

// only for developing tests (GTest needed)!

void Classify_mainwindow::check_params(){

    qDebug()<< "Target files: " ;
    for (auto target : Classify_mainwindow::target_files){

        qDebug()<< QString::fromStdString(target.string()) ;
    }

    qDebug()<< "Deplete files: " ;
    for (auto deplete : Classify_mainwindow::deplete_files){

        qDebug()<< QString::fromStdString(deplete.string()) ;
    }

    qDebug()<< "Read files: " ;
    for (auto read : Classify_mainwindow::read_files){

        qDebug()<< QString::fromStdString(read.string()) ;
    }

    qDebug()<< "Output dir: ";
    qDebug()<< QString::fromStdString(Classify_mainwindow::output_dir.string());

    qDebug()<< "error_rate: ";
    qDebug()<< Classify_mainwindow::error_rate;

    qDebug()<< "chunk_length: ";
    qDebug()<< Classify_mainwindow::chunk_length;

    qDebug()<< "max_chunks: ";
    qDebug()<< Classify_mainwindow::max_chunks;

    qDebug()<< "k: ";
    qDebug()<< Classify_mainwindow::k;

    qDebug()<< "threads: ";
    qDebug()<< Classify_mainwindow::threads;

    qDebug()<< "fragment_size: ";
    qDebug()<< Classify_mainwindow::fragment_size;
}






void Classify_mainwindow::on_pushButton_6_clicked()
{
    ui->output_window->QTextEdit::clear();
}

