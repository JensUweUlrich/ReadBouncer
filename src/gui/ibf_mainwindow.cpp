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
    this->hide();
    QWidget *parent = this->parentWidget();
    parent->show();
}


void IBF_mainwindow::on_pushButton_3_clicked()
{
    close();

}

// output dir
void IBF_mainwindow::on_pushButton_4_clicked()
{
    QString output_log_file =QFileDialog::getExistingDirectory(this, "Select log and output directory");
    IBF_mainwindow::output_dir = output_log_file.toStdString();
    IBF_mainwindow::output_dir.make_preferred();

    if(!std::filesystem::is_directory(IBF_mainwindow::output_dir) || !std::filesystem::exists(IBF_mainwindow::output_dir)){

        std::filesystem::create_directory(IBF_mainwindow::output_dir);
    }
}

// reference files
void IBF_mainwindow::on_pushButton_5_clicked()
{
    QStringList referenceFilesQ = QFileDialog::getOpenFileNames(this, "Select reference files");

    for (QString reference : referenceFilesQ){

        QFile file(reference);

        if(!std::filesystem::exists(std::filesystem::path(reference.toStdString()))){

             QString msg = "The following reference file does not exist: ";
             msg.push_back(reference);
             QMessageBox::warning(this , "Warning" , msg);

          }else if(!file.open(QIODevice::ReadOnly)){

              QString msg = "Cannot open the selected file: ";
              msg.push_back(reference);
              QMessageBox::warning(this , "Warning" , msg);
              return;

           }else{

            IBF_mainwindow::reference_files.emplace_back((std::filesystem::path(reference.toStdString())).make_preferred());

        }
    }
}

// K-mer size
void IBF_mainwindow::on_spinBox_3_valueChanged(int arg1)
{
        IBF_mainwindow::k = arg1;

}


void IBF_mainwindow::on_spinBox_4_valueChanged(int arg1)
{
    IBF_mainwindow::threads = arg1;

    if (IBF_mainwindow::threads < 1){

        IBF_mainwindow::threads = 1;
    }
}


void IBF_mainwindow::on_spinBox_5_valueChanged(int arg1)
{
    IBF_mainwindow::filter_size = arg1;
}


void IBF_mainwindow::on_spinBox_6_valueChanged(int arg1)
{
    IBF_mainwindow::fragment_size = arg1;
}


void IBF_mainwindow::slot_control_std()
{
    QDebugStream* out_cerr = new QDebugStream(std::cerr, ui->output_window);
    QDebugStream* out_cout = new QDebugStream(std::cout, ui->output_window);

}

void IBF_mainwindow::clean_results()
{
    ui->output_window->QTextEdit::clear();
}


// Clear old results in the case of building many IBF's using the same dialog window
void IBF_mainwindow::on_pushButton_6_clicked()
{

    ui->output_window->QTextEdit::clear();

}

// plainTextEdit

void IBF_mainwindow::plainTextEditChange(QString &text){

    //ui->plainTextEdit->setPlainText(text);
}


