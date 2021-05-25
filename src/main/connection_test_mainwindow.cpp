#include "connection_test_mainwindow.h"
#include "ui_connection_test_mainwindow.h"

connection_test_mainwindow::connection_test_mainwindow(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::connection_test_mainwindow)
{
    ui->setupUi(this);
}

connection_test_mainwindow::~connection_test_mainwindow()
{
    delete ui;
}

void connection_test_mainwindow::third_party(){
    QDebugStream* test_1 = new QDebugStream(std::cerr, ui->text1);
    QDebugStream* test = new QDebugStream(std::cout, ui->text1);
}

void connection_test_mainwindow::on_pushButton_3_clicked()
{
    close();
}


void connection_test_mainwindow::on_spinBox_3_valueChanged(int port_)
{
    connection_test_mainwindow::port = port_;
}

void connection_test_mainwindow::on_lineEdit_3_textChanged(const QString &host_)
{
    connection_test_mainwindow::host = host_;
    host_name = connection_test_mainwindow::host.toLocal8Bit().constData();
}

void connection_test_mainwindow::on_lineEdit_2_textChanged(const QString &device_)
{
    connection_test_mainwindow::device = device_;
    device_name = connection_test_mainwindow::device.toLocal8Bit().constData();
}

void connection_test_mainwindow::on_pushButton_7_clicked()
{
     QString host_v {"Host IP address: " +host+"\n"};
     QString port_v {"MinKNOW communication port: " + QString::number(port) +"\n"};
     QString device_v {"Device or Flowcell name: " + device+"\n"};
     QString unBlockAll  {"Unblock all live reads: " +QString::number(unblock_all)+"\n"};
     QString check {""};
     check.append(host_v + port_v + device_v + unBlockAll);
     QMessageBox::information(this, "Connection Test Arguments", check);


}

/*void connection_test_mainwindow::on_pushButton_2_clicked()
{
    connection_test_mainwindow::unblock_all = true;
}*/


void connection_test_mainwindow::on_checkBox_toggled(bool checked)
{
    connection_test_mainwindow::unblock_all = checked;
}

void connection_test_mainwindow::clearResults(){
    ui->text1->QTextEdit::clear();
}
