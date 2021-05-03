#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QFile>          // for opening the files
#include <QFileDialog>    // for pop up
#include <QTextStream>    // for reading text from file
#include <QMessageBox>    // for opening a message box i.e type of user errors
#include <QString>
#include<QtPrintSupport/QPrinter>  // to use the funcionality of printer
#include<QtPrintSupport/QPrintDialog>  // to open a dialog box and let the user choose things from it
#include<QtWidgets>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_pushButton_clicked();

    void on_action_triggered();

    void on_actionOpen_triggered();

    void on_actionExit_triggered();

    void on_actionNanoLive_triggered();

    void on_actionClassify_triggered();

    void on_actionLive_deplete_triggered();

    void on_actionConnection_test_triggered();


private:
    Ui::MainWindow *ui;
     QString currentFile = " ";
    //std::string currentFile = "";
};

#endif // MAINWINDOW_H
