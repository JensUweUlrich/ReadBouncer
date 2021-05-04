#ifndef IBF_MAINWINDOW_H
#define IBF_MAINWINDOW_H

#include <QDialog>

namespace Ui {
class IBF_mainwindow;
}

class IBF_mainwindow : public QDialog
{
    Q_OBJECT

public:
    explicit IBF_mainwindow(QWidget *parent = nullptr);
    ~IBF_mainwindow();

private slots:
    void on_pushButton_2_clicked();

private:
    Ui::IBF_mainwindow *ui;
};

#endif // IBF_MAINWINDOW_H
