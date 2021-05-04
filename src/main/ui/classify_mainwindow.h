#ifndef CLASSIFY_MAINWINDOW_H
#define CLASSIFY_MAINWINDOW_H

#include <QDialog>

namespace Ui {
class Classify_mainwindow;
}

class Classify_mainwindow : public QDialog
{
    Q_OBJECT

public:
    explicit Classify_mainwindow(QWidget *parent = nullptr);
    ~Classify_mainwindow();

private slots:
    void on_pushButton_clicked();

private:
    Ui::Classify_mainwindow *ui;
};

#endif // CLASSIFY_MAINWINDOW_H
