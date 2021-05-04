#ifndef LIVE_DEPLETE_MAINWINDOW_H
#define LIVE_DEPLETE_MAINWINDOW_H

#include <QDialog>

namespace Ui {
class live_deplete_mainwindow;
}

class live_deplete_mainwindow : public QDialog
{
    Q_OBJECT

public:
    explicit live_deplete_mainwindow(QWidget *parent = nullptr);
    ~live_deplete_mainwindow();

private slots:
    void on_pushButton_clicked();

private:
    Ui::live_deplete_mainwindow *ui;
};

#endif // LIVE_DEPLETE_MAINWINDOW_H
