#ifndef CONNECTION_TEST_MAINWINDOW_H
#define CONNECTION_TEST_MAINWINDOW_H

#include <QDialog>

namespace Ui {
class connection_test_mainwindow;
}

class connection_test_mainwindow : public QDialog
{
    Q_OBJECT

public:
    explicit connection_test_mainwindow(QWidget *parent = nullptr);
    ~connection_test_mainwindow();

private:
    Ui::connection_test_mainwindow *ui;
};

#endif // CONNECTION_TEST_MAINWINDOW_H
