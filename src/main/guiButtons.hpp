#ifndef GUIBUTTONS_HPP
#define GUIBUTTONS_HPP

#include <NanoLiveExceptions.hpp>
#include <SafeQueue.hpp>
#include <SafeMap.hpp>
#include <SafeSet.hpp>
#include <StopClock.hpp>
#include "mainwindow.h"
#include "ibf_mainwindow.h"
#include "classify_mainwindow.h"
#include "connection_test_mainwindow.h"
#include "live_deplete_mainwindow.h"

/**
    pushButton to create an IBF
    buildIBF_qt() in ibfbuild.hpp
*/
/*
void IBF_mainwindow::on_pushButton_3_clicked()
{
    StopClock NanoLiveTime;

    QMessageBox::StandardButton ask;
    ask = QMessageBox::question(this, "Building IBF", "The IBF will be built, are you sure?",
          QMessageBox::Yes|QMessageBox::No);

    if (ask == QMessageBox::Yes) {
     third_party();
     //proBar();
    NanoLiveTime.start();
    buildIBF_qt(k, 0, f, t, ref_file_Name, output_file_Name);

    NanoLiveTime.stop();
    size_t peakSize = getPeakRSS();
    int peakSizeMByte = (int)(peakSize / (1024 * 1024));
    std::cout<<"--------------------------------------------------------------"<<std::endl;
    std::cout << "Real time : " << NanoLiveTime.elapsed() << " sec" << std::endl;
    std::cout << "CPU time  : " << cputime() << " sec" << std::endl;
    std::cout << "Peak RSS  : " << peakSizeMByte << " MByte" << std::endl;
    std::cout<<"--------------------------------------------------------------"<<std::endl;
    QMessageBox::StandardButton ask_1;
    ask_1 = QMessageBox::question(this, "Clear Results", "Do you want to keep the last results in the output window?",
            QMessageBox::Yes|QMessageBox::No);
    if (ask_1 == QMessageBox::Yes) {}

     if (ask_1 == QMessageBox::No) {clearResults();}

    }
  }*/

/**
    pushButton to classify reads
   classify_reads_qt() in classify.hpp
*//*
void Classify_mainwindow::on_pushButton_4_clicked()
{
    StopClock NanoLiveTime;
    if (ibf_deplete_file.length() < 1 && ibf_target_file.length() < 1)
    {
    QMessageBox::information(this, "Warning", "Please provide an IBF file for depletion and/or  an IBF file for targeted sequencing. ");
           return;
    }
           QMessageBox::StandardButton ask;
     ask = QMessageBox::question(this, "Classify", "Nanopore Reads will be classified, are you sure?",
           QMessageBox::Yes|QMessageBox::No);

    if (ask == QMessageBox::Yes) {
    //std::string classified {"class_new.fasta"};
    //std::string unclassified {"unclass_new.fasta"};
    third_party();
     NanoLiveTime.start();
    read_classify_parser_qt classify_parser{ibf_deplete_file_name, ibf_target_file_name, read_file_name, classified,
                       unclassified, false, false, kmer_significance, error_rate, threads,
                       preLen, max_chunks, false};
    classify_reads_qt(classify_parser);
    NanoLiveTime.stop();
    size_t peakSize = getPeakRSS();
    int peakSizeMByte = (int)(peakSize / (1024 * 1024));
    std::cout<<"--------------------------------------------------------------"<<std::endl;
    std::cout << "Real time : " << NanoLiveTime.elapsed() << " sec" << std::endl;
    std::cout << "CPU time  : " << cputime() << " sec" << std::endl;
    std::cout << "Peak RSS  : " << peakSizeMByte << " MByte" << std::endl;
    std::cout<<"--------------------------------------------------------------"<<std::endl;
    QMessageBox::StandardButton ask_1;
    ask_1 = QMessageBox::question(this, "Clear Results", "Do you want to keep the last results in the output window?",
            QMessageBox::Yes|QMessageBox::No);

    if (ask_1 == QMessageBox::Yes) {}

    if (ask_1 == QMessageBox::No) {clearResults();}

             //close();
    } else {}
}*/


/**
    pushButton to test connection

*/
/*
void connection_test_mainwindow::on_pushButton_4_clicked()
    {
        StopClock NanoLiveTime;
        QMessageBox::StandardButton ask;
         ask = QMessageBox::question(this, "Test Connection", "Test connection to a working MinKNOW instance, Are you sure?",
                                       QMessageBox::Yes|QMessageBox::No);

         if (ask == QMessageBox::Yes) {
            third_party();
            NanoLiveTime.start();
            //test_connection_qt(host_name, port, device_name, unblock_all);
            NanoLiveTime.stop();
            size_t peakSize = getPeakRSS();
            int peakSizeMByte = (int)(peakSize / (1024 * 1024));
            std::cout<<"--------------------------------------------------------------"<<std::endl;
            std::cout << "Real time : " << NanoLiveTime.elapsed() << " sec" << std::endl;
            std::cout << "CPU time  : " << cputime() << " sec" << std::endl;
            std::cout << "Peak RSS  : " << peakSizeMByte << " MByte" << std::endl;
            std::cout<<"--------------------------------------------------------------"<<std::endl;
            QMessageBox::StandardButton ask_1;
             ask_1 = QMessageBox::question(this, "Clear Results", "Do you want to keep the last results in the output window?",
                                           QMessageBox::Yes|QMessageBox::No);

             if (ask_1 == QMessageBox::Yes) {}

             if (ask_1 == QMessageBox::No) {clearResults();}
         } else {}
    }
*/
#endif // GUIBUTTONS_HPP
