GUI with Qt6.0.3 for ReadBouncer How to use:

Deadline for Readme: 29.06.2021

-------------------------------------------------------------------------------------------------------------------------
**Requirements for Qt (All this kits are automatically installed with Qt)**: 

* CDB Debugger 
* MSVC 2019  x86_64
* QMake 
* CMake tool for QT 6.0.1
* UI Designer 

-------------------------------------------------------------------------------------------------------------------------
**Installing Qt (Windows + Linux)**: 
* Download Qt installer from (https://www.qt.io/download), Downloads for open source users (free version).
* Create a free account for Qt users. 
* Install Qt 6.0.3 or 6.1.2 with all components such as MSVC 2019 64-bit. 
* Install from the same window the CDB Debugger, to debug the application from Qt. 
-------------------------------------------------------------------------------------------------------------------------
**Run ReadBouncer with Qt**:
* After installing Qt you are now able to open Qt-Creator. 
* Open Qt Creator and go to Tools -> Options.
* From Kits choose the kit with MSVC2019. 
* In Debugger set in CDB `Break on: C++ exception`.
* From Qt Creator choose File -> Open File or Project -> select the CMakeLists.txt from `ReadBouncer/src/CMakeLists.txt`.
* While Qt start building ReadBouncer, stop the process and go to `Projects` under Debug in Qt-Creator main window. 
* Check that Qt-Creator is building the ReadBouncer with the MSVC2019 kit. 
* Under CMake in `Projects` change the build directory to the `ReadBouncer/build`.
* Make sure that the initial CMake parameters are:  
   
  `-GVisual Studio 16 2019`  
  `-DCMAKE_BUILD_TYPE:String=RelWithDebInfo`  
  `-DQT_QMAKE_EXECUTABLE:STRING=%{Qt:qmakeExecutable}`   
  `-DCMAKE_PREFIX_PATH:STRING=%{Qt:QT_INSTALL_PREFIX}`  
  `-DCMAKE_C_COMPILER:STRING=%{Compiler:Executable:C}`  
  `-DCMAKE_CXX_COMPILER:STRING=%{Compiler:Executable:Cxx}`  
   
* Under Build Steps: `cmake.exe --build . --target ALL_BUILD`  **NOT**  `cmake.exe --build . --target all`.
* Under Run change the executable to: `C:\NanoLive\build\main\Debug\NanoLive.exe` and the working directory should be: `C:\NanoLive\build\main\Debug\`.
* Just save the changes (**Ctrl+S**) and Qt should rebuild ReadBouncer again automatically.
* After finishing the build step (need about 40-60 minutes), just click on run (**Ctrl+R**) and the GUI will open.  
-------------------------------------------------------------------------------------------------------------------------
**Changing the GUI**:
* For changing you can open the `.ui` files such as `MainWindow.ui`.
* The libraries are designed to generate the files in the right order, sothat you don't need to open any other files.  
-------------------------------------------------------------------------------------------------------------------------
**Deploy Executable**:
* To run the ReadBouncer.exe from `C:\NanoLive\build\main\Debug\NanoLive.exe` you have to add some Qt files, which allow to deploy the executable from there, this libraries should do the Job:   
  - Copy the directory `C:\QT\QT5.12.10\6.0.3\msvc2019_64\plugins` to `C:\NanoLive\build\main\Debug\`      
  - Copy the files: `Qt5Core.dll` `Qt5Cored.dll` `Qt5Guid.dll` `Qt5Widgetsd.dll`   
  - Copy all weights files from `ReadBouncer\weights` to  `C:\NanoLive\build\main\Debug\` .   
  - Copy the `ReadBouncer/src/rpc-certs/ca.crt` to `C:\NanoLive\build\main\Debug\` .    
-------------------------------------------------------------------------------------------------------------------------
**Note**:
* While using c++ exception, **avoid using** `throw` because of CDB-debugger, using `throw` as exception will stop the software, without using `throw` Qt and CDB debugger will do the job also. 
