# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 16:33:45 2022

@author: g05296ar
"""

import sys
from PyQt5 import (
    QtWidgets, 
    QtCore, 
    QtGui,
    )
from PyQt5.QtWidgets import (
    QMainWindow, 
    QApplication, 
    QPushButton, 
    QWidget, 
    QAction, 
    QTabWidget,
    QVBoxLayout,
    QFileDialog
    )
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import pyqtSlot
from PyQt5.uic import loadUi

from PyQt5.QtWidgets import QMessageBox
from PyQt5.QtCore import *

from string import Template
import os

from Range_Widget import Range
from Peak_Widget import Peak
from Output_Widget import Output

import cpf
from cpf.settings import settings

# class logger():
#     def __init__():
# ######################################
#         logger = logging.getLogger(__name__)
#         logger.setLevel(logging.INFO)
#         #formatter = logging.Formatter('%(asctime)s | %(levelname)s | %(message)s')
        
#         stdout_handler = logging.StreamHandler(sys.stdout)
#         #stdout_handler.setLevel(logging.DEBUG)
#         #stdout_handler.setFormatter(formatter)
        
#         file_handler = logging.FileHandler('logs.log')
#         #file_handler.setLevel(logging.DEBUG)
#         #file_handler.setFormatter(formatter)
        
        
#         logger.addHandler(file_handler)
#         logger.addHandler(stdout_handler)

#########################################







class App(QMainWindow):

    def __init__(self):
        super(App, self).__init__()
        loadUi("MainWindow.ui", self)
   
        set_cl = cpf.settings.settings()
   # Set Widgets Height 
        self.Main_Tab.setMinimumHeight(40);
        self.Directory.setMinimumHeight(40);
        self.Basename.setMinimumHeight(40);
        self.Extension.setMinimumHeight(40);
        self.Start_Num.setMinimumHeight(40);
        self.End_Num.setMinimumHeight(40);
        self.Num_Digit.setMinimumHeight(40);
        self.Step.setMinimumHeight(40);
        
        self.Calib_Type.setMinimumHeight(40);
        self.Calib_Detect.setMinimumHeight(40);
        self.Calib_Param.setMinimumHeight(40);
        self.Calib_Mask.setMinimumHeight(40);
        self.Calib_Pixels.setMinimumHeight(40);
        self.Calib_Data.setMinimumHeight(40);
        
        self.fit_bin_type.setMinimumHeight(40);
        self.fit_per_bin.setMinimumHeight(40);
        self.fit_number_bins.setMinimumHeight(40);
        
        self.bg_1.setMinimumHeight(40);
        self.bg_2.setMinimumHeight(40);
        self.ds_1.setMinimumHeight(40);
        self.ds_2.setMinimumHeight(40);
        self.h_1.setMinimumHeight(40);
        self.h_2.setMinimumHeight(40);
        self.pro_1.setMinimumHeight(40);
        self.pro_2.setMinimumHeight(40);
        self.wdt_1.setMinimumHeight(40);
        self.wdt_2.setMinimumHeight(40);
        
        self.AziBins.setMinimumHeight(40);
        self.Output_Dir_1.setMinimumHeight(40);
        self.Output_Dir_2.setMinimumHeight(40);
        
        self.Console_output.setReadOnly(True)
        self.Console_output.textCursor
        
        
        # Function Call aginast button click event   
        self.AddRange_Btn.clicked.connect(self.InsertRange)
        self.RemoveRange_Btn.clicked.connect(self.RemoveRange)
        
        self.Add_Output_Btn.clicked.connect(self.InsertOut)
        self.Remove_Output_Btn.clicked.connect(self.RemoveOut)
        self.Select_Directory_1.clicked.connect(lambda: select_Data_Dir())
        self.Select_Directory_2.clicked.connect(lambda: selectTarget())
        self.Save_input_Btn.clicked.connect(lambda: Insert_Button())
        
        self.input_file_path = None   # global variable
        self.Load_input_Btn.clicked.connect(lambda: Load_Inputs())
        
        self.validate_Btn.clicked.connect(lambda: Validate_Inputs())
        
        self.Run_Range_Btn.clicked.connect(lambda: Run_Range())
        
        self.Run_initial_Btn.clicked.connect(lambda: Run_Initial_Position())
        
        self.Execute_Btn.clicked.connect(lambda: Execute_Fits())
        
        self.Make_Output_Btn.clicked.connect(lambda: Make_Outputs())
        
        
        def Make_Outputs():
            if self.input_file_path == None:
                mess = QMessageBox()
                mess.setIcon(QMessageBox.Warning)
                mess.setText("Please Load input file")
                mess.setStandardButtons(QMessageBox.Ok)
                mess.setWindowTitle("MessageBox")
                returnValue = mess.exec_()
            else:
                cpf.XRD_FitPattern.initiate('./Example1-Fe/BCC1_input_Dioptas')
                text=open('employee.log').read()
                self.DisplayOutput.setText(text)
        
        
        ################################### Run Execute_Fits Button ######################
        
        def Execute_Fits():
            if self.input_file_path == None:
                mess = QMessageBox()
                mess.setIcon(QMessageBox.Warning)
                mess.setText("Please Load input file")
                mess.setStandardButtons(QMessageBox.Ok)
                mess.setWindowTitle("MessageBox")
                returnValue = mess.exec_()
            else:
                cpf.XRD_FitPattern.initiate('./Example1-Fe/BCC1_input_Dioptas')
                text=open('employee.log').read()
                self.DisplayOutput.setText(text)
                
                
        
        ################################### Run Range Button ######################
        def Run_Initial_Position():
            if self.input_file_path == None:
                mess = QMessageBox()
                mess.setIcon(QMessageBox.Warning)
                mess.setText("Please Load input file")
                mess.setStandardButtons(QMessageBox.Ok)
                mess.setWindowTitle("MessageBox")
                returnValue = mess.exec_()
            else:
                cpf.XRD_FitPattern.initiate('./Example1-Fe/BCC1_input_Dioptas')
                text=open('employee.log').read()
                self.DisplayOutput.setText(text)
        
        
        
        
        
        
        
        ################################### Run Range Button ######################
        
        def Run_Range():
            if self.input_file_path == None:
                mess = QMessageBox()
                mess.setIcon(QMessageBox.Warning)
                mess.setText("Please Load input file")
                mess.setStandardButtons(QMessageBox.Ok)
                mess.setWindowTitle("MessageBox")
                returnValue = mess.exec_()
            else:
                cpf.XRD_FitPattern.initiate('./Example1-Fe/BCC1_input_Dioptas')
                text=open('employee.log').read()
                self.DisplayOutput.setText(text)
        
        


        
        ################################### Validate Input Button ######################
        def Validate_Inputs():
            # if self.input_file_path == None:
            #     mess = QMessageBox()
            #     mess.setIcon(QMessageBox.Warning)
            #     mess.setText("Please Load input file")
            #     mess.setStandardButtons(QMessageBox.Ok)
            #     mess.setWindowTitle("MessageBox")
            #     returnValue = mess.exec_()
            # else:
                with open("logs.log",'a') as file:
             # read the lines
                   file.close()
                cpf.XRD_FitPattern.initiate('./Example1-Fe/BCC1_input_Dioptas')
                text=open('logs.log').read()
                self.Console_output.setText(text)
        
        
        
        
        ################################### Load Input Button ######################
        def Load_Inputs():
            
            # When button is clicked clear all data in /result directory
            mydir = "./results"
            filelist = [f for f in os.listdir(mydir) if f.endswith(".png") ]
            for f in filelist:
                os.remove(os.path.join(mydir, f))
             
                
            # Open file dialog to choose input file
            fname= QFileDialog.getOpenFileName(self, "Load Input File", "", "Python Files (*.py)")
            
            if fname:
                self.input_file_path = f"{fname[0]}"
                print(self.input_file_path)
            
            
            set_cl.populate(settings_file=(f"{self.input_file_path}"))
            self.Directory.setText(set_cl.datafile_directory)
            self.Basename.setText(set_cl.datafile_basename)
            # self.Extension.setText(40);
            # self.Start_Num.setText(40);
            # self.End_Num.setText(40);
            # self.Num_Digit.setText(40);
            # self.Step.setText(40);
            
            self.Calib_Type.setWhatsThis(set_cl.calibration_type);
            self.Calib_Detect.setText(set_cl.calibration_detector);
            self.Calib_Param.setText(set_cl.calibration_parameters); 
            self.Calib_Mask.setText(set_cl.calibration_mask);
            self.Calib_Pixels.setText(str(set_cl.calibration_pixel_size));
            self.Calib_Data.setText(set_cl.calibration_data);
            
            self.fit_bin_type.setText(str(set_cl.fit_bin_type))
            self.fit_per_bin.setText(str(set_cl.fit_per_bin))
            self.fit_number_bins.setText(str(set_cl.fit_number_bins))
            
               
            with open("logs.log",'w') as file:
         # read the lines
               file.close()
            cpf.XRD_FitPattern.initiate('./Example1-Fe/BCC1_input_Dioptas')
            text=open('logs.log').read()
            self.Console_output.setText(text)
                
                
                
     ################################### Insert Input Button ######################
        
        def Insert_Button():  
        # Get text from textboxes
            directory = self.Directory.text()
            basename = self.Basename.text()
            extension = self.Extension.text()
            start_num = self.Start_Num.text()
            end_num = self.End_Num.text()
            num_digit = self.Num_Digit.text()
            step = self.Step.text()
            
            calib_type = self.Calib_Type.currentText()
            calib_detect = self.Calib_Detect.text()
            calib_data = self.Calib_Data.text()
            calib_param = self.Calib_Param.text()
            calib_mask = self.Calib_Mask.text()
            calib_pixel = self.Calib_Pixels.text()
            
            Output_Dir_1 = self.Output_Dir_1.text()
            
            
            bg_1 = self.bg_1.text()
            bg_2 = self.bg_2.text()
            ds_1 = self.ds_1.text()
            ds_2 = self.ds_2.text()
            h_1 = self.h_1.text()
            h_2 = self.h_2.text()
            pro_1 = self.pro_1.text()
            pro_2 = self.pro_2.text()
            wdt_1 =self.wdt_1.text()
            wdt_2 =self.wdt_2.text()

            AziBins = self.AziBins.text()
            
            
            ### create file           
                
            data = "\
 datafile_directory = '$dr/'         \n\
 datafile_Basename  = '$bn'     \n\
 datafile_Ending    = '$ext'                   \n\
 datafile_StartNum  = $sn                          \n\
 datafile_EndNum    = $en                           \n\
 datafile_NumDigit  = $nd                       \n\
 datafile_Step      = $increment                          \n\
# Calibration and masking.                                 \n\
 Calib_type   = \"$ct\"                      \n\
 Calib_detector = '$cdetect'                  \n\
 Calib_data     = datafile_directory + '$cadat'       \n\
 Calib_param    = datafile_directory + '$cparam'                \n\
 Calib_mask     = datafile_directory + '$cmask'          \n\
 Calib_pixels = $cpixel                            \n\
# Number of bins for initial fitting.                  \n\
 AziBins = $AziBins                             \n\
# Output settings                                 \n\
 Output_directory   = '$Output_Dir/'                    \n\
 fit_bounds = {                                \n\
   \"background\": ['$bg_1', '$bg_2'],                  \n\
   \"d-space\":    ['$ds_1', '$ds_2'],                     \n\
   \"height\":     [ $h_1,   '$h_2'],                     \n\
   \"profile\":    [$pro_1, $pro_2 ],                           \n\
   \"width\":      ['$wdt_1',  '$wdt_2'],        \n\
   }  \n\
 fit_orders = ["
            
            temp_obj = Template(data)
            after_replacing = temp_obj.substitute(dr=directory, bn=basename,ext=extension, 
                                                  sn=start_num, en=end_num, nd=num_digit, 
                                                  increment=step, ct=calib_type,cdetect=calib_detect,
                                                  cadat=calib_data,cparam=calib_param,
                                                  cmask=calib_mask, cpixel=calib_pixel,AziBins = AziBins,
                                                  bg_1= bg_1, bg_2 = bg_2, ds_1= ds_1,
                                                  ds_2= ds_2, h_1=h_1, h_2= h_2, pro_1= pro_1,
                                                  pro_2= pro_2, wdt_1= wdt_1, wdt_2= wdt_2,
                                                  Output_Dir= Output_Dir_1)
            
        # Show Message Box 
            if directory =='' or basename =='' or extension =='' or start_num =='' or end_num =='' or  calib_type=='' or calib_param=='':
                      mess = QMessageBox()
                      mess.setIcon(QMessageBox.Warning)
                      mess.setText("Insert All Mandaory Details")
                      mess.setStandardButtons(QMessageBox.Ok)
                      mess.setWindowTitle("MessageBox")
                      returnValue = mess.exec_()
                      if os.path.exists("template_3.py"):
                         os.remove("template_3.py")
                      
            else:
                      file1 = open('template_3.py', 'w')
                      file1.writelines((after_replacing))
                      file1.close()
                      #Reset_Button()                                  # After inserting details clear text boxes
                      mess = QMessageBox()
                      mess.setText("Success")
                      mess.setIcon(QMessageBox.Information)
                      mess.setStandardButtons(QMessageBox.Ok)
                      mess.setWindowTitle("MessageBox")
                      returnValue = mess.exec_()
            
        self.clickcount = 0 
        self.clickcountout = 0 
        
    ############## Select Location 1 ##################
        def select_Data_Dir():
                dialog = QtWidgets.QFileDialog()
                if self.Directory.text():
                    dialog.setDirectory(self.Directory.text())

                dialog.setFileMode(dialog.Directory)

                # we cannot use the native dialog, because we need control over the UI
                options = dialog.Options(dialog.DontUseNativeDialog | dialog.ShowDirsOnly)
                dialog.setOptions(options)

                def checkLineEdit(path):
                    if not path:
                        return
                    if path.endswith(QtCore.QDir.separator()):
                        return checkLineEdit(path.rstrip(QtCore.QDir.separator()))
                    path = QtCore.QFileInfo(path)
                    if path.exists() or QtCore.QFileInfo(path.absolutePath()).exists():
                        button.setEnabled(True)
                        return True

                # get the "Open" button in the dialog
                button = dialog.findChild(QtWidgets.QDialogButtonBox).button(
                    QtWidgets.QDialogButtonBox.Open)

                # get the line edit used for the path
                lineEdit = dialog.findChild(QtWidgets.QLineEdit)
                lineEdit.textChanged.connect(checkLineEdit)

                # override the existing accept() method, otherwise selectedFiles() will 
                # complain about selecting a non existing path
                def accept():
                    if checkLineEdit(lineEdit.text()):
                        # if the path is acceptable, call the base accept() implementation
                        QtWidgets.QDialog.accept(dialog)
                dialog.accept = accept

                if dialog.exec_() and dialog.selectedFiles():
                    path = QtCore.QFileInfo(dialog.selectedFiles()[0]).absoluteFilePath()
                    self.Directory.setText(path)
                    
    ############# Select Location2 #################       
        def selectTarget():
                dialog = QtWidgets.QFileDialog()

                if self.Output_Dir_1.text():
                    dialog.setDirectory(self.Output_Dir_1.text())

                dialog.setFileMode(dialog.Directory)

                # we cannot use the native dialog, because we need control over the UI
                options = dialog.Options(dialog.DontUseNativeDialog | dialog.ShowDirsOnly)
                dialog.setOptions(options)
               
                def checkLineEdit(path):
                    if not path:
                        return
                    if path.endswith(QtCore.QDir.separator()):
                        return checkLineEdit(path.rstrip(QtCore.QDir.separator()))
                    path = QtCore.QFileInfo(path)
                    if path.exists() or QtCore.QFileInfo(path.absolutePath()).exists():
                        button.setEnabled(True)
                        return True

                # get the "Open" button in the dialog
                button = dialog.findChild(QtWidgets.QDialogButtonBox).button(
                    QtWidgets.QDialogButtonBox.Open)

                # get the line edit used for the path
                lineEdit = dialog.findChild(QtWidgets.QLineEdit)
                lineEdit.textChanged.connect(checkLineEdit)
                

                # override the existing accept() method, otherwise selectedFiles() will 
                # complain about selecting a non existing path
                def accept():
                    if checkLineEdit(lineEdit.text()):
                        # if the path is acceptable, call the base accept() implementation
                        QtWidgets.QDialog.accept(dialog)
                dialog.accept = accept

                if dialog.exec_() and dialog.selectedFiles():
                    path = QtCore.QFileInfo(dialog.selectedFiles()[0]).absoluteFilePath()
                    self.Output_Dir_1.setText(path)
                    self.Output_Dir_2.setText(path)
                    
      
                
      
      ################################### Add Range Button ######################   
        
    def InsertRange(self):
        self.clickcount += 1 
        self.Range_Tab.addTab(Range() , QIcon("Location of the icon"),"Range "+str(self.clickcount))
        
        
        
      
     ################################### Remove Range Button ###################### 
        
    def RemoveRange(self):
        self.Range_Tab.removeTab(self.Range_Tab.currentIndex())
        
        
        
     ################################### Add Output Button ######################        
    def InsertOut(self):
        self.clickcountout += 1 
        self.Output_Tab.addTab(Output() , QIcon("Location of the icon"),"Output "+str(self.clickcountout))
    
     
    
    ################################### Remove Output Button ######################
    def RemoveOut(self):
        self.Output_Tab.removeTab(self.Output_Tab.currentIndex())
    

# Main
app = QApplication(sys.argv)
mainwindow=App()
widget= QtWidgets.QStackedWidget()
widget.addWidget(mainwindow)
widget.show()

try: 
    sys.exit(app.exec_())
except:
    print("Exiting")