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

from matplotlib_qt import matplotlib_qt
from matplotlib_auto import matplotlib_inline


class MainWindow(QMainWindow):
    
    def __init__(self):
        super(MainWindow, self).__init__()
        loadUi("MainWindow.ui", self)
        self.gui_layout()
        self.set_cl = cpf.settings.settings()
        self.input_file_path = None   
        self.clickcount = 0 
        self.clickcountout = 0 
        
    def gui_layout(self):
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
        self.cascade_bin_type.setMinimumHeight(40);
        self.cascade_per_bin.setMinimumHeight(40);
        self.cascade_number_bins.setMinimumHeight(40);
        self.cascade_track.setMinimumHeight(40);
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
        self.AddRange_Btn.clicked.connect(self.Insert_Range)
        self.RemoveRange_Btn.clicked.connect(self.Remove_Range)
        self.Add_Output_Btn.clicked.connect(self.Insert_Output)
        self.Remove_Output_Btn.clicked.connect(self.Remove_Output)
        self.Select_Directory_1.clicked.connect(self.select_Data_Dir)
        self.Select_Directory_2.clicked.connect(self.selectTarget)
        self.Save_input_Btn.clicked.connect(self.Insert_Button)
        self.Load_input_Btn.clicked.connect(self.Load_Inputs)
        self.validate_Btn.clicked.connect(self.Validate_Inputs)
        self.Run_Range_Btn.clicked.connect(self.Run_Range)
        self.Run_initial_Btn.clicked.connect(self.Run_Initial_Position)
        self.Execute_Btn.clicked.connect(self.Execute_Fits)
        self.Make_Output_Btn.clicked.connect(self.Make_Outputs)
    
    def Load_Inputs(self):
        fname= QFileDialog.getOpenFileName(self, "Load Input File", "..\\", "Python Files (*.py)")
        if fname:
            self.input_file_path = f"{fname[0]}"
        self.set_cl.populate(settings_file=(f"{self.input_file_path}"))
        self.Directory.setText(self.set_cl.datafile_directory)
        self.Basename.setText(self.set_cl.datafile_basename)
        # self.Extension.setText(40);
        # self.Start_Num.setText(40);
        # self.End_Num.setText(40);
        # self.Num_Digit.setText(40);
        # self.Step.setText(40);
        self.Calib_Type.setCurrentText(f"{self.set_cl.calibration_type}");
        self.Calib_Detect.setText(self.set_cl.calibration_detector);
        self.Calib_Param.setText(self.set_cl.calibration_parameters); 
        self.Calib_Mask.setText(self.set_cl.calibration_mask);
        self.Calib_Pixels.setText(str(self.set_cl.calibration_pixel_size));
        self.Calib_Data.setText(self.set_cl.calibration_data);
        self.cascade_bin_type.setText(str(self.set_cl.cascade_bin_type))
        self.cascade_per_bin.setText(str(self.set_cl.cascade_per_bin))
        self.cascade_number_bins.setText(str(self.set_cl.cascade_number_bins))
        self.cascade_track.setText(str(self.set_cl.cascade_track))
        self.bg_1.setText(self.set_cl.fit_bounds.get("background")[0])
        self.bg_2.setText(self.set_cl.fit_bounds.get("background")[1])
        self.ds_1.setText(self.set_cl.fit_bounds.get("d-space")[0])
        self.ds_2.setText(self.set_cl.fit_bounds.get("d-space")[1])
        self.h_1.setText(str(self.set_cl.fit_bounds.get("height")[0]))
        self.h_2.setText(str(self.set_cl.fit_bounds.get("height")[1]))
        self.pro_1.setText(str(self.set_cl.fit_bounds.get("profile")[0]))
        self.pro_2.setText(str(self.set_cl.fit_bounds.get("profile")[1]))
        self.wdt_1.setText(self.set_cl.fit_bounds.get("width")[0])
        self.wdt_2.setText(self.set_cl.fit_bounds.get("width")[1])
        self.Output_Dir_1.setText(self.set_cl.output_directory)
        
        self.range_length = len(self.set_cl.fit_orders)
        
        for remove_range in range(0, self.Range_Tab.count()):
            self.Remove_Range()
        
        for range_tab in range(0, self.range_length):
            self.range_object = Range()
            self.Range_Tab.addTab(self.range_object , QIcon(""),"Range")
            self.range_object.Range_min.setText(str(self.set_cl.fit_orders[range_tab].get("range")[0]))
            self.range_object.Range_max.setText(str(self.set_cl.fit_orders[range_tab].get("range")[1]))
            self.range_object.Range_Background_Val.setText(str(self.set_cl.fit_orders[range_tab].get("background")))
            self.range_object.Intensity_max.setText(str(self.set_cl.fit_orders[range_tab].get("Imax")))   
            self.range_object.Intensity_min.setText(str(self.set_cl.fit_orders[range_tab].get("Imin")))   
            self.range_object.Peak_Pos_Selection.setPlainText(str(self.set_cl.fit_orders[range_tab].get("PeakPositionSelection")))
            self.peak_length = len(self.set_cl.fit_orders[range_tab]["peak"])
            
            for peak_tab in range(0, self.peak_length):
                self.peak_object = Peak()
                self.range_object.Peak_Tab.addTab(self.peak_object , QIcon(""),"Peak ")#+str(self.clickcount))
                self.peak_object.phase_peak.setText(str(self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("phase")))
                self.peak_object.hkl.setText(str(self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("hkl")))
                self.peak_object.d_space_peak.setText(str(self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("d-space")))
                #self.peak_object.d_space_type.setText()   
                #self.peak_object.dspace_fixed.setText(str(self.set_cl.fit_orders[range_tab].get("d-space")))
                self.peak_object.height_peak.setText(str(self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("height")))
                #self.peak_object.height_peak_type.setText()
                #self.peak_object.height_fixed.setText(str(self.set_cl.fit_orders[range_tab].get("height")))
                self.peak_object.profile_peak.setText(str(self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("profile")))
                #self.peak_object.profile_peak_type.setText()
                self.peak_object.profile_fixed.setText(str(self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("profile_fixed")))
                self.peak_object.width_peak.setText(str(self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("width")))
                #self.peak_object.width_peak_type.setText()
                #self.peak_object.width_fixed.setText(str(self.set_cl.fit_orders[range_tab].get("height")))
                self.peak_object.symmetry_peak.setText(str(self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("symmetry")))
                
        with open("../logs/logs.log",'w') as file:
           file.close()
        cpf.XRD_FitPattern.initiate(f"{self.input_file_path}")
        text=open('../logs/logs.log').read()
        self.Console_output.setText(text)
        
    def Validate_Inputs(self):
            with open("../logs/logs.log",'a') as file:
               file.close()
            cpf.XRD_FitPattern.initiate(f"{self.input_file_path}")
            text=open('../logs/logs.log').read()
            self.Console_output.setText(text)
                   
    def Run_Range(self):
        self.matplotlib_inline = matplotlib_inline()
        if self.input_file_path == None:
            mess = QMessageBox()
            mess.setIcon(QMessageBox.Warning)
            mess.setText("Please Load input file")
            mess.setStandardButtons(QMessageBox.Ok)
            mess.setWindowTitle("MessageBox")
            returnValue = mess.exec_()
        else:
            cpf.XRD_FitPattern.set_range(f"{self.input_file_path}", save_all=True)
            text=open('../logs/logs.log').read()
            self.Console_output.setText(text)
   
    def Run_Initial_Position(self):
        self.matplotlib_qt = matplotlib_qt()
        if self.input_file_path == None:
            mess = QMessageBox()
            mess.setIcon(QMessageBox.Warning)
            mess.setText("Please Load input file")
            mess.setStandardButtons(QMessageBox.Ok)
            mess.setWindowTitle("MessageBox")
            returnValue = mess.exec_()
        else:
            cpf.XRD_FitPattern.initial_peak_position(f"{self.input_file_path}", save_all=True)
            text=open('../logs/logs.log').read()
            self.Console_output.setText(text)
   
    def Execute_Fits(self):
        if self.input_file_path == None:
            mess = QMessageBox()
            mess.setIcon(QMessageBox.Warning)
            mess.setText("Please Load input file")
            mess.setStandardButtons(QMessageBox.Ok)
            mess.setWindowTitle("MessageBox")
            returnValue = mess.exec_()
        else:
            cpf.XRD_FitPattern.execute(f"{self.input_file_path}", save_all=True)
            text=open('../logs/logs.log').read()
            self.Console_output.setText(text)
 
    def Make_Outputs(self):
        if self.input_file_path == None:
            mess = QMessageBox()
            mess.setIcon(QMessageBox.Warning)
            mess.setText("Please Load input file")
            mess.setStandardButtons(QMessageBox.Ok)
            mess.setWindowTitle("MessageBox")
            returnValue = mess.exec_()
        else:
            cpf.XRD_FitPattern.write_output(f"{self.input_file_path}", save_all=True)
            text=open('../logs/logs.log').read()
            self.Console_output.setText(text)

    def Insert_Button(self):  
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
                      mess = QMessageBox()
                      mess.setText("Success")
                      mess.setIcon(QMessageBox.Information)
                      mess.setStandardButtons(QMessageBox.Ok)
                      mess.setWindowTitle("MessageBox")
                      returnValue = mess.exec_()

    def select_Data_Dir(self):
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

    def selectTarget(self):
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
  
            def accept():
                    if checkLineEdit(lineEdit.text()):
                        # if the path is acceptable, call the base accept() implementation
                        QtWidgets.QDialog.accept(dialog)
            dialog.accept = accept

            if dialog.exec_() and dialog.selectedFiles():
                    path = QtCore.QFileInfo(dialog.selectedFiles()[0]).absoluteFilePath()
                    self.Output_Dir_1.setText(path)
                    self.Output_Dir_2.setText(path)

    def Insert_Range(self):
        self.clickcount += 1 
        self.Range_Tab.addTab(Range() , QIcon("Location of the icon"),"Range "+str(self.clickcount))

    def Remove_Range(self):
        self.Range_Tab.removeTab(self.Range_Tab.currentIndex())
    
    def Insert_Output(self):
        self.clickcountout += 1 
        self.Output_Tab.addTab(Output() , QIcon("Location of the icon"),"Output "+str(self.clickcountout))

    def Remove_Output(self):
        self.Output_Tab.removeTab(self.Output_Tab.currentIndex())

if __name__=='__main__':
    app = QApplication(sys.argv)
    mainwindow = MainWindow()
    widget = QtWidgets.QStackedWidget()
    widget.addWidget(mainwindow)
    widget.show() 
    try: 
        sys.exit(app.exec_())
    except:
        os._exit(00)
