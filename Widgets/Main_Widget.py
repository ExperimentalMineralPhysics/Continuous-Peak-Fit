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
        loadUi("Main_Widget.ui", self)
        self.setWindowIcon(QtGui.QIcon('logo.png'))
        self.gui_layout()
        self.set_cl = cpf.settings.settings()
        self.input_file_path = None  
    
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
        
    @pyqtSlot()
    def Load_Inputs(self):
        fname= QFileDialog.getOpenFileName(self, "Load Input File", "..\\", "Python Files (*.py)")
        if fname:
            self.input_file_path = f"{fname[0]}"
        self.set_cl.populate(settings_file=(f"{self.input_file_path}"))
        self.Directory.setText(self.set_cl.datafile_directory)
        self.Basename.setText(self.set_cl.datafile_basename)
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
        self.Output_Dir_2.setText(self.Output_Dir_1.text())
        #self.AziBins.setText(self.set_cl.AziBins)  AttributeError: 'settings' object has no attribute 'AziBins'
        # Signals #
        self.Directory.editingFinished.connect(self.Dir_Pressed)
        self.Basename.editingFinished.connect(self.Basename_Pressed)
        self.Calib_Detect.editingFinished.connect(self.Calib_Detect_Pressed)
        self.Calib_Param.editingFinished.connect(self.Calib_Param_Pressed)
        self.Calib_Mask.editingFinished.connect(self.Calib_Mask_Pressed)
        self.Calib_Pixels.editingFinished.connect(self.Calib_Pixels_Pressed)
        self.Calib_Data.editingFinished.connect(self.Calib_Data_Pressed)
        self.cascade_bin_type.editingFinished.connect(self.Cascade_bin_type_Pressed)
        self.cascade_per_bin.editingFinished.connect(self.Cascade_per_bin_Pressed)
        self.cascade_number_bins.editingFinished.connect(self.Cascade_number_bins_Pressed)
        self.cascade_track.editingFinished.connect(self.Cascade_track_Pressed)
        self.bg_1.editingFinished.connect(self.Bg_1_Pressed)
        self.bg_2.editingFinished.connect(self.Bg_2_Pressed)
        self.ds_1.editingFinished.connect(self.Ds_1_Pressed)
        self.ds_2.editingFinished.connect(self.Ds_2_Pressed)
        self.h_1.editingFinished.connect(self.H_1_Pressed)
        self.h_2.editingFinished.connect(self.H_2_Pressed)
        self.pro_1.editingFinished.connect(self.Pro_1_Pressed)
        self.pro_2.editingFinished.connect(self.Pro_2_Pressed)
        self.wdt_1.editingFinished.connect(self.Wdt_1_Pressed)
        self.wdt_2.editingFinished.connect(self.Wdt_2_Pressed)
        self.Output_Dir_1.editingFinished.connect(self.Output_Dir_1_Pressed)
        self.range_length = len(self.set_cl.fit_orders)
        # Remove existing range_tabs whenever new file is loaded
        for remove_range in range(0, self.Range_Tab.count()):
            self.Remove_Range()
        # Manage range_tab data population
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
            
            if self.set_cl.fit_orders[range_tab].get("background-type")== None:
                self.range_object.bg_fixed_checkbox.setChecked(False)
            else:
                self.range_object.bg_fixed_checkbox.setChecked(True)
                self.b_type = str(self.set_cl.fit_orders[range_tab].get("background-type"))
                self.range_object.Background_Type.setCurrentText(f"{self.b_type}")
            # Manage peak_tab data population
            for peak_tab in range(0, self.peak_length):
                self.peak_object = Peak()
                self.range_object.Peak_Tab.addTab(self.peak_object , QIcon(""),"Peak")
                self.peak_object.phase_peak.setText(str(self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("phase")))
                self.peak_object.hkl.setText(str(self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("hkl")))
                self.peak_object.d_space_peak.setText(str(self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("d-space")))
                self.ds_type = str(self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("d-space-type"))
                self.peak_object.d_space_type.setCurrentText(f"{self.ds_type}")
                self.h_type = str(self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("height-type"))
                self.peak_object.height_peak_type.setCurrentText(f"{self.h_type}")
                self.p_type = str(self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("profile-type"))
                self.peak_object.profile_peak_type.setCurrentText(f"{self.p_type}")
                self.w_type = str(self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("width-type"))
                self.peak_object.width_peak_type.setCurrentText(f"{self.w_type}")
                self.peak_object.height_peak.setText(str(self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("height")))
                self.peak_object.profile_peak.setText(str(self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("profile")))
                self.peak_object.profile_fixed.setText(str(self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("profile_fixed")))
                self.peak_object.width_peak.setText(str(self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("width")))
                self.peak_object.symmetry_peak.setText(str(self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("symmetry")))

                if self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("profile_fixed")== None:
                    self.peak_object.profile_checkBox.setChecked(False)
                else:
                    self.peak_object.dspace_fixed.setText(str(self.set_cl.fit_orders[range_tab].get("profile_fixed")))
                    self.peak_object.profile_checkBox.setChecked(True)
                    
                if self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("width_fixed")== None:
                    self.peak_object.width_checkBox.setChecked(False)
                else:
                    self.peak_object.width_fixed.setText(str(self.set_cl.fit_orders[range_tab].get("width_fixed")))
                    self.peak_object.width_checkBox.setChecked(True)
                    
                if self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("height_fixed")== None:
                    self.peak_object.height_checkBox.setChecked(False)
                else:
                    self.peak_object.width_fixed.setText(str(self.set_cl.fit_orders[range_tab].get("height_fixed")))
                    self.peak_object.height_checkBox.setChecked(True)
                    
                if self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("d-space_fixed")== None:
                    self.peak_object.dspace_checkBox.setChecked(False)
                else:
                    self.peak_object.width_fixed.setText(str(self.set_cl.fit_orders[range_tab].get("d-space_fixed")))
                    self.peak_object.dspace_checkBox.setChecked(True)
        
        # Remove existing output_tabs whenever new file is loaded
        for remove_range in range(0, self.Output_Tab.count()):
            self.Remove_Output()
        # Manage output_tab data population    
        self.output_type =  len(self.set_cl.output_types)
        for output in range (0, self.output_type):
            self.output_object = Output()
            self.Output_Tab.addTab(self.output_object , QIcon(""),"Output")

            if 'Polydefix' == self.set_cl.output_types[output]:
                    self.output_object.Output_Type_comboBox.setCurrentText('WritePolydefix')
                # Optional Params
                    
                    for wid in range (len(cpf.output_formatters.WritePolydefix.Requirements()[1])):
                        lineEdit = cpf.output_formatters.WritePolydefix.Requirements()[1][wid]
                        if lineEdit == 'Output_ElasticProperties':
                            self.output_object.lineEdit.setText(str(self.set_cl.output_settings["Output_ElasticProperties"]))
                        if lineEdit == 'phase':
                            self.output_object.lineEdit.setText(str(self.set_cl.output_settings["phase"]))
                        if lineEdit == 'datafile_StartNum':
                            self.output_object.lineEdit.setText(str(self.set_cl.output_settings["datafile_StartNum"]))
                        if lineEdit == 'datafile_EndNum':
                            self.output_object.lineEdit.setText(str(self.set_cl.output_settings["datafile_EndNum"]))
                        if lineEdit == 'datafile_NumDigit':
                            self.output_object.lineEdit.setText(str(self.set_cl.output_settings["datafile_NumDigit"]))
                        if lineEdit == 'Output_NumAziWrite':
                            self.output_object.lineEdit.setText(str(self.set_cl.output_settings["Output_NumAziWrite"]))      
                        
                # Required Params
                    for wid in range (len(cpf.output_formatters.WritePolydefix.Requirements()[0])):
                        lineEdit = cpf.output_formatters.WritePolydefix.Requirements()[0][wid]
                        self.output_object.lineEdit.setText(cpf.output_formatters.WritePolydefix.Requirements()[0][wid])        
                                 
            elif 'CoefficientTable' == self.set_cl.output_types[output]:
                    self.output_object.Output_Type_comboBox.setCurrentText('WriteCoefficientTable')
                
                # Optional Params                    
                    for wid in range (len(cpf.output_formatters.WriteCoefficientTable.Requirements()[1])):
                        lineEdit = cpf.output_formatters.WriteCoefficientTable.Requirements()[1][wid]
                        if lineEdit == 'coefs_vals_write':
                            self.output_object.lineEdit.setText(str(self.set_cl.output_settings["coefs_vals_write"]))
                        
                # Required Params                        
                    for wid in range (len(cpf.output_formatters.WriteCoefficientTable.Requirements()[0])):
                        lineEdit = cpf.output_formatters.WriteCoefficientTable.Requirements()[0][wid]
                        self.output_object.lineEdit.setText(cpf.output_formatters.WriteCoefficientTable.Requirements()[0][wid])
                        
            elif 'DifferentialStrain' == self.set_cl.output_types[output]:
                    self.output_object.Output_Type_comboBox.setCurrentText('WriteDifferentialStrain')
                
                # Optional Params                    
                    for wid in range (len(cpf.output_formatters.WriteDifferentialStrain.Requirements()[1])):
                            lineEdit = cpf.output_formatters.WriteDifferentialStrain.Requirements()[1][wid]
                            if lineEdit == 'Output_ElasticProperties':
                                self.output_object.lineEdit.setText(str(self.set_cl.output_settings["Output_ElasticProperties"]))
                            elif lineEdit == 'phase':
                                self.output_object.lineEdit.setText(str(self.set_cl.output_settings["phase"]))
                            elif lineEdit == 'datafile_StartNum':
                                self.output_object.lineEdit.setText(str(self.set_cl.output_settings["datafile_StartNum"]))
                            elif lineEdit == 'datafile_EndNum':
                                self.output_object.lineEdit.setText(str(self.set_cl.output_settings["datafile_EndNum"]))
                            elif lineEdit == 'datafile_NumDigit':
                                self.output_object.lineEdit.setText(str(self.set_cl.output_settings["datafile_NumDigit"]))
                            elif lineEdit == 'Output_NumAziWrite':
                                self.output_object.lineEdit.setText(str(self.set_cl.output_settings["Output_NumAziWrite"]))
                
                # Required Params                          
                    for wid in range (len(cpf.output_formatters.WriteDifferentialStrain.Requirements()[0])):
                        lineEdit = cpf.output_formatters.WriteDifferentialStrain.Requirements()[0][wid]
                        self.output_object.lineEdit.setText(cpf.output_formatters.WriteDifferentialStrain.Requirements()[0][wid])
                        
            elif 'MultiFit' == self.set_cl.output_types[output]:
                    self.output_object.Output_Type_comboBox.setCurrentText('WriteMultiFit')
    
                # Optional Params                    
                    for wid in range (len(cpf.output_formatters.WriteMultiFit.Requirements()[1])):
                        lineEdit = cpf.output_formatters.WriteMultiFit.Requirements()[1][wid]
                        if lineEdit == 'Output_ElasticProperties':
                            self.output_object.lineEdit.setText(str(self.set_cl.output_settings["Output_ElasticProperties"]))
                        elif lineEdit == 'phase':
                            self.output_object.lineEdit.setText(str(self.set_cl.output_settings["phase"]))
                        elif lineEdit == 'datafile_StartNum':
                            self.output_object.lineEdit.setText(str(self.set_cl.output_settings["datafile_StartNum"]))
                        elif lineEdit == 'datafile_EndNum':
                            self.output_object.lineEdit.setText(str(self.set_cl.output_settings["datafile_EndNum"]))
                        elif lineEdit == 'datafile_NumDigit':
                            self.output_object.lineEdit.setText(str(self.set_cl.output_settings["datafile_NumDigit"]))
                        elif lineEdit == 'Output_NumAziWrite':
                            self.output_object.lineEdit.setText(str(self.set_cl.output_settings["Output_NumAziWrite"]))
                
                # Required Params                        
                    for wid in range (len(cpf.output_formatters.WriteMultiFit.Requirements()[0])):
                        lineEdit = cpf.output_formatters.WriteMultiFit.Requirements()[0][wid]
                        self.output_object.lineEdit.setText(cpf.output_formatters.WriteMultiFit.Requirements()[0][wid])
                        
            elif 'PolydefixED' == self.set_cl.output_types[output]:
                    self.output_object.Output_Type_comboBox.setCurrentText('WritePolydefixED')
                
                # Optional Params                    
                    for wid in range (len(cpf.output_formatters.WritePolydefixED.Requirements()[1])):
                      
                            lineEdit = cpf.output_formatters.WriteMultiFit.Requirements()[1][wid]
                            if lineEdit == 'ElasticProperties':
                                self.output_object.lineEdit.setText(str(self.set_cl.output_settings["ElasticProperties"]))
                            elif lineEdit == 'phase':
                                self.output_object.lineEdit.setText(str(self.set_cl.output_settings["phase"]))
                            elif lineEdit == 'tc':
                                self.output_object.lineEdit.setText(str(self.set_cl.output_settings["tc"]))
                            elif lineEdit == 'Output_tc':
                                self.output_object.lineEdit.setText(str(self.set_cl.output_settings["Output_tc"]))
                            elif lineEdit == 'Output_TemperaturePower':
                                self.output_object.lineEdit.setText(str(self.set_cl.output_settings["Output_TemperaturePower"]))
                            
                # Required Params                        
                    for wid in range (len(cpf.output_formatters.WritePolydefixED.Requirements()[0])):
                        lineEdit = cpf.output_formatters.WritePolydefixED.Requirements()[0][wid]
                        self.output_object.lineEdit.setText(cpf.output_formatters.WritePolydefixED.Requirements()[0][wid])
                        
            else:
                childcount = self.gridLayout_3.count()
                if childcount >=1:
                    for i in range (0,childcount):
                        item = self.gridLayout_3.itemAt(i)
                        widget = item.widget()
                        widget.deleteLater()
                childcount2 = self.gridLayout_2.count()
                if childcount2 >=1:
                    for i in range (0,childcount2):
                        item = self.gridLayout_2.itemAt(i)
                        widget = item.widget()
                        widget.deleteLater()
       
        with open("../logs/logs.log",'w') as file:
           file.close()
        cpf.XRD_FitPattern.initiate(f"{self.input_file_path}")
        text=open('../logs/logs.log').read()
        self.Console_output.setText(text)
        
        self.Directory.setCursorPosition(0);
        self.Basename.setCursorPosition(0);
        self.Extension.setCursorPosition(0);
        self.Start_Num.setCursorPosition(0);
        self.End_Num.setCursorPosition(0);
        self.Num_Digit.setCursorPosition(0);
        self.Step.setCursorPosition(0);
        self.Calib_Detect.setCursorPosition(0);
        self.Calib_Param.setCursorPosition(0);
        self.Calib_Mask.setCursorPosition(0);
        self.Calib_Pixels.setCursorPosition(0);
        self.Calib_Data.setCursorPosition(0);
        self.cascade_per_bin.setCursorPosition(0);
        self.cascade_number_bins.setCursorPosition(0);
        self.cascade_track.setCursorPosition(0);
        self.bg_1.setCursorPosition(0);
        self.bg_2.setCursorPosition(0);
        self.ds_1.setCursorPosition(0);
        self.ds_2.setCursorPosition(0);
        self.h_1.setCursorPosition(0);
        self.h_2.setCursorPosition(0);
        self.pro_1.setCursorPosition(0);
        self.pro_2.setCursorPosition(0);
        self.wdt_1.setCursorPosition(0);
        self.wdt_2.setCursorPosition(0);
        self.AziBins.setCursorPosition(0);
        self.Output_Dir_1.setCursorPosition(0);
        self.Output_Dir_2.setCursorPosition(0);
    
    @pyqtSlot()
    def Validate_Inputs(self):
        if self.input_file_path == None:
            mess = QMessageBox()
            mess.setIcon(QMessageBox.Warning)
            mess.setText("Please Load input file")
            mess.setStandardButtons(QMessageBox.Ok)
            mess.setWindowTitle("MessageBox")
            returnValue = mess.exec_()
        else:
            with open("../logs/logs.log",'a') as file:
               file.close()
            cpf.XRD_FitPattern.initiate(f"{self.input_file_path}")
            text=open('../logs/logs.log').read()
            self.Console_output.setText(text)
            
    @pyqtSlot()             
    def Run_Range(self):
        #self.matplotlib_inline = matplotlib_inline()
        self.matplotlib_qt = matplotlib_qt()
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
    
    @pyqtSlot()
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
    
    @pyqtSlot()
    def Execute_Fits(self):
        self.matplotlib_qt = matplotlib_qt()
        if self.input_file_path == None:
            mess = QMessageBox()
            mess.setIcon(QMessageBox.Warning)
            mess.setText("Please Load input file")
            mess.setStandardButtons(QMessageBox.Ok)
            mess.setWindowTitle("MessageBox")
            returnValue = mess.exec_()
        else:
            cpf.XRD_FitPattern.execute(f"{self.input_file_path}", save_all=True, parallel = False)
            text=open('../logs/logs.log').read()
            # FIX ME Simon: parallel = True doesn't work on Adina's Windows computer
            self.Console_output.setText(text)
    
    @pyqtSlot()
    def Make_Outputs(self):
        self.matplotlib_qt = matplotlib_qt()
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
    
    @pyqtSlot()
    def Insert_Button(self):  
    # Variables not defined in settings class
        # start_num = self.Start_Num.text()
        # end_num = self.End_Num.text()
        # num_digit = self.Num_Digit.text()
        # step = self.Step.text()
     
        self.set_cl.datafile_directory = self.Directory.text()
        self.set_cl.datafile_basename = self.Basename.text()
        self.set_cl.calibration_type = self.Calib_Type.currentText()
        self.set_cl.calibration_detector = self.Calib_Detect.text()
        self.set_cl.calibration_parameters = self.Calib_Param.text()
        self.set_cl.calibration_mask = self.Calib_Mask.text();
        self.set_cl.calibration_pixel_size = self.Calib_Pixels.text()   
        self.set_cl.calibration_data = self.Calib_Data.text()
        self.set_cl.cascade_bin_type = self.cascade_bin_type.text()
        self.set_cl.cascade_per_bin = self.cascade_per_bin.text()
        self.set_cl.cascade_number_bins = self.cascade_number_bins.text()
        self.set_cl.cascade_track = self.cascade_track.text()
        
        self.set_cl.fit_bounds["background"][0] = self.bg_1.text()
        self.set_cl.fit_bounds["background"][1] = self.bg_2.text()
        self.set_cl.fit_bounds["d-space"][0] = self.ds_1.text()
        self.set_cl.fit_bounds["d-space"][1] = self.ds_2.text()
        self.set_cl.fit_bounds["height"][0] = self.h_1.text()
        self.set_cl.fit_bounds["height"][1] = self.h_2.text()
        self.set_cl.fit_bounds["profile"][0] = self.pro_1.text()
        self.set_cl.fit_bounds["profile"][1] = self.pro_2.text()
        self.set_cl.fit_bounds["width"][0] = self.wdt_1.text()
        self.set_cl.fit_bounds["width"][1] = self.wdt_2.text()
        self.set_cl.output_directory = self.Output_Dir_1.text()
        self.range_length = self.Range_Tab.count()
        
        # Save range_tab data 
        for range_tab in range(self.range_length):
            self.range_object = Range()
            
            self.set_cl.fit_orders[range_tab]["range"][0] = self.range_object.Range_min.text()
            self.set_cl.fit_orders[range_tab]["range"][1] = self.range_object.Range_max.text()
            self.set_cl.fit_orders[range_tab]["background"] = self.range_object.Range_Background_Val.text()
            
            self.set_cl.fit_orders[range_tab]["Imax"] = self.range_object.Intensity_max.text() 
            self.set_cl.fit_orders[range_tab]["Imin"] = self.range_object.Intensity_min.text()
            self.set_cl.fit_orders[range_tab]["PeakPositionSelection"] = self.range_object.Peak_Pos_Selection.toPlainText()
            
            self.peak_length = self.range_object.Peak_Tab.count()
            
            # Save peak_tab data 
            for peak_tab in range(self.peak_length):
                self.peak_object = Peak()
               
                self.set_cl.fit_orders[range_tab]["peak"][peak_tab]["phase"] = self.peak_object.phase_peak.text()
                self.set_cl.fit_orders[range_tab]["peak"][peak_tab]["hkl"] = self.peak_object.hkl.text()
                self.set_cl.fit_orders[range_tab]["peak"][peak_tab]["d-space"] = self.peak_object.d_space_peak.text()
                self.set_cl.fit_orders[range_tab]["peak"][peak_tab]["d-space-type"] = self.peak_object.d_space_type.currentText()
                
                self.set_cl.fit_orders[range_tab]["peak"][peak_tab]["height-type"] = self.peak_object.height_peak_type.currentText()
                self.set_cl.fit_orders[range_tab]["peak"][peak_tab]["profile-type"] = self.peak_object.profile_peak_type.currentText()
                self.set_cl.fit_orders[range_tab]["peak"][peak_tab]["width-type"] = self.peak_object.width_peak_type.currentText()
                
                self.set_cl.fit_orders[range_tab]["peak"][peak_tab]["height"] = self.peak_object.height_peak.text()
                self.set_cl.fit_orders[range_tab]["peak"][peak_tab]["profile"] = self.peak_object.profile_peak.text()
                self.set_cl.fit_orders[range_tab]["peak"][peak_tab]["profile_fixed"] = self.peak_object.profile_fixed.text()
                
                self.set_cl.fit_orders[range_tab]["peak"][peak_tab]["width"] = self.peak_object.width_peak.text()
                self.set_cl.fit_orders[range_tab]["peak"][peak_tab]["symmetry"] = self.peak_object.symmetry_peak.text()

                if  self.peak_object.profile_checkBox.isChecked:
                    self.peak_object.profile_fixed.setEnabled(True)
                    self.set_cl.fit_orders[range_tab]["peak"][peak_tab]["profile_fixed"]= self.peak_object.profile_fixed.text() 
                else:
                    self.peak_object.profile_fixed.setEnabled(False)
                    self.set_cl.fit_orders[range_tab]["peak"][peak_tab]["profile_fixed"] = None           
                    
                if  self.peak_object.width_checkBox.isChecked: 
                    self.peak_object.width_fixed.setEnabled(True)
                    self.set_cl.fit_orders[range_tab]["peak"][peak_tab]["width_fixed"] = self.peak_object.width_fixed.text()  
                else:
                    self.peak_object.width_fixed.setEnabled(False)
                    self.set_cl.fit_orders[range_tab]["peak"][peak_tab]["width_fixed"] = None
                    
                if  self.peak_object.height_checkBox.isChecked:
                    self.peak_object.height_fixed.setEnabled(True)
                    self.set_cl.fit_orders[range_tab]["peak"][peak_tab]["height_fixed"] = self.peak_object.height_fixed.text()
                else:
                    self.peak_object.height_fixed.setEnabled(False)
                    self.set_cl.fit_orders[range_tab]["peak"][peak_tab]["height_fixed"] = None  
                    
                if  self.peak_object.dspace_checkBox.isChecked:
                    self.peak_object.dspace_fixed.setEnabled(True)
                    self.set_cl.fit_orders[range_tab]["peak"][peak_tab]["d-space_fixed"]= self.peak_object.dspace_fixed.text()
                else:
                    self.peak_object.dspace_fixed.setEnabled(False)
                    self.set_cl.fit_orders[range_tab]["peak"][peak_tab]["d-space_fixed"]= None
        
        self.output_length = self.Output_Tab.count()
        
        # Save output_tab data 
        for output_tab in range(0, self.output_length):
            self.output_obj = Output()
            if self.output_obj.Output_Type_comboBox.currentText() =='WriteCoefficientTable':
                self.set_cl.output_types[output_tab] = self.output_obj.Output_Type_comboBox.currentText()
                
                if self.output_obj.coefs_vals_write.text() == '':
                    self.set_cl.output_settings["coefs_vals_write"] = None
                else:
                    self.set_cl.output_settings["coefs_vals_write"] = self.output_obj.coefs_vals_write.text()
            
            elif self.output_obj.Output_Type_comboBox.currentText() =='WriteDifferentialStrain':
                self.set_cl.output_types[output_tab] = self.output_obj.Output_Type_comboBox.currentText()
             
            elif self.output_obj.Output_Type_comboBox.currentText() =='WriteMultiFit':
                self.set_cl.output_types[output_tab] = self.output_obj.Output_Type_comboBox.currentText()
                
                if self.output_obj.Output_NumAziWrite.text() == '':
                    self.set_cl.output_settings["Output_NumAziWrite"] = None
                else:
                    self.set_cl.output_settings["Output_NumAziWrite"] = self.output_obj.Output_NumAziWrite.text() 
           
            elif self.output_obj.Output_Type_comboBox.currentText() =='WritePolydefix':
                self.set_cl.output_types[output_tab] = self.output_obj.Output_Type_comboBox.currentText()
                
                if self.output_obj.Output_NumAziWrite.text() == '':
                    self.set_cl.output_settings["Output_NumAziWrite"] = None
                else:
                    self.set_cl.output_settings["Output_NumAziWrite"] = self.output_obj.Output_NumAziWrite.text()
                    
                if self.output_obj.Output_ElasticProperties.text() == '':
                    self.set_cl.output_settings["Output_ElasticProperties"] = None
                else:
                    self.set_cl.output_settings["Output_ElasticProperties"] = self.output_obj.Output_ElasticProperties.text() 
                
                if self.output_obj.tc.text() == '':
                    self.set_cl.output_settings["tc"] = None
                else:
                    self.set_cl.output_settings["tc"] = self.output_obj.tc.text() 
                
                if self.output_obj.phase.text() == '':
                    self.set_cl.output_settings["phase"] = None
                else:
                    self.set_cl.output_settings["phase"] = self.output_obj.phase.text() 
                    
                if self.output_obj.datafile_StartNum.text() == '':
                     self.set_cl.output_settings["datafile_StartNum"] = None
                else:
                     self.set_cl.output_settings["datafile_StartNum"] = self.output_obj.datafile_StartNum.text() 
                
                if self.output_obj.datafile_EndNum.text() == '':
                     self.set_cl.output_settings["datafile_EndNum"] = None
                else:
                     self.set_cl.output_settings["datafile_EndNum"] = self.output_obj.datafile_EndNum.text() 
                     
                if self.output_obj.datafile_NumDigit.text() == '':
                     self.set_cl.output_settings["datafile_NumDigit"] = None
                else:
                     self.set_cl.output_settings["datafile_NumDigit"] = self.output_obj.datafile_NumDigit.text() 
                
            elif self.output_obj.Output_Type_comboBox.currentText() =='WritePolydefixED':
                self.set_cl.output_types[output_tab] = self.output_obj.Output_Type_comboBox.currentText()
                
                if self.output_obj.tc.text() == '':
                    self.set_cl.output_settings["tc"] = None
                else:
                    self.set_cl.output_settings["tc"] = self.output_obj.tc.text() 
                
                if self.output_obj.phase.text() == '':
                    self.set_cl.output_settings["phase"] = None
                else:
                    self.set_cl.output_settings["phase"] = self.output_obj.phase.text()
                
                if self.output_obj.Output_TemperaturePower.text() == '':
                     self.set_cl.output_settings["Output_TemperaturePower"] = None
                else:
                     self.set_cl.output_settings["Output_TemperaturePower"] = self.output_obj.Output_TemperaturePower.text() 
                     
                if self.output_obj.Output_tc.text() == '':
                     self.set_cl.output_settings["Output_tc"] = None
                else:
                     self.set_cl.output_settings["Output_tc"] = self.output_obj.Output_tc.text() 
                
                if self.output_obj.ElasticProperties.text() == '':
                     self.set_cl.output_settings["ElasticProperties"] = None
                else:
                     self.set_cl.output_settings["ElasticProperties"] = self.output_obj.ElasticProperties.text() 

        self.set_cl.save_settings()
        
    @pyqtSlot()
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

                # complain about selecting a non existing path
            def accept():
                    if checkLineEdit(lineEdit.text()):
                        # if the path is acceptable, call the base accept() implementation
                        QtWidgets.QDialog.accept(dialog)
            dialog.accept = accept

            if dialog.exec_() and dialog.selectedFiles():
                    path = QtCore.QFileInfo(dialog.selectedFiles()[0]).absoluteFilePath()
                    self.Directory.setText(path)
            self.Directory.setCursorPosition(0);
    
    @pyqtSlot()
    def selectTarget(self):
            dialog = QtWidgets.QFileDialog()

            if self.Output_Dir_1.text():
                    dialog.setDirectory(self.Output_Dir_1.text())

            dialog.setFileMode(dialog.Directory)
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

            button = dialog.findChild(QtWidgets.QDialogButtonBox).button(
                    QtWidgets.QDialogButtonBox.Open)

            lineEdit = dialog.findChild(QtWidgets.QLineEdit)
            lineEdit.textChanged.connect(checkLineEdit)
  
            def accept():
                    if checkLineEdit(lineEdit.text()):
                        QtWidgets.QDialog.accept(dialog)
            dialog.accept = accept

            if dialog.exec_() and dialog.selectedFiles():
                    path = QtCore.QFileInfo(dialog.selectedFiles()[0]).absoluteFilePath()
                    self.Output_Dir_1.setText(path)
                    self.Output_Dir_2.setText(path)
            self.Output_Dir_1.setCursorPosition(0);
            self.Output_Dir_2.setCursorPosition(0);
    
    @pyqtSlot()
    def Insert_Range(self):
        self.Range_Tab.addTab(Range() , QIcon(""),"Range")
    
    @pyqtSlot()
    def Remove_Range(self):
        self.Range_Tab.removeTab(self.Range_Tab.currentIndex())
    
    @pyqtSlot()
    def Insert_Output(self):
        self.Output_Tab.addTab(Output() , QIcon(""),"Output")
        
    @pyqtSlot()
    def Remove_Output(self):
        self.Output_Tab.removeTab(self.Output_Tab.currentIndex())
        
    @pyqtSlot()
    def Dir_Pressed(self):
        datafile_directory = self.Directory.text()
        self.set_cl.datafile_directory = datafile_directory
        
    @pyqtSlot()
    def Basename_Pressed(self):
        datafile_basename = self.Basename.text()
        self.set_cl.datafile_basename = datafile_basename
        
    @pyqtSlot()
    def Calib_Detect_Pressed(self):
        calibration_detector = self.Calib_Detect.text()
        self.set_cl.calibration_detector = calibration_detector
        
    @pyqtSlot()
    def Calib_Param_Pressed(self):
        calibration_parameters = self.Calib_Param.text()
        self.set_cl.calibration_parameters = calibration_parameters
        
    @pyqtSlot()
    def Calib_Mask_Pressed(self):
        calibration_mask = self.Calib_Mask.text()
        self.set_cl.calibration_mask = calibration_mask
    
    @pyqtSlot()
    def Calib_Pixels_Pressed(self):
        calibration_pixel_size = self.Calib_Pixels.text()
        self.set_cl.calibration_pixel_size = calibration_pixel_size
    
    @pyqtSlot()
    def Calib_Data_Pressed(self):
        calibration_data = self.Calib_Data.text()
        self.set_cl.calibration_pixel_size = calibration_data
    
    @pyqtSlot()
    def Cascade_bin_type_Pressed(self):
        cascade_bin_type = self.cascade_bin_type.text()
        self.set_cl.cascade_bin_type = cascade_bin_type
    
    @pyqtSlot()
    def Cascade_per_bin_Pressed(self):
        cascade_per_bin = self.cascade_per_bin.text()
        self.set_cl.cascade_per_bin = cascade_per_bin
    
    @pyqtSlot()
    def Cascade_number_bins_Pressed(self):
        cascade_number_bins = self.cascade_number_bins.text()
        self.set_cl.cascade_number_bins = cascade_number_bins
    
    @pyqtSlot()
    def Cascade_track_Pressed(self):
        cascade_track = self.cascade_track.text()
        self.set_cl.cascade_track = cascade_track
    
    @pyqtSlot()
    def Bg_1_Pressed(self):
        bg_1 = self.bg_1.text()
        self.set_cl.fit_bounds.get("background")[0] = bg_1
    
    @pyqtSlot()
    def Bg_2_Pressed(self):
        bg_2 = self.bg_2.text()
        self.set_cl.fit_bounds.get("background")[1] = bg_2
        
    @pyqtSlot()
    def Ds_1_Pressed(self):
        ds_1 = self.ds_1.text()
        self.set_cl.fit_bounds.get("d-space")[0] = ds_1
    
    @pyqtSlot()
    def Ds_2_Pressed(self):
        ds_2 = self.ds_2.text()
        self.set_cl.fit_bounds.get("d-space")[1] = ds_2
    
    @pyqtSlot()
    def H_1_Pressed(self):
        h_1 = self.h_1.text()
        self.set_cl.fit_bounds.get("height")[0] = h_1
    
    @pyqtSlot()
    def H_2_Pressed(self):
        h_2 = self.h_2.text()
        self.set_cl.fit_bounds.get("height")[1] = h_2
        
    @pyqtSlot()
    def Pro_1_Pressed(self):
        pro_1 = self.pro_1.text()
        self.set_cl.fit_bounds.get("profile")[0] = pro_1
        
    @pyqtSlot()
    def Pro_2_Pressed(self):
        pro_2 = self.pro_2.text()
        self.set_cl.fit_bounds.get("profile")[1] = pro_2
    
    @pyqtSlot()
    def Wdt_1_Pressed(self):
        wdt_1 = self.wdt_1.text()
        self.set_cl.fit_bounds.get("width")[0] = wdt_1
        
    @pyqtSlot()
    def Wdt_2_Pressed(self):
        wdt_2 = self.wdt_2.text()
        self.set_cl.fit_bounds.get("width")[1] = wdt_2
        
    @pyqtSlot()
    def Output_Dir_1_Pressed(self):
        Output_Dir_1 = self.Output_Dir_1.text()
        self.set_cl.output_directory = Output_Dir_1
        
    def closeEvent(self,event):
        os._exit(00)

def main():
    app = QApplication(sys.argv)
    mainwindow = MainWindow()
    mainwindow.show() 
    try:        
      sys.exit(app.exec_())
    except:
      os._exit(00)

if __name__=='__main__':
    main()
    
