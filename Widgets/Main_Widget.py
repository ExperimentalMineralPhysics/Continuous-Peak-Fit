import sys
from PyQt6 import (
    QtWidgets, 
    QtCore, 
    QtGui
    )
from PyQt6.QtWidgets import (
    QMainWindow, 
    QApplication,  
    QFileDialog
    )
from PyQt6.QtGui import QIcon
from PyQt6.QtCore import pyqtSlot
from PyQt6.uic import loadUi
from PyQt6.QtWidgets import QMessageBox

from Widgets.Range_Widget import Range
from Widgets.Peak_Widget import Peak
from Widgets.Output_Widget import Output

import cpf
from cpf.settings import settings

from matplotlib_qt import matplotlib_qt
from matplotlib_auto import matplotlib_inline

import os


class Main_Widget(QMainWindow):
    
    def __init__(self):
        super(Main_Widget, self).__init__()
        loadUi("Widgets/Main_Widget.ui", self)
        self.setWindowIcon(QtGui.QIcon('logo.png'))
        self.gui_layout()
        
        # Define global variables
        self.set_cl = cpf.settings.settings()
        self.input_file_path = None  
        
        self.range_list = []
        self.output_list = []
        
    def gui_layout(self):
        self.Main_Tab.setMinimumHeight(30);
        self.Directory.setMinimumHeight(30);
        
        self.Basename.setMinimumHeight(30);
        self.Extension.setMinimumHeight(30);
        
        self.Start_Num.setMinimumHeight(30);
        self.End_Num.setMinimumHeight(30);
        
        self.Num_Digit.setMinimumHeight(30);
        self.Step.setMinimumHeight(30);
        
        self.Calib_Type.setMinimumHeight(30);
        self.Calib_Detect.setMinimumHeight(30);
        
        self.Calib_Param.setMinimumHeight(30);
        self.Calib_Mask.setMinimumHeight(30);
        
        self.Calib_Pixels.setMinimumHeight(30);
        self.Calib_Data.setMinimumHeight(30);
        
        self.cascade_bin_type.setMinimumHeight(30);
        self.cascade_per_bin.setMinimumHeight(30);
        
        self.cascade_number_bins.setMinimumHeight(30);
        self.cascade_track.setMinimumHeight(30);
        
        self.bg_1.setMinimumHeight(30);
        self.bg_2.setMinimumHeight(30);
        
        self.ds_1.setMinimumHeight(30);
        self.ds_2.setMinimumHeight(30);
        
        self.h_1.setMinimumHeight(30);
        self.h_2.setMinimumHeight(30);
        
        self.pro_1.setMinimumHeight(30);
        self.pro_2.setMinimumHeight(30);
        
        self.wdt_1.setMinimumHeight(30);
        self.wdt_2.setMinimumHeight(30);
        
        self.AziBins.setMinimumHeight(30);
        self.Output_Dir_1.setMinimumHeight(30);
        
        self.Output_Dir_2.setMinimumHeight(30);
        self.Console_output.setReadOnly(True)
        
        self.AddRange_Btn.clicked.connect(self.Insert_Range)
        self.RemoveRange_Btn.clicked.connect(self.Remove_Range)
        
        self.Add_Output_Btn.clicked.connect(self.Insert_Output)
        self.Remove_Output_Btn.clicked.connect(self.Remove_Output)
        
        self.Select_Directory_1.clicked.connect(self.select_Data_Dir)
        self.Select_Directory_2.clicked.connect(self.select_Output_Dir)
        
        self.Save_input_Btn.clicked.connect(self.Insert_Button)
        self.Load_input_Btn.clicked.connect(self.Load_Inputs)
        
        self.validate_Btn.clicked.connect(self.Validate_Inputs)
        self.Run_Range_Btn.clicked.connect(self.Run_Range)
        
        self.Run_initial_Btn.clicked.connect(self.Run_Initial_Position)
        self.Execute_Btn.clicked.connect(self.Execute_Fits)
        
        self.Make_Output_Btn.clicked.connect(self.Make_Outputs)
        
    @pyqtSlot()
    def Load_Inputs(self):
        '''
        whenever new file is loaded 
        clear settings' class existing data
        otherwise it will show values already 
        stored in variables by repacing None
        with values

        ''' 
        self.set_cl.datafile_directory = "."
        self.set_cl.datafile_basename = None
  
        self.set_cl.datafile_number = 0
        self.set_cl.image_list = None
        self.set_cl.image_number = 0
        
        self.set_cl.datafile_preprocess = None
        self.set_cl.calibration_type = None
        
        self.set_cl.calibration_detector = None
        self.set_cl.calibration_parameters = None
        
        self.set_cl.calibration_mask = None
        self.set_cl.calibration_pixel_size = None 
        
        self.set_cl.calibration_data = None
        self.set_cl.cascade_bin_type = None
        
        self.set_cl.cascade_per_bin = None
        self.set_cl.cascade_number_bins = None
    
        self.set_cl.cascade_track = False
        
        self.set_cl.output_directory = "."
        self.set_cl.fit_orders = None
        
        self.set_cl.fit_bin_type = None
        self.set_cl.fit_per_bin = None
        
        self.set_cl.fit_number_bins = None
        self.set_cl.fit_orders = None
        
        self.set_cl.fit_bounds = {
            "background": ["0.95*min", "1.05*max"],
            "d-space": ["min", "max"],
            "height": [0, "1.05*max"],
            "profile": [0, 1],
            "width": ["range/(ndata)", "range/2"],
        }

        self.set_cl.fit_track = False
        self.set_cl.fit_propagate = True
        
        self.set_cl.output_types = None
        self.set_cl.output_settings = dict()
        
        self.set_cl.subfit_file_position = None
        self.set_cl.subfit_filename = None
        
        self.set_cl.subfit_order_position = None
        self.set_cl.subfit_orders = None

        fname= QFileDialog.getOpenFileName(self, "Load Input File", "./", "Python Files (*.py)")
        
        if fname:
            self.input_file_path = f"{fname[0]}"
        self.set_cl.populate(settings_file=(f"{self.input_file_path}"))
        
        # update new data
        self.Directory.setText(self.set_cl.datafile_directory)
        self.Basename.setText(self.set_cl.datafile_basename)
        
        self.Calib_Type.setCurrentText(f"{self.set_cl.calibration_type}")
        self.Calib_Detect.setText(self.set_cl.calibration_detector)
        
        self.Calib_Param.setText(self.set_cl.calibration_parameters)
        self.Calib_Mask.setText(self.set_cl.calibration_mask)
        
        self.Calib_Pixels.setText(str(self.set_cl.calibration_pixel_size))   
        self.Calib_Data.setText(self.set_cl.calibration_data)
        
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
            
            range_object = Range()
            
            self.Range_Tab.addTab(range_object , QIcon(""),"Range")
            range_object.Range_min.setText(str(self.set_cl.fit_orders[range_tab].get("range")[0]))
            
            range_object.Range_max.setText(str(self.set_cl.fit_orders[range_tab].get("range")[1]))
            
            range_object.Range_Background_Val.setText(str(self.set_cl.fit_orders[range_tab].get("background")))
            range_object.Intensity_max.setText(str(self.set_cl.fit_orders[range_tab].get("Imax")))   
            
            range_object.Intensity_min.setText(str(self.set_cl.fit_orders[range_tab].get("Imin")))   
            range_object.Peak_Pos_Selection.setPlainText(str(self.set_cl.fit_orders[range_tab].get("PeakPositionSelection")))
            
            if self.set_cl.fit_orders[range_tab].get("background-type")== None:
                range_object.bg_fixed_checkbox.setChecked(False)
            else:
                range_object.bg_fixed_checkbox.setChecked(True)
                self.b_type = str(self.set_cl.fit_orders[range_tab].get("background-type"))
                range_object.Background_Type.setCurrentText(f"{self.b_type}")

            self.range_list.append(range_object) 

            self.peak_length = len(self.set_cl.fit_orders[range_tab]["peak"])
            
            # Manage peak_tab data population
            for peak_tab in range(0, self.peak_length):
                
                peak_object = Peak()
                
                range_object.Peak_Tab.addTab(peak_object , QIcon(""),"Peak")
                peak_object.phase_peak.setText(str(self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("phase")))
                
                peak_object.hkl.setText(str(self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("hkl")))
                peak_object.d_space_peak.setText(str(self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("d-space")))
                
                self.ds_type = str(self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("d-space-type"))
                
                peak_object.d_space_type.setCurrentText(f"{self.ds_type}")
                self.h_type = str(self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("height-type"))
                
                peak_object.height_peak_type.setCurrentText(f"{self.h_type}")
               
                self.p_type = str(self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("profile-type"))
                peak_object.profile_peak_type.setCurrentText(f"{self.p_type}")
                
                self.w_type = str(self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("width-type"))
               
                peak_object.width_peak_type.setCurrentText(f"{self.w_type}")
                peak_object.height_peak.setText(str(self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("height")))
               
                peak_object.profile_peak.setText(str(self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("profile")))
                peak_object.profile_fixed.setText(str(self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("profile_fixed")))
                
                peak_object.width_peak.setText(str(self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("width")))
                peak_object.symmetry_peak.setText(str(self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("symmetry")))

                if self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("profile_fixed")== None:
                    peak_object.profile_checkBox.setChecked(False)
                else:
                    peak_object.dspace_fixed.setText(str(self.set_cl.fit_orders[range_tab].get("profile_fixed")))
                    peak_object.profile_checkBox.setChecked(True)
                    
                if self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("width_fixed")== None:
                    peak_object.width_checkBox.setChecked(False)
                else:
                    peak_object.width_fixed.setText(str(self.set_cl.fit_orders[range_tab].get("width_fixed")))
                    peak_object.width_checkBox.setChecked(True)
                    
                if self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("height_fixed")== None:
                    peak_object.height_checkBox.setChecked(False)
                else:
                    peak_object.width_fixed.setText(str(self.set_cl.fit_orders[range_tab].get("height_fixed")))
                    peak_object.height_checkBox.setChecked(True)
                    
                if self.set_cl.fit_orders[range_tab].get("peak")[peak_tab].get("d-space_fixed")== None:
                    peak_object.dspace_checkBox.setChecked(False)
                else:
                    peak_object.width_fixed.setText(str(self.set_cl.fit_orders[range_tab].get("d-space_fixed")))
                    peak_object.dspace_checkBox.setChecked(True)
               
                range_object.peak_list.append(peak_object)
                    
        # Remove existing output_tabs whenever new file is loaded
        for remove_output in range(0, self.Output_Tab.count()):
            self.Remove_Output()
        
        # Manage output_tab data population    
        self.output_type =  len(self.set_cl.output_types)
        
        for output in range (0, self.output_type):
            
            output_object = Output()
            
            self.Output_Tab.addTab(output_object , QIcon(""),"Output")

            if 'Polydefix' == self.set_cl.output_types[output]:
                    
                    output_object.Output_Type_comboBox.setCurrentText('WritePolydefix')
               
                # Optional Params
                    for wid in output_object.WritePolydefix_optional_list:
                        
                        lineEdit = wid
                        
                        if 'Output_ElasticProperties' in self.set_cl.output_settings:
                            
                            if lineEdit.objectName() == 'Output_ElasticProperties':
                               lineEdit.setText(str(self.set_cl.output_settings["Output_ElasticProperties"]))
                        
                        if 'datafile_StartNum' in self.set_cl.output_settings:
                            
                            if lineEdit.objectName() == 'datafile_StartNum':
                               lineEdit.setText(str(self.set_cl.output_settings["datafile_StartNum"]))
                        
                        if 'datafile_EndNum' in self.set_cl.output_settings:
                            
                            if lineEdit.objectName() == 'datafile_EndNum':
                               lineEdit.setText(str(self.set_cl.output_settings["datafile_EndNum"]))
                        
                        if 'datafile_NumDigit' in self.set_cl.output_settings:
                            
                            if lineEdit.objectName() == 'datafile_NumDigit':
                               lineEdit.setText(str(self.set_cl.output_settings["datafile_NumDigit"]))
                        
                        if 'Output_NumAziWrite' in self.set_cl.output_settings:
                            
                            if lineEdit.objectName() =='Output_NumAziWrite':
                               lineEdit.setText(str(self.set_cl.output_settings["Output_NumAziWrite"]))
                       
                        if 'phase' in self.set_cl.output_settings:
                            
                            if lineEdit.objectName() == 'phase':
                               lineEdit.setText(str(self.set_cl.output_settings["phase"]))
                        
                # Required Params
                    for wid in output_object.WritePolydefix_required_list:
                        
                        list_len = len(output_object.WritePolydefix_required_list)
                        
                        for i in range(list_len):
                           
                            if wid.objectName() == cpf.output_formatters.WritePolydefix.Requirements()[0][i]:
                               wid.setText(cpf.output_formatters.WritePolydefix.Requirements()[0][i])       
                                 
            elif 'CoefficientTable' == self.set_cl.output_types[output]:
                    
                    output_object.Output_Type_comboBox.setCurrentText('WriteCoefficientTable')
                
                # Optional Params                    
                    for wid in output_object.WriteCoefficientTable_optional_list:
                        
                        lineEdit = wid
                        
                        if 'coefs_vals_write' in self.set_cl.output_settings:
                            
                            if lineEdit.objectName() == 'coefs_vals_write':
                                lineEdit.setText(str(self.set_cl.output_settings["coefs_vals_write"]))
                        
                # Required Params                        
                    for wid in output_object.WriteCoefficientTable_required_list:
                        
                        list_len = len(output_object.WriteCoefficientTable_required_list)
                        
                        for i in range(list_len):
                           
                            if wid.objectName() == cpf.output_formatters.WriteCoefficientTable.Requirements()[0][i]:
                               wid.setText(cpf.output_formatters.WriteCoefficientTable.Requirements()[0][i])
                        
            elif 'DifferentialStrain' == self.set_cl.output_types[output]:
                    
                    output_object.Output_Type_comboBox.setCurrentText('WriteDifferentialStrain')
                
                # Optional Params                    
                    for wid in output_object.WriteDifferentialStrain_optional_list:
                            
                            lineEdit = wid
                            
                            if 'Output_ElasticProperties' in self.set_cl.output_settings:
                                
                                if lineEdit.objectName() == 'Output_ElasticProperties':
                                    lineEdit.setText(str(self.set_cl.output_settings["Output_ElasticProperties"]))
                           
                            if 'phase' in self.set_cl.output_settings:
                                
                                if lineEdit.objectName() == 'phase':
                                    lineEdit.setText(str(self.set_cl.output_settings["phase"]))
                           
                            if 'datafile_StartNum' in self.set_cl.output_settings:
                                
                                if lineEdit.objectName() == 'datafile_StartNum':
                                    lineEdit.setText(str(self.set_cl.output_settings["datafile_StartNum"]))
                            
                            if 'datafile_EndNum' in self.set_cl.output_settings:
                               
                                if lineEdit.objectName() == 'datafile_EndNum':
                                    lineEdit.setText(str(self.set_cl.output_settings["datafile_EndNum"]))
                            
                            if 'datafile_NumDigit' in self.set_cl.output_settings:
                               
                                if lineEdit.objectName() == 'datafile_NumDigit':
                                    lineEdit.setText(str(self.set_cl.output_settings["datafile_NumDigit"]))
                            
                            if 'Output_NumAziWrite' in self.set_cl.output_settings:
                                
                                if lineEdit.objectName() == 'Output_NumAziWrite':
                                    lineEdit.setText(str(self.set_cl.output_settings["Output_NumAziWrite"]))
                
                # Required Params                          
                    for wid in output_object.WriteDifferentialStrain_required_list:
                        
                        list_len = len(output_object.WriteDifferentialStrain_required_list)
                        
                        for i in range(list_len):
                           
                            if wid.objectName() == cpf.output_formatters.WriteDifferentialStrain.Requirements()[0][i]:
                               wid.setText(cpf.output_formatters.WriteDifferentialStrain.Requirements()[0][i])  
                        
            elif 'MultiFit' == self.set_cl.output_types[output]:
                    
                    output_object.Output_Type_comboBox.setCurrentText('WriteMultiFit')
    
                # Optional Params                    
                    for wid in output_object.WriteMultiFit_optional_list:
                        
                        lineEdit = wid
                        
                        if 'Output_ElasticProperties' in self.set_cl.output_settings:
                            
                            if lineEdit.objectName() == 'Output_ElasticProperties':
                                lineEdit.setText(str(self.set_cl.output_settings["Output_ElasticProperties"]))
                        
                        if 'phase' in self.set_cl.output_settings:    
                            
                            if lineEdit.objectName() == 'phase':
                                lineEdit.setText(str(self.set_cl.output_settings["phase"]))
                       
                        if 'datafile_StartNum' in self.set_cl.output_settings:    
                          
                            if lineEdit.objectName() == 'datafile_StartNum':
                                lineEdit.setText(str(self.set_cl.output_settings["datafile_StartNum"]))
                        
                        if 'datafile_EndNum' in self.set_cl.output_settings:        
                          
                            if lineEdit.objectName() == 'datafile_EndNum':
                                lineEdit.setText(str(self.set_cl.output_settings["datafile_EndNum"]))
                        
                        if 'datafile_NumDigit' in self.set_cl.output_settings:     
                           
                            if lineEdit.objectName() == 'datafile_NumDigit':
                                lineEdit.setText(str(self.set_cl.output_settings["datafile_NumDigit"]))
                        
                        if 'Output_NumAziWrite' in self.set_cl.output_settings:   
                           
                            if lineEdit.objectName() == 'Output_NumAziWrite':
                                lineEdit.setText(str(self.set_cl.output_settings["Output_NumAziWrite"]))
                    
                # Required Params                        
                    for wid in output_object.WriteMultiFit_required_list:
                        
                        list_len = len(output_object.WriteMultiFit_required_list)
                        
                        for i in range(list_len):
                           
                            if wid.objectName() == cpf.output_formatters.WriteMultiFit.Requirements()[0][i]:
                               wid.setText(cpf.output_formatters.WriteMultiFit.Requirements()[0][i])
                        
            elif 'PolydefixED' == self.set_cl.output_types[output]:
                    
                    output_object.Output_Type_comboBox.setCurrentText('WritePolydefixED')
                
                # Optional Params                    
                    for wid in output_object.WritePolydefixED_optional_list:
                            
                            lineEdit = wid
                            
                            if 'ElasticProperties' in self.set_cl.output_settings:
                              
                                if lineEdit.objectName() == 'ElasticProperties':
                                    lineEdit.setText(str(self.set_cl.output_settings["ElasticProperties"]))
                            
                            if 'phase' in self.set_cl.output_settings:    
                               
                                if lineEdit.objectName() == 'phase':
                                    lineEdit.setText(str(self.set_cl.output_settings["phase"]))
                           
                            if 'tc' in self.set_cl.output_settings:   
                              
                                if lineEdit.objectName() == 'tc':
                                    lineEdit.setText(str(self.set_cl.output_settings["tc"]))
                           
                            if 'Output_tc' in self.set_cl.output_settings:     
                                
                                if lineEdit.objectName() == 'Output_tc':
                                    lineEdit.setText(str(self.set_cl.output_settings["Output_tc"]))
                            
                            if 'Output_TemperaturePower' in self.set_cl.output_settings:       
                               
                                if lineEdit.objectName() == 'Output_TemperaturePower':
                                    lineEdit.setText(str(self.set_cl.output_settings["Output_TemperaturePower"]))
                                
                # Required Params                        
                    for wid in output_object.WritePolydefixED_required_list:
                       
                        list_len = len(output_object.WritePolydefixED_required_list)
                        
                        for i in range(list_len):
                           
                            if wid.objectName() == cpf.output_formatters.WritePolydefixED.Requirements()[0][i]:
                               wid.setText(cpf.output_formatters.WritePolydefixED.Requirements()[0][i])
                        
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
               
            self.output_list.append(output_object) 
            
            # clear log file data whenever new file is loaded
            file_to_delete = open(settings.log_file,'w')
            file_to_delete.close()
        
        # set cursor position at the start of text box after filling it
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
            mess.setWindowIcon(QtGui.QIcon('error.JFIF'))
            mess.setIcon(QMessageBox.Icon.Warning)
            mess.setText("Please Load input file")
            mess.setStandardButtons(QMessageBox.StandardButton.Ok)
            mess.setWindowTitle("ALERT")
            returnValue = mess.exec()
        else:
            with open(settings.log_file,'a') as file:
               file.close()
            cpf.XRD_FitPattern.initiate(f"{self.input_file_path}")
            text=open(settings.log_file).read()
            self.Console_output.setText(text)
            
    @pyqtSlot()             
    def Run_Range(self):
        #self.matplotlib_inline = matplotlib_inline()
        self.matplotlib_qt = matplotlib_qt()
        if self.input_file_path == None:
            mess = QMessageBox()
            mess.setWindowIcon(QtGui.QIcon('error.JFIF'))
            mess.setIcon(QMessageBox.Icon.Warning)
            mess.setText("Please Load input file")
            mess.setStandardButtons(QMessageBox.StandardButton.Ok)
            mess.setWindowTitle("ALERT")
            returnValue = mess.exec()
        else:
            cpf.XRD_FitPattern.set_range(f"{self.input_file_path}", save_all=True)
            text=open(settings.log_file).read()
            self.Console_output.setText(text)
    
    @pyqtSlot()
    def Run_Initial_Position(self):
        self.matplotlib_qt = matplotlib_qt()
        QApplication.processEvents()
        if self.input_file_path == None:
            mess = QMessageBox()
            mess.setWindowIcon(QtGui.QIcon('error.JFIF'))
            mess.setIcon(QMessageBox.Icon.Warning)
            mess.setText("Please Load input file")
            mess.setStandardButtons(QMessageBox.StandardButton.Ok)
            mess.setWindowTitle("ALERT")
            returnValue = mess.exec()
        else:
            cpf.XRD_FitPattern.initial_peak_position(f"{self.input_file_path}", save_all=True)
            text=open(settings.log_file).read()
            self.Console_output.setText(text)
            self.matplotlib_inline = matplotlib_inline()
    
    @pyqtSlot()
    def Execute_Fits(self):
        self.matplotlib_qt = matplotlib_qt()
        if self.input_file_path == None:
            mess = QMessageBox()
            mess.setWindowIcon(QtGui.QIcon('error.JFIF'))
            mess.setIcon(QMessageBox.Icon.Warning)
            mess.setText("Please Load input file")
            mess.setStandardButtons(QMessageBox.StandardButton.Ok)
            mess.setWindowTitle("ALERT")
            returnValue = mess.exec()
        else:
            cpf.XRD_FitPattern.execute(f"{self.input_file_path}", save_all=True, parallel = False)
            text=open(settings.log_file).read()
            # FIX ME Simon: parallel = True doesn't work on Adina's Windows computer
            self.Console_output.setText(text)
    
    @pyqtSlot()
    def Make_Outputs(self):
        self.matplotlib_qt = matplotlib_qt()
        if self.input_file_path == None:
            mess = QMessageBox()
            mess.setWindowIcon(QtGui.QIcon('error.JFIF'))
            mess.setIcon(QMessageBox.Icon.Warning)
            mess.setText("Please Load input file")
            mess.setStandardButtons(QMessageBox.StandardButton.Ok)
            mess.setWindowTitle("ALERT")
            returnValue = mess.exec()
        else:
            cpf.XRD_FitPattern.write_output(f"{self.input_file_path}", save_all=True)
            text=open(settings.log_file).read()
            self.Console_output.setText(text)
    
    @pyqtSlot()
    def Insert_Button(self):  
        ''' 
        Variables not defined in settings class
        start_num = self.Start_Num.text()
        end_num = self.End_Num.text()
        num_digit = self.Num_Digit.text()
        step = self.Step.text() '''
        
        # store widgets' data into setting class variables
        self.set_cl.datafile_directory = self.Directory.text()
        self.set_cl.datafile_basename = self.Basename.text()
        
        if self.Calib_Type.currentText() == 'Select Type':
            self.set_cl.calibration_type = None
        else:
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
        
        # Save range_tab data 
        for range_object in range (len(self.range_list)):
            
            self.set_cl.fit_orders[range_object]["range"][0] = self.range_list[range_object].Range_min.text()
            self.set_cl.fit_orders[range_object]["range"][1] = self.range_list[range_object].Range_max.text()
            
            self.set_cl.fit_orders[range_object]["background"] = self.range_list[range_object].Range_Background_Val.text()
            
            if self.range_list[range_object].Background_Type.currentText() == 'Select Type':
                self.set_cl.fit_orders[range_object]["background-type"] = ''
            else:
                self.set_cl.fit_orders[range_object]["background-type"] = self.range_list[range_object].Background_Type.currentText()
            
            self.set_cl.fit_orders[range_object]["Imax"] = self.range_list[range_object].Intensity_max.text() 
            self.set_cl.fit_orders[range_object]["Imin"] = self.range_list[range_object].Intensity_min.text()
            
            self.set_cl.fit_orders[range_object]["PeakPositionSelection"] = self.range_list[range_object].Peak_Pos_Selection.toPlainText()
            
            # Save peak_tab data 
            for peak_object in self.range_list[range_object].peak_list:
                
                self.peak_indices = self.range_list[range_object].peak_list.index(peak_object)
               
                self.set_cl.fit_orders[range_object]["peak"][self.peak_indices]["phase"] = peak_object.phase_peak.text()
                self.set_cl.fit_orders[range_object]["peak"][self.peak_indices]["hkl"] = peak_object.hkl.text()
                
                self.set_cl.fit_orders[range_object]["peak"][self.peak_indices]["d-space"] = peak_object.d_space_peak.text()
                
                if peak_object.d_space_type.currentText() == 'Select Type':
                    self.set_cl.fit_orders[range_object]["peak"][self.peak_indices]["d-space-type"] = ''
                else:
                    self.set_cl.fit_orders[range_object]["peak"][self.peak_indices]["d-space-type"] = peak_object.d_space_type.currentText()
                
                if peak_object.height_peak_type.currentText() == 'Select Type':
                    self.set_cl.fit_orders[range_object]["peak"][self.peak_indices]["height-type"] = ''  
                else: 
                    self.set_cl.fit_orders[range_object]["peak"][self.peak_indices]["height-type"] = peak_object.height_peak_type.currentText()
                
                if peak_object.profile_peak_type.currentText() == 'Select Type':
                    self.set_cl.fit_orders[range_object]["peak"][self.peak_indices]["profile-type"] = ''
                else:
                    self.set_cl.fit_orders[range_object]["peak"][self.peak_indices]["profile-type"] = peak_object.profile_peak_type.currentText()
                
                if peak_object.width_peak_type.currentText() == 'Select Type':
                    self.set_cl.fit_orders[range_object]["peak"][self.peak_indices]["width-type"] = ''
                else:    
                    self.set_cl.fit_orders[range_object]["peak"][self.peak_indices]["width-type"] = peak_object.width_peak_type.currentText()
                
                self.set_cl.fit_orders[range_object]["peak"][self.peak_indices]["height"] = peak_object.height_peak.text()
                self.set_cl.fit_orders[range_object]["peak"][self.peak_indices]["profile"] = peak_object.profile_peak.text()
                
                self.set_cl.fit_orders[range_object]["peak"][self.peak_indices]["profile_fixed"] = peak_object.profile_fixed.text()
                
                self.set_cl.fit_orders[range_object]["peak"][self.peak_indices]["width"] = peak_object.width_peak.text()
                self.set_cl.fit_orders[range_object]["peak"][self.peak_indices]["symmetry"] = peak_object.symmetry_peak.text()

                if  peak_object.profile_checkBox.isChecked:
                    peak_object.profile_fixed.setEnabled(True)
                    self.set_cl.fit_orders[range_object]["peak"][self.peak_indices]["profile_fixed"]= peak_object.profile_fixed.text() 
                else:
                    peak_object.profile_fixed.setEnabled(False)
                    self.set_cl.fit_orders[range_object]["peak"][self.peak_indices]["profile_fixed"] = None           
                    
                if  peak_object.width_checkBox.isChecked: 
                    peak_object.width_fixed.setEnabled(True)
                    self.set_cl.fit_orders[range_object]["peak"][self.peak_indices]["width_fixed"] = peak_object.width_fixed.text()  
                else:
                    peak_object.width_fixed.setEnabled(False)
                    self.set_cl.fit_orders[range_object]["peak"][self.peak_indices]["width_fixed"] = None
                    
                if  peak_object.height_checkBox.isChecked:
                    peak_object.height_fixed.setEnabled(True)
                    self.set_cl.fit_orders[range_object]["peak"][self.peak_indices]["height_fixed"] = peak_object.height_fixed.text()
                else:
                    peak_object.height_fixed.setEnabled(False)
                    self.set_cl.fit_orders[range_object]["peak"][self.peak_indices]["height_fixed"] = None  
                    
                if  peak_object.dspace_checkBox.isChecked:
                    peak_object.dspace_fixed.setEnabled(True)
                    self.set_cl.fit_orders[range_object]["peak"][self.peak_indices]["d-space_fixed"]= peak_object.dspace_fixed.text()
                else:
                    peak_object.dspace_fixed.setEnabled(False)
                    self.set_cl.fit_orders[range_object]["peak"][self.peak_indices]["d-space_fixed"]= None
        
        # Save output_tab data 
        for output_object in self.output_list:
            
            indices = self.output_list.index(output_object)
            
            if output_object.Output_Type_comboBox.currentText() =='WriteCoefficientTable':
                
                self.set_cl.output_types.append('CoefficientTable')
                
                for i in output_object.WriteCoefficientTable_optional_list:
                    
                    if i.objectName() == 'coefs_vals_write':
                        
                        if i.text()=='':
                            self.set_cl.output_settings["coefs_vals_write"] = None
                        else:
                            self.set_cl.output_settings["coefs_vals_write"] = i.text()
            
            elif output_object.Output_Type_comboBox.currentText() =='WriteDifferentialStrain':
                
                self.set_cl.output_types.append('DifferentialStrain')
             
            elif output_object.Output_Type_comboBox.currentText() =='WriteMultiFit':
                
                self.set_cl.output_types.append('MultiFit')
                
                for i in output_object.WriteMultiFit_optional_list:
                   
                    if i.objectName() == 'Output_NumAziWrite':
                        
                        if i.text() == '':
                            self.set_cl.output_settings["Output_NumAziWrite"] = None
                        else:
                            self.set_cl.output_settings["Output_NumAziWrite"] = i.text() 
           
            elif output_object.Output_Type_comboBox.currentText() =='WritePolydefix':
                
                self.set_cl.output_types.append('Polydefix')
                
                for i in output_object.WritePolydefix_optional_list:
                    
                    if i.objectName() == 'Output_NumAziWrite':
                        
                        if i.text() == '':
                            self.set_cl.output_settings["Output_NumAziWrite"] = None
                        else:
                            self.set_cl.output_settings["Output_NumAziWrite"] = i.text()
                            
                    if i.objectName() == 'Output_ElasticProperties':
                       
                        if i.text() == '': 
                            self.set_cl.output_settings["Output_ElasticProperties"] = None
                        else:
                            self.set_cl.output_settings["Output_ElasticProperties"] = i.text() 
                    
                    if i.objectName() == 'tc':   
                        
                        if i.text() == '':
                            self.set_cl.output_settings["tc"] = None
                        else:
                            self.set_cl.output_settings["tc"] = i.text() 
                    
                    if i.objectName() == 'phase':     
                        
                        if i.text() == '':
                            self.set_cl.output_settings["phase"] = None
                        else:
                            self.set_cl.output_settings["phase"] = i.text() 
                            
                    if i.objectName() == 'datafile_StartNum':  
                        
                        if i.text() == '':
                             self.set_cl.output_settings["datafile_StartNum"] = None
                        else:
                             self.set_cl.output_settings["datafile_StartNum"] = i.text() 
                    
                    if i.objectName() == 'datafile_EndNum':     
                        
                        if i.text() == '':
                             self.set_cl.output_settings["datafile_EndNum"] = None
                        else:
                             self.set_cl.output_settings["datafile_EndNum"] = i.text() 
                    
                    if i.objectName() == 'datafile_NumDigit':          
                        
                        if i.text() == '':
                             self.set_cl.output_settings["datafile_NumDigit"] = None
                        else:
                             self.set_cl.output_settings["datafile_NumDigit"] = i.text() 
                
            elif output_object.Output_Type_comboBox.currentText() =='WritePolydefixED':
                
                self.set_cl.output_types.append('PolydefixED')
                
                for i in output_object.WritePolydefixED_optional_list:
                    
                    if i.objectName() == 'tc':
                       
                        if i.text() == '':
                            self.set_cl.output_settings["tc"] = None
                        else:
                            self.set_cl.output_settings["tc"] = i.text() 
                    
                    if i.objectName() == 'phase':    
                        
                        if i.text() == '':
                            self.set_cl.output_settings["phase"] = None
                        else:
                            self.set_cl.output_settings["phase"] = i.text()
                    
                    if i.objectName() == 'Output_TemperaturePower':     
                        
                        if i.text() == '':
                             self.set_cl.output_settings["Output_TemperaturePower"] = None
                        else:
                             self.set_cl.output_settings["Output_TemperaturePower"] = i.text() 
                    
                    if i.objectName() == 'Output_tc':           
                        
                        if i.text() == '':
                             self.set_cl.output_settings["Output_tc"] = None
                        else:
                             self.set_cl.output_settings["Output_tc"] = i.text() 
                    
                    if i.objectName() == 'ElasticProperties':    
                        
                        if i.text() == '':
                             self.set_cl.output_settings["ElasticProperties"] = None
                        else:
                             self.set_cl.output_settings["ElasticProperties"] = i.text() 
            
            elif output_object.Output_Type_comboBox.currentText() =='Select Type':
                mess = QMessageBox()
                mess.setWindowIcon(QtGui.QIcon('error.JFIF'))
                mess.setIcon(QMessageBox.Icon.Warning)
                mess.setText("Please Select Output Type")
                mess.setStandardButtons(QMessageBox.StandardButton.Ok)
                mess.setWindowTitle("ALERT")
                returnValue = mess.exec()
                
        self.set_cl.save_settings()
        
    @pyqtSlot()
    def select_Data_Dir(self):
        dialog = QFileDialog.getExistingDirectory(self, 'Select Directory')
        if dialog:
                self.Directory.setText(dialog)
        self.Directory.setCursorPosition(0);
    
    @pyqtSlot()
    def select_Output_Dir(self):
        dialog = QFileDialog.getExistingDirectory(self, 'Select Directory')
        if dialog:
                self.Output_Dir_1.setText(dialog)
                self.Output_Dir_2.setText(dialog)

        self.Output_Dir_1.setCursorPosition(0);
        self.Output_Dir_2.setCursorPosition(0);
    
    @pyqtSlot()
    def Insert_Range(self):
        range_object = Range()
        self.range_list.append(range_object)
        #print(len(self.range_list))
        self.Range_Tab.addTab(range_object , QIcon(""),"Range")
       
    
    @pyqtSlot()
    def Remove_Range(self):
        if len(self.range_list)>=1:
            self.range_list.pop(self.Range_Tab.currentIndex())
            self.Range_Tab.removeTab(self.Range_Tab.currentIndex())
        else:
            mess = QMessageBox()
            mess.setWindowIcon(QtGui.QIcon('error.JFIF'))
            mess.setIcon(QMessageBox.Icon.Warning)
            mess.setText("No Ranges to remove")
            mess.setStandardButtons(QMessageBox.StandardButton.Ok)
            mess.setWindowTitle("ERROR")
            returnValue = mess.exec()

    @pyqtSlot()
    def Insert_Output(self):
        output_object = Output()
        self.output_list.append(output_object)
        self.Output_Tab.addTab(output_object, QIcon(""),"Output")
        
    @pyqtSlot()
    def Remove_Output(self):
        if len(self.output_list)>=1:
            self.output_list.pop(self.Output_Tab.currentIndex())
            self.Output_Tab.removeTab(self.Output_Tab.currentIndex())
        else:
            mess = QMessageBox()
            mess.setWindowIcon(QtGui.QIcon('error.JFIF'))
            mess.setIcon(QMessageBox.Icon.Warning)
            mess.setText("No Outputs to remove")
            mess.setStandardButtons(QMessageBox.StandardButton.Ok)
            mess.setWindowTitle("ERROR")
            returnValue = mess.exec()
        
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
        

