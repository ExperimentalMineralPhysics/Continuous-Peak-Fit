from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.uic import loadUi
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import (
    QMainWindow, 
    QApplication, 
    QPushButton, 
    QWidget, 
    QAction, 
    QTabWidget,
    QVBoxLayout
)
from PyQt5.QtGui import QIcon

import sys
from string import Template

class Peak(QWidget):
    
    def __init__(self):
        super(Peak, self).__init__()
        loadUi("Peak_Widget.ui", self)
        self.peak_layout()
        
    def peak_layout(self):
        self.phase_peak.setMinimumHeight(40);
        self.hkl.setMinimumHeight(40);
        self.d_space_peak.setMinimumHeight(40);
        self.height_peak.setMinimumHeight(40);
        self.profile_peak.setMinimumHeight(40);
        self.width_peak.setMinimumHeight(40);
        self.symmetry_peak.setMinimumHeight(40);
        self.dspace_fixed.setMinimumHeight(40);
        self.height_fixed.setMinimumHeight(40);
        self.profile_fixed.setMinimumHeight(40);
        self.width_fixed.setMinimumHeight(40);
        self.d_space_type.setMinimumHeight(40);
        self.height_peak_type.setMinimumHeight(40);
        self.profile_peak_type.setMinimumHeight(40);
        self.width_peak_type.setMinimumHeight(40);
        self.Save_Peak_Btn.clicked.connect(self.Save_Peaks)
        self.clickcount = 1
        
    def Save_Peaks(self):
        phase_peak = self.phase_peak.text()
        hkl = self.hkl.text()
        d_space_peak = self.d_space_peak.text()
        d_space_type = self.d_space_type.currentText()
        dspace_fixed = self.dspace_fixed.text()
        height_peak = self.height_peak.text()
        height_peak_type = self.height_peak_type.currentText()
        height_fixed =self.height_fixed.text()
        profile_peak = self.profile_peak.text()
        profile_peak_type = self.profile_peak_type.currentText()
        profile_fixed = self.profile_fixed.text()
        width_peak = self.width_peak.text()
        width_peak_type =self.width_peak_type.currentText()
        width_fixed = self.width_fixed.text()
        symmetry_peak = self.symmetry_peak.text()
        
        if (self.clickcount == 1) :
        
            data =  " \n\
               \"peak\": [{                                              \n\
                          \"phase\": \"$phase_peak\",                               \n\
                          \"hkl\": '$hkl',                                        \n\
                          \"d-space\": $d_space_peak,                                     \n\
                          \"d-space-type\": \"$d_space_type\",                 \n\
                          \"d-space_fixed\": $dspace_fixed,                      \n\
                          \"height\": $height_peak,                                       \n\
                          \"height-type\": \"$height_peak_type\",                           \n\
                          \"height_fixed\": $height_fixed,   \n\
                          \"profile\": $profile_peak,                                     \n\
                          \"profile-type\": $profile_peak_type       \n\
                          \"profile_fixed\": $profile_fixed,                             \n\
                          \"width\": $width_peak,                                       \n\
                          \"width-type\": $width_peak_type,       \n\
                          \"width_fixed\": $width_fixed ,                \n\
                          \"symmetry\": $symmetry_peak                                     \n\
                        },"
            temp_obj = Template(data)
            after_replacing = temp_obj.substitute(phase_peak=phase_peak, 
                                                  hkl=hkl,
                                                  d_space_peak= d_space_peak, 
                                                  d_space_type= d_space_type,
                                                  dspace_fixed=dspace_fixed,
                                                  height_peak=height_peak,
                                                  height_peak_type =height_peak_type,
                                                  height_fixed =height_fixed,
                                                  profile_peak= profile_peak,
                                                  profile_peak_type  = profile_peak_type,
                                                  profile_fixed = profile_fixed,
                                                  width_peak =width_peak,
                                                  width_peak_type= width_peak_type,
                                                  width_fixed  =width_fixed ,
                                                  symmetry_peak = symmetry_peak
                                                  )
            file1 = open('input.py', 'a')
            file1.writelines((after_replacing))
            file1.close()
        else:
            
            data =       "\n\ {                                              \n\
                               \"phase\": \"$phase_peak\",                               \n\
                               \"hkl\": '$hkl',                                        \n\
                               \"d-space\": $d_space_peak,                                     \n\
                               \"d-space-type\": \"$d_space_type\",                 \n\
                               \"d-space_fixed\": $dspace_fixed,                      \n\
                               \"height\": $height_peak,                                       \n\
                               \"height-type\": \"$height_peak_type\",                           \n\
                               \"height_fixed\": $height_fixed,   \n\
                               \"profile\": $profile_peak,                                     \n\
                               \"profile-type\": $profile_peak_type       \n\
                               \"profile_fixed\": $profile_fixed,                             \n\
                               \"width\": $width_peak,                                       \n\
                               \"width-type\": $width_peak_type,       \n\
                               \"width_fixed\": $width_fixed ,                \n\
                               \"symmetry\": $symmetry_peak                                     \n\
                             },"
           
    
            
            temp_obj = Template(data)
            after_replacing = temp_obj.substitute(phase_peak=phase_peak, 
                                                  hkl=hkl,
                                                  d_space_peak= d_space_peak, 
                                                  d_space_type = d_space_type,
                                                  dspace_fixed=dspace_fixed,
                                                  height_peak=height_peak,
                                                  height_peak_type =height_peak_type,
                                                  height_fixed =height_fixed,
                                                  profile_peak= profile_peak,
                                                  profile_peak_type  = profile_peak_type,
                                                  profile_fixed = profile_fixed,
                                                  width_peak =width_peak,
                                                  width_peak_type= width_peak_type,
                                                  width_fixed  =width_fixed ,
                                                  symmetry_peak = symmetry_peak
                                                  )
            file1 = open('Input.py', 'a')
            file1.writelines((after_replacing))
            file1.close()
        self.clickcount +=1
        self.Save_Peak_Btn.setDisabled(True)
