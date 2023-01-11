from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.uic import loadUi
from PyQt6 import QtWidgets
from PyQt6.QtWidgets import (
    QMainWindow, 
    QApplication, 
    QPushButton, 
    QWidget, 
    QTabWidget,
    QVBoxLayout
)
from PyQt6.QtGui import (QIcon, QAction)

import sys
from string import Template

class Peak(QWidget):
    
    def __init__(self):
        super(Peak, self).__init__()
        loadUi("Widgets/Peak_Widget.ui", self)
        self.peak_layout()
        self.profile_fixed.setEnabled(False)
        self.width_fixed.setEnabled(False)
        self.height_fixed.setEnabled(False)
        self.dspace_fixed.setEnabled(False)
        
        self.profile_checkBox.stateChanged.connect(lambda:profilestatechange())
        self.width_checkBox.stateChanged.connect(lambda:widthstatechange())
        self.height_checkBox.stateChanged.connect(lambda:heghtstatechange())
        self.dspace_checkBox.stateChanged.connect(lambda:dspacestatechange())
        
        def profilestatechange():
            if self.profile_checkBox.isChecked()==True:
               self.profile_fixed.setEnabled(True)
            else:
               self.profile_fixed.setEnabled(False)
          
        def widthstatechange():
            if self.width_checkBox.isChecked()==True:
               self.width_fixed.setEnabled(True)
            else:
               self.width_fixed.setEnabled(False)
          
        def heghtstatechange():
            if self.height_checkBox.isChecked()==True:
               self.height_fixed.setEnabled(True)
            else:
               self.height_fixed.setEnabled(False)
          
        def dspacestatechange():
            if self.dspace_checkBox.isChecked()==True:
               self.dspace_fixed.setEnabled(True)
            else:
               self.dspace_fixed.setEnabled(False)
          
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
       