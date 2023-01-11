from PyQt6.QtGui import QIcon
from PyQt6.uic import loadUi
from PyQt6 import QtWidgets
from PyQt6 import (
    QtCore, 
    QtGui, 
    QtWidgets,
    )
from PyQt6.QtWidgets import (
    QMainWindow, 
    QApplication, 
    QPushButton, 
    QWidget, 
    QTabWidget,
    QVBoxLayout,
    )

import sys
from string import Template

from Peak_Widget import Peak

class Range(QWidget):
   
    def __init__(self):
        super(Range, self).__init__()
        loadUi("Range_Widget.ui", self)
        self.range_layout()
        self.bg_fixed_lineEdit.setEnabled(False)
        self.bg_fixed_checkbox.stateChanged.connect(lambda:bg_fixed_checkbox())
        
        self.peak_list = []
        
        def bg_fixed_checkbox():
            if self.bg_fixed_checkbox.isChecked()==True:
               self.bg_fixed_lineEdit.setEnabled(True)
            else:
               self.bg_fixed_lineEdit.setEnabled(False)
    
    def range_layout(self):
        self.Add_Peak_Btn.clicked.connect(self.Insert_Peak)
        self.Remove_Peak_Btn.clicked.connect(self.Remove_Peak)
        self.clickcount = 0 
        self.Range_min.setMinimumHeight(30);
        self.Range_max.setMinimumHeight(30);    
        self.Intensity_min.setMinimumHeight(30);  
        self.Intensity_max.setMinimumHeight(30);
        self.Range_Background_Val.setMinimumHeight(30);  
        self.Background_Type.setMinimumHeight(30);
        self.bg_fixed_lineEdit.setMinimumHeight(30);
    
    def Insert_Peak(self):
        peak_object = Peak()
        self.Peak_Tab.addTab(peak_object , QIcon("Location of the icon"),"Peak")
        self.peak_list.append(peak_object)
     
    def Remove_Peak(self):
        self.peak_list.pop(self.Peak_Tab.currentIndex())
        self.Peak_Tab.removeTab(self.Peak_Tab.currentIndex())
