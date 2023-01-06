from PyQt5.QtGui import QIcon
from PyQt5.uic import loadUi
from PyQt5 import QtWidgets
from PyQt5 import (
    QtCore, 
    QtGui, 
    QtWidgets,
    )
from PyQt5.QtWidgets import (
    QMainWindow, 
    QApplication, 
    QPushButton, 
    QWidget, QAction, 
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
        self.Range_min.setMinimumHeight(40);
        self.Range_max.setMinimumHeight(40);    
        self.Intensity_min.setMinimumHeight(40);  
        self.Intensity_max.setMinimumHeight(40);
        self.Range_Background_Val.setMinimumHeight(40);  
        self.Background_Type.setMinimumHeight(40);
        self.bg_fixed_lineEdit.setMinimumHeight(40);
    
    def Insert_Peak(self):
        peak_object = Peak()
        self.Peak_Tab.addTab(peak_object , QIcon("Location of the icon"),"Peak")
        self.peak_list.append(peak_object)
     
    def Remove_Peak(self):
        self.peak_list.pop(self.Peak_Tab.currentIndex())
        self.Peak_Tab.removeTab(self.Peak_Tab.currentIndex())
