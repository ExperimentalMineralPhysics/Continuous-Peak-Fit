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
        
    def range_layout(self):
        self.Peak_Tab.currentChanged.connect(self.tabChanged)
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
        self.Save_Range_Btn.clicked.connect(self.Save_Ranges)

    def Save_Ranges(self):
        Range_min= self.Range_min.text()
        Range_max= self.Range_max.text()
        Intensity_min = self.Intensity_min.text()
        Intensity_max = self.Intensity_max.text()
        Range_Background_Val = self.Range_Background_Val.text()
        Background_Type = self.Background_Type.currentText()
        Peak_Pos_Selection = self.Peak_Pos_Selection.toPlainText()
        bg_fixed_lineEdit = self.bg_fixed_lineEdit.text()
        
        data = " \n\
             {                                                          \n\
               \"range\": [[$Range_min , $Range_max ]],                              \n\
               \"background\": [$Range_Background_Val],                                  \n\
               \"background-type\": \"$Background_Type\",                      \n\
               \"background-fixed\": [$bg_fixed_lineEdit],        \n\
               \"Imin\": $Intensity_min,                       \n\
               \"imax\": $Intensity_max ,                            \n\
               \"PeakPositionSelection\": [$Peak_Pos_Selection],      \
    "
        temp_obj = Template(data)
        after_replacing = temp_obj.substitute(Range_min=Range_min, Range_max=Range_max,
                                              Intensity_min=Intensity_min,Intensity_max=Intensity_max,
                                              Range_Background_Val=Range_Background_Val,
                                              Background_Type=Background_Type,
                                              bg_fixed_lineEdit=bg_fixed_lineEdit,
                                              Peak_Pos_Selection=Peak_Pos_Selection)
        file1 = open('Input.py', 'a')
        file1.writelines((after_replacing))
        file1.close()
        self.bg_fixed_lineEdit.setVisible(False)
        self.bg_fixed_checkbox.stateChanged.connect(self.onActivated)
        self.Save_Range_Btn.setDisabled(True)
        
    def onActivated(self):
        if self.bg_fixed_checkbox.isChecked()==True:
            self.bg_fixed_lineEdit.setVisible(True)
            bg_fixed_lineEdit= self.bg_fixed_lineEdit.text()
            
        else:
            self.bg_fixed_lineEdit.setVisible(False)    
        
    def Insert_Peak(self):
        self.Peak_Tab.addTab(Peak() , QIcon("Location of the icon"),"Peak")
     
    def Remove_Peak(self):
        self.Peak_Tab.removeTab(self.Peak_Tab.currentIndex())
        
