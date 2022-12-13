import sys
from PyQt5.QtGui import QIcon
from PyQt5.uic import loadUi
from PyQt5 import QtWidgets
from PyQt5 import (
    QtCore, 
    QtGui, 
    QtWidgets
)
from PyQt5.QtWidgets import (
    QMainWindow, 
    QApplication, 
    QPushButton, 
    QWidget, 
    QAction, 
    QTabWidget,
    QVBoxLayout
)


import cpf
from cpf.settings import settings


class Output(QWidget):
    def __init__(self):
        super(Output, self).__init__()
        loadUi("Output_Widget.ui", self)
        self.Output_Type_comboBox.setMinimumHeight(40);
        self.Output_Type_comboBox.currentTextChanged.connect(self.on_combobox_changed)

    def on_combobox_changed(self):
        n = 0
        if self.Output_Type_comboBox.currentText() =='WritePolydefix':
                self.req_item1 = len(cpf.output_formatters.WritePolydefix.Requirements()[1])
                childcount = self.gridLayout_3.count()
                if childcount >=1:
                    for i in range (0,childcount):
                        item = self.gridLayout_3.itemAt(i)
                        widget = item.widget()
                        widget.deleteLater()
                for wid in range (0, self.req_item1):
                    self.lineEdit = QtWidgets.QLineEdit(self.groupBox)
                    self.lineEdit.setObjectName("lineEdit"+str(n))
                    self.lineEdit.setMinimumHeight(40);
                    self.gridLayout_3.addWidget(self.lineEdit, n, 2, 1, 1)
                    self.label_2 = QtWidgets.QLabel(self.groupBox)
                    self.label_2.setObjectName("label_2")
                    self.gridLayout_3.addWidget(self.label_2, n, 0, 1, 1)
                    self.label_2.setText(cpf.output_formatters.WritePolydefix.Requirements()[1][wid])
                    n+=1

                self.req_item11 = len(cpf.output_formatters.WritePolydefix.Requirements()[0])
                childcount2 = self.gridLayout_2.count()
                if childcount2 >=1:
                    for i in range (0,childcount2):
                        item = self.gridLayout_2.itemAt(i)
                        widget = item.widget()
                        widget.deleteLater()
                for wid in range (0, self.req_item11):
                    self.lineEdit = QtWidgets.QLineEdit(self.groupBox_2)
                    self.lineEdit.setObjectName("lineEdit"+str(n))
                    self.lineEdit.setMinimumHeight(40);
                    self.gridLayout_2.addWidget(self.lineEdit, n, 2, 1, 1)
                    self.label_2 = QtWidgets.QLabel(self.groupBox)
                    self.label_2.setObjectName("label_2")
                    self.gridLayout_3.addWidget(self.label_2, n, 0, 1, 1)
                    self.label_2.setText(cpf.output_formatters.WritePolydefix.Requirements()[1][wid])
                    n+=1          
        elif self.Output_Type_comboBox.currentText() =='WriteCoefficientTable':
                self.req_item2 = len(cpf.output_formatters.WriteCoefficientTable.Requirements()[1])
                childcount = self.gridLayout_3.count()
                if childcount >=1:
                    for i in range (0,childcount):
                        item = self.gridLayout_3.itemAt(i)
                        widget = item.widget()
                        widget.deleteLater()
                for wid in range (0, self.req_item2):
                    self.lineEdit = QtWidgets.QLineEdit(self.groupBox)
                    self.lineEdit.setObjectName("lineEdit"+str(n))
                    self.lineEdit.setMinimumHeight(40);
                    self.gridLayout_3.addWidget(self.lineEdit, n, 2, 1, 1)
                    self.label_2 = QtWidgets.QLabel(self.groupBox)
                    self.label_2.setObjectName("label_2")
                    self.label_2.setText(cpf.output_formatters.WriteCoefficientTable.Requirements()[1][wid])
                    self.gridLayout_3.addWidget(self.label_2, n, 0, 1, 1)
                    n+=1
                    
                self.req_item12 = len(cpf.output_formatters.WriteCoefficientTable.Requirements()[0])
                childcount2 = self.gridLayout_2.count()
                if childcount2 >=1:
                    for i in range (0,childcount2):
                        item = self.gridLayout_2.itemAt(i)
                        widget = item.widget()
                        widget.deleteLater() 
                for wid in range (0, self.req_item12):
                    self.lineEdit = QtWidgets.QLineEdit(self.groupBox_2)
                    self.lineEdit.setObjectName("lineEdit"+str(n))
                    self.lineEdit.setMinimumHeight(40);
                    self.gridLayout_2.addWidget(self.lineEdit, n, 2, 1, 1)
                    self.label_2 = QtWidgets.QLabel(self.groupBox)
                    self.label_2.setObjectName("label_2")
                    self.label_2.setText(cpf.output_formatters.WriteCoefficientTable.Requirements()[1][wid])
                    n+=1
        elif self.Output_Type_comboBox.currentText() =='WriteDifferentialStrain':
                self.req_item3 = len(cpf.output_formatters.WriteDifferentialStrain.Requirements()[1])
                childcount = self.gridLayout_3.count()
                if childcount >=1:
                    for i in range (0,childcount):
                        item = self.gridLayout_3.itemAt(i)
                        widget = item.widget()
                        widget.deleteLater()
                for wid in range (0, self.req_item3):
                    self.lineEdit = QtWidgets.QLineEdit(self.groupBox)
                    self.lineEdit.setObjectName("lineEdit"+str(n))
                    self.lineEdit.setMinimumHeight(40);
                    self.gridLayout_3.addWidget(self.lineEdit, n, 2, 1, 1)
                    self.label_2 = QtWidgets.QLabel(self.groupBox)
                    self.label_2.setObjectName("label_2")
                    self.label_2.setText(cpf.output_formatters.WriteDifferentialStrain.Requirements()[1][wid])
                    self.gridLayout_3.addWidget(self.label_2, n, 0, 1, 1)
                    n+=1
                    
                self.req_item13 = len(cpf.output_formatters.WriteDifferentialStrain.Requirements()[0])
                childcount2 = self.gridLayout_2.count()
                if childcount2 >=1:
                    for i in range (0,childcount2):
                        item = self.gridLayout_2.itemAt(i)
                        widget = item.widget()
                        widget.deleteLater() 
                for wid in range (0, self.req_item13):
                    self.lineEdit = QtWidgets.QLineEdit(self.groupBox_2)
                    self.lineEdit.setObjectName("lineEdit"+str(n))
                    self.lineEdit.setMinimumHeight(40);
                    self.gridLayout_2.addWidget(self.lineEdit, n, 2, 1, 1)
                    self.label_2 = QtWidgets.QLabel(self.groupBox)
                    self.label_2.setObjectName("label_2")
                    self.label_2.setText(cpf.output_formatters.WriteDifferentialStrain.Requirements()[1][wid])
                    self.gridLayout_3.addWidget(self.label_2, n, 0, 1, 1)
                    n+=1
        elif self.Output_Type_comboBox.currentText() =='WriteMultiFit':
                self.req_item4 = len(cpf.output_formatters.WriteMultiFit.Requirements()[1])
                childcount = self.gridLayout_3.count()
                if childcount >=1:
                    for i in range (0,childcount):
                        item = self.gridLayout_3.itemAt(i)
                        widget = item.widget()
                        widget.deleteLater()
                for wid in range (0, self.req_item4):
                    self.lineEdit = QtWidgets.QLineEdit(self.groupBox)
                    self.lineEdit.setObjectName("lineEdit"+str(n))
                    self.lineEdit.setMinimumHeight(40);
                    self.gridLayout_3.addWidget(self.lineEdit, n, 2, 1, 1)
                    self.label_2 = QtWidgets.QLabel(self.groupBox)
                    self.label_2.setObjectName("label_2")
                    self.label_2.setText(cpf.output_formatters.WriteMultiFit.Requirements()[1][wid])
                    self.gridLayout_3.addWidget(self.label_2, n, 0, 1, 1)
                    n+=1
                    
                self.req_item14 = len(cpf.output_formatters.WriteMultiFit.Requirements()[0])
                childcount2 = self.gridLayout_2.count()
                if childcount2 >=1:
                    for i in range (0,childcount2):
                        item = self.gridLayout_2.itemAt(i)
                        widget = item.widget()
                        widget.deleteLater() 
                for wid in range (0, self.req_item14):
                    self.lineEdit = QtWidgets.QLineEdit(self.groupBox_2)
                    self.lineEdit.setObjectName("lineEdit"+str(n))
                    self.lineEdit.setMinimumHeight(40);
                    self.gridLayout_2.addWidget(self.lineEdit, n, 0, 1, 1)
                    n+=1
        elif self.Output_Type_comboBox.currentText() =='WritePolydefixED':
                self.req_item5 = len(cpf.output_formatters.WritePolydefixED.Requirements()[1])
                childcount = self.gridLayout_3.count()
                if childcount >=1:
                    for i in range (0,childcount):
                        item = self.gridLayout_3.itemAt(i)
                        widget = item.widget()
                        widget.deleteLater()
                for wid in range (0, self.req_item5):
                    self.lineEdit = QtWidgets.QLineEdit(self.groupBox)
                    self.lineEdit.setObjectName("lineEdit"+str(n))
                    self.lineEdit.setMinimumHeight(40);
                    self.gridLayout_3.addWidget(self.lineEdit, n, 2, 1, 1)
                    self.label_2 = QtWidgets.QLabel(self.groupBox)
                    self.label_2.setObjectName("label_2")
                    self.label_2.setText(cpf.output_formatters.WritePolydefixED.Requirements()[1][wid])
                    self.gridLayout_3.addWidget(self.label_2, n, 0, 1, 1)
                    n+=1        
                    
                self.req_item15 = len(cpf.output_formatters.WritePolydefixED.Requirements()[0])
                childcount2 = self.gridLayout_2.count()
                if childcount2 >=1:
                    for i in range (0,childcount2):
                        item = self.gridLayout_2.itemAt(i)
                        widget = item.widget()
                        widget.deleteLater() 
                for wid in range (0, self.req_item15):
                    self.lineEdit = QtWidgets.QLineEdit(self.groupBox_2)
                    self.lineEdit.setObjectName("lineEdit"+str(n))
                    self.lineEdit.setMinimumHeight(40);
                    self.gridLayout_2.addWidget(self.lineEdit, n, 2, 1, 1)
                    self.label_2 = QtWidgets.QLabel(self.groupBox)
                    self.label_2.setObjectName("label_2")
                    self.label_2.setText(cpf.output_formatters.WritePolydefixED.Requirements()[1][wid])
                    self.gridLayout_3.addWidget(self.label_2, n, 0, 1, 1)
                    n+=1      
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
