from PyQt5.uic import loadUi
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QWidget

import cpf
from cpf.settings import settings


class Output(QWidget):
    def __init__(self):
        super(Output, self).__init__()
        loadUi("Output_Widget.ui", self)
        self.Output_Type_comboBox.setMinimumHeight(30);
        self.Output_Type_comboBox.currentTextChanged.connect(self.on_combobox_changed)
        
        self.WritePolydefix_optional_list = []
        self.WriteCoefficientTable_optional_list = []
        self.WriteDifferentialStrain_optional_list = []
        self.WriteMultiFit_optional_list = []
        self.WritePolydefixED_optional_list = []
        
        self.WritePolydefix_required_list = []
        self.WriteCoefficientTable_required_list = []
        self.WriteDifferentialStrain_required_list = []
        self.WriteMultiFit_required_list = []
        self.WritePolydefixED_required_list = []

    def on_combobox_changed(self):
        n = 0
        if self.Output_Type_comboBox.currentText() =='WritePolydefix':
            # Optional
                self.req_item1 = len(cpf.output_formatters.WritePolydefix.Requirements()[1])
                childcount = self.gridLayout_3.count()
                
                if childcount >=1:
                   
                    for i in range (0,childcount):
                        item = self.gridLayout_3.itemAt(i)
                        widget = item.widget()
                        widget.deleteLater()
                
                for wid in range(self.req_item1):
                  
                    self.lineEdit = QtWidgets.QLineEdit(self.groupBox)
                    self.lineEdit.setObjectName(cpf.output_formatters.WritePolydefix.Requirements()[1][wid])
                    self.lineEdit.setMinimumHeight(30);
                    self.gridLayout_3.addWidget(self.lineEdit, n, 2, 1, 1)
                    self.label_2 = QtWidgets.QLabel(self.groupBox)
                    self.label_2.setObjectName("label_2")
                    self.gridLayout_3.addWidget(self.label_2, n, 0, 1, 1)
                    self.label_2.setText(cpf.output_formatters.WritePolydefix.Requirements()[1][wid])
                    n+=1
                    self.WritePolydefix_optional_list.append(self.lineEdit)
         # Required
                self.req_item11 = len(cpf.output_formatters.WritePolydefix.Requirements()[0])
                
                childcount2 = self.gridLayout_2.count()
                
                if childcount2 >=1:
                   
                    for i in range (0,childcount2):
                        item = self.gridLayout_2.itemAt(i)
                        widget = item.widget()
                        widget.deleteLater()
                
                for wid in range (0, self.req_item11):
                    
                    self.lineEdit = QtWidgets.QLineEdit(self.groupBox_2)
                    self.lineEdit.setObjectName(cpf.output_formatters.WritePolydefix.Requirements()[0][wid])
                    self.lineEdit.setMinimumHeight(30);
                    self.gridLayout_2.addWidget(self.lineEdit, n, 2, 1, 1)
                    self.label_2 = QtWidgets.QLabel(self.groupBox)
                    self.label_2.setObjectName("label_2")
                    self.gridLayout_2.addWidget(self.label_2, n, 0, 1, 1)
                    self.label_2.setText(cpf.output_formatters.WritePolydefix.Requirements()[0][wid])
                    n+=1     
                    self.WritePolydefix_required_list.append(self.lineEdit)
                    
        elif self.Output_Type_comboBox.currentText() =='WriteCoefficientTable':
        # Optional
                self.req_item2 = len(cpf.output_formatters.WriteCoefficientTable.Requirements()[1])
                childcount = self.gridLayout_3.count()
                
                if childcount >=1:
                    
                    for i in range (0,childcount):
                        item = self.gridLayout_3.itemAt(i)
                        widget = item.widget()
                        widget.deleteLater()
                
                for wid in range (0, self.req_item2):
                   
                    self.lineEdit = QtWidgets.QLineEdit(self.groupBox)
                    self.lineEdit.setObjectName(cpf.output_formatters.WriteCoefficientTable.Requirements()[1][wid])
                    self.lineEdit.setMinimumHeight(30);
                    self.gridLayout_3.addWidget(self.lineEdit, n, 2, 1, 1)
                    self.label_2 = QtWidgets.QLabel(self.groupBox)
                    self.label_2.setObjectName("label_2")
                    self.label_2.setText(cpf.output_formatters.WriteCoefficientTable.Requirements()[1][wid])
                    self.gridLayout_3.addWidget(self.label_2, n, 0, 1, 1)
                    n+=1
                    self.WriteCoefficientTable_optional_list.append(self.lineEdit)
        
        # Required            
                self.req_item12 = len(cpf.output_formatters.WriteCoefficientTable.Requirements()[0])
                childcount2 = self.gridLayout_2.count()
                
                if childcount2 >=1:
                   
                    for i in range (0,childcount2):
                        item = self.gridLayout_2.itemAt(i)
                        widget = item.widget()
                        widget.deleteLater() 
                
                for wid in range (0, self.req_item12):
                    self.lineEdit = QtWidgets.QLineEdit(self.groupBox_2)
                    self.lineEdit.setObjectName(cpf.output_formatters.WriteCoefficientTable.Requirements()[0][wid])
                    self.lineEdit.setMinimumHeight(30);
                    self.gridLayout_2.addWidget(self.lineEdit, n, 2, 1, 1)
                    self.label_2 = QtWidgets.QLabel(self.groupBox)
                    self.label_2.setObjectName("label_2")
                    self.gridLayout_2.addWidget(self.label_2, n, 0, 1, 1)
                    self.label_2.setText(cpf.output_formatters.WriteCoefficientTable.Requirements()[0][wid])
                    n+=1
                    self.WriteCoefficientTable_required_list.append(self.lineEdit)
                    
        elif self.Output_Type_comboBox.currentText() =='WriteDifferentialStrain':
        # Optional            
                self.req_item3 = len(cpf.output_formatters.WriteDifferentialStrain.Requirements()[1])
                childcount = self.gridLayout_3.count()
                
                if childcount >=1:
                    
                    for i in range (0,childcount):
                        item = self.gridLayout_3.itemAt(i)
                        widget = item.widget()
                        widget.deleteLater()
                
                for wid in range (0, self.req_item3):
                    self.lineEdit = QtWidgets.QLineEdit(self.groupBox)
                    self.lineEdit.setObjectName(cpf.output_formatters.WriteDifferentialStrain.Requirements()[1][wid])
                    self.lineEdit.setMinimumHeight(30);
                    self.gridLayout_3.addWidget(self.lineEdit, n, 2, 1, 1)
                    self.label_2 = QtWidgets.QLabel(self.groupBox)
                    self.label_2.setObjectName("label_2")
                    self.label_2.setText(cpf.output_formatters.WriteDifferentialStrain.Requirements()[1][wid])
                    self.gridLayout_3.addWidget(self.label_2, n, 0, 1, 1)
                    n+=1
                    self.WriteDifferentialStrain_optional_list.append(self.lineEdit)
        #Required            
                self.req_item13 = len(cpf.output_formatters.WriteDifferentialStrain.Requirements()[0])
                childcount2 = self.gridLayout_2.count()
                
                if childcount2 >=1:
                    
                    for i in range (0,childcount2):
                        item = self.gridLayout_2.itemAt(i)
                        widget = item.widget()
                        widget.deleteLater() 
                
                for wid in range (0, self.req_item13):
                    self.lineEdit = QtWidgets.QLineEdit(self.groupBox_2)
                    self.lineEdit.setObjectName(cpf.output_formatters.WriteDifferentialStrain.Requirements()[0][wid])
                    self.lineEdit.setMinimumHeight(30);
                    self.gridLayout_2.addWidget(self.lineEdit, n, 2, 1, 1)
                    self.label_2 = QtWidgets.QLabel(self.groupBox)
                    self.label_2.setObjectName("label_2")
                    self.label_2.setText(cpf.output_formatters.WriteDifferentialStrain.Requirements()[0][wid])
                    self.gridLayout_2.addWidget(self.label_2, n, 0, 1, 1)
                    n+=1
                    self.WriteDifferentialStrain_required_list.append(self.linEdit)
                    
        elif self.Output_Type_comboBox.currentText() =='WriteMultiFit':
        # Optional            
                self.req_item4 = len(cpf.output_formatters.WriteMultiFit.Requirements()[1])
                childcount = self.gridLayout_3.count()
                
                if childcount >=1:
                    
                    for i in range (0,childcount):
                        item = self.gridLayout_3.itemAt(i)
                        widget = item.widget()
                        widget.deleteLater()
                
                for wid in range (0, self.req_item4):
                    self.lineEdit = QtWidgets.QLineEdit(self.groupBox)
                    self.lineEdit.setObjectName(cpf.output_formatters.WriteMultiFit.Requirements()[1][wid])
                    self.lineEdit.setMinimumHeight(30);
                    self.gridLayout_3.addWidget(self.lineEdit, n, 2, 1, 1)
                    self.label_2 = QtWidgets.QLabel(self.groupBox)
                    self.label_2.setObjectName("label_2")
                    self.label_2.setText(cpf.output_formatters.WriteMultiFit.Requirements()[1][wid])
                    self.gridLayout_3.addWidget(self.label_2, n, 0, 1, 1)
                    n+=1
                    self.WriteMultiFit_optional_list.append(self.lineEdit)
         #Required           
                self.req_item14 = len(cpf.output_formatters.WriteMultiFit.Requirements()[0])
                childcount2 = self.gridLayout_2.count()
                
                if childcount2 >=1:
                    
                    for i in range (0,childcount2):
                        item = self.gridLayout_2.itemAt(i)
                        widget = item.widget()
                        widget.deleteLater() 
                
                for wid in range (0, self.req_item14):
                    self.lineEdit = QtWidgets.QLineEdit(self.groupBox_2)
                    self.lineEdit.setObjectName(cpf.output_formatters.WriteMultiFit.Requirements()[0])
                    self.lineEdit.setMinimumHeight(30);
                    self.gridLayout_2.addWidget(self.lineEdit, n, 0, 1, 1)
                    self.label_2 = QtWidgets.QLabel(self.groupBox)
                    self.label_2.setObjectName("label_2")
                    self.label_2.setText(cpf.output_formatters.WriteMultiFit.Requirements()[0][wid])
                    self.gridLayout_2.addWidget(self.label_2, n, 0, 1, 1)
                    n+=1
                    self.WriteMultiFit_required_list.append(self.linEdit)
                    
        elif self.Output_Type_comboBox.currentText() =='WritePolydefixED':
        # Optional
                self.req_item5 = len(cpf.output_formatters.WritePolydefixED.Requirements()[1])
                childcount = self.gridLayout_3.count()
                
                if childcount >=1:
                    
                    for i in range (0,childcount):
                        item = self.gridLayout_3.itemAt(i)
                        widget = item.widget()
                        widget.deleteLater()
                
                for wid in range (0, self.req_item5):
                    self.lineEdit = QtWidgets.QLineEdit(self.groupBox)
                    self.lineEdit.setObjectName(cpf.output_formatters.WritePolydefixED.Requirements()[1][wid])
                    self.lineEdit.setMinimumHeight(30);
                    self.gridLayout_3.addWidget(self.lineEdit, n, 2, 1, 1)
                    self.label_2 = QtWidgets.QLabel(self.groupBox)
                    self.label_2.setObjectName("label_2")
                    self.label_2.setText(cpf.output_formatters.WritePolydefixED.Requirements()[1][wid])
                    self.gridLayout_3.addWidget(self.label_2, n, 0, 1, 1)
                    n+=1       
                    self.WritePolydefixED_optional_list.append(self.lineEdit)
         #Required           
                self.req_item15 = len(cpf.output_formatters.WritePolydefixED.Requirements()[0])
                childcount2 = self.gridLayout_2.count()
                
                if childcount2 >=1:
                   
                    for i in range (0,childcount2):
                        item = self.gridLayout_2.itemAt(i)
                        widget = item.widget()
                        widget.deleteLater() 
                
                for wid in range (0, self.req_item15):
                    self.lineEdit = QtWidgets.QLineEdit(self.groupBox_2)
                    self.lineEdit.setObjectName(cpf.output_formatters.WritePolydefixED.Requirements()[0][wid])
                    self.lineEdit.setMinimumHeight(30);
                    self.gridLayout_2.addWidget(self.lineEdit, n, 2, 1, 1)
                    self.label_2 = QtWidgets.QLabel(self.groupBox)
                    self.label_2.setObjectName("label_2")
                    self.label_2.setText(cpf.output_formatters.WritePolydefixED.Requirements()[0][wid])
                    self.gridLayout_2.addWidget(self.label_2, n, 0, 1, 1)
                    n+=1      
                    self.WritePolydefixED_required_list.append(self.linEdit)
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
