# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 16:35:15 2022

@author: g05296ar
"""

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.uic import loadUi
import sys
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QMainWindow, QApplication, QPushButton, QWidget, QAction, QTabWidget,QVBoxLayout
from PyQt5.QtGui import QIcon
from string import Template

class Output(QWidget):
    def __init__(self):
        super(Output, self).__init__()
        loadUi("Output_Widget.ui", self)
        self.Output_Type_comboBox.setMinimumHeight(40);
        self.Output_NumAziWrite.setMinimumHeight(40);
        self.Phase_Output.setMinimumHeight(40);
        self.Output_Elastic_Properties.setMinimumHeight(40);
        
        self.datafile_directory.setMinimumHeight(40);