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
        