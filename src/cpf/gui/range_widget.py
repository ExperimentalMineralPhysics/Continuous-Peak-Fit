from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtGui import QDoubleValidator, QIcon, QIntValidator
from PyQt6.QtWidgets import QMessageBox, QWidget
from PyQt6.uic.load_ui import loadUi

import cpf
from cpf.gui.peak_widget import Peak


class Range(QWidget):
    def __init__(self):
        super(Range, self).__init__()
        loadUi("Widgets/cpf.gui.range_widget.ui", self)
        self.range_layout()
        self.bg_fixed_lineEdit.setEnabled(False)
        self.bg_fixed_checkbox.stateChanged.connect(lambda: bg_fixed_checkbox())
        self.set_cl = (
            cpf.settings.Settings()
        )  # Import the actual Settings object instead of creating a copy
        self.peak_list = []

        def bg_fixed_checkbox():
            if self.bg_fixed_checkbox.isChecked() == True:
                self.bg_fixed_lineEdit.setEnabled(True)
            else:
                self.bg_fixed_lineEdit.setEnabled(False)

    def range_layout(self):
        self.Add_Peak_Btn.clicked.connect(self.Insert_Peak)
        self.Remove_Peak_Btn.clicked.connect(self.Remove_Peak)
        self.clickcount = 0
        self.Range_min.setMinimumHeight(30)
        self.Range_max.setMinimumHeight(30)
        self.Intensity_min.setMinimumHeight(30)
        self.Intensity_max.setMinimumHeight(30)
        self.Range_Background_Val.setMinimumHeight(30)
        self.Background_Type.setMinimumHeight(30)
        self.bg_fixed_lineEdit.setMinimumHeight(30)
        self.Range_min.setValidator(QDoubleValidator())
        self.Range_max.setValidator(QDoubleValidator())

    def Insert_Peak(self):
        peak_object = Peak()
        self.Peak_Tab.addTab(peak_object, QIcon(""), "Peak")
        self.peak_list.append(peak_object)

    def Remove_Peak(self):
        if len(self.peak_list) >= 1:
            self.peak_list.pop(self.Peak_Tab.currentIndex())
            self.Peak_Tab.removeTab(self.Peak_Tab.currentIndex())

        else:
            mess = QMessageBox()
            mess.setIcon(QMessageBox.Icon.Warning)
            mess.setText("No Peaks to remove")
            mess.setStandardButtons(QMessageBox.StandardButton.Ok)
            mess.setWindowTitle("WARNING")
            returnValue = mess.exec()
