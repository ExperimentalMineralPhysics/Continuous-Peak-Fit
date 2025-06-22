"""
Main file with which to open and run the app
"""

import sys

from PyQt6.QtWidgets import QApplication

# Add our sub directories to the module search path
sys.path.append("Widgets")
sys.path.append("cpf")

from cpf.gui.main_widget import MainWidget


def main():
    # Instantiate the QT application
    app = QApplication(sys.argv)
    main_widget = MainWidget()
    main_widget.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
