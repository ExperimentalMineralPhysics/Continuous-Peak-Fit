import sys
from PyQt6.QtWidgets import QApplication
import os

# Add our sub directories to the module search path
sys.path.append('Widgets')
sys.path.append('cpf')

from Main_Widget import Main_Widget

def main():
    # Instantiate the QT application
    app = QApplication(sys.argv)
    main_widget = Main_Widget()
    main_widget.show()
    sys.exit(app.exec())
  
	  
if __name__ == '__main__':
	main()