import sys
from PyQt6.QtWidgets import QApplication
import os

# Add our sub directories to the module search path
sys.path.append('Widgets')
sys.path.append('cpf')

from Main_Widget import Main_Widget

def main():
    # Instantiate the QT application
    try: 
        app = QApplication(sys.argv)
        main_widget = Main_Widget()
        main_widget.showMaximized()
        sys.exit(app.exec())
    except Exception as e:
        print("ERROR:", e)
        os._exit(00)
  
	  
if __name__ == '__main__':
	main()