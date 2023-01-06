import sys
from PyQt5.QtWidgets import QApplication

# Add our sub directories to the module search path
sys.path.append('Widgets')
sys.path.append('cpf')

from Main_Widget import MainWindow

def main():
	# Instantiate the QT application
    app = QApplication(sys.argv)
    mainwindow = MainWindow()
    mainwindow.show() 
    try:        
      sys.exit(app.exec_())
    except:
      os._exit(00)
	  
if __name__ == '__main__':
	main()