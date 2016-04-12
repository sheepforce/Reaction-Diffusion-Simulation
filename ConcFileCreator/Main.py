from Data.Calc import Calc
#from Data.GUI import GUI
from PyQt4.QtCore import *
from PyQt4.QtGui import *
import sys
import sip

class GUI(QWidget):
    
    app = QApplication(sys.argv)
    
    MINXSIZE = 750
    MINYSIZE = 750
    
    screenx = 0
    screeny = 0
    
    base = QWidget
    
    pipeArea = QScrollArea
    
    makePipe = QPushButton
    savePipe = QPushButton
    
    subsIn = QPlainTextEdit
    pathIn = QPlainTextEdit
    initIn = QPlainTextEdit
    rIn = QPlainTextEdit
    zIn = QPlainTextEdit
    
    concFields = []
    
    isSolventSel = QCheckBox
    
    z = 0
    
    
    
    def __init__(self):
        QWidget.__init__(self)
         
        
        #determine screen dimensions
        screen = self.app.desktop().screenGeometry()
        self.screenx, self.screeny = screen.width(), screen.height()
        
        #set a screen-sized scrollable window as base for other gui objects
        scrollbase = QScrollArea() 
        scrollbase.resize(self.screenx, self.screeny) 
        #scrollbase.resize(self.MINXSIZE, self.MINYSIZE)
        
        #add a frame to the scrollable window with minimum MINXSIZE x MINYSIZE and maximum screen size dimension
        self.base = QWidget()
        self.base.resize(self.screenx, self.screeny)
       # base.resize(scrollbase.size())
        if self.screenx < self.MINXSIZE:
            self.base.resize(self.MINXSIZE, self.base.height())
        if self.screeny < self.MINYSIZE:
            self.base.resize(self.base.width(), self.MINYSIZE)
        scrollbase.setWidget(self.base)
        scrollbase.setWindowTitle("Concentration File Creator")
        
        #setup the labels, input fields and buttons
        subsLabel = QLabel("substance name:", self.base)
        subsLabel.move(10, 10)
        subsLabel.resize(150, 30)
        
        self.subsIn = QPlainTextEdit(self.base)
        self.subsIn.move(150, 10)
        self.subsIn.setToolTip("Enter the compound's name. It will be used as the file name.")
        self.subsIn.resize(self.base.width() - 160, 30)
                
        self.pathIn = QPlainTextEdit("enter path to target directory", self.base)
        self.pathIn.move(10, 50)
        self.pathIn.setToolTip("enter path to target directory")
        self.pathIn.resize(self.base.width() - 20, 30)
        
        isSolvLabel = QLabel("is solvent:", self.base)
        isSolvLabel.move(10, 90)
        isSolvLabel.resize(150, 30)
        
        self.isSolventSel = QCheckBox(self.base)
        self.isSolventSel.setToolTip("A solvents starting concentration has everywhere in the pipe the same value as in the first layers central grid point.")
        self.isSolventSel.move(160, 95)
        
        initLabel = QLabel("initial value:", self.base)
        initLabel.move(10, 130)
        initLabel.resize(150, 30)
        
        self.initIn = QPlainTextEdit(self.base)
        self.initIn.move(160, 130)
        self.initIn.setToolTip("Enter an initial concentration. Each grid point will be filled with this concentration, but can still be edited. Use '.', not ','.")
        self.initIn.resize(70, 30)
        
        rLabel = QLabel("grid points in radius:", self.base)
        rLabel.move(10, 170)
        rLabel.resize(150, 30)
        self.rIn = QPlainTextEdit(self.base)
        self.rIn.move(160, 170)
        self.rIn.setToolTip("Enter pipe radius. There will be 2r+1 grid points in x and y direction. Has to be integer.")
        self.rIn.resize(70, 30)
        
        zLabel = QLabel("grid points in z:", self.base)
        zLabel.move(10, 210)
        zLabel.resize(150, 30)
        
        self.zIn = QPlainTextEdit(self.base)
        self.zIn.move(160, 210)
        self.zIn.setToolTip("Enter pipe length. There will be z grid points in z direction. Has to be integer.")
        self.zIn.resize(70, 30)
        
        self.makePipe = QPushButton("build pipe", self.base)
        self.makePipe.move(10, 250)
        self.makePipe.setToolTip("Generate pipe from xy and z grid values.")
        self.makePipe.resize(110, 30)
        self.makePipe.clicked.connect(self.drawPipe)
        
        self.savePipe = QPushButton("save pipe", self.base)
        self.savePipe.move(120, 250)
        self.savePipe.setToolTip("Save concentration file generated from pipe.")
        self.savePipe.resize(110, 30)
        self.savePipe.clicked.connect(self.savePipeToFile)
        
        scrollbase.show() 
        sys.exit(self.app.exec_())
    
    
        
    def drawPipe(self):

        
        self.concFields = []
        self.pipeArea = QScrollArea(self.base)
        self.pipeArea.move(260, 90)
        self.pipeArea.resize(self.base.width()-270, self.base.height()-100)
        
        try:
            r = int(self.rIn.toPlainText())
            self.z = int(self.zIn.toPlainText())
            init = float(self.initIn.toPlainText())
            
        except Exception:
            error = QPlainTextEdit("Error while generating pipe. Check the values of r, z and initial value.", self.pipeArea)
            error.move(0,0)
            error.resize(self.pipeArea.size())
            self.pipeArea.show()
            return
        
        frame = QWidget(self.pipeArea)
        frame.move(0,0)
        frame.resize(10+(2*r+1)*60, 10+(2*r+1)*30)
        self.pipeArea.setAutoFillBackground(True)
        self.pipeArea.setWidget(frame)
        calc = Calc()
        pipe = calc.calculatePipe(r)
        for i in xrange(r*2+1):
            self.concFields.append([])
            for j in xrange(r*2+1):
                if pipe[i][j]:
                    self.concFields[i].append(QPlainTextEdit)
                    self.concFields[i][j] = QPlainTextEdit(str(init), frame)
                    self.concFields[i][j].move(60 * j + 10, 10 + 30 * i)
                    self.concFields[i][j].resize(50, 30)
                else:
                    self.concFields[i].append(False)
        self.pipeArea.show()
    
    
    def savePipeToFile(self):
        calc = Calc()
        firstLayer = []
        grid = len(self.concFields)
        for i in xrange(grid):
            firstLayer.append([])
            for j in xrange(grid):
                if self.concFields[i][j] == False:
                    firstLayer[i].append(False)
                else:
                    firstLayer[i].append(str(self.concFields[i][j].toPlainText()))
        
        pipeOut = calc.generatePipe(firstLayer, int(self.zIn.toPlainText()), self.isSolventSel.isChecked())
        pathToOut = str(self.pathIn.toPlainText())+str(self.subsIn.toPlainText())+".dat"
        
        self.z = int(self.zIn.toPlainText())
        f = open(pathToOut, "w")
        f.write("Steps output:1" + "\n")
        f.write("X:" + str(len(firstLayer)) + "\n")
        f.write("Y:" + str(len(firstLayer)) +"\n")
        f.write("Z:" + str(self.z) + "\n")
        f.write("\n")
       
        for i in xrange(self.z):
            for j in xrange(grid):
                for k in xrange(grid):
                    f.write(str(pipeOut[i][j][k]) + "\n")
        f.close()
        
        
             
           
        
        
        
        
        
        
        
        
        



        
        

if __name__ == '__main__':    
    test = GUI()