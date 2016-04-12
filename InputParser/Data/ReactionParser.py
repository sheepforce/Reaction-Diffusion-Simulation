'''
Created on 22.01.2016

@author: maximilian
'''

class ReactionParser():
    
    subTag = '{'
    subOffTag = '}'
    reactionArrow = "->"
    subVec = []
    eduMat = []
    proMat = []
    
    
    def __init__(self, matrix):
        i = 0
        while i < len(matrix):
            self.parseLine(matrix[i])
            i = i + 1
        print("Done. Adjusting Matrices.")
        self.adjustLengthOfReactandMatrices()
        
    def parseLine(self, equation):
        position = 0
        coeVec = self.generateCoeVec()                                  #generiere coeVec fuer Edukte                   
        if equation[0:2] == self.reactionArrow:                         #Reaktionspfeil steht am Anfang der Gleichung (Zufluss)
                self.eduMat.append(coeVec)                              #-> haenge Koeffizienten leer an Eduktmatrix an
                coeVec = self.generateCoeVec()                          #-> setze coeVec fuer Produkte zurueck
                position = position + 2                                 #-> Reaktionspfeil ueberspringen
                
        while position < len(equation):                                 #wenn Ende des Strings erreicht ist abbrechen
            factor = "1"                                                #default-Wert fuer stoech. Koeff.
            substance = ""                                              #Substanzname
            
            #Bestimmung des stoechiometrischen Koeffizienten
            #-----------------------------------------------
            range = 0                                           #counter fuer Laenge des relevanten Bereichs
            while equation[position+range] != self.subTag:      
                range = range+1                                 #bestimme Laenge des Stringteils, in dem der Wert fuer den stoech. Koeff. steht
            if range != 0:
                factor = equation[position:position+range]      #und speichere ihn
            position = position + range                         #setze dann den Positionszeiger hinter den ausgewerteten Stringteil
            
            
            #Bestimmung des Substanznamens
            #-----------------------------            
            stringRange = 0                                             #counter fuer Laenge des relevanten Bereichs, RESET
            while equation[position+stringRange] != self.subOffTag:   
                stringRange = stringRange+1                             #bestimme Laenge des Stringteils, in dem der Wert fuer den Substanznamen. steht
            substance = equation[position+1:position+stringRange]       #und speichere ihn
            position = position + stringRange + 1                       #setze dann den Positionszeiger hinter den ausgewerteten Stringteil
            if substance not in self.subVec:                            #neue Substanz: neue Eintraege in subVec und coeVec
                self.subVec.append(substance)
                coeVec.append("") 
            coeVec[self.subVec.index(substance)] =  factor        
            if equation[position:position+2] == self.reactionArrow:     #Reaktionspfeil erreicht 
                self.eduMat.append(coeVec)                              #-> haenge Koeffizienten an Eduktmatrix an
                coeVec = self.generateCoeVec()                          #-> setze coeVec fuer Produkte zurueck
                position = position + 2                                 #-> Reaktionspfeil ueberspringen
            else:
                position = position + 1                                 #Reaktionspfeil nicht erreicht -> naechstes Zeichen ist + und wird uebersprungen
        
            
        self.proMat.append(coeVec)
         
    def generateCoeVec(self):           #generiert Koeffizientenvektor fuer eine Reaktion
        coeVec = []                     #fuer jeden bisher gefundenen Stoff gibt es einen Eintrag, dessen Index mit dem im subVec uebereinstimmt
        i = 0
        while i < len(self.subVec):
            coeVec.append("0")
            i = i + 1
        return coeVec
    
    def adjustLengthOfReactandMatrices(self):                       #sorge dafuer, dass alle Zeilen in einer Matrix die gleiche Laenge haben
        i = 0
        while i < len(self.proMat):                                 #hier koennte auch self.eduMat stehen; beide haben fuer eine Gleichung eine Zeile
            while len(self.eduMat[i]) < len(self.proMat[-1]):       #laengste Zeile ist systembedingt immer die letzte der Produkmatrix
                self.eduMat[i].append("0")                          #die kuerzeren Zeilen werden mit Nullen aufgefuellt, sowohl in der Eduktmatrix...
            while len(self.proMat[i]) < len(self.proMat[-1]):
                self.proMat[i].append("0")                          #...als auch in der Produktmatrix    
            i = i + 1                                               #und zwar in jeder Zeile der Matrizen     
            
    
    def getSubVec(self):
        return self.subVec   
    
    def getEduMat(self):
        return self.eduMat    
    
    def getProMat(self):
        return self.proMat
                  
  
      

