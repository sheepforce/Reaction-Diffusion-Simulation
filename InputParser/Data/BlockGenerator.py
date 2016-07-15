'''
Created on 22.01.2016

@author: maximilian
'''



class BlockGenerator:
   
    commentTag = "//"
    startTags = ["***useExistingFortranInput",
                 "***zgrid",
                 "***xygridradius",
                 "***reactions",
                 "***reactionRateConstants",
                 "***diffusionCoefficients",
                 "***inflow",
                 "***flowSpeed",
                 "***integrationStepwidth",
                 "***integrationIntervall",
                 "***concFileRepository",
                 "***pipeLength",
                 "***pipeRadius",
                 "***integrationmethod",
                 "***linearFlowSpeed"
                 ]
    endTag = "***end"
    clean = []
    foundBlocks = []
    blocks = []
    minimumInfo = False
   
    def __init__(self, sourceFile):
        self.blocks.append(self.extractExecCommand(sourceFile))
        self.foundBlocks.append("***exec")
        self.clean = self.cleanData(sourceFile)     #entferne stoerende Zeilen aus den Rohdaten
        self.generateBlocks()
        completedata = "***reactions" in self.foundBlocks and "***reactionRateConstants" in self.foundBlocks 
        self.minimumInfo = completedata and len(self.blocks[self.foundBlocks.index("***reactions")]) == len(self.blocks[self.foundBlocks.index("***reactionRateConstants")])
    
    def cleanData(self, data):
        i = 0
        while i < len(data):
            data[i] = data[i].replace(" ","") #entferne saemtliche Leerzeichen
            data[i] = data[i].replace("\n","")#entferne Zeilenumbrueche
            if data[i][0:len(self.commentTag)] == self.commentTag:  #mache Kommentare zu Leerzeilen
                data[i] = "" 
            i = i + 1    
        while "" in data:       #entferne Leerzeilen
            data.remove("")                        
        return data 
    
    def generateBlocks(self):
        storage = []                                        #Zwischenspeicher fuer den gerade bearbeiteten Block
        i = 0
        while i < len(self.clean):
            if self.clean[i] in self.startTags:             #suche nach einem Blockanfang
                self.foundBlocks.append(self.clean[i])      #notiere, welche Blocks in welcher Reihenfolge gefunden wurden
                i += 1
                while self.clean[i] != self.endTag:         #lade alles zwischen Anfangs- und Endtag in den storage
                    storage.append(self.clean[i])
                    i += 1
                self.blocks.append(storage)                 #und wenn der Endtag erreicht ist den storage in die Blockliste
                storage = []                                # dann setze den storage zurueck
            i += 1
    
    def getConcentrations(self):
        if "***concFileRepository" in self.foundBlocks:
            return self.blocks[self.foundBlocks.index("***concFileRepository")]  
        else:
            return [""]

    def getBlockByName(self, name):
        return self.blocks[self.foundBlocks.index(name)]
    
    def extractExecCommand(self, data):
        i = data.index("***exec\n") + 1
        data[i-1] = ""
        command = ""
        while data[i] != self.endTag + "\n":
            data[i] = data[i].replace("\n","")#entferne Zeilenumbrueche
            command = command + data[i]
            data[i] = ""
            i += 1
        data[i] = ""
        return command
                    
                
    
    
 
    
    
       