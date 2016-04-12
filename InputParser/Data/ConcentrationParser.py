'''
Created on 22.01.2016

@author: maximilian
'''
from copy import deepcopy

class ConcentrationParser():
    
        
    subTag = '{'
    subOffTag = '}'
    conVec = []

    def __init__(self, substances, concentrations):
        self.conVec = self.parse(substances, concentrations)
        self.conVec = self.addZeroConcentrations(substances, self.conVec)
        
    def parse(self, subs, cons):        # kopiere die Substanzliste und ersetze die Namen durch die Konzentrationen, wenn bekannt
        i = 0
        ret = deepcopy(subs)
        while i < len(cons):
            position = 0
            while cons[i][position] != self.subTag:
                position +=1
            if cons[i][position+1:-1] in subs:
                ret[subs.index(cons[i][position+1:-1])] = cons[i][0:position]
            i += 1
        return ret
       
    def addZeroConcentrations(self, subs, cons):    #ersetze die unbekannten Konzentrationen (dort stehen noch die Namen in der Liste) durch 0 als Startkonzentration
        i = 0
        while i < len(subs):
            if subs[i] in cons:
                cons[cons.index(subs[i])] = "0"
            i += 1 
        return cons
    
    def getConVec(self):
        return self.conVec
       
 