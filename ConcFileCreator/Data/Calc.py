'''
Created on 03.04.2016

@author: maximilian
'''

from copy import deepcopy

class Calc:


    def __init__(self):
        pass      
        

    def calculatePipe(self, r): #calculate which points in a grid with xyGrid*xyGrid-size are in the circular pipe with same radius
        pipe = []
        for i in xrange(r*2+1):
            pipe.append([])
            for j in xrange(r*2+1):
                x = i - r
                y = j - r
                if x * x + y * y <= r * r:  #function for circle with center at 0|0; all points inside the circle are OK -> less or equal (<=)
                    pipe[i].append(True)
                else:
                    pipe[i].append(False)
        return pipe     #return boolean matrix for gridpoints; true = inside pipe, false = not inside pipe
    
    def generatePipe(self, firstLayer, z, isSolvent):
        pipe = []
        r = len(firstLayer)
        c=0
        while c < z:
            pipe.append(deepcopy(firstLayer))
            c+=1
        for i in xrange(z):
            for j in xrange(r):
                for k in xrange(r):
                    if firstLayer[j][k] == False:
                        pipe[i][j][k] = 1000    #points outside the pipe have a concentration of 1000 to mark them
                    else:
                        if i != 0:
                            if isSolvent == False:       #if it is solvent, the concentration in every layer is equal to the concentration in the first layer
                                pipe[i][j][k] = 0
                                
        return pipe         #return calculated pipe with xy-radius and z length
        
        
                    
        
        
       