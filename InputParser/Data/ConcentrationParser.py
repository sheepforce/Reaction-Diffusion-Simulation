'''
Created on 22.01.2016

@author: maximilian
'''

class ConcentrationParser():

    def __init__(self):
        pass
        
    @staticmethod
    def parse(subs, consRepPath, r, z):
        ret = []
        i = 1
        try:
            if consRepPath == "":
                raise Exception
            import fileinput  
            for line in fileinput.input([consRepPath + subs + '.dat']):
                if i == 1:
                    print("reading data from" + consRepPath + subs + ".dat")
                    steps = int(line.replace(" ","").replace("Stepsoutput:", ""))
                if i > 5 + (2*r+1)*(2*r+1)*z *(steps-1):
                    ret.append(line.replace("\n",""))
                i+=1
        except Exception:
            print("no concentration data found for " + subs + ", using default")
            for i in xrange((2*r+1)*(2*r+1)*z):
                ret.append("0")
        print("parsing of " + subs + "s starting concentration finished")
        return ret
    
    
            
                
                        
 