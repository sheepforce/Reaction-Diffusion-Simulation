'''
Created on 22.01.2016

@author: maximilian
'''
from Data.ReactionParser import ReactionParser
from Data.BlockGenerator import BlockGenerator
from Data.ConcentrationParser import ConcentrationParser
import argparse

if __name__ == "__main__":
    fortranInput = ""
    options = argparse.ArgumentParser()
    options.add_argument("-f","--file", help="enter inputfile (path)")
    args = options.parse_args()
    
    print ("Starting program")
    source = args.file
    if source == None:
        source = raw_input("Enter path of inputfile:")
    print ("Reading inputfile...")
    
    try:                                                                    #inpufile einlesen
        sourceFile = open(source).readlines()
    except Exception:
        print ("Error while reading inputfile. File not existing?")
        raw_input("Press ENTER to exit program.")
        exit()
    print ("Done. Starting to analyse data...")
    
    try:                                                                    #Bloecke finden und pruefen
        makeBlocks = BlockGenerator(sourceFile)
    except Exception:
        print ("Error while analyzing data. Is an input block without end tag?")
        raw_input("Press ENTER to exit program.")
        exit()
    
    if "***pathToExistingFortranInput" in makeBlocks.foundBlocks:
        fortranInput = makeBlocks.getBlockByName("***pathToExistingFortranInput")[0]
        print("Done. Using probably existing fortran input.")
    else:
        print ("Done. Checking minimum information...")
        if not makeBlocks.minimumInfo:
            print ("Not enough data given. Make sure to have one rate constant for each reaction and at least the blocks ***reactions, ***reaction rate constants and ***starting concentrations")
            raw_input("Press ENTER to exit program.")
            exit()
        print ("Done. Start parsing reactions.")
       
        try:                                                                    #Reaktionsgleichungen parsen
            parseReactions = ReactionParser(makeBlocks.getBlockByName("***reactions"))
        except Exception:
            print ("Error while parsing data. Check your reaction-equations in the input file.")
            raw_input("Press ENTER to exit program.")
            exit()
        print("Done. Start parsing starting concentrations.")       
    
        concentrations = []                                                  #Anfangskonzentrationen parsen
        for i in xrange(len(parseReactions.getSubVec())):
            concentrations.append(ConcentrationParser.parse(parseReactions.getSubVec()[i], makeBlocks.getConcentrations()[0], int(makeBlocks.getBlockByName("***zgrid")[0]), int(makeBlocks.getBlockByName("***zgrid")[0])))
        print("Done. Start generating input for simulation program.")  
    
        #Generierung des inputfile fuer fortran
        #pipeLength
        out = str('%e' % float(makeBlocks.getBlockByName("***pipeLength")[0])).replace("e", "d").replace("+","")
        if out[out.index("d")+1] == "0" and out[out.index("d")+1] != "-":
            print(out[:out.index("d")+1] + out[out.index("d")+2:] + " ! pipeLength")
        elif out[out.find("d-")+2] == "0":
            print(out[:out.index("d")+2] + out[out.index("d")+3:] + " ! pipeLength")
        else:
            print(out + " ! pipeLength")
    
        #pipeRadius
        out = str('%e' % float(makeBlocks.getBlockByName("***pipeRadius")[0])).replace("e", "d").replace("+","")
        if out[out.index("d")+1] == "0" and out[out.index("d")+1] != "-":
            print(out[:out.index("d")+1] + out[out.index("d")+2:] + " ! pipeRadius")
        elif out[out.find("d-")+2] == "0":
            print(out[:out.index("d")+2] + out[out.index("d")+3:] + " ! pipeRadius")
        else:
            print(out + " ! pipeRadius")
            
        #xygridradius, int
        
        #zgrid, int
        
        #flowSpeed
        
        #linearFlowSpeed
        
        #integrationStepwidth
        
        #
            
    
    
    #Ausgabe zu Testzwecken. Man sieht ja sonst nichts...
    #----------------------------------------------------
    print("")
    print("Testausgabe")
    print("-----------")
    print("")
    print("Eduktmatrix:")
    print (parseReactions.getSubVec())
    i = 0 
    while i < len(parseReactions.getEduMat()):
        print (parseReactions.getEduMat()[i])
        i += 1
    print ("")
    print("Produktmatrix:")
    print (parseReactions.getSubVec())
    i = 0 
    while i < len(parseReactions.getProMat()):
        print (parseReactions.getProMat()[i])
        i += 1
    print ("")
    print("Anfangskonzentrationen:")
    print (parseReactions.getSubVec())
    print (concentrations)
    print("")
    print("Anfangsgeschwindigkeiten:")
    i = 0 
    while i < len(makeBlocks.getBlockByName("***reactionRateConstants")):
        print (str(i+1)+ " "+ makeBlocks.getBlockByName("***reactionRateConstants")[i])
        i += 1    
        
        
        