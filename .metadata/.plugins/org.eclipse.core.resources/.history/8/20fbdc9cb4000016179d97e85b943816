'''
Created on 22.01.2016

@author: maximilian
'''
from Data.ReactionParser import ReactionParser
from Data.BlockGenerator import BlockGenerator
from Data.ConcentrationParser import ConcentrationParser
import argparse

if __name__ == "__main__":
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
    print ("Done. Checking minimum information...")
    if not makeBlocks.minimumInfo:
        print ("Not enough data given. Make sure to have one rate constant for each reaction and at least the blocks ***reactions, ***reaction rate constants and ***starting concentrations")
        raw_input("Press ENTER to exit program.")
        exit()
    print ("Done. Start parsing reactions.")
       
    try:                                                                    #Reaktionsgleichungen parsen
        parseReactions = ReactionParser(makeBlocks.getReactions())
    except Exception:
        print ("Error while parsing data. Check your reaction-equations in the input file.")
        raw_input("Press ENTER to exit program.")
        exit()
    print("Done. Start parsing starting concentrations.")       
    
    
    
    try:                                                                    #Anfangskonzentrationen parsen
       parseConcentrations = ConcentrationParser(parseReactions.getSubVec(), makeBlocks.getConcentrations())
    except Exception:
        print ("Error while parsing data. Check your starting concentrations in the input file.")
        raw_input("Press ENTER to exit program.")
        exit()
    print("Done. Start something else...")  
    
    
    
    
    
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
    print (parseConcentrations.getConVec())
    print("")
    print("Anfangsgeschwindigkeiten:")
    i = 0 
    while i < len(makeBlocks.getRateConstants()):
        print (str(i+1)+ " "+ makeBlocks.getRateConstants()[i])
        i += 1    
        
        
        