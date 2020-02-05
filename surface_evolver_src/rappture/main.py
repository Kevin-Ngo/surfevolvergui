'''
    Project: GUI for the Surface Evolver - Polycrystalline Grain Growth - main file
    Author: Kevin Ngo
    Program: Purdue University - Network for Computational Nanotechnology
    Date: Summer 2018
    Collaborators: Lucas Robinson
    Faculty Advisor: John Blendell
    Acknowledgements: Jean Taylor and Ken Brakke
'''

import Rappture
import sys
import numpy as np
from Rappture.tools import executeCommand as execute
import os
from surfEvolverGUIFunctions import *  # This is a module that contains all the functions that are used to process the *.txt files created Surface Evolver containing information about the grains during the simulation

# uncomment these to redirect stdout and stderr
# to files for debugging.
# sys.stderr = open('debug.err', 'w')
# sys.stdout = open('debug.out', 'w')

# open the XML file containing the run parameters
io = Rappture.PyXml(sys.argv[1])

# spit out progress messages as you go along...
Rappture.Utils.progress(0, "Starting...")  # Starting the program

#########################################################
# Get input values from Rappture
#########################################################

# get input value for input.integer(numberOfGrains)
numberOfGrains = int(io['input.integer(numberOfGrains).current'].value)

# get input value for input.choice(typeOfSimulation)
typeOfSimulation = io['input.choice(typeOfSimulation).current'].value

# get input value for input.boolean(enableColor)
# returns value as string "yes" or "no"
enableColor = io['input.boolean(enableColor).current'].value == 'yes'

# get input value for input.integer(numberOfIterations)
numberOfIterations = int(io['input.integer(numberOfIterations).current'].value)

#########################################################
#  Add your code here for the main body of your program
#########################################################

#########################################################
# ***** EDIT PATH(s) HERE - start ***** #

# Path to main rappture file (directory when svn checkout)
rappturePath = os.path.dirname(os.path.realpath(__file__))
toolPath = os.path.dirname(rappturePath)

# Path to store the commands to append and also the images generated
pathToData = os.path.join(toolPath,'data/')

# Path to where the output files will be placed (directory where text files containing data for output are placed)
pathToOutputDirectory = os.path.join(toolPath,'rappture/')

# Path to the vor2fe.exe
pathToVor = os.path.join(toolPath,'bin/')

# Path to the graph data report (actual file for Area/Number of grains/Energy(curvature))
pathToOutputFileX = "xxx-report.txt"

# Path to the graph data report (actual file for Number of sides per grain)
pathToOutputFileSide = "sides.txt"

# Path to the *.fe file (append commands to this file)
pathToFe = os.getcwd() + "/" +  str(numberOfGrains) + "grains.fe"

# ***** EDIT PATH(s) HERE  - end *****#
#########################################################

Rappture.Utils.progress(25, "Loading data...")

# (Execute command) - Generate *.fe file by using the vor2fe.exe. The format will be "(numberOfGrains)grains.fe" so that it can be used with the Surface Evolver
command = [(pathToVor + "vor2fe"), '-s101', ('-n' + str(numberOfGrains))]
exitStatus, stdOutput, stdError = execute(command)

# "stdOutput" contains the commands generated by the vor2fe file
# Write the texts generated by the vor2fe executable to a new *.fe file
file = open(pathToFe, 'w+')

# If not 3phase-iso you need to append these at the beginning of the *.fe file before append the vor2fe generated vertices, edges, etc.
mm4Commands = None
if typeOfSimulation == '3phase-mm4-cos':
    mm4Commands = "mm4-cos-Commands.txt"
elif typeOfSimulation == '3phase-mm4-octagon':
    mm4Commands = "mm4-octagon-Commands.txt"

# If it is not 3phase-iso it will append the commands
if typeOfSimulation != '3phase-iso':
    mm4cosCommands = open(pathToData + mm4Commands, "r")
    commandsToAppend = mm4cosCommands.readlines()
    for line in commandsToAppend:
        file.write(line)

# Either way you need to append the vor2fe generated data
file.write(stdOutput)
file.close()

# Opening the "(numberOfGrains)grains.fe" so that commands for simulating/reporting output will be in the *.fe file
# The commands that are appended to the *.fe file are in the "data" directory - each .txt file contains specific commands for a given simulation
commandsForFeFile = open(pathToData + typeOfSimulation + ".txt", "r")
commandsToAppend = commandsForFeFile.readlines()
commandsForFeFile.close()

# Append the commands to the *.fe file used (just created earlier)
feFile = open(pathToFe, "a")
for line in commandsToAppend:
    feFile.write(line)

# -----#########################################################-----#
# EDIT THE *.fe FILE BASED ON RAPPTURE OPTIONS #

# Turn on color before exporting the "before" images
if enableColor:
    feFile.write("show facet where 1;\n")

# Path to the image before
pathToImageBefore = "grainsBefore.ps"

# Call function - from surfEvolverGUIFunctions - to append a command to generate a post script file to the desired path
appendPostScriptCommand(feFile, pathToImageBefore)

# Toggle the color off so that the simulation runs faster
if enableColor:
    feFile.write("show facet where 0;\n")

# Write number of iterations for simulation
feFile.write("gogo " + str(numberOfIterations) + "\n")

# Show colors in pictures if it is wanted
if enableColor:
    feFile.write("show facet where 1;\n")

# EDIT THE *.fe FILE BASED ON RAPPTURE OPTIONS #
# -----#########################################################-----#

# Make a POST-SCRIPT file of the grains after the simulation
pathToImageAfter = "grainsAfter.ps"

# Call function - from surfEvolverGUIFunctions - to append the command to generate a post script file to the desired path
appendPostScriptCommand(feFile, pathToImageAfter)

# Append command to exit the program to the *.fe file
feFile.write("q\nq\n")

# End of editing the *.fe file
feFile.close()

Rappture.Utils.progress(50, "Simulating growth... (this may take a while)")

# Run the Surface Evolver and tell it to execute the script found at "pathToFe"
# (Execute command) - run the surface evolver and pass the *.fe file containing the commands for the Surf. Evolver
command = ['evolver', pathToFe]
exitStatus, stdOutput, stdError = execute(command)

# Change post-script files to *.jpg so that it can be used with Rappture (Rappture doesn't support *.ps files but it does support *.jpg)
pathToImageBeforePng = "grainBefore.png"
pathToImageAfterPng = "grainAfter.png"

# (Execute command) - pass arguments to convert the image files *.ps->*.png
command = ['convert', pathToImageBefore, pathToImageBeforePng]
exitStatus, stdOutput, stdError = execute(command)
command = ['convert', pathToImageAfter, pathToImageAfterPng]
exitStatus, stdOutput, stdError = execute(command)


#########################################################
# Save output values back to Rappture
#########################################################

Rappture.Utils.progress(75, "Generating graphs... (this may take a while)")

# output the image
io.put('output.image(grainImageBefore).current', pathToImageBeforePng, type='file')

# output the image
io.put('output.image(grainImageAfter).current', pathToImageAfterPng, type='file')

# Fetch data from output file, all of these variables are lists containing corresponding x and y coordinates on a curve, all parameters are vs. number of Surface Evolver iterations of "g"
# "getGeneralData()" is a function defined in - surfEvolverGUIFunctions -
numberOfGreenGrains, numberOfRedGrains, numberOfWhiteGrains, totalNumberOfGrains, grainIterations, greenArea, greenIterations, redArea, redIterations, whiteArea, whiteIterations, avgEnergyPerGrain = getGeneralData()

# Output the data for the curve: Area vs. Number of Surface Evolver Iterations "g"
# grouping them together puts it on one curve
# set the color for each curve
io['output.curve(greenAreaCurve).component.xy'] = (greenIterations, greenArea)
io['output.curve(greenAreaCurve).about.style'] = '-color green'
io['output.curve(redAreaCurve).component.xy'] = (redIterations, redArea)
io['output.curve(redAreaCurve).about.style'] = '-color red'
io['output.curve(whiteAreaCurve).component.xy'] = (whiteIterations, whiteArea)
io['output.curve(whiteAreaCurve).about.style'] = '-color blue'

# Output the data for the curve: Grains vs. Number of Surface Evolver Iterations "g"
# grouping them together puts it on one curve
# set the color for each curve
io['output.curve(numberOfGreenGrainsCurve).component.xy'] = (grainIterations, numberOfGreenGrains)
io['output.curve(numberOfGreenGrainsCurve).about.style'] = '-color green'
io['output.curve(numberOfRedGrainsCurve).component.xy'] = (grainIterations, numberOfRedGrains)
io['output.curve(numberOfRedGrainsCurve).about.style'] = '-color red'
io['output.curve(numberOfWhiteGrainsCurve).component.xy'] = (grainIterations, numberOfWhiteGrains)
io['output.curve(numberOfWhiteGrainsCurve).about.style'] = '-color blue'
io['output.curve(totalNumberOfGrainsCurve).component.xy'] = (grainIterations, totalNumberOfGrains)
io['output.curve(totalNumberOfGrainsCurve).about.style'] = '-color orange'

# Fetch data from output file, all of these variables are lists containing corresponding x and y coordinates on a curve, all parameters are vs. number of Surface Evolver iterations of "g"
# "getSideData()" is a function defined in - surfEvolverGUIFunctions -
avgSidesGreen, avgSidesRed, avgSidesWhite, avgSidesTotal, avgIterations = getSideData()

# Output the data for the curve: Avg Sides per Grain vs. Number of Surface Evolver Iterations "g"
io['output.curve(avgSidesGreenCurve).component.xy'] = (avgIterations, avgSidesGreen)
io['output.curve(avgSidesGreenCurve).about.style'] = '-color green'
io['output.curve(avgSidesRedCurve).component.xy'] = (avgIterations, avgSidesRed)
io['output.curve(avgSidesRedCurve).about.style'] = '-color red'
io['output.curve(avgSidesWhiteCurve).component.xy'] = (avgIterations, avgSidesWhite)
io['output.curve(avgSidesWhiteCurve).about.style'] = '-color blue'
io['output.curve(avgSidesTotalCurve).component.xy'] = (avgIterations, avgSidesTotal)
io['output.curve(avgSidesTotalCurve).about.style'] = '-color orange'

# Output the data for the curve: Avg Energy per Grain vs. Number of Iterations "g"
io['output.curve(avgCurvaturePerGrainCurve).component.xy'] = (grainIterations, avgEnergyPerGrain)
io['output.curve(avgCurvaturePerGrainCurve).about.style'] = '-color blue'

# Clean up the files in the system - e.g. (numberOfGrains)grains.fe, xxx-report.txt, and sides.txt
command = ['rm', pathToFe, pathToOutputFileX, pathToOutputFileSide]
exitStatus, stdOutput, stdError = execute(command)

Rappture.Utils.progress(100, "Done")

io.close()
sys.exit()