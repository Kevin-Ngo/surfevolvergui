'''
    Project: GUI for the Surface Evolver - Polycrystalline Grain Growth - functions module
    Author: Kevin Ngo
    Program: Purdue University - Network for Computational Nanotechnology
    Date: Summer 2018
    Collaborators: Lucas Robinson
    Faculty Advisor: John Blendell
    Acknowledgements: Jean Taylor and Ken Brakke
'''


# This module is to store all functions that are used in the Surface Evolver GUI main.py
# Put all functions to process data or manipulate the *.fe file in this module so that in the main.py you just call functions/run the simulation


# Function to append the commands to the *.fe to generate a ps file
def appendPostScriptCommand(feFile, desiredPathForPS):
    feFile.write("P\n")
    feFile.write("3\n")
    feFile.write("y\n")
    feFile.write("\n")
    feFile.write(desiredPathForPS + "\n")


# Function to read the output text file, xxx (contains area, number of grains and energy), and process the data by inserting it into lists that can then be outputted to Rappture
def getGeneralData():
    # Open file for reading
    outputFile = open("xxx-report.txt", 'r')
    
    # Read all the lines in the file and store it into a List named lines
    lines = outputFile.readlines()
    
    # These are all "x" & "y" lists for the graph (time always being the x)
    numberOfRedGrains = []
    numberOfGreenGrains = []
    numberOfWhiteGrains = []
    totalNumberOfGrains = []
    grainIterations = []
    redArea = []
    greenArea = []
    whiteArea = []
    redIterations = []
    greenIterations = []
    whiteIterations = []
    totalArea = 0
    currentIteration = None
    totalGrains = None
    avgEnergyPerGrain = []
    
    # See output file (xxx-report.txt) to see the pattern of how to process the information
    for line in lines:
        line = line.split() ## split each "line" into it's individual strings
        if line[0] == "Iterations:":
            currentIteration = float(line[1])
            numberOfRedGrains.append(int(line[3]))
            numberOfGreenGrains.append(int(line[5]))
            numberOfWhiteGrains.append(int(line[7]))
            totalGrains = int(line[3]) + int(line[5]) + int(line[7])
            totalNumberOfGrains.append(totalGrains)
            grainIterations.append(currentIteration)
        elif line[0] == "Total":
            totalEnergy = float(line[2])
            avgEnergyPGrain = totalEnergy / totalGrains
            avgEnergyPerGrain.append(avgEnergyPGrain)
        # can use grainTime because every calculation per grain & grain times are the same
        elif line[0] == "Green":
            if totalArea > 0:
                whiteArea.append(totalArea)
                whiteIterations.append(currentIteration)
                totalArea = 0
        elif line[0] == "Red":
            greenArea.append(totalArea)
            greenIterations.append(currentIteration)
            totalArea = 0
        elif line[0] == "White":
            redArea.append(totalArea)
            redIterations.append(currentIteration)
            totalArea = 0
        else:
            lowerEnd = line[0]
            upperEnd = line[2]
            frequency = line[3]
            midPoint = (float(lowerEnd) + float(upperEnd)) / 2
            totalArea = totalArea + (midPoint * int(frequency))
    return numberOfGreenGrains, numberOfRedGrains, numberOfWhiteGrains, totalNumberOfGrains, grainIterations, greenArea, greenIterations, redArea, redIterations, whiteArea, whiteIterations, avgEnergyPerGrain


def getSideData():
    # Open file for reading
    outputFile = open("sides.txt")
    
    # Grab lines
    lines = outputFile.readlines()
    
    # Variables in this function
    avgSidesRed = []
    avgSidesGreen = []
    avgSidesWhite = []
    avgIterations = []
    numberOfGrains = None
    totalSides = 0
    average = 0.0
    avgSidesTotal = []
    totalNumberOfGrains = 0
    totalSidesForAllColors = 0
    averageTotal = 0.0
    
    for line in lines:
        line = line.split()
        if line[0] == "Iterations:":
            avgIterations.append(float(line[1]))
            if totalSidesForAllColors > 0:
                averageTotal = float(totalSidesForAllColors) / totalNumberOfGrains
                avgSidesTotal.append(averageTotal)
                totalSidesForAllColors = 0
                totalNumberOfGrains = 0
        elif line[0] == "red:":
            # If at red this means, we were processing white (except first pass)
            if totalSides > 0:
                average = float(totalSides) / numberOfGrains
                avgSidesWhite.append(average)
                average = 0
                totalSides = 0
            if numberOfGrains == 0:
                avgSidesWhite.append(0)
            numberOfGrains = int(line[1])
        elif line[0] == "green:":
            # If at green this means, we were processing sides for red
            if numberOfGrains == 0:
                avgSidesRed.append(0)
            else:
                average = float(totalSides) / numberOfGrains
                avgSidesRed.append(average)
                average = 0
                totalSides = 0
            numberOfGrains = int(line[1])
        elif line[0] == "white:":
            # If at white this means, we were processing sides for green
            if numberOfGrains == 0:
                aavgSidesGreen.append(0)
            else:
                average = float(totalSides) / numberOfGrains
                avgSidesGreen.append(average)
                average = 0
                totalSides = 0
            numberOfGrains = int(line[1])
        else:
            totalSides += int(line[0])
            totalSidesForAllColors += int(line[0])
            totalNumberOfGrains += 1

    # append and calculate last values (after the file reaches end still need to calculate values)
    averageTotal = float(totalSidesForAllColors) / totalNumberOfGrains
    avgSidesTotal.append(averageTotal)

    return avgSidesGreen, avgSidesRed, avgSidesWhite, avgSidesTotal, avgIterations
