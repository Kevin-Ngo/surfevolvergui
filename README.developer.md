# README for *GUI for Surface Evolver* Developers

**Use subversion! (NCN URE Students)**

At Purdue University, specifically nanoHUB, Subversion is how they keep code organized and accessible in the future. Subversion is similar to GitHub (allows you to control versions and is also a safe place to back-up your code), but in the case that you are not familiar with it, click [here](https://www.thegeekstuff.com/2011/04/svn-command-examples/) for a brief overview. If you want to checkout the code to a workspace, run the following command.

~~~~~
svn checkout https://nanohub.org/tools/surfevolvergui/svn/trunk surfevolvergui
~~~~~

## A Quick Rappture Overview
When you start developing with Rappture, you'll see that it is all ran on the Linux Bash. The "documentation" can be found [here](https://nanohub.org/infrastructure/rappture/wiki/Documentation); however, in my experience, the online YouTube series for Development with Rappture and the "examples" page were far more useful. This series can be watched [here](https://youtu.be/2g7lgOr8SJ4) and for the examples, if you have a nanoHUB account with permissions to view the page, it can be found [here](https://nanohub.org/infrastructure/rappture/browser/trunk/examples/zoo?order=name).

### Commands:
~~~~
rappture
~~~~
1. "**rappture**": invokes the "tool.xml" file that is in the current directory.
 
~~~~
rappture -builder
~~~~
2. "**rappture -builder**": invokes the "tool.xml" file, but in "building/editing" mode. This is useful for adding features to the GUI. *Note: When you exit/save the builder,* **__DO NOT CHECK THE BOX TO SAVE THE \*.PY FILE, THE DEFAULT PATH WILL OVERRIDE THE MAIN.PY WITH A BAREBONES SKELETON PROGRAM__**

## A Quick Surface Evolver Overview (For Grain Growth)
When you start the *Surface Evolver*, it expects an input which is the path to an \*.fe file. In the case for modeling the growth of grains, this initial \*.fe file is created through the executable "vor2fe.exe". The vor2fe.exe takes arguments for how many grains to create, which it then distributes the 2D plane into that many "grains" (facets). On the shell, it is invoked using the following command. *Note: "N" represents the number of grains that you want, put any number.*

~~~~
./vor2fe -s101 -nN > "someText.txt"
~~~~

After creating the "grains", it marks the grains by color, randomly, then dissolves edges between facets of the same color. Lastly, it also creates a function that will iterate the surface (move to lower energy) and report qualities about the surface such as area, time, number of sides, etc. The function is called "gogo" and it is important because it is how you simulate growth and also collect data into a text file.

The *Surface Evolver* is an interactive program, meaning you can continuously pass arguments to it, such as "GOGO" (user-defined).

**For *Surface Evolver* commands, view the documentation, which can be found [here](http://facstaff.susqu.edu/brakke/evolver/html/evolver.htm). It is dense; however, very helpful! You will most likely be reading this the majority of the time. Refer to the documentation for all *Surface Evolver* commands**

### Commands (User-Defined):
~~~~
gogo N
~~~~
1. "**gogo**": is a function that when ran in the *Surface Evolver*, will cause the grains to grow for "N" iterations. In addition to simulating the grain growth, "gogo" is the main function where data is "told" to be collected. Typically when you want to create a feature to record a new type of data not already being collected, create a function, then add it to "gogo". "gogo" currently outputs data from the surface to text files named: "xxx-report.txt" and "sides.txt".

<img src="/docs/gogoDemo.png"/>

2. "**reportSize**": is a helper function that is used inside "gogo". This function will output the number of sides each grains has at the time that the function is called to a file called "sides.txt".

3. "**reportEnergyArea**": is a helper function that is used inside "gogo". This function function will output the area of each grain as well as the energy of the system at the time that it is called to a file called "xxx-report.txt".

## The *GUI for Surface Evolver* Files
To run/edit the *GUI for the Surface Evolver* there are many files used.
1. **main.py**: The *GUI for the Surface Evolver* is backended with Python. This main.py uses modules/libraries that are needed for Rappture and are already imported. To make computational changes in the tool, edit this file and also the reportAndSImulateCommands.txt file.

**Note: When running the *GUI for the Surface Evolver* in a local workspace, edit the path for "pathToData" and "pathToOutputDirectory". The "pathToData" is a string containing the address of the directory where data to run the program correctly is stored, such as what commands to append to the \*.fe file and also the images of grains to output. The "pathToOutputDirectory" is a string containing the address of the directory where the collected data from the *Surface Evolver* is outputted.**

<img src="/docs/changePaths.png"/>

**Here is the general flow of the script:
Invoke the Rappture Interface &rarr; Assign in-line variables to the user-parameters &rarr; Generate an \*.fe file &rarr; Append commands to the newly generated \*.fe file &rarr; Append a command to capture the image before simulating &rarr; Append "gogo -N" based on parameters &rarr; Run simulation &rarr; Capture image after &rarr; Extract information from the text files that were created by *Surface Evolver* &rarr; Output the information to Rappture**

2. **tool.xml**: This is the main file used with Rappture. Typically you do not have to edit this file. When you make a change to the tool using "rappture -builder", simply save the changes and check the box to update the "tool.xml".

3. **reportAndSimulateCommands.txt**: This is a text file that is appended to a *.fe file. This text file should be updated accordingly if you want to add new features. Define a new function here then add it to the "gogo" function to run it.

## Final Notes

This project was worked on by Kevin Ngo in collaboration with Lucas Robinson and Dr. John Blendell in the Summer of 2018 at Purdue University for the NCN URE program.
