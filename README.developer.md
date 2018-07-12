# README for *GUI for Surface Evolver* Developers

**Use subversion! (NCN URE Student)** At Purdue University, specifically nanoHUB, it is how they keep code organized and accessible in the future. Subversion is similar to GitHub (allows you to control versions and is also a safe place to back-up your code), but incase you are not familiar with it, click [here](https://www.thegeekstuff.com/2011/04/svn-command-examples/) for a brief overview. If you want to checkout the code to a workspace, run the following command.

~~~~~
svn checkout https://nanohub.org/tools/surfevolvergui/svn/trunk surfevolvergui
~~~~~

## A Quick Rappture Overview
When you start developing with Rappture, you'll see that it is all ran on the Linux Bash.

### Commands:
~~~~
rappture
~~~~
1. "**rappture**": runs the "tool.xml" file that is in the current directory.
 
~~~~
rappture -builder
~~~~
2. "**rappture -builder**": opens up the "tool.xml" file, but in "building" mode. This is useful for adding features on to the GUI. *Note: When you exit/save the builder,* **__DO NOT CHECK THE BOX TO SAVE THE \*.PY FILE, THIS WILL OVERRIDE THE MAIN.PY WITH A BAREBONES SKELETON PROGRAM__**

## A Quick Surface Evolver Overview (For Grain Growth)
When you start the *Surface Evolver*, it expects an input which is the path to an \*.fe file. In the case for modeling the growth of grains, this initial \*.fe file is created through the executable "vor2fe.exe". The vor2fe.exe takes arguments for how many grains to create, which it then distributes the 2D plane into that many "grains" (facets). On the shell, it is invoked using the following command. *Note: "N" represents the number of grains that you want, put any number.*

~~~~
./vor2fe -s101 -nN > "someText.txt"
~~~~

After creating the "grains", it marks the grains by color, randomly, then dissolves edges between facets of the same color. Lastly, it also creates a function that will iterate the surface (move to lower energy) and report qualities about the surface such as area, time, number of sides, etc. The function is called "gogo" and it is important because it is how you simulate growth and also collect data into a text file.

The *Surface Evolver* is an interactive program, meaning you can continuously pass arguments to it, such as "GOGO" (user-defined).

**For *Surface Evolver* commands (most important commands), view the documentation, which can be found [here](http://facstaff.susqu.edu/brakke/evolver/html/evolver.htm). It is dense; however, very helpful! Refer to this for all *Surface Evolver* commands**

### Commands (User-Defined):
~~~~
gogo N
~~~~
1. "**gogo**": is a function that when ran in the *Surface Evolver*, will cause the grains to grow for "N" iterations. In addition to simulating the grain growth, "gogo" is the main function where data is "told" to be collected. Typically when you want to create a feature to record a new type of data not already being collected, create a function, then add it to "gogo". "gogo" currently outputs data from the surface to text files named: "xxx-report.txt" and "sides.txt".

<img src="gogoDemo.png"\>

2. "**reportSize**": is a helper function that is used inside "gogo". This function will output the number of sides each grains has at the time that the function is called to a file called "sides.txt".

3. "**reportEnergyArea**": is a helper function that is used inside "gogo". This function function will output the area of each grain as well as the energy of the system at the time that it is called to a file called "xxx-report.txt".


