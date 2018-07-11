# GUI for the Surface Evolver
The *GUI for the Surface Evolver* is an online simulation tool written by Kevin Ngo in collaboration with Dr. John Blendell and Lucas Robinson at Purdue University. The original program, named *Surface Evolver*, not to be confused with the *GUI for the Surface Evolver*, is an "interactive program for the modelling of surfaces shaped by various forces and constraints" that is written by Ken Brakke. The *GUI for the Surface Evolver* is made through the *Rappture toolkit* and simplifies the process of simulating the growth of crystal grains growing in 2D by offering users an abstraction from the command line as well as a simple way to extract useful data from the simulation. *Note: This version of the src is not the most up to date, if you intend to add features, visit https://nanoHUB.org/tools/surfevolvergui and check the most updated version first*

## Getting Started
### The Easy Way
To use the tool, you do not need to install it, simply visit https://nanoHUB.org/tools/surfevolvergui, create an account and launch the tool. *Note: If you have trouble accessing the tool after creating an account - submit a ticket and request to use the tool.* 
<img src="docs/launch.png"/>

### The Hard Way
In the case that you want to run the program locally (running the Surface Evolver in some instances is computationally expensive so running the program locally often will speed up the process of simulation) or edit the source, the Rappture toolkit must be installed and it can be found at https://nanohub.org/infrastructure/rappture/wiki/rappture_install_linux. The website offers installation guides for Linux/macOS systems as well as documentation for the Rappture interface.

## Starting the *GUI for the Surface Evolver*
If you created a nanoHUB account you can simply click "Launch Tool" to use it.
If you manually installed the Rappture toolkit, clone this directory and using the shell navigate to the directory named "rappture". Once in this directory type "rappture" to run the tool. The command "rappture" will invoke the "tool.xml" file in the directory.
<img src="docs/directory.png" width="500" height="250"/>
