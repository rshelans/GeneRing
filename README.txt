##########################################################################
##################Chromatin Ring Analysis Python Suite####################
##########################################################################

#####################
##Table Of Contents##
#####################
TraceAnalyzer.py
	Program for analyzing Groups of Trace files.
Trace.py
	Module which knows how to do analysis on raw trace files.
Chromatin.py
	Module which represents individual chromatin rings and ring 
	manipulation.
TraceModel.py
	Module which communicates between Trace.py and TraceView.py.
	Is used by the TraceAnalyzer.py
TraceView.py
	Module which is the front end of the TraceAnalyzer.py.

####################
####Instructions####
####################
The following is a set of instrcutions explaining how to do analysis of
chromatin molecules using the TraceAnalyzer.py.

PROLOGUE:

If python3.5 is not installed on your computer instal python5.5 (Check).
If not installed, 
I recommend installing python3.5 using the open data science platform. 

ANACONDA
https://www.continuum.io

I recommend this because *this* suite maintains a couple dependencies
which are taken care of by installing using ANACONDA. Make life easy.

Check to make sure the correct version of python is installed.
run either CMD or terminal (windows/mac). Type python hit enter. 
The first line that is printed after hit enter should be, *like*:
Python 3.5.*  | Anaconda *.*.* (*) ******************************


1) Follow all the tracing guidlines and trace a number of Molcules.
Record the Image Resolution. For convinience put them in the same folder.

2) Open up a Command Prompt set the working directory to the directory
where this suite is located. (using command cd) Alternatively, Add
this suite to your computers *path*.

3) Run the TraceAnalyzer.py by typing TraceAnalyzer.py or
python Traceanalyzer.py remember to give the program the perameters
it needs. Usage Below.

usage: TraceAnalyzer.py [-h] [-v] [-t THRESHOLD] [-l LENGTH] [-s SMOOTH]
                        [-u USER] [-o OUT_PATH] [-p] [-i IMGRES]
                        files [files ...]

Analyze Chromatin Ring Trace Files.

positional arguments:
  files                 The list of trace files you wished to be analyzed.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -t THRESHOLD, --threshold THRESHOLD
                        The threshold distance required for two points on
                        opposing strands to be considered as part of a linker.
  -l LENGTH, --length LENGTH
                        The expected length in basepairs for the DNA molecule.
  -s SMOOTH, --smooth SMOOTH
                        The # of coordinates to include in a sliding window of
                        the average distance between points on opposing
                        strands.
  -u USER, --user USER  The name of the person who completed the trace and is
                        using this software: for record keeping.
  -o OUT_PATH, --out_path OUT_PATH
                        The name of the folder you wish to save the output to.
  -p, --plotless
  -i IMGRES, --imgres IMGRES
                        The image resolution of the raw image used to make the
                        trace.

Assuming I all of my trace files are 
	1  ) By the same tracer(I am the Tracer: OMEGA).
	2  ) Have the Same Image Resolution (ex.  0.3 nm/au)
	3  ) In The directory(folder)  C:\Users\Robert\trace_files\
	3.5) Each trace file is names something like stuff001.txt.
	4  ) I want to use a default threshold of 4.
	5  ) The molecules are expected to be 2246bp long.
	6  ) I want to use a smooth factor of 10.
	7  ) I want to save all the output files to the same directory
											 where the traces are.
	8  ) I want it to autocompute all traces it can without 
												plotting


TraceAnalyzer.py  C:\Users\Robert\trace_files\*.txt -t 4 -l 2246 -s 10  -i 0.3 -u OMEGA -p 

Notice: These are the default threshold length and smooth settings. The following
does the same.

TraceAnalyzer.py  C:\Users\Robert\trace_files\*.txt -u OMEGA -i 0.3 -p

I reccomend the following work flow. Run Without plotting the -p parameter.
Then rerun w/o the -p parameter and *FIX* all the molecules it could not figure out.
Using either the seek slider bar or the next and prev buttons. Remember to save. 
Any edited molecules. The save indicater will appear after save is complete.
Incorrect Molecule fits will be printed to the terminal screen in the -p setting. 
However, I recommend looking at every trace at least once. You must hit save to save 
molecule trace

#########################
######HOW TO FIX#########
#########################
4) How to *FIX* molecules.
You have  a number of metrics/plots/toggles at your disposal.

HOW TO TELL A MOLECULE IS WRONG
EASY mode:
-If no segment calls appear with bp sizes in the bottum right above the buttons
-if BUBBLE CALL INFO contains None or NAN the molecule is wrong.
-If in the trace window right (rainbow one) the two strands are different colors.
-If In the left plot Linker/Bubble (red/blue) the segments dont match color 
-If In the left plot Linker/Bubble the midpoint(Large black point) is not close to the
													actual midpoint of the molecule.
-If the distance plot does not look symmetrical.
HARD mode:
-If the Image resolution does not match the predicted Image Resolution within
										tolerance (probably a couple nms).
-If the steadyness or resolution is above tolerance (~50 bp^2).
-If a particular traces metrics look vastly different then everyother from the 
															same set.

HOW TO FIX A INCORRECT MOLECULE CALL
Typically, The easiest way to fix a molecule is either increase or decrease the threshold to remove or add linker regions that where called incorrectly. Use the
Threshold slider and observe both the distance plot
(blue line plot with red and cyan perpendicular lines). Ensure that every linker
region has a mate in the opposing strand. When this looks good.
You can either set the midpoint yourself or click find midpoint to automate the process. Do this iteratively until a molecule does not *look* wrong anymore.

**recommendation: Set the threshold lower until  all linker look called. Setting it higher often leads to skipped linker regions.

**Recommendation: If molecule is too hard just skip it.

DONT FORGET TO CLICK SAVED AND WAIT FOR IT TO SAY SAVED.  


5) Well Done you have created molecule data. Now go an do some analysis!

######################
##Metric Definitions##
######################

Length- length of molecule in bp
EdgeIgnore- number of coordinates from the midpoint that are ignored.
Treshold- the minimum distance in AU a coordinate needs to be to be considered a 
																		linker.
Smooth- the size of a sliding average window in coordinates
					(Smooths distances between coordinates to prevent jumps).
Midpoint- coordinate which is considered the very middle of the molecule
Image Resolution- Resolution of the TEM Image Trace was made from.
Pred. Image Res.- Guess at the resolution of the TEM iage the Trace was made from.
Trace Scale     - Base Pairs per coordinate distance of trace.
Trace Resolution- A measure of the speed of the Trace.
Trace Steadiness- A measure of how consistant trace speed was.
segments- Number of transitions from linker to bubble in the molecule.
Symmetry Score- How close the midpoint is to the midpoint in terms of distance of the 
															molecule.
Linker Score- The mean square distance between start and end sites for linkers.
Region Score- The mean square difference in size between strands of a region.

#######################
##Option Descriptions##
#######################
SLIDERS
seek- Change the current molecule to any molecule reachable by index.
smooth- Sets the smooth parameter.
threshold- sets the threshold parameter.
midpoint- sets the midpoint of the current molecule.
BUTTONS
find midpoint- given the parameters attempts to find the midpoint of the current 
																		molecule.
<<prev<<- sets the current molecule to the previous molecule in index.
>>next>>- sets the current molecule to the next molecule in index
save    - saves the current molecule if it can with the current settings.






