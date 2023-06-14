
                        CONTACT version "trunk"

Overview
========

CONTACT is an advanced simulation program for the study of three-dimensional
frictional contact problems, such as occur between wheel and rail, in roller
bearings and in offset printing devices. It implements the theories for
rolling contact by Dr.ir. E.A.H. Vollebregt and Prof.dr.ir. J.J. Kalker of
Delft University of Technology. 

CONTACT is provided by Edwin Vollebregt, through his company Vtech CMCC.
 * An "open source version" is provided for the full CONTACT program and
   CONTACT library version. This uses source the Apache v2.0 license as shown
   in "LICENSE.txt". 
 * A "convenience version" is provided at a service fee: with additional
   quality checking, documentation, with an automated installer, and
   detailed support.
See the section "Software products" at the website www.cmcc.nl.

CONTACT is intended for concentrated contact problems:
 * from given profiles, determine where contact occurs;
 * identifying the size and shape of contact areas and the pressures acting;
 * accounting for creepage and friction, determine the tangential surface
   stresses (tractions), and regions with adhesion and (micro-)slip;
 * calculation of elastic displacements and subsurface stresses in the bodies'
   interiors;
 * for stationary and instationary, rolling and sliding problems;
 * with Hertzian and non-Hertzian, "smooth edged" geometries, such that the
   contact is "concentrated" in a small, almost flat part of the bodies' 
   surfaces;
 * for massive homogeneous bodies of elastic and viscoelastic materials.

CONTACT aims to be the worlds fastest detailed contact model: more detailed 
than fast approximations such as analytical formulae or the FASTSIM approach,
and faster than more complete models such as based on Finite Elements.


Getting started
===============

The CONTACT program itself is an old-fashioned text-based application,
reading a plain input-file, computing the cases, producing its various
outputs. This is not a big deal once you are on the way, but may be
confusing at first. To alleviate this we created a simple CONTACT GUI
(graphical user interface).

A. Using the CONTACT GUI.
-------------------------
The CONTACT GUI is a small wrapper around the CONTACT program. On Windows it
is started via its desktop icon or the Windows Start menu. On Linux you use
the script `start_gui.sh'.

To use the GUI with the examples you may use the following steps:

 0. Make sure that CONTACT is installed appropriately, see section
    `Installation' below.

 1. Start the CONTACT GUI (on Linux: by typing `start_gui.sh' in a terminal
    window). 

 2. Press the "Select Input File" button. Navigate to the directory where 
    the examples reside, e.g. "My Documents\CONTACT\examples-v22.2", and 
    select the input file "cattaneo.inp".

 3. Start the simultation by pressing the button "Run Experiment".
    This will start the CONTACT program with experiment name "cattaneo".
    This example is described in the User Guide, Section 5.1.

 4. Shut down the CONTACT GUI, by selecting menu item File -> Quit.

 5. Start Matlab, within Matlab, change to the examples directory.

 6. Within Matlab, make sure that the directory \Program Files\Vtech CMCC\
    contact_v22.2\matlab directory is added to the Matlab
    search path, see the section `Installation' below.

 7. Read the CONTACT results for the Cattaneo example into Matlab with the
    command "s=loadcase('cattaneo');". Check the description of the example
    in the User Guide for further information on the Matlab scripts provided.


B. Using the old-fashioned command-line interface.
--------------------------------------------------
As an alternative to the GUI described above, you may also use a command 
prompt window and execute MS-DOS commands to use the program. For instance, 
to work with the examples you may use the following steps:

 0. Make sure that CONTACT is installed appropriately, see section
    `Installation' below.

 1. Open a command prompt window (e.g. Start -> Windows System -> Command
    Prompt).

 2. Within this window change to the directory where the examples reside:
    cd "My Documents\CONTACT\examples-v22.2 (adjust this command to the 
    location where you installed the software).

 3. Run the CONTACT program, using command "<CONTACTDIR>\contact.exe 2 cattaneo"
    (replace <CONTACTDIR> with the location where you installed the software).
    This will start the program in batch-mode (2) using experiment name
    "cattaneo". This example is described in the User Guide.

 4. Start Matlab, within Matlab, change to the examples directory.

 5. Within Matlab, make sure that the contact_v22.2\matlab directory is
    added to the Matlab search path, see section `Installation' below.

 6. Read the CONTACT results for the Cattaneo example into Matlab with the
    command "s=loadcase('cattaneo');". Check the description of the example
    in the User Guide for further information on the Matlab scripts provided.

Some tricks to permanently set the path variables and simplify these steps
are given below.


Dependencies
============

Usage of the CONTACT GUI requires that the Java runtime environment be
installed (www.java.com). Further the file extension ".jar" should be
associated to the program javaw.exe in the Java installation directory.

Usage of the plot programs requires a license to the (commercial) Matlab 
package.

For using the CONTACT GUI the environment variable CONTACTDIR must be
defined, pointing to the folder <CONTACTDIR> where CONTACT is installed.
You may want to add the bin-directory to your search path as well. 

 * On Windows, add permanent settings:

        Control Panel -> System -> Advanced -> Environment variables -> Add
        CONTACTDIR = <CONTACTDIR>
        PATH = %CONTACTDIR%\bin;%PATH%

 * On Linux, in bash, sh, ksh shells, add in your .profile or .bashrc:

        export CONTACTDIR="<CONTACTDIR>"
        export PATH=${CONTACTDIR}/bin:${PATH}
        export LD_LIBRARY_PATH=${CONTACTDIR}/bin/linux64:${LD_LIBRARY_PATH}

 * On Linux, in tcsh and csh shells, add in your .cshrc file:

        setenv CONTACTDIR "<CONTACTDIR>"
        setenv PATH ${CONTACTDIR}/bin:${PATH}
        setenv LD_LIBRARY_PATH ${CONTACTDIR}/bin/linux64:${LD_LIBRARY_PATH}

You may want to add the matlab-directory to the Matlab search path as well.

 * within one Matlab session, enter the command

        % addpath('<CONTACTDIR>/matlab');

   (Replace <CONTACTDIR> with the actual path to the CONTACT installation.)

 * this command may be put in a file "startup.m" in your work-directory,
   in your personal overall startup-file (MATLAB/startup.m in your
   Documents folder) or in a system-wide configuration-file
   (<MATLAB>/toolbox/local/pathdef.m).


Troubleshooting
===============

In case of problems with this software you may contact us using the feedback
form at the website www.cmcc.nl, or by sending an e-mail to info(at)cmcc.nl.

