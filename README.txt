##############################################################################
#                                                                            #
#                         NIFTYREC TOMOGRAPHY TOOLBOX                        #
#                                                                            #
#                     http://niftyrec.scienceontheweb.net                    #
#                                                                            #
##############################################################################

------------------
1 PACKAGE CONTENTS
------------------

NiftyRec is a software toolbox for Tomographic image reconstruction. 
NiftyRec is written in C and computationally intensive functions have a 
GPU accelerated version based on NVidia CUDA. 
NiftyRec includes a Matlab Toolbox and a Python Package that access 
the low level routines, hiding the complexity of the GPU accelerated 
algorithms. 

##############################################################################

---------------
2 RELEASE NOTES
---------------

This repository is a porting of NiftyRec 2.3.2 as it is available on it website.
Some changes have been made to the source code, so that it is no longer backward
compatible, but it now fully supports the latest stable version of occiput 
(in its current stable version available via pip) and the beta of occiput-suite
now under development at https://github.com/mscipio/occiput-suite.

This new version of NiftyRec will be called 3.0.0, as to underline the dramatic changes
wrt the previous release.

--------------------

For the sake of completeness, and to ackowledge the work of the original developers,
the original changelog will follow:


-- NiftyRec 1.0 (Rel. Nov 2010)-- 
    * Rotation-based projection and backprojection with depth dependent point spread function
    * GPU accelerated version working ok (non accelerated version not working)
    * mex-based Matlab Toolbox
    * Matlab inline documentation
    * Matlab functions for reconstruction

-- NiftyRec 1.1 (Rel. Jan 2011) -- 
    * Python interface
    * Documentation: Programmer's manual

-- NiftyRec 1.2 (Rel. Mar 2011) -- 
    * Bug Fixes

-- NiftyRec 1.3 (Rel. May 2011) --  
    * Embedded NiftyReg
    * Demo

-- NiftyRec 1.4 (Rel. Sep 2011) --
    * Ray-cating based Projector for Transmission Tomography
    * Graphical User Interface
    * Documentation update

-- NiftyRec 1.5 (Rel. Jan 2012) --
    * Faster projection and back-projection (separable convolution) 
    * Fast Fisher Information Estimation
    * Several bug fixes and improved Matlab demos and inline documentation

-- NiftyRec 1.6 (Rel. May 2012) --
    * Command line reconstruction tool for SPECT 
    * 3D Slicer4 plugin (based on Slicer execution model)
    * Independent from Matlab Toolboxes (now only requires basic Matlab install) 
    * Emission tomography handles non-cubic volumes 
    * Emission tomography handles 2D 
    * Simind interface: import and reconstruct Simind simulations 
    * PET-Sorteo interface: import and reconstruct PET-Sorteo simulations 
    * Matlab demos for Simind and PET-Sorteo reconstruction 
    * Barebone CTest based unit testing
    * NiftySeg bindings 
    * Barebone NiftyReg bindings 

-- NiftyRec 1.6.5 (Rel. Jan 2013) --
    * Fixed bug with multiple rotation angles in emission tomography

-- NiftyRec 1.6.6 (Rel. Feb 2013) --
    * Compatible with CUDA 5.0

-- NiftyRec 1.6.7 (Rel. Mar. 2013) --
    * Comaptible with CUDA compute capability 3.0
    * Fixed bugs in the Python interface 
    * Fixed bugs in the Python GUI example 
    * Added simple python example 
    * Added sample dataset for the Python GUI example

-- NiftyRec 1.6.8 (Rel. Mar. 2013) -- 
    * Fixed compilation issues under Windows

######################################################################################

-----------------------------------
4 BUILD AND INSTALL FROM SOURCE 

See the Wiki on the website of 
NiftyRec for detailed instructions
-----------------------------------

NiftyRec is based on the CMake cross-platform build system. 
As a design choice NiftyRec does not have any external dependencies when 
compiled with the basic options. Optional dependencies are the CUDA runtime 
libraries, Matlab (mex and mx) and the Python interpreter. 
CMake simplifies the build process across the supported platforms. 
Further details are in the Programming Manual of NiftyRec.

--Linux and MAC

   Install ccmake or cmake-gui. Download and uncompress source. cd to the 
   project main directory, create here a folder 'build', cd to that folder 
   and run cmake: 

   ..$ mkdir build
   ..$ cd build
   ..$ ccmake ..   (or cmake-gui ..)

   Select options, set the BUILD_TYPE to Release or to Debug and set all the required 
   paths for additional dependencies if you selected any of the options. 
   Configure and Generate. Now quit ccmake/cmake-gui and build
  
   ..$ make build

   In order to create an installation package with CPack run make with option 'package'

   ..$ make package

   In order to install NiftyRec in the system run make with option 'install'

   ..$ sudo make install  

   or install the package created with 'make package'.

--Windows

   Install cmake. Download and uncompress source. Open the source directory 
   with Windows Explorer and create here a new folder 'build'. Launch CMake 
   and open CMakeLists.txt from the main directory of NiftyRec.
   Select options, set the BUILD_TYPE to Release or to Debug and set all the required 
   paths for additional dependencies if you selected any of the options. 
   Configure and Generate. Browse the 'build' directory and double click 
   on the Visual Studio project file. Click Compile button in Visual Studio. 
   Create self-extracting installer by compiling the corresponding target in Visual 
   Studio. 

##############################################################################

-----------------------------------
5 USE

See the Wiki on the website of 
NiftyRec for detailed instructions
-----------------------------------

--Matlab Toolbox
   Launch Matlab. Add path to NiftyRec Toolbox, e.g.:

   >> addpath '/usr/local/niftyrec/matlab'

   Or add permanently by clicking on File->Add path. 
   The path NiftyRec Toolbox is set as an option in CMake. 
   It defaults to '/usr/local/niftyrec/matlab' in Linux and MAC OS, in  Windows C:/ProgramFiles/NiftyRec
   Open Matlab help and click on NiftyRec Toolbox 
   to visualize the documentation of NiftyRec Toolbox. Execute and edit et_osem_demo.m to get started. 

--Python Extension

   The NiftyRec Python module is installed amongst the site-packages 
   by the binary installers and by 'make install'. 
   Open the Python interpreter and import NiftyRec

   >>> from NiftyRec import NiftyRec
   
--C API

   Libraries and headers are installed by the binary installer and by 'make install'. 
   Their location is an option in CMake at build time. It defaults to '/usr/local/lib' 
   and '/usr/local/include' under Linux and Mac OS and to 'C:/ProgramFiles/NiftyRec/lib' 
   and 'C:/ProgramFiles/NiftyRec/include' under Windows.
 
##############################################################################

---------
6 LICENSE
---------

Copyright (c) 2009-2012, University College London, United-Kingdom
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.
Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

Neither the name of the University College London nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
THE POSSIBILITY OF SUCH DAMAGE.

##############################################################################

---------
7 CONTACT
---------
Please contact Stefano Pedemonte (s.pedemonte@cs.ucl.ac.uk)

##############################################################################

------------
8 REFERENCES
------------
[1] S.Pedemonte, A.Bousse, K.Erlandsson, M.Modat, S.Arridge, B.F.Hutton, 
S.Ourselin. GPU Accelerated Rotation-Based Emission Tomography Reconstruction. 
NSS-MIC 2010.

##############################################################################
##############################################################################


