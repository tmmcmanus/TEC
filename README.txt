████████╗███████╗ ██████╗
╚══██╔══╝██╔════╝██╔════╝
   ██║   █████╗  ██║     
   ██║   ██╔══╝  ██║     
   ██║   ███████╗╚██████╗
   ╚═╝   ╚══════╝ ╚═════╝
                         

TEC (Transformation Electromagnetics Code) solves for the material parameters called for by various transformation optics (TO)devices.
There are two different source files found here thus far, along with a number of sample surfaces and sample volumes.  

The first: TEC_2D can be used to create anisotropic surface wave cloaks and illusion devices using any of the structured meshes found in /Sample_Surfaces.  The second: TEC_3D is a very basic code which effectively compares elemental volumes of a virtual and physical space in an attempt to test a Volumetric quasi-conformal transformation optics (VQCTO) technique.  The meshed volumes that TEC_3D operates on can be found in /Sample_Volumes.

Lastly, as these are MATLAB files, they do require a version of MATLAB to run (TEC_2D and TEC_3D has been tested on MATLAB 2015a, but they should also run on earlier versions of MATLAB as well).  Hopefully in the future these source files will be converted into .py.  
