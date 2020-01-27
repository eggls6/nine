# NINE 
Numerical Integrators for N-body problems written by S.Eggl

## Current version: 2.1                                                                  
                                     
                                                                                         
## REFERENCE PUBLICATION:                                                                            

   2010  Eggl,S. & Dvorak, R.                                                              
   An Introduction to Common Numerical Integration Codes Used in Dynamical Astronomy       
   Lecture Notes in Physics, Volume 790, p. 431-477                                          

## CHANGE LOG:

changes since NINE 1.9
Migrated to Github
Input file format changed  
Yarkovsky integrators deprecated                                                                                                                                                                 
   											    
changes since NINE 1.8
								    
bug in output fixed that caused skipping some outputs or producing multiple instances   
of the same output									    	
											    
changes since NINE 1.7								    
   											    
bug in Raeih fixed (Gauss Radau with EIH equations)                                     
bug in Binary barycentric coordinates fixed						    
											    
chagnes since NINE 1.6.5								    
											    
backward integration added								   
oe velocity bug corrected								    
Raeih forces rewritten, NICE model disabled for EIH 				    
											    
changes since NINE 1.6                                                                  
                                                                                             
fixed hco+hel output bug                                                                    
fixed minor OpenMP bugs                                                                     
fixed close encounter velocity bug for massive massive close encounters                  
fixed time output format bug                                                            
modified Yarkovsky effect: spin orientation stays orthogonal to orbit                   
-> max diurnal effect, changed yarkovky input in config.inn                             
corrected bug in OMP version of Rarya integrator                                                                 
                                                                                           
changes since NINE 1.5 

Exploratory integrators LieSY and LieST dropped                                           
New Yarkovsky module for relativistic Gauss-Radau (6)                                    
Possibility to add names in the config.inn file                                         
                           					                            
changes since NINE 1.4

Starting time included                                                                  
Option for close encounter file output added                                            
											      
changes since NINE 1.3.1 

Relativistic Gauss Radau finished (included Nice Model, massless Particles, OpenMP)     
Relativistic Gauss Radau with Yarkovsky effect added (Rarya)	                    
Improvement of Newton Raphson in module transform 		                             
					                            
changes since NINE 1.3   

Merger flags adapted -> no merging between massless bodies				    
                        								                         
changes since NINE 1.2   

Lie Series SW: massless bodies added to stepping control                                
                         								    
changes since NINE 1.1

OpenMP statements improved                        					    
Switch on/off: progress report, stop when particle beyond cutoff                        
Memory consuming bug in Bulirsch Stoer and Radau force calculation fixed                
Output bug (skipping steps) in Yoshida and Candy fixed                                  
Mass added in binary ouptut files                                                        


## Integrators/Propagators:

CANDY SYMPLECTIC  (symplectic, 4th order)
YOSHIDA SYMMETRIC SYMPLECTIC (8th order)
LIE SERIES INTEGRATOR  with adaptive stepsize control (LieSW)
BULIRSCH STOER  with adaptive stepsize control (BulSt)
GAUSS RADAU with adaptive stepsize control (Radau)
Relativistic (Sol Sys approx) GAUSS RADAU with adaptive stepsize control  
Relativistic (EIH) GAUSS RADAU with adaptive stepsize control  


## REQUIREMENTS: 
		fortran 90 compiler (best: gfortran, no guarantees for ifort, though should work if later than version 11.0)
		make 

## INSTALLATION:

1.) open MAKEFILE

2.) change fortran COMPILERNAME to the compiler of your choice (F=gfortran/F=g95/F=whatever)
       
(2.a) optionally change Compilerflags (CFLAGS=-O -Wall) WARNING: DO NOT USE -static together with -OpenMP !

(2.b) choose folder of the source .f90 files (default: ./src/) 

3.) exit MAKEFILE

4.)type MAKE in shell

5.) open CONFIG.INN to customize input

6.) execute NINE.EXE

(6.a) optionally the MAKE clean command will delete all created object files
(6.b) MAKE veryclean cleans up all modules in src as well
 
This code is free, and open source, so there is  
NO SUPPORT OR WARRANTY  WHATSOEVER 

hint: MAKE CLEAN often solves some mysterious problems...

hint: if -100 is an output in the .hel file this means, that the number would have been NAN. 

hint: some features are not available for all integrators, i.e. auto planetoid generation input for Yarkovsky integrators.
 
hint: minimum merging radius will let particles merge for sure, regardless of their Hill's radii, maximum merging radius will limit the merging to the distance specified
