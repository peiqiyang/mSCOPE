# mSCOPE

mSCOPE is an integrated multi-layer model of vegetation reflectance, photosynthesis, fluorescence, temperature and energy balance.


## Brief introduction

This version of the mSCOPE model (mSCOPE_v1_beta) is based upon the SCOPE model (v1.61). It simulates the light interaction and energy balance of vertically heterogeneous canopies. 
mSCOPE keeps the same model structure and output of SCOPE, but uses a different solution for the radiative transfer of incident and emitted radiation in vegetation canopies.  
Questions related to mSCOPE, please contact p.yang@utwente.nl; or peiqiyangweb@gmail.com (Peiqi Yang). 

### Prerequisites

The model is written in Matlab codes. It was developed and tested using Matlab 2017. Using older versions may lead to some problems.
```
For example, the matrix multiplication used in RTMo_m.m and RTMf_m.m. 
To do the multiplication of two matrices with the dimension of 60x2162 (A) and 60x1 (B), in the older versions, one may need to use a loop 
for j=1:60
C(j,:) = A(j,:)*B(j)
end
While, in Matlab 2017 or later versions, it is simply computed as A*B;
```
The input data are structured in a excel spreadsheet. The users are expected to have the Microsoft installed. 
If there is not, please contact p.yang@utwente.nl. This issue can be easily fixed by giving alternative  inpute options. 

## What will you receive when you download the model
In the main folder, you will find 
- manuals
			- SCOPE manuals
			- SCOPE and mSCOPE presetations
			
- inputdata
			- input
				-  soil_spectrum (reflectance spectra of several typic soils)
					soilnew.txt. 
					soil_field.txt.
					Note: the users can also replace with their own soil spectra in the model. 
					
				-  radiationdata
					Esun_.dat. 					Incoming direct solar light spectrum
					Esky_.dat. 					Incomging direct diffuse light spectrum
					FLEX-S3_std.atm. 		T18 system 
					T18 system: Verhoef, W., van der Tol, C., & Middleton, E. M. (2018). Hyperspectral radiative transfer modeling to explore the combined retrieval of biophysical parameters and canopy fluorescence from FLEX–Sentinel-3 tandem mission multi-sensor data. Remote sensing of environment, 204, 942-963.
			
				-  fluspect_parameters  (absorption coefficients of constituents of a leaf, and fluroescence emission basic spectra)
					optipar_fluspect.txt
					Optipar_fluspect_2014.txt
			
				-  directional 				(for directional simulation, the angles in the files will be simulated)
					brdf_angles.dat. 	oversampling in the hot spot position
					brdf_angles_no_oversampling.dat
						
				-  dataset for_verification
		
			- measured  					(you will not need these files for the use of mSCOPE, because it is for time series simulation which is not available yet in mSCOPE)
				-	dataset for_verification
- output 
			Note: the output is saved in the directory 
			
- mSCOPE_code_v1
			-	Input_data.xlsx. 		You change the input parameters in this file
			-  mSCOPE.m. 				The main function 
			
			- RTMs. 							leaf and canopy radiative transfer models
				-	fluspect_b.  			simulating leaf reflectance, transmittance and fluorescence emission matrices
				-	RTMo_m.    			canopy radiative transfer in the solar domain. 
				- 	RTMf_m. 				canopy radiative transfer model for fluorescence 
				- 	RTMt. 						canopy radiative transfer for emitted thermal radiation	
				
			- Supporting.  
				-	Brightness_T.  		converting radiant emittance(energy per time per area) to temperature for backbody by inverting Stefan–Boltzmann law	
				-  calczenithangle. 	computing solar position based on the location and time
				- 	e2phot.					calculating the number of moles of photons corresponding to E Joules of energy of wavelength lambda
				-  ephoton. 			 	calculating the energy content (J) of 1 photon of wavelength lambda (m)
				
			- IO. 								reading inputdata, exporting and ploting (optional) simulation results
			
			- Fluxes 							Computing photosynthesis, latent and sensible heat, leaf temperature 


##  Summary of the main changes in mSCOPE
a. input_mSCOPE.m 
It reads the vertical profiles of leaf optical properties (e.g. Cab, Cw) from input_data.xlsx in spratsheet 'mSCOPE'
The input_mSCOPE is called in the main function mSCOPE.m L61 before executing fluspect_mSCOPE and canopy RTMs. 

b. fluspect_mSCOPE.m 
It runs fluspect_b for different layers to obtain leaf reflectance, transmittance, Mb and Mf
fluspect_mSCOPE is called in the main function mSCOPE.m L249

c. RTMo_m.m 
It is a replacement of the RTMo.m in SCOPE. Many changes have been made here. 
RTMo_m.m is called in the main function mSCOPE.m L279    

d. RTMf_m.m
It is a replacement of the RTMf.m in SCOPE. Many changes have been made here. 
RTMf_m is called in the main function mSCOPE.m L279   

## References
The origial SCOPE paper 

Van der Tol, C., Verhoef, W., Timmermans, J., Verhoef, A., & Su, Z. (2009). An integrated model of soil-canopy spectral radiances, photosynthesis, fluorescence, temperature and energy balance. Biogeosciences, 6(12), 3109-3129. 

The mSCOPE paper

Yang, P., Verhoef, W., & Van Der Tol, C. (2017). The mSCOPE model: A simple adaptation to the SCOPE model to describe reflectance, fluorescence and photosynthesis of vertically heterogeneous canopies. Remote sensing of environment, 201, 1-11.

## Authors

Peiqi Yang (p.yang@utwente.nl; peiqiyangweb@gmail.com)

Wout Verhoef  (w.verhoef@utwente.nl)

Christiaan van der Tol (c.vandertol@utwente.nl)

## License
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.

