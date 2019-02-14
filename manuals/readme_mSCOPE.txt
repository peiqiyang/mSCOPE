mSCOPE: An integrated multi-layer model of vegetation reflectance, photosynthesis, fluorescence, themperature and energy balance
 
1.0 Introduction 
This version of the mSCOPE model (mSCOPE_v1_beta) is based upon the SCOPE model (v1.61). 
It simulates the light interaction and energy balance of vertically heterogeneous canopies. 
mSCOPE keeps the same structure, output of SCOPE, but it uses adifferent solution for the radiative transfer in vegetation canopies.  
For the use of mSCOPE, please contact p.yang@utwente.nl (Peiqi Yang). 

The spreadsheet input_data.xlsx is the main input file.
mSCOPE is the main function.

2.0. Structure of the model
The 

3.0 Summary of the main changes in mSCOPE and the usage.
a. input_mSCOPE.m 
It reads the vertical profiles of leaf optical properties (e.g. Cab, Cw) from input_data.xlsx in spratsheet 'mSCOPE'
The input_mSCOPE is called in the main function mSCOPE.m L61 before excuting fluspect_mSCOPE and canopy RTMs. 
b. fluspect_mSCOPE.m 
It runs fluspect_b for different layers to obtain leaf reflectance, transmittance, Mb and Mf
fluspect_mSCOPE is called in the main function mSCOPE.m L249
c. RTMo_m.m 
It is a replacement of the RTMo.m in SCOPE. Many changes have been made here. 
RTMo_m.m is called in the main function mSCOPE.m L279    
d. RTMf_m.m
It is a replacement of the RTMf.m in SCOPE. Many changes have been made here. 
RTMf_m is called in the main function mSCOPE.m L279   


4.0 References 
The origial SCOPE paper 
Van der Tol, C., Verhoef, W., Timmermans, J., Verhoef, A., & Su, Z. (2009). An integrated model of soil-canopy spectral radiances, photosynthesis, fluorescence, temperature and energy balance. Biogeosciences, 6(12), 3109-3129. 
The mSCOPE paper
Yang, P., Verhoef, W., & Van Der Tol, C. (2017). The mSCOPE model: A simple adaptation to the SCOPE model to describe reflectance, fluorescence and photosynthesis of vertically heterogeneous canopies. Remote sensing of environment, 201, 1-11.

