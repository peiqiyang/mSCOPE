# mSCOPE

mSCOPE is an integrated multi-layer model of vegetation reflectance, photosynthesis, fluorescence, temperature and energy balance.

## Authors

Peiqi Yang (p.yang@utwente.nl; peiqiyangweb@gmail.com)

Wout Verhoef  (w.verhoef@utwente.nl)

Christiaan van der Tol (c.vandertol@utwente.nl)

## Brief introduction

The SCOPE model assumes that vegetation canopies are vertically homogeneous and horizontally infinite, as its radiative transfer routines are based on the classical 1-D SAIL model. However, in reality, canopies generally exhibit large vertical heterogeneity of both biophysical and biochemical properties. The development of mSCOPE is to incorporate the vertical variations of vegetation properties. Therefore, the model can be considered as a 2-D model since it does not consider the horizontal variations. It is noted that mSCOPE works for homogenous canopies as well by just setting all the layers identical or using one layer. The orignial paper on mSCOPE can be found [here](https://www.sciencedirect.com/science/article/pii/S0034425717303954), and the original paper on SCOPE can be found [here](https://www.biogeosciences.net/6/3109/2009/bg-6-3109-2009.html) 

This version of the mSCOPE model (mSCOPE_v1_beta) is based upon the SCOPE model (v1.61). It simulates the light interaction and energy balance of vertically heterogeneous canopies. mSCOPE keeps the same model structure and output of SCOPE, but uses a different solution for the radiative transfer of incident and emitted radiation in vegetation canopies. In mSCOPE, we obtain the reflectance of multi-layer vegetation canopies by using the adding method from bottom to top and then simulate the fluxes profiles using the 'peeling method' from top to bottom. 

Detailed information on SCOPE and mSCOPE can be found from [here](https://scope-model.readthedocs.io/en/latest/) or in the manuals. We expect the users of mSCOPE have already known about the SCOPE model, which can be found from [here](https://github.com/Christiaanvandertol/SCOPE). In the near future, we will update mSCOPE to the latest SCOPE version. 

If you have questions related to mSCOPE, you may contact Peiqi Yang (p.yang@utwente.nl;peiqiyangweb@gmail.com). 

## Prerequisites

The model was tested using **Matlab 2017a**. Using older versions may lead to some problems, for example, the matrix multiplication used in RTMo_m.m and RTMf_m.m. 

```
In the older Matlab versions (e.g. 2013b), multiplying a 60x2162 matrix (A) and a 60x1 matrix (B), one may need to use a loop as follows. 
for j=1:60
C(j,:) = A(j,:)*B(j)
end
In Matlab 2017 or later versions, it is simply computed as A*B, which was used in the model.

```
The users are expected to have the **Microsoft** installed. The input data are structured in an excel spreadsheet. 
If there is not Microsoft, please contact p.yang@utwente.nl. This issue can be easily fixed by giving alternative input options. 

##  Summary of the main changes in mSCOPE

**input_mSCOPE.m** 
It reads the vertical profiles of leaf optical properties (e.g. Cab, Cw) from input_data.xlsx in spratsheet 'mSCOPE'. The input_mSCOPE is called in the main function mSCOPE.m L61 before executing fluspect_mSCOPE and canopy RTMs. 

**fluspect_mSCOPE.m** 
It runs fluspect_b for different layers to obtain leaf reflectance, transmittance, Mb and Mf
fluspect_mSCOPE is called in the main function mSCOPE.m L249

**RTMo_m.m** 
It is a replacement of the RTMo.m in SCOPE. RTMo_m.m is called in the main function mSCOPE.m L279. The main changes are as follows.  

__Thin layer reflectances and transmittances__
 In SCOPE reflectance and transmittances of a vegetation layer are computed in Appendix A in the mSCOPE paper, but in the mSCOPE model, they are directly computed from scattering and extinction coefficients.  
 ``` 
tau_ss = repmat(1-k.*iLAI,nl,1);
tau_dd = 1-a.*iLAI;
tau_sd = sf.*iLAI;
rho_sd = sb.*iLAI;
rho_dd = sigb.*iLAI;  
```
__Effective transmittance and reflectance of a layer bounded beneath a surface__
They are new in the mSCOPE, because the use of the adding method. They are computed from bottom to top of the canopy. 
```
for j=60:-1:1 % from bottom to top. note 61 the background. 1 is the top of canopy.
Xss(j)      = tau_ss(j);
dnorm       = 1-rho_dd(j,:).*R_dd(j+1,:);
Xsd(j,:)    = (tau_sd(j,:)+tau_ss(j).*R_sd(j+1,:).*rho_dd(j,:))./dnorm;
Xdd(j,:)    = tau_dd(j,:)./dnorm;
R_sd(j,:)   = rho_sd(j,:)+tau_dd(j,:).*(R_sd(j+1,:).*Xss(j)+R_dd(j+1,:).*Xsd(j,:));
R_dd(j,:)   = rho_dd(j,:)+tau_dd(j,:).*R_dd(j+1,:).*Xdd(j,:);
end
``` 
__Flux profiles__
They are computed from bottom to top of the canopy. 
``` 
% Eq. 19 in mSCOPE paper
Es_(1,:)       = Esun_;
Emin_(1,:)     = Esky_;
for j=1:60 % from top to bottom
Es_(j+1,:)    =   Xss(j).*Es_(j,:);
Emin_(j+1,:)  =   Xsd(j,:).*Es_(j,:)+Xdd(j,:).*Emin_(j,:);
Eplu_(j,:)    =   R_sd(j,:).*Es_(j,:)+R_dd(j,:).*Emin_(j,:);
end
``` 

**RTMf_m.m** 
It is a replacement of the RTMf.m in SCOPE. RTMf_m is called in the main function mSCOPE.m L279.

## Directory Structure

```
mSCOPE
|-- GNU_General_Public_Licence.txt
|-- Readme.md
|-- inputdata
|   |-- input
|   |   |-- dataset for_verification
|   |   |-- directional
|   |   |   |-- brdf_angles.dat
|   |   |   |-- brdf_angles2.dat
|   |   |   `-- brdf_angles_no_oversampling.dat
|   |   |-- fluspect_parameters
|   |   |   |-- Optipar_fluspect_2014.txt
|   |   |   |-- optipar_fluspect.txt
|   |   |   `-- optipar_fluspect_v152.txt
|   |   |-- radiationdata
|   |   |   |-- Esky_.dat
|   |   |   |-- Esun_.dat
|   |   |   |-- FLEX-S3_std.atm
|   |   |   |-- FLEX-S3_std.tp5
|   |   |   `-- wl.txt
|   |   `-- soil_spectrum
|   |       |-- soil_field.txt
|   |       `-- soilnew.txt
|   `-- measured
|       `-- dataset for_verification
|-- mSCOPE_code_v1
|   |-- Fluxes
|   |   |-- Soil_Inertia0.m
|   |   |-- Soil_Inertia1.m
|   |   |-- biochemical.m
|   |   |-- biochemical_MD12.m
|   |   |-- calculate_vert_profiles.m
|   |   |-- ebal.m              
|   |   |-- heatfluxes.m
|   |   |-- resistances.m
|   |   |-- soil_respiration.m
|   |   `-- zo_and_d.m
|   |-- IO
|   |   |-- MakeDaily.m
|   |   |-- aggreg.m
|   |   |-- assignvarnames.m
|   |   |-- calc_brdf.m
|   |   |-- count.m
|   |   |-- create_output_files.m
|   |   |-- initialize_output_structures.m
|   |   |-- input_mSCOPE.m
|   |   |-- load_timeseries.m
|   |   |-- output_data.m
|   |   |-- output_verification.m
|   |   |-- plot_directional_figure4_function.m
|   |   |-- plots.m
|   |   |-- readStructFromExcel.m
|   |   |-- resizefigure.m
|   |   `-- select_input.m
|   |-- RTMs
|   |   |-- RTMf_m.m
|   |   |-- RTMo_m.m
|   |   |-- RTMt_planck.m
|   |   |-- RTMt_sb.m
|   |   |-- RTMt_sb_new.m
|   |   |-- fluspect_bcar.m
|   |   |-- fluspect_mSCOPE.m
|   |   `-- leafangles.m
|   |-- Supporting
|   |   |-- Brightness_T.m
|   |   |-- Planck.m
|   |   |-- Sint.m
|   |   |-- calc_rssrbs.m
|   |   |-- calczenithangle.m
|   |   |-- define_bands.m
|   |   |-- define_constants.m
|   |   |-- e2phot.m
|   |   |-- ephoton.m
|   |   |-- meanleaf.m
|   |   |-- progressbar.m
|   |   |-- satvap.m
|   |   `-- vangenuchten.m
|   |-- input_data.xlsx
|   `-- mSCOPE.m
|-- manuals
|   |-- SCOPE_presentation.pdf
|   |-- SCOPE_v161 User manual.docx
|   |-- SCOPE_v161 User_manual.pdf
|   |-- mSCOPE_presentation.pdf
|   `-- readme_mSCOPE.txt
|-- output
```
## License
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.
