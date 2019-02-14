
function [rad,profiles] = RTMf_m(spectral,rad,soil,leafopt,canopy,gap,angles,profiles)

% function 'RTMf' calculates the spectrum of fluorescent radiance in the
% observer's direction and also the TOC spectral hemispherical upward Fs flux
%
% Authors:  Wout Verhoef and Christiaan van der Tol (tol@itc.nl)
% Date:     12 Dec 2007
% Update:   26 Aug 2008 CvdT        small correction to matrices
%           07 Nov 2008 CvdT        changed layout
% Update:   19 Mar 2009 CvdT        major corrections: lines 95-96,
%                                   101-107, and 119-120.
% Update:    7 Apr 2009 WV & CvdT   major correction: lines 89-90, azimuth
%                                   dependence was not there in previous verions (implicit assumption of
%                                   azimuth(solar-viewing) = 0). This has been corrected
% Update:   May-June 2012 WV & CvdT Add calculation of hemispherical Fs
%                                   fluxes 
% Update:   Jan-Feb 2013 WV         Inputs and outputs via structures for
%                                   SCOPE Version 1.40
% Update:   Aug-Oct 2016 PY         Re-write the calculation of emitted SIF
%                                   of each layer. It doesnt use loop at
%                                   all. with the function bsxfun, the
%                                   calculation is much faster
% Update:   Oct 2017-Feb 2018 PY    Re-write the RTM of fluorescence                              


% Table of contents of the function:
%   0       preparations
%       0.0     globals
%       0.1     initialisations
%       0.2     geometric quantities
%       0.3     solar irradiance factor and ext. in obs dir for all leaf angle/azumith classes
%   1.0     calculation of fluorescence flux in observation direction
%
% Usage: [rad] = RTMfH(spectral,rad,soil,leafopt,canopy,gap,angles,profiles)
%
% The input and output are structures. These structures are further
% specified in a readme file.
%
% Input:
%   spectral    information about wavelengths and resolutions
%   rad         a large number of radiative fluxes: spectrally distributed
%               and integrated, and canopy radiative transfer coefficients.
%   soil        soil properties
%   leafopt     leaf optical properties
%   canopy      canopy properties (such as LAI and height)
%   gap         probabilities of direct light penetration and viewing
%   angles      viewing and observation angles
%   profiles      vertical profiles of fluxes 
%
% Output:
%   rad         a large number of radiative fluxes: spectrally distributed
%               and integrated, and canopy radiative transfer coefficients.
%               Here, fluorescence fluxes are added
%% 0.0 globals

global constants

%% 0.1 initialisations
wlS          = spectral.wlS';       % SCOPE wavelengths, make column vectors
wlF          = spectral.wlF';       % Fluorescence wavelengths
wlE          = spectral.wlE';       % Excitation wavelengths
[dummy,iwlfi]    = intersect(wlS,wlE); %#ok<ASGLU>
[dummy,iwlfo]    = intersect(wlS,wlF); %#ok<ASGLU>
nf           = length(iwlfo);
ne           = length(iwlfi);
nl           = canopy.nlayers;
LAI          = canopy.LAI;
litab        = canopy.litab;
lazitab      = canopy.lazitab;
lidf         = canopy.lidf;
nlazi        = length(lazitab);         % azumith angle
nlinc        = length(litab);           % inclination
nlori        = nlinc * nlazi;           % total number of leaf orientations
layers       = 1:nl;

Ps           = gap.Ps;
Po           = gap.Po;
Pso          = gap.Pso;

Qso         = (Pso(layers) + Pso(layers+1))/2;
Qs          = (Ps(layers)  + Ps(layers+1))/2;
Qo          = (Po(layers)  + Po(layers+1))/2;

Qsho        = Qo - Qso;                % hot spot
etah        = zeros(nl,1);            % shade leaf
etau        = zeros(nlinc,nlazi,nl);    % sunlit leaf

[Mb,Mf]             = deal(zeros(nf,ne,nl));
[MpluEmin   ,...
    MpluEplu   , ...
    MminEmin   , ...
    MminEplu]       = deal(zeros(nf,nl));

% emitted fluorscence initiation
piLem_              = zeros(nf,nl,2);

% for speed-up the calculation only uses wavelength i and wavelength o part of the spectrum
Esunf_             = rad.Esun_(iwlfi);
Eminf_             = rad.Emin_(:,iwlfi)';          % transpose into [nwlfo,nl] matrix
Epluf_             = rad.Eplu_(:,iwlfi)';
iLAI               = LAI/nl;                       % LAI of a layer        [1]

Xdd         = rad.Xdd(:,iwlfo);
rho_dd      = rad.rho_dd(:,iwlfo);
R_dd        = rad.R_dd(:,iwlfo);
tau_dd      = rad.tau_dd(:,iwlfo);
vb          = rad.vb(:,iwlfo);
vf          = rad.vf(:,iwlfo);
%% 0.2 geometric quantities
MbI                 = leafopt.MbI;
MbII                = leafopt.MbII;
MfI                 = leafopt.MfI;
MfII                = leafopt.MfII;

% geometric factors
deg2rad             = constants.deg2rad;
tto                 = angles.tto;
tts                 = angles.tts;
psi                 = angles.psi;
rs                  = soil.refl(iwlfo,:);           % [nwlfo]     soil reflectance
cos_tto             = cos(tto*deg2rad);             % cos observation zenith angle
sin_tto             = sin(tto*deg2rad);             % sin observation zenith angle

cos_tts             = cos(tts*deg2rad);             % cos solar angle
sin_tts             = sin(tts*deg2rad);             % sin solar angle

cos_ttli            = cos(litab*deg2rad);           % cos leaf inclinaation angles
sin_ttli            = sin(litab*deg2rad);           % sin leaf inclinaation angles
cos_phils           = cos(lazitab*deg2rad);         % cos leaf azimuth angles rel. to sun azi
cos_philo           = cos((lazitab-psi)*deg2rad);   % cos leaf azimuth angles rel. to viewing azi

%% 0.3 geometric factors for all leaf angle/azumith classes
cds                 = cos_ttli*cos_tts*ones(1,36) + sin_ttli*sin_tts*cos_phils;  % [nli,nlazi]
cdo                 = cos_ttli*cos_tto*ones(1,36) + sin_ttli*sin_tto*cos_philo;  % [nli,nlazi]
fs                  = cds/cos_tts;                                               % [nli,nlazi]
absfs               = abs(fs);                                                   % [nli,nlazi]
fo                  = cdo/cos_tto;                                               % [nli,nlazi]
absfo               = abs(fo);                                                   % [nli,nlazi]
fsfo                = fs.*fo;                                                    % [nli,nlazi]
absfsfo             = abs(fsfo);                                                 % [nli,nlazi]
foctl               = fo.*(cos_ttli*ones(1,36));                                 % [nli,nlazi]
fsctl               = fs.*(cos_ttli*ones(1,36));                                 % [nli,nlazi]
ctl2                = cos_ttli.^2*ones(1,36);                                    % [nli,nlazi]

% reshape all the variables
absfs               = reshape(absfs,nlori,1);                                    % [nlori,1]
absfo               = reshape(absfo,nlori,1);                                    % [nlori,1]
fsfo                = reshape(fsfo,nlori,1);                                     % [nlori,1]
absfsfo             = reshape(absfsfo,nlori,1);                                  % [nlori,1]
foctl               = reshape(foctl,nlori,1);                                    % [nlori,1]
fsctl               = reshape(fsctl,nlori,1);                                    % [nlori,1]
ctl2                = reshape(ctl2,nlori,1);                                     % [nlori,1]

%% 1.0 calculation of fluorescence flux in observation direction

% fluorescence efficiencies from ebal, after default fqe has been applied

etahi = profiles.etah;
etaui = profiles.etau;

% fluorescence matrices and efficiencies for PSI and PSII
for PS = 2:-1:1
    [U,Fmin_,Fplu_] =deal(zeros(61,211));
    switch PS 
        case 1, Mb = MbI;  Mf = MfI;    etah(:) = 1;          etau(:) = 1;
        case 2, Mb = MbII; Mf = MfII;   etah(:) = etahi(:);   etau(:) = etaui(:);
    end
    
    Mplu = 0.5*(Mb+Mf);    % [nwlfo,nwlfi]
    Mmin = 0.5*(Mb-Mf);    % [nwlfo,nwlfi]
    
    % in-products: we convert incoming radiation to a fluorescence spectrum using the matrices.
    % resolution assumed is 1 nm
    % [nf,nl] =[nf,ne,nl] * [ne, nl]
    % I have to use a loop here. 
    
    for j = 1:nl
        MpluEmin(:,j)   = Mplu(:,:,j)* Eminf_(:,j);          % [nwlfo,nl+1]
        MpluEplu(:,j)   = Mplu(:,:,j)* Epluf_(:,j);          % [nwlfo,nl+1]
        MminEmin(:,j)   = Mmin(:,:,j)* Eminf_(:,j);          % [nwlfo,nl+1]
        MminEplu(:,j)   = Mmin(:,:,j)* Epluf_(:,j);          % [nwlfo,nl+1]
        MpluEsun(:,j)   = Mplu(:,:,j)* Esunf_;               % integration by inproduct
        MminEsun(:,j)   = Mmin(:,:,j)* Esunf_;               % integration by inproduct
    end
    %% I am trying to not use any loop
    laz= 1/36;
    etau_lidf = bsxfun(@times,reshape(etau,nlori,nl),repmat(lidf*laz,36,1));     %[nlori,nl]
    etah_lidf = bsxfun(@times,repmat(etah,1,nlori)',repmat(lidf*laz,36,1));
    
    wfEs      =  bsxfun(@times,sum(bsxfun(@times,etau_lidf,absfsfo)),MpluEsun) +...
        bsxfun(@times,sum(bsxfun(@times,etau_lidf,fsfo)),MminEsun);
    
    sfEs     =  bsxfun(@times,sum(bsxfun(@times,etau_lidf,absfs)),MpluEsun) -...
        bsxfun(@times,sum(bsxfun(@times,etau_lidf,fsctl)),MminEsun);
    
    sbEs     =  bsxfun(@times,sum(bsxfun(@times,etau_lidf,absfs)),MpluEsun) +...
        bsxfun(@times,sum(bsxfun(@times,etau_lidf,fsctl)),MminEsun);
    
    vfEplu_h  = bsxfun(@times,sum(bsxfun(@times,etah_lidf,absfo)),MpluEplu) -...
        bsxfun(@times,sum(bsxfun(@times,etah_lidf,foctl)),MminEplu);
    
   vfEplu_u  = bsxfun(@times,sum(bsxfun(@times,etau_lidf,absfo)),MpluEplu) -...
    bsxfun(@times,sum(bsxfun(@times,etau_lidf,foctl)),MminEplu);
    
    vbEmin_h  = bsxfun(@times,sum(bsxfun(@times,etah_lidf,absfo)),MpluEmin) +...
        bsxfun(@times,sum(bsxfun(@times,etah_lidf,foctl)),MminEmin);
    
    vbEmin_u  = bsxfun(@times,sum(bsxfun(@times,etau_lidf,absfo)),MpluEmin) +...
        bsxfun(@times,sum(bsxfun(@times,etau_lidf,foctl)),MminEmin);
    
    sigfEmin_h  = bsxfun(@times,sum(etah_lidf),MpluEmin) -...
        bsxfun(@times,sum(bsxfun(@times,etah_lidf,ctl2)),MminEmin);
    
    sigfEmin_u  = bsxfun(@times,sum(etau_lidf),MpluEmin) -...
        bsxfun(@times,sum(bsxfun(@times,etau_lidf,ctl2)),MminEmin);
    
    sigbEmin_h  = bsxfun(@times,sum(etah_lidf),MpluEmin) +...
        bsxfun(@times,sum(bsxfun(@times,etah_lidf,ctl2)),MminEmin);
    
    sigbEmin_u  = bsxfun(@times,sum(etau_lidf),MpluEmin) +...
        bsxfun(@times,sum(bsxfun(@times,etau_lidf,ctl2)),MminEmin);
    
    sigfEplu_h  = bsxfun(@times,sum(etah_lidf),MpluEplu) -...
        bsxfun(@times,sum(bsxfun(@times,etah_lidf,ctl2)),MminEplu);
    
    sigfEplu_u  = bsxfun(@times,sum(etau_lidf),MpluEplu) -...
     bsxfun(@times,sum(bsxfun(@times,etau_lidf,ctl2)),MminEplu);
    
    sigbEplu_h  = bsxfun(@times,sum(etah_lidf),MpluEplu) +...
        bsxfun(@times,sum(bsxfun(@times,etah_lidf,ctl2)),MminEplu);
    sigbEplu_u  = bsxfun(@times,sum(etau_lidf),MpluEplu) +...
        bsxfun(@times,sum(bsxfun(@times,etau_lidf,ctl2)),MminEplu);
%   Emitted fluorescence  
    piLs        =   wfEs+vfEplu_u+vbEmin_u;         % sunlit for each layer
    piLd        =   vbEmin_h+vfEplu_h;              % shade leaf for each layer
    Fsmin       =   sfEs+sigfEmin_u+sigbEplu_u;     % Eq. 29a for sunlit leaf
    Fsplu       =   sbEs+sigbEmin_u+sigfEplu_u;     % Eq. 29b for sunlit leaf
    Fdmin       =   sigfEmin_h+sigbEplu_h;          % Eq. 29a for shade leaf
    Fdplu       =   sigbEmin_h+sigfEplu_h;          % Eq. 29b for shade leaf    
    Femmin      =   iLAI*bsxfun(@times,Qs', Fsmin) +iLAI* bsxfun(@times,(1-Qs)',Fdmin);
    Femplu      =   iLAI*bsxfun(@times,Qs', Fsplu) +iLAI*bsxfun(@times,(1-Qs)',Fdplu);
    Femo        =   iLAI*bsxfun(@times,Qs',piLs) + iLAI*bsxfun(@times,(1-Qs)',piLd);
   

    piLem_(:,:,PS)  =   Femo;
%      for j=1:60      % from bottom to top   it is different in SCOPE where 1 is top, 60 is bottom
%         Y(j,:)  =(rho_dd(j,:).*U(j,:)+Femmin(:,j)')./(1-rho_dd(j,:).*R_dd(j,:));
%         U(j+1,:) =tau_dd(j,:).*(R_dd(j,:).*Y(j,:)+U(j,:))+Femplu(:,j)';
%     end 
    for j=60:-1:1      % from bottom to top   
        Y(j,:)  =(rho_dd(j,:).*U(j+1,:)+Femmin(:,j)')./(1-rho_dd(j,:).*R_dd(j+1,:));
        U(j,:) =tau_dd(j,:).*(R_dd(j+1,:).*Y(j,:)+U(j+1,:))+Femplu(:,j)';  
    end 
    
    for j=1:60          % from top to bottom
        Fmin_(j+1,:)  = Xdd(j,:).*Fmin_(j,:)+Y(j,:);
        Fplu_(j+1,:)  = R_dd(j+1,:).*Fmin_(j+1,:)+U(j+1,:);
    end 
        piLo1(:,PS)     = iLAI*Pso(1:nl)'*piLs';
        piLo2(:,PS)     = iLAI*(Po(1:nl)-Pso(1:nl))'*piLd';
        piLo3(:,PS)     = iLAI*(Po(1:nl)'*(vb.*Fmin_(layers,:) + vf.*Fplu_(layers,:))); 
        piLo4(:,PS)     = rs .* Fmin_(61,:)' * Po(61);     
end 
piLtot      = piLo1 + piLo2 + piLo3 + piLo4;
LoF_        = piLtot/pi;  
Fhem_       = Fplu_; 
rad.LoF_    = LoF_(:,1)  + LoF_(:,2); 
rad.LoF1_   = LoF_(:,1);
rad.LoF2_   = LoF_(:,2);
rad.Fhem_   = Fhem_(:,1) + Fhem_(:,2);

% profiles.fluorescence   = Fiprofile(:,1) + Fiprofile(:,2);
profiles.fluorescence   =   0;
profiles.fluorescenceEm = zeros(nl,1);
for i = 1:nl
    profiles.fluorescenceEm(i) = 0.001 * Sint(sum(piLem_(:,i,:),3)',spectral.wlF)';
end
rad.Eoutf = 0.001 * Sint(sum(Fhem_,1)',spectral.wlF);


