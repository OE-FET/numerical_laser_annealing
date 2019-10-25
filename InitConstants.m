function [  ] = InitConstants( )
%INITCONSTANTS Summary of this function goes here
%   Detailed explanation goes here


polar = 0;

%cured IZO
dieIn4Zn1O = 4+2.1*1i;
%non cured
dieIn4Zn1O = 1.9+0.56*1i;
dieIn2Zn1O = 1.85+1.1i;

%Real gold 4 (https://www.clear.rice.edu/elec603/spring2008/Selecting_a_Paper_files/Phys%20%20Rev%20%20B%201972%20Johnson.pdf)

dieGold = -1 + 1i* 4.451;

dieITO = 3.75 + 2*1i;
dieSiN3 = 5.15 + 0.045*1i;
dieGlass = 4.7 + 0.0001*1i;
dieTiO = 5.76 + 3.5*1i;
dieAlN = 4.812+ 0.005*1i;
dieSi = -10.5 + 12.4*1i; %http://www.pveducation.org/pvcdrom/materials/optical-properties-of-silicon
dieSiO2 = 2.2801;
diePEN = 2.6;
diePoly  = 2.6;
%dieGlass = 3.3856 + 0.0001*1i;
dieAl2O3 = 3.3856;





%specific heat Cp J/kg/k
cpAir = 0.718*1000;
cpGlass = 730;
cpPEN = 1200;
cpPoly = 1200;

cpITO = 360;
cpIZO = 480;
cpGold = 0.13*1000*1.15;
cpSi = 0.7*1000;

cpSiN3 = 	(673+1100)/2;
cpTiO = (683+697)/2;
cpAlN = 730; 
cpAl2O3 = 880;

%note: aluminum oxide has a high heat capacity of 880 J/kg/k (for ZnO Heat
%capacity of ZnO with cubic structure at high temperatures)
%(http://www.sciencedirect.com/science/article/pii/S0038109806007496)
%(https://www.filmetrics.com/refractive-index-database/Si3N4/SiN-SiON-Silicon-Nitride)
%'All-Optical Control of a Single Plasmonic Nanoantenna?ITO Hybrid' -CpITO




%thermal conductivity - k- W/m/k
kAir = 0.25;
kITO = 3.95;
kGlass = 2;
kSiO2ThinFilm = 1.3;
kIZO = kITO;
kGold = 314;
kPEN = 0.26; % typical values between 0.14 to 0.24 for similar  http://www.professionalplastics.com/professionalplastics/ThermalPropertiesofPlasticMaterials.pdf
kPoly = 0.38;
kSiN3 = (10+43)/2;
kTiO = (4.8+11.8)/2;
kAlN = (140)/2;
kSi = 130;
kAl2O3 = 35;

%Density -p- kg/m3
pAir = 1;
pITO = 7140;
pAl2O3 = 3.95*1000;
pGlass = 2500;
pPEN = 1360;
pPoly = 1450;
pGold = 19300;
pSiN3 = 3.44*1000;
pTiO = 4.23*1000;
pIZO = pITO;
pAlN = 3.26*1000;
pSi= 2.3290*1000;

save('SimMaterialConstants')



end

