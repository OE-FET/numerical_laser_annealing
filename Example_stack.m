function Example_stack(fileName)
load('SimMaterialConstants');
%Order of simulations to run
order = [1,2,3,4];
%Number of simulations overall
numGeo = 4;

% Thickness of each layer in stack
Thick(1,1:8) = [180*10^-9,20*10^-9, 7*10^-9,30*10^-9, 26*10^-9,300*10^-9,50*10^-9,10000*10^-9];
Thick(2,1:8) = [180*10^-9,20*10^-9, 7*10^-9,30*10^-9, 26*10^-9,300*10^-9,50*10^-9,10000*10^-9];
Thick(3,1:7) = [180*10^-9,20*10^-9, 7*10^-9,30*10^-9,300*10^-9,50*10^-9,10000*10^-9];
Thick(4,1:6) = [180*10^-9,20*10^-9, 7*10^-9,300*10^-9,50*10^-9,10000*10^-9];

% Finite element resolution of each layer in stack
Res(1,1:8) = [2*10^-9, 2*10^-9, 1*10^-9,2*10^-9,2*10^-9,10*10^-9,10*10^-9,200*10^-9];
Res(2,1:8) = [2*10^-9, 2*10^-9, 1*10^-9,2*10^-9,2*10^-9,10*10^-9,10*10^-9,200*10^-9];
Res(3,1:7) = [2*10^-9, 2*10^-9, 1*10^-9,2*10^-9,10*10^-9,10*10^-9,(10000/63)*10^-9];
Res(4,1:6) = [2*10^-9, 2*10^-9, 1*10^-9,10*10^-9,10*10^-9,(10000/78)*10^-9];

% Name of each layer in stack
nameLayer(1,1:8) = {'air','air','In4ZnO','Au','TiO2','SiO2','PEN','PEN'};
nameLayer(2,1:8) = {'air','air','In4ZnO','ITO','TiO2','SiO2','PEN','PEN'};
nameLayer(3,1:7) = {'air','air','In4ZnO','TiO2','SiO2','PEN','PEN'};
nameLayer(4,1:6) = {'air','air','In4ZnO','SiO2','PEN','PEN'};

% Number of layers in each layer in the stack
NumLayer = [8,8,7,6];

% window length of simulation (m)
winLen = 40;

% Total laser pulse fluence
Fluence = [0.01,0.02,.03,0.04,.05,.06,.07]./0.0001;

% dielectric of each layer in stack
Eps(1,1:8) = [1,1, dieIn4Zn1O, dieGold,dieTiO,dieSiO2,diePEN,diePEN];
Eps(2,1:8) = [1,1, dieIn4Zn1O, dieITO,dieTiO,dieSiO2,diePEN,diePEN];
Eps(3,1:7) = [1,1, dieIn4Zn1O,diePEN,dieTiO,dieSiO2,diePEN];
Eps(4,1:6) = [1,1, dieIn4Zn1O,diePEN,dieSiO2,diePEN];

% magnetic constant of each layer in stack
Mew(1,1:8) = [1,1,1,1,1,1,1,1];
Mew(2,1:8) = [1,1,1,1,1,1,1,1];
Mew(3,1:7) = [1,1,1,1,1,1,1];
Mew(4,1:6) = [1,1,1,1,1,1];


%note: aluminum oxide has a high heat capacity of 880 J/kg/k (for ZnO Heat
%capacity of ZnO with cubic structure at high temperatures)
%(http://www.sciencedirect.com/science/article/pii/S0038109806007496)
%(https://www.filmetrics.com/refractive-index-database/Si3N4/SiN-SiON-Silicon-Nitride)
%'All-Optical Control of a Single Plasmonic Nanoantenna?ITO Hybrid' -CpITO

%specific heat Cp J/kg/k
Cp(1,1:8) =  [cpAir,cpAir, cpIZO,cpGold,cpTiO,cpGlass,cpPEN,cpPEN];
Cp(2,1:8) =  [cpAir,cpAir, cpIZO,cpITO,cpTiO,cpGlass, cpPEN,cpPEN];
Cp(3,1:7) =  [cpAir,cpAir, cpIZO,cpTiO,cpGlass, cpPEN,cpPEN];
Cp(4,1:6) =  [cpAir,cpAir, cpIZO,cpGlass,cpPEN,cpPEN];



%thermal conductivity - k- W/m/k
kcond(1,1:8) =  [kAir,kAir, kIZO,kGold,kTiO,kGlass,kPEN,kPEN];
kcond(2,1:8) =  [kAir,kAir, kIZO,kITO,kTiO,kGlass,kPEN,kPEN];
kcond(3,1:7) =  [kAir,kAir, kIZO,kTiO,kGlass,kPEN,kPEN];
kcond(4,1:6) =  [kAir,kAir, kIZO,kGlass,kPEN,kPEN];

%Density -p- kg/m3
density(1,1:8)= [pAir,pAir,pIZO,pGold,pTiO,pGlass,pPEN,pPEN];
density(2,1:8)= [pAir,pAir,pIZO,pITO,pTiO,pGlass,pPEN,pPEN];
density(3,1:7)= [pAir,pAir,pIZO,pTiO,pGlass,pPEN,pPEN];
density(4,1:6)= [pAir,pAir,pIZO,pGlass,pPEN,pPEN];


%All laser information
%pulse window (m)
window = 20;
% center waveleng (m)
lambda = 248*10^-9;
% width of pulse (s)
PulseWidth = 25*10^-9;

save(strcat(fileName,'_init.mat'))

end
