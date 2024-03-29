function [  ] = Run_EM( fileName )
%RUN_EM Runs a function to solve all relavent EM parameters befor solving
%diffusion equation numerically

% loads appropriate initialization file
load(strcat(fileName,'_init.mat'));

% thickness of each layer - thick(layer)
syms t z;
%Determine crossAR
intedAR(1,1) =t-t;
% for each Fluence
for(indexFluence = 1:length(Fluence))
    disp(['Started to cal ',num2str(Fluence(indexFluence)),' J/m^2'])
    for(indexGeo = 1:numGeo)
        
        % Solved for amplitudes (AR) of EM waves propagating through stack. 
        [AR(indexFluence,indexGeo,1:NumLayer(indexGeo)),tempE,winL] = genAR(Eps(indexGeo,1:NumLayer(indexGeo)), NumLayer(indexGeo), Thick(indexGeo,1:NumLayer(indexGeo)), Mew(indexGeo,1:NumLayer(indexGeo)),Fluence(indexFluence),0);
        E(indexFluence,indexGeo,1:NumLayer(indexGeo)) = tempE(1,1:end);
        
        
    end
end
save(strcat(fileName,'_init_wEM.mat'));
end

