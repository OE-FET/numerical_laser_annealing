function [] = Run_diffusion_eqn(fileName)
% Runs the finite element thermal simulations given a symbollicaly solved
% EM feild. 

% loads EM and stack
load(strcat(fileName,'_init_wEM.mat'));

delete('tempSimulation.mat');
save('tempSimulation');
clearvars -except fileName Fluence;
for(indexF = 1:length(Fluence))
    indexF
    % simulates thermal difussion equation for a given fluence
    Simulate(indexF ,'tempSimulation');
    
    
end
load('tempSimulation');
pause(10);

delete(strcat(fileName,'_outPut.mat'));
pause(10);
save(strcat(fileName,'_outPut.mat'));
disp('saved');
end