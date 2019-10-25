function [] = diffusion_eqn( indexF,fileName)

% load relavent data from saved state
load(fileName);

% goes through stored geometry layers
for(indexGeo = 1:numGeo)
    for(indexLayer =2:NumLayer(indexGeo))
        ARfun{indexGeo,indexLayer} = matlabFunction(AR(indexF,indexGeo,indexLayer));
    end
end

% Initializes the temperature distribtuion of 
T = genT(NumLayer,Thick,Res,order);

% sets up heat capacity matrix based on stack
CpMat = genCp(T,NumLayer,Thick,Res,order,Cp);

%Sets up density matrix of stack
DMat = genDen(T,NumLayer,Thick,Res,order,density);

%Sets up thermal conductivty matrix of stack
KMat = genKcond(T,NumLayer,Thick,Res,order,kcond);

% total time of simulation
totalTime = winLen/(299792458);


% number of times the data is presented to user. 
numberOfInterations = 100;

% Step per iteration 
step = 1000000;

% dt cacluatoin
dt = totalTime/numberOfInterations;

% generates location matrix. (simulation is 1D, Vx will be deprecated)
[Vx, Ymat,LayerPos] = genXYVectorMat(NumLayer,Thick,Res,order,numGeo);

% Presented figure #
count = 0;



% to date simulated run length
SimRunLength = 1;

% Displays only the first geometry
displayGeo = 1;

% max temp achieved in  displayGeo simulation
maxTemp = 0;

% iteration index of maixmum temperature
countMax = 1;

% loops through time intervals
for(presentTime = -winLen/299792458/2+dt:dt:winLen/2/299792458+dt)
    tic;
    

    % converts EM field to an energy abosrobed matrix (QMat), Also returns
    % spatial length of each finite element
    [QMat,ResMat] = genQMat2(T,NumLayer,Thick,Res,order,density,presentTime, presentTime+dt,ARfun,numGeo,Eps);
    
    dTcond = 0;
    % change in temperature for unit of energy absorbed
    dTAbsorb = (QMat./ResMat./(10^-18))./DMat./CpMat;
    
    
    tempTemp = T+dTAbsorb/step+dTcond;
    countMe = 1;
    
    
    
    
    for(temptime=presentTime:dt/step:presentTime+dt-dt/step)
        % resets boundary conditions 
        [tempTemp] = setBC(tempTemp);
        %  tempreature temperature in time step 
        tempTemp = tempTemp+dTcond+dTAbsorb/step;
        % heat transfer due to tempearture
        [dQcond] = condT2(tempTemp,KMat,Vx,Ymat,ResMat);
        % tempeature change due to heat transfer
        dTcond = (dQcond.*dt/step)./DMat./CpMat./ResMat./(10^-18);
        
        
        if(mod(countMe,100000)==1)
            for(indexGeo = 1:numGeo )
                maxTemp(indexGeo,countMax) = max(max(tempTemp(indexGeo,1:end)));
            end
            
            
            
            
            figure(1);
            subplot(2,1,1);
            plot(tempTemp(displayGeo,1:end)-293);
            ylabel('T - 293');
            xlabel('Distance');
            figure(1);
            subplot(2,1,2);
            
            
            plot(maxTemp(displayGeo,1:end)-293);
            ylabel('Max Temp');
            xlabel('Count');
            figure(2);
            subplot(2,1,1);
            plot(dTAbsorb(displayGeo,1:end));
            ylabel('dTAbsorb');
            xlabel('Distance');
            figure(2);
            subplot(2,1,2);
            plot(dTcond(displayGeo,1:end).*step);
            ylabel('dTCond');
            xlabel('Distance');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
            %{
            figure(2);
            subplot(4,1,1);
            plot(tempTemp(3,1:end)-293);
            ylabel('T - 293');
            xlabel('Distance');
            subplot(4,1,2);
            plot(dTcond(3,1:end));
            ylabel('dTCond');
            xlabel('Distance');
            subplot(4,1,3);
            plot(dTAbsorb(3,1:end));
            ylabel('dTAbsorb');
            xlabel('Distance');
            subplot(4,1,4);
            plot(maxTemp(2,1:end)-293);
            ylabel('Max Temp');
            xlabel('Count');
            
            %}
            countMax = countMax +1;
            pause(1/60);
        end
        
        countMe = countMe +1;
        
    end
    
    
    
    % resets boundary conditions to room temperature
    T = setBC(tempTemp);
    
    figure(1);
    subplot(2,1,1);
    plot(tempTemp(displayGeo,1:end)-293);
    ylabel('T - 293');
    xlabel('Distance');
    figure(1);
    subplot(2,1,2);
    plot(maxTemp(displayGeo,1:end)-293);
    ylabel('Max Temp');
    xlabel('Count');
    figure(2);
    subplot(2,1,1);
    plot(dTAbsorb(displayGeo,1:end));
    ylabel('dTAbsorb');
    xlabel('Distance');
    figure(2);
    subplot(2,1,2);
    plot(dTcond(displayGeo,1:end).*step);
    ylabel('dTCond');
    xlabel('Distance');
    
    count = count +1;
    
    for(indexGeo = 1:numGeo )
        saveTemp(indexF,indexGeo,count, 1:length(tempTemp(2,1:end))) = (tempTemp(indexGeo,1:end));
        savedTcond(indexF,indexGeo,count, 1:length(dTcond(2,1:end))) = (dTcond(indexGeo,1:end)).*step;
        saveddTAbsorb(indexF,indexGeo,count, 1:length(dTcond(2,1:end))) = (dTAbsorb(indexGeo,1:end));
        savedmaxTemp(indexF,indexGeo,count, 1:length(maxTemp(2,1:end))) = (maxTemp(indexGeo,1:end));
        
    end
    [presentTime] = clock();
    estTime = toc;
    
    hours = round(mod((presentTime(4)*60*60+presentTime(5)*60+presentTime(6)+estTime*(numberOfInterations-SimRunLength))/60/60,12));
    mins = round(mod(mod((presentTime(4)*60*60+presentTime(5)*60+presentTime(6)+estTime*(numberOfInterations-SimRunLength))/60/60,12)*60,60));
    disp(['Run Fluence ',num2str(indexF),' of ',num2str(length(Fluence)),' (',num2str(Fluence(indexF)),' J/m2). Simulation is %',num2str(SimRunLength/(numberOfInterations+1)*100), ' complete. Estimated time of Completion is ', num2str(hours), ':',num2str(mins)]);
    SimRunLength = SimRunLength+1;
    
    
    
end





pause(10);
clear indexF;
delete('tempSimulation.mat');
save('tempSimulation');

end

function [T] = genT(NumLayer,Thick,Res,order)
% Initializes the temperature distribtuion of 

% Default temp (room temp)
T =293;

for(indexX = 1:length(order))
    start = 0;
    for(indexLayer = NumLayer(order(indexX)):-1:1)
        
        for(indexY = start+1:start+round(Thick(order(indexX),indexLayer)/Res(order(indexX),indexLayer)))
            T(indexX,indexY) = 293;
        end
        start = start+round(Thick(order(indexX),indexLayer)/Res(order(indexX),indexLayer));
    end
    
end

T(1:end,1:end) = 293;

end

function [Xvect, Yvect] = genXYVector(NumLayer,Thick,Res,order,numGeo)


Xvect =(1:length(order))*10^(-9);



for(indexX = 1:numGeo)
    start = 0;
    start2 = 0;
    for(indexLayer = NumLayer(indexX):-1:1)
        
        zStep = round(Thick(indexX,indexLayer)/Res(indexX,indexLayer));
        for(indexY = start+1:start+zStep)
            %indexY
            zStop = start2+Res(indexX,indexLayer)*(indexY-start);
            zStart = start2+Res(indexX,indexLayer)*(indexY-1-start);
            Yvect(indexY) = (zStop+zStart)/2;
            
            temp = indexY;
            
        end
        start = indexY;
        start2 = Thick(indexX,indexLayer)+start2;
    end
    
end


end
function [Xvect, Ymat,LayerPos] = genXYVectorMat(NumLayer,Thick,Res,order,numGeo)


Xvect =(1:length(order))*10^(-9);
LayerPos =0;


for(indexX = 1:numGeo)
    start = 0;
    start2 = 0;
    for(indexLayer = NumLayer(indexX):-1:1)
        
        zStep = round(Thick(indexX,indexLayer)/Res(indexX,indexLayer));
        for(indexY = start+1:start+zStep)
            %indexY
            zStop = start2+Res(indexX,indexLayer)*(indexY-start);
            zStart = start2+Res(indexX,indexLayer)*(indexY-1-start);
            Ymat(indexX,indexY) = (zStop+zStart)/2;
            
            temp = indexY;
            
        end
        start = indexY;
        start2 = Thick(indexX,indexLayer)+start2;
        LayerPos(indexX,indexLayer) = indexY;
    end
    
end


end
function [XSpace,YSpace] = genXYSpace(Xvect, Yvect)
XSpace = 0;
ySpace = 0;
for(index = 1:length(Xvect)-1)
    XSpace(index) = (Xvect(index+1)-Xvect(index));
end
XSpace(length(Xvect))= XSpace(length(Xvect)-1);
for(index = 1:length(Yvect)-1)
    YSpace(index) = -Yvect(index+1)+Yvect(index);
end
YSpace(length(Yvect))= YSpace(length(Yvect)-1);
YSpace=-YSpace;
end

function [KMat] = genKcond(T,NumLayer,Thick,Res,order,kcond)



KMat(1:length(T(1:end,1)),1:length(T(1,1:end))) =kcond(1,1);

for(indexX = 1:length(order))
    start = 0;
    for(indexLayer = NumLayer(order(indexX)):-1:1)
        
        for(indexY = start+1:start+round(Thick(order(indexX),indexLayer)/Res(order(indexX),indexLayer)))
            KMat(indexX,indexY) = kcond(order(indexX),indexLayer);
        end
        start = start+round(Thick(order(indexX),indexLayer)/Res(order(indexX),indexLayer));
    end
    
end



end

function [CpMat] = genCp(T,NumLayer,Thick,Res,order,Cp)



CpMat(1:length(T(1:end,1)),1:length(T(1,1:end))) =Cp(1,1);

for(indexX = 1:length(order))
    start = 0;
    for(indexLayer = NumLayer(order(indexX)):-1:1)
        
        for(indexY = start+1:start+round(Thick(order(indexX),indexLayer)/Res(order(indexX),indexLayer)))
            CpMat(indexX,indexY) = Cp(order(indexX),indexLayer);
        end
        start = start+round(Thick(order(indexX),indexLayer)/Res(order(indexX),indexLayer));
    end
    
end



end
function [DMat] = genDen(T,NumLayer,Thick,Res,order,density)



DMat(1:length(T(1:end,1)),1:length(T(1,1:end))) =density(1,1);

for(indexX = 1:length(order))
    start = 0;
    for(indexLayer = NumLayer(order(indexX)):-1:1)
        
        for(indexY = start+1:start+round(Thick(order(indexX),indexLayer)/Res(order(indexX),indexLayer)))
            DMat(indexX,indexY) = density(order(indexX),indexLayer);
        end
        start = start+round(Thick(order(indexX),indexLayer)/Res(order(indexX),indexLayer));
    end
    
end



end



function [QMat,ResMat] = genQMat(T,NumLayer,Thick,Res,order,density,timeStart, timeEnd,Qfun,numGeo)

area = 1.0*10^-9*10^-9;

QMat(1:length(T(1:end,1)),1:length(T(1,1:end))) =0.0;
QMatTemp(1,1:length(T(1,1:end))) = 0.0;

for(indexX = 1:numGeo)
    start = 0;
    for(indexLayer = NumLayer(indexX):-1:2)
        count = 1;
        zStep = round(Thick(indexX,indexLayer)/Res(indexX,indexLayer));
        for(indexY = start+1:start+zStep)
            %indexY
            zStop = Thick(indexX,indexLayer) - (count-1)*Res(indexX,indexLayer);
            zStart = Thick(indexX,indexLayer) - count*Res(indexX,indexLayer);
            Qtemp = (Qfun{indexX,indexLayer}(timeStart,zStart)+Qfun{indexX,indexLayer}(timeEnd,zStop));
            Qtemp = Qtemp - (Qfun{indexX,indexLayer}(timeEnd,zStart)- Qfun{indexX,indexLayer}(timeStart,zStop));
            QMatTemp(indexX,indexY) = real(Qtemp*area);
            count = count +1;
            
        end
        start = start+round(Thick(indexX,indexLayer)/Res(indexX,indexLayer));
    end
    
end
for(indexX = 1:numGeo)
    start = 0;
    for(indexLayer = NumLayer(indexX):-1:1)
        zStep = round(Thick(indexX,indexLayer)/Res(indexX,indexLayer));
        for(indexY = start+1:start+zStep)
            tempRes(indexX,indexY) = Res(indexX,indexLayer);
        end
        start = start+round(Thick(indexX,indexLayer)/Res(indexX,indexLayer));
    end
end


for(indexX = 1:length(order))
    for(indexY = 1:length(QMatTemp(1,1:end)))
        QMat(indexX,indexY) = QMatTemp(order(indexX),indexY);
        ResMat(indexX,indexY) = tempRes(order(indexX),indexY);
    end
end


end

function [xFun,yFun] = GradFun(Mat,Vx,Vy)

xFun = Mat.*0;
yFun = Mat.*0;
for(index = 1:length(Mat(1,1:end)))
    xFun(1:end,index) = gradient(Mat(1:end,index)',Vx)';
end

for(index = 1:length(Mat(1:end,1)))
    xFun(index,1:end) = gradient(Mat(index,1:end),Vy);
end
end


function [T] = setBC(T)

T(1:end,1) = 293;
T(1:end,end) = 293;
%T(1,1:end) = T(2,1:end);
%T(end,1:end) = T(end-1,1:end);
end

function [dQcond] = condT(T,K,Vx,Vy)
%{dQcond = T-T;

[xFun,yFun] = GradFun(T,Vx,Vy);
Kx = 0;
Ky = 0;
for(indexX = 1:length(T(1:end,1))-1)
    for(indexY = 1:length(T(1,1:end))-1)
        
        Kx(indexX,indexY) = (K(indexX+1,indexY)+K(indexX,indexY))/2;
        Ky(indexX,indexY) = (K(indexX,indexY)+K(indexX,indexY+1))/2;
        
    end
end
Kx(end+1, 1:end+1) = K(end,1:end);
Kx(1:end, end) = K(1:end,end);
Ky(end+1, 1:end+1) = K(end,1:end);
Ky(1:end, end) = K(1:end,end);
xFun = Kx.*xFun;
yFun = Ky.*yFun;

[xFun,~] = GradFun(xFun,Vx,Vy);
[~,yFun] = GradFun(yFun,Vx,Vy);

dQcond = xFun+yFun;


dQcond(1,1:end) = 0;
dQcond(end,1:end) = 0;
dQcond(1:end,1) = dQcond(1:end,2);
dQcond(1:end,end) = dQcond(1:end,end-1);
end
function [dQcond] = condT2(T,K,Vx,Ymat,ResMat)
%{dQcond = T-T;


dQcond = T-T;
for(indexX = 1:length(T(1:end,1)))
    for(indexY = 2:length(T(1,1:end))-1)
        
        limitBC = 1;
        
        
        
        if(K(indexX,indexY+1)==K(indexX,indexY))
            dQcond(indexX,indexY) = dQcond(indexX,indexY) + (T(indexX,indexY+1)-T(indexX,indexY))*(K(indexX,indexY))/abs((Ymat(indexX,indexY)-Ymat(indexX,indexY+1)))*10^-9*10^-9;
        else
            dQcond(indexX,indexY) = dQcond(indexX,indexY) + (T(indexX,indexY+1)-T(indexX,indexY))*min(K(indexX,indexY+1),K(indexX,indexY))*limitBC/abs((Ymat(indexX,indexY)-Ymat(indexX,indexY+1)))*10^-9*10^-9;
        end
        
        
        if(K(indexX,indexY-1)==K(indexX,indexY))
            dQcond(indexX,indexY) = dQcond(indexX,indexY) + (T(indexX,indexY-1)-T(indexX,indexY))*((K(indexX,indexY)))/abs((Ymat(indexX,indexY)-Ymat(indexX,indexY-1)))*10^-9*10^-9;
        else
            dQcond(indexX,indexY) = dQcond(indexX,indexY) + (T(indexX,indexY-1)-T(indexX,indexY))*(min(K(indexX,indexY-1),K(indexX,indexY)))*limitBC/abs((Ymat(indexX,indexY)-Ymat(indexX,indexY-1)))*10^-9*10^-9;
        end
    end
end
dQcond(1:end,1) = 0;
dQcond(1:end,end) = 0;
%dQcond(1:end,2) = 0;
%dQcond(1:end,end-1) = 0;
%dQcond(2,1:end) = dQcond(3,1:end);
%dQcond(end-1,1:end) = dQcond(end-2,1:end);
%dQcond(1,1:end) = dQcond(2,1:end);
%dQcond(end,1:end) = dQcond(end-1,1:end);

end

function [QMat,ResMat] = genQMat2(T,NumLayer,Thick,Res,order,density,timeStart, timeEnd,ARfun,numGeo,Eps)

area = 1.0*10^-9*10^-9;

QMat(1:length(T(1:end,1)),1:length(T(1,1:end))) =0.0;
QMatTemp(1,1:length(T(1,1:end))) = 0.0;

for(indexX = 1:numGeo)
    start = 0;
    for(indexLayer = NumLayer(indexX):-1:2)
        count = 1;
        zStep = round(Thick(indexX,indexLayer)/Res(indexX,indexLayer));
        
        for(indexY = start+1:start+zStep)
            %indexY
            zStop = Thick(indexX,indexLayer) - (count-1)*Res(indexX,indexLayer);
            zStart = Thick(indexX,indexLayer) - count*Res(indexX,indexLayer);
            if(not( isreal(  Eps(indexX,indexLayer))))
                QMatTemp(indexX,indexY) = real(integral2(ARfun{indexX,indexLayer},timeStart,timeEnd,zStart,zStop)*area);
            else
                QMatTemp(indexX,indexY) = 0;
            end
            count = count +1;
            
        end
        start = start+round(Thick(indexX,indexLayer)/Res(indexX,indexLayer));
    end
    
end
for(indexX = 1:numGeo)
    start = 0;
    for(indexLayer = NumLayer(indexX):-1:1)
        zStep = round(Thick(indexX,indexLayer)/Res(indexX,indexLayer));
        for(indexY = start+1:start+zStep)
            tempRes(indexX,indexY) = Res(indexX,indexLayer);
        end
        start = start+round(Thick(indexX,indexLayer)/Res(indexX,indexLayer));
    end
end


for(indexX = 1:length(order))
    for(indexY = 1:length(QMatTemp(1,1:end)))
        QMat(indexX,indexY) = QMatTemp(order(indexX),indexY);
        ResMat(indexX,indexY) = tempRes(order(indexX),indexY);
    end
end


end
