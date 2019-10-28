function [AR,E,winL] = genAR(Eps, numlayers, thick, mew,fluence,polar)
% Solved for absorption rate (AR) of EM waves propagating through stack. 

display('Calculating Real and imagine compontents of the refractive index')
for(index = 1:numlayers)
    n_Er = sqrt((abs(Eps(index))+real(Eps(index)))/2)*sqrt(mew(index));
    k_Er = sqrt((abs(Eps(index)) - real(Eps(index)))/2)*sqrt(mew(index));
    ComplexIndex(index) = n_Er + 1i*k_Er;
end


display('Generating wave equations for the laser pulse')
[AxVal,AyVal, fVal, KVal,numFFTComp,winL] = genPulse(lambda, PulseWidth,window,polar, fluence);

display('Solving for specific reflection coefficents using incident laser')
solA = SolveA(numlayers, thick);

display('Generating Ex of waveforms')
% generate E Field
tempEx = genE(solA,numlayers,numFFTComp, AxVal, fVal, KVal, ComplexIndex, Eps,0);
Ex = CheckE(tempEx,numlayers,numFFTComp, AxVal, fVal, KVal, ComplexIndex, Eps,thick);

display('Generating Ey of waveforms')
tempEy = genE(solA,numlayers,numFFTComp, AyVal, fVal, KVal, ComplexIndex, Eps,0);
Ey = CheckE(tempEy,numlayers,numFFTComp, AyVal, fVal, KVal, ComplexIndex, Eps,thick);

E = [Ex;Ey];

% Displahy E Field
%plotE(Ex,thick, numlayers);
display('Generating electromagnetic wave parameters D,B and H')
[D,B,H] = GenDBH(E,Eps,mew,numlayers);
%[~,~,AR,~, ~] = poynting(E,D,B,H,numlayers);

display('Generating poynting vector components')
[U,Sz,AR,dUt, dSz] = poynting(E,D,B,H,numlayers)
plotDouble(E(1,1:end),AR./2.3./10^4,thick, numlayers,winL);
%plots Electrogmatic wave interaction
end


function [solA  ] = SolveA(numlayers, thick)
%% solved for the Amplitudes of reflected EM waves 

% reduces precision. 
% see ppt for parameter info. 
digits(16);

solA = sym('solA',[1 numlayers*2-1]);
A = sym('a',[1 numlayers*2-1]);
K = sym('k',[1 numlayers]);
syms z;
syms z1 z2;


for index = 1: numlayers-1
    
    Z = thick(index);
    if(index == numlayers - 1)
        Eqns(index) = eq( A(2*(index-1)+1)*exp(i*K(index)*z1) + A(2*(index-1)+2)*exp(-i*K(index)*z1),A(2*(index-1)+3)*exp(i*K(index+1)*z2));
    else
        Eqns(index) = eq( A(2*(index-1)+1)*exp(i*K(index)*z1) + A(2*(index-1)+2)*exp(-i*K(index)*z1),A(2*(index-1)+3)*exp(i*K(index+1)*z2) + A(2*(index-1)+4)*exp(-i*K(index+1)*z2));
    end
    
    diffEqns(index) = diff(Eqns(index),z1)+diff(Eqns(index),z2);
    diffEqns(index) = subs(diffEqns(index),'z1',Z);
    Eqns(index) = subs(Eqns(index),'z1',Z);
    diffEqns(index) = subs(diffEqns(index),'z2',0);
    Eqns(index) = subs(Eqns(index),'z2',0);
    %set z to sum(thick(1:index))
end
BCEqns = [eq(A(1),1),Eqns,diffEqns];
solA = solve(BCEqns,A);

%solA = extractfield(solA,char(A(1)))
end

function [D,B,H] = GenDBH(E,Eps, mew, numlayers)
% generates displacement field (D), magnetic field (B) and magnetizing
% field (H)

digits(10);
D = [Eps;Eps].*E*8.854187817*10^-12;

%[Ex;Ey]ea
syms z t
for(indexLayers = 1:numlayers)
    tempBx = - diff( E(2,indexLayers),z);
    tempBy =   diff( E(1,indexLayers),z);
    B(1,indexLayers) = -int(tempBx,t);
    B(2,indexLayers) = -int(tempBy,t);
end

H = 1/(4*pi()*10^-7)*[1./mew;1./mew].*B;

end

function [U,Sz,AR,dUt, dSz] = poynting(E,D,B,H,numlayers)
% Geneartes poynting vector coeffecients from EM fields and electromagnetic
% energy density
syms t z;
digits(10);
Utemp = (1/2*(REAL(E).*REAL(D)+REAL(B).*REAL(H)));
U = ((Utemp(1,1:end)+Utemp(2,1:end)));
for(indexLayers = 1:numlayers)
    %U x V  = E X H
    %tempSz = Ux*Vy-Uy*Vx
    Sz(indexLayers) = (((REAL(E(1,indexLayers))*REAL(H(2,indexLayers)) - REAL(E(2,indexLayers))*REAL(H(1,indexLayers)))));
end
dUt = ((diff(U,t)));
dSz = ((diff(Sz,z)));
AR = ((-diff(U,t)-diff(Sz,z)));


end

function [realVal] = REAL(complexVal)
% calculates real component of EM field.
syms x t y z k a;
sym('z','real');
sym('t','real');
realVal = 1/2*(complexVal+conj(complexVal));

end

function [AxVal,AyVal, fVal, kVal,numFFTComp,L] = genPulse(lambda, PulseWidth,window,polar,fluence)
%{ Generate a phenotosecound laser pulse

%}
Kcenter = 2*pi()/(lambda);


%FWHM = 2*sqrt(2*ln(2))*sigma
%sigma = FWHM/(2*sqrt(2*ln(2)))

PulseWidth=299792458*PulseWidth;

sigmaEnergy = PulseWidth/(2*sqrt(2*log(2)));
sigmaP = sigmaEnergy*sqrt(2);
syms z n;
pulse = cos(Kcenter*z)/(sigmaP*sqrt(2*pi()))*exp(-1/2*(z/sigmaP)^2);

%p = matlabFunction(pulse);
%plot(-10*lambda:lambda/10:10*lambda,p(-10*lambda:lambda/10:10*lambda));


ncenter = round(window/lambda)*2;
L = lambda*ncenter/2;
tempCn = sqrt(pi())*sigmaP/2/L*exp(-(n-ncenter)^2/(2*L/(pi()*sigmaP))^2);
Cn = matlabFunction(tempCn);
%plot(ncenter-100:1:ncenter+100,Cn(ncenter-100:1:ncenter+100),'o');


AxVal = Cn(ncenter);
kVal = ncenter*pi()/L;
fVal = 299792458*kVal/2/pi();

offSet = 1;
while( abs(Cn(ncenter+offSet)/Cn(ncenter)) >=1/300)
    AxVal((offSet-1)*2+2) = Cn(ncenter+offSet);
    kVal((offSet-1)*2+2) = (ncenter+offSet)*pi()/L;
    fVal((offSet-1)*2+2)  = 299792458*kVal((offSet-1)*2+2)/2/pi();
    
    AxVal((offSet-1)*2+3) = Cn(ncenter-offSet);
    kVal((offSet-1)*2+3) = (ncenter-offSet)*pi()/L;
    fVal((offSet-1)*2+3)  = 299792458*kVal((offSet-1)*2+3)/2/pi();
    offSet = offSet+1;
end

numFFTComp = length(AxVal);
%normalization

if(polar == 1)
    AyVal = AxVal*1i;
else
    AyVal = zeros(1,numFFTComp);
end


%Normalize Pulse
solA = sym('solA',[1 3]);
A = sym('a',[1 3]);
solA = solve([eq(A(1),1),eq(A(2),0),eq(A(3),0)],A);

tempEx = genE(solA,2,numFFTComp, AxVal, fVal, kVal, [1,1], [1,1],0);
Ex = CheckE(tempEx,2,numFFTComp, AxVal, fVal, kVal, [1,1], [1,1],[1,1]);

tempEy = genE(solA,2,numFFTComp, AyVal, fVal, kVal, [1,1], [1,1],0);
Ey = CheckE(tempEy,2,numFFTComp, AyVal, fVal, kVal, [1,1], [1,1],[1,1]);

E = [Ex;Ey];

% Displahy E Field
%plotE(Ex,thick, numlayers);
tempE(1,1) = E(1,1);
tempE(2,1) = E(2,1);
E = tempE;
[D,B,H] = GenDBH(E,1, 1, 1);
[U,Sz,AR,dUt, dSz] = poynting(E,D,B,H,1);
syms t z;
sym('t','real');
sym('z','real');


%intUatfunctionTot =0;
%num = 100;
%for(index = 1:num)
Utotal = simplify(expand(vpa(U,30)));
%  intagrates total energy absorbed over a distance z
Utotal = quickInt(Utotal,z);
Utz = matlabFunction(Utotal);
Uint = Utz(80*10^-9,L)- Utz(80*10^-9,-L);




%Polarization
if(polar == 1)
    %Uint Typically  4.994884504730976e-11
    AxVal = AxVal/ sqrt(Uint)*sqrt(fluence);
    AyVal = AxVal*1i;
else
    %Unit Typically 2.497442252055493e-11
    AyVal = zeros(1,numFFTComp);
    AxVal = AxVal/ sqrt(Uint)*sqrt(fluence);
end



%intUatfunctionTot = intUatfunctionTot+integral(Uatfunction,-L/2+L/num*(index-1),-L/2+L/num*(index),'RelTol',0,'AbsTol',1e-15);
%end

%%%%%%%%%%% To Check it is normalized
%Polarization
%{
if(polar == 1)
    AxVal = AxVal/ sqrt(2.497442252055493e-11)*sqrt(fluence)/sqrt(2);
    AyVal = AxVal*1i;
else
    AyVal = zeros(1,numFFTComp);
    AxVal = AxVal/ sqrt(2.497442252055493e-11)*sqrt(fluence);
end
%}
%%%%


end



function [E] = genE(solA,numlayers, numFFTComp, AVal, fVal, KVal, ComplexIndex, Eps, infiniteConductivity)
% generate total E fields 

A = sym('a',[1 numlayers*2-1]);
K = sym('k',[1 numlayers]);
digits(10);
syms t z;
E=t-t;
for(layerIndex = 1:numlayers)
    E(layerIndex,1)=t-t;
    for(FFTIndex = 1:numFFTComp)
        
        if(layerIndex == numlayers)
            
            solAR = subs(extractfield(solA,char(A(2*(layerIndex-1)+1))),A(1),1);
            
            for(layerIndex2 = 1:numlayers)
                solAR = subs(solAR,K(layerIndex2),KVal(FFTIndex)*ComplexIndex(layerIndex2));
            end
            
            
            solAR = eval(solAR);
            
            
            if(and(or(isnan(solAR),isinf(solAR)), infiniteConductivity))
                E(layerIndex,FFTIndex) = 0*exp(i*KVal(FFTIndex)*ComplexIndex(layerIndex)*z-1i*2*pi()*fVal(FFTIndex)*t);
            else
                E(layerIndex,FFTIndex) =  AVal(FFTIndex)*solAR*exp(i*KVal(FFTIndex)*ComplexIndex(layerIndex)*z-1i*2*pi()*fVal(FFTIndex)*t);
            end
        else
            solAR = subs(extractfield(solA,char(A(2*(layerIndex-1)+1))),A(1),1);
            solAL = subs(extractfield(solA,char(A(2*(layerIndex-1)+2))),A(1),1);
            for(layerIndex2 = 1:numlayers)
                solAR = subs(solAR,K(layerIndex2),KVal(FFTIndex)*ComplexIndex(layerIndex2));
                solAL = subs(solAL,K(layerIndex2),KVal(FFTIndex)*ComplexIndex(layerIndex2));
            end
            
            solAR = eval(solAR);
            solAL = eval(solAL);
            
            
            
            E(layerIndex,FFTIndex) =  AVal(FFTIndex)*solAR*exp(i*KVal(FFTIndex)*ComplexIndex(layerIndex)*z-1i*2*pi()*fVal(FFTIndex)*t) + AVal(FFTIndex)*solAL*exp(-i*KVal(FFTIndex)*ComplexIndex(layerIndex)*z-1i*2*pi()*fVal(FFTIndex)*t);
        end
        
    end
    
end


end

function[EChecked] = CheckE(E,numlayers,numFFTComp, AVal, fVal, KVal, ComplexIndex, Eps,thick)
%Checks the solved electromagnetic field to make sure it is slef consistent
syms t z;
EChecked = t-t;
for(FFTIndex = 1:numFFTComp)
    
    
    checkFail = sum( E(1:end,FFTIndex));
    if(or(isnan(checkFail),isinf(checkFail)))
        E(1:end,FFTIndex) = 0;
        
        tempNumLayers = numlayers-1;
        stillNaN = 1;
        while(stillNaN)
            
            tempNumLayers;
            tempAVal = AVal(FFTIndex);
            tempKVal = KVal(1:tempNumLayers);
            tempfVal = fVal(FFTIndex);
            tempComplexIndex = ComplexIndex(1:tempNumLayers);
            tempEps = Eps(1:tempNumLayers);
            tempthick = thick(1:tempNumLayers);
            tempE = genE(SolveA(tempNumLayers, tempthick),tempNumLayers,1, tempAVal, tempfVal, tempKVal, tempComplexIndex, tempEps,1);
            checkFail = sum(sum( tempE));
            if(or(isnan(checkFail),isinf(checkFail)))
                stillNaN = 1;
                
                if(tempNumLayers == 2)
                    stillNaN =0;
                    failure = 'solution has failed to converge. tempNumLayers == 1'
                    E(1:end,FFTIndex) = -NaN;
                end
            else
                E(1:tempNumLayers,FFTIndex) = tempE(1:end,1);
                
                for(index2 = tempNumLayers+1:numlayers)
                    E(index2,FFTIndex) = 0;
                end
                stillNaN =0;
            end
            tempNumLayers = tempNumLayers-1;
        end
    end
    
end
for(indexLayer = 1:numlayers)
    
    EChecked(indexLayer) = sum( E(indexLayer,1:end));
    
end


end

function [] = plotSingle(E,thick, numlayers,winL)

% function for plotting a symbol equation (Ex. U) over the stack
% Used for debuggig

fo = 1209677419354830;
%hold on;
syms t z;
resolution = 10000;
resolution = sum(thick(1:numlayers))/resolution;
timeIndex = 0;

%Generate Functions
for(layerIndex =  1:numlayers)
    EFunctions{layerIndex} =  matlabFunction([t,z,real(E(layerIndex))]);
end
ymax = 0;
ymin = 0;
EplotSave = 0;
XplotSave = 0;
for(time = 0:1/fo/100:1/fo*100)
    
    Xplot = 0;
    Eplot = 0;
    EplotTemp =0;
    XplotTemp = 0;
    XplotTemp =0:resolution:sum(thick(1));
    EplotTemp =  EFunctions{1}(time,XplotTemp);
    Eplot = [Eplot,EplotTemp(2+length(XplotTemp):end)];
    Xplot = [Xplot,XplotTemp];
    for(layerIndex = 2:numlayers)
        
        EplotTemp =0;
        XplotTemp = 0;
        XplotTemp =sum(thick(1:layerIndex-1)):resolution:sum(thick(1:layerIndex));
        EplotTemp =  EFunctions{layerIndex}(time,XplotTemp-XplotTemp(1));
        if(length(EplotTemp((length(XplotTemp)+2):end))==1)
            Eplot = [Eplot,zeros(1,length(XplotTemp))];
        else
            Eplot = [Eplot,EplotTemp(2+length(XplotTemp):end)];
        end
        
        Xplot = [Xplot,XplotTemp];
        
    end
    
    timeIndex = timeIndex +1;
    EplotSave(timeIndex,1:length(Eplot)) = (Eplot);
    XplotSave(timeIndex,1:length(Xplot)) = Xplot;
    plot(Xplot,Eplot);
    
    ymax = max(max(Eplot),ymax);
    ymin = min(min(Eplot),ymin);
    hold on;
    for(indexBC = 1:numlayers-1)
        plot([sum(thick(1:indexBC)),sum(thick(1:indexBC))],[-2*ymax,2*ymax]);
    end
    hold off;
    %plot(Xplot,log(abs(Eplot)));
    axis([0 sum(thick(1:numlayers))/10 ymin ymax]);
    %pause(1/60);
    pause(1/60);
    'plot'
    
    
    
    
end

for(playBack = 1:timeIndex)
    
    plot(XplotSave(playBack,1:end),EplotSave(playBack,1:end));
    axis([0 sum(thick(1:numlayers)) -3 3]);
    pause();
    
end
%hold on

%hold off

%time
%pause(1/10);
%Eplot =

end

function [] = plotDouble(A,B,thick, numlayers,winL)

% function for plotting two symbol equations (Ex. U and E) over the stack
% Used for debuggig

fo = 1209677419354830;
%hold on;
syms t z;
resolution = 10000;
resolution = sum(thick(1:numlayers))/resolution;
timeIndex = 0;

%Generate Functions
for(layerIndex =  1:numlayers)
    AFunctions{layerIndex} =  matlabFunction([t,z,real(A(layerIndex))]);
end
for(layerIndex =  1:numlayers)
    BFunctions{layerIndex} =  matlabFunction([t,z,real(B(layerIndex))]);
end
ymax = 0;
Bymax = 0;
ymin = 0;
EplotSave = 0;
XplotSave = 0;
count = 1;

for(time = 0:1/fo/100:1/fo*100)
    Xplot = 0;
    Aplot = 0;
    Bplot = 0;
    AplotTemp =0;
    BplotTemp =0;
    XplotTemp = 0;
    XplotTemp =0:resolution:sum(thick(1));
    AplotTemp =  AFunctions{1}(time,XplotTemp);
    BplotTemp =  BFunctions{1}(time,XplotTemp);
    Aplot = [Aplot,AplotTemp(2+length(XplotTemp):end)];
    Bplot = [Bplot,BplotTemp(2+length(XplotTemp):end)];
    Xplot = [Xplot,XplotTemp];
    for(layerIndex = 2:numlayers)
        
        AplotTemp =0;
        XplotTemp = 0;
        XplotTemp =sum(thick(1:layerIndex-1)):resolution:sum(thick(1:layerIndex));
        AplotTemp =  AFunctions{layerIndex}(time,XplotTemp-XplotTemp(1));
        if(length(AplotTemp((length(XplotTemp)+2):end))==1)
            Aplot = [Aplot,zeros(1,length(XplotTemp))];
        else
            Aplot = [Aplot,AplotTemp(2+length(XplotTemp):end)];
        end
        
        BplotTemp =0;
        BplotTemp =  BFunctions{layerIndex}(time,XplotTemp-XplotTemp(1));
        if(length(BplotTemp((length(XplotTemp)+2):end))==1)
            Bplot = [Bplot,zeros(1,length(XplotTemp))];
        else
            Bplot = [Bplot,BplotTemp(2+length(XplotTemp):end)];
        end
        
        
        
        Xplot = [Xplot,XplotTemp];
        
    end
    
    timeIndex = timeIndex +1;
    EplotSave(timeIndex,1:length(Aplot)) = (Aplot);
    XplotSave(timeIndex,1:length(Xplot)) = Xplot;
    plot(Xplot,Aplot);
    
    ymax = max([max(max(Aplot)),ymax]);
    ymin = min([min(min(Aplot)),ymin]);
    
    Bymax = max(max(max(abs(Bplot))),Bymax);
    hold on;
    plot(Xplot,Bplot./Bymax*ymax);
    for(indexBC = 1:numlayers-1)
        plot([sum(thick(1:indexBC)),sum(thick(1:indexBC))],[-2*ymax,2*ymax]);
    end
    hold off;
    %plot(Xplot,log(abs(Eplot)));
     axis([0 sum(thick(1:numlayers))/10 ymin ymax]);
    pause(1/60);
    'plot'
    Asave(count) = max(max(Aplot));
    Bsave(count) = max(max(Bplot));
    count = count +1;
    
    
end


%{
for(playBack = 1:timeIndex)
    
    plot(XplotSave(playBack,1:end),EplotSave(playBack,1:end));
    axis([0 sum(thick(1:numlayers)) -3 3]);
    pause(1/60);
    
end
   %hold on
    
    %hold off
    
    %time
    %pause(1/10);
    %Eplot =
%}
end

