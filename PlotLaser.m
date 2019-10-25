function [  ] = PlotLaser( indexFAgain, fileName)
close all



%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
load(fileName);
indexF =indexFAgain;
%numGeo =2;
count = length(savedTcond(1,1,1:end,1));

for(indexGeo = 1:numGeo)
    T = saveTemp(indexF,indexGeo,round(count/2*1.25), 1:end);
    dTcond =  savedTcond(indexF,indexGeo,round(count/2*1.25), 1:end) ;
    dTAbsorb = saveddTAbsorb(indexF,indexGeo,round(count/2*1.25), 1:end);
    maxTemp =  savedmaxTemp(indexF,indexGeo,end, 1:end);
    
    temp = 0;
    for(indexT = 1:length(T(1:end)))
        temp(indexT) = T(indexT);
    end
    T = temp;
    
    temp = 0;
    for(indexT = 1:length(dTcond(1:end)))
        temp(indexT) = dTcond(indexT);
    end
    dTcond = temp;
    
    temp = 0;
    for(indexT = 1:length(dTAbsorb(1:end)))
        temp(indexT) = dTAbsorb(indexT);
    end
    dTAbsorb = temp;
    
    temp = 0;
    for(indexT = 1:length(maxTemp(1:end)))
        temp(indexT) = maxTemp(indexT);
    end
    maxTemp = temp;
    
    figure((indexGeo-1)*2+1);
    subplot(2,1,1);
    plot(T(1:end)-293);
    LayerPlot(indexGeo, LayerPos,max(max(T(1:end)-293)),0);
    title(['Fluence ', num2str(Fluence(indexFAgain)/10000), 'mj/cm2 -  Geometry ' , num2str(numGeo)]);
    ylabel('T - 293');
    xlabel('Distance');
    
    figure((indexGeo-1)*2+1);
    subplot(2,1,2);
    plot(maxTemp(1:end)-293);
    %LayerPlot(indexGeo, LayerPos,max(max(maxTemp(1:end)-293)));
    ylabel('Max Temp');
    xlabel('Count');
    
    
    figure((indexGeo-1)*2+2);
    
    subplot(2,1,1);
    plot(dTAbsorb(1:end));
    LayerPlot(indexGeo, LayerPos,max(max(dTAbsorb(1:end))),min(min(dTAbsorb(1:end))));
    title(['Fluence ', num2str(Fluence(indexFAgain)/10000), 'mj/cm2 - Geometry ' , num2str(numGeo)]);
    ylabel('dTAbsorb');
    xlabel('Distance');
    
    figure((indexGeo-1)*2+2);
    subplot(2,1,2);
    plot(dTcond(1:end));
    LayerPlot(indexGeo, LayerPos,max(max(dTcond(1:end))),min(min(dTcond(1:end))));
    ylabel('dTCond');
    xlabel('Distance');
    pause(1/60);
    
    
    
    
end


[minLayerPos,maxLayerPos] = layerLimits(NumLayer,order,Thick,Res);
for(indexGeo = 1:numGeo)
    for(indexY = 1:NumLayer(indexGeo))
        for(fluFlu = 1:length(Fluence))
            %saveTemp(indexF,indexGeo,count, 1:length(tempTemp(2,1:end)))
            fluLayerMaxTemp(indexGeo,indexY,fluFlu) = max(max(max(saveTemp(fluFlu,indexGeo,1:end,minLayerPos(indexGeo,indexY):maxLayerPos(indexGeo,indexY)))));
            
            [~,~,indexI,~] = ind2sub(size(saveTemp(fluFlu,indexGeo,1:end,minLayerPos(indexGeo,indexY):maxLayerPos(indexGeo,indexY))), find(saveTemp(fluFlu,indexGeo,1:end,minLayerPos(indexGeo,indexY):maxLayerPos(indexGeo,indexY)) == fluLayerMaxTemp(indexGeo,indexY,fluFlu),1));
            
            maxTempAreaPlot(fluFlu,indexGeo,indexY) = indexI;
            
            
            
        end
    end
end

for(indexGeo = 1:numGeo)
    for(indexY = 1:NumLayer(indexGeo))
        for(fluFlu = 1:length(Fluence))
            for(timeStep = 1:length(saveTemp(1,1,1:end,1)))
                
                %saveTemp(indexF,indexGeo,count, 1:length(tempTemp(2,1:end)))
                MaxTempTime(fluFlu,indexGeo,indexY,timeStep) = max(max(max(saveTemp(fluFlu,indexGeo,timeStep,minLayerPos(indexGeo,indexY):maxLayerPos(indexGeo,indexY)))));
                
                
            end
        end
    end
end
%plot(-293+reshape(MaxTempTime(indexFAgain,indexGeo,indexY,1:end),1,length(MaxTempTime(1,1,3,1:end))))

temp = indexGeo;
indexLayer = 3;

for(tempIndexGeo = 1:numGeo)
    indexGeo = tempIndexGeo;
    figure((temp-1)*2+2+(indexGeo-1)*3+1);
    title(['Temperature profile of geometry ',num2str(indexGeo),' when the ' , char(nameLayer(indexGeo,indexLayer)), ' layer is at its maximum temperature.'  ]);
    %what Geo want to plot
    
    %Plot maximum tempearture of which layer?
    
    
    [~, Ymat,~] = genXYVectorMat(NumLayer,Thick,Res,order,numGeo);
    Ymat=Ymat.*10^9-Thick(indexGeo,NumLayer(indexGeo)).*10^9;
    
    hold on;
    
    
    for(fluFlu = 1:length(Fluence))
        maxTemp =  saveTemp(fluFlu,indexGeo,maxTempAreaPlot(fluFlu,indexGeo,indexLayer), 1:end);
        plot(Ymat(indexGeo,1:end),reshape(maxTemp(1:end),1,length(maxTemp))-293,'LineWidth',1);
        
    end
    ylabel('Temperature (°C)');
    xlabel('Distance (nm)');
    bounds = ylim();
    ylim([min(bounds),max(bounds)*1.1]);
    bounds = ylim();
    LayerPlotRes(indexGeo, LayerPos,bounds(2),bounds(1), Ymat,NumLayer,Res*10^9,nameLayer);
    printLegend = {strcat(num2str(Fluence(1)/10),'mj/cm2')};
    
    for(fluFlu = 2:length(Fluence))
        printLegend = [printLegend,{strcat(num2str(Fluence(fluFlu)/10),'mj/cm2')}];
    end
    legend(printLegend);
    
    xlim([ 0, sum(Thick(indexGeo,1:end)).*10^9- Thick(indexGeo,1).*10^9-Thick(indexGeo,NumLayer(indexGeo)).*10^9]);
    hold off;
    
    %{

 subplot(1,1,1)
 hold on;
 plot(reshape(fluLayerMaxTemp(1,1,1:end),1,length(Fluence)));
 plot(reshape(fluLayerMaxTemp(1,2,1:end),1,length(Fluence)));
 plot(reshape(fluLayerMaxTemp(1,3,1:end),1,length(Fluence)));
 plot(reshape(fluLayerMaxTemp(1,4,1:end),1,length(Fluence)));
 plot(reshape(fluLayerMaxTemp(1,5,1:end),1,length(Fluence)));
 
 hold off;
    %}
    figure((temp-1)*2+2+(indexGeo-1)*3+2);
    hold on;
    
    for(indexY = 2:NumLayer(indexGeo)-1)
        plot((0:length(MaxTempTime(1,1,3,1:end))-1)*dt/10^-9,-293+reshape(MaxTempTime(indexFAgain,indexGeo,indexY,1:end),1,length(MaxTempTime(1,1,3,1:end))));
    end
    xlim([0,(length(MaxTempTime(1,1,3,1:end))-1)*dt/10^-9])
    ylabel('Temperature (°C)');
    xlabel('Time (ns)');
    hold off;
    legend(nameLayer(indexGeo,2:(NumLayer(indexGeo)-1)));
    title(['Maximum temperature of each layer in geometry ',num2str(indexGeo),' with incident fluence of ' ,num2str(Fluence(indexFAgain)/10000), ' mj/cm2.' ]);
    
    figure((temp-1)*2+2+(indexGeo-1)*3+3);
    
        hold on;
    for(indexFluFlu = 1:length(Fluence(1:end)))
        plot((0:length(MaxTempTime(1,1,3,1:end))-1)*dt/10^-9,-293+reshape(MaxTempTime(indexFluFlu,indexGeo,indexLayer,1:end),1,length(MaxTempTime(1,1,3,1:end))));
    end
    xlim([0,(length(MaxTempTime(1,1,3,1:end))-1)*dt/10^-9])
    ylabel('Temperature (°C)');
    xlabel('Time (ns)');
    hold off;
    legend(printLegend);
    title(['Maximum temperature of each layer in geometry ',num2str(indexGeo),' with incident fluence of ' ,num2str(Fluence(indexFluFlu)/10000), ' mj/cm2.' ]);
    
end


end


function [  ] = LayerPlot(indexGeo, LayerPos,maxVal,minVal )
hold on;
for(indexLayer = 2:length(LayerPos(1,1:end)))
    plot([LayerPos(indexGeo,indexLayer),LayerPos(indexGeo,indexLayer)], [minVal*1.5,maxVal*1.5],'--');
end
ylim([minVal,maxVal]);
hold off;
end

function [  ] = LayerPlotRes(indexGeo, LayerPos,maxVal,minVal, Ymat,NumLayer,Res,nameLayer)
hold on;
graphContent = get(gca);
colours = graphContent.ColorOrder();
for(indexLayer = 2:NumLayer(indexGeo))
    %plot([Ymat(indexGeo,LayerPos(indexGeo,indexLayer))+Res(indexGeo,indexLayer)/2,Ymat(indexGeo,LayerPos(indexGeo,indexLayer))+Res(indexGeo,indexLayer)/2], [minVal*1.5,maxVal*1.5],'--');
    plot([Ymat(indexGeo,LayerPos(indexGeo,indexLayer))+Res(indexGeo,indexLayer)/2,Ymat(indexGeo,LayerPos(indexGeo,indexLayer))+Res(indexGeo,indexLayer)/2], [minVal*1.5,maxVal*1.5],'--','Color','black','LineWidth',1);
    graphContent = get(gca);
    
    %text((Ymat(indexGeo,LayerPos(indexGeo,indexLayer))+Res(indexGeo,indexLayer)/2+Ymat(indexGeo,LayerPos(indexGeo,indexLayer-1))+Res(indexGeo,indexLayer-1)/2)/2,maxVal-maxVal/20,nameLayer(indexGeo,indexLayer-1),'HorizontalAlignment','center','Color',colours(mod(length(graphContent.Children())-(indexLayer-2)+6,7)+1,1:end));
    text((Ymat(indexGeo,LayerPos(indexGeo,indexLayer))+Res(indexGeo,indexLayer)/2+Ymat(indexGeo,LayerPos(indexGeo,indexLayer-1))+Res(indexGeo,indexLayer-1)/2)/2,maxVal-maxVal/20,strcat('',nameLayer(indexGeo,indexLayer-1)),'Rotation',90,'HorizontalAlignment','center','Color','black');
    
    
end
ylim([minVal,maxVal]);
hold off;
end

function [minLayerPos,maxLayerPos] = layerLimits(NumLayer,order,Thick,Res)

minLayerPos(1,1) = 1;
maxLayerPos(1,1) = 2;
for(indexX = 1:length(order))
    start = 0;
    for(indexLayer = NumLayer(order(indexX)):-1:1)
        
        minLayerPos(indexX,indexLayer) = start+1;
        maxLayerPos(indexX,indexLayer) = start+round(Thick(order(indexX),indexLayer)/Res(order(indexX),indexLayer));
        
        
        start = start+round(Thick(order(indexX),indexLayer)/Res(order(indexX),indexLayer));
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
