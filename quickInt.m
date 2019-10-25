function [ outInt ] = quickInt( integrand ,IntOver)
% Symbolically ingerates over different layers  in integrand
% IntOver ( variable, like z)
for(index = 1:length(integrand))
    
    % direclty intergrating of z forw whatever reason is prone to errors
    % when given multiple examples of z in an equation. oneInt breaks them 
    % into seperable terms. 

    
    outInt(index) = oneInt(integrand(index), IntOver);
    
    
end
end

function [outInt] = oneInt(integrand, IntOver)
% Seperates out and intergrates each instance of z. removes inherent matlab
% errors. 

syms z t;
integrand = expand(integrand);
intChar = char(integrand);
lengChar = length(char(integrand));
leftBrac = strfind(intChar ,'(');
rightBrac = strfind(intChar ,')');
negSign = strfind(intChar ,' - ');
posSign = strfind(intChar ,' + ');

lengRight = length(rightBrac);
lengLeft = length(leftBrac);
lengNeg = length(negSign);
lengPos = length(posSign);
indexPos =1;
indexNeg = 1;
indexRight = 1;
indexLeft = 1;
countLeft = 1;
countRight = 1;
intComponent(1:round(length(strfind(intChar ,'z'))/2*1.25)) =  t-t;
count = 1;
positon =1;

RunMe = 1;

while(RunMe)
    if(indexLeft<lengLeft)
        if(and(indexPos<=lengPos ,indexNeg<=lengNeg ))
            if(posSign(indexPos) < negSign(indexNeg))
                if(posSign(indexPos)< leftBrac(indexLeft+1))
                    intComponent(count) =simplify(eval( intChar(positon:posSign(indexPos))));
                    positon = posSign(indexPos)+1;
                    intComponent(count);
                    count = count +1;
                    indexPos = indexPos +1;
                    indexLeft = indexLeft+1;
                    indexRight = indexRight+1;
                elseif(posSign(indexPos)>rightBrac(indexRight+1))
                    indexLeft = indexLeft+1;
                    indexRight = indexRight+1;
                else
                    indexPos = indexPos +1;
                    indexLeft = indexLeft+1;
                    indexRight = indexRight+1;
                end
            else
                if(negSign(indexNeg)< leftBrac(indexLeft+1))
                    intComponent(count) =simplify(eval( intChar(positon:negSign(indexNeg))));
                    positon = negSign(indexNeg)+1;
                    intComponent(count);
                    count = count +1;
                    indexNeg = indexNeg +1;
                    indexLeft = indexLeft+1;
                    indexRight = indexRight+1;
                elseif(negSign(indexNeg)>rightBrac(indexRight+1))
                    indexLeft = indexLeft+1;
                    indexRight = indexRight+1;
                else
                    indexNeg = indexNeg +1;
                    indexLeft = indexLeft+1;
                    indexRight = indexRight+1;
                end
            end
        elseif(indexPos>lengPos)
            if(negSign(indexNeg)< leftBrac(indexLeft+1))
                intComponent(count) =simplify(eval( intChar(positon:negSign(indexNeg))));
                positon = negSign(indexNeg)+1;
                intComponent(count);
                count = count +1;
                indexNeg = indexNeg +1;
                indexLeft = indexLeft+1;
                indexRight = indexRight+1;
            elseif(negSign(indexNeg)>rightBrac(indexRight+1))
                indexLeft = indexLeft+1;
                indexRight = indexRight+1;
            else
                indexNeg = indexNeg +1;
                indexLeft = indexLeft+1;
                indexRight = indexRight+1;
            end
        elseif(indexNeg>lengNeg)
            if(posSign(indexPos)< leftBrac(indexLeft+1))
                intComponent(count) =simplify(eval( intChar(positon:posSign(indexPos))));
                positon = posSign(indexPos)+1;
                intComponent(count);
                count = count +1;
                indexPos = indexPos +1;
                indexLeft = indexLeft+1;
                indexRight = indexRight+1;
            elseif(posSign(indexPos)>rightBrac(indexRight+1))
                indexLeft = indexLeft+1;
                indexRight = indexRight+1;
            else
                indexPos = indexPos +1;
                indexLeft = indexLeft+1;
                indexRight = indexRight+1;
            end
        end
    else
        if(length(intChar(positon:end))>1)
            intComponent(count) =simplify(eval(intChar(positon:end)));
        end
        
        RunMe = 0;
    end
    
    %if( mod(positon ,round(lengChar/1000)))
    %positon/lengChar;
    %end
    
end

intTemp = intComponent;

for(indexVar = 1:length(IntOver))
    
    intTemp = int(intTemp,IntOver(indexVar));
    
end
outInt = sum(intTemp);

end