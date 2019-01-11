function [length]=calcFracsLength(fl)
    %{
    Compute total length of fractures
    fl=[x1 y1 x2 y2;
        x1 y1 x2 y2]
    
    Author:Bin Wang
    Date: Nov.21.2018
    %}
    
    
    d=[diff(fl(:,1:2:end),1,2) diff(fl(:,2:2:end),1,2)];
    length=sum(sqrt(sum(d.*d,2)));
    
end