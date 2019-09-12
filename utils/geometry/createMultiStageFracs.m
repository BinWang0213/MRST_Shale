function [fl,xy_wells]=createMultiStageFracs(NumFracs,Frac_Spacing,Frac_halfLength,Frac_StartXY)
%{
Create multi-stage hydraulic fractures and its center well location

Arguments
---------
fl        --  NFx4 array of fracture segment length 
xy_wells  --  NFx2 array of center well location

Author:Bin Wang
Date: Nov.21.2018
%}

fl=[];
xy_wells=[];
for fi=1:NumFracs
    xy_center=[Frac_StartXY(1)+(fi-1)*Frac_Spacing Frac_StartXY(2)+Frac_halfLength];
    [fl_new,~]=createPlanarFrac(xy_center,Frac_halfLength);
    
    fl=[fl; fl_new];% fractures lines in [x1 y1 x2 y2] format.
    xy_wells=[xy_wells; xy_center];
end

end