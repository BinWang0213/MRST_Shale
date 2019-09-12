function [fl,xy_wells]=createMultiStageNonPlanarFracs...
    (NumFracs,Frac_Spacing,Frac_halfLength,Frac_StartXY,varargin)
%{
Create multi-stage non-planar hydraulic fractures
The total length is constrained by planar case

Arguments
---------
fl        --  NFx4 array of fracture segment length 
xy_wells  --  NFx2 array of center well location

Author:Bin Wang
Date: Dec.21.2018
%}

opt = struct('Severity',0.1,...
             'ShapeFactor',1.14,...
             'YRange',1.1);
opt = merge_options(opt, varargin{:});
    
    
TargetLength=NumFracs*2*Frac_halfLength;
fl=[];
xy_wells=[];

%Side Non-planar frac
for fi=1:NumFracs/2
    xy_center=[Frac_StartXY(1)+(2*fi-1-1)*Frac_Spacing Frac_StartXY(2)+Frac_halfLength];
    [fl_new,len]=createNonPlanarFrac(xy_center,opt.Severity,...
        Frac_halfLength*opt.YRange,'random','ShapeFactor',opt.ShapeFactor);
    TargetLength=TargetLength-len;
    fl=[fl; fl_new];% fractures lines in [x1 y1 x2 y2] format.
    xy_wells=[xy_wells; xy_center];
end

%Odd number of fractures final frac
if(mod(NumFracs,2)==1)
    fi=NumFracs-1;
    xy_center=[Frac_StartXY(1)+fi*Frac_Spacing Frac_StartXY(2)+Frac_halfLength];
    [fl_new,len]=createNonPlanarFrac(xy_center,0.3,...
        Frac_halfLength*opt.YRange,'random','ShapeFactor',opt.ShapeFactor);
    TargetLength=TargetLength-len;
    fl=[fl; fl_new];% fractures lines in [x1 y1 x2 y2] format.
    xy_wells=[xy_wells; xy_center];
end

%Center planar frac
NumFrac_left=(NumFracs-size(xy_wells,1));
for fi=1:NumFrac_left
    xy_center=[Frac_StartXY(1)+(2*fi-1)*Frac_Spacing Frac_StartXY(2)+Frac_halfLength];
    [fl_new,len]=createPlanarFrac(xy_center,TargetLength/NumFrac_left/2);
    fl=[fl; fl_new];% fractures lines in [x1 y1 x2 y2] format.
    xy_wells=[xy_wells; xy_center];
end
TargetLength=TargetLength-len*fi;



end