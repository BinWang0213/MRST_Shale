function [IJ_wells]=markWellCell(xy_wells,G)
%{
Locate the well idx for background tensor grid

Author:Bin Wang
Date: Nov.21.2018
%}

[NX,NY]=deal(G.cartDims(1),G.cartDims(2));

%Row order in tensor grid
G_coordX=G.nodes.coords(1:NX+1,1); 
G_coordY=G.nodes.coords(1:NX+1:end,2);

IJ_wells=zeros('like',xy_wells);

NumWells=size(xy_wells,1);
for wi=1:NumWells
    I=find(G_coordX> xy_wells(wi,1)); I=min(I)-1;
    J=find(G_coordY> xy_wells(wi,2)); J=min(J)-1;
    IJ_wells(wi,1)=I;
    IJ_wells(wi,2)=J;
end


end