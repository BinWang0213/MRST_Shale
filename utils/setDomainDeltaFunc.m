function [G,deltaM,deltaHF,deltaNF]=setDomainDeltaFunc(G)
%{
Set up domain delta function for EDFM

Author:Bin Wang
Date: Nov.21.2018
%}

deltaM=[];deltaHF=[];deltaNF=[];

if(~isfield(G,'Matrix')) %This is a Explicit Case
    NumTotalCells=G.cells.num;
    G.FracCellMask=zeros(NumTotalCells,1);
    cellInx = sub2ind(G.cartDims,G.FracCell.I,G.FracCell.J);
    G.FracCellMask(cellInx)=1;

    deltaM=~G.FracCellMask;
else
    G.FracCellMask=findFracCells(G); %All frac cells
    
    G.NFCellMask=findFracCells(G,'NFs',true); %Natural frac cells
    G.HFCellMask=findFracCells(G,'HFs',true); %Hydraulic frac cells

    deltaM=~G.FracCellMask;
    deltaHF=G.HFCellMask;
    deltaNF=G.NFCellMask;
end



end