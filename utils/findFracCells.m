function [cellInx]=findFracCells(G,varargin)
%{
Find fracture cell idx in EDFM

Author:Bin Wang
Date: Nov.21.2018
%}
opt = struct('NFs', false, ...
             'HFs', false);  % Compatibility only

opt = merge_options(opt, varargin{:});


NumMatrixCells=G.Matrix.cells.num; NumTotalCells=G.cells.num;


if(opt.NFs==false && opt.HFs==false)
    cellInx=ones(NumTotalCells,1);
    cellInx(1:NumMatrixCells)=0;
else
    NumHFCells=0;
    for i=1:G.NumHFs
        name=strcat('Frac',num2str(i));
        NumHFCells=NumHFCells+G.FracGrid.(name).cells.num;
    end
    NumHFCells;
    NumNFCells=NumTotalCells-NumMatrixCells-NumHFCells;

    cellInx=zeros(NumTotalCells,1);
    if(opt.NFs==false)%Hydraulic Fracs
        cellInx(NumMatrixCells+1:NumMatrixCells+NumHFCells)=1;
    end
    if(opt.HFs==false)%Natural Fracs
        cellInx(end-NumNFCells+1:end)=1;
    end
end



end