function [cellInx]=markFracWellCell(xy_wells,G)
%{
Locate the well idx for fracture cell in EDFM

Author:Bin Wang
Date: Nov.21.2018
%}

%Collect all nodes along each fracture
NumFracs=length(fieldnames(G.FracGrid));
fields = fieldnames(G.FracGrid);

cell_center=[];
for i = 1:numel(fields)
  frac=G.FracGrid.(fields{i});
  %Check node pattern
  sum1=sum(frac.nodes.coords(2,:)-frac.nodes.coords(1,:));
  sum2=sum(frac.nodes.coords(3,:)-frac.nodes.coords(2,:));
  if(abs(sum1-sum2)<1e-8)
      endpts=frac.nodes.coords(1:end/2,:);
  else
      endpts=frac.nodes.coords(1:2:end,:);
  end
  centerpts=0.5 * (endpts(1:end-1,:) + endpts(2:end,:));
  cell_center=[cell_center; centerpts];
end

%Find the cell id based on nearest pts
Idx=knnsearch(cell_center,xy_wells,'K',1);

cellInx=Idx+G.Matrix.cells.num;

end