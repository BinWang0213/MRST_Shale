function [G]=ExplicitFracGrid(physdim,...
    NumFracs,Frac_Spacing,Frac_halfLength,Frac_StartXY,varargin)
%{
Create rectilinear grid with log-refinement near the vertical fracture and
well

Mesh zone distribution
|----Out--|-frac_spacing-|....|-frac_spacing-|--Out----|

Arguments
---------
NX_FracRefine --  Number of refinement for a domain of a |frac_spaing|
NX_OutRefine  --  Number of refinement for domain out of frac zone
FracCell_IJK  --  Frac Cell index (i,j,k)
    
Author:Bin Wang
Date: Dec.21.2018
%}
%clf;

opt = struct('NX_FracRefine',3,'NX_OutRefine',-1,'FracCellSize',0.01*ft);  % Compatibility only
opt = merge_options(opt, varargin{:});

n_refinement=opt.NX_FracRefine;
if(opt.NX_OutRefine~=-1)
    n_refinement2=opt.NX_OutRefine;
else
    n_refinement2=n_refinement;
end

%Reservoir dimension
[DimX,DimY]=deal(physdim(1),physdim(2));
ratio=DimX/DimY;
%Background mesh pts for tensor grid
FracCell_I=[];
FracCell_J=[];

% X-direction for fracture
x = [];
%Left uniform space
x_uniform=symmetric_linspace(0,Frac_StartXY(1)-Frac_Spacing/2,n_refinement2,0);
x=[x x_uniform];
CellSize_X=x(2)-x(1);

%Log-scale refinement for fractures
center_cell_width=opt.FracCellSize; %Fracture aperature
x_refinement=symmetric_logspace(0,Frac_Spacing,n_refinement,center_cell_width);
if(opt.NX_FracRefine==3) %Special refinement for benchmakr with CMG
    DX_CMG=[0 4.8419    0.1531    0.0100    0.1531    4.8419].*ft;
    x_refinement=cumsum(DX_CMG);
end
for fi=1:NumFracs
    StartPts=Frac_StartXY(1)-Frac_Spacing/2+(fi-1)*Frac_Spacing;
    x=[x(1:end-1) x_refinement+StartPts];
    FracCell_I=[FracCell_I; numel(x)-numel(x_refinement)/2];
end

%Right unifrom space
x_uniform=symmetric_linspace(x(end),physdim(1),n_refinement2,0);
x=[x(1:end-1) x_uniform];

% Y-direction for center well
y_refinement_space=Frac_halfLength/2.0;
y_center=Frac_StartXY(2)+Frac_halfLength;

y=[];
%bottom uniform space
CellSize_Y=CellSize_X/ratio;
n_refinement2=DimY/CellSize_Y/5.7;
y_uniform=symmetric_linspace(0,y_center-Frac_halfLength,n_refinement2,0);
y=[y y_uniform];
CellSize_Y=y(end)-y(end-1);

%Log-scale refinement for well
n_refinement=round(2*Frac_halfLength/CellSize_Y/1.5);
y_refinement=symmetric_linspace(y_center-Frac_halfLength,...
    y_center+Frac_halfLength,...
    n_refinement,center_cell_width);
y=[y(1:end-1) y_refinement];
FracCell_J=[numel(y)-numel(y_refinement)+1:numel(y)-1]';
CellSize_Y=y(end)-y(end-1);

%Top uniform space
y_uniform=symmetric_linspace(y_center+Frac_halfLength,DimY,n_refinement2,0);
y=[y(1:end-1) y_uniform];
CellSize_Y=y(end)-y(end-1);

%% Genetrate Grid and plot it
if(numel(physdim)==2)
    G = tensorGrid(x, y);
else
    G = tensorGrid(x, y,[0 physdim(3)]);
end
G.FracCell.I=[];
for fi=1:NumFracs
    G.FracCell.I=[G.FracCell.I; repmat(FracCell_I(fi),numel(FracCell_J),1)];
end
G.FracCell.J=FracCell_J;
if(numel(physdim)==3)
    G.FracCell.K=repmat(1,G.cartDims(1)*G.cartDims(2),1);
end

%plotGrid(G);
end


function [x_space]=symmetric_linspace(a,b,n_refinement,center_space)
x_refine=(linspace(0,log10(11),n_refinement)-1)/10;

center=(a+b)/2;
x_space_right=scaledata(x_refine,center+center_space/2,b);
%x_space_right=rescale(x_refine,center+center_space/2,b);
x_space_left=-fliplr(x_space_right-center) + (center);

if(center_space==0.0)
    x_space_right=x_space_right(2:end);
end

x_space=[x_space_left x_space_right];
end

function [x_space]=symmetric_logspace(a,b,n_refinement,center_space)
%     Log space pattern
% |    | || center || |    |
% a                        b

%Create a log space range from (0,1)
x_refine=(logspace(0,log10(11),n_refinement)-1)/10;

center=(a+b)/2;
%x_space_right=rescale(x_refine,center+center_space/2,b);
x_space_right=scaledata(x_refine,center+center_space/2,b);
x_space_left=-fliplr(x_space_right-center) + (center);

if(center_space==0.0)
    x_space_right=x_space_right(2:end);
end

x_space=[x_space_left x_space_right];
end

function dataout = scaledata(datain,minval,maxval)
%
% Program to scale the values of a matrix from a user specified minimum to a user specified maximum
%
% Usage:
% outputData = scaleData(inputData,minVal,maxVal);
%
% Example:
% a = [1 2 3 4 5];
% a_out = scaledata(a,0,1);
% 
% Output obtained: 
%            0    0.1111    0.2222    0.3333    0.4444
%       0.5556    0.6667    0.7778    0.8889    1.0000
%
% Program written by:
% Aniruddha Kembhavi, July 11, 2007
dataout = datain - min(datain(:));
dataout = (dataout/range(dataout(:)))*(maxval-minval);
dataout = dataout + minval;
end