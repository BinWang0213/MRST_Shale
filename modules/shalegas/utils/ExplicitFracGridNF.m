function [G,fl]=ExplicitFracGridNF(physdim,...
    NumFracs,Frac_Spacing,Frac_halfLength,Frac_StartXY,varargin)
%{
Create rectilinear grid with log-refinement near the 
vertical hydraulic fracture, horizontal natural fracture and well

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

opt = struct('NX_FracRefine',3, ...
             'NX_OutRefine',-1,...
             'FracCellSize',0.01*ft,...
             'FracCellSize_Y',-1,...
             'NY_OutRefine',-1,...
             'NumNFs',-1,...
             'NFCellSize',-1,...
             'NY_Refine',-1,...
             'NY_LogRefine',true,...
             'NF_FracRefine',-1,...
             'NF_Spacing',-1,...
             'NF_Length',-1,...
             'NF_RepeatPatternSpace',-1,...
             'NF_StartXY',-1);
opt = merge_options(opt, varargin{:});

n_refinement=opt.NX_FracRefine;
if(opt.NX_OutRefine~=-1)
    n_refinement2=opt.NX_OutRefine;
else
    n_refinement2=n_refinement;
end

if(opt.FracCellSize_Y==-1)
    FracCellSize_Y=opt.FracCellSize;
else
    FracCellSize_Y=opt.FracCellSize_Y;
end

if(opt.NY_OutRefine==-1), n_refinement2_y=n_refinement2;
else, n_refinement2_y=opt.NY_OutRefine;end

%Reservoir dimension
[DimX,DimY]=deal(physdim(1),physdim(2));
ratio=DimX/DimY;
%Background mesh pts for tensor grid
fl=[];
FracCell_I=[];
FracCell_J=[];

% X-direction for fracture
x = [];
%Left uniform space
x_uniform=symmetric_linspace(0,Frac_StartXY(1)-Frac_Spacing/2,n_refinement2,0);
x=[x x_uniform];

%Log-scale refinement for fractures
center_cell_width=opt.FracCellSize; %Fracture aperature
x_refinement=symmetric_logspace(0,Frac_Spacing,n_refinement,center_cell_width);

for fi=1:NumFracs
    StartPts=Frac_StartXY(1)-Frac_Spacing/2+(fi-1)*Frac_Spacing;
    x=[x(1:end-1) x_refinement+StartPts];
    FracCell_I=[FracCell_I; numel(x)-numel(x_refinement)/2];
    
    %Collect hydraulic fracture line info
    fl=[fl;...
        StartPts+Frac_Spacing/2 Frac_StartXY(2)...
        StartPts+Frac_Spacing/2 Frac_StartXY(2)+2*Frac_halfLength];
end
CellSize_X=x(2)-x(1)

%Right unifrom space
x_uniform=symmetric_linspace(x(end),physdim(1),n_refinement2,0);
if(numel(x_uniform)>0) 
    x=[x(1:end-1) x_uniform];
end

% Y-direction for center well
y_refinement_space=Frac_halfLength/2.0;
y_center=Frac_StartXY(2)+Frac_halfLength;

y=[];
%bottom uniform space
%Option 1. Match with cell size in X
CellSize_Y=CellSize_X/ratio;
%n_refinement2=DimY/CellSize_Y/5.7;
n_refinement2=DimY/CellSize_Y/6.5; 
%Option2. Defined by user
n_refinement2=n_refinement2_y;
y_uniform=symmetric_linspace(0,y_center-Frac_halfLength,n_refinement2,0);
y=[y y_uniform];
CellSize_Y=y(end)-y(end-1);

%Log-scale refinement for well
if(opt.NumNFs==-1)
    center_cell_width=FracCellSize_Y; %Well cell centerlizer
    n_refinement=round(2*Frac_halfLength/CellSize_Y/1.5);
    y_refinement=symmetric_linspace(y_center-Frac_halfLength,...
        y_center+Frac_halfLength,...
        n_refinement,center_cell_width);
    if(numel(y_refinement)>0) 
        y=[y(1:end-1) y_refinement];
    end
    FracCell_J=[numel(y)-numel(y_refinement)+1:numel(y)-1]';
    CellSize_Y=y(end)-y(end-1)
end

if(opt.NumNFs>0)
    %Unifrom refine Y
    y_refinement=symmetric_linspace(0,opt.NF_Spacing,opt.NY_Refine,FracCellSize_Y);
    %Log refine
    if(opt.NY_LogRefine==true)
        y_refinement=symmetric_logspace(0,opt.NF_Spacing,opt.NY_Refine,FracCellSize_Y);
    end

    %Log refine single fracture cell
    center_cell_width=opt.NFCellSize; %Well cell centerlizer
    y_frac_cell_refinement=symmetric_logspace(0,FracCellSize_Y,opt.NF_FracRefine,center_cell_width);
    y_refinement=[y_refinement(1:numel(y_refinement)/2-1) ...
                  y_frac_cell_refinement+y_refinement(numel(y_refinement)/2) ...
                  y_refinement(numel(y_refinement)/2+2:end)];
    
    for fi=1:opt.NumNFs
        StartPts=opt.NF_StartXY(2)+(fi-1)*opt.NF_Spacing;
        y=[y(1:end-1) y_refinement+StartPts];
        FracCell_J=[FracCell_J; numel(y)-numel(y_refinement)/2];
        
        %Collect natural fracture line info
        fl=[fl;...
        opt.NF_StartXY(1) StartPts+opt.NF_Spacing/2 ...
        opt.NF_StartXY(1)+opt.NF_Length StartPts+opt.NF_Spacing/2 ];
        if(opt.NF_RepeatPatternSpace~=-1)
            space=opt.NF_RepeatPatternSpace+opt.NF_Length;
            fl=[fl;...
            space+opt.NF_StartXY(1) StartPts+opt.NF_Spacing/2 ...
            space+opt.NF_StartXY(1)+opt.NF_Length StartPts+opt.NF_Spacing/2 ];
        end
    end    
    
    %Repeat NF on other side
    if(opt.NF_RepeatPatternSpace~=-1),opt.NumNFs=opt.NumNFs*2; end
end

%Top uniform space
y_uniform=symmetric_linspace(y_center+Frac_halfLength,DimY,n_refinement2,0);
if(numel(y_uniform)>0) 
    y=[y(1:end-1) y_uniform];
end
CellSize_Y=y(end)-y(end-1);

%% Genetrate Grid and plot it
if(numel(physdim)==2)
    G = tensorGrid(x, y);
else
    G = tensorGrid(x, y,[0 physdim(3)]);
end
G.FracCell.I=[];
G.FracCell.J=[];
for fi=1:NumFracs
    G.FracCell.I=[G.FracCell.I; repmat(FracCell_I(fi),numel(FracCell_J),1)];
    G.FracCell.J=[G.FracCell.J; FracCell_J];
end
if(numel(physdim)==3)
    G.FracCell.K=repmat(1,G.cartDims(1)*G.cartDims(2),1);
end

%Find Frac Cell index
G=markFracCells(G,NumFracs,opt.NumNFs,fl);
%plotGrid(G);
end


function [G]=markFracCells(G,NumFracs,NumNFs,fl)

[NX,NY]=deal(G.cartDims(1),G.cartDims(2));
G_coordX=G.nodes.coords(1:NX+1,1); 
G_coordY=G.nodes.coords(1:NX+1:(NX+1)*(NY+1),2);

avg_DX=mean(diff(G_coordX))/2;
avg_DY=mean(diff(G_coordY))/2;
tol=1e-5;

FracCell_I=[];
FracCell_J=[];

%Hydraulic fractures
for fi=1:NumFracs
    pts_on_frac=[fl(fi,1)+tol fl(fi,2)+tol; fl(fi,3)+tol fl(fi,4)-tol];
    [IJ_cells]=markCellbyXY(pts_on_frac,G);
    J_low=min(IJ_cells(:,2));J_up=max(IJ_cells(:,2));
    FracCell_J=[FracCell_J  J_low:J_up];
    FracCell_I=[FracCell_I; repmat(IJ_cells(1,1),J_up-J_low+1,1)];
end
FracCell_I=FracCell_I';
FracCell_J=FracCell_J';

%Natural fractures
NFStartIdx=numel(FracCell_I)+1;
if(NumNFs>0)
    for fi=1:NumNFs
        pts_on_frac=[fl(NumFracs+fi,1)+tol fl(NumFracs+fi,2)+tol;...
            fl(NumFracs+fi,3)-tol fl(NumFracs+fi,4)+tol];
        [IJ_cells]=markCellbyXY(pts_on_frac,G);
        I_low=min(IJ_cells(:,1));I_up=max(IJ_cells(:,1));
        FracCell_I=[FracCell_I  I_low:I_up];
        FracCell_J=[FracCell_J; repmat(IJ_cells(1,2),I_up-I_low+1,1)];
    end
end

G.FracCell.I=FracCell_I;
G.FracCell.J=FracCell_J';
G.FracCell.NFStartIdx=NFStartIdx;
if(numel(G.cartDims)==3)
    G.FracCell.K=repmat(1,G.cartDims(1)*G.cartDims(2),1);
end

%debug
%{
plotFracGeo([1990*ft 1990*ft 150*ft],fl,[],'FigSize',600);
field=repmat(0,NX,NY);
cellInx = sub2ind(G.cartDims, FracCell_I', FracCell_J);
field(cellInx)=1;
plotCellData (G , field(:),'FaceAlpha',0.5);
colorbar ('horiz'); view (2); axis equal tight ;
%}

end

function [xy_centers]=splitline(line,n)
%Find the mid point of subdivided cells along a line segment 
%https://www.mathworks.com/matlabcentral/answers/218445-how-can-i-divide-a-line-into-equal-parts-and-divide-them-into-small-circles-of-equal-radius

x=[line(1) line(3)];
y=[line(2) line(4)];

n;  % #segments
deltax = (x(2) - x(1))/n;
deltay = (y(2) - y(1))/n;
halfdx = deltax/2;
halfdy = deltay/2;
xcents = x(1) + halfdx + (0:n-1)*deltax;
ycents = y(1) + halfdy + (0:n-1)*deltay;
%plot(x, y, 'bs', xcents, ycents, 'r*');
xy_centers=[xcents' ycents'];
end

function [x_space]=symmetric_linspace(a,b,n_refinement,center_space)

%Speical case handle, interval size = center cell size or no refinement
if (abs(a-b)<1e-9)
    x_space=[];
elseif(abs(abs(a-b)-abs(center_space))<1e-9 || n_refinement<=1)
    x_space=[a b];
else
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

end

function [x_space]=symmetric_logspace(a,b,n_refinement,center_space)
%     Log space pattern
% |    | || center || |    |
% a                        b

%Speical case handle, interval size = center cell size or no refinement
if (abs(a-b)<1e-9)
   x_space=[];
elseif(abs(abs(a-b)-abs(center_space))<1e-9 || n_refinement<=1)
   x_space=[a b];
else
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