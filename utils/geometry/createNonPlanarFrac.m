function [fl,length]=createNonPlanarFrac(xy_center,degree,v_offset,type,varargin)
%{
Create a non-planar frac based on x^2 func

func:
(x-a)=-c*(y-b)^2

Arguments
---------
xy_center -- (a,b) center of non-planar frac
degree    --  c    curvture of curve (0.1-1 gives good shape)
v_offset  --  y(-v_offset,v_offset)vertical range of curve
type -- 'left' left curve
             'right' right curve
             'random' random curve

Author:Bin Wang
Date: Nov.21.2018
%}

opt = struct('ShapeFactor',1.14);  % Compatibility only
opt = merge_options(opt, varargin{:});

NumPts=3;
a=xy_center(1);
b=xy_center(2);
d=v_offset;

%Y=linspace(b-v_offset,b,NumPts);
%Y=[Y(1:end-1) linspace(b,b+v_offset,NumPts)];
Y=linspace(b-v_offset,b+v_offset,NumPts*2);%Well needs to be put on cell center


if(strcmp(type,'left'))
    c=degree*0.1;
    X=-1*c.*(Y-b).^2+a;
elseif(strcmp(type,'right'))
    c=degree*0.1;
    X=1*c.*(Y-b).^2+a;
else
    c=degree*v_offset;
    random_magn=(1.5-(-1.5)).*rand + (-1.5);
    X=random_magn*c.*sin(opt.ShapeFactor.*(Y-b))+a;
end


line=[X' Y'];
%Change [Pts1 Pts2 Pts3]->[Pts1 Pts2 Pts2 Pts3]
fl=[line(1:end-1,:) line(2:end,:)];

%Compute total length of line
length=calcFracsLength(fl);

%plot(X,Y,'bo-','LineWidth', 1);
%axis equal tight;
%hold on;

end