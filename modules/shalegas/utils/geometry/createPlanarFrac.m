function [fl,length]=createPlanarFrac(xy_center,half_length)
%{
Create a planar frac based on x^2 func

Arguments
---------
xy_center -- (a,b) center of planar frac

Author:Bin Wang
Date: Nov.21.2018
%}


X=[xy_center(1) xy_center(1)];
Y=[xy_center(2)-half_length xy_center(2)+half_length];
line=[X' Y'];
length=2*half_length;

%Change [Pts1 Pts2 Pts3]->[Pts1 Pts2 Pts2 Pts3]
fl=[line(1:end-1,:) line(2:end,:)];


%plot(X,Y,'bo-','LineWidth', 1);
%axis equal tight;
%hold on;

end