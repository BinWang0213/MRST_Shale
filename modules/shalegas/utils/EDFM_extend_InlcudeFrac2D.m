function [fracp,ratio]=EDFM_extend_InlcudeFrac2D(edges,frac_endp)
%{
Locate the well idx for background tensor grid

Author:Bin Wang
Date: Nov.21.2018
%}

clf;
%Draw Cell
xx=[edges(:,1:2);edges(:,3:4)]; xx=unique(xx,'rows');
x=xx(:,1);y=xx(:,2);
cx = mean(x);cy = mean(y);
a = atan2(y - cy, x - cx);
[~, order] = sort(a);
x=x(order);y=y(order);
x=[x;x(1)];y=[y;y(1)];
plot(x,y,'bo-');
hold on;

%Draw Frac
plot(frac_endp(1:2:3),frac_endp(2:2:4),'r^-');
%hold off;
%pause(0.1);

% A vector along the ray from pto1 to pto2...
v=frac_endp(3:4)-frac_endp(1:2);
frac_length=pdist([frac_endp(3:4);frac_endp(1:2)],'euclidean');
box_size=pdist([[x(1) y(1)];[x(3) y(3)]],'euclidean');
factor_distance = box_size/frac_length;
new_pts1=frac_endp(1:2)-factor_distance*2*v;
new_pts2=frac_endp(3:4)+factor_distance*2*v;
plot([new_pts1(1) new_pts2(1)],[new_pts1(2) new_pts2(2)],'m^-');


%Find intersection pts
extended_endp=[new_pts1,new_pts2];
out = lineSegmentIntersect(extended_endp,edges);
xi = out.intMatrixX(out.intAdjacencyMatrix).';
yi = out.intMatrixY(out.intAdjacencyMatrix).';
new_frac_length=pdist([xi,yi],'euclidean');

%Draw Intersection Pts
scatter(xi,yi,'g');
hold off;

pause(0.5);

fracp = [xi,yi];
[~,ia] = unique(roundsd(fracp,14),'rows');
fracp = fracp(ia,:);
ratio=frac_length/new_frac_length;


end