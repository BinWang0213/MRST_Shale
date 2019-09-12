function drawCICalc(edges,xi,yi,frac_endp)
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
%Draw Intersection Pts
scatter(xi,yi,'k');
%Draw Frac
plot(frac_endp(1:2:3),frac_endp(2:2:4),'r^-');
%hold off;
pause(0.1);

end