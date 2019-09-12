function CompareLinePres(y_center,G,Gr,states,states_r,step)
%Currently, only suport extract center horizontal line y=center pts

figure;
set(gcf,'color','w');

I=[1:Gr.cartDims(1)]';
J=int16(repmat(Gr.cartDims(1)/2, [Gr.cartDims(1), 1]));
cellInx = int16(sub2ind(Gr.cartDims, I, J));

x=Gr.cells.centroids(1:Gr.cartDims(1),1);
p=convertTo(states_r{step}.pressure(cellInx),barsa);
plot(x,p,'-b','DisplayName','Fully resolved solution','LineWidth',2);
hold on;

I=[1:G.cartDims(1)]';
J=repmat(G.cartDims(1)/2, [G.cartDims(1), 1]);
cellInx = int16(sub2ind(G.cartDims, I, J));


x=G.cells.centroids(1:G.cartDims(1),1);
p=convertTo(states{step}.pressure(cellInx),barsa);

cellFracInx=[];
for i=1:length(x)
    Inx=markFracWellCell([x(i) y_center],G); %4.5 is the center 
    cellFracInx=[cellFracInx;Inx];
end
[C, ia, ic] = unique(cellFracInx);
if(ia(2)-ia(1)>1)
   ia(1)=ia(2)-1;
end
p_frac=convertTo(states{step}.pressure(C),barsa);
p(ia)=p_frac;

plot(x,p,'--r','DisplayName','EDFM','LineWidth',2);
set(gca,'FontSize',35);

xlabel('x [m]')
ylabel('Pressure [Bar]')
title('Pressure along y=4.5','FontSize',30)
xlim([int16(min(x)) int16(max(x))]);
legend('Location','north');
hold off;

end