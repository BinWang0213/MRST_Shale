function ComparePresContour(fl,G,Gr,states,states_r,step)

figure
pres=convertTo(states_r{step}.pressure, barsa);
contour(reshape(Gr.cells.centroids(:,1),Gr.cartDims),...
   reshape(Gr.cells.centroids(:,2),Gr.cartDims),...
   reshape(pres,Gr.cartDims),10);
hold on

nc = prod(G.cartDims);
pres=convertTo(states{step}.pressure, barsa);
contour(reshape(G.cells.centroids(1:nc,1),G.cartDims),...
   reshape(G.cells.centroids(1:nc,2),G.cartDims),...
   reshape(pres(1:nc),G.cartDims),10,'--','LineWidth',2);
line(fl(:,1:2:3)',fl(:,2:2:4)','Color','r','LineWidth',0.5);
hold off;
axis equal tight;
legend('Fully resolved','Embedded fractures');
set(gca,'FontSize',20);
end