function ComparePresSurf(fl,G,Gr,states,states_r,step)

figure

subplot(1,2,1);
set(gcf,'color','w');
pres=convertTo(states_r{step}.pressure, barsa);
plotCellData(Gr, pres,'EdgeColor','none')
line(fl(:,1:2:3)',fl(:,2:2:4)','Color','r','LineWidth',0.5);
colormap jet(25)
view(0, 90); colorbar; axis equal tight, cx = caxis();
title('Fully resolved solution')
set(gca,'FontSize',20);


subplot(1,2,2);
set(gcf,'color','w');
nc = prod(G.cartDims);
pres=convertTo(states{step}.pressure, barsa);
plotCellData(G, pres,'EdgeColor','none')
line(fl(:,1:2:3)',fl(:,2:2:4)','Color','r','LineWidth',0.5);
colormap jet(25)
view(0, 90); colorbar; axis equal tight, caxis(cx);
title('EDFM');
set(gca,'FontSize',20);

end