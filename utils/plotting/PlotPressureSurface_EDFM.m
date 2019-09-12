figure(222);
steps = [10, numSteps/4, numSteps/2, numSteps];
for i = 1 : numel(steps)
   subplot(2,2,i)
   plotCellData(G0, convertTo(states_e{steps(i)}.pressure(1:G0.cells.num, 1), barsa),...
       'EdgeColor', repmat(0.5, [1, 3]));
   line(fl(:,1:2:3)',fl(:,2:2:4)',1e-3*ones(2,size(fl,1)),'Color','r','LineWidth',0.5);
   view(2);

   caxis([40, 160])
   axis tight off;

   text(200, 170, -8, ...
        sprintf('%.1f days', convertTo(steps(i)*dt, day)), 'FontSize', 14)
end

h = colorbar('South', 'Position', [0.1, 0.05, 0.8, .025]);
colormap(jet);