figure(223);
steps = [10, numSteps/4, numSteps/2, numSteps];
for i = 1 : numel(steps)
   subplot(2,2,i)
   plotCellData(G, convertTo(states{steps(i)}.pressure, barsa), ...
                'EdgeColor', repmat(0.5, [1, 3]));
   view(2);

   caxis([40, 160])
   axis tight off;
   
   text(200, 170, -8, ...
        sprintf('%.1f days', convertTo(steps(i)*dt, day)), 'FontSize', 14)
end

h = colorbar('South', 'Position', [0.1, 0.05, 0.8, .025]);
colormap(jet);