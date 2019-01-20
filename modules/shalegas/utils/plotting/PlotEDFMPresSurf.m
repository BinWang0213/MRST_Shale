function PlotEDFMPresSurf(fl,G,states,step,varargin)

opt = struct('ColorLim',[]);  % Compatibility only
opt = merge_options(opt, varargin{:});

figure('rend','painters','pos',[10 10 1000 600]);
set(gcf,'color','w');

nc = prod(G.cartDims);
pres=convertTo(states{step}.pressure, mega*Pascal);
%pres=convertTo(states{step}.pressure, psia);
[min(pres) max(pres)]
%plotCellData(G, pres,'EdgeColor','none')
line(fl(:,1:2:3)',fl(:,2:2:4)','Color','k','LineWidth',1);
colormap jet(64);
hcb = colorbar;
%view(0, 90); colorbar; axis equal tight, caxis([min(pres) max(pres)]);
view(0, 90); colorbar; axis equal tight;
if(~isempty(opt.ColorLim))
    caxis(opt.ColorLim);
else
    caxis([min(pres) max(pres)]);
end
set(gca,'FontSize',20);
%title(sprintf('EDFM %dx%d', G.cartDims(1),G.cartDims(2)))
%title(sprintf('EDFM %dx%d', G.cartDims(1),G.cartDims(2)))
grid2vtk(['Pressure_',num2str(step),'.vtk'],G,'Pressure',pres);
end