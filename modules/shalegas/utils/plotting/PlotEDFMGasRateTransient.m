function PlotEDFMGasRateTransient(time_list,ws,varargin)
%When run 2D model, we need scale production rate by formation thickness
opt = struct('Formation_thickness', 1, ...
             'Reference_data', 'filename.csv',...
             'LogLog',1);  % Compatibility only
opt = merge_options(opt, varargin{:});

h=opt.Formation_thickness;
fname=strcat(pwd,opt.Reference_data);

ws_mat=cell2mat(ws);
figure('rend','painters','pos',[10 10 800 600]);
set(gcf,'color','w');

y0=convertTo(-[ws_mat(:).qGs]*h, meter^3/day);
plot(time_list,y0 ,'r.','MarkerSize',25)
hold off;

set(gca,'FontSize',20);
xlabel('Time [Days]')
%xlabel('Time [Years]')
ylabel('Gas Production Rate [m^3/day]')
xlim([1e-4 1e4]);
ylim([1e2 1e5]);
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
grid on;
set(gca,'YMinorGrid','off')
set(gca,'XMinorGrid','off')

GasRateProEDFM_NaturalFrac250=y0;
time_list_Base=time_list;

end