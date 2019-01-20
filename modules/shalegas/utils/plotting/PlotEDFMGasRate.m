function PlotEDFMGasRate(time_list,ws,varargin)
%When run 2D model, we need scale production rate by formation thickness
opt = struct('Formation_thickness', 1, ...
             'Reference_data', '',...
             'Reference_data2', '',...
             'CumPlot',0,...
             'XLog',0,...
             'YLog',0,...
             'Xlim',[],...
             'Ylim',[],...
             'XUnit',day,...
             'YUnit',meter^3/day,...
             'LogLog',0);  % Compatibility only
opt = merge_options(opt, varargin{:});

h=opt.Formation_thickness;

ws_mat=cell2mat(ws);
figure('rend','painters','pos',[10 10 1000 600]);
set(gcf,'color','w');

time=convertTo(time_list*day,opt.XUnit);
PROs=convertTo(-[ws_mat(:).qWs]*h, opt.YUnit);
PROr=convertTo(-[ws_mat(:).qWr]*h, opt.YUnit);

plot(time,PROs,'b-', 'LineWidth', 2)
%plot(time_list, convertTo(-[ws_mat(:).qWs]*h, 1e5*ft^3/day),'b-', 'LineWidth', 2)
hold on;
if(~isempty(opt.Reference_data))
    %fname=strcat(pwd,opt.Reference_data);
    fname=opt.Reference_data;
    M = csvread(fname,1);
    %plot(M(:,1),M(:,2),'r-.', 'LineWidth', 2)
    plot(M(:,1),M(:,2),'ro', 'LineWidth', 2)
    legend('MRST-Shale','Reference Solution')
end
if(~isempty(opt.Reference_data2))
    %fname=strcat(pwd,opt.Reference_data2);
    fname=opt.Reference_data2;
    M = csvread(fname,1);
    plot(M(:,1),M(:,2),'ko', 'LineWidth', 2)
    %scatter(M(:,1),M(:,3),50,'r','LineWidth',1.5);
    %scatter(M(:,1)/365,M(:,2)./1e6,50,'r','LineWidth',1.5);
    legend('MRST-Shale','Field Data','Field Data2')
end
hold off;



set(gca,'FontSize',25);
xlabel('Time [Days]')
ylabel('Gas Production Rate [m^3/day]')
if(~isempty(opt.Xlim))
    xlim(opt.Xlim);
end
if(~isempty(opt.Ylim))
    ylim(opt.Ylim);
end
if(opt.LogLog==1)
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
end
if(opt.XLog==1)
   set(gca, 'XScale', 'log')
end
if(opt.YLog==1)
   set(gca, 'YScale', 'log')
end


if(opt.CumPlot==1)%Plot Cumulaive production plot
    figure('rend','painters','pos',[10 10 1000 600]);
    set(gcf,'color','w');
    
    cumPro=cumtrapz(time_list, PROs);
    plot(time/365,cumPro/1e6,'b-', 'LineWidth', 2)
    hold on;
    if(~isempty(opt.Reference_data))
        cumPro_ref=cumtrapz(M(:,1),M(:,2));
        plot(M(:,1)/365,cumPro_ref/1e6,'ro', 'LineWidth', 2)
        legend('MRST-Shale','Reference Solution')
    end
    hold off;
    
    Error=(cumPro(end)-cumPro_ref(end))/cumPro_ref(end)
    
    set(gca,'FontSize',25);
    xlabel('Time [Years]')
    ylabel('Cum Production Rate [m^3/day]')
end

Case4_Nonplanar_NFs_Geomechs_PROS_EDFM=convertTo(-[ws_mat(:).qWs]*h, meter^3/day);
cumPro=cumtrapz(time_list, convertTo(-[ws_mat(:).qWs]*h, meter^3/day));
CumulativeProduction=cumPro(end)
%[mean(M(15:end,2)./PROs(15:end)') mean(M(15:end,3)./PROr(15:end)')]
csvwrite('LGR150New.csv',[time PROs'],1);
end