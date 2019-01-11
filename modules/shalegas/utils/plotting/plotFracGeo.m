function plotFracGeo(physdim,fl,xy_wells,varargin)
%When run 2D model, we need scale production rate by formation thickness
opt = struct('NaturalFrac', [],'FigSize',400,'Title','Fracture Map'); 
opt = merge_options(opt, varargin{:});

figure('rend','painters','pos',[10 10 opt.FigSize*(physdim(1)/physdim(2)) opt.FigSize]);
set(gcf,'color','w');

%Plot boundary of Domain
X=[0 physdim(1) physdim(1) 0 0];
Y=[0 0 physdim(2) physdim(2) 0];

plot(X,Y,'k-','LineWidth',0.5);
hold on;

%Plot fracture
NumLines=size(fl,1);
NumHFs=NumLines;

%Plot Natural fracture
if(~isempty(opt.NaturalFrac))
    fnm=opt.NaturalFrac;
    NumNFs=size(fnm,1);
    NumHFs=NumLines-NumNFs;
    for i=1:NumNFs
        x = [fnm(i,1);fnm(i,3)];
        y = [fnm(i,2);fnm(i,4)];
        line(x,y,'Color',[0.5 0.5 0.5],'LineWidth',1.5);
        %text(mean(x),mean(y),num2str(i)); % show line numbering
    end
end

%Plot Hydraulic fractures
% Color Table
% https://www.cimat.mx/~max/InformaticaAplicadaII_2017/bibliografia/MATLAB_reference.html
for i=1:NumHFs
    x = [fl(i,1);fl(i,3)];
    y = [fl(i,2);fl(i,4)];
    line(x,y,'Color','b','LineWidth',5);
    %text(mean(x),mean(y),num2str(i)); % show line numbering
end

%Plot wells
if(~isempty(xy_wells))
    scatter(xy_wells(:,1),xy_wells(:,2),35,'k','filled')
end

%Figure improve
xlim([0 physdim(1)])
ylim([0 physdim(2)])
set(gca,'FontSize',20);
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [0 0] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XGrid'       , 'off'      )
%  'XColor'      , [.3 .3 .3], ...
%  'YColor'      , [.3 .3 .3], ...
%  'YTick'       , 0:500:2500, ...
%  'LineWidth'   , 1         );
grid on;
title(opt.Title,'FontWeight','normal')


end