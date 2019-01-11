ws_mat=cell2mat(ws_e);
hp = loglog(...
   time_list(:), ...
   convertTo(-[ws_mat(:).qWs], meter^3/day),'DisplayName','EDFM');
hold on;
set(hp, 'LineStyle', '-', 'Marker', 'o', 'LineWidth', 1);

xlabel('time [days]');
ylabel('rate [m^3/day]');
legend;