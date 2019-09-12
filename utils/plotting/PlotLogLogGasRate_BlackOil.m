clf
ws_mat=cell2mat(ws);
hp = loglog(...
   time_list(:), ...
   convertTo(-[ws_mat(:).qGs], meter^3/day),'DisplayName','BlackOil');
hold on;
set(hp, 'LineStyle', '-', 'Marker', 'o', 'LineWidth', 1);

xlabel('time [days]');
ylabel('rate [m^3/day]');