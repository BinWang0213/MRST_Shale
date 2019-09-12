%clf
hp = loglog(...
   convertTo([sol(2:end).time], day), ...
   convertTo(-[sol(2:end).qS], meter^3/day),'DisplayName','Explicit');

set(hp, 'LineStyle', '-', 'Marker', 'o', 'LineWidth', 1);

xlabel('time [days]');
ylabel('rate [m^3/day]');
legend;