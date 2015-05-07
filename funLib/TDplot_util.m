function TDplot_util(voltageX_BD2,currentX_BD2,xAxis)

xAxisType = 'Linear';
stringXAxis='[ns]';
xAxisNorm = 1e-9;

set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 25)

% Change default text fonts.
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 25)

figure,
plot(xAxis/xAxisNorm,voltageX_BD2,'b','linewidth',2);hold on;
xlabel(stringXAxis);
ylabel('|Vin| [V]');
xlim([xAxis(1) xAxis(end)]/xAxisNorm);
ylim([0 1])
set(gca, 'XScale',xAxisType);
hold off;

figure,
plot(xAxis/xAxisNorm,currentX_BD2,'r','linewidth',2);hold on;
xlabel(stringXAxis);
ylabel('|I_R| [A]');
xlim([xAxis(1) xAxis(end)]/xAxisNorm);
set(gca, 'XScale',xAxisType);
hold off;

figure,
plot(xAxis/xAxisNorm,voltageX_BD2./currentX_BD2,'k','linewidth',2);hold on;
xlabel(stringXAxis);
ylabel('|v(t)/i(t)| [\Omega]');
xlim([xAxis(1) xAxis(end)]/xAxisNorm);
ylim([0 10])
set(gca, 'XScale',xAxisType);
hold off;

end