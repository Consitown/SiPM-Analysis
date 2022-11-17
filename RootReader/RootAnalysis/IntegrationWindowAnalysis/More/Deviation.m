%x=FWValues.x;
%y=FWValues.y;
%z=FWValues.z;
%h=FWValues.h;

x=CorrectionValuesFormatted.x;
y=CorrectionValuesFormatted.y;
z=CorrectionValuesFormatted.z;
h=CorrectionValuesFormatted.h;


ymean = mean(y);

cn = 500;                                        % Number Of Colors (Scale AsAppropriate)
cm = colormap(jet(cn));                                         % Define Colormap
figure(2)
scatter(x, y, 100,h,'filled')
ylabel('f_W') 
xlabel('run number') 
%axis([0 38 -4 4])
hold on;
errorbar(x, y, z, 'LineStyle','none');
%yline(0.731,'LineWidth',3);
yline(ymean,'LineWidth',3);

set(gca,'FontSize',25)
grid on
