
cn = 500;                                        % Number Of Colors (Scale AsAppropriate)
cm = colormap(jet(cn));                                         % Define Colormap
figure(2)
scatter(Data.Cut, Data.DistMean, 100, 'filled')
hold on;
scatter(Data.Cut, Data.CSWMean , 100, 'filled','d')
hold on;
errorbar(Data.Cut, Data.CSWMean, Data.ErrDistMean, 'LineStyle','none');
%yline(0.6866,'LineWidth',3);
set(gca,'FontSize',25)
grid on
legend({'f_W^{DIST}','f_W^{CSW}'},'Location','best', 'Fontsize',30, 'Orientation','horizontal')
ylabel('f_W') 
xlabel('threshold [mV]') 
