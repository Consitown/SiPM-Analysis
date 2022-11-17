x=FinalEfficiencies.x;
y=FinalEfficiencies.y;
z=FinalEfficiencies.z3;

[xi,yi] = meshgrid(-400:1:400, -600:1:600);
zi = griddata(x,y,z,xi,yi,'cubic');
h=surf(xi,yi,zi);

c=colorbar
set( c, 'YTick', 90:0.2:100 )
c.FontSize=20;
set(get(c,'Title'),'String','efficiency[%]')
view(0,90)
ylim([-600 600])
set(h,'LineStyle','none')
grid on
daspect([1 1 1])



%EXTRAPOLATION

%figure

%F = scatteredInterpolant(x,y,z,"linear","linear");
%zi = F(xi,yi);
%zi(zi > 100.00) = 100.00;
%surf(xi,yi,zi, 'EdgeAlpha', 0)
%c=colorbar
%set( c, 'YTick', 80:0.2:100 )
%c.FontSize=20;
%set(get(c,'Title'),'String','efficiency[%]')
%view(0,90)
%grid on
%daspect([1 1 1])


%figure

%F = scatteredInterpolant(x,y,z,"linear","linear");
%zi = F(xi,yi);
%[M,c] = contour3(xi,yi,zi,20);
%c.LineWidth = 5;
%grid on
%daspect([1 1 1])