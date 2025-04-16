figure

subplot(1,2,1)
y_fci1=-9.608760262340200597e+02;
y_fci2=-9.608218870699513445e+02;

h1=plot(1-x1,y3-y1,'o','linewidth',2,'Markersize',12);
col1=get(h1,'color');
set(h1,'MarkerFaceColor',col1);
set(gca,'fontsize',20)
set(gca,'linewidth',2)

hold on
h2=plot(1-x2,y4-y2,'o','linewidth',2,'Markersize',12);
col2=get(h2,'color');
set(h2,'MarkerFaceColor',col2);

h3=plot(1-x1,yy3-yy1,'o','linewidth',2,'Markersize',12);
col3=get(h3,'color');
set(h3,'MarkerFaceColor',col3);

h4=plot(1-x2,yy4-yy2,'o','linewidth',2,'Markersize',12);
col4=get(h4,'color');
set(h4,'MarkerFaceColor',col4);

h5=plot([1e-3 0.1],[y_fci2-y_fci1 y_fci2-y_fci1],'g--','linewidth',2,'Markersize',12);

legend('Barrier-HAA','Barrier-HEA','Barrier-HAA-ZNE','Barrier-HEA-ZNE','Barrier-FCI')
legend('boxoff')

xlabel('Depolarizing error','interpreter','latex')
ylabel('Energy Barrier (Ha)','interpreter','latex')

subplot(1,2,2)
x_sq = [1e-3,0.1,0.1,1e-3];   % 逆时针遍历每个点的x值
y_sq = [1e-5,1e-5,1.6e-3,1.6e-3];   % 逆时针遍历每个点的y值
fill(x_sq,y_sq,'p');   % 填充函数
hold on

h6=loglog(1-x1,abs((y3-y1)-(y_fci2-y_fci1)),'o','Color',col1,'linewidth',2,'Markersize',12);
set(h6,'MarkerFaceColor',col1);

hold on
h7=loglog(1-x2,abs((y4-y2)-(y_fci2-y_fci1)),'o','Color',col2,'linewidth',2,'Markersize',12);
set(h7,'MarkerFaceColor',col2);

h8=loglog(1-x1,abs((yy3-yy1)-(y_fci2-y_fci1)),'o','Color',col3,'linewidth',2,'Markersize',12);
set(h8,'MarkerFaceColor',col3);

h9=loglog(1-x2,abs((yy4-yy2)-(y_fci2-y_fci1)),'o','Color',col4,'linewidth',2,'Markersize',12);
set(h9,'MarkerFaceColor',col4);

set(gca,'fontsize',20)
set(gca,'linewidth',2)

xlabel('Depolarizing error','interpreter','latex')
ylabel('Absoluste Energy error (Ha)','interpreter','latex')
