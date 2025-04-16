figure

subplot(2,2,1)

y_fci1=-9.608760262340200597e+02;
y_fci2=-9.608218870699513445e+02;

h1=plot(1-x1,y1,'o','linewidth',2,'Markersize',12);
col1=get(h1,'color');
set(h1,'MarkerFaceColor',col1);
set(gca,'fontsize',20)
set(gca,'linewidth',2)

hold on
h2=plot(1-x2,y2,'o','linewidth',2,'Markersize',12);
col2=get(h2,'color');
set(h2,'MarkerFaceColor',col2);

h3=plot(1-x1,yy1,'o','linewidth',2,'Markersize',12);
col3=get(h3,'color');
set(h3,'MarkerFaceColor',col3);

h4=plot(1-x1,yy2,'o','linewidth',2,'Markersize',12);
col4=get(h4,'color');
set(h4,'MarkerFaceColor',col4);

legend('HAA','HEA','HAA-ZNE','HEA-ZNE')
legend('boxoff')

xlabel('Depolarizing error','interpreter','latex')
ylabel('Energy (Ha)','interpreter','latex')

subplot(2,2,2)
x_sq = [1e-3,0.1,0.1,1e-3];   % 逆时针遍历每个点的x值
y_sq = [1e-5,1e-5,1.6e-3,1.6e-3];   % 逆时针遍历每个点的y值
fill(x_sq,y_sq,'p');   % 填充函数
hold on

h5=loglog(1-x1,abs(y1-y_fci1),'o','Color',col1,'linewidth',2,'Markersize',12);
set(h5,'MarkerFaceColor',col1);

hold on
h6=loglog(1-x2,abs(y2-y_fci1),'o','Color',col2,'linewidth',2,'Markersize',12);
set(h6,'MarkerFaceColor',col2);

h7=loglog(1-x1,abs(yy1-y_fci1),'o','Color',col3,'linewidth',2,'Markersize',12);
set(h7,'MarkerFaceColor',col3);

h8=loglog(1-x2,abs(yy2-y_fci1),'o','Color',col4,'linewidth',2,'Markersize',12);
set(h8,'MarkerFaceColor',col4);

set(gca,'fontsize',20)
set(gca,'linewidth',2)

xlabel('Depolarizing error','interpreter','latex')
ylabel('Absoluste Energy error (Ha)','interpreter','latex')

subplot(2,2,3)
h9=plot(1-x3,y3,'o','Color',col1,'linewidth',2,'Markersize',12);
set(h9,'MarkerFaceColor',col1);

hold on
h10=plot(1-x4,y4,'o','Color',col2,'linewidth',2,'Markersize',12);
set(h10,'MarkerFaceColor',col2);

h11=plot(1-x3,yy3,'o','Color',col3,'linewidth',2,'Markersize',12);
set(h11,'MarkerFaceColor',col3);

h12=plot(1-x4,yy4,'o','Color',col4,'linewidth',2,'Markersize',12);
set(h12,'MarkerFaceColor',col4);

set(gca,'fontsize',20)
set(gca,'linewidth',2)

xlabel('Depolarizing error','interpreter','latex')
ylabel('Energy (Ha)','interpreter','latex')

subplot(2,2,4)
x_sq = [1e-3,0.1,0.1,1e-3];   % 逆时针遍历每个点的x值
y_sq = [1e-5,1e-5,1.6e-3,1.6e-3];   % 逆时针遍历每个点的y值
fill(x_sq,y_sq,'p');   % 填充函数
hold on

h13=loglog(1-x3,abs(y3-y_fci2),'o','Color',col1,'linewidth',2,'Markersize',12);
set(h13,'MarkerFaceColor',col1);

hold on
h14=loglog(1-x4,abs(y4-y_fci2),'o','Color',col2,'linewidth',2,'Markersize',12);
set(h14,'MarkerFaceColor',col2);

h15=loglog(1-x3,abs(yy3-y_fci2),'o','Color',col3,'linewidth',2,'Markersize',12);
set(h15,'MarkerFaceColor',col3);

h16=loglog(1-x4,abs(yy4-y_fci2),'o','Color',col4,'linewidth',2,'Markersize',12);
set(h16,'MarkerFaceColor',col4);

set(gca,'fontsize',20)
set(gca,'linewidth',2)

xlabel('Depolarizing error','interpreter','latex')
ylabel('Absoluste Energy error (Ha)','interpreter','latex')