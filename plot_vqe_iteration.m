figure

y_fci1=-9.608760262340200597e+02;
y_fci2=-9.608126606978445352e+02;

nn1=length(y1);
x1=1:nn1;

nn2=length(y2);
x2=1:nn2;

nn3=length(y3);
x3=1:nn3;

nn4=length(y4);
x4=1:nn4;

subplot(2,2,1)
h1=plot(x1,y1,'o','linewidth',2,'Markersize',12);
col1=get(h1,'color');
set(h1,'MarkerFaceColor',col1);

hold on
h2=plot(x2,y2,'o','linewidth',2,'Markersize',12);
col2=get(h2,'color');
set(h2,'MarkerFaceColor',col2);

nn=max(nn1,nn2);
h3=plot([1 nn],[y_fci1,y_fci1],'g--','linewidth',2);

legend([h1,h2],{'HAA','HEA'})
legend('boxoff')

set(gca,'fontsize',20)
set(gca,'linewidth',2)

xlabel('Iteration','interpreter','latex')
ylabel('Energy (Ha)','interpreter','latex')

subplot(2,2,2)
x_sq = [1,nn,nn,1];   % 逆时针遍历每个点的x值
y_sq = [1e-5,1e-5,1.6e-3,1.6e-3];   % 逆时针遍历每个点的y值
fill(x_sq,y_sq,'p');   % 填充函数
hold on

h4=plot(x1,y1-y_fci1,'o','linewidth',2,'Color',col1,'Markersize',12);
set(h4,'MarkerFaceColor',col1);

h5=plot(x2,y2-y_fci1,'o','linewidth',2,'Color',col2,'Markersize',12);
set(h5,'MarkerFaceColor',col2);

set(gca,'fontsize',20)
set(gca,'linewidth',2)

xlabel('Iteration','interpreter','latex')
ylabel('Absolute energy error(Ha)','interpreter','latex')

subplot(2,2,3)
h6=plot(x3,y3,'o','linewidth',2,'Color',col1,'Markersize',12);
set(h6,'MarkerFaceColor',col1);

hold on
h7=plot(x4,y4,'o','linewidth',2,'Color',col2,'Markersize',12);
set(h7,'MarkerFaceColor',col2);

nn=max(nn3,nn4);
h8=plot([1 nn],[y_fci2,y_fci2],'g--','linewidth',2);

set(gca,'fontsize',20)
set(gca,'linewidth',2)

xlabel('Iteration','interpreter','latex')
ylabel('Energy (Ha)','interpreter','latex')

subplot(2,2,4)
x_sq = [1,nn,nn,1];   % 逆时针遍历每个点的x值
y_sq = [1e-5,1e-5,1.6e-3,1.6e-3];   % 逆时针遍历每个点的y值
fill(x_sq,y_sq,'p');   % 填充函数
hold on

h9=plot(x3,y3-y_fci2,'o','linewidth',2,'Color',col1,'Markersize',12);
set(h9,'MarkerFaceColor',col1);

h10=plot(x4,y4-y_fci2,'o','linewidth',2,'Color',col2,'Markersize',12);
set(h5,'MarkerFaceColor',col2);

set(gca,'fontsize',20)
set(gca,'linewidth',2)

xlabel('Iteration','interpreter','latex')
ylabel('Absolute energy error(Ha)','interpreter','latex')