figure

y_fci1=-9.608760262340200597e+02;
y_fci2=-9.608126606978445352e+02;

subplot(1,3,1)
xx1=1-x1;
yy1=(y1-y_fci1)*627.503;
yy2=(y2-y_fci1)*627.503;
h1=plot(xx1(1:4:end),yy1(1:4:end),'o','linewidth',2,'Markersize',12);
hold on
h2=plot(xx1(1:4:end),yy2(1:4:end),'o','linewidth',2,'Markersize',12);

legend([h1,h2],{'HAA','HEA'})
legend('boxoff')

set(gca,'fontsize',20)
set(gca,'linewidth',2)

xlabel('Depolarizing error','interpreter','latex')
ylabel('Absolute energy error (kcal/mol)','interpreter','latex')

subplot(1,3,2)
yy3=(y3-y_fci2)*627.503;
yy4=(y4-y_fci2)*627.503;
h3=plot(xx1(1:4:end),yy3(1:4:end),'o','linewidth',2,'Markersize',12);
hold on
h4=plot(xx1(1:4:end),yy4(1:4:end),'o','linewidth',2,'Markersize',12);

set(gca,'fontsize',20)
set(gca,'linewidth',2)

xlabel('Depolarizing error','interpreter','latex')
ylabel('Absolute energy error (kcal/mol)','interpreter','latex')

subplot(1,3,3)
y_barrier=y_fci2-y_fci1;
yy5=abs(y3-y1-y_barrier)*627.503;
yy6=abs(y4-y2-y_barrier)*627.503;
%yy5=abs(yy1+yy3);
%yy6=abs(yy2+yy4);
h5=plot(xx1(1:4:end),yy5(1:4:end),'o','linewidth',2,'Markersize',12);
hold on
h6=plot(xx1(1:4:end),yy6(1:4:end),'o','linewidth',2,'Markersize',12);

set(gca,'fontsize',20)
set(gca,'linewidth',2)

xlabel('Depolarizing error','interpreter','latex')
ylabel('Absolute energy error (kcal/mol)','interpreter','latex')
