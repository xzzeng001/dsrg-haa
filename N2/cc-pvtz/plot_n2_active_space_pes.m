figure

subplot(1,2,1)

% for the AVAS method
h1=plot(x1,ya1,'o','linewidth',2,'Markersize',15);
col1=get(h1,'color');
set(h1,'MarkerFaceColor',col1);
hold on

h2=plot(x1,ya2,'o','color',col1,'linewidth',2,'Markersize',15);
set(h2,'MarkerFaceColor',col1);
h3=plot(x1,ya3,'o','color',col1,'linewidth',2,'Markersize',15);
set(h3,'MarkerFaceColor',col1);
h4=plot(x1,ya4,'o','color',col1,'linewidth',2,'Markersize',15);
set(h4,'MarkerFaceColor',col1);
h5=plot(x1,ya5,'o','color',col1,'linewidth',2,'Markersize',15);
set(h5,'MarkerFaceColor',col1);
h6=plot(x1,ya6,'o','color',col1,'linewidth',2,'Markersize',15);
set(h6,'MarkerFaceColor',col1);

h7=plot(xr1,yr1,'o','color',col1,'linewidth',2,'Markersize',15);
set(h7,'MarkerFaceColor',col1);

% for the entropy-based method
h8=plot(x1,ye1,'s','color',[0.4660 0.6740 0.1880],'linewidth',2,'Markersize',15);
col2=get(h8,'color');
%set(h8,'MarkerFaceColor',col2);

h9=plot(x1,ye2,'s','color',col2,'linewidth',2,'Markersize',15);
%set(h9,'MarkerFaceColor',col2);
h10=plot(x1,ye3,'s','color',col2,'linewidth',2,'Markersize',15);
%set(h10,'MarkerFaceColor',col2);
h11=plot(x1,ye4,'s','color',col2,'linewidth',2,'Markersize',15);
%set(h11,'MarkerFaceColor',col2);

h12=plot(xr2,yr2,'s','color',col2,'linewidth',2,'Markersize',15);
%set(h12,'MarkerFaceColor',col2);
h13=plot(xr3,yr3,'s','color',col2,'linewidth',2,'Markersize',15);
%set(h13,'MarkerFaceColor',col2);
xr4=[0.6,0.8,1.4];
yr4=[8,9,8];
h14=plot(xr4,yr4,'s','color',col2,'linewidth',2,'Markersize',15);
%set(h14,'MarkerFaceColor',col2);
xr5=[0.6,1.4];
yr5=[9,9];
h15=plot(xr5,yr5,'s','color',col2,'linewidth',2,'Markersize',15);
%set(h15,'MarkerFaceColor',col2);

% for the orbital correlation
h16=plot(x1,yc1,'p','color',[0.3010 0.7450 0.9330],'linewidth',2,'Markersize',15);
col3=get(h16,'color');
set(h16,'MarkerFaceColor',col3);

h17=plot(x1,yc2,'p','color',col3,'linewidth',2,'Markersize',15);
set(h17,'MarkerFaceColor',col3);
h18=plot(x1,yc3,'p','color',col3,'linewidth',2,'Markersize',15);
set(h18,'MarkerFaceColor',col3);
h19=plot(x1,yc4,'p','color',col3,'linewidth',2,'Markersize',15);
set(h19,'MarkerFaceColor',col3);

h20=plot(xr6,yr6,'p','color',col3,'linewidth',2,'Markersize',15);
set(h20,'MarkerFaceColor',col3);
xr7=[2.4,2.8,3.2,3.6,4.0,5.0,6.0];
yr7=[9,9,9,9,9,9,9];
h21=plot(xr7,yr7,'p','color',col3,'linewidth',2,'Markersize',15);
set(h21,'MarkerFaceColor',col3);

%h22=plot(0.6,15,'p','color',col3,'linewidth',2,'Markersize',15);
%set(h22,'MarkerFaceColor',col3);

set(gca,'fontsize',20)
set(gca,'linewidth',2)

legend([h1,h8,h16],{'AVAS','autoCAS','This work'})
legend('boxoff')

xlabel('R$_{N-N} (\AA)$','interpreter','latex')
ylabel('Active space index','interpreter','latex')


subplot(1,2,2)

h30=plot(x1,yya1,'-o','color',col1,'linewidth',2,'Markersize',12);
set(h30,'MarkerFaceColor',col1);

hold on
h31=plot(x1,yya2,'-o','color',col1,'linewidth',2,'Markersize',15);
set(h31,'MarkerFaceColor',col1);
h32=plot(x1,yya3,'-+','color',col1,'linewidth',2,'Markersize',12);
set(h32,'MarkerFaceColor',col1);
hold on

h33=plot(x1,yye1,'-s','color',col2,'linewidth',2,'Markersize',12);
set(h33,'MarkerFaceColor',col2);
h34=plot(x1,yye2,'-s','color',col2,'linewidth',2,'Markersize',15);
set(h34,'MarkerFaceColor',col2);
h35=plot(x1,yye3,'-+','color',col2,'linewidth',2,'Markersize',12);
set(h35,'MarkerFaceColor',col2);

h36=plot(x1,yyc1,'-p','color',col3,'linewidth',2,'Markersize',12);
set(h36,'MarkerFaceColor',col3);
h37=plot(x1,yyc2,'-p','color',col3,'linewidth',2,'Markersize',15);
set(h37,'MarkerFaceColor',col3);
h38=plot(x1,yyc3,'-+','color',col3,'linewidth',2,'Markersize',12);
set(h38,'MarkerFaceColor',col3);

legend([h30,h31,h32],{'CASCI','CASSCF','CASSCF-NEVPT2'})
legend('boxoff')

set(gca,'fontsize',20)
set(gca,'linewidth',2)

xlabel('R$_{N-N} (\AA)$','interpreter','latex')
ylabel('Energy (Ha)','interpreter','latex')
