figure

subplot(2,2,1)

y_fci1=-9.608760262340200597e+02;
y_fci2=-9.608126606978445352e+02;
nn=length(y1);
x1=1:nn;
nn2=length(y2);
x2=1:nn2;
nn3=length(y3);
x3=1:nn3;
nn4=length(y4);
x4=1:nn4;

h1=semilogx(x1,y1,'-o','linewidth',2,'Markersize',12);
col1=get(h1,'color');
set(h1,'MarkerFaceColor',col1);
set(gca,'fontsize',20)
set(gca,'linewidth',2)

hold on
h2=semilogx(x2,y2,'-o','linewidth',2,'Markersize',12);
col2=get(h2,'color');
set(h2,'MarkerFaceColor',col2);

legend('HAA','HEA')
legend('boxoff')

xlabel('Iterations','interpreter','latex')
ylabel('Energy (Ha)','interpreter','latex')

subplot(2,2,2)
x_sq = [1,1e4,1e4,1];   % 逆时针遍历每个点的x值
y_sq = [1e-5,1e-5,1.6e-3,1.6e-3];   % 逆时针遍历每个点的y值
fill(x_sq,y_sq,'p');   % 填充函数
hold on

h3=loglog(x1,abs(y1-y_fci1),'-o','Color',col1,'linewidth',2,'Markersize',12);
set(h3,'MarkerFaceColor',col1);

hold on
h4=loglog(x2,abs(y2-y_fci1),'-o','Color',col2,'linewidth',2,'Markersize',12);
set(h4,'MarkerFaceColor',col2);

set(gca,'fontsize',20)
set(gca,'linewidth',2)

xlabel('Iterations','interpreter','latex')
ylabel('Absoluste Energy error (Ha)','interpreter','latex')

subplot(2,2,3)
h5=semilogx(x3,y3,'-o','Color',col1,'linewidth',2,'Markersize',12);
set(h5,'MarkerFaceColor',col1);

hold on
h6=semilogx(x4,y4,'-o','Color',col2,'linewidth',2,'Markersize',12);
set(h6,'MarkerFaceColor',col2);

set(gca,'fontsize',20)
set(gca,'linewidth',2)

%legend('UHF','CISD','VQE','NNQS','HF+NNQS','CISD+NNQS','VQE+NNQS','FCI')
%legend('boxoff')

%legend('VQE-Simulator','FCI')
%legend('boxoff')
xlabel('Iterations','interpreter','latex')
ylabel('Energy (Ha)','interpreter','latex')

subplot(2,2,4)
x_sq = [1,1e5,1e5,1];   % 逆时针遍历每个点的x值
y_sq = [1e-5,1e-5,1.6e-3,1.6e-3];   % 逆时针遍历每个点的y值
fill(x_sq,y_sq,'p');   % 填充函数
hold on

h7=loglog(x3,abs(y3-y_fci2),'-o','Color',col1,'linewidth',2,'Markersize',12);
set(h7,'MarkerFaceColor',col1);

hold on
h8=loglog(x4,abs(y4-y_fci2),'-o','Color',col2,'linewidth',2,'Markersize',12);
set(h7,'MarkerFaceColor',col2);

set(gca,'fontsize',20)
set(gca,'linewidth',2)

xlabel('Iterations','interpreter','latex')
ylabel('Absoluste Energy error (Ha)','interpreter','latex')