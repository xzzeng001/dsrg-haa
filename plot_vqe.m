figure

subplot(2,2,1)

h1=plot(x1,y1,'-o','linewidth',2,'Markersize',12);
col1=get(h1,'color');
set(h1,'MarkerFaceColor',col1);
set(gca,'fontsize',20)
set(gca,'linewidth',2)

%legend('UHF','CISD','VQE','NNQS','HF+NNQS','CISD+NNQS','VQE+NNQS','FCI')
%legend('boxoff')

%legend('VQE-Simulator','FCI')
%legend('boxoff')
xlabel('Iterations','interpreter','latex')
ylabel('Energy (Ha)','interpreter','latex')

subplot(2,2,2)
h2=plot(x1,abs(z1),'-o','linewidth',2,'Markersize',12);
col2=get(h2,'color');
set(h2,'MarkerFaceColor',col2);

set(gca,'fontsize',20)
set(gca,'linewidth',2)

xlabel('Iterations','interpreter','latex')
ylabel('Absoluste Energy error (Ha)','interpreter','latex')

subplot(2,2,3)

h3=plot(x2,y2,'-o','linewidth',2,'Markersize',12);
col3=get(h3,'color');
set(h3,'MarkerFaceColor',col3);
set(gca,'fontsize',20)
set(gca,'linewidth',2)

%legend('UHF','CISD','VQE','NNQS','HF+NNQS','CISD+NNQS','VQE+NNQS','FCI')
%legend('boxoff')

%legend('VQE-Simulator','FCI')
%legend('boxoff')
xlabel('Iterations','interpreter','latex')
ylabel('Energy (Ha)','interpreter','latex')

subplot(2,2,4)
h4=plot(x2,abs(z2),'-o','linewidth',2,'Markersize',12);
col4=get(h4,'color');
set(h4,'MarkerFaceColor',col4);
set(gca,'fontsize',20)
set(gca,'linewidth',2)

xlabel('Iterations','interpreter','latex')
ylabel('Absoluste Energy error (Ha)','interpreter','latex')