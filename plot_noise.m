figure

subplot(1,2,1)
X=categorical({'HF', 'CCSD', 'CISD', 'CASCI', 'CASSCF', 'This work'});
X=reordercats(X,{'HF', 'CCSD', 'CISD', 'CASCI', 'CASSCF', 'This work'});

bar(X(1:end-1),yy1(1:end-1))
hold on
bar(X(end),yy1(end))

ax = gca;
for i = 1:length(yy1)
text(i, yy1(i)+1, num2str(yy1(i),'%.2f'), 'HorizontalAlignment', 'center','fontsize',15); % 在每个柱子上显示数值
end

hold on
plot([0 9],[34.5 34.5],'g--','linewidth',2)

set(gca,'fontsize',20)
set(gca,'linewidth',2)

ylabel('Energy (kcal/mol)','interpreter','latex')


y_fci1=-9.608760262340200597e+02;
y_fci2=-9.608126606978445352e+02;
subplot(1,2,2)
y_barrier=y_fci2-y_fci1;
yy5=abs(y3-y1-y_barrier)*627.503;
yy6=abs(y4-y2-y_barrier)*627.503;
h1=plot(xx1(1:4:end),yy5(1:4:end),'o','linewidth',2,'Markersize',12);
hold on
h2=plot(xx1(1:4:end),yy6(1:4:end),'o','linewidth',2,'Markersize',12);

legend([h1,h2],{'HAA','HEA'})
legend('boxoff')

set(gca,'fontsize',20)
set(gca,'linewidth',2)

xlabel('Depolarizing error','interpreter','latex')
ylabel('Absolute energy error (kcal/mol)','interpreter','latex')
