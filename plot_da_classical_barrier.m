figure

X=categorical({'HF', 'CCSD', 'CASCI', 'CASSCF', 'This work'});
X=reordercats(X,{'HF', 'CCSD', 'CASCI', 'CASSCF', 'This work'});

bar(X(1:end-1),yy1(1:end-1))
hold on
bar(X(end),yy1(end))

ax = gca;
for i = 1:length(yy1)
text(i, yy1(i)+1, num2str(yy1(i),'%.2f'), 'HorizontalAlignment', 'center','fontsize',15); % 在每个柱子上显示数值
end

hold on
plot([0 6],[34.5 34.5],'g--','linewidth',2)

set(gca,'fontsize',20)
set(gca,'linewidth',2)

ylabel('Energy (kcal/mol)','interpreter','latex')

