
figure; 
hold on
for i = 1:17
    for j = 1:17
        if confusion_matrix(i,j) ~= 0
            scatter(i,j,confusion_matrix(i,j)*50,'k','filled')
        end
    end
end

set(gca,'fontsize',14)
xticks(1:4:17)
xticklabels(Flist(1:4:17))
xlim([0 18])
ylim([0 18])
yticks(1:4:17)
yticklabels(Flist(1:4:17))

plot([0 18],[0,18],':k','linewidth',2)
