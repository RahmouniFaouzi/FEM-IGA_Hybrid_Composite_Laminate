function Visualisation_STR_DISP(stress, disp, node, ply, elementV, stressComp)

PlyIn = ['Ply : 0',num2str(ply)];
dir = {'Stress in X direction',
       'Stress in Y direction',
       'Stress in XY direction'};
for i = 1 : stressComp
    figure
    clf
    plot_field(node,elementV,'Q4',stress(:,:,i,ply));
    hold on
    colorbar
    title({dir{i} ; PlyIn})
    axis off
end
hold off
%
figure
clf
plot_field(node,elementV,'Q4',disp(:,:,1));
hold on
colorbar
title({'Displacement in x direction';PlyIn})
axis off
hold off
%
figure
clf
plot_field(node,elementV,'Q4',disp(:,:,2));
hold on
colorbar
title({'Displacement in y direction';PlyIn})
axis off
hold off

% opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
% exportfig(gcf,fileName,opts)
end