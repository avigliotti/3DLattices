clear
%% unit cell model definition

L = 1;  r0 = 1e-2; t = r0/5;
Es = 7e7; nu = 0.3; Gm = Es/2/(1+nu);
Area = pi*r0^2;
Ixx = r0^4/12; Iyy = r0^4/12; Izz = 2*r0^4/12;

nodes = [-(1/2),1/2+1/sqrt(2),-(1/2)-sqrt(2);-(1/2),1/2+1/sqrt(2),1/2+sqrt(2);-(1/2),1/(2-2*sqrt(2)),-(1/2)-sqrt(2);-(1/2),1/(2-2*sqrt(2)),1/2+sqrt(2);-(1/2),-(1/2)-sqrt(2),1/2+1/sqrt(2);-(1/2),-(1/2)-sqrt(2),1/(2-2*sqrt(2));-(1/2),1/2+sqrt(2),1/2+1/sqrt(2);-(1/2),1/2+sqrt(2),1/(2-2*sqrt(2));1/2,1/2+1/sqrt(2),-(1/2)-sqrt(2);1/2,1/2+1/sqrt(2),1/2+sqrt(2);1/2,1/(2-2*sqrt(2)),-(1/2)-sqrt(2);1/2,1/(2-2*sqrt(2)),1/2+sqrt(2);1/2,-(1/2)-sqrt(2),1/2+1/sqrt(2);1/2,-(1/2)-sqrt(2),1/(2-2*sqrt(2));1/2,1/2+sqrt(2),1/2+1/sqrt(2);1/2,1/2+sqrt(2),1/(2-2*sqrt(2));1/2+1/sqrt(2),-(1/2),-(1/2)-sqrt(2);1/2+1/sqrt(2),-(1/2),1/2+sqrt(2);1/2+1/sqrt(2),1/2,-(1/2)-sqrt(2);1/2+1/sqrt(2),1/2,1/2+sqrt(2);1/2+1/sqrt(2),-(1/2)-sqrt(2),-(1/2);1/2+1/sqrt(2),-(1/2)-sqrt(2),1/2;1/2+1/sqrt(2),1/2+sqrt(2),-(1/2);1/2+1/sqrt(2),1/2+sqrt(2),1/2;1/(2-2*sqrt(2)),-(1/2),-(1/2)-sqrt(2);1/(2-2*sqrt(2)),-(1/2),1/2+sqrt(2);1/(2-2*sqrt(2)),1/2,-(1/2)-sqrt(2);1/(2-2*sqrt(2)),1/2,1/2+sqrt(2);1/(2-2*sqrt(2)),-(1/2)-sqrt(2),-(1/2);1/(2-2*sqrt(2)),-(1/2)-sqrt(2),1/2;1/(2-2*sqrt(2)),1/2+sqrt(2),-(1/2);1/(2-2*sqrt(2)),1/2+sqrt(2),1/2;-(1/2)-sqrt(2),-(1/2),1/2+1/sqrt(2);-(1/2)-sqrt(2),-(1/2),1/(2-2*sqrt(2));-(1/2)-sqrt(2),1/2,1/2+1/sqrt(2);-(1/2)-sqrt(2),1/2,1/(2-2*sqrt(2));-(1/2)-sqrt(2),1/2+1/sqrt(2),-(1/2);-(1/2)-sqrt(2),1/2+1/sqrt(2),1/2;-(1/2)-sqrt(2),1/(2-2*sqrt(2)),-(1/2);-(1/2)-sqrt(2),1/(2-2*sqrt(2)),1/2;1/2+sqrt(2),-(1/2),1/2+1/sqrt(2);1/2+sqrt(2),-(1/2),1/(2-2*sqrt(2));1/2+sqrt(2),1/2,1/2+1/sqrt(2);1/2+sqrt(2),1/2,1/(2-2*sqrt(2));1/2+sqrt(2),1/2+1/sqrt(2),-(1/2);1/2+sqrt(2),1/2+1/sqrt(2),1/2;1/2+sqrt(2),1/(2-2*sqrt(2)),-(1/2);1/2+sqrt(2),1/(2-2*sqrt(2)),1/2]*L;
beams = struct('nodes', [1, 8;1, 9;1, 27;2, 7;2, 10;2, 28;3, 6;3, 11;3, 25;4, 5;4, 12;4, 26;5, 13;5, 30;6, 14;6, 29;7, 15;7, 32;8, 16;8, 31;9, 16;9, 19;10, 15;10, 20;11, 14;11, 17;12, 13;12, 18;13, 22;14, 21;15, 24;16, 23;17, 19;17, 42;18, 20;18, 41;19, 44;20, 43;21, 22;21, 47;22, 48;23, 24;23, 45;24, 46;25, 27;25, 34;26, 28;26, 33;27, 36;28, 35;29, 30;29, 39;30, 40;31, 32;31, 37;32, 38;33, 35;33, 40;34, 36;34, 39;35, 38;36, 37;37, 38;39, 40;41, 43;41, 48;42, 44;42, 47;43, 46;44, 45;45, 46;47, 48]);

mat = struct('E', Es, 'nu', nu);
prop.beams = struct('A', Area, ...
    'Ixx', Ixx, 'Iyy', Iyy, 'Izz', Izz);
model = struct('nodes', nodes, 'beams',  beams, ...
    'mat', mat, 'prop', prop);
%% periodic directions
dirs = [10,14; 18,34; 12,16];
a1 = (model.nodes(dirs(1,2),:)-model.nodes(dirs(1,1), :))';
a2 = (model.nodes(dirs(2,2),:)-model.nodes(dirs(2,1), :))';
a3 = (model.nodes(dirs(3,2),:)-model.nodes(dirs(3,1), :))';
%% plot unit cell
figure(4); clf

plot3DModel(model, [],...
    struct('nodeslabels', false, 'FontSize', 18, ...
    'plotundef', true, 'undeflinewidth', 1, 'undeflinestyle', '-', ...
    'undeflinecolor', 'b'));

plot3DModel(model, [],...
    struct('nodeslabels', true, 'FontSize', 18, ...
    'plotundef', true, 'undeflinewidth', 6, 'undeflinestyle', '-'));
xlabel('x'); ylabel('y'); zlabel('z');

daspect([1, 1, 1]); view(3)
light('Position',[1 1 1],'Style','infinite');
set(gca, 'Visible', 'off')
myaxis = axis;
shg
p0 = -[1,1,1]/10;
cg = mean(nodes);
set(quiver3(cg(1), cg(2), cg(3), a1(1), a1(2), a1(3), 1), ...
    'linewidth', 2, 'color', 'b');
set(text(cg(1)+a1(1), cg(2)+a1(2), cg(3)+a1(3), 'a$_1$'), ...
    'fontSize', 18, 'interpreter', 'latex');

set(quiver3(cg(1), cg(2), cg(3), a2(1), a2(2), a2(3), 1), ...
    'linewidth', 2, 'color', 'b');
set(text(cg(1)+a2(1), cg(2)+a2(2), cg(3)+a2(3), 'a$_2$'), ...
    'fontSize', 18, 'interpreter', 'latex');
set(quiver3(cg(1), cg(2), cg(3), a3(1), a3(2), a3(3), 1), ...
    'linewidth', 2, 'color', 'b');
set(text(cg(1)+a3(1), cg(2)+a3(2), cg(3)+a3(3), 'a$_3$'), ...
    'fontSize', 18, 'interpreter', 'latex');

camlight; lighting phong; material metal;
alpha (0.25);
%% find material properties
[Keps, Vol, Vol0, deps] = Find3DMatProp(model, a1, a2, a3);

fprintf ('the material volume of the unit cell is : %.4f\n', Vol);
fprintf ('the volume of the unit cell is          : %.4f\n', Vol0);
fprintf ('the relative density of the lattice is  : %.4f\n', Vol/Vol0);
fprintf ('the lattice stiffnes matrix is : \n');
disp (Keps);
%% plot pure deformation modes of the unit cell

nnodes = size(model.nodes,1);
stitles = {{'$\epsilon_{11} = 1, \epsilon_{22} = 0, \epsilon_{33} = 0$', ...
    '$\gamma_{12} = 0, \gamma_{23} = 0, \gamma_{31} = 0$'}, ...
    {'$\epsilon_{11} = 0, \epsilon_{22} = 1, \epsilon_{33} = 0$', ...
    '$\gamma_{12} = 0, \gamma_{23} = 0, \gamma_{31} = 0$'}, ...
    {'$\epsilon_{11} = 0, \epsilon_{22} = 0, \epsilon_{33} = 1$', ...
    '$\gamma_{12} = 0, \gamma_{23} = 0, \gamma_{31} = 0$'}, ...
    {'$\epsilon_{11} = 0, \epsilon_{22} = 0, \epsilon_{33} = 0$', ...
    '$\gamma_{12} = 0, \gamma_{23} = 1, \gamma_{31} = 0$'}, ...
    {'$\epsilon_{11} = 0, \epsilon_{22} = 0, \epsilon_{33} = 0$', ...
    '$\gamma_{12} = 0, \gamma_{23} = 0, \gamma_{31} = 1$'}, ...
    {'$\epsilon_{11} = 0, \epsilon_{22} = 0, \epsilon_{33} = 0$', ...
    '$\gamma_{12} = 1, \gamma_{23} = 0, \gamma_{31} = 0$'}};

figure(2); clf
set(gcf, 'position', [20, 400, 1200, 800]);
for nMode = 1:6
    subplot(2, 3, nMode); set(gca, 'fontSize', 14);
    
    disp0 = reshape(deps(:,nMode), 6, nnodes)';
    plot3DModel(model, disp0,...
        struct('nodeslabels', false, 'FontSize', 18, ...
        'plotundef', true, 'undeflinewidth', 2, 'undeflinestyle', '-', ...
        'plotdef', true, 'linewidth', 4, 'maxdisp', 0.3,...
        'drawplates', false, 'NPoints', 20, ...
        'drawfaces', true, 'facecolor', 'g'));
    set(title(stitles{nMode}), 'interpreter', 'latex');
    set(xlabel('x$_1$'), 'interpreter', 'latex');
    set(ylabel('x$_2$'), 'interpreter', 'latex');
    set(zlabel('x$_3$'), 'interpreter', 'latex');
    set(gca, 'xtick', [], 'ytick', [], 'ztick', []);
    view (-13.5, 12);
end
shg

