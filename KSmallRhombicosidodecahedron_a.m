clear
%% unit cell model definition

L = 1;  r0 = 1e-2;
Es = 7e7; nu = 0.3; Gm = Es/2/(1+nu);
Area = pi*r0^2;
Ixx = r0^4/12; Iyy = r0^4/12; Izz = 2*r0^4/12;

nodes = [-1/2, -1/2, -1 - sqrt(5)/2;-1/2, -1/2, (2 + sqrt(5))/2;-1/2, 1/2, -1 - sqrt(5)/2;-1/2, 1/2, (2 + sqrt(5))/2;-1/2, -1 - sqrt(5)/2, -1/2;-1/2, -1 - sqrt(5)/2, 1/2;-1/2, (2 + sqrt(5))/2, -1/2;-1/2, (2 + sqrt(5))/2, 1/2;0, (-3 + sqrt(5))^(-1), (-5 - sqrt(5))/4;0, (-3 + sqrt(5))^(-1), (5 + sqrt(5))/4;0, (3 + sqrt(5))/4, (-5 - sqrt(5))/4;0, (3 + sqrt(5))/4, (5 + sqrt(5))/4;1/2, -1/2, -1 - sqrt(5)/2;1/2, -1/2, (2 + sqrt(5))/2;1/2, 1/2, -1 - sqrt(5)/2;1/2, 1/2, (2 + sqrt(5))/2;1/2, -1 - sqrt(5)/2, -1/2;1/2, -1 - sqrt(5)/2, 1/2;1/2, (2 + sqrt(5))/2, -1/2;1/2, (2 + sqrt(5))/2, 1/2;(-5 - sqrt(5))/4, 0, (-3 + sqrt(5))^(-1);(-5 - sqrt(5))/4, 0, (3 + sqrt(5))/4;(-1 - sqrt(5))/4, (-1 - sqrt(5))/2, (-3 + sqrt(5))^(-1);(-1 - sqrt(5))/4, (-1 - sqrt(5))/2, (3 + sqrt(5))/4;(-1 - sqrt(5))/4, (1 + sqrt(5))/2, (-3 + sqrt(5))^(-1);(-1 - sqrt(5))/4, (1 + sqrt(5))/2, (3 + sqrt(5))/4;(-1 - sqrt(5))/2, (-3 + sqrt(5))^(-1), (-1 - sqrt(5))/4;(-1 - sqrt(5))/2, (-3 + sqrt(5))^(-1), (1 + sqrt(5))/4;(-1 - sqrt(5))/2, (3 + sqrt(5))/4, (-1 - sqrt(5))/4;(-1 - sqrt(5))/2, (3 + sqrt(5))/4, (1 + sqrt(5))/4;-1 - sqrt(5)/2, -1/2, -1/2;-1 - sqrt(5)/2, -1/2, 1/2;-1 - sqrt(5)/2, 1/2, -1/2;-1 - sqrt(5)/2, 1/2, 1/2;(-3 + sqrt(5))^(-1), (-5 - sqrt(5))/4, 0;(-3 + sqrt(5))^(-1), (-1 - sqrt(5))/4, (-1 - sqrt(5))/2;(-3 + sqrt(5))^(-1), (-1 - sqrt(5))/4, (1 + sqrt(5))/2;(-3 + sqrt(5))^(-1), (1 + sqrt(5))/4, (-1 - sqrt(5))/2;(-3 + sqrt(5))^(-1), (1 + sqrt(5))/4, (1 + sqrt(5))/2;(-3 + sqrt(5))^(-1), (5 + sqrt(5))/4, 0;(1 + sqrt(5))/4, (-1 - sqrt(5))/2, (-3 + sqrt(5))^(-1);(1 + sqrt(5))/4, (-1 - sqrt(5))/2, (3 + sqrt(5))/4;(1 + sqrt(5))/4, (1 + sqrt(5))/2, (-3 + sqrt(5))^(-1);(1 + sqrt(5))/4, (1 + sqrt(5))/2, (3 + sqrt(5))/4;(1 + sqrt(5))/2, (-3 + sqrt(5))^(-1), (-1 - sqrt(5))/4;(1 + sqrt(5))/2, (-3 + sqrt(5))^(-1), (1 + sqrt(5))/4;(1 + sqrt(5))/2, (3 + sqrt(5))/4, (-1 - sqrt(5))/4;(1 + sqrt(5))/2, (3 + sqrt(5))/4, (1 + sqrt(5))/4;(2 + sqrt(5))/2, -1/2, -1/2;(2 + sqrt(5))/2, -1/2, 1/2;(2 + sqrt(5))/2, 1/2, -1/2;(2 + sqrt(5))/2, 1/2, 1/2;(3 + sqrt(5))/4, (-5 - sqrt(5))/4, 0;(3 + sqrt(5))/4, (-1 - sqrt(5))/4, (-1 - sqrt(5))/2;(3 + sqrt(5))/4, (-1 - sqrt(5))/4, (1 + sqrt(5))/2;(3 + sqrt(5))/4, (1 + sqrt(5))/4, (-1 - sqrt(5))/2;(3 + sqrt(5))/4, (1 + sqrt(5))/4, (1 + sqrt(5))/2;(3 + sqrt(5))/4, (5 + sqrt(5))/4, 0;(5 + sqrt(5))/4, 0, (-3 + sqrt(5))^(-1);(5 + sqrt(5))/4, 0, (3 + sqrt(5))/4]*L;
beams = struct('nodes', [1, 3;1, 9;1, 13;1, 36;2, 4;2, 10;2, 14;2, 37;3, 11;3, 15;3, 38;4, 12;4, 16;4, 39;5, 6;5, 17;5, 23;5, 35;6, 18;6, 24;6, 35;7, 8;7, 19;7, 25;7, 40;8, 20;8, 26;8, 40;9, 13;9, 23;9, 41;10, 14;10, 24;10, 42;11, 15;11, 25;11, 43;12, 16;12, 26;12, 44;13, 15;13, 54;14, 16;14, 55;15, 56;16, 57;17, 18;17, 41;17, 53;18, 42;18, 53;19, 20;19, 43;19, 58;20, 44;20, 58;21, 31;21, 33;21, 36;21, 38;22, 32;22, 34;22, 37;22, 39;23, 27;23, 36;24, 28;24, 37;25, 29;25, 38;26, 30;26, 39;27, 31;27, 35;27, 36;28, 32;28, 35;28, 37;29, 33;29, 38;29, 40;30, 34;30, 39;30, 40;31, 32;31, 33;32, 34;33, 34;41, 45;41, 54;42, 46;42, 55;43, 47;43, 56;44, 48;44, 57;45, 49;45, 53;45, 54;46, 50;46, 53;46, 55;47, 51;47, 56;47, 58;48, 52;48, 57;48, 58;49, 50;49, 51;49, 59;50, 52;50, 60;51, 52;51, 59;52, 60;54, 59;55, 60;56, 59;57, 60]);

mat = struct('E', Es, 'nu', nu);
prop.beams = struct('A', Area, ...
    'Ixx', Ixx, 'Iyy', Iyy, 'Izz', Izz);
model = struct('nodes', nodes, 'beams',  beams, ...
    'mat', mat, 'prop', prop);
%% periodic directions
dirs = [13,14; 50,32; 26,53];
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
