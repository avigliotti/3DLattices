clear
%% unit cell model definition
L = 1;  r0 = 1e-2; t = r0/10;
Es = 7e7; nu = 0.3; Gm = Es/2/(1+nu);
Area = pi*r0^2;
Ixx = r0^4/12; Iyy = r0^4/12; Izz = 2*r0^4/12;

nodes= [-(3/2), -(1/2), 0; -(3/2), 1/2, 0; -1, -1, -(1/sqrt(2));...
    -1, -1, 1/sqrt(2); -1, 1, -(1/sqrt(2)); -1, 1, 1/sqrt(2);...
    -(1/2), -(3/2), 0; -(1/2), -(1/2), -sqrt(2);...
    -(1/2), -(1/2), sqrt(2); -(1/2), 1/2, -sqrt(2);...
    -(1/2), 1/2, sqrt(2);-(1/2), 3/2, 0;...
    1/2, -(3/2), 0;1/2, -(1/2), -sqrt(2);...
    1/2, -(1/2), sqrt(2);1/2, 1/2, -sqrt(2);...
    1/2, 1/2, sqrt(2);1/2, 3/2, 0;1, -1, -(1/sqrt(2));...
    1, -1, 1/sqrt(2);1, 1, -(1/sqrt(2));...
    1, 1, 1/sqrt(2);3/2, -(1/2), 0;3/2, 1/2, 0]*L;

beams = struct('nodes', [1, 2; 1, 3; 1, 4; 2, 5; 2, 6; 3, 7;...
    3, 8; 4, 7; 4, 9; 5, 10; 5, 12; 6, 11; 6, 12; 7, 13;...
    8, 10; 8, 14; 9, 11; 9, 15; 10, 16; 11, 17; 12, 18;...
    13, 19; 13, 20; 14, 16; 14, 19; 15, 17; 15, 20; 16, 21;...
    17, 22; 18, 21; 18, 22; 19, 23; 20, 23; 21, 24; 22, 24; 23, 24]);

e1 = (nodes(20,:)-nodes(06,:)); e1=e1/norm(e1);
e2 = (nodes(22,:)-nodes(04,:)); e2=e2/norm(e2);
e3 = cross(e1,e2);
R = [e1;e2;e3];
nodes = nodes*R';

mat = struct('E', Es, 'nu', nu);
prop.beams = struct('A', Area, ...
    'Ixx', Ixx, 'Iyy', Iyy, 'Izz', Izz);
model = struct('nodes', nodes, 'beams',  beams, ...
    'mat', mat, 'prop', prop);
%% periodic directions
dirs = [19,22; 20,21; 05,22];
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

