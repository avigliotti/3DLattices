
%% primitive unit, geometry and mesh

L0 = 10;  thk = L0/20;
Es = 1e4; nu = 0.3;

pv=[0, 0; 1, 0;
    1, 1; 0, 1;
    0, 0]*L0;

figure(1); clf;
[pcoords, tnodes] = distmesh2d(@dpoly, @huniform, L0/9, ...
    [0, 0; L0, L0], pv, pv);
%% boundary conditions
nnodes = size(pcoords, 1);
dTol = 1e-6;
pcoords = fix(pcoords/dTol)*dTol;

xmin = min(pcoords(:,1));
xmax = max(pcoords(:,1));

id_leftside = find(pcoords(:,1) == xmin);
id_rightside = find(pcoords(:,1) == xmax);

Disp  	= nan(nnodes, 6);
Force   = zeros(nnodes, 6);

Disp(:, 6)           = 0;       % remove oop rotation DoF
Disp(id_leftside, :) = 0;       % left side is clamped

Force(id_rightside, 3) = -1;    % a vertical force is applied to the right side

ifree = reshape(isnan(Disp.'), ...
    nnodes*6, 1);               % indexes of free DoFs
icnst = ~ifree;
%% model

nodes = [pcoords, zeros(size(pcoords, 1), 1)];
tria = struct('nodes', tnodes);

mat = struct('E', Es, 'nu', nu);
prop.tria = struct('t', thk);
model = struct('nodes', nodes, 'tria',  tria, ...
    'mat', mat, 'prop', prop);
%% find stiffeness matrix and solve

[K, Vol] = Solve3DModelLin(model);

u = reshape(Disp.', nnodes*6, 1);
f = reshape(Force.', nnodes*6, 1);

u(ifree) = K(ifree, ifree)\f(ifree) - K(ifree, icnst)*u(icnst);
f(icnst) = K(icnst, ifree)*u(ifree) + K(icnst, icnst)*u(icnst);

Disp = reshape(u, 6, nnodes).';
Force = reshape(f, 6, nnodes).';
%% plot deformed structure
figure(2); clf;
set(gcf, 'position', [200, 100, 700, 850]);

subplot(2,1,1);
plot3DModel(model, Disp,...
    struct('nodeslabels', false, 'FontSize', 18, ...
    'plotundef', true, 'undeflinewidth', 2, 'undeflinestyle', '-', ...
    'plotdef', false, 'linewidth', 4, 'maxdisp', L0/5,...
    'drawplates', true, 'NPoints', 20, ...
    'drawfaces', true, 'facecolor', 'g'));
set(xlabel('x$_1$'), 'interpreter', 'latex');
set(ylabel('x$_2$'), 'interpreter', 'latex');
set(zlabel('x$_3$'), 'interpreter', 'latex');
axis tight;
view(3);

subplot(2,1,2);
plot3DModel(model, Disp,...
    struct('nodeslabels', false, 'FontSize', 18, ...
    'plotundef', false, 'undeflinewidth', 2, 'undeflinestyle', '-', ...
    'plotdef', true, 'linewidth', 4, 'maxdisp', L0/5,...
    'drawplates', true, 'NPoints', 20, ...
    'drawfaces', true, 'facecolor', 'g'));
set(xlabel('x$_1$'), 'interpreter', 'latex');
set(ylabel('x$_2$'), 'interpreter', 'latex');
set(zlabel('x$_3$'), 'interpreter', 'latex');
axis tight;
view(3);

disp(mean(abs(Disp(id_rightside, 3))));
