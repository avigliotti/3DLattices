% plot3DModel(model, [Disp], [opt])
%
function plot3DModel(varargin)

%% preamble and defaults

if nargin ==0, error('empty argument list'); end
model = varargin{1};

Disp = [];
opt = [];

if nargin >=2, Disp = varargin{2}; end
if isempty(Disp), Disp=zeros(size(model.nodes)); end;
if nargin >=3, opt = varargin{3}; end
if ~isfield(opt, 'scalefactor'),    opt.scalefactor = 1; end
if ~isfield(opt, 'plotundef'),      opt.plotundef = true; end
if ~isfield(opt, 'plotdef'),        opt.plotdef = false; end
if ~isfield(opt, 'linewidth'),      opt.linewidth = 1; end
if ~isfield(opt, 'linestyle'),      opt.linestyle = '-'; end
if ~isfield(opt, 'linecolor'),      opt.linecolor = 'b'; end
if ~isfield(opt, 'undeflinestyle'), opt.undeflinestyle = '-'; end
if ~isfield(opt, 'undeflinewidth'), opt.undeflinewidth = 1; end
if ~isfield(opt, 'undeflinecolor'), opt.undeflinecolor = 'r'; end
if ~isfield(opt, 'NPoints'),        opt.NPoints = 2; end
if ~isfield(opt, 'marker'),         opt.marker = '.'; end
if ~isfield(opt, 'markersize'),     opt.markersize = 3; end
if ~isfield(opt, 'FontSize'),       opt.FontSize = 9; end
if ~isfield(opt, 'nodeslabels'),    opt.nodeslabels = false; end
if ~isfield(opt, 'barslabels'),     opt.barslabels = false; end
if ~isfield(opt, 'drawforces'),     opt.drawforces = false; end
if ~isfield(opt, 'drawfaces'),      opt.drawfaces = true; end
if ~isfield(opt, 'drawbeams'),      opt.drawbeams = true; end
if ~isfield(opt, 'drawplates'),     opt.drawplates = true; end
if ~isfield(opt, 'FaceAlpha'),      opt.FaceAlpha = 1; end
if ~isfield(opt, 'FaceLineStyle'),  opt.FaceLineStyle = '-'; end
if ~isfield(opt, 'platesedgecolor'), opt.platesedgecolor = 'k'; end
if ~isfield(opt, 'platesedgewidth'), opt.platesedgewidth = 0.5; end
if ~isfield(opt, 'plateslinestyle'), opt.plateslinestyle = '-'; end

if isfield(opt, 'maxdisp')
    if ~all(all(Disp==0))
        opt.scalefactor = opt.maxdisp/max(max(abs(Disp(:, 1:3))));
    end
end

Disp = Disp*opt.scalefactor;
washold = ishold;
% if (washold==0), clf; end
hold on;
%% plot tria
if isfield(model, 'tria') && (opt.drawplates),
    nPlates = size(model.tria.nodes, 1);
    if isfield(opt, 'platescolor')
        opt.platescolor = opt.platescolor(:);
        if ischar(opt.platescolor),
            for kk=1:nPlates,
                ElNodes = model.tria.nodes(kk,:);
                set(patch(...
                    model.nodes(ElNodes,1) + Disp(ElNodes,1),...
                    model.nodes(ElNodes,2) + Disp(ElNodes,2),...
                    model.nodes(ElNodes,3) + Disp(ElNodes,3),...
                    opt.platescolor), ...
                    'LineStyle', '-');
            end
        else
            colors = (opt.platescolor-min(opt.platescolor))/...
                (max(opt.platescolor)-min(opt.platescolor))*0.8+0.1;
            colors = [colors, colors, colors];
            for kk=1:nPlates,
                ElNodes = model.tria.nodes(kk,:);
                patch(...
                    model.nodes(ElNodes,1) + Disp(ElNodes,1),...
                    model.nodes(ElNodes,2) + Disp(ElNodes,2),...
                    model.nodes(ElNodes,3) + Disp(ElNodes,3),...
                    colors(kk));
            end
        end
    else
        if opt.plotdef
            platescolor = sqrt(Disp(:,1).^2 + ...
                Disp(:,2).^2 + ...
                Disp(:,3).^2);
            set(patch('Faces', model.tria.nodes, ...
                'Vertices', model.nodes+Disp(:, 1:3)), ...
                'FaceColor', 'interp', ...
                'CDataMapping', 'scaled', ...
                'LineStyle', opt.plateslinestyle, ...
                'EdgeColor', opt.platesedgecolor, ...
                'FaceAlpha', opt.FaceAlpha, ...
                'FaceVertexCData', platescolor, ...
                'LineWidth', opt.platesedgewidth);
        end
        if opt.plotundef
            set(patch('Faces', model.tria.nodes, ...
                'Vertices', model.nodes), ...
                'FaceColor', 'c', ...
                'CDataMapping', 'scaled', ...
                'LineStyle', opt.plateslinestyle, ...
                'EdgeColor', opt.platesedgecolor, ...
                'FaceAlpha', opt.FaceAlpha, ...
                'LineWidth', opt.platesedgewidth);
        end        
    end
end
%% plot faces
if isfield(model, 'faces') && (opt.drawfaces),
    nfaces = length(model.faces.nodes);
    % opt.color = 'c';
    for kk=1:nfaces,
        face=model.faces.nodes{kk};
        if isfield(opt, 'facecolor')
            set(patch(model.nodes(face, 1)+Disp(face,1), ...
                model.nodes(face, 2)+Disp(face,2), ...
                model.nodes(face, 3)+Disp(face,3), ...
                opt.facecolor), ...
                'LineStyle', opt.FaceLineStyle, 'FaceAlpha', opt.FaceAlpha);
        else
            set(patch(model.nodes(face, 1)+Disp(face,1), ...
                model.nodes(face, 2)+Disp(face,2), ...
                model.nodes(face, 3)+Disp(face,3), ...
                sqrt(Disp(face,1).^2+Disp(face,2).^2)), ...
                'LineStyle', opt.FaceLineStyle, 'FaceAlpha', opt.FaceAlpha);
        end
    end
end

%% plot beams
if isfield(model, 'beams') && (opt.drawbeams),
    nBeams = size(model.beams.nodes, 1);
    csi = linspace(-1, 1, opt.NPoints);
    rl = zeros(3, opt.NPoints);
    sl = zeros(3, opt.NPoints);
    if opt.NPoints >2,
        for kk=1:nBeams,
            p0 = model.nodes(model.beams.nodes(kk,:),:);
            s0 = Disp(model.beams.nodes(kk,:),:);
            [T0, L0] = calcRotMatrix(p0);
            T = blkdiag(T0, T0);
            s0l = (T*s0')';
            r0l = [0, L0; 0, 0; 0, 0];
            for ii = 1:opt.NPoints,
                [Nu, Nv, Nw] = shapefunc(csi(ii), L0);
                rl(:, ii) = r0l*Nu;
                sl(:, ii) = [...
                    [s0l(1,1), s0l(2,1)]*Nu;...
                    [s0l(1,[2,6]), s0l(2,[2,6])]*Nv;...
                    [s0l(1,[3,5]), s0l(2,[3,5])]*Nw];
            end
            r = (T0')*rl;
            s = (T0')*sl;
            if (opt.plotundef),
                set(plot3(p0(1,1)+r(1,:),...
                    p0(1,2)+r(2,:),...
                    p0(1,3)+r(3,:)), ...
                    'LineStyle', opt.undeflinestyle,...
                    'LineWidth', opt.undeflinewidth, ...
                    'Color', opt.undeflinecolor);
            end
            if (opt.plotdef),
                set(plot3(p0(1,1)+r(1,:)+s(1,:),...
                    p0(1,2)+r(2,:)+s(2,:),...
                    p0(1,3)+r(3,:)+s(3,:)), ...
                    'Color', opt.linecolor , 'LineWidth', opt.linewidth, ...
                    'LineStyle', opt.linestyle);
            end
        end
    else
        for kk=1:nBeams,
            p0 = model.nodes(model.beams.nodes(kk,:),:);
            s0 = Disp(model.beams.nodes(kk,:),:);
            if (opt.plotundef),
                set(plot3(p0(:,1), p0(:,2), p0(:,3)), ...
                    'LineStyle', opt.undeflinestyle, ...
                    'LineWidth', opt.undeflinewidth, ...
                    'Color', opt.undeflinecolor);
            end
            if (opt.plotdef),
                set(plot3(p0(:,1)+s0(:,1), ...
                    p0(:,2)+s0(:,2), ...
                    p0(:,3)+s0(:,3)), ...
                    'Color', opt.linecolor , 'LineWidth', opt.linewidth, ...
                    'LineStyle', opt.linestyle);
            end
        end
    end
end
%% node labels
if opt.nodeslabels
    if opt.plotundef
        for kk=1:size(model.nodes, 1)
            set(text(model.nodes(kk, 1), ...
                model.nodes(kk, 2), ...
                model.nodes(kk, 3), ...
                sprintf('%i', kk)),'FontSize', opt.FontSize);
        end
    end
    if opt.plotdef
        for kk=1:size(model.nodes, 1)
            set(text(model.nodes(kk, 1)+Disp(kk, 1), ...
                model.nodes(kk, 2)+Disp(kk, 2), ...
                model.nodes(kk, 3)+Disp(kk, 3), ...
                sprintf('%i', kk)),'FontSize', opt.FontSize, ...
                'EraseMode', 'xor');
        end
    end
end
%% plot node markers
if opt.plotundef
    set(plot3 (model.nodes(:,1), model.nodes(:,2), model.nodes(:,3)), ...
        'marker', opt.marker, ...
        'color', opt.undeflinecolor, ...
        'markersize', opt.markersize, ...
        'markerfacecolor', opt.undeflinecolor, ...
        'linestyle', 'none');
end
if opt.plotdef
    set(plot3 (model.nodes(:,1) + Disp(:,1), ...
        model.nodes(:,2) + Disp(:,2), ...
        model.nodes(:,3) + Disp(:,3)), ...
        'marker', opt.marker, ...
        'color', opt.linecolor, ...
        'markersize', opt.markersize, ...
        'markerfacecolor', opt.linecolor, ...
        'linestyle', 'none');
end
%% barslabels labels
if opt.barslabels
    for kk=1:size(model.beams.nodes, 1)
        r0 = mean([model.nodes(model.beams.nodes(kk,1),:);
            model.nodes(model.beams.nodes(kk,2),:)]);
        set(text(r0(1), r0(2), r0(3), ...
            sprintf('%i', kk)),'FontSize', opt.FontSize*0.7, ...
            'FontAngle', 'italic', 'EraseMode', 'xor');
    end
end

daspect([1, 1, 1]);
box on;
% camlight; lighting phong; material metal;
if (washold), hold on; end
% drawnow;

end
%% function [Nu, Nv] = shapefunc(xi, L0)
function [Nu, Nv, Nw] = shapefunc(xi, L0)

Nu = [(1-xi)/2; (1+xi)/2];
Nv = [(1-xi).^2.*(2+xi)/4;
    L0*(1-xi).^2.*(1+xi)/8;
    (1+xi).^2.*(2-xi)/4;
    -L0*(1+xi).^2.*(1-xi)/8];
Nw = [(1-xi).^2.*(2+xi)/4;
    -L0*(1-xi).^2.*(1+xi)/8;
    (1+xi).^2.*(2-xi)/4;
    L0*(1+xi).^2.*(1-xi)/8];
end
%% calcRotMatrix
function [T0, L0] = calcRotMatrix(r0)

deltar = r0(2,:)-r0(1,:);
L0 = norm(deltar);
e1 = deltar(:)/L0;
L12 = norm(e1(1:2));
if (L12==0),
    e2 = [0; 1; 0];
else
    e2 = [-e1(2); e1(1); 0]/L12;
end
e3 = cross(e1, e2);
T0 = [e1, e2, e3]';
end
