%% matprop = Find3DMatProp(model, {[a1, a2, a3], [dirs]})
%
function [Keps, Vol, Vol0, deps] = Find3DMatProp(model, varargin)
    
    %% preamble
    if nargin == 2,
        dirs = varargin{1};
        a1 = (model.nodes(dirs(1,2),:)-model.nodes(dirs(1,1), :))';
        a2 = (model.nodes(dirs(2,2),:)-model.nodes(dirs(2,1), :))';
        a3 = (model.nodes(dirs(3,2),:)-model.nodes(dirs(3,1), :))';
        Vol0= abs(det([a1, a2, a3]));            % unit cell volume
        [~, Eqs, C, RBeams, RFaces, Rtria] = ...
            UnitCellCheck3D(model, a1, a2, a3, ...
            struct('Echo', false));
    elseif nargin == 4,
        a1   = varargin{1};
        a2   = varargin{2};
        a3   = varargin{3};
        Vol0 = abs(det([a1, a2, a3]));            % unit cell volume
        [~, Eqs, C, RBeams, RFaces, Rtria] = ...
            UnitCellCheck3D(model, a1, a2, a3, ...
            struct('Echo', false));
    end
    %% find Eqs and reduce model to minimum form
    
    if ~isempty(RBeams),
        model.beams.nodes(RBeams,:) = [];       % remove duplicated beams
    end
    if ~isempty(RFaces),
        model.faces.nodes(RFaces) = [];         % remove duplicated faces
    end
    if ~isempty(Rtria),
        model.tria.nodes(Rtria,:) = [];     % remove duplicated faces
    end
    %% calc UC stiff Matrix
    [Kuc, Vol] = Solve3DModelLin(model);
    
    %% find stiffness matrix as a function of Deltaa
    [B0, B1] = makeB0andBa(Eqs, C, 3);
    nsp = null([B0'*Kuc*B0, B0'*Kuc*B1], 'r'); % this is more accurate than pinv
    D0 = nsp(1:end-9,end-8:end);
    
    dDeltaa = B0*D0+B1;
    KDeltaa = dDeltaa'*Kuc*dDeltaa;
    
    %% find stiffness matrix as a function of epsilon
    %       e11,    e22,    e33,    e23,  e31,  e12
    Beps = [bee(a1); bee(a2); bee(a3)];
    
    Keps = 1/Vol*Beps'*KDeltaa*Beps;
    deps = dDeltaa*Beps;
    
end
%% UnitCellModelCheck3D(model, a1, a2, a3, varargin)
function [minTopN, Eqs, C, RBeams, RFaces, RTria] = ...
        UnitCellCheck3D(model, a1, a2, a3, varargin)
    %% test
    Lstar = max([norm(a1), norm(a2), norm(a3)]);
    model.nodes = myfix(model.nodes/Lstar, 1e-8);
    a1 = myfix(a1/Lstar, 1e-8);
    a2 = myfix(a2/Lstar, 1e-8);
    a3 = myfix(a3/Lstar, 1e-8);
    %% check rank
    if rank([a1,a2,a3]) < 3,
        error('periodic directions are not independent');
    end
    %% preamble
    opt = [];
    N1 = 1;
    N2 = 1;
    N3 = 1;
    Echo = false;
    if (nargin >=5), opt = varargin{1}; end
    if isfield(opt, 'N1'), N1 = opt.N1; end
    if isfield(opt, 'N2'), N2 = opt.N2; end
    if isfield(opt, 'N3'), N1 = opt.N3; end
    if isfield(opt, 'Echo'), Echo = opt.Echo; end
    if (Echo), clc; end
    if ~isfield(model,'tria'), model.tria.nodes = []; end
    if ~isfield(model,'faces'), model.faces.nodes = []; end
    
    nodes = model.nodes; nNodes = size(nodes,1);
    
    if isfield(model, 'beams')
        beams = model.beams.nodes; 

    else
        beams = [];
    end
    
    nBeams = size(beams,1);
    faces = model.faces.nodes; nFaces = length(faces);
    tria = model.tria.nodes; nPlates = size(tria, 1);
    %% total topology for nodes
    
    [TotTopN, minTopN] = entityCheck(nodes, a1, a2, a3, N1, N2, N3);
    if (Echo),
        for kk=1:size(TotTopN,1),
            fprintf(...
                '%%  nodes(%02i) = nodes(%02i)%+i*a1%+i*a2%+i*a3;\n', ...
                TotTopN(kk,1), TotTopN(kk,2), TotTopN(kk,3), TotTopN(kk,4), TotTopN(kk,5));
        end
    end
    %% Build equilibrium equations
    
    if ~isempty(minTopN),
        nEqs = max(minTopN(:,2));
        Eqs = zeros(nEqs, nNodes);
        for kk = 1:nEqs,
            idxx = find(minTopN(:,2)==kk);
            Eqs(kk,kk) = 1;
            Eqs(kk,minTopN(idxx,1)) = ones(size(idxx));
        end
        
        % remove duplicate or incomplete equations
        for kk=1:nNodes,
            idxx = find(Eqs(:,kk));
            for jj=idxx',
                Eqs(idxx(1),:) = Eqs(idxx(1),:) | Eqs(jj,:);
            end
            Eqs(idxx(2:end),:)=[];
        end
        
        % add equations for the nodes that do not have correspondents on the unit
        % cell
        idxx = find(sum(Eqs,1)==0);
        for kk=idxx,
            Eqs(end+1,kk) = 1;  %#ok<*EMGRO>
        end
        
        % remove null equations
        normEqs = [];
        for kk=1:size(Eqs,1),
            normEqs(kk) = norm(Eqs(kk,:)); %#ok<AGROW>
        end
        Eqs(normEqs==0,:) = [];
    else
        nEqs = nNodes;
        Eqs = eye(nNodes);
    end
    %% conections matrix
    C = zeros(nEqs, 3);     %  number of independent periodic directions
    if ~isempty(minTopN),
        for kk = 1:nNodes,
            idxx = find(minTopN(:,1) == kk);
            if ~isempty(idxx),
                C(kk, :) = minTopN(idxx, 3:5);
            end
        end
    end
    %% look for repeated beams
    RBeams = [];
    cgBeams = zeros(nBeams,3);
    for kk=1:nBeams,
        cgBeams(kk,:) = mean(nodes(beams(kk,:),:));
    end
    [~, minTopB] = entityCheck(cgBeams, a1, a2, a3, N1, N2, N3);
    
    if ~isempty(minTopB),
        idxx = union(minTopB(:,2),[]);
        % corresponding Beams
        BBeams = cell(length(idxx),1);
        for kk=idxx',
            idrp = minTopB(:,2)==kk;
            BBeams{kk} = [kk, minTopB(idrp,1)'];
        end
        for kk=1:length(BBeams),
            dummy = BBeams{kk};
            RBeams = union(RBeams,  dummy(2:end));
        end
    else
        RBeams=[];
    end
    %% look for repeated faces
    RFaces= [];
    cgFaces = zeros(nFaces,3);
    for kk=1:nFaces,
        FaceNodes = faces{kk};
        cgFaces(kk,:) = mean(nodes(FaceNodes,:));
    end
    [~, minTopF] = entityCheck(cgFaces, a1, a2, a3, N1, N2, N3);
    
    if ~isempty(minTopF),
        idxx = union(minTopF(:,2),[]);
        % corresponding faces
        FFaces = cell(length(idxx),1);
        for kk=idxx',
            idrp = minTopF(:,2)==kk;
            FFaces{kk} = [kk, minTopF(idrp,1)'];
        end
        for kk=1:length(FFaces),
            dummy = FFaces{kk};
            RFaces = union(RFaces,  dummy(2:end));
        end
    else
        RFaces = [];
    end
    %% look for repeated tria
    cgplates = zeros(nPlates,3);
    for kk=1:nPlates,
        PlateNodes = tria(kk,:);
        cgplates(kk,:) = mean(nodes(PlateNodes,:));
    end
    [~, minTopF] = entityCheck(cgplates, a1, a2, a3, N1, N2, N3);
    
    if ~isempty(minTopF),
        idxx = union(minTopF(:,2),[]);
        % corresponding tria
        FPlates = cell(length(idxx),1);
        for kk=idxx',
            idrp = minTopF(:,2)==kk;
            FPlates{kk} = [kk, minTopF(idrp,1)'];
        end
        RTria= [];
        for kk=1:length(FPlates),
            dummy = FPlates{kk};
            RTria = union(RTria,  dummy(2:end));
        end
    else
        RTria = [];
    end
end
%% [TotTop, minTop] = entityCheck(entities, a1, a2, a3, N1, N2, N3)
function [TotTop, minTop] = entityCheck(entities, a1, a2, a3, N1, N2, N3)
    
    TotTop=[];
    nEntities = size(entities,1);
    
    %% find topology correspondaces
    for nn = nEntities:-1:1,
        count = 0;
        for ll = -N3:N3,
            for kk = -N2:N2,
                for hh = -N1:N1,
                    count = count +1;
                    dist = myfix(...
                        (entities(:, 1)-(entities(nn, 1)+hh*a1(1)+kk*a2(1)+ll*a3(1))).^2 + ...
                        (entities(:, 2)-(entities(nn, 2)+hh*a1(2)+kk*a2(2)+ll*a3(2))).^2 + ...
                        (entities(:, 3)-(entities(nn, 3)+hh*a1(3)+kk*a2(3)+ll*a3(3))).^2);
                    idx=find(dist == 0);
                    if (~isempty(idx)),
                        for ii=idx',
                            if (ii ~= nn),
                                TotTop(end+1,:) = [ii, nn, hh, kk, ll]; %#ok<AGROW>
                            end
                        end
                    end
                end
            end
        end
    end
    
    %% reduce topology correspondances to minimum
    TotTop = TotTop(end:-1:1,:);
    if ~isempty(TotTop),
        minTop = TotTop(1,:);
        for kk = 2:size(TotTop, 1),
            if ~any(TotTop(kk,2) == TotTop(1:kk,1));
                minTop = [minTop; TotTop(kk,:)]; %#ok<AGROW>
            end
        end
    else
        minTop = [];
    end
    
end
%% M = myfix(M, varargin)
function M = myfix(M, varargin)
    
    if (nargin <=1)
        M = fix(M/1e-8)*1e-8;
    else
        myeps = varargin{1};
        M = fix(M/myeps)*myeps;
    end
end
%% [B0, Ba] = makeB0andBa(Eqs, C, varargin)
function [B0, Ba] = makeB0andBa(Eqs, C, varargin)
    
    nDoF  = 6;   % number of DoFs per node
    naDoF = 6;   % number of DoFs per periodic direction
    if nargin>=3, naDoF = varargin{1}; end
    
    [nEqs, nNodes] = size(Eqs);
    
    A0 = zeros(nEqs*nDoF, nNodes*nDoF);
    for kk=1:nEqs,
        for jj=1:size(Eqs,2),
            A0(1+(kk-1)*nDoF:kk*nDoF, 1+(jj-1)*nDoF:jj*nDoF) = ...
                eye(nDoF)*Eqs(kk,jj);
        end
    end
    B0 = A0';
    
    nas = size(C,2);    % number of periodic directions, 3 for 3D
    Ba = zeros(nNodes*nDoF, nas*naDoF);
    for kk=1:nas,
        idxs = find(C(:, kk));
        for jj = idxs',
            Ba((jj-1)*nDoF+1:(jj-1)*nDoF+naDoF,...
                (kk-1)*naDoF+1:(kk-1)*naDoF+naDoF) = eye(naDoF)*C(jj,kk);
        end
    end
end
%% bee=bee(dr)
function bee=bee(dr)
    %       e11,    e22,    e33,    e23,        e31,        e12
    bee = [ dr(1),  0,      0,      0,          dr(3)/2,    dr(2)/2;
        0,      dr(2),  0,      dr(3)/2,    0,          dr(1)/2;
        0,      0,      dr(3),  dr(2)/2,    dr(1)/2,    0];
end
%% function [K, Vol] = Solve3DModelLin
function [K, Vol] = Solve3DModelLin(model)
    %% preamble and defaults
    %
    nDoF = 6;    
    nNodes = size(model.nodes, 1);    
    
    if isfield(model, 'beams')
        nBeams = size(model.beams.nodes, 1);
    else
        nBeams = 0;
    end
    %
    if isfield(model, 'tria')
        nTria = size(model.tria.nodes, 1);
    else
        nTria = 0;
    end    
    %% assemble stiffness matrix
    %    
    K = zeros(nDoF*nNodes);
    Vol = 0;
    
    E     = model.mat.E;
    nu    = model.mat.nu;
    G     = E/2/(1+nu);
    
    % Assemble beams
    if nBeams>0
        for kk=1:nBeams,
            [Kelkk, L] = makeBeamKMatrix(model.nodes(model.beams.nodes(kk,:), :), ...
                E, G, ...
                model.prop.beams.A, ...
                model.prop.beams.Ixx,...
                model.prop.beams.Iyy, ...
                model.prop.beams.Izz);
            
            indexg = [model.beams.nodes(kk,1)*nDoF-[5, 4, 3, 2, 1, 0], ...
                model.beams.nodes(kk,2)*nDoF-[5, 4, 3, 2, 1, 0]];
            K(indexg, indexg) = K(indexg, indexg) + Kelkk;
            Vol = Vol +  model.prop.beams.A*L;
        end
    end
    
    % asseble tria
    if nTria>0
        DMat = E/(1-nu^2)* ...
            [1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2]; % plane stress rigidity mtx
        thk = model.prop.tria.t;
        
        for kk=1:nTria
            idnodes = model.tria.nodes(kk,:);
            nodeskk = model.nodes(idnodes, :);
            
            [KShell, Areakk] = makeTriaShellKMatrix(nodeskk, DMat*thk);
            KPlate = makeTriaPlateKMatrix(nodeskk, DMat*thk^3/12);
            
            indexg = [model.tria.nodes(kk,1)*nDoF-[5, 4, 3, 2, 1, 0], ...
                model.tria.nodes(kk,2)*nDoF-[5, 4, 3, 2, 1, 0], ...
                model.tria.nodes(kk,3)*nDoF-[5, 4, 3, 2, 1, 0]];
            K(indexg, indexg) = K(indexg, indexg) + KShell + KPlate;
            
            Vol = Vol +  thk*Areakk;
        end
    end
    %% return
end
%% functions for beams

function varargout = makeBeamKMatrix(p0, E, G, A, Ixx, Iyy, Izz, varargin)
    [T0, L] = calcRotMatrix(p0);                % calc. element cosines
    T = blkdiag(T0, T0, T0, T0);                % assemble DoF rotation matrix
    
    Kell = makeBeamKlMatrix(E,G,L,A,Ixx,Iyy,Izz);
    
    if nargin == 7
        varargout(1) = {T'*Kell*T};
        varargout(2) = {L};
    else
        Disp = varargin{1};
        dispg = [Disp(1,:), Disp(2,:)]';    % nodal DoFs in the global ref sys
        displ = T*dispg;                    % nodal DoFs in the local ref sys
        Fil = Kell*displ;                       % forces in the local ref sys
        Uax = 1/2*Fil([3,9])'*displ([3,9]);     % axial strain energy
        Ubd = 1/2*Fil([1,2,4,5,7,8,10,11])'*...
            displ([1,2,4,5,7,8,10,11]);         % bending strain energy
        Ut = 1/2*Fil([6, 12])'*displ([6,12]);   % torsion strain energy
        varargout = {[Uax, Ubd, Ut]};
    end
end

function Kell = makeBeamKlMatrix(E,G,L,A,Ixx,Iyy,Izz)
    % G = E/2/(1+nu);
    EIxx = E*Ixx; EIyy = E*Iyy; GIzz = G*Izz; EA = E*A;
    Kell = reshape([EIyy.*1.0./L.^3.*1.2e1,0.0,0.0,0.0,...
        EIyy.*1.0./L.^2.*6.0,0.0,EIyy.*1.0./L.^3.*-1.2e1,0.0,...
        0.0,0.0,EIyy.*1.0./L.^2.*6.0,0.0,0.0,EIxx.*1.0./L.^3.*1.2e1,...
        0.0,EIxx.*1.0./L.^2.*-6.0,0.0,0.0,0.0,EIxx.*1.0./L.^3.*-1.2e1,...
        0.0,EIxx.*1.0./L.^2.*-6.0,0.0,0.0,0.0,0.0,EA./L,0.0,0.0,0.0,...
        0.0,0.0,-EA./L,0.0,0.0,0.0,0.0,EIxx.*1.0./L.^2.*-6.0,0.0,...
        (EIxx.*4.0)./L,0.0,0.0,0.0,EIxx.*1.0./L.^2.*6.0,0.0,...
        (EIxx.*2.0)./L,0.0,0.0,EIyy.*1.0./L.^2.*6.0,0.0,0.0,0.0,...
        (EIyy.*4.0)./L,0.0,EIyy.*1.0./L.^2.*-6.0,0.0,0.0,0.0,...
        (EIyy.*2.0)./L,0.0,0.0,0.0,0.0,0.0,0.0,GIzz./L,0.0,0.0,0.0,...
        0.0,0.0,-GIzz./L,EIyy.*1.0./L.^3.*-1.2e1,0.0,0.0,0.0,...
        EIyy.*1.0./L.^2.*-6.0,0.0,EIyy.*1.0./L.^3.*1.2e1,0.0,0.0,...
        0.0,EIyy.*1.0./L.^2.*-6.0,0.0,0.0,EIxx.*1.0./L.^3.*-1.2e1,...
        0.0,EIxx.*1.0./L.^2.*6.0,0.0,0.0,0.0,EIxx.*1.0./L.^3.*1.2e1,...
        0.0,EIxx.*1.0./L.^2.*6.0,0.0,0.0,0.0,0.0,-EA./L,0.0,0.0,0.0,0.0,...
        0.0,EA./L,0.0,0.0,0.0,0.0,EIxx.*1.0./L.^2.*-6.0,0.0,...
        (EIxx.*2.0)./L,0.0,0.0,0.0,EIxx.*1.0./L.^2.*6.0,0.0,...
        (EIxx.*4.0)./L,0.0,0.0,EIyy.*1.0./L.^2.*6.0,0.0,0.0,0.0,...
        (EIyy.*2.0)./L,0.0,EIyy.*1.0./L.^2.*-6.0,0.0,0.0,0.0,...
        (EIyy.*4.0)./L,0.0,0.0,0.0,0.0,0.0,0.0,-GIzz./L,0.0,0.0,...
        0.0,0.0,0.0,GIzz./L],[12,12]);
end

function [T0, L0] = calcRotMatrix(r0)
    deltar = r0(2,:)-r0(1,:);
    L0 = norm(deltar);
    e3 = deltar/L0;
    if abs(e3(3))==1,
        e2=[1,0,0];
    else
        e2 = cross(e3, [0,0,1]);
        e2 = e2/norm(e2);
    end
    e1 = cross(e2,e3);
    T0 = [e1; e2; e3];
end
%% functions for tria
function [KShell, Area] = makeTriaShellKMatrix(p0g, Kmat, varargin)
    
    T0 = calcPlateRotMatrix(p0g);
    T = blkdiag(T0, T0, T0, T0, T0, T0);
    p0 = p0g*T0';
    [B0, Area] = makeShellBMatrix(p0);
    
    idxl = [1, 2, 7, 8, 13, 14];
    Kell = Area*B0'*Kmat*B0;
    
    Kel0 = zeros(18);
    Kel0(idxl, idxl) = Kell;
    Kel0 = T'*Kel0*T;
    KShell = Kel0;    
end
%
function [B, A] = makeShellBMatrix(p0)
    indexu = (1:2:6);	indexv = (2:2:6);
    [A, Nx, Ny] = calcNvectors(p0);
    Nux = zeros(6,1);   Nuy = zeros(6,1);
    Nvx = zeros(6,1);   Nvy = zeros(6,1);
    Nux(indexu) = Nx;	Nuy(indexu) = Ny;
    Nvx(indexv) = Nx;	Nvy(indexv) = Ny;
    B = [Nux'; Nvy'; Nuy'+Nvx'];
end
function [A, Nx, Ny] = calcNvectors(p0)
    x0 = p0(:,1)';	y0 = p0(:,2)';
    A = abs(det([x0(2:3)-x0(1);y0(2:3)-y0(1)]))/2;
    Nx = [(y0(2) - y0(3)) / (x0(1) * y0(2) - x0(1) * y0(3) + x0(2) * y0(3) - x0(2) * y0(1) + x0(3) * y0(1) - x0(3) * y0(2)) (y0(3) - y0(1)) / (x0(1) * y0(2) - x0(1) * y0(3) + x0(2) * y0(3) - x0(2) * y0(1) + x0(3) * y0(1) - x0(3) * y0(2)) (y0(1) - y0(2)) / (x0(1) * y0(2) - x0(1) * y0(3) + x0(2) * y0(3) - x0(2) * y0(1) + x0(3) * y0(1) - x0(3) * y0(2))]';
    Ny = [(-x0(2) + x0(3)) / (x0(1) * y0(2) - x0(1) * y0(3) + x0(2) * y0(3) - x0(2) * y0(1) + x0(3) * y0(1) - x0(3) * y0(2)) (x0(1) - x0(3)) / (x0(1) * y0(2) - x0(1) * y0(3) + x0(2) * y0(3) - x0(2) * y0(1) + x0(3) * y0(1) - x0(3) * y0(2)) (-x0(1) + x0(2)) / (x0(1) * y0(2) - x0(1) * y0(3) + x0(2) * y0(3) - x0(2) * y0(1) + x0(3) * y0(1) - x0(3) * y0(2))]';
end
function KPlate = makeTriaPlateKMatrix(p0g, Mat, varargin)
    T0 = calcPlateRotMatrix(p0g);
    T = blkdiag(T0, T0, T0, T0, T0, T0);
    p0 = p0g*T0';
    x1 = p0(1,1); x2 = p0(2,1); x3 = p0(3,1);
    y1 = p0(1,2); y2 = p0(2,2); y3 = p0(3,2);
    
    % A = (-x2*y1 + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3)/2;
    % use the following that come from the simplification
    % is not the area of the tringle
    A = (-x2*y1 + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3);
    J2 = [(y2 - y3)^2, (y1 - y3)^2, (y1 - y2)^2, 2*(y1 - y3)*(y3 - y2),...
        -2*(y1 - y2)*(y1 - y3), 2*(y1 - y2)*(y2 - y3);
        (x2 - x3)^2, (x1 - x3)^2, (x1 - x2)^2, 2*(x1 - x3)*(x3 - x2),...
        -2*(x1 - x2)*(x1 - x3), 2*(x1 - x2)*(x2 - x3);
        (-(x2 - x3))*(y2 - y3), (-(x1 - x3))*(y1 - y3), ...
        (-(x1 - x2))*(y1 - y2),...
        (-x3)*(y1 + y2 - 2*y3) + x2*(y1 - y3) + x1*(y2 - y3), ...
        x3*(y2 - y1) + x1*(2*y1 - y2 - y3) + x2*(y3 - y1),...
        x3*(y1 - y2) - x2*(y1 - 2*y2 + y3) + x1*(y3 - y2)]/A^2;
    J2(3,:) = J2(3,:)*2;    % last row of the Jacobian squared is D[w,{x,x}]
    % but we need 2D[w,{x,x}]
    
    B12 = J2*calcPlateHessian([0.5, 0.5, 0], p0);
    B23 = J2*calcPlateHessian([0, 0.5, 0.5], p0);
    B31 = J2*calcPlateHessian([0.5, 0, 0.5], p0);
    
    idxl = [3, 4, 5, 9, 10, 11, 15, 16, 17];	%indeces of local DoFs
    
    Kell = (B12'*Mat*B12 + ...
        B23'*Mat*B23 + ...
        B31'*Mat*B31)*A/2/3;
    % A is twice the area
    % /3 is the integation formula
    
    Kel0 = zeros(18);
    Kel0(idxl,idxl) = Kell;
    KPlate = T'*Kel0*T;
end
function H = calcPlateHessian(phi, p0)
    % this implents th BCIZ plate element
    
    phi1 = phi(1); phi2 = phi(2); phi3 = phi(3);
    x1 = p0(1,1); x2 = p0(2,1); x3 = p0(3,1);
    y1 = p0(1,2); y2 = p0(2,2); y3 = p0(3,2);
    
    H = [6*(phi1+phi2+phi3), 2*(phi2*(y2-y1)+phi3*(y3-y1)), ...
        2*(phi2*(x1-x2)+phi3*(x1-x3)), 0, 0, 0, 0, 0, 0;
        0, 0, 0, 6*(phi1+phi2+phi3), 2*(phi1*(y1-y2)+phi3*(y3-y2)), ...
        2*phi1*(x2-x1)+2*phi3*(x2-x3), 0, 0, 0;
        0, 0, 0, 0, 0, 0, 6*(phi1+phi2+phi3), 2*(phi1*(y1-y3)+phi2*(y2-y3)), ...
        2*phi1*(x3-x1)+2*phi2*(x3-x2);
        2*(3*phi1+phi3), 1/2*(phi3*(-2*y1+y2+y3)-4*phi1*(y1-y2)), ...
        2*phi1*(x1-x2)+1/2*phi3*(2*x1-x2-x3), 2*(3*phi2+phi3), ...
        1/2*(4*phi2*(y1-y2)+phi3*(y1-2*y2+y3)), ...
        -2*phi2*(x1-x2)-1/2*phi3*(x1-2*x2+x3), 2*phi3, ...
        1/2*phi3*(y1+y2-2*y3), -(1/2)*phi3*(x1+x2-2*x3);
        2*phi1, 1/2*phi1*(-2*y1+y2+y3), 1/2*phi1*(2*x1-x2-x3), ...
        2*(phi1+3*phi2), 1/2*(phi1*(y1-2*y2+y3)+4*phi2*(y3-y2)), ...
        2*phi2*(x2-x3)-1/2*phi1*(x1-2*x2+x3), 2*(phi1+3*phi3), ...
        1/2*(phi1*(y1+y2-2*y3)+4*phi3*(y2-y3)), ...
        2*phi3*(x3-x2)-1/2*phi1*(x1+x2-2*x3);
        2*(3*phi1+phi2), 1/2*(phi2*(-2*y1+y2+y3)-4*phi1*(y1-y3)), ...
        2*phi1*(x1-x3)+1/2*phi2*(2*x1-x2-x3), 2*phi2, 1/2*phi2*(y1-2*y2+y3), ...
        -(1/2)*phi2*(x1-2*x2+x3), 2*(phi2+3*phi3), ...
        1/2*(phi2*(y1+y2-2*y3)+4*phi3*(y1-y3)), ...
        2*phi3*(x3-x1)-1/2*phi2*(x1+x2-2*x3)];
end
function [hr, hs] = calchvectors(r,s)
    % called by makeQuadShellKMatrix
    hr = [s-1, 1-s, s, -s]';
    hs = [r-1, -r, r, 1-r]';
end
function T0 = calcPlateRotMatrix(p0)
    e1 = p0(2,:)-p0(1,:); e1 = e1/norm(e1);
    e3 = cross(e1, p0(3,:)-p0(1,:)); e3 = e3/norm(e3);
    e2 = cross(e3,e1);
    T0 = [e1', e2', e3']';
end


