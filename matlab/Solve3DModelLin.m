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
    
    
    % asseble tria
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


