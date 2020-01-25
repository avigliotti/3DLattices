(* ::Package:: *)

(* ::Section:: *)
(*Calc3DMatProp*)


BeginPackage["Calc3DMatProp`"]


(* ::Subsection:: *)
(*Misc Funcs*)


BlockDiagonalMatrix[b:{__?MatrixQ}] := Block[{r, c, n = Length[b], i, j}, 
    {r, c} = Transpose[Dimensions /@ b]; ArrayFlatten[
      Table[If[i == j, b[[i]], ConstantArray[0, {r[[i]], c[[j]]}]], {i, n}, 
       {j, n}]]]; 


MakeB0andBe[Eqs_, nodes_] := Block[{jj, kk, idxi, A0, nEqs, 
     nDoF, nNodes, Be, bee, idxt, nOnes, idxs}, 
    nDoF = 6;
bee[dr_] := {{dr[[1]], 0, 0, 0, dr[[3]]/2, 
        dr[[2]]/2}, {0, dr[[2]], 0, dr[[3]]/2, 0, dr[[1]]/2}, 
       {0, 0, dr[[3]], dr[[2]]/2, dr[[1]]/2, 0}}; 
nNodes = Dimensions[nodes][[1]]; 
     nEqs = Dimensions[Eqs][[1]]; 
     A0 = ConstantArray[0, {nDoF*nEqs, nDoF*nNodes}]; 
     idxi = Range[nDoF] - 1;
For[kk = 1, kk <= nEqs, kk++,
For[jj = 1, jj <= nNodes, jj++, 
       A0[[kk*nDoF - idxi,jj*nDoF - idxi]] = 
        IdentityMatrix[nDoF]*Eqs[[kk,jj]]]]; 
     Be = ConstantArray[0, {nDoF*nNodes, nDoF}]; 
     idxt = {5, 4, 3}; For[kk = 1, kk <= nEqs, kk++, 
      idxs = Position[Eqs[[kk]], 1]; nOnes = Dimensions[idxs][[
         1]]; If[nOnes > 1, For[jj = 2, jj <= nOnes, jj++, 
         Be[[nDoF*idxs[[jj,1]] - idxt]] = 
          bee[nodes[[idxs[[jj,1]]]] - nodes[[idxs[[1,1]]]]]]]]; 
     Return[{Transpose[A0], Be}]]; 


PlateRotationMatrix[r0_] := Block[{e1, e2, e3}, 
    e1 = Simplify[(r0[[2]] - r0[[1]])/Norm[r0[[2]] - r0[[1]]]]; 
     e3 = Cross[e1, r0[[3]] - r0[[1]]]; e3 = e3/Norm[e3]; e2 = Cross[e3, e1]; 
     Return[{e1, e2, e3}]]; 


BeamRotationMatrix[deltar_] := Block[{L0, e1, e2, e3, L12, T0}, 
    L0 = Norm[deltar]; e1 = Simplify[deltar/L0, L0 > 0]; 
     If[Simplify[e1[[3]] == 1], e2 = {1, 0, 0}, e2 = Cross[e1, {0, 0, 1}]; 
       e2 = Simplify[e2/Norm[e2], L0 > 0]]; e3 = Simplify[Cross[e1, e2], L0 > 0]; 
     T0 = {e1, e2, e3}; Return[{T0, L0}]]; 


CalcMohrCircle[\[Sigma]_] := Block[{C, R}, C = (\[Sigma][[1]] + \[Sigma][[2]])/2; 
    R = Sqrt[((\[Sigma][[1]] - \[Sigma][[2]])/2)^2 + \[Sigma][[3]]]; Return[Simplify[{C - R, C + R}]]];


CalcVMStress[\[Sigma]_] := Block[{\[Sigma]1, \[Sigma]2}, {\[Sigma]1, \[Sigma]2} = CalcMohrCircle[\[Sigma]]; 
     Return[Simplify[\[Sigma]1^2 + \[Sigma]2^2 - \[Sigma]1 \[Sigma]2]]]; 


Calc\[Tau]maxStress[\[Sigma]_] := Block[{\[Sigma]1, \[Sigma]2}, {\[Sigma]1, \[Sigma]2} = CalcMohrCircle[\[Sigma]]; 
     Return[Simplify[Abs[(\[Sigma]2 - \[Sigma]1)/2]]]]; 


AddNode[node_, nodes_]:=Block[{}, Return[Join[nodes,node]]];


(* ::Section:: *)
(*Finite Element Funcs*)


CalcModelVolume[nodes_, beams_, shells_, beamprop_, shellprop_] := 
   Block[{nbeams, nshells, kk, Vol, A, Ix, Iy, Iz, Em, \[Nu], t}, 
    {A, Ix, Iy, Iz} = beamprop; {Em, \[Nu], t} = shellprop; 
     nbeams = Dimensions[beams][[1]]; nshells = Dimensions[shells][[1]]; Vol = 0; 
     For[kk = 1, kk <= nbeams, kk++, 
      Vol += A*Simplify[Norm[nodes[[beams[[kk,2]]]] - nodes[[beams[[kk,1]]]]]]; ]; 
     For[kk = 1, kk <= nshells, kk++, 
      Vol += (1/2)*t*Simplify[Norm[Cross[nodes[[shells[[kk,2]]]] - 
            nodes[[shells[[kk,1]]]], nodes[[shells[[kk,3]]]] - 
            nodes[[shells[[kk,1]]]]]]]];
Return[Simplify[Vol]]]; 


Make3DStiffMatrix[nodes_, beams_, shells_, beamat_, beamprop_, shellprop_, flags_] := 
   Block[{idxg, Kfrm, Kel, nbeams, nnodes, kk, ii, jj, nshells, beamon, shellon, 
     plateon}, {beamon, shellon, plateon} = flags; nbeams = Dimensions[beams][[1]]; 
     nnodes = Dimensions[nodes][[1]]; nshells = Dimensions[shells][[1]]; 
     Kfrm = ConstantArray[0, {6*nnodes, 6*nnodes}]; 
     If[beamon == 1, For[kk = 1, kk <= nbeams, kk++, 
       idxg = Join[6*beams[[kk]][[1]] - Range[5, 0, -1], 6*beams[[kk]][[2]] - 
           Range[5, 0, -1]]; Kel = Beam3DStiffness[nodes[[beams[[kk]][[2]]]] - 
           nodes[[beams[[kk]][[1]]]], beamat, beamprop];
For[ii = 1, ii <= 12, ii++, 
         For[jj = 1, jj <= 12, jj++, Kfrm[[idxg[[ii]]]][[idxg[[jj]]]] += 
           Kel[[ii]][[jj]]]]]];
If[shellon == 1, For[kk = 1, kk <= nshells, kk++, 
       idxg = Join[6*shells[[kk]][[1]] - Range[5, 0, -1], 6*shells[[kk]][[2]] - 
           Range[5, 0, -1], 6*shells[[kk]][[3]] - Range[5, 0, -1]];
Kel = Shell3DStiffness[{nodes[[shells[[kk]][[1]]]], 
           nodes[[shells[[kk]][[2]]]], nodes[[shells[[kk]][[3]]]]}, shellprop];
For[ii = 1, ii <= 18, ii++, For[jj = 1, jj <= 18, jj++, 
          Kfrm[[idxg[[ii]]]][[idxg[[jj]]]] += Kel[[ii]][[jj]]]]]];
If[plateon == 1, For[kk = 1, kk <= nshells, kk++, 
       idxg = Join[6*shells[[kk]][[1]] - Range[5, 0, -1], 6*shells[[kk]][[2]] - 
           Range[5, 0, -1], 6*shells[[kk]][[3]] - Range[5, 0, -1]];
Kel = plateon*Plate3DStiffness[{nodes[[shells[[kk]][[1]]]], 
            nodes[[shells[[kk]][[2]]]], nodes[[shells[[kk]][[3]]]]}, shellprop];
For[ii = 1, ii <= 18, ii++, For[jj = 1, jj <= 18, jj++, 
          Kfrm[[idxg[[ii]]]][[idxg[[jj]]]] += Kel[[ii]][[jj]]]]]]; 
     Return[Simplify[Kfrm]]]; 


(* ::Subsection:: *)
(*Beams*)


Beam3DStiffness[r12_, mprop_, fprop_] := Block[{Em, Gm, A, Ix, Iy, Iz, Kww, Kvv, Kuu, 
     L0, Kel, R0, R}, {R0, L0} = BeamRotationMatrix[r12]; 
     R = BlockDiagonalMatrix[{R0, R0, R0, R0}]; {Em, Gm} = mprop; 
     {A, Ix, Iy, Iz} = fprop; Kuu = {{1, -1}, {-1, 1}}; 
     Kww = {{12/L0^3, -6/L0^2, -12/L0^3, -6/L0^2}, {-6/L0^2, 4/L0, 6/L0^2, 2/L0}, 
       {-12/L0^3, 6/L0^2, 12/L0^3, 6/L0^2}, {-6/L0^2, 2/L0, 6/L0^2, 4/L0}}; 
     Kvv = {{12/L0^3, 6/L0^2, -12/L0^3, 6/L0^2}, {6/L0^2, 4/L0, -6/L0^2, 2/L0}, 
       {-12/L0^3, -6/L0^2, 12/L0^3, -6/L0^2}, {6/L0^2, 2/L0, -6/L0^2, 4/L0}}; 
     Kel = ConstantArray[0, {12, 12}]; Kel[[{1, 7},{1, 7}]] = ((Em*A)/L0)*Kuu; 
     Kel[[{4, 10},{4, 10}]] = ((Gm*Ix)/L0)*Kuu; Kel[[{2, 6, 8, 12},{2, 6, 8, 12}]] = 
      Em*Iz*Kvv; Kel[[{3, 5, 9, 11},{3, 5, 9, 11}]] = Em*Iy*Kww; 
     Return[Simplify[Transpose[R] . Kel . R, L0 > 0]]]; 


Frame3DStiffMatrix[nodes_, elem_, mprop_, beamprop_, t_] := 
   Block[{idxg, Kfrm, Kel, nelem, nnodes, kk, ii, jj}, nelem = Dimensions[elem][[1]]; 
     nnodes = Dimensions[nodes][[1]]; Kfrm = ConstantArray[0, {6*nnodes, 6*nnodes}]; 
     For[kk = 1, kk <= nelem, kk++, idxg = Join[6*elem[[kk]][[1]] - Range[5, 0, -1], 
         6*elem[[kk]][[2]] - Range[5, 0, -1]]; 
       Kel = Beam3DStiffness[nodes[[elem[[kk]][[2]]]] - nodes[[elem[[kk]][[1]]]], 
         mprop, beamprop]; For[ii = 1, ii <= 12, ii++, For[jj = 1, jj <= 12, jj++, 
         Kfrm[[idxg[[ii]]]][[idxg[[jj]]]] += Kel[[ii]][[jj]]]]]; 
     Return[Simplify[Kfrm]]]; 


Beam3DLocForces[r12_, dispg_, mprop_, fprop_] := 
   Block[{Em, Gm, A, Ix, Iy, Iz, Kww, Kvv, Kuu, L0, Kel, R0, R}, 
    {R0, L0} = BeamRotationMatrix[r12]; R = BlockDiagonalMatrix[{R0, R0, R0, R0}]; 
     {Em, Gm} = mprop; {A, Ix, Iy, Iz} = fprop; Kuu = {{1, -1}, {-1, 1}}; 
     Kww = {{12/L0^3, -6/L0^2, -12/L0^3, -6/L0^2}, {-6/L0^2, 4/L0, 6/L0^2, 2/L0}, 
       {-12/L0^3, 6/L0^2, 12/L0^3, 6/L0^2}, {-6/L0^2, 2/L0, 6/L0^2, 4/L0}}; 
     Kvv = {{12/L0^3, 6/L0^2, -12/L0^3, 6/L0^2}, {6/L0^2, 4/L0, -6/L0^2, 2/L0}, 
       {-12/L0^3, -6/L0^2, 12/L0^3, -6/L0^2}, {6/L0^2, 2/L0, -6/L0^2, 4/L0}}; 
     Kel = ConstantArray[0, {12, 12}]; Kel[[{1, 7},{1, 7}]] = ((Em*A)/L0)*Kuu; 
     Kel[[{4, 10},{4, 10}]] = ((Gm*Ix)/L0)*Kuu; Kel[[{2, 6, 8, 12},{2, 6, 8, 12}]] = 
      Em*Iz*Kvv; Kel[[{3, 5, 9, 11},{3, 5, 9, 11}]] = Em*Iy*Kww; 
     Return[Kel . R . dispg]]; 


GetBeamStress[nbeam_, disp_, nodes_, beams_, mprop_, fprop_] := 
   Block[{r12, idxg, elforces, T0, L0, A, Ix, Iy, Iz}, {A, Ix, Iy, Iz} = fprop; 
     r12 = nodes[[beams[[nbeam,2]]]] - nodes[[beams[[nbeam,1]]]]; 
     idxg = Join[6*beams[[nbeam]][[1]] - Range[5, 0, -1], 6*beams[[nbeam]][[2]] - Range[5, 0, -1]]; 
     elforces = Beam3DLocForces[r12, disp[[idxg]], mprop, fprop]; 
     Return[Simplify[{elforces[[7]]/A, elforces[[11]]/Iy, elforces[[12]]/Iz, 
        elforces[[10]]/Ix}]]; ]; 


GetBeamBuckl[nbeam_, disp_, nodes_, beams_, mprop_, fprop_] := 
   Block[{r12, idxg, elforces, T0, A, Ix, Iy, Iz, L0, \[Sigma]crz, \[Sigma]cry, Em,Gm}, {Em, Gm} = mprop;
    {A, Ix, Iy, Iz} = fprop; r12 = nodes[[beams[[nbeam,2]]]] - nodes[[beams[[nbeam,1]]]]; 
     L0=Simplify[Norm[r12]];
idxg = Join[6*beams[[nbeam]][[1]] - Range[5, 0, -1], 
       6*beams[[nbeam]][[2]] - Range[5, 0, -1]]; elforces = Beam3DLocForces[r12, disp[[idxg]], 
       mprop, fprop]; \[Sigma]cry = (Pi*Iy)/(A*L0^2); \[Sigma]crz = (Pi*Iz)/(A*L0^2); 
     Return[Simplify[{elforces[[7]]/A/\[Sigma]cry,elforces[[7]]/A/\[Sigma]crz}]]; ]; 


GetBeamForces[nbeam_, disp_, nodes_, beams_, mprop_, fprop_] := 
   Block[{r12, idxg, elforces, T0, L0, A, Ix, Iy, Iz}, {A, Ix, Iy, Iz} = fprop; 
     r12 = nodes[[beams[[nbeam,2]]]] - nodes[[beams[[nbeam,1]]]]; 
     idxg = Join[6*beams[[nbeam]][[1]] - Range[5, 0, -1], 6*beams[[nbeam]][[2]] - Range[5, 0, -1]]; 
     elforces = Beam3DLocForces[r12, disp[[idxg]], mprop, fprop]; 
     Return[Simplify[elforces[[Range[7,12,1]]]]]; ]; 


(* ::Subsection:: *)
(*Plates*)


Plate3DStiffness[r0g_, prop_] := Block[{R, R0, r0l, Kll, idxx, K}, 
    R0 = PlateRotationMatrix[r0g]; R = BlockDiagonalMatrix[
       {R0, R0, R0, R0, R0, R0}]; r0l = Simplify[r0g . Transpose[R0]]; 
     K = ConstantArray[0, {18, 18}]; idxx = {3, 4, 5, 9, 10, 11, 15, 16, 17}; 
     K[[idxx,idxx]] = MakePlateKll[r0l, prop]; 
     Return[Simplify[Transpose[R] . K . R]]; ]; 


MakePlateKll[r0l_, prop_] := Block[{x1, y1, z1, x2, y2, z2, x3, y3, z3, x, y, t, 
     Em, \[Nu], \[CapitalPhi]0, \[CapitalPhi], B, \[Kappa], A, J2, H\[CapitalPhi], Fk, Kmat}, 
    {Em, \[Nu], t} = prop; {{x1, y1, z1}, {x2, y2, z2}, {x3, y3, z3}} = r0l; 
     A = (1/2)*Norm[Cross[r0l[[2]] - r0l[[1]], r0l[[3]] - r0l[[1]]]]; 
     \[Kappa] = Simplify[((t^3*Em)/(1 - \[Nu]^2)/12)*{{1, \[Nu], 0}, {\[Nu], 1, 0}, 
         {0, 0, (1 - \[Nu])/2}}]; 
     J2 = {{(y2 - y3)^2, (y1 - y3)^2, (y1 - y2)^2, 2*(y1 - y3)*(-y2 + y3), 
         -2*(y1 - y2)*(y1 - y3), 2*(y1 - y2)*(y2 - y3)}, {(x2 - x3)^2, (x1 - x3)^2, 
         (x1 - x2)^2, 2*(x1 - x3)*(-x2 + x3), -2*(x1 - x2)*(x1 - x3), 
         2*(x1 - x2)*(x2 - x3)}, {-((x2 - x3)*(y2 - y3)), -((x1 - x3)*(y1 - y3)), 
         -((x1 - x2)*(y1 - y2)), -(x3*(y1 + y2 - 2*y3)) + x2*(y1 - y3) + 
          x1*(y2 - y3), x3*(-y1 + y2) + x1*(2*y1 - y2 - y3) + x2*(-y1 + y3), 
         x3*(y1 - y2) - x2*(y1 - 2*y2 + y3) + x1*(-y2 + y3)}}*(1/A^2); 
     J2[[3]] = 2*J2[[3]]; H\[CapitalPhi][\[Phi]1_, \[Phi]2_, \[Phi]3_] := 
      {{6*(\[Phi]1 + \[Phi]2 + \[Phi]3), 2*(y2*\[Phi]2 + y3*\[Phi]3 - y1*(\[Phi]2 + \[Phi]3)), 
        2*(-(x2*\[Phi]2) - x3*\[Phi]3 + x1*(\[Phi]2 + \[Phi]3)), 0, 0, 0, 0, 0, 0}, 
       {0, 0, 0, 6*(\[Phi]1 + \[Phi]2 + \[Phi]3), 2*(y1*\[Phi]1 + y3*\[Phi]3 - y2*(\[Phi]1 + \[Phi]3)), 
        2*(-(x1*\[Phi]1) - x3*\[Phi]3 + x2*(\[Phi]1 + \[Phi]3)), 0, 0, 0}, 
       {0, 0, 0, 0, 0, 0, 6*(\[Phi]1 + \[Phi]2 + \[Phi]3), 
        2*(y1*\[Phi]1 + y2*\[Phi]2 - y3*(\[Phi]1 + \[Phi]2)), 
        2*(-(x1*\[Phi]1) - x2*\[Phi]2 + x3*(\[Phi]1 + \[Phi]2))}, {2*(3*\[Phi]1 + \[Phi]3), 
        (y3*\[Phi]3 - 2*y1*(2*\[Phi]1 + \[Phi]3) + y2*(4*\[Phi]1 + \[Phi]3))/2, 
        x1*(2*\[Phi]1 + \[Phi]3) + (-(x3*\[Phi]3) - x2*(4*\[Phi]1 + \[Phi]3))/2, 2*(3*\[Phi]2 + \[Phi]3), 
        (y3*\[Phi]3 - 2*y2*(2*\[Phi]2 + \[Phi]3) + y1*(4*\[Phi]2 + \[Phi]3))/2, 
        -(x3*\[Phi]3)/2 + x2*(2*\[Phi]2 + \[Phi]3) - (x1*(4*\[Phi]2 + \[Phi]3))/2, 2*\[Phi]3, 
        ((y1 + y2 - 2*y3)*\[Phi]3)/2, -((x1 + x2 - 2*x3)*\[Phi]3)/2}, 
       {2*\[Phi]1, ((-2*y1 + y2 + y3)*\[Phi]1)/2, ((2*x1 - x2 - x3)*\[Phi]1)/2, 
        2*(\[Phi]1 + 3*\[Phi]2), (y1*\[Phi]1 - 2*y2*(\[Phi]1 + 2*\[Phi]2) + y3*(\[Phi]1 + 4*\[Phi]2))/2, 
        -(x1*\[Phi]1)/2 + x2*(\[Phi]1 + 2*\[Phi]2) - (x3*(\[Phi]1 + 4*\[Phi]2))/2, 
        2*(\[Phi]1 + 3*\[Phi]3), (y1*\[Phi]1 - 2*y3*(\[Phi]1 + 2*\[Phi]3) + y2*(\[Phi]1 + 4*\[Phi]3))/2, 
        -(x1*\[Phi]1)/2 + x3*(\[Phi]1 + 2*\[Phi]3) - (x2*(\[Phi]1 + 4*\[Phi]3))/2}, 
       {2*(3*\[Phi]1 + \[Phi]2), (y2*\[Phi]2 - 2*y1*(2*\[Phi]1 + \[Phi]2) + y3*(4*\[Phi]1 + \[Phi]2))/2, 
        x1*(2*\[Phi]1 + \[Phi]2) + (-(x2*\[Phi]2) - x3*(4*\[Phi]1 + \[Phi]2))/2, 2*\[Phi]2, 
        ((y1 - 2*y2 + y3)*\[Phi]2)/2, -((x1 - 2*x2 + x3)*\[Phi]2)/2, 2*(\[Phi]2 + 3*\[Phi]3), 
        (y2*\[Phi]2 - 2*y3*(\[Phi]2 + 2*\[Phi]3) + y1*(\[Phi]2 + 4*\[Phi]3))/2, 
        -(x2*\[Phi]2)/2 + x3*(\[Phi]2 + 2*\[Phi]3) - (x1*(\[Phi]2 + 4*\[Phi]3))/2}}; 
     B[\[Phi]1_, \[Phi]2_, \[Phi]3_] := J2 . H\[CapitalPhi][\[Phi]1, \[Phi]2, \[Phi]3]; 
     Fk[\[Phi]1_, \[Phi]2_, \[Phi]3_] := Transpose[B[\[Phi]1, \[Phi]2, \[Phi]3]] . \[Kappa] . 
       B[\[Phi]1, \[Phi]2, \[Phi]3]; Kmat = Simplify[(A/3)*(Fk[1/2, 1/2, 0] + 
         Fk[0, 1/2, 1/2] + Fk[1/2, 0, 1/2])]; Return[Kmat]]; 


MakePlateBMatrix[r0l_, \[Phi]_]:=Block[{x1, y1, z1, x2, y2, z2, x3, y3, z3, x, y, t, 
     Em, \[Nu], \[CapitalPhi]0, \[CapitalPhi], B, \[Kappa], A, J2, H\[CapitalPhi], \[Phi]1,\[Phi]2,\[Phi]3}, 
{\[Phi]1,\[Phi]2,\[Phi]3}=\[Phi];
    {{x1, y1, z1}, {x2, y2, z2}, {x3, y3, z3}} = r0l; 
     J2 = {{(y2 - y3)^2, (y1 - y3)^2, (y1 - y2)^2, 2*(y1 - y3)*(-y2 + y3), 
         -2*(y1 - y2)*(y1 - y3), 2*(y1 - y2)*(y2 - y3)}, {(x2 - x3)^2, (x1 - x3)^2, 
         (x1 - x2)^2, 2*(x1 - x3)*(-x2 + x3), -2*(x1 - x2)*(x1 - x3), 
         2*(x1 - x2)*(x2 - x3)}, {-((x2 - x3)*(y2 - y3)), -((x1 - x3)*(y1 - y3)), 
         -((x1 - x2)*(y1 - y2)), -(x3*(y1 + y2 - 2*y3)) + x2*(y1 - y3) + 
          x1*(y2 - y3), x3*(-y1 + y2) + x1*(2*y1 - y2 - y3) + x2*(-y1 + y3), 
         x3*(y1 - y2) - x2*(y1 - 2*y2 + y3) + x1*(-y2 + y3)}}*(1/A^2); 
     J2[[3]] = 2*J2[[3]];
H\[CapitalPhi]= 
      {{6*(\[Phi]1 + \[Phi]2 + \[Phi]3), 2*(y2*\[Phi]2 + y3*\[Phi]3 - y1*(\[Phi]2 + \[Phi]3)), 
        2*(-(x2*\[Phi]2) - x3*\[Phi]3 + x1*(\[Phi]2 + \[Phi]3)), 0, 0, 0, 0, 0, 0}, 
       {0, 0, 0, 6*(\[Phi]1 + \[Phi]2 + \[Phi]3), 2*(y1*\[Phi]1 + y3*\[Phi]3 - y2*(\[Phi]1 + \[Phi]3)), 
        2*(-(x1*\[Phi]1) - x3*\[Phi]3 + x2*(\[Phi]1 + \[Phi]3)), 0, 0, 0}, 
       {0, 0, 0, 0, 0, 0, 6*(\[Phi]1 + \[Phi]2 + \[Phi]3), 
        2*(y1*\[Phi]1 + y2*\[Phi]2 - y3*(\[Phi]1 + \[Phi]2)), 
        2*(-(x1*\[Phi]1) - x2*\[Phi]2 + x3*(\[Phi]1 + \[Phi]2))}, {2*(3*\[Phi]1 + \[Phi]3), 
        (y3*\[Phi]3 - 2*y1*(2*\[Phi]1 + \[Phi]3) + y2*(4*\[Phi]1 + \[Phi]3))/2, 
        x1*(2*\[Phi]1 + \[Phi]3) + (-(x3*\[Phi]3) - x2*(4*\[Phi]1 + \[Phi]3))/2, 2*(3*\[Phi]2 + \[Phi]3), 
        (y3*\[Phi]3 - 2*y2*(2*\[Phi]2 + \[Phi]3) + y1*(4*\[Phi]2 + \[Phi]3))/2, 
        -(x3*\[Phi]3)/2 + x2*(2*\[Phi]2 + \[Phi]3) - (x1*(4*\[Phi]2 + \[Phi]3))/2, 2*\[Phi]3, 
        ((y1 + y2 - 2*y3)*\[Phi]3)/2, -((x1 + x2 - 2*x3)*\[Phi]3)/2}, 
       {2*\[Phi]1, ((-2*y1 + y2 + y3)*\[Phi]1)/2, ((2*x1 - x2 - x3)*\[Phi]1)/2, 
        2*(\[Phi]1 + 3*\[Phi]2), (y1*\[Phi]1 - 2*y2*(\[Phi]1 + 2*\[Phi]2) + y3*(\[Phi]1 + 4*\[Phi]2))/2, 
        -(x1*\[Phi]1)/2 + x2*(\[Phi]1 + 2*\[Phi]2) - (x3*(\[Phi]1 + 4*\[Phi]2))/2, 
        2*(\[Phi]1 + 3*\[Phi]3), (y1*\[Phi]1 - 2*y3*(\[Phi]1 + 2*\[Phi]3) + y2*(\[Phi]1 + 4*\[Phi]3))/2, 
        -(x1*\[Phi]1)/2 + x3*(\[Phi]1 + 2*\[Phi]3) - (x2*(\[Phi]1 + 4*\[Phi]3))/2}, 
       {2*(3*\[Phi]1 + \[Phi]2), (y2*\[Phi]2 - 2*y1*(2*\[Phi]1 + \[Phi]2) + y3*(4*\[Phi]1 + \[Phi]2))/2, 
        x1*(2*\[Phi]1 + \[Phi]2) + (-(x2*\[Phi]2) - x3*(4*\[Phi]1 + \[Phi]2))/2, 2*\[Phi]2, 
        ((y1 - 2*y2 + y3)*\[Phi]2)/2, -((x1 - 2*x2 + x3)*\[Phi]2)/2, 2*(\[Phi]2 + 3*\[Phi]3), 
        (y2*\[Phi]2 - 2*y3*(\[Phi]2 + 2*\[Phi]3) + y1*(\[Phi]2 + 4*\[Phi]3))/2, 
        -(x2*\[Phi]2)/2 + x3*(\[Phi]2 + 2*\[Phi]3) - (x1*(\[Phi]2 + 4*\[Phi]3))/2}}; 
     Return[Simplify[J2 . H\[CapitalPhi]]]]; 


GetPlateStrain[nplate_, disp_, nodes_, shells_] := 
   Block[{idxg, R, R0, r0l, r0g, dispg, displ, Kel, B, idxx}, 
    idxx = {3, 4, 5, 9, 10, 11, 15, 16, 17}; 
     idxg = Join[6*shells[[nplate]][[1]] - Range[5, 0, -1], 
       6*shells[[nplate]][[2]] - Range[5, 0, -1], 6*shells[[nplate]][[3]] - 
        Range[5, 0, -1]]; r0g = nodes[[shells[[nplate]]]]; 
     R0 = PlateRotationMatrix[r0g]; R = BlockDiagonalMatrix[
       {R0, R0, R0, R0, R0, R0}]; r0l = Simplify[r0g . Transpose[R0]]; 
     dispg = disp[[idxg]]; displ = (R . dispg)[[idxx]]; 
     Return[{MakePlateBMatrix[r0l, {1, 0, 0}] . displ, 
       MakePlateBMatrix[r0l, {0, 1, 0}] . displ, MakePlateBMatrix[r0l, {0, 0, 1}] . 
        displ}]; ]; 


GetPlateStress[nplate_, disp_, nodes_, shells_, prop_] := 
   Block[{\[Epsilon],\[Kappa],Em,\[Nu],t}, \[Epsilon] = GetPlateStrain[nplate, disp, nodes, shells]; 
     {Em, \[Nu],t} = prop; \[Kappa] = Simplify[(Em/(1 - \[Nu]^2))*
        {{1, \[Nu], 0}, {\[Nu], 1, 0}, {0, 0, (1 - \[Nu])/2}}]; 
     Return[Simplify[t/2 \[Epsilon] . \[Kappa]]]]; 


GetPlateVMStress[nplate_, disp_, nodes_, shells_, prop_] := 
  Block[{\[Sigma]}, \[Sigma] = GetPlateStress[nplate, disp, nodes, shells, prop]; 
    Return[{CalcVMStress[\[Sigma][[1]]], CalcVMStress[\[Sigma][[2]]], CalcVMStress[\[Sigma][[3]]]}]]


GetPlate\[Tau]maxStress[nplate_, disp_, nodes_, shells_, prop_] := 
  Block[{\[Sigma]}, \[Sigma] = GetPlateStress[nplate, disp, nodes, shells, prop]; 
    Return[{Calc\[Tau]maxStress[\[Sigma][[1]]], Calc\[Tau]maxStress[\[Sigma][[2]]], Calc\[Tau]maxStress[\[Sigma][[3]]]}]]


GetPlateMainStress[nplate_, disp_, nodes_, shells_, prop_] := 
  Block[{\[Sigma]}, \[Sigma] = GetPlateStress[nplate, disp, nodes, shells, prop]; 
    Return[{CalcMohrCircle[\[Sigma][[1]]], CalcMohrCircle[\[Sigma][[2]]], CalcMohrCircle[\[Sigma][[3]]]}]]


(* ::Subsection::Closed:: *)
(*Shells*)


MakeShellBMatrix[r0l_]:=Block[{x1, y1, z1, x2, y2, z2, x3, y3, z3, x, y,\[CapitalPhi]0, \[CapitalPhi], B},
     {{x1, y1, z1}, {x2, y2, z2}, {x3, y3, z3}} = r0l;
     \[CapitalPhi]0[x_, y_] := {1, x, y}; \[CapitalPhi][x_, y_] := 
      Simplify[\[CapitalPhi]0[x, y] . Inverse[{\[CapitalPhi]0[x1, y1], \[CapitalPhi]0[x2, y2], \[CapitalPhi]0[x3, y3]}]]; 
     B = ConstantArray[0, {3, 6}]; B[[1,{1, 3, 5}]] = Simplify[D[\[CapitalPhi][x, y], x]]; 
     B[[2,{2, 4, 6}]] = Simplify[D[\[CapitalPhi][x, y], y]]; B[[3,{1, 3, 5}]] = 
      Simplify[D[\[CapitalPhi][x, y], y]]; B[[3,{2, 4, 6}]] = Simplify[D[\[CapitalPhi][x, y], x]]; 
     Return[FullSimplify[B]]];


MakeShellKll[r0l_, prop_] := Block[{x1, y1, z1, x2, y2, z2, x3, y3, z3, x, y, t, 
     Em, \[Nu], \[CapitalPhi]0, \[CapitalPhi], B, Kmat, Area}, {Em, \[Nu], t} = prop; 
     {{x1, y1, z1}, {x2, y2, z2}, {x3, y3, z3}} = r0l; 
     Area = (1/2)*Norm[Cross[r0l[[2]] - r0l[[1]], r0l[[3]] - r0l[[1]]]]; 
     \[CapitalPhi]0[x_, y_] := {1, x, y}; \[CapitalPhi][x_, y_] := 
      Simplify[\[CapitalPhi]0[x, y] . Inverse[{\[CapitalPhi]0[x1, y1], \[CapitalPhi]0[x2, y2], \[CapitalPhi]0[x3, y3]}]]; 
     B = ConstantArray[0, {3, 6}]; B[[1,{1, 3, 5}]] = Simplify[D[\[CapitalPhi][x, y], x]]; 
     B[[2,{2, 4, 6}]] = Simplify[D[\[CapitalPhi][x, y], y]]; B[[3,{1, 3, 5}]] = 
      Simplify[D[\[CapitalPhi][x, y], y]]; B[[3,{2, 4, 6}]] = Simplify[D[\[CapitalPhi][x, y], x]]; 
     Kmat = Simplify[Area*((t*Em)/(1 - \[Nu]^2))*{{1, \[Nu], 0}, {\[Nu], 1, 0}, 
         {0, 0, (1 - \[Nu])/2}}]; Return[FullSimplify[Transpose[B] . Kmat . B]]]; 


Shell3DStiffness[r0g_, prop_] := Block[{R, R0, r0l, Kll, idxx, K}, 
    R0 = PlateRotationMatrix[r0g]; R = BlockDiagonalMatrix[
       {R0, R0, R0, R0, R0, R0}]; r0l = Simplify[r0g . Transpose[R0]]; 
     K = ConstantArray[0, {18, 18}]; idxx = {1, 2, 7, 8, 13, 14}; 
     K[[idxx,idxx]] = MakeShellKll[r0l, prop]; 
     Return[Simplify[Transpose[R] . K . R]]]; 


GetShellStrain[nplate_,disp_,nodes_,shells_]:=Block[{idxg,R,R0,r0l,r0g,dispg,displ,Kel,B,idxx},
idxx = {1, 2, 7, 8, 13, 14}; 
idxg = Join[6*shells[[nplate]][[1]] - Range[5, 0, -1], 6*shells[[nplate]][[2]] - 
           Range[5, 0, -1], 6*shells[[nplate]][[3]] - Range[5, 0, -1]];
r0g = nodes[[shells[[nplate]]]];
R0 = PlateRotationMatrix[r0g];
R = BlockDiagonalMatrix[{R0, R0, R0, R0, R0, R0}];
r0l = Simplify[r0g.Transpose[R0]];
dispg=disp[[idxg]];displ=R.dispg;
B= MakeShellBMatrix[r0l];
 Return[Simplify[B.displ[[idxx]]]];];


GetShellStress[nplate_, disp_, nodes_, shells_, prop_] := 
  Block[{\[Epsilon], Em, \[Nu], Kmat}, 
   {Em, \[Nu]} = prop; \[Epsilon]=GetShellStrain[nplate, disp, nodes, shells];
Kmat = (Em/(1 - \[Nu]^2))*{{1, \[Nu], 0}, {\[Nu], 1, 0}, 
       {0, 0, (1 - \[Nu])/2}};
 Return[Simplify[\[Epsilon].Kmat]];];


GetShellVMStress[nplate_, disp_, nodes_, shells_, prop_] := 
  Block[{\[Sigma]}, \[Sigma] = GetShellStress[nplate, disp, nodes, shells, prop]; 
    Return[CalcVMStress[\[Sigma]]]];


GetShell\[Tau]maxStress[nplate_, disp_, nodes_, shells_, prop_] := 
  Block[{\[Sigma]}, \[Sigma] = GetShellStress[nplate, disp, nodes, shells, prop]; 
    Return[Calc\[Tau]maxStress[\[Sigma]]]];


GetShellMainStress[nplate_, disp_, nodes_, shells_, prop_] := 
  Block[{\[Sigma]}, \[Sigma] = GetShellStress[nplate, disp, nodes, shells, prop]; 
    Return[CalcMohrCircle[\[Sigma]]]];


(* ::Subsection:: *)
(*Graphic Funcs*)


DrawUnitCell[nodes_, beams_, plates_] := 
   Block[{ucline, kk, nbeams, nfaces, ucfaces}, nbeams = Dimensions[beams][[1]]; 
      nfaces = Dimensions[plates][[1]]; ucline = {Line[nodes[[beams[[1]]]]]}; 
      For[kk = 2, kk <= nbeams, kk++, ucline = Append[ucline, 
            Line[nodes[[beams[[kk]]]]]]]; ucfaces = {Polygon[nodes[[plates[[1]]]]]}; 
      For[kk = 2, kk <= nfaces, kk++, ucfaces = Append[ucfaces, 
            Polygon[nodes[[plates[[kk]]]]]]]; Graphics3D[{EdgeForm[], Yellow, 
          Opacity[0.25], ucfaces, Thickness[0.05], Red, Opacity[1], ucline}, 
        {Axes -> False, Boxed -> False, AspectRatio -> 1}]]


EndPackage[]
