(* Mathematica code *)
path[zeta_,p_,q_] := Module[{angle},
    angle = ArcCos[p.q];
    (Sin[angle-zeta]*p + Sin[zeta]*q)/Sin[angle]
];
geodesic[p_, q_] := ParametricPlot3D[path[zeta,p,q], {zeta, 0, ArcCos[p.q]}];
pts = Normalize/@RandomReal[{-1,1},{3,3}];
(* Ensure that the vertices come in the right order *)
If[Det[pts]<0, pts = {pts[[1]],pts[[3]],pts[[2]]}];
(* Normals to planes containing pairs of vertices plus the origin *)
n1 = Cross[pts[[2]], pts[[3]]];
n2 = Cross[pts[[3]], pts[[1]]];
n3 = Cross[pts[[1]], pts[[2]]];
(* Unit sphere plot *)
sp = Graphics3D[{Opacity[0.5], Sphere[]}];
(* Geodesic curve plots *)
g12p = geodesic[pts[[1]], pts[[2]]];
g23p = geodesic[pts[[2]], pts[[3]]];
g31p = geodesic[pts[[3]], pts[[1]]];
(* Vertex plots *)
pp = ListPointPlot3D[pts, PlotStyle->PointSize[Medium]];
(* Spherical triangle described implicitly and a plot of it *)
st = ImplicitRegion[x^2+y^2+z^2==1 && {x,y,z}.n1 >= 0 && {x,y,z}.n2 >=0 && {x,y,z}.n3 >=0, {x,y,z}];
dst = DiscretizeRegion[st];
Show[sp,g12p,g23p,g31p,dst,pp]
(* Area of the spherical triangle via Mathematica's adaptive numerical integration *)
area1 = Integrate[1, Element[{x,y,z}, st]];
(* Area of the spherical triangle via Gauss' formula *)
angle[o_,p_,q_]:= Module[{po, pq},
    po = Normalize[(o-p) - ((o-p).p)*p];
    pq = Normalize[(q-p) - ((q-p).p)*p];
    ArcCos[po.pq]
];
angle1 = angle[pts[[3]],pts[[1]],pts[[2]]];
angle2 = angle[pts[[1]],pts[[2]],pts[[3]]];
angle3 = angle[pts[[2]],pts[[3]],pts[[1]]];
area2 = angle1 + angle2 + angle3 - Pi;
{area1, area2}
(* Integral of first Cartesian coordinate via Mathematica's adaptive numerical integration *)
Ix1 = Integrate[x, Element[{x,y,z}, st]];
(* Integral of first Cartesian coordinate using my own formula *)
IxTerm[p_,q_] := Module[{cosine},
    cosine = p.q;
    1/2*ArcCos[cosine]/Sqrt[1-cosine^2]*(q[[3]]*p[[2]]-q[[2]]*p[[3]])
];
Ix2 = IxTerm[pts[[1]],pts[[2]]] + IxTerm[pts[[2]],pts[[3]]] + IxTerm[pts[[3]],pts[[1]]];
{Ix1, Ix2}
