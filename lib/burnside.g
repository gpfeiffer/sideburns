#############################################################################
##
##  burnside.g
##
##  towards the double burnside ring...
##


###############################################################################
##
##  MackeyReps(L, M)
##
##  Given $L$ a subgroup of $G \times H$ and $M$ a subgroup
##  of $H \times K$,
##  this function computes the Mackey product of $L$ and $M$
##  as a linear combination
##  over double coset representatives of the top of the
##  second and first section of $L$
##  and $M$ respectively in $H$ of subgroups in $G \times K$.
##
##  given two triples L, M, find double coset reps for Mackey formula
##
MackeyReps:= function(L, M)
    local   H,  emb,  P2,  P1,  D;

    H:= Sections(L)[2]!.G;
    emb:= Embedding(ProductGroup(L), 2);
    P2:= TopSec(Sections(L)[2]);
    P1:= TopSec(Sections(M)[1]);
    D:= List(DoubleCosets(H, P2, P1), Representative);

    return List(D, h-> CompositionTriples(L^(h^emb), M));
end;


############################################################################
##
##  BasisDoubleBurnsideRing(G)
##
##  Constructs the right regular representation of the double Burnside
##  algebra of $G$ over the rationals using
##  the set of conjugacy classes  of subgroups of $G \times G$ as a basis.
##
BasisDoubleBurnsideRing:= function(G)
    local   classes,  qqq,  iii,  mats,  i,  mat,  j,  row,  sub,  k;

    qqq:= Concatenation(List(TriplesDirectProduct(G, G), x-> Flat(x.trip)));
    classes:= ConjugacyClassesSubgroupsDirectProduct(G, G);
    iii:= [1..Length(qqq)];

    # mackey
    mats:= [];
    for i in iii do
        mat:= [];  Add(mats, mat);
        for j in iii do
            row:= 0 * iii;  Add(mat, row);
            for sub in Collected(MackeyReps(qqq[j], qqq[i])) do
                k:= PositionProperty(classes, c-> sub[1] in c);
                row[k]:= row[k] + sub[2];
            od;
        od;
        Print(".\c");
    od;
    Print("\n");

    return rec(
               algebra:= AlgebraWithOne(Rationals, mats, "basis"),
               basis:= mats,
               );
end;


#############################################################################
##
##  GelfandModule( G )
##
##  compute a basis for the Gelfand module ...
##

##  FIXME: use all sections!
GelfandModule:= function(G)
    local   secs,  mors;
    secs:= RepresentativesSections(G);
    mors:= Concatenation(List(secs, x-> AllMorphismsSection(x)));
    return List(mors, x-> Opposite(TripleMorphism(x)));
end;