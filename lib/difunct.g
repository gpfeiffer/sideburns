#############################################################################
##
##  Difunctional Relations
##

#############################################################################
##
##  Sub-Partitions
##
##  For now, a subpartition (of [1..n]) is simply a list of disjoint lists,
##  or rather a set of disjoint subsets of [1..n].
##
##  We will use *three* different representations:
##
##  + parts: the disjoint subsets (lacks knowledge on n)
##
##  + points: an image list of the binary relation
##
##  + equivalence: a GAP binary relation
##

#############################################################################
##
##  SetComposition
##
##  Given a composition lambda of m, construct our favourite set partition
##  of type lambda.
##
SetComposition:= function(lambda)
    local   sum,  list,  l;
    sum:= 0;  list:= [];
    for l in lambda do
        Add(list, sum + [1..l]);
        sum:= sum + l;
    od;
    return list;
end;

TypeParts:= function(list)
    return Reversed(SortedList(List(list, Length)));
end;

#############################################################################
##
##  how to turn a parts into a points into a binary relation
##
PointsParts:= function(n, parts)
    local   points,  part,  i;
    points:= List([1..n], i-> []);
    for part in parts do
        for i in part do
            points[i]:= part;
        od;
    od;
    return points;
end;

PartsPoints:= function(points)
    return Difference(points, [[]]);
end;

EquivalencePoints:= BinaryRelationOnPoints;

PointsEquivalence:= Successors;

EquivalenceParts:= function(n, parts)
    return EquivalencePoints(PointsParts(n, parts));
end;

PartsEquivalence:= function(bin)
    return PartsPoints(PointsEquivalence(bin));
end;

#############################################################################
##
##  all equivalences on n points
##
AllEquivalencesType:= function(n, k)
    local   list,  sym,  m,  lambda;
    list:= [];
    sym:= SymmetricGroup(n);
    for m in [0..n] do
        for lambda in Partitions(m, k) do
            Append(list, Orbit(sym, SetComposition(lambda), OnSetsSets));
        od;
    od;
    return List(list, x-> EquivalenceParts(n, x));
end;

AllEquivalences:= function(n)
    return Concatenation(List([0..n], k-> AllEquivalencesType(n, k)));
end;


RandomEquivalence:= function(n)
    local   lambda,  parts;
    lambda:= Random(Union(List([0..n], Partitions)));
    parts:= OnSetsSets(SetComposition(lambda), Random(SymmetricGroup(n)));
    return EquivalenceParts(n, parts);
end;

#############################################################################
##
##  restriction to a subset
##
RestrictedParts:= function(set, parts)
    return Difference(List(parts, part-> Intersection(part, set)), [[]]);
end;

#############################################################################
##
##  Support
##
##  the *support* of a partition is the union of its parts
##
SupportParts:= Union;
SupportPoints:= Union;
SupportEquivalence:= bin -> SupportPoints(PointsEquivalence(bin));

LargestPointParts:= function(parts)
    if parts = [] then  return 1;  fi;
    return Maximum(List(parts, Maximum));
end;

#############################################################################
IntersectionParts:= function(l, r)
    local   n,  ptsL,  ptsR;
    n:= Maximum(List([l, r], LargestPointParts));
    ptsL:= PointsParts(n, l);
    ptsR:= PointsParts(n, r);
    return Difference(List([1..n], i-> Intersection(ptsL[i], ptsR[i])), [[]]);
end;


#############################################################################
##
##  Butterfly
##
ButterflyParts:= function(l, r)
    local   suppL,  suppR,  n,  relL,  relR,  fly;
    suppL:= SupportParts(l);
    suppR:= SupportParts(r);
    n:= Maximum(List([l, r], LargestPointParts));
    relL:= EquivalenceParts(n, RestrictedParts(suppR, l));
    relR:= EquivalenceParts(n, RestrictedParts(suppL, r));
    fly:= TransitiveClosureBinaryRelation(relL + relR);
    return PartsEquivalence(fly);
end;


#############################################################################
##
##  DescentMapPartitions
##
##  check condition [x]' \supseteq [x] \cap P' for all x \in P'
##  and then return the connecting binary relation ...
##  So here hi and lo are binary relations representing partitions.
##
DescentMapEquivalences:= function(hi, lo)
    local   partsHi,  partsLo,  suppHi,  suppLo,  part,  n;

    partsHi:= PartsEquivalence(hi);
    partsLo:= PartsEquivalence(lo);
    suppHi:= SupportParts(partsHi);
    suppLo:= SupportParts(partsLo);

    # check
    for part in partsLo do
        if not IsSubset(part, Intersection(part^hi, suppLo)) then
            return fail;
        fi;
    od;

    #  build relation
    n:= Length(Successors(hi));
    return Difunction(n, partsLo, List(partsLo, x-> x^hi), ());
end;


#############################################################################
##
##  RightStar
##
RightStar:= function(all, r)
    local   n,  mat,  i;
    n:= Length(all);
    mat:= NullMat(n, n);
    for i in [1..n] do
        mat[i][Position(all, ButterflyParts(all[i], r))]:= 1;
    od;
    return mat;
end;


#############################################################################
##
##  Now for the difuncs ...
##
##  A difunc is three things: two subpartitions, and a permutation
##  describing a bijection between them.
##
##  We'll try and implement them as GAP binary relations
##
Difunction:= function(n, parts1, parts2, perm)
    local   points,  p,  i;
    points:= List([1..n], i-> []);
    for p in [1..Length(parts1)] do
        for i in parts1[p] do
            points[i]:= parts2[p^perm];
        od;
    od;
    return BinaryRelationOnPoints(points);
end;

RandomDifunction:= function(n, m)
    local   sn,  sm,  p;
    sn:= SymmetricGroup(n);
    sm:= SymmetricGroup(m);
    p:= List([1, 2], i->
     OnSetsSets(SetComposition(Random(RestrictedPartitions(n, [1..n], m))), Random(sn)));
    return Difunction(n, p[1], p[2], Random(sm));
end;

OppositeDifunction:= function(dif)
    return BinaryRelationOnPoints(Successors(InverseGeneralMapping(dif)));
end;

#############################################################################
##
##  PartsDifunction
##
##  A difunction is a bijection between two equivalences.  Extract them!
##
RightParts:= function(dif)
    return Difference(Successors(dif), [[]]);
end;

LeftParts:= function(dif)
    return RightParts(OppositeDifunction(dif));
end;

PartsDifunction:= function(dif)
    return [LeftParts(dif), RightParts(dif)];
end;

CoRestrictionDifunction:= function(dif, parts)
    local   n,  hi,  lo,  des;
    n:= Length(Successors(dif));
    hi:= EquivalenceParts(n, RightParts(dif));
    lo:= EquivalenceParts(n, parts);
    des:= DescentMapEquivalences(hi, lo);
    return dif * des;
end;

RestrictionDifunction:= function(dif, parts)
    return OppositeDifunction(CoRestrictionDifunction(OppositeDifunction(dif), parts));
end;

ProductDifunctions:= function(l, r)
    local   fly;
    fly:= ButterflyParts(RightParts(l), LeftParts(r));
    return CoRestrictionDifunction(l, fly) * RestrictionDifunction(r, fly);
end;


AllDifunctionsType:= function(n, k)
    local   eq,  sym,  all;
    eq:= List(AllEquivalencesType(n, k), PartsEquivalence);
    sym:= SymmetricGroup(k);
    all:= List(eq, x-> List(eq, y-> List(sym, a-> Difunction(n, x, y, a))));
    return Concatenation(Concatenation(all));
end;

AllDifunctions:= function(n)
    return Concatenation(List([0..n], k-> AllDifunctionsType(n, k)));
end;

RightStar:= function(all, r)
    local   n,  mat,  i;
    n:= Length(all);
    mat:= NullMat(n, n);
    for i in [1..n] do
        mat[i][Position(all, ProductDifunctions(all[i], r))]:= 1;
    od;
    return mat;
end;

#############################################################################
##
##  Test
##
parts1:= [[1],[3,4]];
parts2:= [[1],[3],[4]];
parts3:= [[1,2],[3],[4,5]];
parts4:= [[1,2],[3],[4,5],[6]];
e1:= EquivalenceParts(6, parts1);
e2:= EquivalenceParts(6, parts2);
e3:= EquivalenceParts(6, parts3);
e4:= EquivalenceParts(6, parts4);
