#############################################################################
##
##  boucindex.g
##

##  The Bouc index of a section P/K is the number
##
##    m_{P,K} = \frac1{|P|} \sum_{U: KU = P} \mu(P, U) |U|
##
##  Note that KU = P implies |P||K \cap U| = |K||U|, so
##
##    m_{P,K} = \frac1{|K|} \sum_{U: KU = P} \mu(P, U) |K \cap U|
##
##  * Does this have a name in lattice theory?
##
##  * Are there efficient ways to compute it?
##
##  The following is a brute force implementation.
##
BoucIndexSection:= function(section)
    local   P,  K,  subs,  zeta,  moebius,  elts,  list,  xxx,  m,  i,
            sub;

    P:= TopSec(section);
    K:= BotSec(section);

    # find all subgroups of P that supplement K
    subs:= Concatenation(List(ConjugacyClassesSubgroups(P), Elements));
    zeta:= List(subs, x-> List(subs, y-> Iverson(IsSubgroup(x, y))));
    moebius:= zeta^-1;
    moebius:= moebius[Length(moebius)];
    elts:= Elements(P);
    list:= [];
    xxx:= List(elts, x-> 0);
    m:= 0;
    for i in [1..Length(subs)] do
        sub:= subs[i];
        if ClosureGroup(K, sub) = P then
            Add(list, sub);
            xxx:= xxx + moebius[i] * List(elts, x-> Iverson(x in sub));
            m:= m + moebius[i] * Size(sub);
        fi;
    od;

    return rec(list:= list, m:= m, xxx:= xxx);
end;

BoucIndexNormalSubgroup:= function(P, A)
    local   ct,  ccs,  a,  tom,  mat,  l,  wt,  pc,  mun,  inv,  moe,
            ppp,  cf;

    ct:= CharacterTable(P);
    ccs:= ConjugacyClassesSubgroups(P);
    a:= PositionProperty(ccs, c-> A in c);
    tom:= TableOfMarks(P);
    mat:= MatTom(tom);
    l:= Length(mat);
    wt:= List([1..l], i-> mat[i][i]);
    pc:= PermCharsTom(ct, tom);
    mun:= List([1..l], i-> mat[i]/wt[i]);
    inv:= mun^-1;
    moe:= inv[l];
    ppp:= Filtered([1..l], i-> moe[i] <> 0 and ScalarProduct(pc[i], pc[a]) = 1);
    cf:= Sum(ppp, i-> moe[i] * pc[i] / wt[i]);
#Error();
    return cf;
end;

LambdaNormalSubgroup:= function(P, K)
    local   bin;

    bin:= BoucIndexNormalSubgroup(P, K);
    return Size(K) * ScalarProduct(TrivialCharacter(P), bin);
end;

##  FIXME: the following is wrong.  It probably needs to be weighted
##  with some normalizer indices ....
LambdasSections:= function(G)
    local   iso,  list,  i,  mat,  rep,  nor,  nor1,  sz;

    iso:= IsoClassIncMatSections(G);
    list:= [];
    for i in [1..Length(iso.mats)] do
        mat:= iso.mats[i];
        rep:= iso.reps[i];
        nor:= List(rep, x-> Size(NormalizerSection(x)));
        nor1:= List(nor, x-> 1/x);
        sz:= List(rep, x-> Size(BotSec(x)));
        Append(list, HadamardProduct(nor, (mat^-1 * HadamardProduct(nor1, sz))));
    od;
    return list;
end;
