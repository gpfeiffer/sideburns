#############################################################################
##
##  zgroups.g
##
ZGroup:= function(m, n, r)
    local   A,  a,  B,  b,  alpha,  aut,  phi;

    # check conditions first
    if GcdInt(n, m) > 1 then
        Error("gcd(n, m) must be 1");
    fi;
    if GcdInt(r-1, m) > 1 then
        Error("gcd(r-1, m) must be 1");
    fi;
    if PowerModInt(r, n, m) <> 1 then
        Error("r^n must be 1 mod m");
    fi;


    # construct the group as semidirect product
    A:= CyclicGroup(IsPermGroup, m);
    a:= GeneratorsOfGroup(A)[1];
    B:= CyclicGroup(IsPermGroup, n);
    b:= GeneratorsOfGroup(B)[1];
    alpha:= GroupHomomorphismByImages(A, A, [a], [a^r]);
    aut:= Group(alpha);
    phi:= GroupHomomorphismByImages(B, aut, [b], [alpha]);
    return SemidirectProduct(B, phi, A);
end;


# given n and m, find all possible r
ZParams:= function(m, n)
    local   list,  r;
    list:= [];
    if GcdInt(m, n) > 1 then
        return list;
    fi;
    for r in [1..m-1] do
        if GcdInt(r-1, m) = 1 and PowerMod(r, n, m) = 1 then
            Add(list, r);
        fi;
    od;
    return list;
end;

ZParamsSize:= function(mn)
    local   zzz,  m,  n,  r;
    zzz:= [];
    for m in DivisorsInt(mn) do
        n:= mn/m;
        for r in ZParams(m, n) do
            Add(zzz, [m, n, r]);
        od;
    od;
    return zzz;
end;

ZParamsSizeUnique:= function(mn)
    local   zzz,  m,  n,  new,  ids,  r,  idx;
    zzz:= [];
    for m in DivisorsInt(mn) do
        n:= mn/m;
        new:= [];  ids:= [];
        for r in ZParams(m, n) do
            Add(new, [m, n, r]);
            Add(ids, IdGroup(ZGroup(m, n, r)));
        od;
        idx:= List(Set(ids), x-> Position(ids, x));
        Append(zzz, new{idx});
    od;
    return zzz;
end;

OutZGroupParams:= function(param)
    local   G,  A,  I;
    G:= ZGroup(param[1], param[2], param[3]);
    A:= AutomorphismGroup(G);
    I:= InnerAutomorphismsAutomorphismGroup(A);
    return A/I;
end;

#############################################################################
##
##  compute a couple of regular representations of B(G, G)
##  where G is a Z-group.
##
BaseChangesDoubleBurnsideZ:= function(G)
    local   trips,  reps,  sec1s,  ccs,  GG,  top,  bot,  iso,  tops,
            bots,  isos,  topm,  botm,  isom,  diam,  ll,  N,  tap,
            potm,  i,  bas,  chg,  a,  wt1,  wt1m,  wt2,  wt2m,  b,
            c,  d;

    ##  its conjugacy classes of subgroups
    trips:= Flat(List(TriplesDirectProduct(G, G), x-> x.trip));
    reps:= List(trips, SubgroupTriple);
    sec1s:= List(trips, x-> Sections(x)[1]);
    ccs:= ConjugacyClassesSubgroupsDirectProduct(G, G);
    GG:= DirectProduct(G, G);

    ##  decomposition of mat
    top:= TopClassIncMatDirectProduct(G, G);
    bot:= BotClassIncMatDirectProduct(G, G);
    iso:= IsoClassIncMatDirectProduct(G, G);
    reps:= Concatenation(top.reps);
    tops:= List(ccs, x-> PositionProperty(reps, r-> r in x));
    reps:= Concatenation(bot.reps);
    bots:= List(ccs, x-> PositionProperty(reps, r-> r in x));
    reps:= Concatenation(iso.reps);
    isos:= List(ccs, x-> PositionProperty(reps, r-> r in x));
    topm:= DirectSumMat(top.mats){tops}{tops};
    botm:= DirectSumMat(bot.mats){bots}{bots};
    isom:= DirectSumMat(iso.mats){isos}{isos};
    diam:= List(ccs, Representative);
    diam:= List(diam, x-> Index(Normalizer(GG, x), x));
    diam:= DiagonalMat(diam/Size(G));

    # transpose topm:
    ll:= List(ccs, Size);
    N:= Length(topm);
    tap:= List([1..N], i-> List([1..N], j-> topm[i][j] * ll[j] / ll[i]));
    tap:= TransposedMat(tap);;

    # ... and neutralize noncyclic part:
    potm:= MutableCopyMat(tap);;
    for i in [1..N] do
        if not IsCyclic(AsGroup(Sections(trips[i])[1])) then
            potm[i]{[1..i-1]}:= 0*[1..i-1];
        fi;
    od;

    # Burnside
    bas:= BasisDoubleBurnsideRing(G);

    chg:= diam * botm * isom;;
    a:= RightRegularBaseChange(bas.basis, chg);;

    ###  this is the number of G-conjugates of (P_1,K_1) -> U:
    ###  the index of N_G(P_1, K_1) in G  x  the size of Aut_{\theta_1}(U)
    wt1:= List(sec1s, x-> Index(G, NormalizerSection(x))
           * Size(Conjugators(OneMorphismSection(x))));

    wt1m:= DiagonalMat(wt1);

    ### this is |Aut(P_1)| / |Aut(P_1/K_1)|
    wt2:= List(sec1s, x-> Size(AutomorphismGroup(TopSec(x)))
               / Size(AutomorphismGroup(AsGroup(x))));

    wt2m:= DiagonalMat(wt2);

    chg:= diam * botm * isom * wt2m;;
    b:= RightRegularBaseChange(bas.basis, chg);;

    chg:= diam * botm * isom * wt2m * potm;;
    c:= RightRegularBaseChange(bas.basis, chg);;

    chg:= diam * botm * isom * wt2m * potm * wt1m;;
    d:= RightRegularBaseChange(bas.basis, chg);;

    return rec(
               diam:= diam,
               botm:= botm,
               isom:= isom,
               topm:= topm,
               potm:= potm,
               wt1m:= wt1m,
               wt2m:= wt2m,
               a:= a,
               b:= b,
               c:= c,
               d:= d,
               chg:= chg
               );
end;


#############################################################################
##
##  identify columns and shrink a regular rep'n z accordingly.
##
ShrinkByCols:= function(z, cols)
    local   shrink,  col,  i,  j,  shrink1,  firsts,  seconds,  f,  m,
            h;

    shrink:= z[1]^0;;
    for col in cols do
        i:= col[1];
        for j in col{[2..Length(col)]} do
            shrink[j][i]:= -1;
        od;
    od;
    shrink1:= shrink^-1;;

    firsts:= List(cols, x-> x[1]);;
    seconds:= Difference([1..Length(z)], firsts);;

    f:= List(z, x-> x^shrink1);;
    m:= List(f, x-> x{firsts}{firsts});;
    h:= List(f, x-> x{seconds}{firsts});;

    return rec(m:= m, h:= h);
end;
