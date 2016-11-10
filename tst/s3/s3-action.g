# G:= SymmetricGroup(3);
G:= Image(IsomorphismPermGroup(SmallGroup(6, 1)));
SetName(G, "s3");

# sections cache;
AllSections:= [];
secs:= RepresentativesSections(G);
for grp in Set(List(secs, sec-> OneMorphismSection(sec)!.group)) do
    sss:= RepresentativesSections(grp);
    sss:= Concatenation(List(sss, x-> Orbit(grp, x)));
    Append(AllSections, sss);
od;

idss:= Set(List(secs, IdSection));
secs_by_id:= List(idss, x-> []);
for sec in secs do
    Add(secs_by_id[Position(idss, IdSection(sec))], sec);
od;

secs:= Concatenation(secs_by_id);
secs:= Concatenation(List(secs, x-> Orbit(G, x)));

#  make really all morphisms, one per automorphism.
AllMorphismsSection0:= function(sec)
    local   mor,  U,  A;

    mor:= OneMorphismSection(sec);
    U:= mor!.group;
    
    A:= AutomorphismGroup(U);
##  FIXME: use mor * auto
    return List(Elements(A), d-> Morphism(sec, U, mor!.theta*d));
end;

GelfandModule:= function(G)
    local   secs,  idss,  secs_by_id,  sec,  mors;
    secs:= RepresentativesSections(G);
    idss:= Set(List(secs, IdSection));
    secs_by_id:= List(idss, x-> []);
    for sec in secs do
        Add(secs_by_id[Position(idss, IdSection(sec))], sec);
    od;
    secs:= Concatenation(secs_by_id);
    secs:= Concatenation(List(secs, x-> Orbit(G, x)));
    mors:= Concatenation(List(secs, x-> AllMorphismsSection0(x)));
    return List(mors, x-> Opposite(TripleMorphism(x)));
end;

# a basis of morphisms U -> P/K, represented as triples.
basis:= GelfandModule(G);

# list all subsgroups
## FIXME: conjugation of triples doesn't work???
GG:= DirectProduct(G, G);
trips:= Flat(List(TriplesDirectProduct(G, G), x-> x.trip));
subs:= List(trips, SubgroupTriple);
subs:= Concatenation(List(subs, x-> (Orbit(GG, x))));
trips:= List(subs, x-> TripleSubgroup(GG, x));


GelfandImage:= function(b, t)
    local   co,  pos,  pre;
    co:= CompositionTriples(b, t);
    pos:= Position(AllSections, Sections(co)[1]);
    pre:= TripleMorphism(OneMorphismSection(AllSections[pos]));
    return Position(basis, CompositionTriples(Opposite(pre), co));
end;

GelfandMatrix:= function(basis, t)
    local   mat,  b,  co,  kappa,  pos,  pre,  kappa1,  row;
    mat:= [];
    for b in basis do
        co:= CompositionTriples(b, t);
        kappa:= KappaTriples(b, t);
        if Sections(b)[2] <> Sections(t)[1] then
            kappa:= 0;
        fi;
        pos:= Position(AllSections, Sections(co)[1]);
        pre:= Opposite(TripleMorphism(OneMorphismSection(AllSections[pos])));
        pos:= Position(basis, CompositionTriples(pre, co));
        kappa1:= KappaTriples(pre, co);
        kappa2:= KappaTriples(pre, b);
        row:= List(basis, x-> 0);
        row[pos]:= kappa;
        Add(mat, row);
    od;
    return mat;
end;

# the matrices!
mats:= List(trips, t-> GelfandMatrix(basis, t));

# now apply zeta functions ...
comsec:= List(secs, x-> List(secs, y-> Iverson(IsCompSubsection(x, y))));
smpsec:= List(secs, x-> List(secs, y-> Iverson(TopSec(x) = TopSec(y) and IsCompSubsection(x, y))));
smksec:= List(secs, x-> List(secs, y-> Iverson(BotSec(x) = BotSec(y) and IsCompSubsection(x, y))));
smusec:= smksec^-1 * comsec * smpsec^-1;
dessec:= List(secs, x-> List(secs, y-> Iverson(IsDescendantSection(x, y))));

butts:= List(secs, x-> List(Filtered(secs, y-> IsDescendantSection(x, y)), z-> Butterfly(x, z)));


desmor:= [];
for b in basis do
    pos:= Position(secs, Sections(b)[2]);
    co:= List(butts[pos], x-> CompositionTriples(b, x));
    lis:= List(co, x-> Position(basis, CompositionTriples(Opposite(TripleMorphism(OneMorphismSection(AllSections[Position(AllSections, Sections(x)[1])]))), x)));
    Add(desmor, lis);
od;

smkmor:= List([1..Length(basis)], i-> Filtered(desmor[i], j-> BotSec(Sections(basis[j])[2]) = BotSec(Sections(basis[i])[2])));
smtmor:= List([1..Length(basis)], i-> Filtered(desmor[i], j-> TopSec(Sections(basis[j])[2]) = TopSec(Sections(basis[i])[2])));

desmor:= List(desmor, x-> BlistList([1..Length(basis)], x));
desmor:= List(desmor, x-> List(x, Iverson));

smkmor:= List(smkmor, x-> BlistList([1..Length(basis)], x));
smkmor:= List(smkmor, x-> List(x, Iverson));

smtmor:= List(smtmor, x-> BlistList([1..Length(basis)], x));
smtmor:= List(smtmor, x-> List(x, Iverson));

smumor:= smkmor^-1 * desmor * smtmor^-1;
