##  the group
G:= SymmetricGroup(3);
SetName(G, "S3");
GG:= DirectProduct(G, G);

##  its conjugacy classes of subgroups
trips:= Flat(List(TriplesDirectProduct(G, G), x-> x.trip));
reps:= List(trips, SubgroupTriple);
sec1s:= List(trips, x-> Sections(x)[1]);
ccs:= ConjugacyClassesSubgroupsDirectProduct(G, G);

##  its table of marks (just to be sure).
tom:= TableOfMarks(GG);
mat:= MatTom(tom);

poss:= List(ccs, x-> PositionProperty(ConjugacyClassesSubgroups(GG), c-> Representative(x) in c));

mat:= mat{poss}{poss};

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
diam:= DiagonalMat(diam);
tomm:= diam * botm * isom * topm;

# transpose topm:
ll:= List(ccs, Size);
N:= Length(mat);
tap:= List([1..N], i-> List([1..N], j-> topm[i][j] * ll[j] / ll[i]));
tap:= TransposedMat(tap);

# ... and neutralize noncyclic part:
potm:= MutableCopyMat(tap);
for i in [1..N] do
    if not IsCyclic(AsGroup(Sections(trips[i])[1])) then
        potm[i]{[1..i-1]}:= 0*[1..i-1];
    fi;
od;

# Burnside
bas:= BasisDoubleBurnsideRing(G);

chg:= diam * botm * isom;
a:= RightRegularBaseChange(bas.basis, chg);;

chg:= diam * botm * isom * potm;
b:= RightRegularBaseChange(bas.basis, chg);;

##  this is the number of G-conjugates of (P_1,K_1) -> U:
##  the index of N_G(P_1, K_1) in G  x  the size of Aut_{\theta_1}(U)
wt1:= List(sec1s, x-> Index(G, NormalizerSection(x))
           * Size(Conjugators(OneMorphismSection(x))));

wt1m:= DiagonalMat(wt1/Size(G));

## this is |Aut(P_1)| / |Aut(P_1/K_1)|
wt2:= List(sec1s, x-> Size(AutomorphismGroup(TopSec(x)))
           / Size(AutomorphismGroup(AsGroup(x))));

wt2m:= DiagonalMat(wt2);


#chg:= diam * botm * isom * wgtm * potm;
chg:= diam * botm * isom * wt2m * potm * wt1m;
d:= RightRegularBaseChange(bas.basis, chg);;
new:= chg;

cols:= Concatenation(
               List([1..4], i-> i + 4 * [0..3]),
               List([17..18], i-> i + 2 * [0..1]),
               List([21..22], i-> [i])
               );

shrink:= chg^0;
for col in cols do
    i:= col[1];
    for j in col{[2..Length(col)]} do
        shrink[j][i]:= -1;
    od;
od;
shrink1:= shrink^-1;

firsts:= List(cols, x-> x[1]);
seconds:= Difference([1..Length(chg)], firsts);
poss:= Concatenation(firsts, seconds);

f:= List(d, x-> x^shrink1);;
m:= List(f, x-> x{firsts}{firsts});;
h:= List(f, x-> x{seconds}{firsts});;

invo:= 0 * chg;
ii:= Concatenation(cols);
for i in [1..Length(ii)] do
    invo[i][ii[i]]:= 1;
od;

# marks
mmm:= List(chg, x-> Sum([1..Length(x)], i-> x[i] * m[i]));;

# original marks
org0:= m-> [
  [m[ 1], m[ 2], m[ 3], m[ 4], 0, 0, 0, 0],
  [m[ 5], m[ 6], m[ 7], m[ 8], 0, 0, 0, 0],
  [m[ 9], m[10], m[11], m[12], 0, 0, 0, 0],
  [m[13], m[14], m[15], m[16], 0, 0, 0, 0],
  [0, 0, 0, 0, m[17], m[18], 0, 0],
  [0, 0, 0, 0, m[19], m[20], 0, 0],
  [0, 0, 0, 0, 0, 0, m[21], 0],
  [0, 0, 0, 0, 0, 0, 0, m[22]],
];

org:= function(m)
    local   d,  aaa,  iii,  i;
    d:= 8;
    aaa:= NullMat(d, d);
    iii:= [1..4];
    for i in iii do
        aaa[ i]{iii}:= m{(i-1)*4 + iii};
    od;
    iii:= [1..2];
    aaa[ 5]{4+iii}:= m{16+iii};
    aaa[ 6]{4+iii}:= m{18+iii};
    aaa[7][7]:= m[21];
    aaa[8][8]:= m[22];
    return aaa;
end;


mymarks:= function(m)
    local   d,  aaa,  iii,  i;
    d:= 8;
    aaa:= NullMat(d, d);
    iii:= [1..4];
    for i in [1..3] do
        aaa[ i]{iii}:= m{(i-1)*4 + iii};
    od;
    i:= 4;
    aaa[ 8]{iii}:= m{(i-1)*4 + iii};
    iii:= [1..2];
    aaa[ 5]{4+iii}:= m{16+iii};
    aaa[ 8]{4+iii}:= m{18+iii};
    aaa[7][7]:= m[21];
    aaa[4][4]:= m[22];
    aaa[6][6]:= m[22];
    aaa[8][8]:= m[22];
    return aaa;
end;


ddd:= DiagonalMat([1, 3, 2, 6, 3, 6, 2, 6]/6);
ppp:= ddd^0;
for ij in [
        [5,2], [7,3]
        ] do
    i:= ij[1];  j:= ij[2];
    ppp[i][j]:= 1;
od;