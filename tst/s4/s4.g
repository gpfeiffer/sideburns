##  the group
G:=  SymmetricGroup(4);
SetName(G, "S4");
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

###  this is the number of G-conjugates of (P_1,K_1) -> U:
###  the index of N_G(P_1, K_1) in G  x  the size of Aut_{\theta_1}(U)

wt1:= List(sec1s, x-> Index(G, NormalizerSection(x))
           * Size(Conjugators(OneMorphismSection(x))));

wt1m:= DiagonalMat(wt1);

### this is |Aut(P_1)| / |Aut(P_1/K_1)|
wt2:= List(sec1s, x-> Size(AutomorphismGroup(TopSec(x)))
           / Size(AutomorphismGroup(AsGroup(x))));

wt2m:= DiagonalMat(wt2/Size(G));

chg:= diam * botm * isom * wt2m * potm * wt1m;
c:= RightRegularBaseChange(bas.basis, chg);;

max:= chg^0;;
x:= [[-1,1,1],[1,-1,1],[1,1,-1]];
for p in List([199..209], i -> i + 11 * [0..2]) do
    max{p}{p}:= x;
od;

max1:= chg^0;;
x:= [[-1,1,1],[1,-1,1],[1,1,-1]];
for p in List([0..10], i -> 11*i + [129..131]) do
    max1{p}{p}:= -2*x;
od;

d:= RightRegularBaseChange(c, max1^-1);;

max1a:= chg^0;;
x:= [[-1,1,1],[1,-1,1],[1,1,-1]];
for p in List([0..10], i -> 11*i + [125..127]) do
    max1a{p}{p}:= x;
od;

dix1a:= chg^0;;
x:= [[-1,1,1],[1,-1,1],[1,1,-1]];
for p in Concatenation(List([0..10], i -> 11*i + [125..127])) do
    dix1a[p][p]:= -3;
od;

min:= chg^0;;
one:= [1..121];
two:= [122..242];
kl4:= [248..266];
dih:= [271..272];
min{kl4}{two}:= tap{kl4}{two};; # 2 in V4
min{dih}{kl4}:= tap{dih}{kl4};; # V4 in D8
min{dih}{two}:= tap{dih}{two};; # 2 in D8
min{dih}{one}:= tap{dih}{one};; # 1 in D8


e:= RightRegularBaseChange(d, min);;

new:= chg / max1 * min;

cols:= Concatenation(
               List([1..11], i-> i + 11 * [0..10]),
               List([122..132], i-> i + 11 * [0..10]),
               [[243,245],[244,246]],
               [[247]],
               [[248,251,257],
                [258,252,249],
                [259,253,249],
                [260,253,249],
                [261,254,250],
                [262,255,250],
                [263,255,250],
                [264,256,250],
                [265,256,250],
                [266,254,250],
                ],
               #### ??? ####
               [[267,269],[268,270]],
               [[271],[272]],
               [[273]],
               [[274]]
               );

shrink:= chg^0;;
for col in cols do
    i:= col[1];
    for j in col{[2..Length(col)]} do
        shrink[j][i]:= -1;
    od;
od;
shrink1:= shrink^-1;;

firsts:= List(cols, x-> x[1]);;
seconds:= Difference([1..Length(chg)], firsts);;
poss:= Concatenation(firsts, seconds);;

f:= List(d, x-> x^shrink1);;
m:= List(f, x-> x{firsts}{firsts});;
h:= List(f, x-> x{seconds}{firsts});;
