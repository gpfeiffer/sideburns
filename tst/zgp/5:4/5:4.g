##  the group
C:= CyclicGroup(IsPermGroup, 5);
G:= Normalizer(SymmetricGroup(5), C);
SetName(G, "5:4");
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

wt1m:= DiagonalMat(wt1/Size(G));

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

cols:= Concatenation(
               List([1..6], i-> i + 6 * [0..5]), #1
               List([37..40], i-> i + 4 * [0..3]), #2
               List([53..56], i-> i + 4 * [0..1]),
               List([61..63], i-> [i])
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

z:= d;
f:= List(z, x-> x^shrink1);;
m:= List(f, x-> x{firsts}{firsts});;
h:= List(f, x-> x{seconds}{firsts});;

fou:= chg^0;
for p in [ [53,54], [55,56], [57,58], [59,60] ] do
    fou{p}{p}:= [[1,1],[1,-1]];
od;

e:= RightRegularBaseChange(d, fou);;

