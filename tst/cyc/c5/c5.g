##  the group
n:= 5;
G:=  CyclicGroup(IsPermGroup, n);
SetName(G, Concatenation("C", String(n)));
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
tops:= List(ccs, x-> PositionProperty(top.reps, r-> r in x));
bots:= List(ccs, x-> PositionProperty(bot.reps, r-> r in x));
isos:= List(ccs, x-> PositionProperty(iso.reps, r-> r in x));
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
potm:= TransposedMat(tap);

# Burnside
bas:= BasisDoubleBurnsideRing(G);

chg:= diam * botm * isom;
a:= RightRegularBaseChange(bas.basis, chg);;

chg:= diam * botm * isom * potm;
b:= RightRegularBaseChange(bas.basis, chg);;

# this is |Aut(P_1)| / |Aut(P_1/K_1)|
wt2:= List(sec1s, x-> Size(AutomorphismGroup(TopSec(x)))
           / Size(AutomorphismGroup(AsGroup(x))));

wt2m:= DiagonalMat(wt2/Size(G));

chg:= diam * botm * isom * wt2m * potm;;
c:= RightRegularBaseChange(bas.basis, chg);;

d:= c;

cols:= Concatenation(
               List([1..2], i-> i + 2 * [0..1]),
               List([5..8], i-> [i])
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

fou:= chg^0;
for p in [ [5,6] ] do
    fou{p}{p}:= [[1,1],[1,-1]];
od;

e:= RightRegularBaseChange(d, fou);;
