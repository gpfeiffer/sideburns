##  the group
G:=  DihedralGroup(IsPermGroup, 42);
SetName(G, "D42");
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

###  this is the number of G-conjugates of (P_1,K_1) -> U:
###  the index of N_G(P_1, K_1) in G  x  the size of Aut_{\theta_1}(U)
wt1:= List(sec1s, x-> Index(G, NormalizerSection(x))
           * Size(Conjugators(OneMorphismSection(x))));

wt1m:= DiagonalMat(wt1/Size(G));

### this is |Aut(P_1)| / |Aut(P_1/K_1)|
wt2:= List(sec1s, x-> Size(AutomorphismGroup(TopSec(x)))
           / Size(AutomorphismGroup(AsGroup(x))));

wt2m:= DiagonalMat(wt2);

chg:= diam * botm * isom * wt2m * potm * wt1m;
c:= RightRegularBaseChange(bas.basis, chg);;
#
## this has to do with the radical
min:= chg^0;
min[85]{[28,70]}:= [1,1];
min[86]{[32,72]}:= [1,1];
min[87]{[60,78]}:= [1,1];
min[88]{[64,80]}:= [1,1];

min[101]{[46,75]}:= [1,1];
min[102]{[46,75]}:= [1,1];
min[103]{[46,75]}:= [1,1];
min[104]{[48,76]}:= [1,1];
min[105]{[48,76]}:= [1,1];
min[106]{[48,76]}:= [1,1];
min[107]{[62,79]}:= [1,1];
min[108]{[62,79]}:= [1,1];
min[109]{[62,79]}:= [1,1];
min[110]{[64,80]}:= [1,1];
min[111]{[64,80]}:= [1,1];
min[112]{[64,80]}:= [1,1];

min[119]{[64,80,88,110]}:= [1,1,1,1];
min[120]{[64,80,88,111]}:= [1,1,1,1];
min[121]{[64,80,88,112]}:= [1,1,1,1];
min[122]{[64,80,88,112]}:= [1,1,1,1];
min[123]{[64,80,88,111]}:= [1,1,1,1];
min[124]{[64,80,88,110]}:= [1,1,1,1];

d:= RightRegularBaseChange(c, min^-1);;

new:= chg / min;

cols:= Concatenation(
               List([1..8], i-> i + 8 * [0..7]),
               List([65..68], i-> i + 4 * [0..3]),
               List([81..82], i-> i + 2 * [0..1]),
               List([85..86], i-> i + 2 * [0..1]),
               List([89..94], i-> i + 6 * [0..1]),
               List([101..106], i-> i + 6 * [0..1]),
               List([113..124], i-> [i])
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

