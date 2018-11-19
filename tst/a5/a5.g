##  the group
G:=  AlternatingGroup(5);
SetName(G, "A5");
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

#  this is the number of G-conjugates of (P_1,K_1) -> U:
#  the index of N_G(P_1, K_1) in G  x  the size of Aut_{\theta_1}(U)
wt1:= List(sec1s, x-> Index(G, NormalizerSection(x))
           * Size(Conjugators(OneMorphismSection(x))));

wt1m:= DiagonalMat(wt1/Size(G));;

# this is |Aut(P_1)| / |Aut(P_1/K_1)|
wt2:= List(sec1s, x-> Size(AutomorphismGroup(TopSec(x)))
           / Size(AutomorphismGroup(AsGroup(x))));

wt2m:= DiagonalMat(wt2);

chg:= diam * botm * isom * wt2m * potm * wt1m;
c:= RightRegularBaseChange(bas.basis, chg);;

# some columns (with P = V4) need to be tripled
dia:= chg^0;;
for i in [83,87,91,95] do
    dia[i][i]:= -3;
od;

# reinstate potm between V4 and 2
min:= chg^0;;
min{[103,104]}{[82..97]}:= tap{[103,104]}{[82..97]};

d:= RightRegularBaseChange(c, dia^-1 * min);;

new:= chg / dia * min;

cols:= [
        1 + 9 * [0..8],
        2 + 9 * [0..8],
        3 + 9 * [0..8],
        4 + 9 * [0..8],
        5 + 9 * [0..8],
        6 + 9 * [0..8],
        7 + 9 * [0..8],
        8 + 9 * [0..8],
        9 + 9 * [0..8],
        82 + 4 * [0..3],
        83 + 4 * [0..3],
        84 + 4 * [0..3],
        85 + 4 * [0..3],
        [98, 100],
        [101,99],
        [102,99],
        [103],
        [104],
        [105],
        [106],
        [107],
        [108],
        [109],
        [110],
        [111],
        [112],
        [113]
        ];

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

fou:= chg^0;;
for p in [ [101,102], [103,104], [105,106], [108,109], [110,111], [112,113] ] do
    fou{p}{p}:= [[1,1],[1,-1]]/2;
od;
fou[99][99]:= 1/2;

e:= RightRegularBaseChange(d, fou^-1);;

cols[16]:= [102];

sshrink:= chg^0;
for col in cols do
    i:= col[1];
    for j in col{[2..Length(col)]} do
        sshrink[j][i]:= -1;
    od;
od;
sshrink1:= sshrink^-1;

ff:= List(e, x-> x^sshrink1);;
mm:= List(ff, x-> x{firsts}{firsts});;
hh:= List(ff, x-> x{seconds}{firsts});;
