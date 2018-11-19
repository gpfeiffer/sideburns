##  the group
G:= Group((1,2),(3,4));
SetName(G, "V4");
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

##  I'm not sure why, but changing these from 3 to -1 does the trick!!
#wt1{[28,29]}:= -[1,1];

wt1m:= DiagonalMat(wt1);

# this is |Aut(P_1)| / |Aut(P_1/K_1)|
wt2:= List(sec1s, x-> Size(AutomorphismGroup(TopSec(x)))
           / Size(AutomorphismGroup(AsGroup(x))));

wt2m:= DiagonalMat(wt2/Size(G));

chg:= diam * botm * isom * wt2m * potm * wt1m;;
c:= RightRegularBaseChange(bas.basis, chg);;

m:= [[-1,1,1],[1,-1,1],[1,1,-1]];
max:= chg^0;
for p in List([0..5], k-> k + [44,50,56]) do
    max{p}{p}:= m;
od;
max1:= chg^0;;
for p in List([0..5], k-> 6*k + [29..31]) do
    max1{p}{p}:= m;
od;
max2:= chg^0;;
for p in List([0..5], k-> 6*k + [29..31]) do
    max2{p}{p}:= -3*m;
od;

# rows or cols? In fact, d = d1 !!
d:= RightRegularBaseChange(c, max^-1);;
d1:= RightRegularBaseChange(c, max1^-1);;

# some columns (with P = V4) need to be tripled
dia:= chg^0;;
for i in Concatenation(List([0..5], k-> 6*k + [29..31])) do
    dia[i][i]:= -3;
od;

# cols! dd = d2 !!
d2:= RightRegularBaseChange(d, dia^-1);;
d3:= RightRegularBaseChange(c, max2^-1);;

# reinstate potm between V4 and 2
min:= chg^0;;
min{[62..67]}{[26..61]}:= tap{[62..67]}{[26..61]};

e:= RightRegularBaseChange(d2, min);;

new:= chg / max / dia * min;
new2:= chg / max2 * min;

cols:= Concatenation(
               List([1..5], i-> i + 5 * [0..4]),
               List([26..31], i-> i + 6 * [0..5]),
               List([62..67], i-> [i])
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

f:= List(e, x-> x^shrink1);;
m:= List(f, x-> x{firsts}{firsts});;
h:= List(f, x-> x{seconds}{firsts});;

# marks
mmm:= List(new, x-> Sum([1..Length(x)], i-> x[i] * m[i]));;
mmm2:= List(new2, x-> Sum([1..Length(x)], i-> x[i] * m[i]));;


fous3:= [
 [ 1,  1,  1,  1,  1,  1 ]/6,
 [ 1, -1,  0,  1, -1,  0 ]/3,
 [ 0,  1, -1,  0,  1, -1 ]/3,
 [ 0,  0,  1,  1, -1, -1 ]/3,
 [ 1,  1, -1, -1,  0,  0 ]/3,
 [ 1, -1,  1, -1,  1, -1 ]/6,
];

fou:= chg^0;
for iii in [ [62..67] ] do
    fou{iii}{iii}:= fous3;
od;

ee:= RightRegularBaseChange(e, fou^-1);;

fin:= new / fou;;
fin2:= new2/ fou;;

ff:= List(ee, x-> x^shrink1);;
mm:= List(ff, x-> x{firsts}{firsts});;
hh:= List(ff, x-> x{seconds}{firsts});;

