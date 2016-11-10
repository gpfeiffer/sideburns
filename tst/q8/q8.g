##  the group
G:=  Group((1,2,7,5)(3,6,4,8), (1,4,7,3)(2,6,5,8));  # stab of 11, 10, 9 in M11
SetName(G, "Q8");
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

##  I'm not sure why, but changing these ...
#wt1{[31..36]}:= 1/2*[1,1,1,1,1,1];
#wt1{[65..85]}:= 0*[65..85] + 1/6;
#wt1{[101..106]}:= 1/2*[1,1,1,1,1,1];

#wt1{[65..85]}:= 0*[65..85] - 1/6;
#wt1{[95..100]}:= 1/2*[1,1,1,1,1,1];

wt1m:= DiagonalMat(wt1/Size(G));

# this is |Aut(P_1)| / |Aut(P_1/K_1)|
wt2:= List(sec1s, x-> Size(AutomorphismGroup(TopSec(x)))
           / Size(AutomorphismGroup(AsGroup(x))));

wt2m:= DiagonalMat(wt2);

chg:= diam * botm * isom * wt2m * potm * wt1m;;
c:= RightRegularBaseChange(bas.basis, chg);;

# unwind
max:= chg^0;;
m:= [[-1,1,1],[1,-1,1],[1,1,-1]];
for p in List([65..71], i -> i + 7 * [0..2]) do
    max{p}{p}:= -6*m;
od;

d:= RightRegularBaseChange(c, max^-1);;

di2:= List([1..Length(c)], x-> 1);
di2{[95..100]}:= 2*[1,1,1,1,1,1];
dia2:= DiagonalMat(di2);;
d1:= RightRegularBaseChange(d, dia2^-1);;

# reinstate potm between V4 and 2, and between Q8 and V4, Q8 and 2.
min:= chg^0;;
min{[95..100]}{[37..85]}:= tap{[95..100]}{[37..85]};;
min{[101..106]}{[95..100]}:= tap{[101..106]}{[95..100]}/4;;
min{[101..106]}{[37..85]}:= tap{[101..106]}{[37..85]}/4;;


e:= RightRegularBaseChange(d1, min);;

fous3:= [
 [ 1,  1,  1,  1,  1,  1 ]/6,
 [ 1, -1,  0,  1, -1,  0 ]/3,
 [ 0,  1, -1,  0,  1, -1 ]/3,
 [ 0,  0,  1,  1, -1, -1 ]/3,
 [ 1,  1, -1, -1,  0,  0 ]/3,
 [ 1, -1,  1, -1,  1, -1 ]/6,
];

fou:= chg^0;
for iii in [ [95..100], [101..106] ] do
    fou{iii}{iii}:= fous3;
od;

ee:= RightRegularBaseChange(e, fou^-1);;

cols:= Concatenation(
               List([1..6], i-> i + 6 * [0..5]),
               List([37..43], i-> i + 7 * [0..6]),
               List([86..88], i-> i + 3 * [0..2]),
               List([95..106], i-> [i])
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

ff:= List(ee, x-> x^shrink1);;
mm:= List(ff, x-> x{firsts}{firsts});;
hh:= List(ff, x-> x{seconds}{firsts});;

