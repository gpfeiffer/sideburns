##  the group
G:=  DihedralGroup(IsPermGroup, 14);
SetName(G, "D14");
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

chg:= diam * botm * isom * potm;;
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

chg:= diam * botm * isom * wt2m * potm * wt1m;;
c:= RightRegularBaseChange(bas.basis, chg);;


e:= c;

new:= chg;

fous3:= [
 [ 1,  1,  1,  1,  1,  1 ]/6,
 [ 1, -1,  0,  1, -1,  0 ]/3,
 [ 0,  1, -1,  0,  1, -1 ]/3,
 [ 0,  0,  1,  1, -1, -1 ]/3,
 [ 1,  1, -1, -1,  0,  0 ]/3,
 [ 1, -1,  1, -1,  1, -1 ]/6,
];


fou2:= [
  [ 1,  1 ]/2,
  [ 1, -1 ]/2,
];


#dia2:= List([1..Length(chg)], x-> 1);
#for i in [249, 255, 261, 267, 254, 260, 266, 272] do
#    dia2[i]:= 1/6;
#od;
#for i in [250,251,252,253,256,257,258,259,262,263,264,265,268,269,270,271] do
#    dia2[i]:= 1/3;
#od;
#dia2[283]:= 1/2;
#dia2[284]:= 1/2;
#di2:= DiagonalMat(dia2);;
#
#e:= RightRegularBaseChange(c, (di2 * fou * min * max)^-1);;
#
#fin:= chg^0;
#for iii in List([1..12], i-> 12*i + [98,100]) do
#    fin{iii}{iii}:= fou2;
#od;
#e:= RightRegularBaseChange(e, fin);;
#
#
#dia:= chg^0;
#for i in [6,12,18,24,30,36,39,42,45,51,53] do
#    dia[i][i]:= 3;
#od;
#
#e:= RightRegularBaseChange(d, dia^-1);;
#
#fin:= chg^0;
#fin{[57,58,59]}[53]:= -[1,1,1];
#
#e:= RightRegularBaseChange(e, fin^-1);;
#
#new:= chg / min / dia / fin;

e:= c;

cols:= Concatenation(
               List([1..4], i-> i + 4 * [0..3]), #1
               List([17..18], i-> i + 2 * [0..1]), #2
               List([21..26], i-> [i])
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

f:= List(e, x-> x^shrink1);;
m:= List(f, x-> x{firsts}{firsts});;
h:= List(f, x-> x{seconds}{firsts});;
