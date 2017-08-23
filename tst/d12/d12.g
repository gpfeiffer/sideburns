##  the group
G:=  DihedralGroup(IsPermGroup, 12);
SetName(G, "D12");
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

now:= chg^0;;
now[249]{[ 140, 153, 166 ]}:= [1,1,1];
now[250]{[ 140, 154, 165 ]}:= [1,1,1];
now[251]{[ 141, 154, 164 ]}:= [1,1,1];
now[252]{[ 141, 152, 166 ]}:= [1,1,1];
now[253]{[ 142, 152, 165 ]}:= [1,1,1];
now[254]{[ 142, 153, 164 ]}:= [1,1,1];
now[255]{[ 148, 158, 171 ]}:= [1,1,1];
now[256]{[ 148, 159, 170 ]}:= [1,1,1];
now[257]{[ 146, 159, 172 ]}:= [1,1,1];
now[258]{[ 146, 160, 171 ]}:= [1,1,1];
now[259]{[ 147, 160, 170 ]}:= [1,1,1];
now[260]{[ 147, 158, 172 ]}:= [1,1,1];

now[261]{[ 213, 226, 236 ]}:= [1,1,1];
now[262]{[ 214, 225, 236 ]}:= [1,1,1];
now[263]{[ 214, 224, 237 ]}:= [1,1,1];
now[264]{[ 212, 226, 237 ]}:= [1,1,1];
now[265]{[ 212, 225, 238 ]}:= [1,1,1];
now[266]{[ 213, 224, 238 ]}:= [1,1,1];
now[267]{[ 218, 231, 244 ]}:= [1,1,1];
now[268]{[ 219, 230, 244 ]}:= [1,1,1];
now[269]{[ 219, 232, 242 ]}:= [1,1,1];
now[270]{[ 220, 231, 242 ]}:= [1,1,1];
now[271]{[ 220, 230, 243 ]}:= [1,1,1];
now[272]{[ 218, 232, 243 ]}:= [1,1,1];

now[273][ 192 ]:= 1;
now[274][ 191 ]:= 1;
now[275][ 195 ]:= 1;
now[276][ 180 ]:= 1;
now[277][ 179 ]:= 1;
now[278][ 183 ]:= 1;
now[279][ 228 ]:= 1;
now[280][ 227 ]:= 1;
now[281][ 231 ]:= 1;


#now[279]{[ 97, 228 ]}:= [6,6];
#now[280]{[ 99, 227 ]}:= [6,6];
#now[281]{[ 100, 231 ]}:= [6,6];

#now[283]{[218,244]}:= -2* [1,1];
#now[284]{[220,242]}:= -2* [1,1];


chg:= diam * botm * isom * wt2m * potm * wt1m;;
c:= RightRegularBaseChange(bas.basis, chg);;

# unwind
max:= chg^0;;
m:= [[-1,1,1],[1,-1,1],[1,1,-1]];
for p in List([1..12], i-> 136 + i + 12 * [0..2]) do
    max{p}{p}:= m;
od;

dix:= chg^0;;
for p in Concatenation(List([0..11], i-> 12*i + [104..106])) do
    dix[p][p]:= -3;
od;
for p in Concatenation(List([0..11], i-> 12*i + [110..112])) do
    dix[p][p]:= -3;
od;

dix1:= chg^0;;
for p in Concatenation(List([0..11], i-> 12*i + [104..106])) do
    dix1[p][p]:= -3;
od;

dix2:= chg^0;
for p in Concatenation(List([0..11], i-> 12*i + [110..112])) do
    dix2[p][p]:= -3;
od;

max1:= chg^0;;
m:= [[-1,1,1],[1,-1,1],[1,1,-1]];
for p in List([0..11], i-> 12*i + [104..106]) do
    max1{p}{p}:= -3*m;
od;

#di1:= List([1..Length(chg)], x-> 1);;
#for iii in List([0..11], i-> 12*i + [110..112]) do
#    di1{iii}:= -2*[1,1,1];
#od;
#dia1:= DiagonalMat(di1);;

#dia1:= chg^0;;
#m:= [[-1,1,1],[1,-1,1],[1,1,-1]];
#for iii in List([0..11], i-> 12*i + [110..112]) do
#    dia1{iii}{iii}:= -2*m;
#od;

#dia1:= chg^0;;
#m:= [[-1,1,1],[1,-1,1],[1,1,-1]];
#for iii in List([0..11], i-> i + [209,221,233]) do
#    dia1{iii}{iii}:= -m;
#od;

di1:= List([1..Length(chg)], x-> 1);;
for i in Concatenation([255..260],[267..272]) do
    di1[i]:= 1/3;
od;
dia1:= DiagonalMat(di1);;
d0:= RightRegularBaseChange(c, dia1^-1);;


di2:= List([1..Length(chg)], x-> 1);;
di2{[275,278,281]}:= -2*[1,1,1];
dia2:= DiagonalMat(di2);;

d1:= RightRegularBaseChange(c, max1^-1);;

d:= RightRegularBaseChange(d1, dia2^-1);;

## reinstate potm between V4 and 2
min:= chg^0;;
min{[249..272]}{[101..244]}:= tap{[249..272]}{[101..244]};;

## reinstate some potm between V4 and 2
minx:= chg^0;;
minx{[249..254]}{[101..244]}:= tap{[249..254]}{[101..244]};;
minx{[261..266]}{[101..244]}:= tap{[261..266]}{[101..244]};;

e:= RightRegularBaseChange(d, min);;

#min{[249..254]}[56]:= [1,1,1,1,1,1];
#min{[255..260]}[60]:= [1,1,1,1,1,1];
#min{[261..266]}[96]:= [1,1,1,1,1,1];
#min{[267..272]}[100]:= [1,1,1,1,1,1];
#min[273]{[67, 192]}:= [1,1];
#min[274]{[69, 191]}:= [1,1];
#min[275]{[70, 195]}:= [1,1];
#min[276]{[87, 180]}:= [1,1];
#min[277]{[89, 179]}:= [1,1];
#min[278]{[90, 183]}:= [1,1];
#min[279]{[97, 228]}:= [1,1];
#min[280]{[99, 227]}:= [1,1];
#min[281]{[100, 231]}:= [1,1];
#min[283]{[100, 218, 231, 244, 267, 281]}:= [1,1,1,1,1,1];
#min[284]{[100, 220, 231, 242, 270, 281]}:= [1,1,1,1,1,1];
#
#e:= RightRegularBaseChange(d, min^-1);;

new:= chg/max/min;;

fous3:= [
 [ 1,  1,  1,  1,  1,  1 ]/6,
 [ 1, -1,  0,  1, -1,  0 ]/3,
 [ 0,  1, -1,  0,  1, -1 ]/3,
 [ 0,  0,  1,  1, -1, -1 ]/3,
 [ 1,  1, -1, -1,  0,  0 ]/3,
 [ 1, -1,  1, -1,  1, -1 ]/6,
];

#r:= ER(3);
#fous3:= [
# [ 1,  1,  1,  1,  1,  1 ],
# [ 2, 2, -1, -1, -1, -1 ]/2,
# [ 0, 0, r, r, -r, -r ]/2,
# [ 0, 0, -r, r, r, -r ]/2,
# [ 2, -2, -1, 1, -1, 1 ]/2,
# [ 1, -1,  1, -1,  1, -1 ],
#         ];
#
#fous3:= [
# [ 1,  1,  1,  1,  1,  1 ],
# [ 1, -1, -1/2, 1/2, -1/2, 1/2 ],
# [ 0, 0, -3/4, 3/4, 3/4, -3/4 ],
# [ 0, 0, 1, 1, -1, -1 ],
# [ 1, 1, -1/2, -1/2, -1/2, -1/2 ],
# [ 1, -1,  1, -1,  1, -1 ],
#];


fou2:= [
  [ 1,  1 ]/2,
  [ 1, -1 ]/2,
];

fou:= chg^0;;
iii:= [283,284];
fou{iii}{iii}:= fou2;
for iii in [ [249..254], [255..260], [261..266], [267..272] ] do
    fou{iii}{iii}:= fous3;
od;

foux:= chg^0;;
for iii in List([1..12], i-> i + [208, 232]) do
    foux{iii}{iii}:= fou2;
od;

## c9:= RightRegularBaseChange(c1, wt2m * potm * wt1m / max1 / dia2 * min / fou / foux);;
## 2891, but h <> 0  :-(

# columns!
fou1:= chg^0;;
for iii in List([1..12], i-> 12*i + [98,100]) do
    fou1{iii}{iii}:= fou2;
od;

##c9:= RightRegularBaseChange(c1, wt2m * potm * wt1m / max / dix / dia2 * min / fou / fou1 );;
## 2761

max11:= chg^0;;
m:= [[-1,1,1],[1,-1,1],[1,1,-1]];
for iii in List([1..12], i-> 12*i + [98..100]) do
    max11{iii}{iii}:= m;
od;

max12:= chg^0;;
m:= [[-1,1,1],[1,-1,1],[1,1,-1]];
for iii in List([0..11], i-> i + [209,221,233]) do
    max12{iii}{iii}:= m;
od;


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
               List([1..10], i-> i + 10 * [0..9]),
               List([101..112], i-> i + 12 * [0..11]),
               List([245..246], i-> i + 2 * [0..1]),
               List([249..260], i-> i + 12 * [0..1]),
               List([273..275], i-> i + 3 * [0..2]),
               List([282..284], i-> [i])
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
