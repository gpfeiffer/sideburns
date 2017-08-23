##  the group
G:=  DihedralGroup(IsPermGroup, 8);
SetName(G, "D8");
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

#for i in [200..211] do
#    wt1[i]:= 2;
#od;


wt1m:= DiagonalMat(wt1);;

### this is |Aut(P_1)| / |Aut(P_1/K_1)|
wt2:= List(sec1s, x-> Size(AutomorphismGroup(TopSec(x)))
           / Size(AutomorphismGroup(AsGroup(x))));

wt2m:= DiagonalMat(wt2/Size(G));

### this seems more useful here.
wt3:= List(sec1s, x-> 1);
for i in Concatenation([41..48], [109..119]) do
    wt3[i]:= 2;
od;
wt3m:= DiagonalMat(wt3/Size(G));


chg:= diam * botm * isom * wt2m * potm * wt1m;;
c:= RightRegularBaseChange(bas.basis, chg);;

m:= [[-1,1,1],[1,-1,1],[1,1,-1]];
max:= chg^0;;
for p in List([0..10], k-> 11*k + [73..75]) do
    max{p}{p}:= -2*m;
od;
d:= RightRegularBaseChange(c, max^-1);;

##  leaves 4644 non-0's

di1:= List(chg, x-> 1);
for p in [71,82..181] do
    di1[p]:= -3;
od;
for p in [72,83..182] do
    di1[p]:= -3;
od;
for p in [188,195,201,202,203,190,197,204,205,206] do
    di1[p]:= 1/2;
od;
for p in [191,192,193,198,199,200,207,208,209,210,211,212] do
    di1[p]:= 2/3;
od;
dia1:= DiagonalMat(di1);;
d1:= RightRegularBaseChange(d, dia1^-1);;

min:= chg^0;;
# this is too much too early ...
min{[187..212]}{[65..185]}:= tap{[187..212]}{[65..185]};;
min{[213,214]}{[187..212]}:= tap{[213,214]}{[187..212]};;
min{[213,214]}{[65..185]}:= tap{[213,214]}{[65..185]};;
#min{[213,214]}{[1..64]}:= tap{[213,214]}{[1..64]};
e:= RightRegularBaseChange(d1, min);;


## this has to do with the radical
#min:= chg^0;
#min[93]{[37,70]}:= [1,1];
#min[94]{[40,72]}:= [1,1];
#min[95]{[61,78]}:= [1,1];
#min[96]{[64,80]}:= [1,1];
#
#min[97]{[46,75]}:= [1,1];
#min[98]{[46,75]}:= [1,1];
#min[99]{[48,76]}:= [1,1];
#min[100]{[48,76]}:= [1,1];
#min[101]{[62,79]}:= [1,1];
#min[102]{[62,79]}:= [1,1];
#min[103]{[64,80]}:= [1,1];
#min[104]{[64,80]}:= [1,1];
#
#min[109]{[64,80,96,103]}:= [1,1,1,1];
#min[110]{[64,80,96,104]}:= [1,1,1,1];
#min[111]{[64,80,96,104]}:= [1,1,1,1];
#min[112]{[64,80,96,103]}:= [1,1,1,1];
#

#min[186]{[37,101,137]}:= [1,1,1];
#min[187]{[37,101,137]}:= 2*[1,1,1];
#min[188]{[39,103,138]}:= [1,1,1];
#min[189]{[39,103,138]}:= 2*[1,1,1];
#min[190]{[40,139,140,141]}:= [1,1,1,1];
#min[191]{[40,139,140,141]}:= [1,1,1,1];
#min[192]{[40,139,140,141]}:= [1,1,1,1];

new:= chg / min;

poss2:=
[ [  186,  187,  187,  188,  189,  189,  190,  190,  191,  192,  192,  191 ],
  [  187,  186,  187,  189,  188,  189,  191,  192,  190,  191,  190,  192 ],
  [  187,  187,  186,  189,  189,  188,  192,  191,  192,  190,  191,  190 ],
  [  193,  194,  194,  195,  196,  196,  197,  197,  198,  199,  199,  198 ],
  [  194,  193,  194,  196,  195,  196,  198,  199,  197,  198,  197,  199 ],
  [  194,  194,  193,  196,  196,  195,  199,  198,  199,  197,  198,  197 ],
  [  200,  201,  202,  203,  204,  205,  206,  207,  209,  211,  210,  208 ],
  [  200,  202,  201,  203,  205,  204,  207,  206,  208,  210,  211,  209 ],
  [  201,  200,  202,  204,  203,  205,  209,  210,  206,  208,  207,  211 ],
  [  202,  201,  200,  205,  204,  203,  211,  208,  210,  206,  209,  207 ],
  [  202,  200,  201,  205,  203,  204,  208,  211,  207,  209,  206,  210 ],
  [  201,  202,  200,  204,  205,  203,  210,  209,  211,  207,  208,  206 ] ];

cols:= Concatenation(
               List([1..8], i-> i + 8 * [0..7]),
               List([65..75], i-> i + 11 * [0..10]),
               [[186]],
               [
                [201,187,194],
                [202,188,195],
                [203,188,195],
                [204,189,196],
                [205,190,197],
                [206,190,197],
                [207,191,198],
                [208,191,198],
                [209,192,199],
                [210,192,199],
                [211,193,200],
                [212,193,200]
                ]  ,
               List([213..214], i-> [i])
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
