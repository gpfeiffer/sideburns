##  the group
G:=  AlternatingGroup(4);
SetName(G, "A4");
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

chg:= diam * botm * isom * potm;;
b:= RightRegularBaseChange(bas.basis, chg);;

#  this is the number of G-conjugates of (P_1,K_1) -> U:
#  the index of N_G(P_1, K_1) in G  x  the size of Aut_{\theta_1}(U)
wt1:= List(sec1s, x-> Index(G, NormalizerSection(x))
           * Size(Conjugators(OneMorphismSection(x))));

wt1m:= DiagonalMat(wt1/Size(G));

# this is |Aut(P_1)| / |Aut(P_1/K_1)|
wt2:= List(sec1s, x-> Size(AutomorphismGroup(TopSec(x)))
           / Size(AutomorphismGroup(AsGroup(x))));

wt2m:= DiagonalMat(wt2);

chg:= diam * botm * isom * wt2m * potm * wt1m;;
c:= RightRegularBaseChange(bas.basis, chg);;

# some columns (with P = V4) need to be tripled
dia:= chg^0;;
for i in [27, 29] do
    dia[i][i]:= -3;
od;

# reinstate potm between V4 and 2
min:= chg^0;;
min{[38,39]}{[26..29]}:= tap{[38,39]}{[26..29]};

d1:= RightRegularBaseChange(c, dia^-1);;
d:= RightRegularBaseChange(d1, min);;

new:= chg / dia * min;

cols:= [
        1 + 5 * [0..4],
        2 + 5 * [0..4],
        3 + 5 * [0..4],
        4 + 5 * [0..4],
        5 + 5 * [0..4],
        [26,28],
        [27,29],
        [30,34],
        [31,35],
        [32,36],
        [33,37],
        [38],
        [39],
        [40],
        [41],
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

fou:= chg^0;
for p in [ [30,31], [32,33], [34,35], [36,37], [38,39], [40,41] ] do
    fou{p}{p}:= [[1,1],[1,-1]];
od;

e:= RightRegularBaseChange(d, fou);;

# marks
mmm:= List(new, x-> Sum([1..Length(x)], i-> x[i] * m[i]));;

# original marks
org:= function(m)
    local d, aaa;
    d:= 15;
    aaa:= NullMat(d, d);
    aaa[ 1]{[1..5]}:= m{[1..5]};
    aaa[ 2]{[1..5]}:= m{5+[1..5]};
    aaa[ 3]{[1..5]}:= m{10+[1..5]};
    aaa[ 4]{[1..5]}:= m{15+[1..5]};
    aaa[ 5]{[1..5]}:= m{20+[1..5]};
    aaa[ 6]{5+[1..2]}:= m{25+[1..2]};
    aaa[ 7]{5+[1..2]}:= m{27+[1..2]};
    aaa[ 8]{7+[1..4]}:= m{29+[1..4]};
    aaa[ 9]{7+[1..4]}:= m{[31,30,33,32]};
    aaa[10]{7+[1..4]}:= m{33+[1..4]};
    aaa[11]{7+[1..4]}:= m{[35,34,37,36]};
    aaa[12]{11+[1..2]}:= m{37+[1..2]};
    aaa[13]{11+[1..2]}:= m{[39,38]};
    aaa[14]{13+[1..2]}:= m{39+[1..2]};
    aaa[15]{13+[1..2]}:= m{[41,40]};
    return aaa;
end;

# try and make a mark hom
mymarks0:= function(m)
    local d, aaa;
    d:= 15;
    aaa:= NullMat(d, d);
    aaa[ 1]{[1..5]}:= 1/12 * m{[1..5]};
    aaa[ 2]{[1..5]}:= 1/4 * m{5+[1..5]};
    aaa[ 2][ 2]:= aaa[ 2][ 2] + 1/4 * m[26];
    aaa[ 2][ 4]:= aaa[ 2][ 4] + 1/4 * m[27];
    aaa[ 3]{[1..5]}:= 2/3 * m{10+[1..5]};
    aaa[ 3][ 3]:= aaa[ 3][ 3] + 1/3 * (m[30] + m[31]);
    aaa[ 3][ 5]:= aaa[ 3][ 5] + 1/3 * (m[32] + m[33]);
    aaa[ 4][ 4]:= 1/4 * (m[38] + m[39]);
    aaa[ 5][ 5]:= m[40] + m[41];

    aaa[ 6]{5+[1..2]}:= 1/4 * m{25+[1..2]};
    aaa[ 7]{5+[1..2]}:= 2/4 * m{27+[1..2]};
    aaa[ 7][ 7]:= aaa[ 7][ 7] + 1/4 * (m[38] + m[39]);

    aaa[ 8]{7+[1..4]}:= 1/3 * m{29+[1..4]};
    aaa[ 9]{7+[1..4]}:= 1/3 * m{[31,30,33,32]};
    aaa[10][10]:= m[40];
    aaa[10][11]:= m[41];
    aaa[11][10]:= m[41];
    aaa[11][11]:= m[40];

    aaa[12]{[1..5]}:= 2/4 * m{15+[1..5]};
    aaa[13]{[1..5]}:= 2/4 * m{15+[1..5]};
    aaa[12][ 2]:= aaa[12][ 2] + 2/4 * m[28];
    aaa[13][ 2]:= aaa[13][ 2] + 2/4 * m[28];
    aaa[12][ 4]:= aaa[12][ 4] + 2/4 * m[29];
    aaa[13][ 4]:= aaa[13][ 4] + 2/4 * m[29];
    aaa[12]{11+[1..2]}:= 1/4 * m{37+[1..2]};
    aaa[13]{11+[1..2]}:= 1/4 * m{[39,38]};

    aaa[14]{[1..5]}:= 2/1 * m{20+[1..5]};
    aaa[15]{[1..5]}:= 2/1 * m{20+[1..5]};
    aaa[14][ 3]:= aaa[14][ 3] + m[34] + m[35];
    aaa[14][ 5]:= aaa[14][ 5] + m[36] + m[37];
    aaa[15][ 3]:= aaa[15][ 3] + m[34] + m[35];
    aaa[15][ 5]:= aaa[15][ 5] + m[36] + m[37];
    aaa[14]{7+[1..4]}:= m{33+[1..4]};
    aaa[15]{7+[1..4]}:= m{[35,34,37,36]};
    aaa[14]{13+[1..2]}:= m{39+[1..2]};
    aaa[15]{13+[1..2]}:= m{[41,40]};
    return aaa;
end;

mymarks1:= function(m)
    local d, aaa;
    d:= 15;
    aaa:= NullMat(d, d);
    aaa[ 1]{[1..5]}:= 1/12 * m{[1..5]};
    aaa[ 2]{[1..5]}:= 1/4 * m{5+[1..5]};
    aaa[ 3]{[1..5]}:= 1/3 * m{10+[1..5]};
    aaa[ 4][ 4]:= 1/4 * (m[38] + m[39]);
    aaa[ 5][ 5]:= m[40] + m[41];

    aaa[ 6]{5+[1..2]}:= 1/4 * m{25+[1..2]};
    aaa[ 7]{5+[1..2]}:= 1/4 * m{27+[1..2]};

    aaa[ 8]{7+[1..4]}:= 1/3 * m{29+[1..4]};
    aaa[ 9]{7+[1..4]}:= 1/3 * m{[31,30,33,32]};
    aaa[10][10]:= m[40];
    aaa[10][11]:= m[41];
    aaa[11][10]:= m[41];
    aaa[11][11]:= m[40];

    aaa[12]{[1..5]}:= 1/4 * m{15+[1..5]};
    aaa[13]{[1..5]}:= 1/4 * m{15+[1..5]};
    aaa[12]{11+[1..2]}:= 1/4 * m{37+[1..2]};
    aaa[13]{11+[1..2]}:= 1/4 * m{[39,38]};

    aaa[14]{[1..5]}:= 1/1 * m{20+[1..5]};
    aaa[15]{[1..5]}:= 1/1 * m{20+[1..5]};
    aaa[14]{7+[1..4]}:= m{33+[1..4]};
    aaa[15]{7+[1..4]}:= m{[35,34,37,36]};
    aaa[14]{13+[1..2]}:= m{39+[1..2]};
    aaa[15]{13+[1..2]}:= m{[41,40]};
    return aaa;
end;

mymarks:= function(m)
    local d, aaa;
    d:= 15;
    aaa:= NullMat(d, d);
    aaa[ 1]{[1..5]}:= 1/12 * m{[1..5]};
    aaa[ 2]{[1..5]}:= 1/4 * m{5+[1..5]};
    aaa[ 3]{[1..5]}:= 2/3 * m{10+[1..5]};
    aaa[ 4][ 4]:= 1/4 * (m[38] + m[39]);
    aaa[ 5][ 5]:= m[40] + m[41];

    aaa[ 6]{5+[1..2]}:= 1/4 * m{25+[1..2]};
    aaa[ 7]{5+[1..2]}:= 2/4 * m{27+[1..2]};

    aaa[ 8]{7+[1..4]}:= 1/3 * m{29+[1..4]};
    aaa[ 9]{7+[1..4]}:= 1/3 * m{[31,30,33,32]};
    aaa[10][10]:= m[40];
    aaa[10][11]:= m[41];
    aaa[11][10]:= m[41];
    aaa[11][11]:= m[40];

    aaa[12]{[1..5]}:= 2/4 * m{15+[1..5]};
    aaa[13]{[1..5]}:= 2/4 * m{15+[1..5]};
    aaa[12]{11+[1..2]}:= 1/4 * m{37+[1..2]};
    aaa[13]{11+[1..2]}:= 1/4 * m{[39,38]};

    aaa[14]{[1..5]}:= 2/1 * m{20+[1..5]};
    aaa[15]{[1..5]}:= 2/1 * m{20+[1..5]};
    aaa[14]{7+[1..4]}:= m{33+[1..4]};
    aaa[15]{7+[1..4]}:= m{[35,34,37,36]};
    aaa[14]{13+[1..2]}:= m{39+[1..2]};
    aaa[15]{13+[1..2]}:= m{[41,40]};
    return aaa;
end;
