G:= SymmetricGroup(3);
secs:= RepresentativesSections(G);

GG:= DirectProduct(G, G);
ccs:= ConjugacyClassesSubgroupsDirectProduct(G, G);
subs:= Concatenation(List(ccs, Elements));
tris:= List(subs, x-> TripleSubgroup(GG, x));


cc:= ConjugacyClassesSubgroups(G);
sameTops:= function(x, y)
    return List(Sections(x), TopSec) = List(Sections(y), TopSec);
end;
sameBots:= function(x, y)
    return List(Sections(x), BotSec) = List(Sections(y), BotSec);
end;
    
com:= List(tris, x-> List(tris, y-> Iverson(IsSubgroup(x, y))));
smp:= List(tris, x-> List(tris, y-> Iverson(sameTops(x, y) and IsSubgroup(x, y))));
smk:= List(tris, x-> List(tris, y-> Iverson(sameBots(x, y) and IsSubgroup(x, y))));
smu:= smk^-1 * com * smp^-1;
smt:= TransposedMat(smp);

des:= smk * smu;

mats:= List(tris, x-> RightKappaStar(tris, x));

new:= RightRegularBaseChange(mats, des);

#c:= RightRegularBaseChange(new, smt);;


#  these numbers compare sizes of Conjugators ...
vec:= List(des, x-> 1);
for i in [25..30] do vec[i]:= 2; od;
for i in [31..36] do vec[i]:= 6; od;
for i in [49..52] do vec[i]:= 6; od;
dia:= DiagonalMat(1/6*vec);

c:= RightRegularBaseChange(new, dia * smt);;

# this has to do with the radical
min:= des^0;
for i in [55..60] do min[i][52]:= 1; min[i][36]:= 1; od;

d:= RightRegularBaseChange(c, min^-1);;

chg:= des * dia * smt / min;

cols:= [
        [1, 7, 8, 9, 25, 31],
        [2, 10, 11, 12, 26, 32],
        [3, 13, 14, 15, 27, 33],
        [4, 16, 17, 18, 28, 34],
        [5, 19, 20, 21, 29, 35],
        [6, 22, 23, 24, 30, 36],
        [37, 40, 43, 49],
        [38, 41, 44, 50],
        [39, 42, 45, 51],
        [46, 47, 48, 52],
        [53],
        [54],
        [55],
        [56],
        [57],
        [58],
        [59],
        [60]
        ];

shrink:= des^0;
for col in cols do
    i:= col[1];
    for j in col do
        shrink[j][i]:= -1;
        shrink[i][j]:= 1;
    od;
od;

firsts:= List(cols, x-> x[1]);
seconds:= Difference([1..60], firsts);
poss:= Concatenation(firsts, seconds);

