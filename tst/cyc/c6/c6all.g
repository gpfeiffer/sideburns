G:= Group((1,2,3,4,5,6));

ccs:= ConjugacyClassesSubgroupsDirectProduct(G, G);
subs:= Concatenation(List(ccs, Elements));
GG:= DirectProduct(G, G);
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

c:= RightRegularBaseChange(new, smt);;


#  these numbers compare sizes of Conjugators ...
vec:= List(des, x-> 1);
for i in [9..12] do vec[i]:= 2; od;
for i in [13..16] do vec[i]:= 2; od;
for i in [19..20] do vec[i]:= 2; od;
dia:= DiagonalMat(1/6*vec);

d:= RightRegularBaseChange(new, dia * smt);;

chg:= des * dia * smt;
