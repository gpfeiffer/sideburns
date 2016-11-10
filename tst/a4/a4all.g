G:= AlternatingGroup(4);
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

c:= RightRegularBaseChange(new, smt);;


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

c:= RightRegularBaseChange(c, min^-1);;

chg:= des * dia * smt / min;
