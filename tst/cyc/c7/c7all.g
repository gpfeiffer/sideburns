G:= Group((1,2,3,4,5,6,7));

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

des:= smk * smu * TransposedMat(smp);

mats:= List(tris, x-> RightKappaStar(tris, x));

new:= RightRegularBaseChange(mats, des);

min:= des^0;
min[10][4]:= 5/6;
min[9][4]:= 5/6;
min[8][4]:= 5/6;
min[7][4]:= 5/6;
min[6][4]:= 5/6;
min[5][4]:= 5/6;

c:= RightRegularBaseChange(mats, des/min);

min[3][3]:= 1/6;
min[4][4]:= 1/6;

d:= RightRegularBaseChange(mats, des/min);
