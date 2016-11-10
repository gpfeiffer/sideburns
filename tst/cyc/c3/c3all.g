G:= Group((1,2,3));


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

des:= smk * smu * smt;

mats:= List(tris, x-> RightKappaStar(tris, x));

new:= RightRegularBaseChange(mats, des);

min:= des^0;
min[6][4]:= 1/2;
min[5][4]:= 1/2;


min[3][3]:= 1/2;
min[4][4]:= 1/2;

c:= RightRegularBaseChange(mats, des/min);


# the ultimate solution: matrix units!
vec:= List(mats, x-> 1);
for i in [3,4] do
    vec[i]:= 2;
od;
dia:= DiagonalMat(1/3*vec);

chg:= smk*smu*dia*smt;

d:= RightRegularBaseChange(mats, chg);
