G:= Group((1,2,3,4));

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

des:= smk * smu;

mats:= List(tris, x-> RightKappaStar(tris, x));

new:= RightRegularBaseChange(mats, des);


smt:= TransposedMat(smp);
min:= des^0;

min[10][5]:= 1;
min[11][6]:= 1;
min[12][8]:= 1;
min[13][9]:= 1;

min[14][13]:= 1;
min[15][13]:= 1;
min[13][13]:= 2;
min[12][12]:= 2;

min[14][9]:= 1;
min[15][9]:= 1;


c:= RightRegularBaseChange(new, min);;


# the ultimate solution: matrix units!
vec:= List(mats, x-> 1);
for i in [7,8,9,12,13] do
    vec[i]:= 2;
od;
dia:= DiagonalMat(1/4*vec);

d:= RightRegularBaseChange(mats, des*dia*smt);
