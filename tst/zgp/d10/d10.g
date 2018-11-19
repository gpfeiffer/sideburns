##  the group
G:=  DihedralGroup(IsPermGroup, 10);
SetName(G, "D10");

abcd:= BaseChangesDoubleBurnsideZ(G);

cols:= Concatenation(
               List([1..4], i-> i + 4 * [0..3]), #1
               List([17..18], i-> i + 2 * [0..1]), #2
               List([21..24], i-> [i])
               );

ma:= ShrinkByCols(abcd.a, cols);
md:= ShrinkByCols(abcd.d, cols);

fou:= abcd.chg^0;
for p in [ [21,22], [23,24] ] do
    fou{p}{p}:= [[1,1],[1,-1]];
od;

e:= RightRegularBaseChange(abcd.d, fou);;

me:= ShrinkByCols(e, cols);;
