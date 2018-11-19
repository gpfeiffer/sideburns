##  the group
G:= Group((1,7,6,5,4,3,2), (2,3,5)(4,7,6));
SetName(G, "7:3");

abcd:= BaseChangesDoubleBurnsideZ(G);

cols:= Concatenation(
               List([1..4], i-> i + 4 * [0..3]), #1
               List([17..20], i-> i + 4 * [0..1]), #3
               List([25..28], i-> [i])
               );

md:= ShrinkByCols(abcd.d, cols);

## Fourier transform
fou:= abcd.chg^0;;
for p in [ [17,18], [19,20], [21,22], [23,24], [25,26], [27,28] ] do
    fou{p}{p}:= [[1,1],[1,-1]];
od;

e:= RightRegularBaseChange(abcd.d, fou);;

me:= ShrinkByCols(e, cols);

