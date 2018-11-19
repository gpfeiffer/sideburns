##  the group
G:= DihedralGroup(IsPermGroup, 22);
SetName(G, "D22");

abcd:= BaseChangesDoubleBurnsideZ(G);

cols:= Concatenation(
               List([1..4], i-> i + 4 * [0..3]), #1
               List([17..18], i-> i + 2 * [0..1]), #3
               List([21..30], i-> [i])
               );

md:= ShrinkByCols(abcd.d, cols);

## Fourier transform
fou:= abcd.chg^0;;
fou5:= List([0..4], i-> List([0..4], j-> E(5)^(i*j)));
for p in [ [21..25], [26,27,29,28,30] ] do  # swap 28 and 29 !?
    fou{p}{p}:= fou5;
od;

e:= RightRegularBaseChange(abcd.d, fou);;

me:= ShrinkByCols(e, cols);

