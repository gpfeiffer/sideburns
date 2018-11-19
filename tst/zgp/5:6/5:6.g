##  the group
G:= ZGroup(5, 6, 4);
SetName(G, "5:6");

abcd:= BaseChangesDoubleBurnsideZ(G);

cols:= Concatenation(
               List([1..8], i-> i + 8 * [0..7]), #1
               List([65..68], i-> i + 4 * [0..3]), #2
               List([81..88], i-> i + 8 * [0..3]), #2
               List([113..116], i-> i + 4 * [0..1]), #2
               List([121..124], i-> i + 4 * [0..1]), #2
               List([129..132], i-> i + 4 * [0..1]), #2
               List([93..94], i-> i + 2 * [0..1]), #2
               List([97..100], i-> i + 4 * [0..1]), #2
               List([137..144], i-> [i])
               );

md:= ShrinkByCols(abcd.d, cols);

## Fourier transform
fou:= abcd.chg^0;;

fou2:= List([0..1], i-> List([0..1], j-> E(2)^(i*j)));
fou4:= List([0..3], i-> List([0..3], j-> E(4)^(i*j)));
for p in [
        [85,86], [87,88], [89,90], [91,92],
        [97,98], [99,100], [101,102], [103,104],
        ] do
    fou{p}{p}:= fou2;
od;
#for p in [ [105..108] ] do
#    fou{p}{p}:= KroneckerProduct(fou2, fou2);
#od;
for p in [ [105,107,106,108], [109,110,112,111] ] do
    fou{p}{p}:= fou4;
od;

e:= RightRegularBaseChange(abcd.d, fou);;

me:= ShrinkByCols(e, cols);;
