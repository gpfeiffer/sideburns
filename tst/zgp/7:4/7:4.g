##  the group
G:= ZGroup(7, 4, 6);
SetName(G, "7:4");

abcd:= BaseChangesDoubleBurnsideZ(G);

cols:= Concatenation(
               List([1..6], i-> i + 6 * [0..5]), #1
               List([37..40], i-> i + 4 * [0..3]), #2
               List([53..56], i-> i + 4 * [0..1]),
               List([61..66], i-> i + 6 * [0..1]),
               List([73..84], i-> [i])
               );

md:= ShrinkByCols(abcd.d, cols);;

## Fourier transform
fou:= abcd.chg^0;;
fou3:= List([0..2], i-> List([0..2], j-> E(3)^(i*j)));
fou2:= List([0..1], i-> List([0..1], j-> E(2)^(i*j)));
#fou6:= List([0..5], i-> List([0..5], j-> E(6)^(i*j)));
fou6:= KroneckerProduct(fou2, fou3);
for p in [
        [53,54], [55,56], [57,58], [59,60],
        ] do
    fou{p}{p}:= fou2;
od;
for p in [
        [61,62,63], [64,65,66], [67,68,69], [70,71,72],
        [73,74,75], [76,77,78],
        ] do
    fou{p}{p}:= fou3;
od;
for p in [
        [79..84]
        ] do
    fou{p}{p}:= fou6;
od;

e:= RightRegularBaseChange(abcd.d, fou);;

me:= ShrinkByCols(e, cols);;

# separate ...
lst:= abcd.chg^0;
lst[79][73]:= 1;
lst[80][74]:= 1;
lst[81][75]:= 1;

ee:= RightRegularBaseChange(e, lst);;
mee:= ShrinkByCols(ee, cols);;

# undo fourier transform
dd:= RightRegularBaseChange(ee, fou^-1);;
mdd:= ShrinkByCols(dd, cols);;
