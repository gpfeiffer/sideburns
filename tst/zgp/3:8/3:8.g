##  the group
G:= ZGroup(3, 8, 2);
SetName(G, "3:8");

abcd:= BaseChangesDoubleBurnsideZ(G);

cols:= Concatenation(
               List([1..8], i-> i + 8 * [0..7]), #1
               List([65..70], i-> i + 6 * [0..5]), #2
               List([101..103], i-> i + 3 * [0..2]), #3
               List([110..117], i-> i + 8 * [0..3]), #4
               [[142]],
               List([143..144], i-> i + 2 * [0..1]),
               List([147..154], i-> i + 8 * [0..1]),
               List([163..170], i-> [i])
               );

md:= ShrinkByCols(abcd.d, cols);;

## Fourier transform
fou:= abcd.chg^0;;
#fou4:= List([0..3], i-> List([0..3], j-> E(4)^(i*j)));
#fou4:= fou4^PermutationMat((3,4), 4);
fou2:= List([0..1], i-> List([0..1], j-> E(2)^(i*j)));
for p in [
        [110,111], [112,113], [114,115], [116,117],
        [118,119], [120,121], [122,123], [124,125],
        [126,127], [128,129], [130,131], [132,133],
        [134,135], [136,137], [138,139], [140,141],
        [163,164], [165,166],
        ] do
    fou{p}{p}:= fou2;
od;

fou22:= KroneckerProduct(fou2, fou2);
for p in [
        [147..150], [151..154], [155..158], [159..162],
        [167..170],
        ] do
    fou{p}{p}:= fou22;
od;

e:= RightRegularBaseChange(abcd.d, fou);;

me:= ShrinkByCols(e, cols);;

# separate ...
lst:= abcd.chg^0;
lst[167][163]:= 1;
lst[168][164]:= 1;
lst[167][142]:= 1;
lst[163][142]:= 1;


ee:= RightRegularBaseChange(e, lst);;
mee:= ShrinkByCols(ee, cols);;

dd:= RightRegularBaseChange(ee, fou^-1);;
mdd:= ShrinkByCols(dd, cols);;
