#############################################################################
##
##
##
RightRegularBaseChange:= function(basis, ch)
    local   ch1,  iii;
    ch1:= ch^-1;
    iii:= [1..Length(basis)];
    return List(iii, j-> ch1 * Sum(iii, l-> ch1[j][l] * basis[l]) * ch);
end;


#############################################################################
##
##
supp:= l -> Filtered([1..Length(l)], i-> l[i] <> 0*l[i]);

# represent sparse matrix compactly
compact:= function(a)
    local   map,  i,  x;
    map:= rec();
    for i in [1..Length(a)] do
        x:= supp(a[i]);
        if x <> [] then
            map.(i):= [x, a[i]{x}];
        fi;
    od;
    return map;
end;


#############################################################################
##
##
##
pretty:= function(a)
    local   c,  n,  more,  i;
    c:= compact(a);
    for n in Set(RecFields(c), Int) do
        Print(n, " -> ");
        more:= false;
        for i in [1..Length(c.(n)[1])] do
            if more then Print(" + "); fi;
            if c.(n)[2][i] <> 1 then
                Print(c.(n)[2][i], "*");
            fi;
            Print(c.(n)[1][i]);
            more:= true;
        od;
        Print("\n");
    od;
end;


#############################################################################
##
##  display a sparse matrix with single digit entries ...
##
short:= function(m)
    local l, a;
    for l in m do
        for a in l do
            if a = 0 then
                Print(". \c");
                elif a = -1 then
                  Print("- \c");
            else
                Print(a, " \c");
            fi;
        od;
        Print("\n");
    od;
end;
