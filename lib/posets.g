#############################################################################
##
#W  posets.g                                                        sideBurns
##
#W  Götz Pfeiffer <goetz.pfeiffer@nuigalway.ie>
##
#Y  Copyright (C) 2016  Götz Pfeiffer
##
#Y  This file is part of the GAP 4 _sideBurns_ package.  _sideBurns_ is
#Y  free software, see the license information at the end of this file.
##

#############################################################################
##
##  $G$-Posets.
##
##  This file contains structures and functions for partially ordered sets
##  acted upon by a finite group.
##
##  Of particular interest is the so-called class incidence matrix,
##  that is the incidence matrix of the poset, collapsed under the group
##  action with one row and one column per orbit, counting incidences
##  between orbits.
##

#############################################################################
##
##
##
ClassIncMat:= function(G, reps, action, contains)
    local   mat,  x,  orb;
    mat:= [];
    for x in reps do
        orb:= Orbit(G, x, action);
        #Print("[", Length(orb), "]\c");
        Add(mat, List(reps, lo-> Number(orb, hi-> contains(hi, lo))));
    od;
    return mat;
end;


#############################################################################
##
##  MacNeille Completion
##
##  Suppose a poset is given in the form of its incidence matrix
##  with entries 0 and 1,  compute its completion in terms of cuts.
##
MacNeilleCompletion:= function(poset)
    local   iii,  lower,  upper,  orb,  set,  i,  new;

    # convert poset into lower set
    iii:= [1..Length(poset)];
    lower:= List(iii, x-> Filtered(iii, y-> poset[x][y] = 1));
    upper:= List(iii, x-> Filtered(iii, y-> poset[y][x] = 1));

    # orbit algorithm
    orb:= [iii];
    for set in orb do
        for i in iii do
            new:= Intersection(set, upper[i]);
            if not new in orb then
                Add(orb, new);
            fi;
        od;
    od;

    # filter
    orb:= Filtered(orb, x-> x = Intersection(upper{Intersection(lower{x})}));

    # reduce to minima
    orb:= List(orb, x-> Filtered(x, i-> x<> Union(upper{Difference(x, upper[i])})));

    return orb;
end;


#############################################################################
##
##  License Information.
##
##  _sideBurns_ is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##
##  _sideBurns_ is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with _sideBurns_.  If not, see <http://www.gnu.org/licenses/>.
##
