#############################################################################
##
#W  marks.g                                                         sideBurns
##
#W  by Götz Pfeiffer <goetz.pfeiffer@nuigalway.ie>
##
#Y  Copyright (C) 2016  Götz Pfeiffer
##
#Y  This file is part of the GAP 4 _sideBurns_ package.  _sideBurns_ is
#Y  free software, see the license information at the end of this file.
##

#############################################################################
##
##  Marks of Direct Products.
##
##  This file contains structures and functions for the computation of the
##  table of marks of a direct product of finite groups.
##
##  The table of marks of a direct product is a product of three simpler
##  matrices, each of which can be derived from information about the
##  constituent groups.
##

#############################################################################
##
##  TriplesByTop
##
##  Sort conjugacy classes of subgroups of GxH by top groups of sections.
##  This yields a rectangular arrangement, with rows labelled by the
##  subgroup classes of G, and columns labelled by those of H.
##
##  FIXME: this could probably be more efficient and less repetitive ...
##
TriplesByTop:= function(G, H)
    local   ccG,  ccH,  mat,  tri,  tops,  posG,  posH;
    ccG:= ConjugacyClassesSubgroups(G);
    ccH:= ConjugacyClassesSubgroups(H);
    mat:= List(ccG, x-> List(ccH, y-> []));
    for tri in Flat(List(TriplesDirectProduct(G, H), x-> x.trip)) do
        tops:= List(Sections(tri), TopSec);
        posG:= PositionProperty(ccG, c-> tops[1] in c);
        posH:= PositionProperty(ccH, c-> tops[2] in c);
        Add(mat[posG][posH], tri);
    od;
    return mat;
end;

#############################################################################
##
##  TriplesByBot
##
##  Sort conjugacy classes of subgroups of GxH by bottom groups of sections.
##  This yields a rectangular arrangement, with rows labelled by the
##  subgroup classes of G, and columns labelled by those of H.
##
##  FIXME: this could probably be more efficient and less repetitive ...
##
TriplesByBot:= function(G, H)
    local   ccG,  ccH,  mat,  tri,  bots,  posG,  posH;
    ccG:= ConjugacyClassesSubgroups(G);
    ccH:= ConjugacyClassesSubgroups(H);
    mat:= List(ccG, x-> List(ccH, y-> []));
    for tri in Flat(List(TriplesDirectProduct(G, H), x-> x.trip)) do
        bots:= List(Sections(tri), BotSec);
        posG:= PositionProperty(ccG, c-> bots[1] in c);
        posH:= PositionProperty(ccH, c-> bots[2] in c);
        Add(mat[posG][posH], tri);
    od;
    return mat;
end;

#############################################################################
##
##  Same Top Matrix.
##
##  FIXME: get rid of GxH
##
TopClassIncMatDirectProduct:= function(G, H)
    local   aligned,  subs,  mats,  sub,  secs,  PP,  NN,  GxH,  em;

    # how to align two triples for the same tops
    aligned:= function(GH, em, PP, trip)
        local   tops,  a;
        tops:= List(Sections(trip), TopSec);
        a:= List([1,2], i-> RepresentativeAction(GH[i], tops[i], PP[i]));
        a:= a[1]^em[1] * a[2]^em[2];
        return trip^a;
    end;

    subs:= Concatenation(TriplesByTop(G, H));
    mats:= [];
    for sub in subs do

        # special case first
        if Length(sub) = 1 then
            Add(mats, [[1]]);
        else
            secs:= Sections(sub[1]);
            PP:= List(secs, TopSec);
            NN:= [Normalizer(G, PP[1]), Normalizer(H, PP[2])];

            # make direct product of normalizers in GxH
            GxH:= ProductGroup(sub[1]);
            em:= List([1,2], i-> Embedding(GxH, i));
            NN:= List([1,2], i-> Image(em[i], NN[i]));
            NN:= ClosureGroup(NN[1], NN[2]);

            # align subs to have tops PP
            sub:= List(sub, s-> aligned([G,H], em, PP, s));
            Add(mats, ClassIncMat(NN, sub, OnPoints, IsSubgroup));
        fi;
    od;

    return rec(
      reps:= Concatenation(subs),
      mats:= mats
    );
end;


#############################################################################
##
##  Same Bottom Matrix
##
##  FIXME: get rid of GxH
##
##
BotClassIncMatDirectProduct:= function(G, H)
    local   aligned,  subs,  mats,  sub,  secs,  KK,  NN,  GxH,  em;

    # how to align two triples for the same bots
    aligned:= function(GH, em, KK, trip)
        local   bots,  a;
        bots:= List(Sections(trip), BotSec);
        a:= List([1,2], i-> RepresentativeAction(GH[i], bots[i], KK[i]));
        a:= a[1]^em[1] * a[2]^em[2];
        return trip^a;
    end;

    subs:= Concatenation(TriplesByBot(G, H));
    mats:= [];
    for sub in subs do

        # special case first
        if Length(sub) = 1 then
            Add(mats, [[1]]);
        else
            secs:= Sections(sub[1]);
            KK:= List(secs, BotSec);
            NN:= [Normalizer(G, KK[1]), Normalizer(H, KK[2])];

            # make direct product of normalizers in GxH
            GxH:= ProductGroup(sub[1]);
            em:= List([1,2], i-> Embedding(GxH, i));
            NN:= List([1,2], i-> Image(em[i], NN[i]));
            NN:= ClosureGroup(NN[1], NN[2]);

            # align subs to have tops PP
            sub:= List(sub, s-> aligned([G,H], em, KK, s));
            Add(mats, ClassIncMat(NN, sub, OnPoints, IsSubgroup));
        fi;
    od;

    return rec(
      reps:= Concatenation(subs),
      mats:= mats
    );
end;


#############################################################################
##
##  IsoClassIncMatDirectProduct
##
##  FIXME: get rid of GxH
##
IsoClassIncMatDirectProduct0:= function(G, H)
    local   subs,  mats,  GxH,  sub;

    subs:= List(TriplesDirectProduct(G, H), x-> Flat(x.trip));
    mats:= [];
    GxH:= DirectProduct(G, H);
    for sub in subs do
        Add(mats, ClassIncMat(GxH, sub, OnPoints, IsDescendantTriple));
    od;

    return rec(
      reps:= Concatenation(subs),
      mats:= mats
    );
end;

IsoClassIncMatDirectProduct:= function(G, H)
    local   subs,  mats,  sub,  trips,  dblcs,  morG,  dblG,  rows,
            grp,  morH,  dblH,  cols,  mat,  poss,  inverse,  i;

    subs:= TriplesDirectProduct(G, H);

    mats:= [];
    for sub in subs do
        trips:= Flat(sub.trip);
        dblcs:= Flat(List([1..Length(sub.rows)], i-> List([1..Length(sub.cols)], j-> List(sub.dblc[i][j], x-> rec(i:= i,j:= j, dc:= x)))));


        ## FIXME: special case rows = cols?

        # the rows correspond to G
        morG:= Concatenation(List(sub.rows, x-> AllMorphismsSection(x!.section)));
        dblG:= Flat(List([1..Length(sub.rows)], i-> List(CosetsSection(sub.rows[i]!.section), c-> rec(i:= i, cos:= c))));

        rows:= List(morG, TripleMorphism);
        grp:= ProductGroup(rows[1]);
        rows:= ClassIncMat(grp, rows, OnPoints, IsSubgroup);

        # the cols correspond to H
        morH:= Concatenation(List(sub.cols, x-> AllMorphismsSection(x!.section)));
        dblH:= Flat(List([1..Length(sub.cols)], j-> List(CosetsSection(sub.cols[j]!.section), c-> rec(j:= j, cos:= c))));
        cols:= List(morH, TripleMorphism);
        grp:= ProductGroup(cols[1]);
        cols:= ClassIncMat(grp, cols, OnPoints, IsSubgroup);

        # Kronecker product
        mat:= KroneckerProduct(rows, cols);
        poss:= Flat(List(dblG, x-> List(dblH, y-> PositionProperty(dblcs, dc-> x.i = dc.i and y.j = dc.j and Representative(x.cos)/Representative(y.cos) in dc.dc))));
        inverse:= List(Set(poss), x-> []);
        for i in [1..Length(poss)] do Add(inverse[poss[i]], i); od;
        mat:= List(inverse, x-> Sum(mat{x}));
        poss:= List(inverse, x-> x[1]);
        mat:= List(mat, x-> x{poss});

        Add(mats, mat);
    od;

    return rec(
      reps:= Concatenation(List(subs, x-> Flat(x.trip))),
      mats:= mats
    );
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
