#############################################################################
##
#W  sectlatt.g                                                      sideBurns
##
#W  by Götz Pfeiffer <goetz.pfeiffer@nuigalway.ie>
##
#Y  Copyright (C) 2016  NUI Galway
##
#Y  This file is part of the GAP 4 _sideBurns_ package.  _sideBurns_ is
#Y  free software, see the license information at the end of this file.
##

#############################################################################
##
##  The Sections Lattice.
##
##  This file contains functions for the lattice of sections of a
##  finite group.
##
##  The sections of a finite group $G$ form lattices,
##  in different ways ...
##


#############################################################################
##
##  Iverson( bool )
##
##  converts `true`/`false` into 1/0.
##
Iverson:= function(property)
    if property then return 1; else return 0; fi;
end;

#############################################################################
##
##  The meet in the descendent lattice
##
MeetSections:= function(secL, secR)
    local   G,  top,  one,  two;

    G:= secL!.G;
    if secR!.G <> G then
        Error("Parent groups are not identical");
    fi;

    top:= Intersection(TopSec(secL), TopSec(secR));
    one:= Intersection(TopSec(secL), BotSec(secR));
    two:= Intersection(BotSec(secL), TopSec(secR));

    return Section(G, top, ClosureGroup(one, two));
end;


#############################################################################
##
##  The join in the descendent lattice
##
JoinSections:= function(secL, secR)
    local   G,  top,  bot;

    G:= secL!.G;
    if secR!.G <> G then
        Error("Parent groups are not identical");
    fi;

    top:= ClosureGroup(TopSec(secL), TopSec(secR));
    bot:= Intersection(BotSec(secL), BotSec(secR));

    return Section(G, top, Core(top, bot));
end;


#############################################################################
##
##  componentwise intersection
##
IntersectionSections:= function(secL, secR)
    local   G,  top,  bot;

    G:= secL!.G;
    if secR!.G <> G then
        Error("Parent groups are not identical");
    fi;

    top:= Intersection(TopSec(secL), TopSec(secR));
    bot:= Intersection(BotSec(secL), BotSec(secR));

    return Section(G, top, bot);
end;


#############################################################################
##
##  componentwise closure
##
ClosureSections:= function(secL, secR)
    local   G,  top,  bot;

    G:= secL!.G;
    if secR!.G <> G then
        Error("Parent groups are not identical");
    fi;

    top:= ClosureGroup(TopSec(secL), TopSec(secR));
    bot:= ClosureGroup(BotSec(secL), BotSec(secR));

    return Section(G, top, NormalClosure(top, bot));
end;


#############################################################################
##
##  descending sub: P' <= P and K \cap P' = K \cap K'
##
IsDescendantSection:= function(hi, lo)
    if hi!.G <> lo!.G then
        Error("Parent groups are not identical");
    fi;

    return  IsSubgroup(TopSec(hi), TopSec(lo)) and
            IsSubgroup(BotSec(lo), Intersection(BotSec(hi), TopSec(lo)));
end;


#############################################################################
##
##  squeeze sub: P' <= P and K' >= K
##
##  A subsection is a section of a (quotient group of a) section ...
##
IsSubsection:= function(hi, lo)
    if hi!.G <> lo!.G then
        Error("Parent groups are not identical");
    fi;

    return  IsSubgroup(TopSec(hi), TopSec(lo)) and
            IsSubgroup(BotSec(lo), BotSec(hi));
end;


#############################################################################
##
##  component wise sub: P' <= P and K' <= K
##
IsCompSubsection:= function(hi, lo)
    if hi!.G <> lo!.G then
        Error("Parent groups are not identical");
    fi;

    return  IsSubgroup(TopSec(hi), TopSec(lo)) and
            IsSubgroup(BotSec(hi), BotSec(lo));
end;


#############################################################################
##
##  IntermediateSections( secL, secR )
##
##  Computes two new sections from the given ones.
##  Under suitable circumstances this yields a chain for a particular order.
##
##  FIXME:  compute the corresponding homomorphism
##
IntermediateSections:= function(hi, lo)
    local   G,  K,  P1,  P1K,  int;
    G:= hi!.G;  K:= hi!.K;  P1:= lo!.P;
    P1K:= ClosureGroup(P1, K);
    int:= Intersection(P1, K);
    return [hi, Section(G, P1K, K), Section(G, P1, int), lo];
end;


#############################################################################
##
##  Right regular Matrices for Meet and Intersection
##
RightMeet:= function(secs, secR)
    local   mat,  secL,  meet;
    mat:= [];
    for secL in secs do
        meet:= MeetSections(secL, secR);
        Add(mat, List(secs, x-> Iverson(x = meet)));
    od;
    return mat;
end;

RightIntersection:= function(secs, secR)
    local   mat,  secL,  meet;
    mat:= [];
    for secL in secs do
        meet:= IntersectionSections(secL, secR);
        Add(mat, List(secs, x-> Iverson(x = meet)));
    od;
    return mat;
end;


#############################################################################
##
##  ClassIncMatSections
##
##  the class incidence matrix of componentwise subsections
##
ClassIncMatSections:= function(G)
    local   reps;
    reps:= RepresentativesSections(G);
    return rec(
      reps:= reps,
      mat:= ClassIncMat(G, reps, OnPoints, IsCompSubsection)
    );
end;


#############################################################################
##
##  TopClassIncMatSections
##
##  the class incidence matrix of the same top incidence
##
TopClassIncMatSections:= function(G)
    local   secs,  mats,  sec,  P,  N;

    secs:= SectionsByTop(G);
    mats:= [];
    for sec in secs do
        P:= TopSec(sec[1]);
        N:= Normalizer(G, P);

        # align secs to have top P
        sec:= List(sec, s-> s^RepresentativeAction(G, TopSec(s), P));
        Add(mats, ClassIncMat(N, sec, OnPoints, IsCompSubsection));
    od;


    return rec(
      reps:= Concatenation(secs),
      mats:= mats
    );
end;


#############################################################################
##
##  BotClassIncMatSections
##
##  the class incidence matrix of the same bottom incidence
##
BotClassIncMatSections:= function(G)
    local   secs,  mats,  sec,  K,  N;

    secs:= SectionsByBot(G);
    mats:= [];
    for sec in secs do
        K:= BotSec(sec[1]);
        N:= Normalizer(G, K);

        # align secs to have bot K
        sec:= List(sec, s-> s^RepresentativeAction(G, BotSec(s), K));
        Add(mats, ClassIncMat(N, sec, OnPoints, IsCompSubsection));
    od;

    return rec(
      reps:= Concatenation(secs),
      mats:= mats
    );
end;


#############################################################################
##
##  IsoClassIncMatSection
##
##  the class incidence matrix of the isomorphic sections incidence
##
##  Note that some of the isomorphic sections posets could be further
##  decomposed as a direct sum of posets ...
##
IsoClassIncMatSections:= function(G)
    local   secs,  mats,  sec;
    secs:= List(SectionsByType(G), x-> x.sections);
    mats:= [];
    for sec in secs do
        Add(mats, ClassIncMat(G, sec, OnPoints, IsDescendantSection));
    od;

    return rec(
      reps:= Concatenation(secs),
      mats:= mats
    );
end;


#############################################################################
##
##
##
ProductClassIncMats:= function(G)
    local   secs,  top,  bot,  iso,  tops,  bots,  isos,  t,  b,  i;

    secs:= RepresentativesSections(G);
    top:= TopClassIncMatSections(G);
    bot:= BotClassIncMatSections(G);
    iso:= IsoClassIncMatSections(G);
    tops:= List(secs, x-> Position(top.reps, x));
    bots:= List(secs, x-> Position(bot.reps, x));
    isos:= List(secs, x-> Position(iso.reps, x));
    t:= DirectSumMat(top.mats);
    b:= DirectSumMat(bot.mats);
    i:= DirectSumMat(iso.mats);

    return b{bots}{bots} * i{isos}{isos} * t{tops}{tops};
end;


#############################################################################
##
##  DesClassIncMatSections
##
##  FIXME: for bigger examples it might be faster to compute this as a product
##
DesClassIncMatSections:= function(G)
    return ClassIncMat(G, RepresentativesSections(G),
                   OnPoints, IsDescendantSection);
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
