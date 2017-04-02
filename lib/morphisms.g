#############################################################################
##
#W  morphisms.g                                                     sideBurns
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
##  Morphisms.
##
##  This file contains structures and functions for morphisms of sections of
##  a finite group.
##
##  A _morphism_, in this context is an isomorphism between a section
##  and a standard representative permutation group for the isomorphism
##  class of the section's quotient group.
##

#############################################################################
##
##  Representation and Type
##
DeclareRepresentation("IsMorphism",
        IsComponentObjectRep and IsAttributeStoringRep and IsPermGroup,
        ["section", "group", "theta"]);

MorphismType:=  NewType(CollectionsFamily(PermutationsFamily), IsMorphism);


#############################################################################
##
#C  Morphism( section, group, theta ) . . . . . . . . . . . . .  constructor
##
##  returns a new morphism object...
##
Morphism:= function(section, group, theta)
    local   r;
    r:= rec(section:= section, group:= group, theta:= theta);
    return Objectify(MorphismType, r);
end;


##  print method
InstallMethod(PrintObj, "for morphisms", true, [IsMorphism], 0, function(gt)
    Print("Morphism( ", gt!.section, ", ", gt!.group, ", ", gt!.theta,  " )");
end);

InstallMethod(ViewObj, "for morphisms", true, [IsMorphism], 0, function(gt)
    Print("<morphism>");
end);


#############################################################################
##
##  Automizer( mor )
##
DeclareAttribute("Automizer", IsMorphism);

InstallMethod(Automizer, "for a morphism", [IsMorphism], function(mor)
    local   sec,  P,  map,  map1,  auto,  gens,  inv,  reps,  imgs,
            conj,  hom;

    sec:= mor!.section;
    P:= TopSec(sec);
    map:= NaturalHomomorphism(sec) * mor!.theta;
    map1:= InverseGeneralMapping(map);
    auto:= Automizer(sec);
    inv:= InverseHomomorphism(auto);
    gens:= GeneratorsOfGroup(AsGroup(auto));
    reps:= List(gens, x-> Representative(Images(inv, x)));
    imgs:= List(reps, x-> map1 * ConjugatorAutomorphism(P, x) * map);
    conj:= Group(imgs);
    hom:= GroupHomomorphismByImages(AsGroup(auto), conj, gens, imgs);

    return Morphism(auto, conj, hom);
end);


##############################################################################
##
##  Conjugators( mor )
##
##  The subgroup of the automorphism group of $U$ induced by conjugation.
##
Conjugators:= function(mor)
    return Automizer(mor)!.group;
end;


#############################################################################
##
#F  OneMorphismSection( section )
#A  AllMorphismsSection( section )
##
##  Computes a list of all monomorphisms from a permutation group isomorphic
##  to the quotient of section into the section up to conjugation in G.
##
DeclareAttribute("OneMorphismSection", IsSection);

InstallMethod(OneMorphismSection, "for a section", [IsSection], function(sec)
    local   U,  iso;
    U:= Image(IsomorphismPermGroup(SmallGroup(IdSection(sec))));
    iso:= IsomorphismGroups(AsGroup(sec), U);
    return Morphism(sec, U, iso);
end);

DeclareAttribute("AllMorphismsSection", IsSection);

# FIXME: use this function below.  make it an attribute?
CosetsSection:= function(sec)
    local   mor,  U,  C;
    mor:= OneMorphismSection(sec);
    U:= mor!.group;
    C:= Conjugators(mor);
    return RightCosets(AutomorphismGroup(U), C);
end;


InstallMethod(AllMorphismsSection, "for a section", [IsSection], function(sec)
    local   mor,  U,  C,  A,  cos;

    mor:= OneMorphismSection(sec);
    U:= mor!.group;

    C:= Conjugators(mor);
    A:= AutomorphismGroup(U);
    cos:= List(RightCosets(A, C), Representative);
##  FIXME: use mor * auto
    return List(cos, d-> Morphism(sec, U, mor!.theta*d));
end);


#############################################################################
SubgroupMorphism:= mor -> SubgroupTriple(TripleMorphism(mor));


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
