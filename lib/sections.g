#############################################################################
##
#W  sections.g                                                      sideBurns
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
##  Sections.
##
##  This file contains structures and functions for sections of a
##  finite group.
##
##  A _section_ (or subquotient) of a finite group $G$ with
##  subgroups $P$ and $K$ such that $K$ is normal in $P$.
##


#############################################################################
##
##  Family, Representation, Type.
##
SectionFamily:= NewFamily("SectionFamily", IsObject);

DeclareRepresentation("IsSection",
        IsComponentObjectRep and IsAttributeStoringRep,
        ["G", "P", "K"]);

SectionType:= NewType(SectionFamily, IsSection);


#############################################################################
##
#C  Section( G, P, K ) . . . . . . . . . . . . . . . . . . constructor.
##
##  constructs a new `Section` object with components `G`, `P`, and `K`.
##
Section:= function(G, P, K)
    local   r;
    r:= rec(G:= G, P:= P, K:= K);
    return Objectify(SectionType, r);
end;

#############################################################################
##
##  print method
##
InstallMethod(PrintObj, "for sections", true, [IsSection], 0, function(sec)
    Print("Section( ", sec!.G, ", ", sec!.P, ", ", sec!.K, " )");
end);


##############################################################################
##
#O  TopSec( section )
##
##  gets the top group `P` of the given `section`.
##
DeclareOperation("TopSec", [IsSection]);

InstallMethod(TopSec, "for a section", [IsSection], sec -> sec!.P);

##############################################################################
##
#O  BotSec( section )
##
##  gets the bottom group `K` of the given `section`.
##
DeclareOperation("BotSec", [IsSection]);

InstallMethod(BotSec, "for a section", [IsSection], sec -> sec!.K);

##############################################################################
##
##  section ^ g
##
InstallOtherMethod(\^, "for a section", true, [IsSection, IsMultiplicativeElementWithInverse], 0, function(sec, g)
    return Section((sec!.G)^g, (sec!.P)^g, (sec!.K)^g);
end);

#############################################################################
##
##  sectionL = sectionR
##  sectionL < sectionR
##
InstallMethod(\=, "for sections", true, [IsSection, IsSection], 0, function(sec1, sec2)
    return ForAll(["G", "P", "K"], p-> sec1!.(p) = sec2!.(p));
end);

InstallMethod(\<, "for sections", true, [IsSection, IsSection], 0, function(secL, secR)
    if secL!.G = secR!.G then
        if secL!.P = secR!.P then
            return secL!.K < secR!.K;
        else
            return secL!.P < secR!.P;
        fi;
    else
        return secL!.G < secR!.G;
    fi;
end);


#############################################################################
##
##  Normalizer and homomorphism
##
#A  HomomorphismSection( section )
#A  NarrowHomomorphism( section )
#A  InverseHomomorphism( section )
##
##  returns the natural homomorphism $\alpha: N_G(K) \to N_G(K)/K$,
##  its restriction to $P$,
##  or the generalized mapping which is the inverse of that restriction.
##
DeclareAttribute("HomomorphismSection", IsSection);

InstallMethod(HomomorphismSection, "for a section", [IsSection], function(sec)
    local   N;
    N:= Normalizer(sec!.G, sec!.K);
    return NaturalHomomorphismByNormalSubgroup(N, sec!.K);
end);

##  its restriction to P ...
##
DeclareAttribute("NarrowHomomorphism", IsSection);

InstallMethod(NarrowHomomorphism, "for a section", [IsSection], function(sec)
    local   phi,  P,  gens;

    #    return RestrictedMapping(HomomorphismSection(sec), sec!.P);
    phi:= HomomorphismSection(sec);
    P:= TopSec(sec);
    gens:= GeneratorsOfGroup(P);
    return GroupHomomorphismByImages(P, Image(phi, P),
                   gens, List(gens, x-> x^phi));
end);

##  and its inverse ...
##
DeclareAttribute("InverseHomomorphism", IsSection);

InstallMethod(InverseHomomorphism, "for a section", [IsSection], function(sec)
    return InverseGeneralMapping(NarrowHomomorphism(sec));
end);

##############################################################################
##
#A  NormalizerSection( section )
##
##  The _normalizer_ of a section is the intersection of the normalizers
##  of its top and bottom groups.
##
DeclareAttribute("NormalizerSection", IsSection);

InstallMethod(NormalizerSection, "for a section", [IsSection], function(sec)
    return Intersection(Normalizer(sec!.G, sec!.P), Normalizer(sec!.G, sec!.K));
end);


##############################################################################
##
#A  CentralizerSection( section )
##
##  The _centralizer_ of a section is the preimage in $N_G(K)$ of
##  the centralizer of the quotient group $P/K$ in $N_G(K)/K$.
##
DeclareAttribute("CentralizerSection", IsSection);

InstallMethod(CentralizerSection, "for a section", [IsSection], function(sec)
    local hom;
    hom:= HomomorphismSection(sec);
    return PreImages(hom, Centralizer(Image(hom), Image(hom, sec!.P)));
end);


#############################################################################
##
#A  Automizer( section )
##
##  The automizer of a section of a group $G$ is the section of $G$
##  consisting of its normalizer and its centralizer.
##
DeclareAttribute("Automizer", IsSection);

InstallMethod(Automizer, "for a section", [IsSection], function(sec)
    return Section(sec!.G, NormalizerSection(sec), CentralizerSection(sec));
end);


##############################################################################
##
#A  AsGroup( section )
##
##  Returns the quotient group $P/K$ of a section.
##
DeclareAttribute("AsGroup", IsSection);

InstallMethod(AsGroup, "for a section", [IsSection], function(sec)
    return Image(NarrowHomomorphism(sec));
end);


#############################################################################
##
#A  Size( section )
##
##  The size of a section is the size of its quotient group.
##
DeclareAttribute("Size", IsSection);

InstallMethod(Size, "for a section", [IsSection], function(sec)
    return Size(AsGroup(sec));
end);


#############################################################################
##
#F  IdSection( section )
##
##  The ID of a section is the ID of its quotient group.
##
IdSection:= function(sec)
    return IdGroup(AsGroup(sec));
end;


##############################################################################
##
#A  RepresentativesSections( G )
##
##  computes a list of representatives  of the conjugacy classses of sections
##  of the finite group $G$.
##
DeclareAttribute("RepresentativesSections", IsGroup);

InstallMethod( RepresentativesSections,
        "for a group",
        [IsGroup],
        function(G)
    local   sections,  subs,  K,  N,  phi,  Q,  ccl,  x,  P,  section;

    sections:= [ ];
    subs:= List(ConjugacyClassesSubgroups(G), Representative);
    for K in subs do
        N:= Normalizer(G, K);
        phi:= NaturalHomomorphismByNormalSubgroup(N, K);
        Q:= Image(phi);
        ccl:= List(ConjugacyClassesSubgroups(Q), Representative);
        for x in ccl do
            P:= PreImages(phi, x);
            section:= Section(G, P, K);
            SetHomomorphismSection(section, phi);
            SetAsGroup(section, x);
            Add(sections, section);
        od;
    od;

    return sections;
end);


SectionsByType:= function(G)
    local   types,  recs,  secs,  sec,  type,  pos;
    types:= [];
    recs:= [];
    secs:= RepresentativesSections(G);
    for sec in secs do
        type:= IdSection(sec);
        pos:= Position(types, type);
        if pos = fail then
            Add(types, type);
            Add(recs, rec(type:= type, sections:= [sec]));
        else
            Add(recs[pos].sections, sec);
        fi;
    od;
    Sort(recs, function(a, b) return a.type < b.type; end);
    return recs;
end;


#############################################################################
##
##  KappaSections
##
KappaSections:= function(secL, secR)
    local   G;
    G:= secL!.G;
    if secR!.G <> G then
        Error("Parent groups are not identical");
    fi;

    return Index(G, Intersection(secL!.K, secR!.K))^-1;
end;


#############################################################################
##
##  SectionsByTop( G )
##  SectionsByBot( G )
##  SectionsByMid( G )
##
SectionsByTop:= function(G)
    local   secs,  ccs,  list,  sec;
    secs:= RepresentativesSections(G);
    ccs:= ConjugacyClassesSubgroups(G);
    list:= List(ccs, x-> []);
    for sec in secs do
        Add(list[PositionProperty(ccs, c-> TopSec(sec) in c)], sec);
    od;
    return list;
end;

SectionsByBot:= function(G)
    local   secs,  ccs,  list,  sec;
    secs:= RepresentativesSections(G);
    ccs:= ConjugacyClassesSubgroups(G);
    list:= List(ccs, x-> []);
    for sec in secs do
        Add(list[PositionProperty(ccs, c-> BotSec(sec) in c)], sec);
    od;
    return list;
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
