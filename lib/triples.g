#############################################################################
##
#W  triples.g                                                       sideBurns
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
##  Subgroups as Goursat Triples.
##
##  This file contains structures and functions for Groursat triples of a
##  direct product of finite groups.
##
##  A _Goursat triple_ of a direct product of groups $G \times H$
##  consists of a section of $G$ (the _source_ of the triple),
##  a section of $H$ (its _target_), and an isomorphism between the two
##  sections.
##
##  By a classical result known as Goursat's Lemma, the subgroups of
##  a direct product $G \times H$ are in bijection to its Goursat triples,
##  in a way that is made explicit by the programs in this file.
##


#############################################################################
##
##  Representation and Type
##
DeclareRepresentation("IsTriple",
        IsComponentObjectRep and IsAttributeStoringRep and IsPermGroup,
        ["source", "target", "map"]);

TripleType:=  NewType(CollectionsFamily(PermutationsFamily), IsTriple);


#############################################################################
##
#C  Triple( source, target, map ) . . . . . . . . . . . . . . .  constructor
##
##  a new triple object...
##
##  A Goursat triple is an object with three components:
##  two sections and an isomorphism between them.
##
Triple:= function(source, target, map)
    local   r;
    r:= rec(source:= source, target:= target, map:= map);
    return Objectify(TripleType, r);
end;

##  print method
InstallMethod(PrintObj, "for triples", true, [IsTriple], 0, function(gt)
    Print("Triple( ", gt!.source, ", ", gt!.target, ", ", gt!.map,  " )");
end);

InstallMethod(ViewObj, "for triples", true, [IsTriple], 0, function(gt)
    Print("<triple>");
end);

##############################################################################
##
##  Sections( triple )
##  Source( triple )
##  Target( triple )
##
##  Returns the sections of `triple`.
##
DeclareOperation("Sections", [IsTriple]);

InstallMethod(Sections, "for a triple", [IsTriple], function(gt)
      return [gt!.source, gt!.target];
end);

#DeclareOperation("Source", [IsTriple]);

InstallOtherMethod(Source, "for a triple", [IsTriple], function(gt)
      return gt!.source;
end);

DeclareOperation("Target", [IsTriple]);

InstallMethod(Target, "for a triple", [IsTriple], function(gt)
      return gt!.target;
end);

##############################################################################
##
##  Map( triple )
##  Mapping( triple )
##
##  `Map` returns the isomorphism of `triple`.
##  `Mapping` returns this isomorphism regarded as a general mapping
##  between the top groups of the sections of `triple`.
##
DeclareOperation("Map", [IsTriple]);

InstallMethod(Map, "for a triple", [IsTriple], function(gt)
    return gt!.map;
end);

DeclareOperation("Mapping", [IsTriple]);

InstallMethod(Mapping, "for a triple", [IsTriple], function(gt)
    return NarrowHomomorphism(gt!.source) * gt!.map
           * InverseHomomorphism(gt!.target);
end);


#############################################################################
##
##  Opposite( triple )
##
##  The _opposite_ of a triple $\theta: P_1/K_1 \to P_2/K_2$
##  is the inverse isomorphism $P_2/K_2 /to P_1/K_1$.
##  If `triple` corresponds to a subgroup of $G \times H$ then
##  its opposite corresponds to a subgroup of $H \times G$.
##
DeclareOperation("Opposite", [IsTriple]);

InstallMethod(Opposite, "for a triple", [IsTriple], function(gt)
    return Triple(gt!.target, gt!.source, InverseGeneralMapping(gt!.map));
end);


#############################################################################
##
#A  ProductGroup
##
DeclareAttribute("ProductGroup", IsTriple);

InstallMethod(ProductGroup, "for a triple", [IsTriple], function(gt)
    return DirectProduct(gt!.source!.G, gt!.target!.G);
end);


##############################################################################
##
##  triple ^ g
##
##  FIXME: this method is out of line as it relies on the GxH object ...
##
InstallOtherMethod(\^, "for triples", true, [IsTriple, IsMultiplicativeElementWithInverse], 0, function(gt, element)
    local   GH,  sec,  e,  s,  c,  i,  gens,  map,  imgs;

    GH:= ProductGroup(gt);
    sec:= Sections(gt);
    e:= [];  s:= []; c:= [];
    for i in [1,2] do
        e[i]:= element^Projection(GH, i);
        s[i]:= sec[i]^e[i];
        c[i]:= ConjugatorIsomorphism(sec[i]!.G, e[i]);
    od;

    # P1^c1/K1^c1 -> P1^c1 -> P1 -> P1/K1 -> P2/K2 -> P2 -> P2^c2 -> P2^c2/K2^c2
    gens:= GeneratorsOfGroup(AsGroup(s[1]));
    map:= InverseHomomorphism(s[1]) * c[1]^-1
          * NarrowHomomorphism(sec[1]) * Map(gt);
    imgs:= List(gens, x-> x^map);
    map:= InverseHomomorphism(sec[2]) * c[2]
          * NarrowHomomorphism(s[2]);
    imgs:= List(imgs, x-> x^map);
    map:= GroupHomomorphismByImages(AsGroup(s[1]), AsGroup(s[2]), gens, imgs);

    return Triple(s[1], s[2], map);
end);


##############################################################################
##
##  SubgroupTriple( triple )
##
##  Converts `triple` into a subgroup of $G \times H$.
##
DeclareAttribute("SubgroupTriple", IsTriple);

InstallMethod(SubgroupTriple, "for triples", [IsTriple], function(gt)
    local   sec,  GH,  P,  PP,  K,  KK,  p1,  p2,  hom1,  hom2,  prop;

    sec:= Sections(gt);
    GH:= DirectProduct(sec[1]!.G, sec[2]!.G);

    P:= List([1, 2], i-> Images(Embedding(GH, i), TopSec(sec[i])));
    PP:= ClosureGroup(P[1], P[2]);
    K:= List([1, 2], i-> Images(Embedding(GH, i), BotSec(sec[i])));
    KK:= ClosureGroup(K[1], K[2]);

    p1:= Projection(GH, 1);
    p2:= Projection(GH, 2);
    hom1:= NarrowHomomorphism(sec[1]);
    hom2:= NarrowHomomorphism(sec[2]);
    prop:= gh-> ((gh^p1)^hom1)^Map(gt) = (gh^p2)^hom2;

    return SubgroupProperty(PP, prop, KK);
end);


#############################################################################
##
##  tripleL = tripleR
##  tripleL < tripleR
##
##  GAP cannot compare homomorphisms.  So instead of comparing the components
##  we compare the corresponding subgroups.
##  Note that this does not take the ambient groups into account.
##
InstallMethod(\=, "for triples", true, [IsTriple, IsTriple], 0, function(gtL, gtR)
    return SubgroupTriple(gtL) = SubgroupTriple(gtR);
end);

InstallMethod(\<, "for triples", true, [IsTriple, IsTriple], 0, function(gtL, gtR)
    return SubgroupTriple(gtL) < SubgroupTriple(gtR);
end);


#############################################################################
##
##  IsDescendantTriple
##
##  test whether hi \geq' lo
##
IsDescendantTriple:= function(hi, lo)
  return IsDescendantSection(Sections(hi)[1], Sections(lo)[1])
         and IsDescendantSection(Sections(hi)[2], Sections(lo)[2])
         and IsSubgroup(hi, lo);
end;

#############################################################################
##
##  IsSubgroupSameU
##
IsSubgroupSameU:= function(hi, lo)
    return IsCompSubsectionSameU(Sections(hi)[1], Sections(lo)[1])
           and IsCompSubsectionSameU(Sections(hi)[2], Sections(lo)[2])
           and IsSubgroup(hi, lo);
end;

##############################################################################
##
#A  GeneratorsOfGroup( triple )
##
##  Converts `triple` into a list of generators of the corresponding subgroup
##  of $G \times H$.
##
##  The presence of this attribute makes `triple` behave like a group.
##
DeclareAttribute("GeneratorsOfGroup", IsTriple);

InstallMethod(GeneratorsOfGroup, "for triples", [IsTriple], function(gt)
    return GeneratorsOfGroup(SubgroupTriple(gt));
end);


##############################################################################
##
#F  IdentityTriple( sec )
##
##  Creates a triple with sections equal to `sec`
##  and an identity isomorphism between the quotient groups of the sections.
##
IdentityTriple:= function(sec)
    return Triple(sec, sec, IdentityMapping(AsGroup(sec)));
end;


##############################################################################
##
#F  TripleSubgroup( GxH, L )
##
##  Converts a subgroup $L$  of $G \times H$ into a triple.
##
TripleSubgroup:= function(GH, L)
    local   sections,  list,  i,  p,  e,  P,  K,  map;

    sections:= [];  list:= [];
    for i in [1, 2] do
        p:= Projection(GH, i);
        e:= Embedding(GH, i);
        P:= Image(p, L);
        K:= Image(p, Intersection(Image(e), L));
        sections[i]:= Section(Image(p), P, K);
        list[i]:= List(GeneratorsOfGroup(L), ab ->
                       (ab^p)^NarrowHomomorphism(sections[i]));
    od;
    map:= GroupHomomorphismByImages(AsGroup(sections[1]), AsGroup(sections[2]),
                    list[1], list[2]);

    return Triple(sections[1], sections[2], map);
end;


#############################################################################
##
#F  IdTriple( triple )
##
##  Returns the library number of the quotient group of the sections of
##  `triple`
##
IdTriple:= function(gt)
    return IdSection(Sections(gt)[1]);
end;


#############################################################################
##
##  TripleMorphism( mor )
##
TripleMorphism:= function(mor)
    local   u,  target;
    u:= mor!.group;
    target:= Section(u, u, TrivialSubgroup(u));
    return Triple(mor!.section, target, mor!.theta);
end;


#############################################################################
##
##  TripleMorphisms( morL, morR, d )
##
##  Produces a triple from two morphisms and a double coset rep
##
##  Needs to carefully engineer an isomorphism between the sections.
##
TripleMorphisms:= function(morL, morR, d)
    local   secL,  secR,  map,  gens;
    secL:= AsGroup(morL!.section);
    secR:= AsGroup(morR!.section);
    map:= morL!.theta * d * InverseGeneralMapping(morR!.theta);
    gens:= GeneratorsOfGroup(secL);
    map:= GroupHomomorphismByImages(secL, secR, gens, List(gens, x-> x^map));
    return Triple(morL!.section, morR!.section, map);
end;


#############################################################################
##
##  TriplesMorphisms( morL, morR )
##
##  Up to conjugacy, the triples of $G \times H$ which
##  up to automorphisms of $U$, can be decomposed as `morL` times the
##  opposite of `morR`.   This set corresponds to the conjugacy classes
##  of subgroups of $G \times H$ whose sections are conjugates of
##  those of `morL` and `morR` respectively.
##
# FIXME: use this function below and then move into morphisms.g
DoubleCosetsMorphisms:= function(morL, morR)
    local   U,  C1,  C2;

    U:= morL!.group;
    if morR!.group <> U  then return [];  fi;

    C1:= Conjugators(morL);
    C2:= Conjugators(morR);

    return DoubleCosets(AutomorphismGroup(U), C1, C2);
end;

# FIXME: use the above function below.
TriplesMorphisms:= function(morL, morR)
    local   U,  C1,  C2,  A,  ddd;

    U:= morL!.group;
    if morR!.group <> U  then return [];  fi;

    C1:= Conjugators(morL);
    C2:= Conjugators(morR);
    A:= AutomorphismGroup(U);
    ddd:= List(DoubleCosets(A, C1, C2), Representative);
    return List(ddd, d-> TripleMorphisms(morL, morR, d));
end;


#############################################################################
##
##  MorphismsTriple( triple )
##
##  how to split a triple into a pair of morphisms
##
MorphismsTriple:= function(gt)
    local   mor2,  mor1;
    mor2:= OneMorphismSection(Target(gt));
    mor1:= Morphism(gt!.source, mor2!.group, gt!.map * mor2!.theta);
    return [mor1, mor2];
end;


#############################################################################
##
##  ConjugacyClassesSubgroupsDirectProduct( G, H )
##
##  Computes the conjugacy classes of subgroups of the direct product
##  $G \times H$.
##
##  TriplesDirectProduct is the version that does not depend on an
##  object representing the direct product GxH.
##  It returns a list of records, one for each section type.
##  Each record contains a rectangle, with rows labelled by sections of
##  G, and columns by sections of H, and entries consisting of
##  conjugacy class representatives with those coordinates,
##  corresponding to double cosets of induced automorphism groups.
##
TriplesDirectProduct:= function(G, H)
    local   types,  secsG,  secsH,  rect,  monG,  row,  dbl,  monH;

    types:= [];
    for secsG in SectionsByType(G) do
        for secsH in SectionsByType(H) do
            if secsG.type = secsH.type then
                rect:= rec(type:= secsG.type, trip:= [], dblc:= []);
                rect.rows:= List(secsG.sections, OneMorphismSection);
                rect.cols:= List(secsH.sections, OneMorphismSection);
                Add(types, rect);
                for monG in rect.rows do
                    row:= [];  dbl:= [];
                    Add(rect.trip, row);
                    Add(rect.dblc, dbl);
                    for monH in rect.cols do
                        Add(row, TriplesMorphisms(monG, monH));
                        Add(dbl, DoubleCosetsMorphisms(monG, monH));
                    od;
                od;
            fi;
        od;
    od;
    return types;
end;


ConjugacyClassesSubgroupsDirectProduct:= function(G, H)
    local   GH,  triples,  reps;
    GH:= DirectProduct(G, H);
    triples:= List(TriplesDirectProduct(G, H), x-> x.trip);
    reps:= List(Flat(triples), SubgroupTriple);
    return List(reps, x-> ConjugacyClassSubgroups(GH, x));
end;

#############################################################################
##
##  ConstrictedTriple( triple, section )
##
##  Given a section that is a descendant of the source of `triple`,
##  constrict triple to sec.
##
##  FIXME:  this could be an operation RestrictedMapping ...
##
ConstrictedTriple:= function(gt, sec)
    local   mapping,  P,  K,  target,  inv,  nat,  gens,  imgs,  map;

    # check arguments for compatibility
    if not IsDescendantSection(gt!.source, sec) then
        Error("section must descend from source");
    fi;

    mapping:= Mapping(gt);
    P:= Image(mapping, TopSec(sec));
    K:= Image(mapping, BotSec(sec));
    target:= Section(gt!.target!.G, P, K);
    inv:= InverseHomomorphism(sec);
    nat:= NarrowHomomorphism(target);
    gens:= GeneratorsOfGroup(AsGroup(sec));
    imgs:= List(gens, x-> Representative(Images(mapping, Images(inv, x)))^nat);
    map:= GroupHomomorphismByImages(AsGroup(sec), AsGroup(target), gens, imgs);

    return Triple(sec, target, map);
end;


##############################################################################
##
#F  CompositionTriples( tripleL, tripleR )
##
##  computes the triple corresponding to the relation product of the
##  subgroups corresponding to the given triples.
##
CompositionTriples:= function(gtL, gtR)
    local   meet,  r,  l;
    meet:= MeetSections(gtL!.target, gtR!.source);
    r:= ConstrictedTriple(gtR, meet);
    l:= Opposite(ConstrictedTriple(Opposite(gtL), meet));
    return Triple(l!.source, r!.target, l!.map * r!.map);
end;


#############################################################################
##
#F  ButterflyIsomorphismSections( secL, secR )
##
##  By Zassenhaus' Butterfly Lemma, any two sections of a finite group $G$
##  determine a unique isomorphism between subsections of them.
##
Butterfly:= function(secL, secR)
    return CompositionTriples(IdentityTriple(secL), IdentityTriple(secR));
end;


##############################################################################
##
##  NormalizerTriple( triple )
##
##  Calculates the triple of the normalizer $N_{G \times H}(L)$, given the
##  triple of the subgroup $L$ of $G \times H$.
##
NormalizerTriple:= function(gt)
    local   aut,  U,  iso,  U1,  sec1,  U2,  sec2,  l,  r;
    aut:= List(MorphismsTriple(gt), Automizer);
    U:= ClosureGroup(aut[1]!.group, aut[2]!.group);
    iso:= IsomorphismPermGroup(U);
    U:= Image(iso);
    U1:= Image(iso, aut[1]!.group);
    sec1:= Section(U, U1, TrivialSubgroup(U1));
    U2:= Image(iso, aut[2]!.group);
    sec2:= Section(U, U2, TrivialSubgroup(U2));
    l:= Triple(aut[1]!.section, sec1, aut[1]!.theta * iso);
    r:= Triple(sec2, aut[2]!.section, InverseGeneralMapping(aut[2]!.theta * iso));
    return CompositionTriples(l, r);
end;


#############################################################################
##
##  Right regular Matrices for kappa-composition
##
KappaTriples:= function(gtL, gtR)
    return KappaSections(gtL!.target, gtR!.source);
end;

RightKappaStar:= function(gts, gtR)
    local   mat,  gtL,  comp;
    mat:= [];
    for gtL in gts do
        comp:= CompositionTriples(gtL, gtR);
        Add(mat, KappaTriples(gtL, gtR) * List(gts, x-> Iverson(x = comp)));
    od;
    return mat;
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
