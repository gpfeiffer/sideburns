#############################################################################
##
#W  latex.g                                                         sideBurns
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
##  LaTeX Helpers.
##
##  This file contains functions supporting LaTeX output of complex
##  structures.
##


#############################################################################
Ltx:= rec();

Ltx.display:= function(text)
    return Concatenation("\\[", text, "\\]\n");
end;

Ltx.command:= function(name, args)
    return Concatenation("\\", name, "{", args, "}");
end;

Ltx.open:= function(type, attrs)
    return Concatenation(Ltx.command("begin", type), attrs, "\n");
end;

Ltx.close:= function(type)
    return Concatenation(Ltx.command("end", type), "\n");
end;

Ltx.element:= function(type, attrs, content)
    return Concatenation(Ltx.open(type, attrs), content, Ltx.close(type));
end;

Tikz:= rec();

Tikz.picture:= function(attrs, text)
   return Ltx.element("tikzpicture", attrs, text);
end;


#############################################################################
TriplePositionsDirectProduct:= function(G, H)
    local   cc,  types,  secsG,  secsH,  rect,  monG,  row,  monH,
            tri;

    cc:= ConjugacyClassesSubgroupsDirectProduct(G, H);
    types:= [];
    for secsG in SectionsByType(G) do
        for secsH in SectionsByType(H) do
            if secsG.type = secsH.type then
                rect:= rec(type:= secsG.type, trip:= []);
                rect.rows:= Concatenation(List(secsG.sections, AllMorphismsSection));
                rect.cols:= Concatenation(List(secsH.sections, AllMorphismsSection));
                Add(types, rect);
                for monG in rect.rows do
                    row:= [];
                    Add(rect.trip, row);
                    for monH in rect.cols do
                        tri:= TriplesMorphisms(monG, monH)[1];
                        Add(row, PositionProperty(cc, c-> tri in c));
                    od;
                od;
            fi;
        od;
    od;
    return types;
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
