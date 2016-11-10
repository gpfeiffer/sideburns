#############################################################################
##
#W  read.g                                                          sideBurns
##
#A  Götz Pfeiffer                                 goetz.pfeiffer@nuigalway.ie
##

#############################################################################
##
## read the actual code.
##
for name in [
        "sections",
        "posets",
        "sectlatt",
        "morphisms",
        "triples",
        "marks",
        "burnside",
        "latex",
        "utility",
        ] do
    ReadPackage(Concatenation("sideburns/lib/", name, ".g"));
od;
