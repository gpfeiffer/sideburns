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
        "utility",
        "zgroups",
        "boucindex",
        "difunct",
        "latex",
        ] do
    ReadPackage(Concatenation("sideburns/lib/", name, ".g"));
od;
