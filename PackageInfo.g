#############################################################################
##
#W  PackageInfo.g                                                   sideBurns
##

#############################################################################
SetPackageInfo( rec(

PackageName := "sideBurns",
Subtitle := "Double Burnside Rings and Subgroups of Direct Products",
Version := "0.0.1",
Date := "26/03/2016",

Persons := [
  rec(
    LastName      := "Pfeiffer",
    FirstNames    := "Götz",
    IsAuthor      := true,
    IsMaintainer  := true,
    Email         := "goetz.pfeiffer@nuigalway.ie",
    WWWHome       := "http://schmidt.nuigalway.ie/goetz",
    PostalAddress := Concatenation( [
                     "NUI Galway\n",
                     "University Road\n",
                     "Galway, Ireland" ] ),
    Place         := "Galway",
    Institution   := "NUI Galway"
     ),
],

PackageWWWHome := "http://schmidt.nuigalway.ie/sideburns",

PackageDoc := rec(
  BookName := "sideburns",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0.html",
  PDFFile := "doc/manual.pdf",
  SixFile := "doc/manual.six",
  LongTitle := "Double Burnside Rings and Subgroups of Direct Products",
  Autoload := true
),

Dependencies := rec(
  GAP := ">=4.4",
  NeededOtherPackages := [ ["GAPDoc", ">= 1.1"] ],
  SuggestedOtherPackages := [],
  ExternalConditions := []
),

AvailabilityTest := ReturnTrue,
Autoload := false,
TestFile := "tst/testall.g",
Keywords := [
  "table of marks",
  "direct product",
  "Burnside ring",
  "double Burnside ring",
]
));
