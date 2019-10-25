#############################################################################
##
#W  PackageInfo.g                                                   sideBurns
##

#############################################################################
SetPackageInfo( rec(

PackageName := "sideBurns",
Subtitle := "Double Burnside Rings and Subgroups of Direct Products",
Version := "0.0.2",
Date := "25/10/2019",

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

Status:= "dev",
PackageWWWHome := "https://github.com/gpfeiffer/sideburns",
ArchiveURL := Concatenation( ~.PackageWWWHome, "archive/v", ~.Version),
ArchiveFormats:= ".tar.gz",
README_URL := Concatenation( ~.PackageWWWHome, "README.md" ),
PackageInfoURL := Concatenation( ~.PackageWWWHome, "PackageInfo.g" ),
AbstractHTML:= "<span class='pkgname'>sideBurns</span> provides functionality for computations with Double Burnside Rings.",

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
