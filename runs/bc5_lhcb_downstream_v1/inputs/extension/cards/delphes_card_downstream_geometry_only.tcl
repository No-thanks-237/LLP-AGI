# Software geometry tests only. Not the production kinematic selection.
# Keep identical geometry by loading the formal card and explicitly disabling
# only the two birth thresholds; this file is copied alongside it in the package.
# Delphes overrides Tcl source to resolve relative to the top card directory.
source delphes_card_downstream.tcl
set Downstream::MinBirthEnergy 0.0
set Downstream::MinBirthPT 0.0
