# Fixed geometry approximation. Lengths in m, By in T.
# Full-envelope traversal is the explicit initial approximation, not a measured
# tracking efficiency. There is no UT requirement and no smearing or PID.
# RandomSeed must be set and recorded per run by the execution driver.
set ExecutionPath { Downstream TreeWriter }
module DownstreamPropagator Downstream {
  set InputArray Delphes/stableParticles
  set OutputArray tracks
  # Strict thresholds on the original birth four-vector, BEFORE the dipole.
  # Energy is E, not |p|. Physical event inputs must use HepMC GEV units.
  set MinBirthEnergy 5.0
  set MinBirthPT 0.5
  set By 1.0
  set FieldZMin 3.0
  set FieldZMax 7.0
  set SciFiZFront 7.7
  set SciFiZBack 9.4
  set SciFiHalfWidthX 3.24
  set SciFiHalfWidthY 2.415
  set BeamHoleRadius 0.09
  set AcceptanceMode full_envelope
}
module TreeWriter TreeWriter {
  add Branch Delphes/allParticles Particle GenParticle
  add Branch Downstream/tracks Track Track
}
# LLP 1 < decay z < 2.5 m, 2 < eta < 5, ancestry and event-level selections
# belong to the downstream event analysis. They are not track birth cuts here.
