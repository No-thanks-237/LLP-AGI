
import json, os, hashlib, datetime

RUN = r"E:\AAAwork\LLP-AGI\runs\20261025_1325_UTC_2HDM_lightH_lhcb_downstream"
MG = os.path.join(RUN, "work", "magnus")
SAMP = os.path.join(MG, "samples")

def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()

val = json.load(open(os.path.join(MG, "validation.json")))
vint = json.load(open(os.path.join(MG, "validation_interface.json")))
inst_pik = json.load(open(os.path.join(MG, "jobs", "install_report_pik.json")))
inst_hl = json.load(open(os.path.join(MG, "jobs", "install_report_hadronlib.json")))
smk_sum = json.load(open(os.path.join(MG, "jobs", "install_report_pik.json")))["summary"][0]
norm_p = os.path.join(SAMP, "normalisation.json")
norm = json.load(open(norm_p))

pik_path = os.path.join(SAMP, "minbias_pik.txt.gz")
norm["inel"].update({
    "rows": val["pik"]["rows"],
    "sha256": sha256(pik_path),
    "source_job_id": "c588d44d0bd8ac7a",
    "source_bundle_sha256": "bfaf8806619f18fbdff3d73f84a3e4f9ffb8939c5def972b03130807216bddd6",
    "pythia_sigmaGen_mb": smk_sum["sigmaGen_mb"],
    "pythia_sigmaErr_mb": smk_sum["sigmaErr_mb"],
    "sigmaGen_vs_normalisation_note": "Pythia LO sigmaGen = 78.05 mb, the normalisation uses the measured 79.5 mb (DESIGN 4); the 1.8 % difference is recorded here, not corrected, because it affects a channel that is ~1e-4 of the signal.",
    "events_with_pi_or_K": smk_sum["events_with_pi_or_K"],
    "events_with_pi_or_K_abs_eta_lt6": smk_sum["events_with_pi_or_K_eta6"],
    "rows_abs_eta_lt6": smk_sum["rows_eta6"],
    "rows_per_event_mean": val["pik"]["rows_per_event_mean"],
    "charge_sum_piK": val["pik"]["charge_sum_piK"],
    "generation_seconds": 822,
    "format_note": "10 columns: event pdg px py pz E x y z t (x,y,z in mm, t in mm/c); pdg in {211,-211,321,-321,311,-311}; status 91/92 fragmentation hadrons, no eta cut; K0/K0bar are flavour states -> consumer applies the 50/50 K_L/K_S treatment of DESIGN 3.3",
})
norm["hadron_lib_info"] = {
    "dir": "../hadron_lib",
    "files": {k: {"rows": val["hadron_lib"][k]["rows"],
                  "sha256": inst_hl["outputs"]["hadron_lib_%s.txt.gz" % k]["sha256"],
                  "n_mass_nodes": val["hadron_lib"][k]["n_mass_nodes"]}
              for k in ("cc", "gg", "ss", "tautau")},
    "n_attempt_per_node": 1500,
    "source_job_id": "5b3b8b76029ea0af",
    "note": "generated on the exact per-channel subsets of work/analysis/grids/mass_grid.txt; all 71 open-channel (scan mass, channel) pairs have a node within the 0.6 % consumer tolerance (work/magnus/validation_interface.json); the consumer rotates each event by a random SO(3) matrix as in bc5_lib.HadronLib.",
    "channel_file_map": {"gg": "gg", "ss": "ss", "cc": "cc", "tau": "tautau"},
}
json.dump(norm, open(norm_p, "w"), indent=1)

manifest = {
    "run": "20261025_1325_UTC_2HDM_lightH_lhcb_downstream",
    "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "job_id": "5b3b8b76029ea0af",
    "bundle_sha256": "ff7cbd921a476ae1ee9fce3fa03ab94a56f6e9b07013c428a9bb8a5135d848fd",
    "generator_source_sha256": "ebae345ef7cc76a5776b87e37b27ee25114f8bfe4d5951e98c42bf865f94ad28",
    "generator": "Pythia 8.316, shadlib.cc verbatim from inputs/bc5_reuse",
    "seeds": {"gg": 911, "ss": 912, "cc": 913, "tautau": 914},
    "n_sample_per_node": 1500,
    "channels": {},
}
chans = {"gg": ("mass_list_gg.txt", 911), "ss": ("mass_list_ss.txt", 912),
         "cc": ("mass_list_cc.txt", 913), "tautau": ("mass_list_tau.txt", 914)}
for ch, (mf, seed) in chans.items():
    f = "hadron_lib_%s.txt.gz" % ch
    manifest["channels"][ch] = {
        "file": f, "sha256": inst_hl["outputs"][f]["sha256"], "bytes": inst_hl["outputs"][f]["bytes"],
        "rows": val["hadron_lib"][ch]["rows"],
        "mass_nodes": val["hadron_lib"][ch]["n_mass_nodes"],
        "mass_file": "work/analysis/grids/" + mf,
        "mass_list_sha256": sha256(os.path.join(RUN, "work", "common", mf)),
        "seed": seed,
        "empty_attempt_fraction": val["hadron_lib"][ch]["empty_attempt_fraction"],
        "nonzero_charge_groups": val["hadron_lib"][ch]["nonzero_charge_groups"],
    }
json.dump(manifest, open(os.path.join(MG, "hadron_lib", "manifest.json"), "w"), indent=1)
print(json.dumps(manifest, indent=1)[:2500])
print("normalisation updated")
