# Decay module: provenance and source patches

Task-1 (DESIGN.md section 2). Author program: **Light_scalar_decay**
(github.com/shiggs90/Light_scalar_decay), distributed as `Combined_Code.zip`,
shipped in this run as `inputs/decay_source/` (src, include, makefile, manual.pdf).
The program is the one used by Li et al. arXiv:2212.06186v1 (Ref. [LSDedcay] of that
paper) and by the FORESEE 2HDM_H model; see "Cross-check against the authors' own
tables" in `logs/decay_validation.md`.

The physical model is **not** touched. Only the two adaptations allowed by
DESIGN.md section 2 / prompt_2HDM.md section 3 are applied:

1. **PATCH-1**: define `M_PI` for the Windows/MSYS2 compiler in
   `include/Parameters.h`, guarded by `#ifndef` so that `-DM_PI` on the command
   line is also harmless. The author code uses `M_PI` in `CPevenDecay.cpp`
   (loop functions AQ0H/AQ1H/AQ12H) and `loop_func` in `ToolFunctions.cpp` uses
   `PI` from Parameters.h; MSVC does not predefine M_PI. No numerical constant
   is changed: 3.14159265358979323846 is the double-precision value.
2. **PATCH-2**: the original top-level `main()` (in the authors' `T1Aa.cpp`,
   which is **not** part of the shipped source tree) is replaced by the batch
   driver `src_patched/drv.cpp`. It initialises the model exactly as
   `manual.pdf` / prompt_2HDM.md section 3 prescribe and reads the public API
   of `Model.h`.

`src_patched/src/` is **byte-identical** to `inputs/decay_source/src/`
(`diff -r` empty); only `include/Parameters.h` differs, by the PATCH-1 hunk.
The patched tree also contains `drv.cpp` (new file) and a copy of the original
makefile (`makefile_orig`, unused).

## Initialisation sequence used by the driver

```cpp
BSModel model;
model.set_model(1);        // 2HDM Type-I
model.set_mH(mH);
model.set_mA(600.0);
model.set_mHC(600.0);
model.set_lambdav(0.0);
model.set_tanb(tanb);
model.set_cosba(1.0/tanb);
// read: DecayWidth_HTotal(), DecayWidth_H{gaga,ee,mumu,tautau,PiPi,KK,PiPiPiPi,cc,ss,Gluon}()
```

## Compile

```bash
bash work/decay/build.sh
# g++ (MSYS2 ucrt64, 15.2.0) -std=c++11 -O2 -DM_PI=3.14159265358979323846 \
#   -I src_patched/include -o build/light_scalar_decay.exe \
#   src_patched/drv.cpp src_patched/src/*.cpp src_patched/src/data_base/*.cpp
```

Compile log: `work/decay/build_compile.log` (clean, no warnings on stdout/stderr;
exit 0). Build time < 5 s.

## Hashes

| file | sha256 |
|---|---|
| `work/decay/build/light_scalar_decay.exe` | `5c63124c9b60ead8fd4c68874aa8a5ce916ddada3a2dd9d6ed877ec5ec1c909a` |
| `work/decay/build.sh` | `c0fc7f48ee8704d9b4d2284798f5883d89e41850af280dcd3ce7c24f0df1e04c` |
| `work/decay/src_patched/drv.cpp` | `c761be699d67cf04c557a9bb3b04b39ace8659582c3bbfc09a7b12633e3fcb77` |
| `work/decay/src_patched/include/Parameters.h` (patched) | `772af5bd28ec784e68400d15694c2b8bc463b94061b39da44e28556496a86395` |
| `work/decay/decay_widths.py` | `23f82ab8bf1a29462560bd2b6d1e7527b331dbcbab65e1c6596912ad09e5b908` |
| `work/decay/make_grids.py` | `ce3b2833d796a0eb0406582ce1c981c51fefd07486812bcbc73e21f049b0ce8b` |
| `work/decay/validate_decay.py` | `50b5bd2cb01b804b7afb14e1e39f59538b8466e937274a313eb2dcac4f2d4f29` |
| `work/decay/digitise_paper.py` | `408efebdfa12ae92b09e06ad35926c51007ed78b84bf0db036f072fc463cc52a` |
| `inputs/decay_source/include/Parameters.h` (original) | `5d1ffea04abedbeebadd45c1663231f0e321b3d89a74778795befa0686aaac61` |
| `inputs/decay_source/src/CPevenDecay.cpp` | `b1fd6986c792d11455c7c5fcb02d7f183f2d3b663b120b1d641b8796a3897d69` |
| `inputs/decay_source/src/Model.cpp` | `f493a8de23ccbb2e0a22475d9c3d3504c7ffb4198e74b17df171dea8899adb42` |
| `inputs/decay_source/src/ToolFunctions.cpp` | `9f9fe8990f357106d9d35cd3317876710ad3c70d9f7eef2d0b1b96ccd4b69835` |
| `inputs/decay_source/manual.pdf` | `09e2684b78af50065ca299b6bd44ed3586e5633e3bb057f263e88a1d01d5248e` |

The complete list of every input file of this run (including all 13
`src/data_base/*.cpp` interpolation tables) is `inputs/source.sha256`; the
decay-relevant subset is reproduced in `work/decay/validation.json`.

## Patch diff (only change to the author source)

```diff
--- inputs/decay_source/include/Parameters.h
+++ work/decay/src_patched/include/Parameters.h
@@ -3,6 +3,10 @@
 
 #include <cmath>
 
+#ifndef M_PI
+#define M_PI 3.14159265358979323846  /* PATCH-1 (decay module): define M_PI for the Windows/MSYS2 compiler. No physics change. */
+#endif
+
 #define PI 3.141592654
 
 //All masses are in units of GeV
```
## Other artefacts produced by this module (sha256 in `work/decay/grids_provenance.json`)

| file | content |
|---|---|
| `outputs/data/decay_widths_grid.csv` | 2680 rows (67 masses x 40 tanb), all Gamma_i, BR_i, GammaTot, ctau (canonical and 2.0e-16), mH=2.0 double-count flag |
| `work/common/scan_points_tanb40.txt` | 2680 lines "<mH> <tanb>" (%.17g, exact round-trip) |
| `work/common/{scan_points_tanb40,widths_img_13cols}_no2GeV.txt` | same but on `work/common/scan_mass_grid.txt` (66 masses, mH=2.0 dropped): 2640 rows — recommended for the scan |
| `work/common/widths_img_13cols.txt` | 2680 lines, 13 numbers per point, same order (header '#' lines document the columns) |
| `work/decay/scan_points_tanb40.txt`, `work/decay/widths_img_13cols.txt` | identical copies inside the owner's area |
| `work/decay/validation.json` | full machine-readable validation (sections basic, foresee, paper_figures, validation_table) |
| `outputs/figures/decay_validation.png` | overlay of our program vs the digitised paper figures |
| `work/decay/digitised_paper_points.csv` | every digitised point: figure, panel, channel, mH, digitised value, model value, log10 ratio |
