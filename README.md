Microkinetic Pipeline — README

Last updated: 2025-10-17

What this pipeline does

It collects energetics from your structure folders, builds temperature-dependent Potential Energy Diagrams (PEDs), generates MKM input files with Eyring rates, and assembles combined mechanisms (Kinetics_Combined, Kinetics_All). Optional utilities plot the PED and prepare flux inputs.

# Recommended order (or run 'bash run.sh')
python extracted.py
python ped.py
python temp.py
python input.py
python input_rates.py
python input_combined.py
python all_input.py

# optional
python plot.py
python flux.py

The required folders must be present in the execution directory. Each folder corresponds to a reaction intermediate along the mechanistic pathways and collectively serves as the state-space input for the MKMCXXX package. To add, remove, or modify intermediates, update the Python script `x` (e.g., adjust the species list and pathway definitions) before running the pipeline.

1Vo
1Vo_down
1Vo_sub
2O
2Vo
2Vo_down
CO
CO_1Vo
CO_2O
CO_2Vo
CO_2Vo_up
CO_g
CO_O
CO2
CO2_1Vo
CO2_2Vo
CO2_g
CO2_O
O2
O2_1Vo
O2_2Vo
O2_g
OCOO
OCOO_1Vo
TS_COox_LH1
TS_COox_LH2
TS_COox_MvK1
TS_COox_MvK2
TS_O2diss_LH
TS_O2diss_MvK1
TS_O2diss_MvK2
TS_OCOOdiss_cLH1
TS_OCOOdiss_cMvK1
TS_Omig
temp (can be empty, it will complain if there is nno folder)

The TS folders must contain one transition state having a total.e file with the respective energy:
01
02
03
04_TS
05
06
07
08
09


SCII schema (data flow)
Species/TS (folders)   mech/ (defs)         input/ (placeholders)
         │               │                             │
         ▼               ▼                             │
                extracted.py  ─────────►  ped.py / temp_ped.py
                           │                     │
                           ▼                     ▼
                      data/*_data.txt        temp/ped_<T>.txt
                           │                     │
                           └──────────────┬──────┘
                                          ▼
                                       input.py
                                          ▼
                                    input_rates.py
                                          ▼
                        Kinetics_{Lh, L2, Mix, Mix2, Mvk}
                                 │               │
                                 │               └──► flux.py (optional)
                                 ▼
                           input_combined.py
                                 ▼
                          Kinetics_Combined
                                 ▼
                             all_input.py
                                 ▼
                            Kinetics_All

(Optionally, plot.py reads data/ped.txt to render the PED.)

Script-by-script (purpose & I/O)

1) extracted.py
Reads: mech/species.txt (lists of folders), then per folder the energetics (total.e, optional ZPE/S; TS energies from TS_*/*_TS/total.e).
Writes: data/*_data.txt for each mechanism (lh, l2, mix, mix2, mvk).
Notes: Missing pieces become “−”; those propagate into later steps if not fixed.

2) ped.py (T = 0 K) & temp_ped.py (variable T)
Reads: data/*_data.txt + mech/*_mech.txt.
Writes: data/ped.txt.
Notes: Each mechanism line defines which ZPE/S columns to include; T is 0 in ped.py, adjustable in temp_ped.py.

3) temp.py
Task: Temperature sweep (default: 0…1000 K, step 5 K).
Writes: temp/ped_<T>.txt.
Notes: Ensures temp/ exists; runs temp_ped.py per temperature.

4) input.py
Reads: temp/ped_0.txt (by default), mech/*_mech_kinetics.txt, input/input_placeholder_*.mkm.
Writes: input/{mech}_input.mkm (energies converted to kJ/mol).
Notes: Populates placeholders; does not touch your base input.mkm.

5) input_rates.py
Reads: input/{mech}_input.mkm.
Writes: Kinetics_{Lh|L2|Mix|Mix2|Mvk}/input.mkm with Eyring rate constants inserted.
Default T for Eyring: 423.15 K (editable in the script).
Notes: If placeholder tokens don’t match, rates won’t be replaced.

6) input_combined.py
Reads: all Kinetics_* (without NN branches).
Writes: Kinetics_Combined/input.mkm with unified &compounds, &reactions, and &settings.

7) all_input.py
Reads: all Kinetics_* + optionally nn/Kinetics_* (e.g., Mvknn).
Writes: Kinetics_All/input.mkm.

8) plot.py (optional)
Reads: data/ped.txt (or per-T files).
Writes: a PED figure (PNG) with your color/line-style rules.

9) flux.py (optional)
Task: Produces flux/ inputs from the Kinetics outputs (see script for exact expectations).

Minimal sanity checklist
All species/TS folders exist and contain the expected energetics (total.e; TS in *_TS/total.e).
mech/*_mech.txt and mech/*_mech_kinetics.txt are present and use consistent names/indices.
input/input_placeholder_*.mkm placeholders match the tokens your scripts replace.
After running extracted.py, verify data/*_data.txt do not contain only “−”.

Frequently adjusted parameters
PED temperature: ped.py (T = ...) or the sweep in temp.py (T_start/T_end/T_step).
Eyring temperature: input_rates.py (default T = 423.15 K).
Runs/pressures & grids: in input_combined.py / all_input.py (&settings, &runs blocks).
Units: PED in eV, kinetics in kJ/mol (conversion handled in input.py).

Troubleshooting (short & practical)
“Mix is missing”
mech/species.txt: is the mixed_folders block filled and pointing to real folders?
data/mix_data.txt: does it contain numbers (not just “−”)? If not, folders are empty/missing total.e.
mech/mix_mech.txt & mech/mix_mech_kinetics.txt: correct names/indices?
input/input_placeholder_mix.mkm: placeholder names match what input_rates.py expects?
TS not picked up → does the TS energy live in TS_*/*_TS/total.e?
Rates not replaced → placeholder tokens in the placeholder and the code must match exactly.
Empty temp-PEDs → temp.py ran fine, but if data/*_data.txt had “−”, the PED at those T will be incomplete.
Units confusion → Energies are in eV on the PED side; converted to kJ/mol in input.py for kinetics.

