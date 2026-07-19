# ================================================================
# synthon-TS docked poses - live view
#
# LIVE REFRESH: while a screen is running it rewrites these files
# after every round. To pull in the newest poses, in the PyMOL
# command line just re-run this script:
#     PyMOL>  @view.pml
# (or use the menu: File > Reload All). 'reinitialize' below keeps
# re-loads clean - no duplicated/stacked objects.
# ================================================================
reinitialize
bg_color white
cd /Users/bb/repos/hackathons/ebolathon/runs/glycoprotein_batch_stress/shard_1

# --- receptor ---
load receptor.pdb, receptor
hide everything, receptor
show cartoon, receptor
color grey70, receptor
set cartoon_transparency, 0.15, receptor

# --- docked ligand poses (best score first) ---
pseudoatom box_center, pos=[-44.989, 14.937, -8.254]
zoom box_center, 15

set ray_shadows, 0
# end view.pml
