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
cd /Users/bb/repos/hackathons/ebolathon/runs/glycoprotein_t0r_trial

# --- receptor ---
load receptor.pdb, receptor
hide everything, receptor
show cartoon, receptor
color grey70, receptor
set cartoon_transparency, 0.15, receptor

# --- docked ligand poses (best score first) ---
load rank1_-7.70.pdb, pose_rank1
show sticks, pose_rank1
util.cbag('pose_rank1')
label first pose_rank1, "1: -7.70"
load rank2_-7.90.pdb, pose_rank2
show sticks, pose_rank2
util.cbag('pose_rank2')
label first pose_rank2, "2: -7.90"

# --- pocket context + framing ---
group poses, pose_rank1 pose_rank2
select pocket_res, byres (receptor within 5 of (poses))
show lines, pocket_res
set label_size, 18
set label_color, black
zoom (poses), 6

set ray_shadows, 0
# end view.pml
