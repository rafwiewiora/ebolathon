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
cd /Users/bb/repos/hackathons/ebolathon/runs/glycoprotein_t0r_100

# --- receptor ---
load receptor.pdb, receptor
hide everything, receptor
show cartoon, receptor
color grey70, receptor
set cartoon_transparency, 0.15, receptor

# --- docked ligand poses (best score first) ---
load rank1_-9.20.pdb, pose_rank1
show sticks, pose_rank1
util.cbag('pose_rank1')
label first pose_rank1, "1: -9.20"
load rank2_-8.70.pdb, pose_rank2
show sticks, pose_rank2
util.cbag('pose_rank2')
label first pose_rank2, "2: -8.70"
load rank3_-8.50.pdb, pose_rank3
show sticks, pose_rank3
util.cbag('pose_rank3')
label first pose_rank3, "3: -8.50"
load rank4_-8.40.pdb, pose_rank4
show sticks, pose_rank4
util.cbag('pose_rank4')
label first pose_rank4, "4: -8.40"
load rank5_-8.20.pdb, pose_rank5
show sticks, pose_rank5
util.cbag('pose_rank5')
label first pose_rank5, "5: -8.20"
load rank6_-8.10.pdb, pose_rank6
show sticks, pose_rank6
util.cbag('pose_rank6')
label first pose_rank6, "6: -8.10"
load rank7_-8.00.pdb, pose_rank7
show sticks, pose_rank7
util.cbag('pose_rank7')
label first pose_rank7, "7: -8.00"
load rank8_-8.00.pdb, pose_rank8
show sticks, pose_rank8
util.cbag('pose_rank8')
label first pose_rank8, "8: -8.00"
load rank9_-8.00.pdb, pose_rank9
show sticks, pose_rank9
util.cbag('pose_rank9')
label first pose_rank9, "9: -8.00"
load rank10_-7.90.pdb, pose_rank10
show sticks, pose_rank10
util.cbag('pose_rank10')
label first pose_rank10, "10: -7.90"

# --- pocket context + framing ---
group poses, pose_rank1 pose_rank2 pose_rank3 pose_rank4 pose_rank5 pose_rank6 pose_rank7 pose_rank8 pose_rank9 pose_rank10
select pocket_res, byres (receptor within 5 of (poses))
show lines, pocket_res
set label_size, 18
set label_color, black
zoom (poses), 6

set ray_shadows, 0
# end view.pml
