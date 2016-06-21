import os
import glob
from chimera import runCommand as rc
from chimera import replyobj

#change to folder with data files
#os.chdir("/Users/pett/data")

#em_map = "../mediator_modeling/em_map_files/asturias_mediator_translated.mrc"

cluster_dir = "./kmeans_weight_0_500_1/cluster.0/"
#cluster_dir = "../cluster.0/"

print cluster_dir

pdb_file = cluster_dir + "0.pdb"
view1_matrix = "./view1_matrix"
view2_matrix = "./view2_matrix"
background = "white"


colors = dict(
  #        colo       pdb code    mw
  beta_ANT_MG6beta  = ["#4c4c4c",  None,   64252],
  ANA  = ["#7f7fff",  None,   47718],
  MG7  = ["#ffdd40",  None,   43080],
  CUBf_CUBg  = ["#199999",  "H",    32205],
  TED  = ["#bfff3f",  None,   128796],
  MG8  = ["#8d381c",  "A",    32819],
  Anchor_C345C  = ["#ff7fff",  "I",    25585]
)


def open_mrc(dir, thresh=0.1):
   density_names = glob.glob(dir+"*.mrc")

   for model_num,d in enumerate(density_names):
      print "A"
      model_name = d.replace(".mrc","").replace("med","").replace(dir,"")
      print "B"
      print(d)
      print "C"
      replyobj.status("Processing " + d) # show what file we're working on
      print "D"
      rc("open " + d )
      print "E"
      rc("col " +colors[model_name][0] + " #" + str(model_num))
      print "F"
      if thresh == "byVol":
         volume = str(float(colors[model_name][2]) * 1.21)
         rc("vol #" + str(model_num) + " encloseVolume " + volume)
      else:
         rc("vol #" + str(model_num) + " level " + str(thresh))

   rc("focus")

def display_pdb(pdb, trans):
   rc("open #0 " + pdb)
   for prot, dat in colors.iteritems():
      model = prot.replace("Med","")
      if dat[1]:
        rc("transparency " + str(trans) + " #" + model)
        rc("col " + dat[0] + " #0:." + dat[1])
   try:
      rc("ribscale med #0")
   except:
      print("Ribbon style not found")

def show_em(em_file):
   rc("open #100 " + em_file)
   rc("vol #100 level 0.35 style mesh")


def view1_rotate(view1_matrix):
   rc("matrixset " + view1_matrix)
   rc("focus")

def view2_rotate(view2_matrix):
    rc("matrixset " + view2_matrix)
    rc("focus")

def save_img(file_name):
    rc("windowsize 1024 768")
    rc("focus")
    rc("copy file " + file_name + ".png" + " width 3072" + " height 2304" + " supersample 3" + " raytrace")

def set_background(background):
   rc("background solid " + background)

set_background(background)

open_mrc(cluster_dir)
#open_mrc(cluster_dir, "byVol")

#display_pdb(pdb_file, 70)

#view1_rotate(view1_matrix)
#file_name = ("kmeans_500_1_half1_view1")
#save_img(file_name)


#view2_rotate(view2_matrix)
#file_name = ("kmeans_500_1_half1_view2")
#save_img(file_name)

#show_em(em_map)

#save_img("out")

