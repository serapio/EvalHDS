#! /usr/bin/python
# this script converts the output of HierBayes into a segmented text
# example usage can be seen in runwphb.sh
# @author: lucien@discurs.us

import re, os, sys

hboutdir = '../wpeos_hb-em'
wpdir = '../wpeos_hb'
outdir = '../wpeos_hb-em_fix'
hboutdir, wpdir, outdir = sys.argv[1:]
print hboutdir, wpdir, outdir
if not os.path.exists(outdir): os.mkdir(outdir)

filelist = ['Albert_Einstein.txt.eos.tok']
filelist = os.listdir(hboutdir)

for fname in filelist:
    hbfile = open(os.path.join(hboutdir,fname)).readlines()
    hyp = [line for line in hbfile[-10:] if line.startswith('[[[')]
    if (len(hyp) != 1):
        print "Whoops! Segmentation failed on",fname
        print "hyp:",hyp
        continue
    hyp = hyp[0].strip()
    depth = len(re.match(r'^\[*',hyp).group(0))
    numstr = re.split(r'([0-9]+[ \]\[]*)',hyp)
    nummatch = [re.match(r'(?P<linenum>[0-9]+) ?(?P<ends>( \])*) ?(\[)*$',s) for s in numstr]
    boundaries = [(int(m.group("linenum")),depth - len(m.group("ends"))/2) for m in nummatch if m]
    print fname
    print depth,hyp
    #print boundaries

    if boundaries:
        boundaries = [(0,1)] + boundaries[:-1] + [(boundaries[-1][0],1)]
        boundaries.reverse()
        print boundaries

        txtfile = open(os.path.join(wpdir,fname)).readlines()
        txtfile = [line for line in txtfile if not re.search('===',line)]
        for linenum, depth in boundaries:
            txtfile = txtfile[0:linenum] + ['='*4 +' '+str(depth)+' '+ '='*4 +'\n'] + txtfile[linenum:]
            #print linenum
            #print txtfile[linenum-1:linenum+2]
        open(os.path.join(outdir,fname),'w').writelines(txtfile)
    else:
        print "No boundaries!",fname
