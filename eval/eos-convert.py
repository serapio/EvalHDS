#! python
# This script converts between Eisenstein's format and the format used in the WP corpus
# @author: lucien@discurs.us



import re, os, textwrap

doHBayesFormat = True
doTextWrap = False

breakre = re.compile('(==+) ([^=]+) ==+')

def reformat(filename):
    text = [line.strip() for line in open(filename).readlines() if line.strip()]
    return ['==== 1 ====']+text+['==== 1 ====']
    
def old_reformat(filename):
    fstr = open(filename)
    outlines = ['==== 1 ====']
    for line in fstr:
        line = line.strip()
        m = breakre.search(line)
        if m:
            if not line.startswith('=='):
                prefix = line[:line.index('==')]
                line = line[line.index('=='):]
                outlines = outlines + [prefix]
            depth = len(m.group(1))
            line = ' '.join(['='*4,str(depth),'='*4])
            if outlines[-1].endswith('==') and line.endswith('=='):
                prevrank = int(outlines[-1].split()[-2])
                thisrank = int(line.split()[-2])
                if thisrank < prevrank: outlines[-1] = line
            else:
                outlines = outlines + [line]            
        elif line:
            outlines = outlines + [line]
    outlines = outlines + ['==== 1 ====']
    return outlines

def HierBayesFormat(text):
    breaks = [breakre.search(line) for line in text]
    breakdepths = [int(m.group(2)) for m in breaks if m != None]
    print sorted(dict(zip(breakdepths,[0]*len(breakdepths))).keys())
    maxdepth = max(breakdepths)
    for idx in range(len(text),0,-1):
        line = text[idx-1]
        m = breakre.search(text[idx-1])
        if m:
            depth = int(m.group(2))
            #breaks = [' '.join(['='*4,str(d),'='*4]) for d in range(depth,maxdepth+1)]
            breaks = [' '.join(['='*8,str(d)]) for d in range(depth,maxdepth+1)]
            text = text[:idx-1]+ breaks + text[idx:]
    return text
        
    

#def main():
if True:
    try:
        import psyco
        psyco.full()
    except:
        pass
    
    wpdir = "wpeos"
    textdir = "wpeos_hb"
    if not os.path.exists(textdir): os.mkdir(textdir)

    for fname in os.listdir(wpdir):
            print fname
            filename = os.path.join(wpdir, fname)
            text = reformat(filename)
            if doHBayesFormat:
                text = HierBayesFormat(text)
            if doTextWrap:
                wrapper = textwrap.TextWrapper(width=100,break_long_words=False)
                text = [wrapper.fill(line) for line in text]
            if text:
                outfname = os.path.join(textdir, fname)
                outfile = open(outfname,'w')
                print >>outfile, '\n'.join(text) #.encode('utf8')
                outfile.close()
