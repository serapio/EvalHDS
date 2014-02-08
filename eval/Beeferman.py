#! /usr/bin/python
## Beeferman.py
# @author: lucien@discurs.us

"""
Creates table of errors for hierarchical Beeferman segmentation metric

For usage see ./Beeferman.py --help
"""

import re, sys, os, random, math

from optparse import OptionParser

# global flags and constants
_DEBUG_ = False 	# Also set by option --debug
_END_FLAG_ = "wrap" 	# What to do at end of text. Also set by option --end=
initBoundScore = 1000	# Nominal rank of unlabeled possible boundaries

# boundaries are expected to look like this:  ==== 3 ====
segLineRegex = re.compile(' *=+ *(?P<rank>[0-9]*) *=+ *') 


def meanBeeferman(refData,hypData,useWindowDiff,precisionFactor):
    # compute P_k average over sequence of linearizations

    tmp = {}
    for line,rank in refData.items():
        tmp[rank] = tmp.get(rank,[]) + [(line,rank)]
    refList = sorted(tmp.items(),lambda x,y: cmp(x[0],y[0]))
    if _DEBUG_: print >>sys.stderr, 'ref',len(refList),refList

    # cut off the first two implied boundaries
    hypList = sorted(hypData.items(),lambda x,y: cmp(x[1],y[1]))
    if _DEBUG_: print >>sys.stderr, 'hyp:',len(hypList),hypList
    tmpRefList = sorted(refData.items(),lambda x,y: cmp(x[1],y[1]))
    if _DEBUG_: print >>sys.stderr, 'ref:',len(tmpRefList),tmpRefList
    maxRank = hypList[len(tmpRefList)-1][1]
    hypList = [item for item in hypList if item[1] <= maxRank]
    hypList = hypList[2:] #len(refList)]

    if _DEBUG_:
        print >>sys.stderr, 'maxRank',maxRank,'lists', refList, hypList
    # find degeneracy
    tmpDict = {}
    for line,rank in hypList:
        tmpDict[rank] = 0
    hypDegeneracy = len(hypList) - len(tmpDict)
    if _DEBUG_: print >>sys.stderr, 'hypList:',len(hypList),hypList, 'tmpDict:',len(tmpDict),tmpDict 

    # if there's some degeneracy in the ranks,
    # do some averages over trials with random jitter
    #hypTrials = (hypDegeneracy * 10 +1) 
    hypTrials = int(min(math.sqrt(hypDegeneracy), 30)/ precisionFactor /precisionFactor) + 1
    #hypTrials = (hypDegeneracy + 1)**2
    #hypTrials = int(math.sqrt(len(hypList)*(hypDegeneracy + 1)))
    if _DEBUG_: print >>sys.stderr, 'degeneracy:',hypDegeneracy,'trials:',hypTrials
    total = 0
    hypTmp = {}
    for trialj in range(hypTrials):
        for line in hypData:
            hypTmp[line] = hypData[line] + random.random()*0.1
        hypList = sorted(hypTmp.items(),lambda x,y: cmp(x[1],y[1]))

        total += mean_P_k(refList,hypList,useWindowDiff)
    return total/(hypTrials)

def mean_P_k(refList,hypList,useWindowDiff):
    # average P_k for sampled hyp boundaries

    refList = list(refList) # copy it so we're not changing the original

    # take off the two implied boundaries at beginning and end of text
    previous = refList.pop(0)[1]
    if len(previous)==1:
        previous += refList.pop(0)[1]

    #previous += refList.pop(0)[1]

    total = 0
    lenprev = len(previous)
    while refList:
        rankSet = refList.pop(0)[1]
        R_i = previous + rankSet
        H_i = hypList[:len(R_i)]
#        H_i = hypList 
        if False: #_DEBUG_:
            print >>sys.stderr, 'rankSet:',rankSet, 'boundaries:',len(R_i),R_i,len(H_i),H_i
        total += P_k(R_i,H_i,useWindowDiff) * len(rankSet)
        previous = R_i
    if _DEBUG_: print >>sys.stderr, 'len(refList)',len(refList), 'lenprev',lenprev, len(previous), previous
    return total/(len(previous)-lenprev)

def P_k(R_i,H_i,useWindowDiff):

    refNums = makeSegNumList(dict(R_i))
    K = computeK(refNums)
    refJs = makeJudgementList(refNums, K)
        
    hypNums = makeSegNumList(dict(H_i))
    hypJs = makeJudgementList(hypNums, K)
    if _DEBUG_:
        print >>sys.stderr, 'hypNums',len(hypNums),'hypJs',len(hypJs),'refNums',len(refNums),'refJs',len(refJs)
    
    errs = makeErrorList(hypJs, refJs, useWindowDiff)
    average = sum(errs)/float(len(errs))
    return average
    
def computeK( refSegNums ):
    """Computes the window width K based on the reference segmentation"""
    wordCnt = len(refSegNums)
    segCnt = refSegNums[-1]
    K = 0.5 * wordCnt / segCnt
    return int(K)

def collectBoundaryData( segFilename, doWordRate ):
    """Gather location and rank information for the boundaries in a file"""
    segFile = open( segFilename )
    segData = {}

    atomCnt = 0
    for line in segFile:
        line = line.strip()
        #if _DEBUG_: print >>sys.stderr, line
        m = segLineRegex.match(line)
        if m:
            if _DEBUG_: print >>sys.stderr, m.groups()
            if m.group('rank') == '':
                segData[atomCnt] = 1
            else:
                segData[atomCnt] = int(m.group('rank'))
        else:
            if doWordRate: # calculating word error rate
                for word in line.split():
                    if word != '' and word != '\x00':
                        atomCnt += 1
                        segData[atomCnt] = initBoundScore
            else: # calculating line error rate
                atomCnt += 1    
                segData[atomCnt] = initBoundScore
    
    # implied boundaries at beginning and end
    segData[atomCnt] = 0
    segData[0] = 0
    if _DEBUG_:
        print >>sys.stderr, 'total atoms:',atomCnt
        
    return segData
    
def makeSegNumList( segData ):
    """Compose list of which segment (1st, 2nd, etc) each word is in"""
    wordCnt = max(segData.keys())
    segNums = [0]*wordCnt
    seg = 0
    for idx in range(wordCnt):
        if idx in segData:
            seg += 1
        segNums[idx] = seg
    
    return segNums

    
def makeJudgementList( segNumList, K ):
    """Compose list of number of boundaries within K-word window"""
    # _END_FLAG_ = ('stop','lookback','wrap')[2]
    atomCnt = len(segNumList)
    maxSegNum = segNumList[-1]
    if _DEBUG_: # assert
        if max(segNumList) != maxSegNum:
            raise Exception, "last isn't max!"
    if _END_FLAG_ in ('wrap','lookback'):
        atomCnt = atomCnt - 1 # don't count the last window, since textend boundary doesn't count
    jlist = [0]*atomCnt
    for idx in range(atomCnt):
        if idx+K < atomCnt:
            jlist[idx] =  segNumList[idx+K] - segNumList[idx]
        else: # near the end of the text
            if _END_FLAG_ =='wrap':
                # wrap around tail and head
                #             ( segment boundaries in wrapped portion ) + (segment boundaries in trailing portion)
                jlist[idx] =  segNumList[(idx+K)%atomCnt] - segNumList[0] + maxSegNum - segNumList[idx] 
            elif _END_FLAG_ == 'stop':
                #jlist[idx] =  segNumList[idx+K] - segNumList[idx]
                #if _DEBUG_: print >>sys.stderr, "Not stopping yet",str(idx),str(atomCnt),str(K)
                jlist = jlist[:idx]
                break
            elif _END_FLAG_ =='lookback':
                # look back over top
                jlist[idx] =  segNumList[idx] - segNumList[(idx+K)%atomCnt] 
    if _DEBUG_: print >>sys.stderr, 'jlist atoms:', atomCnt, len(jlist),"_END_FLAG_",_END_FLAG_
    return jlist

def makeErrorList( refJudgements, hypJudgements, useWindowDiff ):
    """Compose binary list of when there is an error"""
    total = len(refJudgements)
    elist = [0]*total
    for idx in range( total ):
        refJ = refJudgements[idx]
        try: hypJ = hypJudgements[idx]
        except Exception, message:
            print idx, len(hypJudgements), len(refJudgements)
            raise Exception, message
        if useWindowDiff:
            #isErr = int(abs(refJ-hypJ))
            isErr = int(refJ!=hypJ)
        else:
            isErr = int((refJ and not hypJ) or (hypJ and not refJ))
        elist[idx] = isErr
    if _DEBUG_:
        print >>sys.stderr, 'P_k:',sum(elist),total,sum(elist)/float(total)
    return elist

def doAverages( errs, K ):
    """Average error rate over K-size bins"""
    aves = []
    while(len(errs) >= K):
        sub = errs[:K]
        errs = errs[K:]
        aves.append( sum(sub)/float(len(sub)) )
    return aves
                
        
def write( scores, showAverages ):
    """Print table to stdout"""

    print 'file ',
    for hyp in scores[0][1]:
        print '\t'+hyp['judge'][:6],
    if showAverages:
        print '\taverage'
    else:
        print

    if showAverages:
        averages = [{'score':0} for idx in range(len(scores[0][1]))]
        for filename,hypList in scores:
            for idx in range(len(hypList)):
                averages[idx]['score'] += hypList[idx]['score']
        for elem in averages:
            elem['score'] = elem['score'] / len(scores)
        scores.append(('average',averages))
    
    for filename,hypList in scores:
        print filename,
        if showAverages:
            average = sum([hyp['score'] for hyp in hypList])/len(hypList)
            hypList.append({'score':average})
        for hyp in hypList:
            print '\t%s' % round(hyp['score'],3),
        print


def parseArgs(args):
    """Parse command line arguments"""
    usage = """
    $ %prog [options] config1 [config2 [...]] > outfile
    
    where config* are configuration files in attribute value format

    ## config file for one text
    
        filename    task-identifier
        gold        path/to/goldstandard

        file        path/to/firstfile
        judge       producer-of-firstfile
        file        path/to/secondfile
        judge       ...
        ..."""
    parser = OptionParser(usage)
    parser.add_option( "-w", "--windowdiff", action="store_true",
                       default=False,
                       help="use Pevzner & Hearst's algorithm" )
    parser.add_option( "-a", "--averages", action="store_true",
                       default=False,
                       help="print judge and file averages" )
    parser.add_option( "-e", "--end", default="wrap",  
                       help="behavior within K atoms of text end: stop, lookback, or wrap [Default: %default]" )
    parser.add_option( "-l", "--line_rate", action="store_false",
                       default=True, dest="word_rate",
                       help="compute line (sentence) error rather than word error rate" )
    parser.add_option( "-H", "--no_hyp_levels", action="store_true",
                       default=False,
                       help="ignore level rankings of marked hypothesized boundaries" )
    parser.add_option( "-R", "--no_ref_levels", action="store_true",
                       default=False,
                       help="ignore level rankings of marked reference boundaries" )
    parser.add_option( "-p", "--precision", default=0.1, type="float",
                       help="precision factor [Default: %default]" )

    parser.add_option( "-D","--debug", action="store_true",
                       default=False,
                       help="turn on debugging output" )
    opts,args = parser.parse_args(args[1:])
    opts.files = args
    return opts

def parseConfig(filename):
    """Parse configuration file in attribute-value format:
        filename    task-identifier
        gold        path/to/goldstandard

        file        path/to/firstfile
        judge       producer-of-firstfile
        file        path/to/secondfile
        judge       ...
        ...
    """
    hyps = []
    for line in open(filename):
        nonComment = line.strip().split('#')[0]
        if nonComment != '':
            attr, val = nonComment.split()
            if attr == 'filename': filename = val
            if attr == 'gold': gold = val
            if attr == 'file':
                hyps.append({'file':val})
            if attr == 'judge':
                hyps[-1]['judge'] = val
    return filename, gold, hyps
                

def main(args=sys.argv):
    options = parseArgs(args)

    global _DEBUG_
    _DEBUG_ = options.debug
    global _END_FLAG_
    _END_FLAG_ = options.end

    scores = []
    for configfile in options.files:
        filename, ref, hyps = parseConfig( configfile )
        print >>sys.stderr, '\n'+filename
        refData = collectBoundaryData( ref, options.word_rate )
        for k,v in refData.items(): 
            if (v == initBoundScore): 
                del refData[k]
        if options.no_ref_levels:
            for k,v in refData.items(): refData[k] = min(1,int(v))

        if _DEBUG_: 
            print >>sys.stderr, 'ref:',refData
            boundaries = sorted(refData.keys())
            lengths = [boundaries[i+1]-boundaries[i] for i in range(len(boundaries)-1)]
            print >>sys.stderr, 'lengths:',' '.join([str(i) for i in lengths])

        for idx in range(len(hyps)):
            print >>sys.stderr, '.',
            hypData = collectBoundaryData( hyps[idx]['file'], options.word_rate )
            if options.no_hyp_levels:
                for k,v in hypData.items():
                    if (v < initBoundScore): hypData[k] = min(1,int(v))

            shortage = len(refData) - len(hypData)
            if shortage > 0:
                print >>sys.stderr, "Oops! This shouldn't happen anymore!"
                if _DEBUG_: 
                    print >>sys.stderr, len(hypData), shortage, max(refData.keys())
                randpts = random.sample(xrange(max(refData.keys())),shortage*5)
                for p in randpts:
                    if p not in hypData:
                        hypData[p] = initBoundScore
            if _DEBUG_:
                print >>sys.stderr, '\nfile:',hyps[idx]['file']
                print >>sys.stderr, 'hyp:',hypData
            hyps[idx]['score'] = meanBeeferman(refData,hypData,options.windowdiff,options.precision)
        scores.append((filename,hyps))

    
    write( scores,options.averages )

if __name__ == '__main__': 
    try:
        import psyco
        psyco.full()
    except:
        pass

    main()
