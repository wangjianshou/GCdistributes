import sys
import collections
from operator import itemgetter
from itertools import dropwhile
import time
import pdb

chrName = {'1':'chr1', '2':'chr2', '3':'chr3', '4':'chr4', '5':'chr5', '6':'chr6',
           '7':'chr7', '8':'chr8', '9':'chr9','10':'chr10', '11':'chr11', '12':'chr12',
           '13':'chr13', '14':'chr14', '15':'chr15', '16':'chr16', '17':'chr17', '18':'chr18',
           '19':'chr19', '20':'chr20', '21':'chr21', '22':'chr22'}

snp = collections.defaultdict(list)
with open(sys.argv[1], "r") as f:
  for line in f:
    i, j = line.split()[0:2]
    j = int(j)
    snp[chrName[i]].append([j, 0, 0, j>100000 and j-100000 or 0, j+100000])
snp = dict(snp)
for i in snp: snp[i].sort(key=itemgetter(0))
sam = sys.argv[2] == '-' and sys.stdin or open(sys.argv[2], 'r')
sam = dropwhile(lambda x: x.startswith('@'), sam)
reads = [[i[2],int(i[3]),i[9]] for i in (line.split('\t') for line in sam)]
reads.sort(key=itemgetter(0,1))
#sam.close()
#t = time.time()
chr = None
for line in reads:
  #chr, pos, seq
  #pdb.set_trace()  #################
  if chr != line[0]:
    try:
      chrsnp = snp[line[0]]
    except KeyError:
      continue
    start, end = 0, len(chrsnp)
    chr, pos = line[0], int(line[1])
  else:
    pos = int(line[1])
    start = midL
    i = 1
    try:
      while (chrsnp[midR+i][0]-chrsnp[midR][0]) < (pos-lastPos):
        i += 1
      end = midR + i
    except IndexError:
      end = len(chrsnp) - 1
  while end >= start:
    #pdb.set_trace()  ##############
    mid = (end+start)>>1
    #print(str(mid))
    snpStart = chrsnp[mid][3]
    snpEnd = chrsnp[mid][4]
    #if max(pos, chrsnp[mid]-100000) >= min(pos+75, chrsnp[mid]+100000):
    if pos > snpEnd:
      start = mid + 1
      continue
    elif pos+75 < snpStart:
      end = mid - 1
      continue
    elif pos < snpStart:
      seq = line[2][snpStart-pos:]
    elif pos+75 > snpEnd:
      seq = line[2][:snpEnd-pos+1]
    else:
      seq = line[2]
    chrsnp[mid][1] += seq.count('G') + seq.count('C')
    chrsnp[mid][2] += len(seq)
    midL = mid - 1
    midR = mid + 1
    while midL >= 0 and max(pos, chrsnp[midL][3]) <= min(pos+75, chrsnp[midL][4]):
      #pdb.set_trace()  ####################
      snpStart = chrsnp[midL][3]
      snpEnd = chrsnp[midL][4]
      if pos < snpStart:
        seq = line[2][snpStart-pos:]
      elif pos+75 > snpEnd:
        seq = line[2][:snpEnd-pos+1]
      else:
        seq = line[2]
      chrsnp[midL][1] += seq.count('G') + seq.count('C')
      chrsnp[midL][2] += len(seq)
      midL -= 1
    midL = 0 if midL < 0 else midL
    while midR <= len(chrsnp)-1 and max(pos, chrsnp[midR][3]) <= min(pos+75, chrsnp[midR][4]):
      #pdb.set_trace()  ###############
      snpStart = chrsnp[midR][3]
      snpEnd = chrsnp[midR][4]
      if pos < snpStart:
        seq = line[2][snpStart-pos:]
      elif pos+75 > snpEnd:
        seq = line[2][:snpEnd-pos+1]
      else:
        seq = line[2]
      chrsnp[midR][1] += seq.count('G') + seq.count('C')
      chrsnp[midR][2] += len(seq)
      midR += 1
    lastPos = pos
    break

with open(sys.argv[3], 'w') as f:
  for i in sorted(chrName, key=lambda x: int(x)):
    chrsnpresult = '\n'.join(['\t'.join((i, '%d'%snppos[0], '%.4f'%(snppos[1]/snppos[2])))
                             for snppos in snp[chrName[i]] if snppos[2] != 0])
    if chrsnpresult.strip()=='':
      continue
    else:
      f.write(chrsnpresult+'\n')

















