# Owen Chapman
# 5/28/2021
# code for converting the cycles found in the _cycles.txt file into
# .bedpe (breakpoints) and
# .bed (sequences)

# TODO Only reads the first cycle.

import copy

class Breakpoint:
    def __init__(self,chr1,start1,chr2,start2,name='',score='',
                strand1='',strand2=''):
        self.chr1=chr1
        self.start1=int(start1)
        self.chr2=chr2
        self.start2=int(start2)
        self.name=name
        self.score=score
        self.strand1=strand1
        self.strand2=strand2
    
    def all(self):
        return [self.chr1,str(self.start1),str(self.start1),
                self.chr2,str(self.start2),str(self.start2),self.name,self.score,self.strand1, self.strand2]
    def __repr__(self):
        s = '\t'.join(self.all())+'\n'
        return s

class Segment:
    def __init__(self,chrom,start,end,name='',strand=''):
        self.chrom=chrom
        self.start=int(start)
        self.end=int(end)
        self.name=name
        self.strand=strand
    
    def all(self):
        return [self.chrom,str(self.start),str(self.end),self.name,self.score,self.strand]
    
    def __str__(self):
        return '\t'.join(self.all())+'\n'
    #def __copy__(self):
    #    return Segment(self.chrom,self.start.self.end,self.name,self.strand)
    
class Cycle:
    '''
    An ordered list of segments.
    '''
    @staticmethod
    def _merge_seqs(seqs):
        '''
        merge adjacent Segments.
        '''
        merged_seqs = []
        curr = None
        for i in range(len(seqs)):
            if curr == None:
                curr = copy.copy(seqs[i])
            next = seqs[(i+1)%len(seqs)]
            if curr.chrom == next.chrom and curr.end+1 == next.start and curr.strand == '+' and next.strand == '+':
                curr.end = next.end
            elif curr.chrom == next.chrom and curr.start == next.end+1 and curr.strand == '-' and next.strand == '-':
                curr.start = next.start
            else:
                merged_seqs.append(curr)
                curr=None
        if curr != None:
            merged_seqs[0]=curr
        return merged_seqs
    
    def __init__(self, filepath):
        ordered_seqs = []
        with open(filepath,'r') as f:
            seqs = []
            for line in f:
                if line.startswith('Segment'):
                    line = line.strip().split('\t')
                    seqs.append(Segment(line[2],line[3],line[4]))
                elif line.startswith('Cycle'):
                    line = line.split(';')[-1]
                    line = line.split('=')[-1]
                    line = line.split(',')
                    line = [(int(n[:-1]),n[-1]) for n in line] #[(47,'+')]
                    for (i,sign) in line:
                        s = copy.copy(seqs[i-1])
                        s.strand = sign
                        ordered_seqs.append(s)
					break
        ordered_seqs = Cycle._merge_seqs(ordered_seqs)
        self.segments = ordered_seqs
    
    def __str__(self):
        r=''
        for s in self.sequences:
            r+=str(s)
        return r
    
    def get_breakpoints(self):
        bps=[]
        for i in range(len(self.segments)):
            curr = self.segments[i]
            next = self.segments[(i+1)%len(self.segments)]
            pt1 = curr.end if curr.strand == '+' else curr.start
            pt2 = next.start if next.strand == '+' else next.end
            bp = Breakpoint(curr.chrom,pt1,next.chrom,pt2,strand1=curr.strand,strand2=next.strand)
            bps.append(bp)
        return bps

    def write_breakpoints(self,filepath):
        '''
        Write breakpoints as .bedpe
        '''
        with open(filepath,'w') as f:
            f.writelines(str(s) for s in self.get_breakpoints())
    
    def write_sequences(self,filepath):
        '''
        Write sequences as .bed
        '''
        with open(filepath,'w') as f:
            f.writelines(str(s) for s in self.sequences)
