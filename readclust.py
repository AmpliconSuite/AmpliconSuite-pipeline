import logging

class dummy_read(object):
    def __init__(self,mrn,mstart,mir,mqn=""):
        self.reference_start = mstart
        self.reference_name = mrn
        self.is_reverse = mir
        self.query_name = mqn
        self.qstart,self.qend = -1,-1
        self.reference_end = mstart
        self.template_length = -1
        self.is_read1 = False
        self.is_read2 = True
        self.mapping_quality = 0

    def get_tags(self):
        return ["Dummy read",]

    def has_tag(self,_):
        return False


#TODO: how to check if two read clusts overlap? - same as other - but make sure ref chroms sorted
class pe_read_clust(object):
    def __init__(self,r1,r2,clustDelta=500):
        self.clust_ID = None
        self.left_reads, self.right_reads = [],[]
        self.size = 0
        self.centroid = (0,0)
        #self.is_foldback = is_foldback
        self.r_IDs = (r1.reference_name,r2.reference_name)
        self.add_pair_to_clust(r1,r2)
        self.total_diff = 0.0
        self.clustDelta = clustDelta

    def add_pair_to_clust(self,r1,r2):
        if (r1.reference_name,r2.reference_name) != self.r_IDs:
            logging.warning("Tried to add pair to cluster located on different chromosome(s)")

        else:
            self.left_reads.append(r1)
            self.right_reads.append(r2)
            self.size+=1
            if self.size > 1:
                fwd_diff = abs(r1.reference_end - self.centroid[0]) + abs(r2.reference_start - self.centroid[1])
                self.total_diff+=fwd_diff

            self.update_centroid()
            self.clust_ID = str(self.centroid)


    def update_centroid(self):
        currL,currR = self.centroid
        prevSS = len(self.left_reads)-1
        wL = prevSS*currL
        wR = prevSS*currR
        meanL = (wL + self.left_reads[-1].reference_end)/len(self.left_reads)
        meanR = (wR + self.right_reads[-1].reference_start)/len(self.right_reads)
        self.centroid = (meanL,meanR)


    # check if read pair (r1, r2) has overlap with centroid of this cluster
    def rp_has_overlap(self,r1,r2):
        tt = (r1.reference_name,r2.reference_name)
        if tt == self.r_IDs or tt == self.r_IDs[::-1]:
            fwd_match = abs(r1.reference_end - self.centroid[0]) < self.clustDelta and \
               abs(r2.reference_start - self.centroid[1]) < self.clustDelta

            rev_match =  abs(r2.reference_end - self.centroid[0]) < self.clustDelta and \
               abs(r1.reference_start - self.centroid[1]) < self.clustDelta

            return fwd_match or rev_match

        return False


    # check if cluster (tc) centroid overlaps with centroid of this cluster
    def clust_has_overlap(self,tc):
        if tc.r_IDs == self.r_IDs or tc.r_IDs == self.r_IDs[::-1]:
            fwd_match = abs(tc.centroid[0] - self.centroid[0]) < self.clustDelta and \
                   abs(tc.centroid[1] - self.centroid[1]) < self.clustDelta

            rev_match = abs(tc.centroid[1] - self.centroid[0]) < self.clustDelta and \
                   abs(tc.centroid[0] - self.centroid[1]) < self.clustDelta

            return fwd_match or rev_match

        return False


    def clust_to_bedpe(self):
        a = self.left_reads[-1]
        b = self.right_reads[-1]
        # return [str(a.reference_name),str(self.centroid[0]),str(b.reference_name),str(self.centroid[1]),
        #                   str(int(self.is_foldback)),str(self.size)]

        return [str(a.reference_name), str(self.centroid[0]), str(b.reference_name), str(self.centroid[1]),
                str(self.size)]

    def clust_to_string(self):
        s = self.clust_ID + " | #read_pairs: " + str(self.size) + "\n"
        s +=" ".join(["ReadNum","isReverse","qstart","qend","refname","ref_start","ref_end","template_length",
                      "qname","tags","\n"])
        for v in zip(self.left_reads,self.right_reads):
            for a in v:
                readno = ""
                if a.is_read1:
                    readno = "Read 1:  "
                elif a.is_read2:
                    readno = "Read 2:  "

                if a.is_reverse:
                    adir = "-"
                else:
                    adir = "+"

                sstart = readno + " " + adir + " "
                s+=sstart
                s+=" ".join([str(x) for x in [a.qstart,a.qend,a.reference_name,a.reference_start,a.reference_end,
                                             a.template_length,a.mapping_quality,a.query_name,a.get_tags()]])
                s+="\n"

            s+="\n"

        s+="\n"
        return s

