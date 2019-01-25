#Daniel Ence
#May 23, 2018

import argparse
import pysam
from scipy.stats import chisquare
import scipy.stats as stats
from numpy import mean
import pandas as pd


def main(args):
    current_GTer = PAV_GTer(args.cov_file, args.ref_position)

    print "Discovering PAVs"
    current_GTer.discover_PAVs(0,"")
    current_GTer.dump_PAVs("test_min_covg_0.bed","bed")
    current_GTer.dump_PAVs("test_min_covg_0.txt")
    current_GTer.dump_PAVs("test_min_covg_0.joinmap4.loc","joinmap")

    #print "Dumping covg table"
    #current_GTer.dump_covg_table(args.intervals, "test_covg_table.txt")


#    current_GTer.discover_PAVs(args.intervals,0)


#    print "Discovering PAVs"
#    current_GTer.discover_PAVs(args.intervals, 1)
#    current_GTer.dump_PAVs("test_mean_covg_1.bed", "bed")

#    print "Discovering PAVs"
#    current_GTer.discover_PAVs(args.intervals, 2)
#    current_GTer.dump_PAVs("test_mean_covg_2.bed", "bed")

#    print "Discovering PAVs"
#    current_GTer.discover_PAVs(args.intervals, 3)
#    current_GTer.dump_PAVs("test_mean_covg_3.bed", "bed")

class PAV_GTer(object):

    #just a dummy constructor for now
    def __init__(self):
        #list of bam file objects or just the bam file names?
        self.__cov_file = ""
        self.__ref_position
        self.__my_PAVs = {}

    def __init__(self,cov_file, reference_position):
        self.__cov_file = cov_file
        self.__ref_position=reference_position
        self.__my_PAVs = {}

    def __PAVs_properly_segregate(self, PAV_array):
        #print "testing for proper segregation"
        #print PAV_presence_absence_array
        #print "Chi-square test is:"
        #print "f_obs:
        #print [PAV_presence_absence_array.count(1), PAV_presence_absence_array.count(0)]
        #chi2_result = self.__PAVs_chi2_test(PAV_presence_absence_array)
        chi2_result = chisquare(f_obs=[PAV_array.count(0),PAV_array.count(1)])
        #print chi2_result
        #print chi2_result
        alpha = 0.05
        if(chi2_result[1] <= alpha):
#            print "less than or equal to alpha"
            return 0
        else:
#            print "greater than alpha"
            return 1

#    def __PAVs_chi2_test(self, PAV_array):

#        expected_table = pd.crosstab(index=pd.DataFrame([1]*(len(PAV_array)/2) + [0]*(len(PAV_array)/2))[0], columns="count")
#        observed_table = pd.crosstab(index=pd.DataFrame(PAV_array)[0], columns="count")

#        chi_squared_stat = (((observed_table-expected_table)**2)/expected_table).sum()
#        print chi_squared_stat

#        crit = stats.chi2.ppf(q=0.95,df=1)
#        print "critical value:" + str(crit)

#        p_value = 1 - stats.chi2.cdf(x=chi_squared_stat, df=1)
#        print "P value" + str(p_value)

#        return chisquare(f_obs=PAV_array)

    #print table of covgs
    def dump_covg_table(self, interval_file, filename):
        print "using this interval file:\t" +filename
        target_intervals = self.get_intervals(interval_file)

        table_file = open(filename,'w')
        table_file.write(self.__make_header_line() + "\n")

        for interval in target_intervals.keys():
            print "in this interval:\t" + interval

            curr_interval_covg_list = []

            for curr_bam in self.__bams:
                curr_bam_file = pysam.AlignmentFile(curr_bam,"rb")
                curr_interval_covg_list.append(
                    self.__mean_covg(curr_bam_file,interval.split("-")[0],target_intervals[interval]))
            curr_interval_covg_list.append(str(mean(curr_interval_covg_list)))
            curr_interval_covg_list.append(str(curr_interval_covg_list.count(0)))


            table_file.write("\t".join( str(x) for x in curr_interval_covg_list) + "\n")


    def discover_PAVs(self, max_read_count, intervals_file):

        if (len(self.__my_PAVs) > 0):
            self.__my_PAVs = {}
        #Default for now is assume that there is no intervals file
        if(intervals_file != None):
            #for every line in the cov_file
            coverage = open(self.__cov_file,'r')
            for line in coverage:
                curr_gene_PAV_list=[]
                curr_cov_line = line.strip().split("\t")
                curr_line_scaf=curr_cov_line[0]
                curr_line_start=curr_cov_line[1]
                curr_line_stop=curr_cov_line[2]
                curr_line_samples = curr_cov_line
                del curr_line_samples[0:3]

                del curr_line_samples[self.__ref_position]

                for sample in curr_line_samples:
                    if(int(sample) <= max_read_count):
                        curr_gene_PAV_list.append(0)
                    else:
                        curr_gene_PAV_list.append(1)
                #print self.__make_PAV_name(curr_line_scaf, (curr_line_start, curr_line_stop), "-") + "\t" + str(
                #    curr_gene_PAV_list)
                if(self.__PAVs_properly_segregate(curr_gene_PAV_list)):
                    self.__my_PAVs.setdefault(
                        self.__make_PAV_name(curr_line_scaf,(curr_line_start, curr_line_stop),"-"),
                        self.__make_PAV_entry(curr_gene_PAV_list))


        #else:
        #    print "in discover_PAVs method"
        #    print "using this file:\t" + interval_file

        #    target_intervals = self.get_intervals(interval_file)

        # make a nxm matrix where n=number of gams
        # m = number of genes (intervals in the #intervals" file
        #for interval in target_intervals.keys():
        #    #print "in this interval:\t" + interval
        #    curr_gene_PAV_list = []
        #    for curr_bam in self.__bams:
        #        curr_bam_file = pysam.AlignmentFile(curr_bam,"rb")

        #        counter = 0
        #        try:
        #            for read in curr_bam_file.fetch(interval.split("-")[0],target_intervals[interval][0], target_intervals[interval][1]):
        #                counter = counter + 1
        #        except ValueError:
        #            counter = 0

        #        if(counter > max_read_count):
        #            curr_gene_PAV_list.append(1)
        #        else:
        #            curr_gene_PAV_list.append(0)

            #print "this is the PAV list:"
            #print curr_gene_PAV_list
        #    if(self.__PAVs_properly_segregate(curr_gene_PAV_list)):
        #        self.__my_PAVs.setdefault(self.__make_PAV_name(interval,target_intervals[interval],"-"),
        #                                  self.__make_PAV_entry(curr_gene_PAV_list))
        #print "This many PAVs:"
        #print len(self.__my_PAVs)

    def __mean_covg(self, bam_file, seqID, interval_tuple):

        sum_covg = 0
        counter = 0
        if(self.__seqID_in_bam_file(bam_file,seqID) == 1):
            for pileupcolumn in bam_file.pileup(seqID,interval_tuple[0], interval_tuple[1]):
                counter = counter + 1
                sum_covg = pileupcolumn.n

            if(counter == 0):
                return 0
            else:
                return float(sum_covg) / float(counter)
        else:
            return 0

    def __seqID_in_bam_file(self,bam_file,seqID):

        match_found = 0

        for curr_ref in bam_file.references:
            if(curr_ref == seqID):
                match_found = 1
        return match_found

    def __make_PAV_entry(self, PAV_array):

        pvalue = chisquare(f_obs=[PAV_array.count(0),PAV_array.count(1)])
        entry_tuple = (PAV_array, pvalue)

        return entry_tuple

    def dump_PAVs(self,filename,format_string=None):

        print "This many PAVs:"
        print len(self.__my_PAVs)

        if(format_string == None):
            self.output_PAVs_table(filename)
        elif(format_string == "bed"):
            self.output_PAVs_bed(filename)
        elif(format_string == "joinmap"):
            #output formatted for joinmap
            self.output_PAVs_joinmap(filename)

    def output_PAVs_joinmap(self,filename):
        curr_file = open(filename,'w')

        #print the header for the joinmap file
        name_string = "name = HAP_PAV"
        curr_file.write(name_string + "\n")
        curr_file.write("popt = HAP\n")
        curr_file.write("nloc = " + str(len(self.__my_PAVs.keys())) + "\n")
        num_indvs = len(self.__my_PAVs[self.__my_PAVs.keys()[0]][0])
        curr_file.write("nind = " + str(num_indvs) + "\n")

        PAV_name_file = open("PAV_name_file.txt",'w')
        counter = 0
        for PAV in self.__my_PAVs.keys():
            locus_name = "PAV_" + str(counter)
            line = locus_name + " {1} (a,b)"
            for indv in self.__my_PAVs[PAV][0]:
                if(indv == 0):
                    line = line + " a "
                else:
                    line = line + " b "
            curr_file.write(line +"\n")
            PAV_name_file.write(PAV + "\t" + locus_name + "\n")
            counter = counter + 1



    def output_PAVs_table(self,filename):
        curr_file = open(filename,'w')
        for PAV in self.__my_PAVs.keys():
            PAV_name_parts = PAV.split("-")
            #print PAV_name_parts[0] + "\t" + str(PAV_name_parts[1]) + \
            #      "\t" + str(PAV_name_parts[2]) + "\t".join(str(x) for x in self.__my_PAVs[PAV][0]) + "\n"
            curr_file.write(PAV_name_parts[0] + "\t" + str(PAV_name_parts[1]) + \
                  "\t" + str(PAV_name_parts[2]) + "\t".join(str(x) for x in self.__my_PAVs[PAV][0]) + "\n")


    def output_PAVs_bed(self,filename):
        curr_file = open(filename,'w')
        for PAV in self.__my_PAVs.keys():
            PAV_name_parts = PAV.split("-")
            print PAV_name_parts[0] + "\t" + str(PAV_name_parts[1]) + "\t" + str(PAV_name_parts[2]) + "\t" + str(self.__my_PAVs[PAV]) + "\n"
            curr_file.write(PAV_name_parts[0] + "\t" + str(PAV_name_parts[1]) + "\t" + str(PAV_name_parts[2]) + "\n")


    def genotype_PAVs(self, bam_list, intervals):

        mean_covgs_dict = get_mean_covgs(bam_list,intervals)

        target_intervals = get_intervals(intervals)
        #for each gene in the intervals list
            #genotype each sample (bam)
                #is there no covg?
                    #then absent
                #is there covg == mean covg?
                    #then present
                #is there covg ~= 1/mean covg?
                    #then het


    def __get_mean_covgs(self, bam_list, intervals):

        #returns a list of the mean covg of each bam over all the defined intervals
        target_intervals = get_intervals(intervals)

        mean_covg_list = {}
        for curr_bam_filename in bam_list:
            curr_bam_covgs = {}
            curr_bam = pysam.AlignmentFile(curr_bam_filename, "rb")
            for gene in target_intervals.keys():
                curr_gene_curr_bam_covg = []
                curr_gene_pileup = curr_bam.pileup(gene, target_intervals[gene][0],target_intervals[1])
                for bp in curr_gene_pileup:
                    curr_gene_curr_bam_covg.append(bp.nsegments)
                curr_bam_covg[self.__make_PAV_name(gene,target_intervals[gene])] = numpy.mean(curr_gene_curr_bam_covg)
            mean_covg_list.setdefault(curr_bam_filename,numpy.mean(curr_bam_covgs))
        return mean_covg_list

    def __make_PAV_name(self, gene_name, interval, delim):
        name_delim = delim
        PAV_name = name_delim.join((gene_name, str(interval[0]), str(interval[1])))
        return PAV_name

    def get_intervals(self, intervals_filename):
        #assume bed file for now. Maybe need to accomodate for gff files too
        interval_dict = {}
        intervals_file = open(intervals_filename,'r')
        for line in intervals_file:
            line_parts = line.strip().split("\t")
            curr_interval_tuple = (int(line_parts[1]),int(line_parts[2]))
            curr_name = self.__make_PAV_name(line_parts[0],curr_interval_tuple,"-")
            interval_dict.setdefault(curr_name, curr_interval_tuple)

        return interval_dict

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cov_file",type=str,help="coverage file generated by bedtools \"multicov\"")
    parser.add_argument("--ref_position",type=int,help="position of the reference sample in the cov_file. zero-based. -1 means last")
    #add optional interval file to analyze a subset of the positions in the cov file
    args = parser.parse_args()
    main(args)

