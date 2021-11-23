#!/usr/bin/env python

import argparse
from pyfaidx import Fasta
from pybedtools import BedTool, Interval
import sys

##########################
# FUNCTION
def gtf_to_exon_bed(file, chrom_size):
    """
    @summary :  store +/- 5bp of exon bounderies in bed object
    @param file : [str] path to gtf to parse
    @param chrom_size : [dict] keys = chrom name : value = length
    return BedTool object of exon boundaries (+/- 5bp)
    """
    intervals = list()
    gff = BedTool(file)

    for feature in gff:
        if feature[2] == 'exon':
            # if exon start very close to the begining of the chromosome
            start = feature.start-4 if feature.start > 4 else 0
            intervals.append(Interval(feature.chrom, start, feature.start+5))
            # if exon stop very close to the end of the chromosome
            stop = feature.stop+4 if feature.stop < chrom_size[feature.chrom]-4 else chrom_size[feature.chrom]
            intervals.append(Interval(feature.chrom, feature.stop-5, stop))
            # case of exon shorter than 6 pb (so start very close to the end of chromosome in vice versa) is not taken into account
            
    return BedTool(intervals)



def process(params):
    """
    @summary : Annot variant, by adding INFO fields :
         - LOCALISATION: PASS or a combination of jnct5nt (near exon junction), homodimer5 ( in homopolymer), cl3snp35nt(in SNP cluster)
         - GTpopCR= X% : percentage of individuals genotyped
            - DPgt05rPopCR= X% : percentage of individuals genotyped with DP>5
    """

    # load indexed fasta genome
    genome = Fasta(params.in_fasta)
    chrom_size = { chr:len(genome[chr]) for chr in genome.keys()}

    # store +/- 5bp of exon bounderies in bed object
    exon_bed = gtf_to_exon_bed(params.in_gtf,chrom_size)

    # input vcf to annotate
    FH_in = open(params.in_vcf, "rt")

    # annotated VCF 
    FH_out = open(params.out_vcf, "wt")
    FH_out_gt = open(params.out_vcf_gt, "wt")

    # GT DP AD tables of filtered variant (GTpopCR >50 && DPgt05rPopCR > 20)    
    FH_gt = open(params.out_gt, "wt")
    FH_dp = open(params.out_dp, "wt")
    FH_ad = open(params.out_ad, "wt")

    first = True
    samples = list()
    nb_samples = 0

    # summary:
    tot_snp = 0 
    FS = 0
    QD = 0
    not_SNP_BiAll = 0
    FILTER_PASS = 0
    cl3snp35nt = 0
    jnct5nt = 0
    homodimer5 = 0
    GTpopCR50 = 0
    DPgt05rPopCR20 = 0
    GTpopCR50_DPgt05rPopCR20 = 0
    LOC_PASS = 0

    # FILTER PASS, LOC PASS, GTpopCR > 50, DPgt05rPopCR > 20
    PERFECT = 0

    # compteur d'ecriture par bloc
    write_variant = list()
    write_variant_gt = list()
    write_table_dp = list()
    write_table_ad = list()
    write_table_gt = list()
    # parse input file
    for line in FH_in:
        if line.startswith("#"):
            FH_out.write(line)
            FH_out_gt.write(line)

            # add INFO annotation in header
            if first:
                FH_out.write("##INFO=<ID=LOCALISATION,Number=.,Type=String,Description=\"is the variant present in a cluster (3 SNP in a window of 35nt) : cl3snp35nt, at 5bp of an exon junction : jnct5nt, in a 5bp homopolymer : homodimer5 or not : PASS \">\n")
                FH_out.write("##INFO=<ID=GTpopCR,Number=A,Type=Float,Description=\"percentage of individuals genotyped \">\n")
                FH_out.write("##INFO=<ID=DPgt05rPopCR,Number=A,Type=Float,Description=\"percentage of individuals genotyped with a DP > 5 \">\n")
                first = False
                FH_out_gt.write("##INFO=<ID=LOCALISATION,Number=.,Type=String,Description=\"is the variant present in a cluster (3 SNP in a window of 35nt) : cl3snp35nt, at 5bp of an exon junction : jnct5nt, in a 5bp homopolymer : homodimer5 or not : PASS \">\n")
                FH_out_gt.write("##INFO=<ID=GTpopCR,Number=A,Type=Float,Description=\"percentage of individuals genotyped \">\n")
                FH_out_gt.write("##INFO=<ID=DPgt05rPopCR,Number=A,Type=Float,Description=\"percentage of individuals genotyped with a DP > 5 \">\n")
                first = False

            # save total number of samples
            if line.startswith("#CHROM"):
                samples = line.split()[9:]
                nb_samples = len(samples)
                FH_gt.write("\t".join(["#chrom", "pos", "REF", "ALT", "DPpop", "jnct5nt", "homodimer5" ,"cl3snp35nt" , "LOCALISATION"]) + "\t" + 
                    "\t".join([s + "_GT" for s in samples]) + "\n")
                FH_dp.write("\t".join(["#chrom", "pos", "REF", "ALT", "DPpop", "jnct5nt", "homodimer5" ,"cl3snp35nt" , "LOCALISATION"]) + "\t" + 
                    "\t".join([s + "_DP" for s in samples]) + "\n")
                FH_ad.write("\t".join(["#chrom", "pos", "REF", "ALT", "DPpop", "jnct5nt", "homodimer5" ,"cl3snp35nt" , "LOCALISATION"]) + "\t" + 
                    "\t".join([s + "_AD_" + a for s in samples for a in ["REF","ALT"] ]) + "\n")
                #~ print ("end parsing header\n")
            continue

        tot_snp += 1
        line_list = line.strip().split("\t")
        chrom = line_list[0]
        pos = int(line_list[1])
        ref = line_list[3]
        alt = line_list[4].split(",")
        filter = line_list[6].split(";")
        info= line_list[7].split(";")
        variant_bed = BedTool(chrom + " " + str(pos - 1) + " " + str(pos), from_string=True)

        ### Update and count Filter criteria
        
        ## check SNP bi-allelic
        if len(ref) != 1 or len(alt) != 1 or len(alt[0]) != 1:
            not_SNP_BiAll += 1
            filter.append("not_SNP_BiAll")
            if "PASS" in filter:
                filter.remove("PASS")
        
        ## check SNP not filtered on FS
        if "FS" in filter:
            FS += 1
        
        ## check SNP not filtered on QD
        if "QD" in filter :
            QD += 1

        ## check localisation annotation (move it on INFO, and count only if FILTER PASS)
        LOCALISATION = list()
        # add cluster localisation annotation
        if "SnpCluster" in filter:
            LOCALISATION.append("cl3snp35nt")
            filter.remove("SnpCluster")
            if len(filter) == 0:
                filter.append("PASS")
            if "PASS" in filter:
                cl3snp35nt += 1
        
        ## count kept SNP
        if "PASS" in filter:
            FILTER_PASS += 1
        else: # ignore call rate computing, writing SNP and GT DP AD tables
            continue
        
        # add junction localisation annotation
        if exon_bed.intersect(variant_bed) :
            LOCALISATION.append("jnct5nt")
            jnct5nt += 1
        
        # add homopolymer localisation annotation
        start = pos - 5 if pos > 5 else 0
        stop = pos + 4 if pos < len(genome[chrom])-4 else len(genome[chrom])
        env_seq = genome[chrom][start : stop ].seq.upper()
            
        if 'A'*5 in env_seq or 'C'*5 in env_seq or 'G'*5 in env_seq or 'T'*5 in env_seq:
            LOCALISATION.append("homodimer5")
            homodimer5 += 1

        # add unambiguous localisation
        if len(LOCALISATION) == 0 :
            LOCALISATION.append("PASS")
            LOC_PASS += 1

        # update INFO with LOCALISATION
        info.append("LOCALISATION="+",".join(LOCALISATION))

        GT_idx = line_list[8].split(":").index("GT")
        DP_idx = line_list[8].split(":").index("DP")
        AD_idx = line_list[8].split(":").index("AD")

        nb_genotype = 0
        nb_genotype_dp5 = 0
        gt_list = list()
        dp_list = list()
        ad_list = list()
        
        # parse genotype to compute call rates
        for idx,call in enumerate(line_list[9:]) :
            call_list = call.split(":")
            genotype = call_list[GT_idx]
            gt_list.append(genotype)
            dp = int(call_list[DP_idx]) if DP_idx < len(call_list) and call_list[DP_idx] != "." else 0
            dp_list.append(str(dp))
            ad = call_list[AD_idx] if AD_idx < len(call_list) else "0,0"
            ad_list.append(ad)
            if genotype != "./." :
                nb_genotype += 1
                if dp > 5:
                    nb_genotype_dp5 += 1

        GTpopCR = round(nb_genotype * 100.0 / nb_samples, 2)
        DPgt05rPopCR = round(nb_genotype_dp5 * 100.0 / nb_samples, 2)

        # update INFO with call rates
        info.append("GTpopCR="+str(GTpopCR))
        info.append("DPgt05rPopCR="+str(DPgt05rPopCR))

        if GTpopCR >= 50 :
            GTpopCR50 += 1

        if DPgt05rPopCR >= 20 :
            DPgt05rPopCR20 += 1

        if "PASS" in filter and "PASS" in LOCALISATION and GTpopCR > 50 and DPgt05rPopCR > 20 :
            PERFECT += 1

        # update line_list
        line_list[6] = ";".join(filter)
        line_list[7] = ";".join(info)
        # ECRIRE LES VARIANTS FILTER PASS, GTpopCR50 & DPgt05rPopCR20 DANS LES TABLES GT DP AD
        if "PASS" in filter:
            write_variant.append("\t".join(line_list) + '\n')
            if len(write_variant) == 10000:
                print ("writing 10000 filtered variants")
                FH_out.write("".join(write_variant))
                write_variant=list()
            if GTpopCR >= 50 and DPgt05rPopCR >= 20 :
                write_variant_gt.append("\t".join(line_list) + '\n')
                if len(write_variant_gt) == 10000:
                    print ("writing 10000 filtered variants on GTCR50 & DPgt5CR20")
                    FH_out_gt.write("".join(write_variant_gt))
                    write_variant_gt=list()
                
                GTpopCR50_DPgt05rPopCR20 += 1
                # recover DP
                dp_all = ""
                for i in info:
                    if i.startswith("DP"):
                        dp_all = i.split("=")[1]
                        break

                common_str=[chrom, str(pos), ref, alt[0], dp_all]
                common_str.append("1") if  "jnct5nt" in LOCALISATION else common_str.append("0")
                common_str.append("1") if  "homodimer5" in LOCALISATION else common_str.append("0")
                common_str.append("1") if  "cl3snp35nt" in LOCALISATION else common_str.append("0")
                common_str.append("PASS") if  "PASS" in LOCALISATION else common_str.append(".")
                
                write_table_gt.append("\t".join(common_str) + "\t" + "\t".join(gt_list) + "\n")
                write_table_dp.append("\t".join(common_str) + "\t" + "\t".join(dp_list) + "\n")
                write_table_ad.append("\t".join(common_str) + "\t" + "\t".join([ad.replace(",","\t") for ad in ad_list]) + "\n")
                
                if len(write_table_gt) == 10000:
                    FH_gt.write("".join(write_table_gt))
                    FH_dp.write("".join(write_table_dp))
                    FH_ad.write("".join(write_table_ad))
                    write_table_dp = list()
                    write_table_ad = list()
                    write_table_gt = list()

    if len(write_variant) > 0:
        FH_out.write("".join(write_variant))
    FH_out.close()
    if len(write_variant_gt) > 0:
        FH_out_gt.write("".join(write_variant_gt))
    FH_out_gt.close()
        
    if len(write_table_gt) > 0:
        FH_gt.write("".join(write_table_gt))
        FH_dp.write("".join(write_table_dp))
        FH_ad.write("".join(write_table_ad))
    FH_gt.close()
    FH_dp.close()
    FH_ad.close()
    
    # ECRIRE SUMMARY
    FILTER_PASS_percent = round(FILTER_PASS * 100.0 / tot_snp , 2) if tot_snp > 0 else 0.00
    GTpopCR50_DPgt05rPopCR20_percent = round(GTpopCR50_DPgt05rPopCR20 * 100.0/ tot_snp,2) if tot_snp > 0 else 0.00
    PERFECT_percent = round(PERFECT * 100.0 / tot_snp , 2) if tot_snp > 0 else 0.00
    FH_summary = open(params.summary, "wt")
    FH_summary.write("# Input vcf file is " + params.in_vcf + "\n" \
    + "\t contains : " + str(tot_snp) + " variants\n" \
    + "# Filtering step\n" \
    + "\t filtered on non biallelic SNP : " + str(not_SNP_BiAll) + "\n" \
    + "\t filtered on FS : " + str(FS) + "\n" \
    + "\t filtered on QD : " + str(QD) + "\n" \
    + "\t kept SNP : " + str(FILTER_PASS) + " (" + str(FILTER_PASS_percent) + "%)\n" \
    + "\t\t those SNP are stored in " + params.out_vcf + "\n" \
    + "\t# filter PASS SNP annotation step\n" \
    +"\t\t localized in SNP cluster : " + str(cl3snp35nt) + "\n" \
    +"\t\t localized near exon junctions : " + str(jnct5nt) + "\n" \
    +"\t\t localized in homopolymers : " + str(homodimer5) + "\n" \
    +"\t\t localized in unbiased region : " + str(LOC_PASS) + "\n" \
    + "\t# computing call rate step\n" \
    + "\t\t number of SNP with a call rate >= 50% : " + str(GTpopCR50) + "\n" \
    + "\t\t number of SNP with a signigicant (GT with DP>5) call rate >= 20% : " + str(DPgt05rPopCR20) + "\n" \
    + "\t\t number of SNP with GTpopCR >= 50% and DPgt05rPopCR >= 20 : " + str(GTpopCR50_DPgt05rPopCR20) + " ("+ str(GTpopCR50_DPgt05rPopCR20_percent) + "%)\n" \
    + "\t\t\t those SNP are described in " + params.out_vcf_gt + ", " + params.out_gt + ", "+ params.out_dp + ", " + params.out_ad + "\n" \
    + "\t\t number of SNP with FILTER PASS, LOCALISATION PASS, GTpopCR >=50%, DPgt05rPopCR >= 20% : " + str(PERFECT) + " (" + str(PERFECT_percent) + "%)\n")
    

##########################
# MAIN

if __name__ == '__main__':
    # Parameters
    parser = argparse.ArgumentParser(description="Annot variant, by adding INFO fields : \n " +
        " - LOCALISATION: PASS or a combination of jnct5nt (near exon junction), homodimer5 ( in homopolymer), cl3snp35nt(in SNP cluster) \n" + 
        " - GTpopCR= X% : percentage of individuals genotyped \n" +
           " - DPgt05rPopCR= X% : percentage of individuals genotyped with DP>5\n")
    #     Inputs
    group_input = parser.add_argument_group( 'Inputs' )
    group_input.add_argument('--in-vcf', required=True, help="The input vcf file to annotate.")
    group_input.add_argument('--in-gtf', required=True, help="The genome annotations GTF file, to tag SNP closed to exon boundaries (+/- 5bp)")
    group_input.add_argument('--in-fasta', required=True, help="The genome Fasta file, to tag SNP in homopolymer ( 5 bp).")
    #     Outputs
    group_output = parser.add_argument_group( 'Outputs' )
    group_output.add_argument('--out-vcf', required=True, help="The annotated vcf file output.")
    group_output.add_argument('--out-vcf-gt', required=True, help="The annotated vcf file output filtered on GTpopCR > 50 and DPgt05rPopCR > 20.")
    group_output.add_argument('--out-gt', required=True, help="A GT table of variants with FILTER PASS, GTpopCR > 50 and DPgt05rPopCR > 20 ")
    group_output.add_argument('--out-dp', required=True, help="A DP table of variants with FILTER PASS, GTpopCR > 50 and DPgt05rPopCR > 20 ")
    group_output.add_argument('--out-ad', required=True, help="A AD table of variants with FILTER PASS, GTpopCR > 50 and DPgt05rPopCR > 20 ")
    group_output.add_argument('--summary', required=True, help="Annotation summary")
    args = parser.parse_args()

    process(args)
