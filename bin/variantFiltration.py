#!/usr/bin/env python

import argparse
import os,sys
import time
import subprocess
from subprocess import Popen, PIPE
import multiprocessing

##########################
# FUNCTION
def select_on_chrom(in_vcf, accepted_chrom, excluded_chrom, selected_vcf):
    
    #~ FH_in = open(in_vcf, encoding='utf8', 'rt')
    # pour FLLL je ne comprends pas pourquoi mais utf8 ne fonctionne pas ... 
    FH_in = open(in_vcf, encoding='iso-8859-1', mode='rt')
    FH_out = open(selected_vcf,'wt')
    
    for line in FH_in:
        if line.startswith('#'):
            FH_out.write(line)
        elif accepted_chrom and line.split()[0] in accepted_chrom :
            FH_out.write(line)
        elif excluded_chrom and line.split()[0] not in excluded_chrom :
            FH_out.write(line)

def split_vcf(in_vcf, nb_cpus, tmp_files ):
    
    # count variant
    cmd = 'grep -vc "#" ' + in_vcf
    p = Popen(cmd, shell=True, universal_newlines=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    nb_var = int(stdout.strip())
    nb_by_file = int(nb_var / nb_cpus)
    
    # write splitted VCF
    tmp_FH = [ open(f,'wt') for f in tmp_files]
    FH_in = open(in_vcf)
    header = list()
    
    c=0
    var=list()

    for line in FH_in:
        if line.startswith('#'):
            header.append(line)
        else:
            if len(header) > 0:
                for t in tmp_FH:
                    t.write("".join(header))
                header.clear()
            c += 1
            var.append(line)
            
            # write variant by batch of 10000
            if len(tmp_FH)>0 and len(var) == min(10000,nb_by_file):
                tmp_FH[0].write(''.join(var))
                var.clear()
                
            # change temp output file
            if len(tmp_FH)>0 and c == nb_by_file:
                c = 0
                tmp_FH[0].close()
                del(tmp_FH[0])

    
    # write remaining var
    if len(var) > 0:
        tmp_FH = open(tmp_files[-1],'a')
        tmp_FH.write(''.join(var))
        var.clear()
        tmp_FH.close()
    
def submit_cmd( cmd, cwd=None):
    """
    @summary: Submits a command on system.
    @param cmd: [list] The command.
    @param stdout_path: The path to the file where the standard outputs will be written.
    @param stderr_path: The path to the file where the error outputs will be written.
    """
    p = Popen(cmd, universal_newlines=True, stdout=PIPE, stderr=PIPE, cwd=cwd)
    stdout, stderr = p.communicate()

    # check error status
    if p.returncode != 0:
        raise Exception( stderr )
        
def process_variantFiltering_and_Annot(in_vcf, in_gtf, in_fasta, out_vcf, out_vcf_gt, out_gt, out_dp, out_ad, out_summary):

    cmd = ['vcf_filter_and_annot.py',  '--in-vcf',  in_vcf, '--in-gtf', in_gtf, '--in-fasta', in_fasta, 
        '--out-vcf', out_vcf, '--out-vcf-gt', out_vcf_gt, '--out-gt', out_gt, '--out-dp', out_dp, '--out-ad', out_ad, '--summary', out_summary]
    print(" ".join(cmd) + "\n")
    submit_cmd( cmd )


def parallel_submission( function, inputs_vcf, in_gtf, in_fasta, outputs_vcf, outputs_vcf_gt, outputs_gt, outputs_dp, outputs_ad, outputs_summary):
    cpu_used = len(inputs_vcf)
    processes = [{'process':None, 'input_vcf':None, 'gtf' : in_gtf, 'fasta':in_fasta,
                  'output_vcf':None, 'output_vcf_gt':None, 'output_gt':None, 'output_dp':None, 'output_ad':None, 'output_summary':None } for idx in range(cpu_used)]
    
    # Launch processes
    for idx in range(cpu_used):
        process_idx = idx % cpu_used
        processes[process_idx]['input_vcf'] = inputs_vcf[idx]
        processes[process_idx]['output_vcf'] = outputs_vcf[idx]
        processes[process_idx]['output_vcf_gt'] = outputs_vcf_gt[idx]
        processes[process_idx]['output_gt'] = outputs_gt[idx]
        processes[process_idx]['output_dp'] = outputs_dp[idx]
        processes[process_idx]['output_ad'] = outputs_ad[idx]
        processes[process_idx]['output_summary'] = outputs_summary[idx]

    for current_process in processes:
        if idx == 0:  # First process is threaded with parent job
            current_process['process'] = threading.Thread(target=function,
                                                          args=(current_process['input_vcf'],current_process['gtf'], current_process['fasta'], current_process['output_vcf'], current_process['output_vcf_gt'], current_process['output_gt'], current_process['output_dp'], current_process['output_ad'], current_process['output_summary']))
        else:  # Others processes are processed on diffrerent CPU
            current_process['process'] = multiprocessing.Process(target=function,
                                                                 args=(current_process['input_vcf'],current_process['gtf'], current_process['fasta'], current_process['output_vcf'], current_process['output_vcf_gt'], current_process['output_gt'], current_process['output_dp'], current_process['output_ad'], current_process['output_summary']))
        current_process['process'].start()
    # Wait processes end
    for current_process in processes:
        current_process['process'].join()
    # Check processes status
    for current_process in processes:
        if issubclass(current_process['process'].__class__, multiprocessing.Process) and current_process['process'].exitcode != 0:
            raise Exception("Error in sub-process execution.")

def merge_vcf(in_list,out):
    
    cmd = ['bcftools','concat','-o',out,'-O','v'] + in_list
    print ("\n#merge temp VCFs\n" + " ".join(cmd) + "\n")
    submit_cmd( cmd  )

def merge_tsv(in_list, out):
    
    print ("\n#merge temp tables\n\tInputs : " + " ".join(in_list) + "\n\tOutputs : " + out + "\n")
    
    FH_out = open(out,'wt')
    first = True
    
    for f in in_list:
        FH_in = open(f,'rt')
        header=FH_in.readline()
        
        if first :
            FH_out.write(header)
            first = False
            
        for line in FH_in:
            FH_out.write(line)
            
def merge_summary(in_list, out):
    
    print ("\n#merge temp summaries\n\tInputs : " + " ".join(in_list) + "\n\tOutputs : " + out + "\n")
    
    dic_results = dict()
    first=in_list.pop(0)
    FH_first = open(first,'rt')
    # prefix = tmp_dir/time_pid_i_
    prefix = os.path.join(os.path.dirname(first), '_'.join(os.path.basename(first).split('_')[0:3])) + '_'
    summary = ''
    first_key = ''
    percent_list = list()
    for line in FH_first:
        if prefix in line:
            summary += line.replace(prefix,'')
            
        elif ':' in line:
            category = line.split(':')[0]
            key = category.upper().strip().replace(' ','_')
            
            # store first key concidered as the total of variant
            if first_key == '':
                first_key = key
            
            # store key where we want to compute percentage
            if '%)' in line:
                percent_list.append(key)
                
            value = int(line.split(':')[1].strip().split()[0])
            
            summary += category + ': ###' + key + '###\n'
            dic_results[key] = value
        
        else:
            summary += line

    for f in in_list:
        FH_in = open(f,'rt')
        for line in FH_in:
            if ':' in line:
                category = line.strip().split(':')[0]
                key = category.upper().strip().replace(' ','_')
                value = int(line.strip().split(':')[1].strip().split()[0])
                dic_results[key] += value

    for k in dic_results:
        if k in percent_list:
            percent = round(dic_results[k]*100/dic_results[first_key],2)
            summary = summary.replace('###' + k + '###' , str(dic_results[k]) + ' (' + str(percent) + '%)')
        else:
            summary = summary.replace('###' + k + '###' , str(dic_results[k]))
            
    FH_out = open(out,'wt')
    FH_out.write(summary)
    FH_out.close()
    
def process(params):
    
    in_vcf=params.in_vcf
    out_dir = os.path.dirname(params.out_vcf)
    temp_dir = os.path.join(out_dir,str(time.time()) + "_" + str(os.getpid()) + '_temp_variantFiltration')
    temp_files = list()
    try:
        # select variant based on chrom
        if params.accepted_chrom is not None or params.excluded_chrom is not None:
            os.mkdir(temp_dir)
            selected_vcf = os.path.join(temp_dir, str(time.time()) + "_" + str(os.getpid()) + "_selectOnChrom_" + os.path.basename(in_vcf))
            temp_files.append(selected_vcf)
            select_on_chrom(in_vcf, params.accepted_chrom, params.excluded_chrom, selected_vcf)
            in_vcf = selected_vcf
            
        if params.nb_cpus == 1:
            process_variantFiltering_and_Annot(in_vcf, params.in_gtf, params.in_fasta, params.out_vcf, params.out_vcf_gt, params.out_gt, params.out_dp, params.out_ad, params.summary)
        else :
            #### create temp_dir
            if not os.path.exists(temp_dir):
                os.mkdir(temp_dir)
            
            #### split VCF in [nb_cpus] files
            prefix = str(time.time()) + "_" + str(os.getpid())
            splitted_in_vcf_files = [ os.path.join(temp_dir, prefix + "_" + str(i) + "_" + os.path.basename(in_vcf)) for i in range(0,params.nb_cpus )]
            temp_files += splitted_in_vcf_files
            split_vcf(in_vcf, params.nb_cpus, splitted_in_vcf_files )

            #### filter_and_annotate
            splitted_out_vcf_files = [ os.path.join(temp_dir, prefix + "_" + str(i) + "_" + os.path.basename(params.out_vcf)) for i in range(0,params.nb_cpus)]
            splitted_out_vcf_gt_files = [ os.path.join(temp_dir, prefix + "_" + str(i) + "_" + os.path.basename(params.out_vcf_gt)) for i in range(0,params.nb_cpus)]
            splitted_out_gt_files = [ os.path.join(temp_dir, prefix + "_" + str(i) + "_" + os.path.basename(params.out_gt)) for i in range(0,params.nb_cpus)]
            splitted_out_dp_files = [ os.path.join(temp_dir, prefix + "_" + str(i) + "_" + os.path.basename(params.out_dp)) for i in range(0,params.nb_cpus)]
            splitted_out_ad_files = [ os.path.join(temp_dir, prefix + "_" + str(i) + "_" + os.path.basename(params.out_ad)) for i in range(0,params.nb_cpus)]
            splitted_out_summary_files = [ os.path.join(temp_dir, prefix + "_" + str(i) + "_" + os.path.basename(params.summary)) for i in range(0,params.nb_cpus)]
            temp_files += splitted_out_vcf_files + splitted_out_vcf_gt_files + splitted_out_gt_files + splitted_out_ad_files + splitted_out_dp_files + splitted_out_summary_files
            
            parallel_submission( process_variantFiltering_and_Annot, splitted_in_vcf_files, params.in_gtf, params.in_fasta, splitted_out_vcf_files, splitted_out_vcf_gt_files, splitted_out_gt_files, splitted_out_dp_files, splitted_out_ad_files, splitted_out_summary_files )
            
            # merge outputs
            merge_vcf(splitted_out_vcf_files, params.out_vcf)
            merge_vcf(splitted_out_vcf_gt_files, params.out_vcf_gt)
            merge_tsv(splitted_out_gt_files, params.out_gt)
            merge_tsv(splitted_out_dp_files, params.out_dp)
            merge_tsv(splitted_out_ad_files, params.out_ad)
            merge_summary(splitted_out_summary_files, params.summary)
        
    finally :
        # remove temp_dir
        if not params.debug and os.path.exists(temp_dir):
            # remove temp files
            for f in temp_files:
                if os.path.exists(f):
                    os.remove(f)
            # remove temp dir
            os.rmdir(temp_dir)
        

##########################
# MAIN

if __name__ == '__main__':
    # Parameters
    parser = argparse.ArgumentParser(description="Annot variant, by adding INFO fields : \n " +
        " - LOCALISATION: PASS or a combination of jnct5nt (near exon junction), homodimer5 ( in homopolymer), cl3snp35nt(in SNP cluster) \n" + 
        " - GTpopCR= X% : percentage of individuals genotyped \n" +
           " - DPgt05rPopCR= X% : percentage of individuals genotyped with DP>5\n")
    parser.add_argument('-n','--nb-cpus', type=int, default=1, help='Number of CPU to paralelise the treatment')
    parser.add_argument('-d','--debug', default=False, action='store_true', help='keep temporary files')
    parser.add_argument('-a','--accepted-chrom', nargs="+", required=False, help='list of chromosome to include')
    parser.add_argument('-e','--excluded-chrom', nargs="+", required=False, help='list of chromosome to exclude')
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

    if args.accepted_chrom is not None and args.excluded_chrom is not None:
        intersect_chrom = [ c for c in args.accepted_chrom if c in args.excluded_chrom]
        if len(intersect_chrom) > 0:
            raise Exception("accepted_chrom and excluded_chrom list must not intersect\n")

    process(args)
