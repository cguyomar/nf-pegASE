params.options = [:]

process FILTER_JUNCTIONS {
    tag "$meta.id"
    
    input:
    tuple val(meta), path(tab)

    output:
    tuple val(meta), path('*.SJ.out.filtered.tab') , emit:sjout

    script:
    """
    # filter junctions on: mitochondria , non canonical, already known, covered by less than 3 uniquely mapped reads in at least 1 sample
    cat $tab | awk '\$6==0 && \$5>0 && \$7>=2' | cut -f1-6 | sort | uniq -c |awk '{if(\$1>1){print \$2,\$3,\$4,\$5,\$6,\$7}}' |sed 's/ /	/g' > ${meta.id}.SJ.out.filtered.tab
    """
}
