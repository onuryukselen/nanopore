$HOSTNAME = ""
params.outdir = 'results'  

$HOSTNAME = "default"
params.DOWNDIR = (params.DOWNDIR) ? params.DOWNDIR : ""
//* params.genome_build =  ""  //* @dropdown @options:"human_hg19, human_hg38, macaque_macFas5, rat_rn6, rat_rn6ens, mousetest_mm10, custom"
//* params.run_Minimap2 =  "yes"  //* @dropdown @options:"yes","no" @show_settings:"Map_Minimap2"
//* params.run_checkAndBuild =  "yes"  //* @dropdown @options:"yes","no" @show_settings:"Build_Minimap2_Index"
//* params.run_Talon =  "yes"  //* @dropdown @options:"yes","no" @show_settings:"talon_initialize_database","talon"


_species = ""
_build = ""
_share = ""
//* autofill
if (params.genome_build == "mousetest_mm10"){
    _species = "mousetest"
    _build = "mm10"
} else if (params.genome_build == "human_hg19"){
    _species = "human"
    _build = "hg19"
} else if (params.genome_build == "human_hg38"){
    _species = "human"
    _build = "hg38"
} else if (params.genome_build == "mouse_mm10"){
    _species = "mouse"
    _build = "mm10"
} else if (params.genome_build == "macaque_macFas5"){
    _species = "macaque"
    _build = "macFas5"
} else if (params.genome_build == "rat_rn6"){
    _species = "rat"
    _build = "rn6"
} else if (params.genome_build == "rat_rn6ens"){
    _species = "rat"
    _build = "rn6ens"
}
if ($HOSTNAME == "default"){
    _share = "${params.DOWNDIR}/genome_data"
    $SINGULARITY_IMAGE = "dolphinnext/nanopore:1.0"
    $CPU  = 1
    $MEMORY = 10
}

if ($HOSTNAME == "fs-bb7510f0"){
    _share = "/mnt/efs/share/genome_data"
    $SINGULARITY_IMAGE = "dolphinnext/nanopore:1.0"
} else if ($HOSTNAME == "192.168.20.150"){
    _share = "/home/botaoliu/share/genome_data"
    $SINGULARITY_IMAGE = "dolphinnext/nanopore:1.0"
} else if ($HOSTNAME == "50.228.141.2"){
    _share = "/share/genome_data"
    $SINGULARITY_IMAGE = "dolphinnext/nanopore:1.0"
} else if ($HOSTNAME == "ghpcc06.umassrc.org"){
    _share = "/share/data/umw_biocore/genome_data"
    $SINGULARITY_IMAGE = "dolphinnext/nanopore:1.0"
    $TIME = 500
    $CPU  = 1
    $MEMORY = 32 
    $QUEUE = "long"
}
if (params.genome_build && $HOSTNAME){
    params.genome ="${_share}/${_species}/${_build}/${_build}.fa"
    params.gtf ="${_share}/${_species}/${_build}/ucsc.gtf"
    params.minimap2_index = "${_share}/${_species}/${_build}/Minimap2Index/genome.mmi"
}
if ($HOSTNAME){
    params.samtools_path = "samtools"
}
//*


if (!params.reads){params.reads = ""} 
if (!params.mate){params.mate = ""} 

Channel
	.fromFilePairs( params.reads , size: params.mate == "single" ? 1 : params.mate == "pair" ? 2 : params.mate == "triple" ? 3 : params.mate == "quadruple" ? 4 : -1 )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.set{g_1_reads_g_0}

Channel.value(params.mate).set{g_2_mate_g_0}

//* params.gtf =  ""  //* @input


process talon_initialize_database {


output:
 file "Database.db"  into g_6_database_g_7
 val genome_build  into g_6_genome_build_g_7
 val platform  into g_6_platform_g_0
 val annotation_name  into g_6_annotation_g_7

when:
(params.run_Talon && (params.run_Talon == "yes")) || !params.run_Talon

script:

custom_genome_build = params.talon_initialize_database.custom_genome_build
annotation = params.talon_initialize_database.annotation
genome_build = custom_genome_build != "" ? custom_genome_build : params.genome_build
annotation_name = annotation != "" ? annotation : params.genome_build
platform = params.talon_initialize_database.platform
"""
talon_initialize_database --f ${params.gtf} --g ${genome_build} --a ${annotation_name} --o Database
"""

}

//* params.genome =  ""  //* @input
//* params.minimap2_index =  ""  //* @input

process Build_Minimap2_Index {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.mmi$/) "Minimap2_Index/$filename"
}


output:
 file "*.mmi"  into g_3_index_g_0

when:
params.run_checkAndBuild == "yes"

script:
params_build_Minimap2 = params.Build_Minimap2_Index.params_build_Minimap2
indexPath = params.minimap2_index.toString()
if (indexPath.contains('/')) {
	indexDir = indexPath.substring(0, indexPath.lastIndexOf('/'))
} else {
	indexDir = "index"
}
"""
if [ ! -e "${params.minimap2_index}" ] ; then
    echo "Minimap2 index not found"
    minimap2 ${params_build_Minimap2} genome.mmi ${params.genome}
    echo "Minimap2 index built."
    mkdir -p ${indexDir} 2> /dev/null || true
    if [ -w "${indexDir}" ]; then
    	cp genome.mmi ${params.minimap2_index}
    fi
else 
	cp ${params.minimap2_index} genome.mmi
fi
"""




}

g_3_index_g_0= g_3_index_g_0.ifEmpty([""]) 
g_6_platform_g_0= g_6_platform_g_0.ifEmpty([""]) 


process Map_Minimap2 {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}_sorted.(bam|bam.bai)$/) "Alignment/$filename"
}

input:
 set val(name), file(reads) from g_1_reads_g_0
 val mate from g_2_mate_g_0
 file minimap_index from g_3_index_g_0
 val platform from g_6_platform_g_0

output:
 file "${name}_sorted.sam"  into g_0_sam_file_g_7
 file "${name}_sorted.{bam,bam.bai}"  into g_0_bam_bai
 file "${name}_config.csv"  into g_0_config_g_7

when:
(params.run_Minimap2 && (params.run_Minimap2 == "yes")) || !params.run_Minimap2

script:
params_Minimap2 = params.Map_Minimap2.params_Minimap2
index = minimap_index.toString().matches("(.*).mmi") ? minimap_index : params.minimap2_index
"""
minimap2 ${params_Minimap2} ${index} $reads > ${name}.sam
echo "Alignment completed."
#Convert SAM to BAM and remove SAM
samtools view -bT ${params.genome} ${name}.sam > ${name}.bam
rm ${name}.sam

#Sort & Index:
samtools sort ${name}.bam -o ${name}_sorted.bam
samtools index ${name}_sorted.bam ${name}_sorted.bam.bai
rm ${name}.bam

# Convert into a single SAM file
samtools view -h ${name}_sorted.bam > ${name}_sorted.sam

# Talon Config line:
# SIRV_Rep1,SIRV,PacBio-Sequel2,labeled/SIRV_rep1_labeled.sam
echo "${name},${name},${platform},${name}_sorted.sam" > ${name}_config.csv

"""


}

//* params.gtf =  ""  //* @input


process talon {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*$/) "Talon/$filename"
}

input:
 file db from g_6_database_g_7
 val genome_build from g_6_genome_build_g_7
 file sams from g_0_sam_file_g_7.collect()
 file configs from g_0_config_g_7.collect()
 val annotation from g_6_annotation_g_7

output:
 file "*"  into g_7_outputDir

script:

MIN_COVERAGE = params.talon.MIN_COVERAGE
MIN_IDENTITY = params.talon.MIN_IDENTITY
THREADS = params.talon.THREADS
"""
awk 'FNR==1 && NR!=1 {  getline; } 1 {print} ' *.csv > conf.csv
# talon writes on db file. Copied for resume functionality.
cp ${db} talon.db
talon --f conf.csv --db talon.db --build ${genome_build} --cov $MIN_COVERAGE --identity $MIN_IDENTITY --o . --threads ${THREADS} > talon.log 2>&1
echo "Talon completed."
#Summarize TALON.
talon_summarize --db talon.db --v --o . > talon_summarize.log 2>&1
echo "talon_summarize completed."

# Export isoform abundance
talon_abundance --db talon.db  -a ${annotation} -b ${genome_build} --o talon_abundance > talon_abundance.log 2>&1
echo "talon_abundance completed."

# Generate GTF
talon_create_GTF --db talon.db -b ${genome_build} -a ${annotation} --observed --o talon_create_GTF > talon_create_GTF.log 2>&1
echo "talon_create_GTF completed."

"""

}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
