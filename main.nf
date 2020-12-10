$HOSTNAME = ""
params.outdir = 'results'  


if (!params.reads){params.reads = ""} 
if (!params.mate){params.mate = ""} 

Channel
	.fromFilePairs( params.reads , size: params.mate == "single" ? 1 : params.mate == "pair" ? 2 : params.mate == "triple" ? 3 : params.mate == "quadruple" ? 4 : -1 )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.set{g_1_reads_g_0}

Channel.value(params.mate).set{g_2_mate_g_0}


process Map_Minimap2 {

input:
 set val(name), file(reads) from g_1_reads_g_0
 val mate from g_2_mate_g_0


# User provided inputs:
# $5: Minimap2 Index of reference genome (e.g., ~/Ref/hg38/genome.mmi)
# $6: Fasta file of reference genome (e.g., ~/Ref/hg38/genome.fa)
# $7: Lowest used Barcode (01)
# $8: Highest used Barcode (08)

#Alignment using minimap2 + generation of sortedn BAM + SAM
#files
for bc in $(seq -w $7 $8)
do

minimap2  $5 $2/Demultiplex/barcode$bc/*.fastq --MD >$2/Alignment/barcode$bc.sam

when:
(params.run_Minimap2 && (params.run_Minimap2 == "yes")) || !params.run_Minimap2

script:
params_Minimap2 = params.Map_Minimap2.params_Minimap2

"""
$runGzip
minimap2 ${params_Minimap2}  --genomeDir ${params.minimap2_index} $reads >barcode${name}.sam
echo "Alignment completed."


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
