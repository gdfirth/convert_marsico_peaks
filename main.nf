#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow{

	// Create channel for bed file inputs
	OQs_K_bed_ch = Channel.fromPath("${projectDir}/../K/*.hits.max.K.w50.25.bed", checkIfExists: true)
	OQs_K_PDS_bed_ch = Channel.fromPath("${projectDir}/../K-PDS/*.hits.max.PDS.w50.35.bed", checkIfExists: true)

	convert_marsico_peaks(OQs_K_bed_ch.collect(), 
						  OQs_K_PDS_bed_ch.collect(), 
						  params.target_genome, 
						  params.original_genome)
}

process convert_marsico_peaks {
	tag "Converting Marsico Peaks"
	conda "envs/env.yml"

	publishDir "${projectDir}/results", mode: 'copy'

	input:
	tuple path(OQs_K_bed_1), path(OQs_K_bed_2)
	tuple path(OQs_PDS_bed_1), path(OQs_PDS_bed_2)
	path target_genome
	path original_genome

	output:
	path "Marsico_OQs_K_combined.gff3" 
	path "Marsico_OQs_K_PDS_combined.gff3"

	script:
	"""
	cat ${OQs_K_bed_1} ${OQs_K_bed_2} | sort-bed - > Marsico_OQs_K_combined.bed
	cat ${OQs_PDS_bed_1} ${OQs_PDS_bed_2} | sort-bed - > Marsico_OQs_K_PDS_combined.bed
	
	# Rename chromosome "TP26M21-2a10.p1k_v5.1" to "TP26M21-2a10.p1k" in both BED files
	sed -i '' 's/TP26M21-2a10.p1k_v5.1/TP26M21-2a10.p1k/g' Marsico_OQs_K_combined.bed
	sed -i '' 's/TP26M21-2a10.p1k_v5.1/TP26M21-2a10.p1k/g' Marsico_OQs_K_PDS_combined.bed
	#NOTE: If you're using Linux, remove the '' after -i in the sed commands above.
	#NOTE: THIS IS A KNOWN DIFFERENCE BETWEEN MACOS SED AND LINUX SED.


	# Convert BED to GFF3 for Marsico_OQs_K_combined.bed-------------------------------
    echo "##gff-version 3" > Marsico_OQs_K_combined.gff3
    
    # Convert BED to GFF3
    awk -v feature_type="${params.desired_feature_name}" 'BEGIN {OFS="\\t"} 
    {
        if (\$0 !~ /^#/ && NF >= 3) {
            # BED is 0-based, GFF is 1-based
            start = \$2 + 1
            end = \$3
            
            # Handle optional fields
            name = (NF >= 4 && \$4 != "") ? \$4 : "feature_" NR
            score = (NF >= 5 && \$5 != "") ? \$5 : "."
            strand = (NF >= 6 && \$6 != "") ? \$6 : "."
            
            # Create attributes
            attributes = "score=" name
            
            print \$1, "BED_conversion", feature_type, start, end, score, strand, ".", attributes
        }
    }' Marsico_OQs_K_combined.bed >> Marsico_OQs_K_combined.gff3

	# Convert BED to GFF3 for Marsico_OQs_K_PDS_combined.bed-------------------------------
	echo "##gff-version 3" > Marsico_OQs_K_PDS_combined.gff3
    
    # Convert BED to GFF3
    awk -v feature_type="${params.desired_feature_name}" 'BEGIN {OFS="\\t"} 
    {
        if (\$0 !~ /^#/ && NF >= 3) {
            # BED is 0-based, GFF is 1-based
            start = \$2 + 1
            end = \$3
            
            # Handle optional fields
            name = (NF >= 4 && \$4 != "") ? \$4 : "feature_" NR
            score = (NF >= 5 && \$5 != "") ? \$5 : "."
            strand = (NF >= 6 && \$6 != "") ? \$6 : "."
            
            # Create attributes
            attributes = "score=" name
            
            print \$1, "BED_conversion", feature_type, start, end, score, strand, ".", attributes
        }
    }' Marsico_OQs_K_PDS_combined.bed >> Marsico_OQs_K_PDS_combined.gff3

	echo "${params.desired_feature_name}" > desired_features.txt

	liftoff -g Marsico_OQs_K_combined.gff3 \
	-o Marsico_OQs_K_combined_v10.gff3 \
	-f desired_features.txt \
	${target_genome} \
	${original_genome}

	liftoff -g Marsico_OQs_K_PDS_combined.gff3 \
	-o Marsico_OQs_K_PDS_combined_v10.gff3 \
	-f desired_features.txt \
	${target_genome} \
	${original_genome}

	"""
}