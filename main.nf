#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow{

	// Create channel for bed file inputs
	OQs_K_bed_plus_ch = Channel.fromPath("${projectDir}/../K/*plus.hits.max.K.w50.25.bed", checkIfExists: true)
	OQs_K_bed_minus_ch = Channel.fromPath("${projectDir}/../K/*minus.hits.max.K.w50.25.bed", checkIfExists: true)
	OQs_K_PDS_bed_plus_ch = Channel.fromPath("${projectDir}/../K-PDS/*plus.hits.max.PDS.w50.35.bed", checkIfExists: true)
	OQs_K_PDS_bed_minus_ch = Channel.fromPath("${projectDir}/../K-PDS/*minus.hits.max.PDS.w50.35.bed", checkIfExists: true)

	convert_marsico_peaks(OQs_K_bed_plus_ch,
						  OQs_K_bed_minus_ch, 
						  OQs_K_PDS_bed_plus_ch,
						  OQs_K_PDS_bed_minus_ch, 
						  params.target_genome, 
						  params.original_genome)
}

process convert_marsico_peaks {
	tag "Converting Marsico Peaks"
	conda "envs/env.yml"

	publishDir "${projectDir}/results", mode: 'copy'

	input:
	path(OQs_K_bed_plus)
	path(OQs_K_bed_minus)
	path(OQs_PDS_bed_plus)
	path(OQs_PDS_bed_minus)
	path target_genome
	path original_genome

	output:
	path "Marsico_OQs_K_combined_v10.bed" 
	path "Marsico_OQs_K_PDS_combined_v10.bed"
	path "Marsico_OQs_K_combined_v10_circos.bed" 
	path "Marsico_OQs_K_PDS_combined_v10_circos.bed"

	script:
	"""
	awk -v s="+" 'BEGIN{OFS="\\t"} \
	/^(track|browser|#)/ {print; next} \
	{ print \$1,\$2,\$3,".", \$4, s }' ${OQs_K_bed_plus} > OQs_K_bed_plus.bed

	awk -v s="-" 'BEGIN{OFS="\\t"} \
	/^(track|browser|#)/ {print; next} \
	{ print \$1,\$2,\$3,".", \$4, s }' ${OQs_K_bed_minus} > OQs_K_bed_minus.bed

	awk -v s="+" 'BEGIN{OFS="\\t"} \
	/^(track|browser|#)/ {print; next} \
	{ print \$1,\$2,\$3,".", \$4, s }' ${OQs_PDS_bed_plus} > OQs_PDS_bed_plus.bed

	awk -v s="-" 'BEGIN{OFS="\\t"} \
	/^(track|browser|#)/ {print; next} \
	{ print \$1,\$2,\$3,".", \$4, s }' ${OQs_PDS_bed_minus} > OQs_PDS_bed_minus.bed

	cat OQs_K_bed_plus.bed OQs_K_bed_minus.bed | sort-bed - > Marsico_OQs_K_combined.bed
	cat OQs_PDS_bed_plus.bed OQs_PDS_bed_minus.bed | sort-bed - > Marsico_OQs_K_PDS_combined.bed
	
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
            attributes = "score=" score
            
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
            attributes = "score=" score
            
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

	#Now convert the lifted GFF3 files into bed files with the structure chrom, start, end, score
	sed -i '' '/^#/d' Marsico_OQs_K_combined_v10.gff3
	sed -i '' '/^#/d' Marsico_OQs_K_PDS_combined_v10.gff3
	awk '{
		split(\$9, attrs, ";");
		score = ".";
		for (i in attrs) {
			if (attrs[i] ~ /^score=/) {
				score = substr(attrs[i], 7);
				break;
			}
		}
		print \$1 "\\t" \$4-1 "\\t" \$5 "\\t" "." "\\t" score "\\t" \$7
	}' "Marsico_OQs_K_combined_v10.gff3" > "Marsico_OQs_K_combined_v10.bed"
	
	awk '{
		split(\$9, attrs, ";");
		score = ".";
		for (i in attrs) {
			if (attrs[i] ~ /^score=/) {
				score = substr(attrs[i], 7);
				break;
			}
		}
		print \$1 "\\t" \$4-1 "\\t" \$5 "\\t" "." "\\t" score "\\t" \$7
	}' "Marsico_OQs_K_PDS_combined_v10.gff3" > "Marsico_OQs_K_PDS_combined_v10.bed"

	awk -F'\\t' 'BEGIN {OFS="\\t"} \$6=="-" {\$5=-\$5} {print \$0}' Marsico_OQs_K_combined_v10.bed > Marsico_temp_K.bed
	awk -F'\\t' 'BEGIN {OFS="\\t"} \$6=="-" {\$5=-\$5} {print \$0}' Marsico_OQs_K_PDS_combined_v10.bed > Marsico_temp_K_PDS.bed
	
	cut -f 1,2,3,5 Marsico_temp_K.bed > Marsico_OQs_K_combined_v10_circos.bed
	cut -f 1,2,3,5 Marsico_temp_K_PDS.bed > Marsico_OQs_K_PDS_combined_v10_circos.bed

	sed -i '' 's/_Tb427v10//g' Marsico_OQs_K_combined_v10_circos.bed
	sed -i '' 's/_Tb427v10//g' Marsico_OQs_K_PDS_combined_v10_circos.bed

	"""
}