nextflow {

    params {
        reads = "data/*_{1,2}.fastq.gz" // Path to paired-end reads
        ref = "genome.fasta"          // Reference genome
        outdir = "results"            // Output directory
        threads = 8                   // Number of threads for tools
        admixture_k = 3               // Number of ancestral populations for ADMIXTURE
    }

    process INDEX_REFERENCE {
        input:
        path ref

        output:
        path "*.bwt", "*.ann", "*.amb", "*.pac", "*.sa"

        script:
        """
        bwa index ${ref}
        """
    }

    process ALIGN_READS {
        input:
        tuple path(reads1), path(reads2)
        path ref
        path index

        output:
        path "*.sam"

        script:
        """
        bwa mem -t ${task.cpus} ${ref} ${reads1} ${reads2} > ${reads1.simpleName.replace('_1', '')}.sam
        """
    }

    process CONVERT_SAM_TO_BAM {
        input:
        path sam

        output:
        path "*.bam"

        script:
        """
        samtools view -@ ${task.cpus} -Sb ${sam} > ${sam.simpleName}.bam
        """
    }

    process SORT_BAM {
        input:
        path bam

        output:
        path "*.sorted.bam"

        script:
        """
        samtools sort -@ ${task.cpus} ${bam} > ${bam.simpleName}.sorted.bam
        """
    }

    process MARK_DUPLICATES {
        input:
        path sorted_bam

        output:
        path "*.dedup.bam", path "*.metrics"

        script:
        """
        gatk MarkDuplicates -I ${sorted_bam} -O ${sorted_bam.simpleName}.dedup.bam -M ${sorted_bam.simpleName}.metrics
        """
    }

    process CALL_VARIANTS {
        input:
        path dedup_bam
        path ref

        output:
        path "*.vcf"

        script:
        """
        gatk HaplotypeCaller -R ${ref} -I ${dedup_bam} -O ${dedup_bam.simpleName}.vcf
        """
    }

    process FILTER_VARIANTS {
        input:
        path vcf

        output:
        path "*.filtered.vcf"

        script:
        """
        gatk VariantFiltration -V ${vcf} -O ${vcf.simpleName}.filtered.vcf --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" --filter-name "basic_snp_filter"
        """
    }

    process POPULATION_GENOMICS {
        input:
        path filtered_vcf

        output:
        path "*.plink"

        script:
        """
        plink --vcf ${filtered_vcf} --make-bed --out ${filtered_vcf.simpleName}
        """
    }

    process PCA_ANALYSIS {
        input:
        path plink_files

        output:
        path "*.pca"

        script:
        """
        plink --bfile ${plink_files.simpleName} --pca 10 --out ${plink_files.simpleName}
        """
    }

    process ADMIXTURE_ANALYSIS {
        input:
        path plink_files

        output:
        path "*.Q", path "*.P"

        script:
        """
        admixture --cv ${plink_files.simpleName}.bed ${params.admixture_k} | tee ${plink_files.simpleName}.log
        """
    }

    process PLOT_ADMIXTURE {
        input:
        path q_file
        path sample_metadata

        output:
        path "admixture_plot.pdf"

        script:
        """
        Rscript -e 'library(pophelper); \
                  qdata <- readQ(files = "${q_file}"); \
                  metadata <- read.table("${sample_metadata}", header = TRUE); \
                  plotQ(qdata, grplab = metadata, exportpath = "admixture_plot.pdf")'
        """
    }

    workflow {
        ref_index = INDEX_REFERENCE(params.ref)

        reads_ch = Channel.fromFilePairs(params.reads, flat: false)
        reads_ch.map { pair ->
            ALIGN_READS(pair[0], pair[1], params.ref, ref_index)
        }
        .set { aligned_sam }

        aligned_sam.map { sam ->
            CONVERT_SAM_TO_BAM(sam)
        }
        .map { bam ->
            SORT_BAM(bam)
        }
        .map { sorted_bam ->
            MARK_DUPLICATES(sorted_bam)
        }
        .map { dedup_bam ->
            CALL_VARIANTS(dedup_bam, params.ref)
        }
        .map { vcf ->
            FILTER_VARIANTS(vcf)
        }
        .map { filtered_vcf ->
            POPULATION_GENOMICS(filtered_vcf)
        }
        .set { plink_files }

        plink_files.map { plink ->
            PCA_ANALYSIS(plink)
        }

        plink_files.map { plink ->
            ADMIXTURE_ANALYSIS(plink)
        }
        .map { q_file ->
            PLOT_ADMIXTURE(q_file, "metadata.txt")
        }
    }
}
