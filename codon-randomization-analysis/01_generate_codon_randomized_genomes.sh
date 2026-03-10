#mamba create -n biopython
#mamba install biopython pandas

mamba activate biopython

mkdir -p processed_data/output_gene_ssrs_codon_resampled_by_genome
mkdir -p processed_data/output_gene_ssrs_actual
mkdir -p processed_data/output_gene_gc_content

cd raw_data/Basel_deduplicated_genbank
perl ../../batch_run.pl -0 "python ../../randomize_genome.py  --minimum-length 1 --randomization-method genome --resamplings 1000 --input #d --output-gc-content ../../processed_data/output_gene_gc_content/#e.csv --output-actual ../../processed_data/output_gene_ssrs_actual/#e.csv --output-resampled ../../processed_data/output_gene_ssrs_codon_resampled_by_genome/#e.csv"


python ../tabulate_SSRs.py  --minimum-length 6 --input Basel_deduplicated_genbank --output ../processed_data/genomewide_ssrs_actual.csv

