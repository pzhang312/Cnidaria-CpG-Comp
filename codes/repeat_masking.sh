#!/bin/bash

/home/each/RepeatScout-1/build_lmer_table -sequence Acropora_digitefera_GCF_000222465.1_unplaced.scaf.fna -freq Acropora_digitefera_GCF_000222465.1_unplaced_output_lmer.frequency
/home/each/RepeatScout-1/RepeatScout -sequence Acropora_digitefera_GCF_000222465.1_unplaced.scaf.fna -output Acropora_digitefera_GCF_000222465.1_unplaced_repeats.fasta -freq Acropora_digitefera_GCF_000222465.1_unplaced_output_lmer.frequency
cat Acropora_digitefera_GCF_000222465.1_unplaced_repeats.fasta | /home/each/RepeatScout-1/filter-stage-1.prl > Acropora_digitefera_GCF_000222465.1_unplaced_repeats_filtered1.fasta
/home/each/RepeatMasker/RepeatMasker -xm -xsmall -gff -s -pa 14 -lib Acropora_digitefera_GCF_000222465.1_unplaced_repeats_filtered1.fasta Acropora_digitefera_GCF_000222465.1_unplaced.scaf.fna
cat Acropora_digitefera_GCF_000222465.1_unplaced_repeats_filtered1.fasta | /home/each/RepeatScout-1/filter-stage-2.prl --cat=Acropora_digitefera_GCF_000222465.1_unplaced.scaf.fna.out > Acropora_digitefera_repeats_filtered2.fasta
/home/each/RepeatMasker/RepeatMasker -xm -xsmall -gff -pa 14 -s -lib Acropora_digitefera_repeats_filtered2.fasta -nolow -norna -no_is -gff Acropora_digitefera_GCF_000222465.1_unplaced.scaf.fna

/home/each/RepeatScout-1/build_lmer_table -sequence Acropora_millepora_GCF_004143615.1_amil_sf_1.1_genomic.fna -freq Acropora_millepora_GCF_004143615.1_output_lmer.frequency
/home/each/RepeatScout-1/RepeatScout -sequence Acropora_millepora_GCF_004143615.1_amil_sf_1.1_genomic.fna -output Acropora_millepora_GCF_004143615.1_repeats.fasta -freq Acropora_millepora_GCF_004143615.1_output_lmer.frequency
cat Acropora_millepora_GCF_004143615.1_repeats.fasta | /home/each/RepeatScout-1/filter-stage-1.prl > Acropora_millepora_GCF_004143615.1_repeats_filtered1.fasta
/home/each/RepeatMasker/RepeatMasker -xm -xsmall -gff -s -pa 14 -lib Acropora_millepora_GCF_004143615.1_repeats_filtered1.fasta Acropora_millepora_GCF_004143615.1_amil_sf_1.1_genomic.fna
cat Acropora_millepora_GCF_004143615.1_repeats_filtered1.fasta | /home/each/RepeatScout-1/filter-stage-2.prl --cat=Acropora_millepora_GCF_004143615.1_amil_sf_1.1_genomic.fna.out > Acropora_millepora_repeats_filtered2.fasta
/home/each/RepeatMasker/RepeatMasker -xm -xsmall -gff -pa 14 -s -lib Acropora_millepora_repeats_filtered2.fasta -nolow -norna -no_is -gff Acropora_millepora_GCF_004143615.1_amil_sf_1.1_genomic.fna

/home/each/RepeatScout-1/build_lmer_table -sequence Actinia_tenebrosa_GCF_009602425.1_ASM960242v1_genomic.fna -freq Actinia_tenebrosa_GCF_009602425.1_output_lmer.frequency
/home/each/RepeatScout-1/RepeatScout -sequence Actinia_tenebrosa_GCF_009602425.1_ASM960242v1_genomic.fna -output Actinia_tenebrosa_GCF_009602425.1_repeats.fasta -freq Actinia_tenebrosa_GCF_009602425.1_output_lmer.frequency
cat Actinia_tenebrosa_GCF_009602425.1_repeats.fasta | /home/each/RepeatScout-1/filter-stage-1.prl > Actinia_tenebrosa_GCF_009602425.1_repeats_filtered1.fasta
/home/each/RepeatMasker/RepeatMasker -xm -xsmall -gff -s -pa 14 -lib Actinia_tenebrosa_GCF_009602425.1_repeats_filtered1.fasta Actinia_tenebrosa_GCF_009602425.1_ASM960242v1_genomic.fna
cat Actinia_tenebrosa_GCF_009602425.1_repeats_filtered1.fasta | /home/each/RepeatScout-1/filter-stage-2.prl --cat=Actinia_tenebrosa_GCF_009602425.1_ASM960242v1_genomic.fna.out > Actinia_tenebrosa_repeats_filtered2.fasta
/home/each/RepeatMasker/RepeatMasker -xm -xsmall -gff -pa 14 -s -lib Actinia_tenebrosa_repeats_filtered2.fasta -nolow -norna -no_is -gff Actinia_tenebrosa_GCF_009602425.1_ASM960242v1_genomic.fna

/home/each/RepeatScout-1/build_lmer_table -sequence Aiptasia_pallida_GCF_001417965.1_genome_1.1_genomic.fna -freq Aiptasia_pallida_GCF_001417965.1_output_lmer.frequency
/home/each/RepeatScout-1/RepeatScout -sequence Aiptasia_pallida_GCF_001417965.1_genome_1.1_genomic.fna -output Aiptasia_pallida_GCF_001417965.1_repeats.fasta -freq Aiptasia_pallida_GCF_001417965.1_output_lmer.frequency
cat Aiptasia_pallida_GCF_001417965.1_repeats.fasta | /home/each/RepeatScout-1/filter-stage-1.prl > Aiptasia_pallida_GCF_001417965.1_repeats_filtered1.fasta
/home/each/RepeatMasker/RepeatMasker -xm -xsmall -gff -s -pa 14 -lib Aiptasia_pallida_GCF_001417965.1_repeats_filtered1.fasta Aiptasia_pallida_GCF_001417965.1_genome_1.1_genomic.fna
cat Aiptasia_pallida_GCF_001417965.1_repeats_filtered1.fasta | /home/each/RepeatScout-1/filter-stage-2.prl --cat=Aiptasia_pallida_GCF_001417965.1_genome_1.1_genomic.fna.out > Aiptasia_pallida_repeats_filtered2.fasta
/home/each/RepeatMasker/RepeatMasker -xm -xsmall -gff -pa 14 -s -lib Aiptasia_pallida_repeats_filtered2.fasta -nolow -norna -no_is -gff Aiptasia_pallida_GCF_001417965.1_genome_1.1_genomic.fna

/home/each/RepeatScout-1/build_lmer_table -sequence GCF_004324835.1_Dendronephthya_gigantea_1.0_genomic.fna -freq GCF_004324835.1_Dendronephthya_gigantea_output_lmer.frequency
/home/each/RepeatScout-1/RepeatScout -sequence GCF_004324835.1_Dendronephthya_gigantea_1.0_genomic.fna -output Dendronephthya_gigantea_repeats.fasta -freq GCF_004324835.1_Dendronephthya_gigantea_output_lmer.frequency
cat Dendronephthya_gigantea_repeats.fasta | /home/each/RepeatScout-1/filter-stage-1.prl > Dendronephthya_gigantea_repeats_filtered1.fasta
/home/each/RepeatMasker/RepeatMasker -xm -xsmall -gff -s -pa 14 -lib Dendronephthya_gigantea_repeats_filtered1.fasta GCF_004324835.1_Dendronephthya_gigantea_1.0_genomic.fna
cat Dendronephthya_gigantea_repeats_filtered1.fasta | /home/each/RepeatScout-1/filter-stage-2.prl --cat=GCF_004324835.1_Dendronephthya_gigantea_1.0_genomic.fna.out > Dendronephthya_gigantea_repeats_filtered2.fasta
/home/each/RepeatMasker/RepeatMasker -xm -xsmall -gff -pa 14 -s -lib Dendronephthya_gigantea_repeats_filtered2.fasta -nolow -norna -no_is -gff GCF_004324835.1_Dendronephthya_gigantea_1.0_genomic.fna

/home/each/RepeatScout-1/build_lmer_table -sequence Henneguya_salminicola_GCA_009887335.1_TAU_Hsal_1_genomic.fna -freq Henneguya_salminicola_output_lmer.frequency
/home/each/RepeatScout-1/RepeatScout -sequence Henneguya_salminicola_GCA_009887335.1_TAU_Hsal_1_genomic.fna -output Henneguya_salminicola_repeats.fasta -freq Henneguya_salminicola_output_lmer.frequency
cat Henneguya_salminicola_repeats.fasta | /home/each/RepeatScout-1/filter-stage-1.prl > Henneguya_salminicola_repeats_filtered1.fasta
/home/each/RepeatMasker/RepeatMasker -xm -xsmall -gff -s -pa 14 -lib Henneguya_salminicola_repeats_filtered1.fasta Henneguya_salminicola_GCA_009887335.1_TAU_Hsal_1_genomic.fna
cat Henneguya_salminicola_repeats_filtered1.fasta | /home/each/RepeatScout-1/filter-stage-2.prl --cat=Henneguya_salminicola_GCA_009887335.1_TAU_Hsal_1_genomic.fna.out > Henneguya_salminicola_repeats_filtered2.fasta
/home/each/RepeatMasker/RepeatMasker -xm -xsmall -gff -pa 14 -s -lib Henneguya_salminicola_repeats_filtered2.fasta -nolow -norna -no_is -gff Henneguya_salminicola_GCA_009887335.1_TAU_Hsal_1_genomic.fna

/home/each/RepeatScout-1/build_lmer_table -sequence Hydra_vulgaris_GCF_000004095.1_unplaced.scaf.fna -freq Hydra_vulgaris_output_lmer.frequency
/home/each/RepeatScout-1/RepeatScout -sequence Hydra_vulgaris_GCF_000004095.1_unplaced.scaf.fna -output Hydra_vulgaris_repeats.fasta -freq Hydra_vulgaris_output_lmer.frequency
cat Hydra_vulgaris_repeats.fasta | /home/each/RepeatScout-1/filter-stage-1.prl > Hydra_vulgaris_repeats_filtered1.fasta
/home/each/RepeatMasker/RepeatMasker -xm -xsmall -gff -s -pa 14 -lib Hydra_vulgaris_repeats_filtered1.fasta Hydra_vulgaris_GCF_000004095.1_unplaced.scaf.fna
cat Hydra_vulgaris_repeats_filtered1.fasta | /home/each/RepeatScout-1/filter-stage-2.prl --cat=Hydra_vulgaris_GCF_000004095.1_unplaced.scaf.fna.out > Hydra_vulgaris_repeats_filtered2.fasta
/home/each/RepeatMasker/RepeatMasker -xm -xsmall -gff -pa 14 -s -lib Hydra_vulgaris_repeats_filtered2.fasta -nolow -norna -no_is -gff Hydra_vulgaris_GCF_000004095.1_unplaced.scaf.fna

/home/each/RepeatScout-1/build_lmer_table -sequence Myxobolus_squamalis_GCA_010108815.1_TAU_Msqu_1_genomic.fna -freq Myxobolus_squamalis_output_lmer.frequency
/home/each/RepeatScout-1/RepeatScout -sequence Myxobolus_squamalis_GCA_010108815.1_TAU_Msqu_1_genomic.fna -output Myxobolus_squamalis_repeats.fasta -freq Myxobolus_squamalis_output_lmer.frequency
cat Myxobolus_squamalis_repeats.fasta | /home/each/RepeatScout-1/filter-stage-1.prl > Myxobolus_squamalis_repeats_filtered1.fasta
/home/each/RepeatMasker/RepeatMasker -xm -xsmall -gff -s -pa 14 -lib Myxobolus_squamalis_repeats_filtered1.fasta Myxobolus_squamalis_GCA_010108815.1_TAU_Msqu_1_genomic.fna
cat Myxobolus_squamalis_repeats_filtered1.fasta | /home/each/RepeatScout-1/filter-stage-2.prl --cat=Myxobolus_squamalis_GCA_010108815.1_TAU_Msqu_1_genomic.fna.out > Myxobolus_squamalis_repeats_filtered2.fasta
/home/each/RepeatMasker/RepeatMasker -xm -xsmall -gff -pa 14 -s -lib Myxobolus_squamalis_repeats_filtered2.fasta -nolow -norna -no_is -gff Myxobolus_squamalis_GCA_010108815.1_TAU_Msqu_1_genomic.fna

/home/each/RepeatScout-1/build_lmer_table -sequence Orbicella_faveolata_GCF_002042975.1_ofav_dov_v1_genomic.fna -freq Orbicella_faveolata_output_lmer.frequency
/home/each/RepeatScout-1/RepeatScout -sequence Orbicella_faveolata_GCF_002042975.1_ofav_dov_v1_genomic.fna -output Orbicella_faveolata_repeats.fasta -freq Orbicella_faveolata_output_lmer.frequency
cat Orbicella_faveolata_repeats.fasta | /home/each/RepeatScout-1/filter-stage-1.prl > Orbicella_faveolata_repeats_filtered1.fasta
/home/each/RepeatMasker/RepeatMasker -xm -xsmall -gff -s -pa 14 -lib Orbicella_faveolata_repeats_filtered1.fasta Orbicella_faveolata_GCF_002042975.1_ofav_dov_v1_genomic.fna
cat Orbicella_faveolata_repeats_filtered1.fasta | /home/each/RepeatScout-1/filter-stage-2.prl --cat=Orbicella_faveolata_GCF_002042975.1_ofav_dov_v1_genomic.fna.out > Orbicella_faveolata_repeats_filtered2.fasta
/home/each/RepeatMasker/RepeatMasker -xm -xsmall -gff -pa 14 -s -lib Orbicella_faveolata_repeats_filtered2.fasta -nolow -norna -no_is -gff Orbicella_faveolata_GCF_002042975.1_ofav_dov_v1_genomic.fna

/home/each/RepeatScout-1/build_lmer_table -sequence Morbakka_virulenta_GCA_003991215.1_MVIv1_genomic.fna -freq Morbakka_virulenta_output_lmer.frequency
/home/each/RepeatScout-1/RepeatScout -sequence Morbakka_virulenta_GCA_003991215.1_MVIv1_genomic.fna -output Morbakka_virulenta_repeats.fasta -freq Morbakka_virulenta_output_lmer.frequency
cat Morbakka_virulenta_repeats.fasta | /home/each/RepeatScout-1/filter-stage-1.prl > Morbakka_virulenta_repeats_filtered1.fasta
/home/each/RepeatMasker/RepeatMasker -xm -xsmall -gff -s -pa 14 -lib Morbakka_virulenta_repeats_filtered1.fasta Morbakka_virulenta_GCA_003991215.1_MVIv1_genomic.fna
cat Morbakka_virulenta_repeats_filtered1.fasta | /home/each/RepeatScout-1/filter-stage-2.prl --cat=Morbakka_virulenta_GCA_003991215.1_MVIv1_genomic.fna.out > Morbakka_virulenta_repeats_filtered2.fasta
/home/each/RepeatMasker/RepeatMasker -xm -xsmall -gff -pa 14 -s -lib Morbakka_virulenta_repeats_filtered2.fasta -nolow -norna -no_is -gff Morbakka_virulenta_GCA_003991215.1_MVIv1_genomic.fna

/home/each/RepeatScout-1/build_lmer_table -sequence Aurelia.Genome_v1.2_11-27-18.fasta -freq Aurelia_coerulia_output_lmer.frequency
/home/each/RepeatScout-1/RepeatScout -sequence Aurelia.Genome_v1.2_11-27-18.fasta -output Aurelia_coerulia_repeats.fasta -freq Aurelia_coerulia_output_lmer.frequency
cat Aurelia_coerulia_repeats.fasta | /home/each/RepeatScout-1/filter-stage-1.prl > Aurelia_coerulia_repeats_filtered1.fasta
/home/each/RepeatMasker/RepeatMasker -xm -xsmall -gff -s -pa 14 -lib Aurelia_coerulia_repeats_filtered1.fasta Aurelia.Genome_v1.2_11-27-18.fasta
cat Aurelia_coerulia_repeats_filtered1.fasta | /home/each/RepeatScout-1/filter-stage-2.prl --cat=Aurelia.Genome_v1.2_11-27-18.fasta.out > Aurelia_coerulia_repeats_filtered2.fasta
/home/each/RepeatMasker/RepeatMasker -xm -xsmall -gff -pa 14 -s -lib Aurelia_coerulia_repeats_filtered2.fasta -nolow -norna -no_is -gff Aurelia.Genome_v1.2_11-27-18.fasta

/home/each/RepeatScout-1/build_lmer_table -sequence Pocillopora_damicornis_GCF_003704095.1_ASM370409v1_genomic.fna -freq Pocillopora_damicornis_output_lmer.frequency
/home/each/RepeatScout-1/RepeatScout -sequence Pocillopora_damicornis_GCF_003704095.1_ASM370409v1_genomic.fna -output Pocillopora_damicornis_repeats.fasta -freq Pocillopora_damicornis_output_lmer.frequency
cat Pocillopora_damicornis_repeats.fasta | /home/each/RepeatScout-1/filter-stage-1.prl > Pocillopora_damicornis_repeats_filtered1.fasta
/home/each/RepeatMasker/RepeatMasker -xm -xsmall -gff -s -pa 14 -lib Pocillopora_damicornis_repeats_filtered1.fasta Pocillopora_damicornis_GCF_003704095.1_ASM370409v1_genomic.fna
cat Pocillopora_damicornis_repeats_filtered1.fasta | /home/each/RepeatScout-1/filter-stage-2.prl --cat=Pocillopora_damicornis_GCF_003704095.1_ASM370409v1_genomic.fna.out > Pocillopora_damicornis_repeats_filtered2.fasta
/home/each/RepeatMasker/RepeatMasker -xm -xsmall -gff -pa 14 -s -lib Pocillopora_damicornis_repeats_filtered2.fasta -nolow -norna -no_is -gff Pocillopora_damicornis_GCF_003704095.1_ASM370409v1_genomic.fna
bedtools subtract -A -a Pocillopora_damicornis_GCF_003704095.1_ASM370409v1_genomic.fna.out.gff -b Pocillopora_damicornis_GCF_003704095.1_ASM370409v1_gene.gff > Pocillopora_damicornis_repeats_not_overlapping_genes.gff
bedtools getfasta -fi Pocillopora_damicornis_GCF_003704095.1_ASM370409v1_genomic.fna -bed Pocillopora_damicornis_repeats_not_overlapping_genes.gff -name -fo Pocillopora_damicornis_repeats_not_overlapping_genes.fasta


/home/each/RepeatScout-1/build_lmer_table -sequence Thelohanellus_kitauei_GCA_000827895.1_ASM82789v1_genomic.fna -freq Thelohanellus_kitauei_output_lmer.frequency
/home/each/RepeatScout-1/RepeatScout -sequence Thelohanellus_kitauei_GCA_000827895.1_ASM82789v1_genomic.fna -output Thelohanellus_kitauei_repeats.fasta -freq Thelohanellus_kitauei_output_lmer.frequency
cat Thelohanellus_kitauei_repeats.fasta | /home/each/RepeatScout-1/filter-stage-1.prl > Thelohanellus_kitauei_repeats_filtered1.fasta
/home/each/RepeatMasker/RepeatMasker -xm -xsmall -gff -s -pa 14 -lib Thelohanellus_kitauei_repeats_filtered1.fasta Thelohanellus_kitauei_GCA_000827895.1_ASM82789v1_genomic.fna
cat Thelohanellus_kitauei_repeats_filtered1.fasta | /home/each/RepeatScout-1/filter-stage-2.prl --cat=Thelohanellus_kitauei_GCA_000827895.1_ASM82789v1_genomic.fna.out > Thelohanellus_kitauei_repeats_filtered2.fasta
/home/each/RepeatMasker/RepeatMasker -xm -xsmall -gff -pa 14 -s -lib Thelohanellus_kitauei_repeats_filtered2.fasta -nolow -norna -no_is -gff Thelohanellus_kitauei_GCA_000827895.1_ASM82789v1_genomic.fna
bedtools subtract -A -a Thelohanellus_kitauei_GCA_000827895.1_ASM82789v1_genomic.fna.out.gff -b Thelohanellus_kitauei_GCA_000827895.1_ASM82789v1_gene.gff > Thelohanellus_kitauei_repeats_not_overlapping_genes.gff
bedtools getfasta -fi Thelohanellus_kitauei_GCA_000827895.1_ASM82789v1_genomic.fna -bed Thelohanellus_kitauei_repeats_not_overlapping_genes.gff -name -fo Thelohanellus_kitauei_repeats_not_overlapping_genes.fasta

/home/each/RepeatScout-1/build_lmer_table -sequence GCF_002571385.1_Stylophora_pistillata_v1_genomic.fna -freq Stylophora_pistillata_output_lmer.frequency
/home/each/RepeatScout-1/RepeatScout -sequence GCF_002571385.1_Stylophora_pistillata_v1_genomic.fna -output Stylophora_pistillata_repeats.fasta -freq Stylophora_pistillata_output_lmer.frequency
cat Stylophora_pistillata_repeats.fasta | /home/each/RepeatScout-1/filter-stage-1.prl > Stylophora_pistillata_repeats_filtered1.fasta
/home/each/RepeatMasker/RepeatMasker -xm -xsmall -gff -s -pa 14 -lib Stylophora_pistillata_repeats_filtered1.fasta GCF_002571385.1_Stylophora_pistillata_v1_genomic.fna
cat Stylophora_pistillata_repeats_filtered1.fasta | /home/each/RepeatScout-1/filter-stage-2.prl --cat=GCF_002571385.1_Stylophora_pistillata_v1_genomic.fna.out > Stylophora_pistillata_repeats_filtered2.fasta
/home/each/RepeatMasker/RepeatMasker -xm -xsmall -gff -pa 14 -s -lib Stylophora_pistillata_repeats_filtered2.fasta -nolow -norna -no_is -gff GCF_002571385.1_Stylophora_pistillata_v1_genomic.fna
bedtools subtract -A -a GCF_002571385.1_Stylophora_pistillata_v1_genomic.fna.out.gff -b GCF_002571385.1_Stylophora_pistillata_v1_gene.gff > Stylophora_pistillata_repeats_not_overlapping_genes.gff
bedtools getfasta -fi GCF_002571385.1_Stylophora_pistillata_v1_genomic.fna -bed Stylophora_pistillata_repeats_not_overlapping_genes.gff -name -fo Stylophora_pistillata_repeats_not_overlapping_genes.fasta

/home/each/RepeatScout-1/build_lmer_table -sequence Nematostella_vectensis_GCF_000209225.1_ASM20922v1_genomic.fna -freq Nematostella_vectensis_output_lmer.frequency
/home/each/RepeatScout-1/RepeatScout -sequence Nematostella_vectensis_GCF_000209225.1_ASM20922v1_genomic.fna -output Nematostella_vectensis_repeats.fasta -freq Nematostella_vectensis_output_lmer.frequency
cat Nematostella_vectensis_repeats.fasta | /home/each/RepeatScout-1/filter-stage-1.prl > Nematostella_vectensis_repeats_filtered1.fasta
/home/each/RepeatMasker/RepeatMasker -xm -xsmall -gff -s -pa 14 -lib Nematostella_vectensis_repeats_filtered1.fasta Nematostella_vectensis_GCF_000209225.1_ASM20922v1_genomic.fna
cat Nematostella_vectensis_repeats_filtered1.fasta | /home/each/RepeatScout-1/filter-stage-2.prl --cat=Nematostella_vectensis_GCF_000209225.1_ASM20922v1_genomic.fna.out > Nematostella_vectensis_repeats_filtered2.fasta
/home/each/RepeatMasker/RepeatMasker -xm -xsmall -gff -pa 14 -s -lib Nematostella_vectensis_repeats_filtered2.fasta -nolow -norna -no_is -gff Nematostella_vectensis_GCF_000209225.1_ASM20922v1_genomic.fna
bedtools subtract -A -a Nematostella_vectensis_GCF_000209225.1_ASM20922v1_genomic.fna.out.gff -b Nematostella_vectensis_GCF_000209225.1_ASM20922v1_gene.gff > Nematostella_vectensis_repeats_not_overlapping_genes.gff
bedtools getfasta -fi Nematostella_vectensis_GCF_000209225.1_ASM20922v1_genomic.fna -bed Nematostella_vectensis_repeats_not_overlapping_genes.gff -name -fo Nematostella_vectensis_repeats_not_overlapping_genes.fasta





































