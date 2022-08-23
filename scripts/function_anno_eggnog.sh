#!/bin/bash
emapper.py  --target_orthologs all --go_evidence non-electronic -m diamond --pident 40 --seed_ortholog_evalue 0.001 --seed_ortholog_score 60 --query_cover 20 --subject_cover 20 --override -i Athaliana_167_protein.fa -o Arabidopsis_thaliana.protein.fa.eggnog --output_dir eggnog/ --temp_dir eggnog_tmp/ --database none --cpu 138 --target_taxa 3699

emapper.py  --target_orthologs all --go_evidence non-electronic -m diamond --pident 40 --seed_ortholog_evalue 0.001 --seed_ortholog_score 60 --query_cover 20 --subject_cover 20 --override -i Aspergillus_fumigatus.protein.fa -o Aspergillus_fumigatus.protein.fa.eggnog --output_dir eggnog/ --temp_dir eggnog_tmp/ --database none --cpu 138 --target_taxa 4751

emapper.py  --target_orthologs all --go_evidence non-electronic -m diamond --pident 40 --seed_ortholog_evalue 0.001 --seed_ortholog_score 60 --query_cover 20 --subject_cover 20 --override -i Coprinopsis_cinerea.proteins.fa -o Coprinopsis_cinerea.proteins.fa.eggnog --output_dir eggnog/ --temp_dir eggnog_tmp/ --database none --cpu 138 --target_taxa 4751

emapper.py  --target_orthologs all --go_evidence non-electronic -m diamond --pident 40 --seed_ortholog_evalue 0.001 --seed_ortholog_score 60 --query_cover 20 --subject_cover 20 --override -i Danio_rerio.proteins.fa -o Danio_rerio.proteins.fa.eggnog --output_dir eggnog/ --temp_dir eggnog_tmp/ --database none --cpu 138 --target_taxa 7898

emapper.py  --target_orthologs all --go_evidence non-electronic -m diamond --pident 40 --seed_ortholog_evalue 0.001 --seed_ortholog_score 60 --query_cover 20 --subject_cover 20 --override -i Drosophila_melanogaster.proteins.fa -o Drosophila_melanogaster.proteins.fa.eggnog --output_dir eggnog/ --temp_dir eggnog_tmp/ --database none --cpu 138 --target_taxa 7214

emapper.py  --target_orthologs all --go_evidence non-electronic -m diamond --pident 40 --seed_ortholog_evalue 0.001 --seed_ortholog_score 60 --query_cover 20 --subject_cover 20 --override -i Fusarium_graminearum.protein.fa -o Fusarium_graminearum.protein.fa.eggnog --output_dir eggnog/ --temp_dir eggnog_tmp/ --database none --cpu 138 --target_taxa 4751

emapper.py  --target_orthologs all --go_evidence non-electronic -m diamond --pident 40 --seed_ortholog_evalue 0.001 --seed_ortholog_score 60 --query_cover 20 --subject_cover 20 --override -i Neurospora_crassa.protein.fa -o Neurospora_crassa.protein.fa.eggnog --output_dir eggnog/ --temp_dir eggnog_tmp/ --database none --cpu 138 --target_taxa 4751

emapper.py  --target_orthologs all --go_evidence non-electronic -m diamond --pident 40 --seed_ortholog_evalue 0.001 --seed_ortholog_score 60 --query_cover 20 --subject_cover 20 --override -i Rhizopus_delemar.protein.fa -o Rhizopus_delemar.protein.fa.eggnog --output_dir eggnog/ --temp_dir eggnog_tmp/ --database none --cpu 138 --target_taxa 4751
