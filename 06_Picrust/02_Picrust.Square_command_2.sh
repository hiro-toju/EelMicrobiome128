wget https://github.com/picrust/picrust2/archive/v2.4.2.tar.gz
tar xvzf  v2.4.2.tar.gz
cd picrust2-2.4.2/

conda env create -f  picrust2-env.yaml
conda activate picrust2
pip install --editable .



conda activate picrust2

picrust2_pipeline.py -s All.ASVs.fasta -i square.th30.txt -o Square.th30 -p 5

cd Square.th30
pathway_pipeline.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -o pathways_out --coverage --intermediate minpath_working -p 5
add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC -o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz
add_descriptions.py -i KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO -o KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz
add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC  -o pathways_out/path_abun_unstrat_descrip.tsv.gz
add_descriptions.py -i pathways_out/path_cov_unstrat.tsv.gz -m METACYC  -o pathways_out/path_cov_unstrat_descrip.tsv.gz

