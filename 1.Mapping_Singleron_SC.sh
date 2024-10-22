#Using CeleScope software to get Gene Expression Matrix
git clone https://github.com/singleron-RD/CeleScope.git

cd CeleScope
conda create -n celescope
conda activate celescope
conda install -y --file conda_pkgs.txt
pip install celescope

cd ref
conda activate celescope
celescope rna mkref \
 --genome_name Sus_scrofa.Sscrofa11.1 \
 --fasta Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
 --gtf Sus_scrofa.Sscrofa11.1.101.gtf

du -sh LC*|awk -v path=`pwd` '{print $2, " ", path"/"$2, " ", $2}' > ./map_file

multi_rna --mapfile ./map_file --genomeDir ./ref --thread 8 --mod shell
