###目的：使用新格元软件将测序数据进行比对###
###代码过程：单细胞组-姚天雄负责整理，完善。 -20221104 ###
###来源：https://github.com/singleron-RD/CeleScope ###

#1.获取文件
git clone https://github.com/singleron-RD/CeleScope.git

#2.创建环境，并安装软件
cd CeleScope
conda create -n celescope_1127
conda activate celescope_1127
conda install -y --file conda_pkgs.txt
pip install celescope

#3.创建参考基因组索引
cd ref
conda activate celescope_1127
celescope rna mkref \
 --genome_name Sus_scrofa.Sscrofa11.1 \
 --fasta Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
 --gtf Sus_scrofa.Sscrofa11.1.101.gtf


#4.举例创建样本表达矩阵
#准备my.mapfile
#格式如下:
# LC210124054 /Business/Jxnd_company/jxndyangbin/SeqData/20210405/LC210124054 LC210124054
# LC210124088 /Business/Jxnd_company/jxndyangbin/SeqData/20210405/LC210124088 LC210124088
# LC210124112 /Business/Jxnd_company/jxndyangbin/SeqData/20210405/LC210124112 LC210124112

#4.1创建（在数据样本目录下）map_file文件
du -sh LC*|awk -v path=`pwd` '{print $2, " ", path"/"$2, " ", $2}' > ../01.out/map_file_1127_14sample

#4.2利用mapfile创建运行脚本
multi_rna --mapfile ./map_file_1127_14sample --genomeDir /home/SCell/ref --thread 8 --mod shell

#5.运行脚本（记得到生成脚本的路径下运行）
conda activate celescope_1127
nohup bash LC210124054.sh > ./LC210124054.sh.out &
