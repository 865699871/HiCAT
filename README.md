[![DOI](https://zenodo.org/badge/558735701.svg)](https://zenodo.org/badge/latestdoi/558735701)
# HiCAT: Hierarchical Centromere structure AnnoTation Tool

Advanced long-read sequencing technologies have revolutionized genome assembly, unlocking the complex region centromere and signaling the new stage in genomics research. The new computing problems generated by these new areas, like centromere annotation problem, required novel bioinformatics methods. Here, we proposed HiCAT, a generalized computational tool based on hierarchical tandem repeat mining (HTRM) method to automatically process centromere annotation.
## Dependencies
Python 3.9.13

Packages  | Version |
--------- | --------|
biopython  | 1.79 |
setuptools  | 61.2.0 |
joblib  | 1.1.0 |
numpy  | 1.22.3 |
pandas  | 1.4.0 |
python-levenshtein  | 0.12.2 |
python-edlib  | 1.3.9 |
networkx  | 2.7.1 |
matplotlib  | 3.5.1 |

StringDecomposer (https://github.com/ablab/stringdecomposer) with version 1.1.2.   
Development environment: Linux  
Development tool: Pycharm  

## Installation and Quick start
#### Source code (g++ version 5.3.1 or higher for stringdecomposer)
```Bash
#install
conda install -y --file requirements.txt
cd ./stringdecomposer && make
#run
python HiCAT.py -i ./testdata/cen21.fa -t ./testdata/AlphaSat.fa
#For more details, please use
python HiCAT.py -h
```
#### Conda 
```Bash
#install
conda install -c shenghangao hicat
#run
hicat -i ./testdata/cen21.fa -t ./testdata/AlphaSat.fa
#For more details, please use
hicat -h
```
### Visualization
HiCAT default visualized the top five HORs with repeat numbers greater than 10 in maximum HOR coverage similarity. 

```Bash
Custom visualization can use visualization.py
-r is HiCAT result directory.
-s is which similarity be visualized.
-sp is the number of top HORs.
-sn is the minimum repeat number of HOR.
```
### Output
The result is in the directory named out.

out_hor.raw.fa contains HOR DNA sequences. 

out_hor.normal.fa contains each normalized HOR DNA sequence. We normalized the raw DNA sequence to one represent HOR. For example, normalized 9_10_(1_2_3_1_2_3)_4_5_6_7_4_8 to 1_2_3_4_5_6_7_4_8_9_10 in CEN21.

out_monomer.fa contains monomer DNA sequences.

pattern_static.xls is the HOR repeat number.

pattern_static.pdf is the HOR repeat number in bar plot.

plot_pattern.pdf is the plot for hierarchical centromere structure annotation

out_top_layer.xls is the annotation in top layer.

out_all_layer.xls is the annotation in all layer. Label "top" represent this region is in top layer. Label "cover" represent this region is covered by a top layer region.

## Contact
If you have any questions, please feel free to contact: gaoxian15002970749@163.com, xfyang@xjtu.edu.cn, kaiye@xjtu.edu.cn





