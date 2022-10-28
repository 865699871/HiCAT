conda create --name HiCAT python=3.9
source activate HiCAT
conda install -y --file requirements.txt
cd ./stringdecomposer
make
