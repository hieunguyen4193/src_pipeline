# Trong-Hieu Nguyen, last modified on 24.08.2022.
current_dir=$(pwd);

path_to_cellranger="/home/hieunguyen/CRC1382/src/src_pipeline/CellRanger/cellranger-7.1.0";
export PATH=${path_to_cellranger}:$PATH;
path_to_save_output="/media/hieunguyen/My Book/GS";
path_to_config_file="/home/hieunguyen/CRC1382/src/single_cell_pipeline_23052023/10k_PBMC_5pv2_nextgem_Chromium_Controller_Multiplex_config.csv";

cellranger multi --id=test_pbmc_10k --csv=${path_to_config_file} --localcores=25;

# rsync -avh --progress ${current_dir}/test_pbmc_10k ${path_to_save_output}/test_pbmc_10k;

# rm -rf ${current_dir}/test_pbmc_10k;

# EOF

