# -- coding: utf-8 --
import os


ddsmdir = "/data/DDSM/out_16/out_benign/"
LMLO_dir = "/home/zju/rzx/dm/normalization/target/000138-LMLO.png"
LCC_dir = "/home/zju/rzx/dm/normalization/target/100151-LCC.png"
RMLO_dir = "/home/zju/rzx/dm/normalization/target/105552-RMLO.png"
RCC_dir = "/home/zju/rzx/dm/normalization/target/105551-RCC.png"


for root, subFolders, file_names in os.walk(ddsmdir):
    for file_name in file_names:
        if ".LEFT_MLO" in file_name:
            lmlo_path = os.path.join(root, file_name)
            cmd = './normalization '+lmlo_path+' '+LMLO_dir
            os.system(cmd)

        if ".LEFT_CC" in file_name:
            lcc_path = os.path.join(root, file_name)
            cmd = './normalization '+lcc_path+' '+LCC_dir
            os.system(cmd)

        if ".RIGHT_MLO" in file_name:
            rmlo_path = os.path.join(root, file_name)
            cmd = './normalization '+rmlo_path+' '+RMLO_dir
            os.system(cmd)

        if ".RIGHT_CC" in file_name:
            rcc_path = os.path.join(root, file_name)
            cmd = './normalization '+rcc_path+' '+RCC_dir
            os.system(cmd)

print('normalization cancer done')