# -- coding: utf-8 --
import os
import csv


ddsmdir = "/data/DDSM/out_16/out_benign/"
LMLO_dir = "/home/zju/rzx/dm/normalization/target/000138-LMLO.png"
LCC_dir = "/home/zju/rzx/dm/normalization/target/100151-LCC.png"
RMLO_dir = "/home/zju/rzx/dm/normalization/target/105552-RMLO.png"
RCC_dir = "/home/zju/rzx/dm/normalization/target/105551-RCC.png"
mammo_dir = "/data/DDSM/firstpeople/presentation/"


# ddsm-normalization
# for root, subFolders, file_names in os.walk(ddsmdir):
#     for file_name in file_names:
#         if ".LEFT_MLO" in file_name:
#             lmlo_path = os.path.join(root, file_name)
#             cmd = './normalization '+lmlo_path+' '+LMLO_dir
#             os.system(cmd)
#
#         if ".LEFT_CC" in file_name:
#             lcc_path = os.path.join(root, file_name)
#             cmd = './normalization '+lcc_path+' '+LCC_dir
#             os.system(cmd)
#
#         if ".RIGHT_MLO" in file_name:
#             rmlo_path = os.path.join(root, file_name)
#             cmd = './normalization '+rmlo_path+' '+RMLO_dir
#             os.system(cmd)
#
#         if ".RIGHT_CC" in file_name:
#             rcc_path = os.path.join(root, file_name)
#             cmd = './normalization '+rcc_path+' '+RCC_dir
#             os.system(cmd)
#
# print('normalization cancer done')

# firstpeople-normalization
counter = 0
view,laterality,filename = [],[],[]
data = []
with open('/home/zju/rzx/dm/metadata/images_crosswalk_mix3_cut.tsv', 'r') as file_crosswalk:
    reader_crosswalk = csv.reader(file_crosswalk, delimiter='\t')
    for row in reader_crosswalk:
        if counter == 0:
            counter += 1
            continue
        counter += 1
        view = row[3]
        laterality = row[4]
        filename = row[5]
        data.append((filename,view, laterality))
# print(X)
print(counter, " lines in images crosswalk file.")

for root, subFolders, file_names in os.walk(mammo_dir):
    for file_name in file_names:
        if ".dcm" in file_name:
            dcm = [x for x in data if file_name in x[0]]
            if dcm != []:
                if dcm[0][2] in 'L' and dcm[0][1] in 'MLO':
                    # print(dcm[0][0], 'LMLO')
                    lmlo_path = os.path.join(root, dcm[0][0])
                    cmd = './normalization ' + lmlo_path + ' ' + LMLO_dir
                    os.system(cmd)
                elif dcm[0][2] in 'L' and dcm[0][1] in 'CC':
                    # print(dcm[0][0],'LCC')
                    lcc_path = os.path.join(root, dcm[0][0])
                    cmd = './normalization ' + lcc_path + ' ' + LCC_dir
                    os.system(cmd)
                elif dcm[0][2] in 'R' and dcm[0][1] in 'MLO':
                    # print(dcm[0][0],'RMLO')
                    rmlo_path = os.path.join(root, dcm[0][0])
                    cmd = './normalization ' + rmlo_path + ' ' + RMLO_dir
                    os.system(cmd)
                elif dcm[0][2] in 'R' and dcm[0][1] in 'CC':
                    # print(dcm[0][0],'RCC')
                    rcc_path = os.path.join(root, dcm[0][0])
                    cmd = './normalization ' + rcc_path + ' ' + RCC_dir
                    os.system(cmd)
                else:
                    print(dcm[0][0])

print('normalization cancer done')