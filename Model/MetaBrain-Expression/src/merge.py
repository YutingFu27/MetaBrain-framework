import os
files=os.listdir('./')
import re
matched_files=[file for file in files if file.endswith('.npy')]
print(matched_files)
matched_files=[file for file in matched_files if file not in ['chrX_chromatin.npy','chrY_chromatin.npy']]
def numerical_sort(value):
    number=re.findall(r'\d+',value)
    if number:
        return int(number[0])
    else:
        return value
matched_files=sorted(matched_files,key=numerical_sort)
matched_files=matched_files+['chrX_chromatin.npy','chrY_chromatin.npy']
import numpy as np
Xreducedall_list=[]
for file in matched_files:
    print(file)
    x=np.load(file)
    print(x.shape)
    Xreducedall_list.append(x)

Xreducedall_all=np.concatenate(Xreducedall_list,axis=0)
Xreducedall_all =Xreducedall_all.reshape(24339,6103*10)
print('done')
print(Xreducedall_all.shape)
np.save('Xreducedall_all.npy',Xreducedall_all)

    
