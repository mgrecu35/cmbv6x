import glob
import os
fs=glob.glob("src_f90/*f90")
fs2=glob.glob("*.o")

for f in []:#fs:
    f1=f[8:-3]
    if f1+'o' not in fs2:
        cmd='mv %s oldStuff\n'%f
        print(cmd)
        os.system(cmd)
        

fs=glob.glob("src_c/*.cc")
for f in fs:
    f1=f[6:-2]
    print(f1)
    if f1+'o' not in fs2:
        cmd='mv %s oldStuff\n'%f
        print(cmd)
        #os.system(cmd)
