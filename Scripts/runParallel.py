from subprocess import Popen, PIPE
import glob

f_list = glob.glob('paramX/*param')
f_list=sorted(f_list)[0:200]
nodes1=[i for i in range(236,246)]
nodes=[]
for i in range(20):
    nodes.extend(nodes1)

cmds_list = [['bpsh', str(i), './combAlg.exe', 'junk ',file_name] for (i,file_name) in zip(nodes,f_list)]


for i in range(20):
    cmds_list1=cmds_list[i*10:i*10+10]
    print(cmds_list1)
    procs_list = [Popen(cmd, stdout=PIPE, stderr=PIPE) for cmd in cmds_list1]
    for proc in procs_list:
        proc.wait()
        
