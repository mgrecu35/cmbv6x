import subprocess
import glob
import time
import sys
#pfiles = ['20140901-S204159-E221431.002897','20140902-S025211-E042443.002901','20140905-S123710-E140941.002954','20140909-S182446-E195717.003020','20140912-S002324-E015555.003055','20140916-S104736-E122006.003124','20140917-S234911-E012144.003148']
#pfiles = ['20140902-S025211-E042443.002901','20140905-S123710-E140941.002954','20140909-S182446-E195717.003020','20140912-S002324-E015555.003055','20140916-S104736-E122006.003124','20140917-S234911-E012144.003148']
paramdir = 'Param_Files'
pfiles = glob.glob(paramdir + '/*')
pfiles.sort()
print pfiles

nodes=['178', '179', '181', '183']

i=0
p=[]
for file in pfiles:
	if file.find("output") > 0: continue
	print 'found file: ' + file
#	if((int(file[-6:]) < 3684) | (int(file[-6:]) > 3716)): continue

	#if(file[31:33] != '04'): continue
	i=i+1
	#print i
	#if(i < 547): continue
	#check that node is available
      	node_inuse = True
      	while(node_inuse):
      		for node in nodes:
      	  		process = subprocess.Popen(['nodeps.sh'],stdin=None,stdout=subprocess.PIPE,stderr=None,shell=True)
      	  		node_inuse=False
      	  		for line in process.stdout:
      	    			#print line
	    			if(line[0:3] == node): 
	      				node_inuse = True
	      				break
	      
      	  		print node, node_inuse
      	  		if(node_inuse == False): break
      
       		if(node_inuse): 
      	  		print 'Waiting for available node'
			for child in p:
				print child.poll()
      	  		time.sleep(60)
  
	print 'Using node '+node
	stdoutputfile = file + '.output'
	fout=open(stdoutputfile, mode = 'w')
	p.append(subprocess.Popen(['bpsh',node,'./combAlg.v5.2.exe', 'junk', file], stdin=None, stdout=fout,stderr=fout,shell=False))
	fout.close()
	#p.append(subprocess.Popen(['bpsh',node,'./combAlg_test.exe', 'junk', 'paramFiles/paramFile_sjm.'+file], stdin=None,stdout=None,stderr=None,shell=False))
