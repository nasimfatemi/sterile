import os
import subprocess
import time
from subprocess import call
count=10
for i in range(0,1000):  
            outputfilename='psterile'+str(i)
            with open(outputfilename,'w') as dest:
                        test1="root -b -q '/home/fatemigh/projects/rpp-jillings/fatemigh/sterile/Likelihood/LikelihoodSterile_modified7.C+("+str(i)+")'"
                        print test1
                        dest.write("#!/bin/bash\n");
                        dest.write("#SBATCH --account=rpp-jillings\n")
                        dest.write("#SBATCH --time=00:40:00\n")
                        dest.write("#SBATCH --job-name="+str(outputfilename)+"\n")
                        dest.write("#SBATCH --error=/home/fatemigh/projects/rpp-jillings/fatemigh/sterile/Likelihood/%x-%j.err"+"\n")
                        #dest.write("#SBATCH --mem-per-cpu=2024M\n")
                        dest.write("source  /home/fatemigh/projects/rpp-jillings/fatemigh/env.sh\n")
                        dest.write(test1)
                        exe="chmod +x  "+str(outputfilename)
                        subprocess.Popen(exe,shell=True)
                        sub="sbatch  /home/fatemigh/projects/rpp-jillings/fatemigh/sterile/Likelihood/"+str(outputfilename)
                        print sub
                        time.sleep(5)
                        subprocess.Popen(sub,shell=True)
 
            

