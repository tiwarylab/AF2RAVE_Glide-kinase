import os

def sbatchpreamble(path,filename,jobname="prodrun",time=5):
    f=open(path+"/"+filename,"w+")
    f.write("#!/usr/bin/env bash\n#SBATCH --partition=gpu\n#SBATCH --gres=gpu:a100:1\n#SBATCH --time=%i:00:00\n#SBATCH --job-name=%s\n"%(time,jobname))
    f.write("module purge\nmodule load cuda\nconda activate openmm\n")
    return f

#def runequil(mutant,afid,predid):
#    f=sbatchpreamble(".","%s_%i_job.sh"%(mutant,predid),mutant+"_eqrun")
#    f.write("\n python /home/xg23//equilprotocol.py -file /home/xg23/scratch/kinase/AF2structs/DDR1_%s_%s/pred_%i.pdb"%(mutant,afid,predid+1))
#    f.close()
#    os.system("sbatch %s_%i_job.sh"%(mutant,predid))
#    
#def basicmd(mutant,equilid):
#    f=sbatchpreamble(".","%s_job.sh"%mutant,mutant+"_prodrun")
#    f.write("\n python /home/xg23/scratch/kinase/pyscripts/basicmd.py -source /home/xg23/scratch/kinase/equilruns/%s/equil_%i/equilNPT_0"%(mutant,equilid))
#    f.close()
#    os.system("sbatch %s_job.sh"%(mutant))

f=sbatchpreamble(".","job_head.sh","test")
f.close()
