import os
cmd = 'f2py model_tools.f90 -h model_tools.pyf -m model_tools --overwrite-signature'
os.system(cmd)
#Create driver library
#cmd = 'f2py model_tools.f90 -lgfortran -lgomp -c model_tools.pyf --fcompiler=gnu95 --f90flags="-w -fopenmp"'
cmd = 'f2py model_tools.f90 -lgfortran -lgomp -c model_tools.pyf --fcompiler=gnu95 --f90flags="-w -fopenmp"'
#cmd = 'f2py model_tools.f90 -lgfortran -lgomp -c model_tools.pyf --fcompiler=gnu95 --f90flags="-w -fopenmp" -L/usr/lib64 -L/lib64'
#cmd = 'f2py model_tools.f90 -c model_tools.pyf'
os.system(cmd)
#Move to the previois directory
#os.system('mv model_tools.*.so ../model_tools.so')
os.system('mv model_tools*so ../model_tools.so')
#os.system('mv model_tools.so ../.')
