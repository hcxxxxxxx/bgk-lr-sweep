"""Cross-check uninstrumented archive main/solver against the audited run."""
from pathlib import Path
import json
import os
import re
import subprocess
import zipfile
import numpy as np

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/gray_reference_platform'
SRC=OUT/'source';NATIVE=OUT/'native_gsis';BUILD=OUT/'native_build'
NATIVE.mkdir(exist_ok=True);BUILD.mkdir(exist_ok=True)
(NATIVE/'res').mkdir(exist_ok=True);(NATIVE/'msh').mkdir(exist_ok=True)
with zipfile.ZipFile('/Users/hanchenxi/Desktop/GSIS_Gray_20261110.zip') as z:
    for name in ['Solvers.f90','Callaway_2D_2V_DG.f90']:
        (BUILD/name).write_bytes(z.read('GSIS_Gray_20261110/'+name))
    (NATIVE/'control.in').write_bytes(z.read('GSIS_Gray_20261110/control.in'))
    (NATIVE/'msh/A1_Nx11_Ny11.msh').write_bytes(z.read('GSIS_Gray_20261110/msh/A1_Nx11_Ny11.msh'))
compiler='/opt/homebrew/bin/gfortran'
with (OUT/'native_build.log').open('w') as log:
    for name in ['Solvers','Callaway_2D_2V_DG']:
        subprocess.run([compiler,'-O3','-fopenmp','-ffree-line-length-none',f'-J{BUILD}',
            f'-I{BUILD}',f'-I{SRC/"bin"}','-c',str(BUILD/f'{name}.f90'),'-o',str(BUILD/f'{name}.o')],
            cwd=BUILD,stdout=log,stderr=subprocess.STDOUT,check=True)
    objects=[p for p in sorted((SRC/'bin').glob('*.o')) if p.stem not in ['Solvers','Callaway_2D_2V_DG']]
    objects += [BUILD/'Solvers.o',BUILD/'Callaway_2D_2V_DG.o']
    subprocess.run([compiler,*map(str,objects),'-L/opt/homebrew/opt/suite-sparse/lib','-lumfpack',
        '-L/opt/homebrew/opt/openblas/lib','-lopenblas','-fopenmp','-o',str(BUILD/'DGACC_native')],
        stdout=log,stderr=subprocess.STDOUT,check=True)
env=dict(os.environ,OMP_THREAD_LIMIT='4',OPENBLAS_NUM_THREADS='1',GFORTRAN_UNBUFFERED_ALL='1')
with (NATIVE/'run.log').open('w') as log:
    subprocess.run([str(BUILD/'DGACC_native')],cwd=NATIVE,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)

def updates(path):
    rows=[]
    for line in path.read_text().splitlines():
        m=re.fullmatch(r'\s*(\d+)\s+([\d.Ee+-]+)\s+([\d.Ee+-]+)\s*',line)
        if m:rows.append([int(m[1]),float(m[2]),float(m[3])])
    return np.array(rows)

native=updates(NATIVE/'run.log');audit=updates(OUT/'gsis/run.log')[:len(native)]
assert native.shape==audit.shape and np.array_equal(native[:,0],audit[:,0])
error=float(np.max(np.abs(native[:,1:]-audit[:,1:])))
assert error<1e-12
result=dict(native_stop_step=int(native[-1,0]),native_final_macro_update=float(native[-1,1]),
    maximum_macro_update_or_mass_difference=error,read_only_audit_verified=True)
(OUT/'native_audit_check.json').write_text(json.dumps(result,indent=2)+'\n');print(result)
