"""Run latest uploaded Gray code with read-only original/coupled residual audits.

The shipped model, mesh, P3 basis, 80x40 angular rule, and tau values remain.
Only diagnostic code and the minimum number of iterations are instrumented.
Uses the archive's own macOS PARDISO/UMFPACK shim and native transport solver.
"""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import argparse
import difflib
import hashlib
import json
import os
import re
import shutil
import subprocess
import time
import zipfile

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = Path('/Users/hanchenxi/Desktop/GSIS_Gray_20261110.zip')
OUTPUT = ROOT/'results/gray_reference_platform'
FILES = ['Global_Parameter','USD_Math','Boundary_Conditions','Spatial_Mesh',
         'Velocity_Mesh','Matrix','Basis_Function','Integration','Velocity_Distribution',
         'Solvers','Synthetic_Acceleration','Out_Put_Result','PARDISO_Interface',
         'Callaway_2D_2V_DG']


def build():
    source = OUTPUT/'source'; source.mkdir(parents=True, exist_ok=True)
    originals = {}
    with zipfile.ZipFile(ARCHIVE) as z:
        for name in [f'{f}.f90' for f in FILES]+['Makefile','control.in']:
            original = z.read(f'GSIS_Gray_20261110/{name}')
            (source/name).write_bytes(original); originals[name] = original
        mesh = z.read('GSIS_Gray_20261110/msh/A1_Nx11_Ny11.msh')
    solver = (source/'Solvers.f90').read_text()
    solver = solver.replace('SUBROUTINE Calculate_Discrete_Transport_Residual(RESIDUAL)',
        'SUBROUTINE Calculate_Discrete_Transport_Residual(RESIDUAL,ABSOLUTE_RESIDUAL,REFERENCE_NORM)',1)
    marker = 'REAL(KIND=DBL), INTENT(OUT) :: RESIDUAL'
    assert solver.count(marker) == 1
    solver = solver.replace(marker,marker+'\nREAL(KIND=DBL), OPTIONAL, INTENT(OUT) :: ABSOLUTE_RESIDUAL,REFERENCE_NORM',1)
    marker = 'RESIDUAL=SQRT(NUMERATOR/MAX(DENOMINATOR,TINY(1.d0)))'
    solver = solver.replace(marker,marker+'\nIF(PRESENT(ABSOLUTE_RESIDUAL)) ABSOLUTE_RESIDUAL=SQRT(NUMERATOR)\n'
        'IF(PRESENT(REFERENCE_NORM)) REFERENCE_NORM=SQRT(DENOMINATOR)',1)
    (source/'Solvers.f90').write_text(solver)
    main = (source/'Callaway_2D_2V_DG.f90').read_text()
    main = main.replace('INTEGER :: ios, STEP, INNER_SWEEP',
        'INTEGER :: ios, STEP, INNER_SWEEP,AUDIT_MINIMUM,AUDIT_EVERY,AUDIT_FIRST_STOP\nCHARACTER(LEN=32) :: AUDIT_ENV')
    marker = 'STEP = 1\nDO'
    assert marker in main
    main = main.replace(marker,"""AUDIT_MINIMUM=0;AUDIT_EVERY=10;AUDIT_FIRST_STOP=-1
CALL GET_ENVIRONMENT_VARIABLE('GRAY_AUDIT_MINIMUM',AUDIT_ENV)
IF(LEN_TRIM(AUDIT_ENV).GT.0) READ(AUDIT_ENV,*)AUDIT_MINIMUM
CALL GET_ENVIRONMENT_VARIABLE('GRAY_AUDIT_EVERY',AUDIT_ENV)
IF(LEN_TRIM(AUDIT_ENV).GT.0) READ(AUDIT_ENV,*)AUDIT_EVERY
OPEN(96,FILE='equation_residual.csv',STATUS='REPLACE',ACTION='WRITE')
WRITE(96,'(A)')'step,macro_update,absolute_original,absolute_coupled,rhs_relative_original,rhs_relative_coupled,'// &
 'initial_relative_original,initial_relative_coupled,relative_source_gap,relative_temperature_gap,common_initial_absolute'
CALL Audit_Record(0,0.d0)
STEP = 1
DO""",1)
    marker = 'IF (RESIDUL.LT.TOL) EXIT'
    assert marker in main
    main = main.replace(marker,"""IF(RESIDUL.LT.TOL .AND. AUDIT_FIRST_STOP.LT.0)THEN
 AUDIT_FIRST_STOP=STEP
 WRITE(*,*)'NATIVE_FIRST_STOP',STEP,RESIDUL
END IF
IF(STEP.EQ.1 .OR. MOD(STEP,AUDIT_EVERY).EQ.0 .OR. STEP.EQ.AUDIT_FIRST_STOP .OR. &
 STEP.EQ.TMAX .OR. (RESIDUL.LT.TOL .AND. STEP.GE.AUDIT_MINIMUM)) CALL Audit_Record(STEP,RESIDUL)
IF (RESIDUL.LT.TOL .AND. STEP.GE.AUDIT_MINIMUM) EXIT""",1)
    main = main.replace('IF (ACCFLAG.EQ.1) CALL Release_PARDISO ()',
        'CLOSE(96)\nIF (ACCFLAG.EQ.1) CALL Release_PARDISO ()',1)
    include = (ROOT/'debug/gray_reference_audit.inc').read_text()
    main = main.replace('END PROGRAM Callaway_2D_2V_DG','CONTAINS\n'+include+'\nEND PROGRAM Callaway_2D_2V_DG')
    (source/'Callaway_2D_2V_DG.f90').write_text(main)
    diff=[]
    for name,original in originals.items():
        actual=(source/name).read_bytes()
        if actual!=original:
            diff.extend(difflib.unified_diff(original.decode().splitlines(True),actual.decode().splitlines(True),
                fromfile=f'archive/{name}',tofile=f'audit/{name}'))
    (OUTPUT/'diagnostic_only.patch').write_text(''.join(diff))
    manifest = dict(archive_sha256=hashlib.sha256(ARCHIVE.read_bytes()).hexdigest(),
        original_source_sha256={n:hashlib.sha256(b).hexdigest() for n,b in originals.items()},
        unchanged_transport=solver.split('SUBROUTINE DG_Solver_VDF ()')[1].split('END SUBROUTINE DG_Solver_VDF')[0]
            == originals['Solvers.f90'].decode().split('SUBROUTINE DG_Solver_VDF ()')[1].split('END SUBROUTINE DG_Solver_VDF')[0],
        unchanged_synthetic_acceleration=(source/'Synthetic_Acceleration.f90').read_bytes()==originals['Synthetic_Acceleration.f90'],
        runtime=dict(OMP_THREAD_LIMIT=4,OPENBLAS_NUM_THREADS=1),mesh_sha256=hashlib.sha256(mesh).hexdigest())
    (OUTPUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    args=['make','FCOMP=/opt/homebrew/bin/gfortran',
          'FLAG=-O3 -fopenmp -ffree-line-length-none -J./bin/ -I./bin/',
          'LIB_LIST=-L/opt/homebrew/opt/suite-sparse/lib -lumfpack -L/opt/homebrew/opt/openblas/lib -lopenblas -fopenmp',
          'PARDISO_OBJ=./bin/PARDISO_Interface.o']
    # Force a complete rebuild when diagnostic code changes within one second.
    if (source/'bin').exists():shutil.rmtree(source/'bin')
    with (OUTPUT/'build.log').open('w') as log:
        subprocess.run(args,cwd=source,stdout=log,stderr=subprocess.STDOUT,check=True)
    return source/'DGACC',originals['control.in'].decode(),mesh


def run(binary,control,mesh,mode,minimum,maximum,every):
    folder=OUTPUT/mode.lower();folder.mkdir(exist_ok=True)
    (folder/'msh').mkdir(exist_ok=True)
    (folder/'res').mkdir(exist_ok=True)
    (folder/'msh/A1_Nx11_Ny11.msh').write_bytes(mesh)
    for key,value in [('ACCFLAG',int(mode=='GSIS')),('TMAX',maximum)]:
        control,count=re.subn(rf'(?m)^(\s*{key}\s*=)[^\n]*',rf'\g<1> {value},',control);assert count==1
    (folder/'control.in').write_text(control)
    env=dict(os.environ,OMP_THREAD_LIMIT='4',OPENBLAS_NUM_THREADS='1',
        GRAY_AUDIT_MINIMUM=str(minimum),GRAY_AUDIT_EVERY=str(every),GFORTRAN_UNBUFFERED_ALL='1')
    start=time.monotonic()
    with (folder/'run.log').open('w') as log:
        result=subprocess.run([str(binary)],cwd=folder,env=env,stdout=log,stderr=subprocess.STDOUT)
    elapsed=time.monotonic()-start
    (folder/'elapsed_seconds.txt').write_text(f'{elapsed:.9g}\n')
    if result.returncode:raise RuntimeError((folder/'run.log').read_text()[-4000:])
    print(f'DONE {mode} {elapsed:.1f}s',flush=True)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--gsis-steps',type=int,default=100)
    parser.add_argument('--cis-max',type=int,default=500)
    parser.add_argument('--mode',choices=['both','CIS','GSIS'],default='both')
    args=parser.parse_args();OUTPUT.mkdir(parents=True,exist_ok=True)
    binary,control,mesh=build();print('BUILT latest uploaded Gray reference',flush=True)
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures=[]
        if args.mode in ('both','GSIS'):futures.append(pool.submit(run,binary,control,mesh,'GSIS',args.gsis_steps,args.gsis_steps,10))
        if args.mode in ('both','CIS'):futures.append(pool.submit(run,binary,control,mesh,'CIS',0,args.cis_max,100))
        for future in futures:future.result()


if __name__=='__main__':main()
