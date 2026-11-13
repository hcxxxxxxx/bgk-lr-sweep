"""Continue the original P4/578-triangle cavity and monitor equation residuals.

Restore both H and the saved mixed collision source, so the GSIS continuation
does not silently replace its source with CIS moments on the first new sweep.
Only the audit driver changes; original solver modules remain byte-identical.
"""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import argparse
import csv
import importlib.util
import json
import os
import re
import shutil
import subprocess
import time

ROOT=Path(__file__).resolve().parents[1]
REFERENCE=ROOT/'reference_runs/nonoscillatory'
ORIGINAL=Path('/Users/hanchenxi/Desktop/NonOscillatory')
OUTPUT=ROOT/'results/reference_cavity_residual_platform'

RESTORE_SOURCE=r'''
SUBROUTINE Audit_Restore_Source()
INTEGER::ai,am,readi,readm,stat
REAL(KIND=DBL)::nodeweight,unusedkin(7),mixed(6),unusedmacro(6)
CALL GET_ENVIRONMENT_VARIABLE('REF_SOURCE_NODES',AUDIT_SOURCE_NODES)
IF(LEN_TRIM(AUDIT_SOURCE_NODES).GT.0)THEN
 OPEN(95,FILE=TRIM(AUDIT_SOURCE_NODES),STATUS='OLD',ACTION='READ')
 DO ai=1,N_TRIS
  DO am=1,NDOF_TRI
   READ(95,*,IOSTAT=stat)readi,readm,nodeweight,unusedkin,mixed,unusedmacro
   IF(stat.NE.0.OR.readi.NE.ai.OR.readm.NE.am) ERROR STOP 'Incompatible saved source nodes'
   IF(ABS(nodeweight-INT_NODFUNC_TRI(am,ai)).GT.1.d-12) ERROR STOP 'Restart mesh mismatch'
   RHO_s(am,ai)=mixed(1);VELX_s(am,ai)=mixed(2);VELY_s(am,ai)=mixed(3)
   TEMP_s(am,ai)=mixed(4);HFX_s(am,ai)=mixed(5);HFY_s(am,ai)=mixed(6)
  END DO
 END DO
 CLOSE(95)
 DO ai=1,N_TRIS
  RHO(ai)=TRIANGLES_INF(ai,1)
  TEMP(ai)=TRIANGLES_INF(ai,1)
  VELX(ai)=0.d0;VELY(ai)=0.d0
  DO am=1,NDOF_TRI
   RHO(ai)=RHO(ai)+RHO_s(am,ai)*INT_NODFUNC_TRI(am,ai)
   TEMP(ai)=TEMP(ai)+TEMP_s(am,ai)*INT_NODFUNC_TRI(am,ai)
   VELX(ai)=VELX(ai)+VELX_s(am,ai)*INT_NODFUNC_TRI(am,ai)
   VELY(ai)=VELY(ai)+VELY_s(am,ai)*INT_NODFUNC_TRI(am,ai)
  END DO
 END DO
 MASS=SUM(RHO)
END IF
IF(LEN_TRIM(AUDIT_RESTART).GT.0)THEN
 RHO_OLD=RHO;TEMP_OLD=TEMP;VELX_OLD=VELX;VELY_OLD=VELY
END IF
END SUBROUTINE Audit_Restore_Source
'''


def build():
    spec=importlib.util.spec_from_file_location('residual_driver',ROOT/'debug/run_reference_residual_platform.py')
    driver=importlib.util.module_from_spec(spec);spec.loader.exec_module(driver)
    driver.OUTPUT=OUTPUT
    binary=driver.build()
    path=OUTPUT/'source/Shakhov_2D_3V_DG.f90'
    text=path.read_text()
    text=text.replace('CHARACTER(LEN=512) :: AUDIT_MODE,AUDIT_RESTART',
                      'CHARACTER(LEN=512) :: AUDIT_MODE,AUDIT_RESTART,AUDIT_SOURCE_NODES,AUDIT_ONLY')
    text=text.replace('CALL Calculate_Macro_Properties ()\nMASS0=MASS',
                      'CALL Calculate_Macro_Properties ()\nCALL Audit_Restore_Source()\nMASS0=MASS',1)
    text=text.replace('CALL Audit_Equation_Residual()\nAUDIT_PREVIOUS=',
                      "CALL Audit_Equation_Residual()\nAUDIT_ONLY=''\nCALL GET_ENVIRONMENT_VARIABLE('REF_AUDIT_ONLY',AUDIT_ONLY)\nIF(TRIM(AUDIT_ONLY).EQ.'1') STOP\nAUDIT_PREVIOUS=",1)
    text=text.replace('END PROGRAM Shakhov_2D_3V_DG',RESTORE_SOURCE+'\nEND PROGRAM Shakhov_2D_3V_DG')
    path.write_text(text)
    # macOS make can compare mtimes at whole-second precision. The first
    # build and this driver edit can land in the same second; force this
    # object to be rebuilt so source restoration is actually in the binary.
    (path.parent/'bin_macos/Shakhov_2D_3V_DG.o').unlink(missing_ok=True)
    with (OUTPUT/'build_cavity.log').open('w') as log:
        subprocess.run(['make','-f','Makefile.macos'],cwd=path.parent,stdout=log,stderr=subprocess.STDOUT,check=True)
    return binary


def prepare(name,steps,tolerance):
    folder=OUTPUT/name;folder.mkdir(parents=True,exist_ok=True)
    control=(ORIGINAL/'control.in').read_text()
    for key,value in [('TMAX',steps),('TOL',f'{tolerance:.16e}')]:
        control,count=re.subn(rf'(?m)^(\s*{key}\s*=)[^\n]*',rf'\g<1> {value},',control)
        assert count==1
    (folder/'control.in').write_text(control)
    shutil.copy2(ORIGINAL/'Rect_A1_Nx17_Ny17.msh',folder/'Rect_A1_Nx17_Ny17.msh')
    return folder


def run(binary,name,mode,steps,tolerance,scale=None,restart=None,source_nodes=None,only=False):
    folder=prepare(name,steps,tolerance)
    env=dict(os.environ,REF_MODE=mode,OPENBLAS_NUM_THREADS='1')
    if scale is not None:env['REF_RESIDUAL_SCALE']=f'{scale:.17e}'
    if restart is not None:env['REF_RESTART']=str(restart)
    if source_nodes is not None:env['REF_SOURCE_NODES']=str(source_nodes)
    if only:env['REF_AUDIT_ONLY']='1'
    start=time.monotonic()
    with (folder/'run.log').open('w') as log:
        result=subprocess.run([str(binary)],cwd=folder,env=env,stdout=log,stderr=subprocess.STDOUT)
    (folder/'elapsed_seconds.txt').write_text(f'{time.monotonic()-start:.9g}\n')
    if result.returncode:
        raise RuntimeError((folder/'run.log').read_text()[-4000:])
    print(f'DONE {name} seconds={time.monotonic()-start:.1f}',flush=True)
    with (folder/'equation_residual.csv').open() as stream:
        return [{key:float(value) for key,value in row.items()} for row in csv.DictReader(stream)]


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--additional-gsis-steps',type=int,default=100)
    args=parser.parse_args();OUTPUT.mkdir(parents=True,exist_ok=True)
    binary=build();print('BUILT original-cavity diagnostic executable',flush=True)
    initial=run(binary,'zero_state','CIS',1,1e-11,only=True)[0]
    scale=initial['absolute_original']
    saved=REFERENCE/'original_config'
    settings=dict(model='2D3V Shakhov',delta_approximately=.1,Kn=8.8623,degree=4,
                  triangles=578,velocity=[24,24,24],stored_velocity=[24,24,12],Pr=2/3,
                  stabilization=[1,1,1,1.875],lid_velocity=1,acceleration=[0,0],
                  wall_type='fully diffuse isothermal',tolerance=1e-11,
                  common_zero_state_equation_residual=scale,
                  cis_restart_total_steps=16,gsis_restart_total_steps=50,
                  additional_gsis_steps=args.additional_gsis_steps,
                  restart_source='H and all six mixed collision-source moments restored for GSIS')
    (OUTPUT/'settings.json').write_text(json.dumps(settings,indent=2)+'\n')
    def cis():
        return run(binary,'cis','CIS',100,1e-11,scale=scale,restart=saved/'cis/final_vdf.bin')
    def gsis():
        return run(binary,'gsis','GSIS',args.additional_gsis_steps,1e-11,scale=scale,
                   restart=saved/'gsis/final_vdf.bin',source_nodes=saved/'gsis/final_nodes.dat')
    with ThreadPoolExecutor(max_workers=2) as pool:
        cf=pool.submit(cis);gf=pool.submit(gsis)
        cis_rows=cf.result();gsis_rows=gf.result()
    rows=[]
    for mode,offset,values in [('CIS',16,cis_rows),('GSIS',50,gsis_rows)]:
        for value in values:
            rows.append(dict(mode=mode,total_step=offset+int(value['step']),**value))
    with (OUTPUT/'summary.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    print('SAVED',OUTPUT/'summary.csv',flush=True)


if __name__=='__main__':main()
