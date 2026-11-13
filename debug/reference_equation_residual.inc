! Read-only diagnostic for the unmodified reference kinetic operator.
! Called after each selected iteration, using CURRENT neighbor/wall data.
! The sole global diagnostic side effect is recomputing RHO_WALL, which the
! reference solver recomputes before every subsequent kinetic sweep anyway.
SUBROUTINE Audit_Equation_Residual()
INTEGER :: ri,rm,rl,rj1,rj2,rj3,rside,rf,rb,rneighbor,rpair,rinfo
INTEGER :: rpivot(NDOF_TRI)
INTEGER :: rscaleios
CHARACTER(LEN=128) :: rscaletext
REAL(KIND=DBL) :: rkin(NDOF_TRI,N_TRIS,6),rmassinv(NDOF_TRI,NDOF_TRI,N_TRIS)
REAL(KIND=DBL) :: rwork(NDOF_TRI),rmat(NDOF_TRI,NDOF_TRI)
REAL(KIND=DBL) :: rbase(NDOF_TRI),rorig(NDOF_TRI),rcoupled(NDOF_TRI),rsource(NDOF_TRI)
REAL(KIND=DBL) :: rgain(NDOF_TRI),rsgain(NDOF_TRI),rwall(NDOF_TRI),rbasis(6)
REAL(KIND=DBL) :: rw,rspeed,rwallnormal,rnorm(3),rvolume,rkmass,rflow,rminrho,rmaxrho
REAL(KIND=DBL),SAVE :: rinitial=-1.d0

CALL Calculate_RHO_WALL()
rkin=0.d0
DO rj3=1,NCZ/2
 DO rj2=1,NCY
  DO rj1=1,NCX
   rw=DC(rj1,rj2,rj3)*VDF_EQ(rj1,rj2,rj3)
   rbasis=[1.d0,CX(rj1),CY(rj2),2.d0*C2(rj1,rj2,rj3)/3.d0-1.d0, &
    CX(rj1)*(C2(rj1,rj2,rj3)-2.5d0),CY(rj2)*(C2(rj1,rj2,rj3)-2.5d0)]
   DO rm=1,6
    rkin(:,:,rm)=rkin(:,:,rm)+rw*rbasis(rm)*VDF(:,:,rj1,rj2,rj3)
   END DO
  END DO
 END DO
END DO

rvolume=SUM(TRIANGLES_INF(:,1));rkmass=rvolume;rflow=0.d0
rminrho=MINVAL(1.d0+rkin(:,:,1));rmaxrho=MAXVAL(1.d0+rkin(:,:,1))
DO ri=1,N_TRIS
 rkmass=rkmass+DOT_PRODUCT(INT_NODFUNC_TRI(:,ri),rkin(:,ri,1))
 rflow=rflow+DOT_PRODUCT(INT_NODFUNC_TRI(:,ri),rkin(:,ri,2))
 rmat=INT_NODFUNC_TRI_TRI(:,:,ri)
 CALL DGETRF(NDOF_TRI,NDOF_TRI,rmat,NDOF_TRI,rpivot,rinfo)
 IF(rinfo.NE.0) ERROR STOP 'Residual monitor mass factorization failed'
 CALL DGETRI(NDOF_TRI,rmat,NDOF_TRI,rpivot,rwork,NDOF_TRI,rinfo)
 IF(rinfo.NE.0) ERROR STOP 'Residual monitor mass inverse failed'
 rmassinv(:,:,ri)=rmat
END DO

rnorm=0.d0
DO rj3=1,NCZ/2
 DO rj2=1,NCY
  DO rj1=1,NCX
   rw=DC(rj1,rj2,rj3)*VDF_EQ(rj1,rj2,rj3)
   DO ri=1,N_TRIS
    rbase=DELTA_REF*MATMUL(INT_NODFUNC_TRI_TRI(:,:,ri),VDF(:,ri,rj1,rj2,rj3)) &
     -CX(rj1)*MATMUL(INT_NODFUNC_TRI_TRI_X(ri,:,:),VDF(:,ri,rj1,rj2,rj3)) &
     -CY(rj2)*MATMUL(INT_NODFUNC_TRI_TRI_Y(ri,:,:),VDF(:,ri,rj1,rj2,rj3)) &
     -2.d0*(ACC(1)*CX(rj1)+ACC(2)*CY(rj2))*INT_NODFUNC_TRI(:,ri)
    DO rside=1,3
     rf=TRIANGLES_TAG(ri,3+rside);rb=FACES_TAG(rf,3)
     rspeed=CX(rj1)*TRIANGLES_INF(ri,2*rside)+CY(rj2)*TRIANGLES_INF(ri,2*rside+1)
     rbase=rbase+MAX(rspeed,0.d0)*MATMUL(INT_NODFUNC_TRI_FC(ri,rside,:,:),VDF(:,ri,rj1,rj2,rj3))
     IF(rspeed.GE.0.d0) CYCLE
     IF(rb.EQ.0)THEN
      rneighbor=FACES_TAG(rf,5)
      IF(rneighbor.EQ.ri) rneighbor=FACES_TAG(rf,6)
      rbase=rbase+rspeed*MATMUL(INT_NODFUNC_TRI_TRI_FC(ri,rside,:,:),VDF(:,rneighbor,rj1,rj2,rj3))
     ELSE IF(BC_TYP(rb).EQ.3)THEN
      rpair=FACES_TAG(rf,4);rneighbor=FACES_TAG(rpair,5)
      IF(rneighbor.EQ.0) rneighbor=FACES_TAG(rpair,6)
      rbase=rbase+rspeed*MATMUL(INT_NODFUNC_TRI_TRI_FC(ri,rside,:,:),VDF(:,rneighbor,rj1,rj2,rj3))
     ELSE IF(BC_TYP(rb).EQ.1)THEN
      rwallnormal=-BC_U(rb)*TRIANGLES_INF(ri,2*rside)-BC_V(rb)*TRIANGLES_INF(ri,2*rside+1)
      rwall=(SQRT(PI)*rwallnormal-.5d0*BC_TEMP(rb)+2.d0*(CX(rj1)*BC_U(rb)+CY(rj2)*BC_V(rb)) &
       +BC_TEMP(rb)*(C2(rj1,rj2,rj3)-1.5d0))*INT_NODFUNC_TRI_TRI_FC(ri,rside,:,1) &
       -2.d0*SQRT(PI)*RHO_WALL(:,rf)
      rbase=rbase+rspeed*rwall
     ELSE
      ERROR STOP 'Residual monitor currently supports walls and periodic boundaries only'
     END IF
    END DO
    rgain=rkin(:,ri,1)+2.d0*(CX(rj1)*rkin(:,ri,2)+CY(rj2)*rkin(:,ri,3)) &
     +(C2(rj1,rj2,rj3)-1.5d0)*rkin(:,ri,4) &
     +4.d0*(1.d0-Pr)/5.d0*(CX(rj1)*rkin(:,ri,5)+CY(rj2)*rkin(:,ri,6))*(C2(rj1,rj2,rj3)-2.5d0)
    rsgain=RHO_s(:,ri)+2.d0*(CX(rj1)*VELX_s(:,ri)+CY(rj2)*VELY_s(:,ri)) &
     +(C2(rj1,rj2,rj3)-1.5d0)*TEMP_s(:,ri) &
     +4.d0*(1.d0-Pr)/5.d0*(CX(rj1)*HFX_s(:,ri)+CY(rj2)*HFY_s(:,ri))*(C2(rj1,rj2,rj3)-2.5d0)
    rorig=rbase-DELTA_REF*MATMUL(INT_NODFUNC_TRI_TRI(:,:,ri),rgain)
    rcoupled=rbase-DELTA_REF*MATMUL(INT_NODFUNC_TRI_TRI(:,:,ri),rsgain)
    rsource=DELTA_REF*MATMUL(INT_NODFUNC_TRI_TRI(:,:,ri),rsgain-rgain)
    IF(MAXVAL(ABS(rorig-rcoupled-rsource)).GT.1.d-10*MAX(1.d0,MAXVAL(ABS(rbase)))) &
     ERROR STOP 'Original/coupled/source residual identity failed'
    rnorm(1)=rnorm(1)+rw*DOT_PRODUCT(rorig,MATMUL(rmassinv(:,:,ri),rorig))
    rnorm(2)=rnorm(2)+rw*DOT_PRODUCT(rcoupled,MATMUL(rmassinv(:,:,ri),rcoupled))
    rnorm(3)=rnorm(3)+rw*DOT_PRODUCT(rsource,MATMUL(rmassinv(:,:,ri),rsource))
   END DO
  END DO
 END DO
END DO
rnorm=SQRT(MAX(rnorm,0.d0)/rvolume)
IF(rinitial.LT.0.d0)THEN
 rinitial=MAX(rnorm(1),1.d-100)
 ! Continuation runs can share the zero-state denominator of a baseline run.
 rscaletext=''
 CALL GET_ENVIRONMENT_VARIABLE('REF_RESIDUAL_SCALE',rscaletext)
 IF(LEN_TRIM(rscaletext).GT.0)THEN
  READ(rscaletext,*,IOSTAT=rscaleios)rinitial
  IF(rscaleios.NE.0.OR.rinitial.LE.0.d0) ERROR STOP 'Invalid common residual normalization'
 END IF
 OPEN(94,FILE='equation_residual.csv',STATUS='REPLACE')
 WRITE(94,'(A)')'step,macro_update,relative_original,relative_coupled,relative_source,absolute_original,absolute_coupled,absolute_source,kinetic_mass,mixed_mass,mean_streamwise_velocity,min_density,max_density'
ELSE
 OPEN(94,FILE='equation_residual.csv',STATUS='OLD',POSITION='APPEND')
END IF
WRITE(94,'(I8,12(",",ES25.16))')STEP,RESIDUL,rnorm/rinitial,rnorm,rkmass,MASS,rflow/rvolume,rminrho,rmaxrho
CLOSE(94)
WRITE(*,'(A,I6,A,3ES13.5)')' EQUATION ',STEP,' original/coupled/source ',rnorm/rinitial
END SUBROUTINE Audit_Equation_Residual
