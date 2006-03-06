##
## $Id$
##
## Gromacs 4.0                         Copyright (c) 1991-2003 
## David van der Spoel, Erik Lindahl
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.
##
## To help us fund GROMACS development, we humbly ask that you cite
## the research papers on the package. Check out http://www.gromacs.org
## 
## And Hey:
## Gnomes, ROck Monsters And Chili Sauce
##



.globl nb_kernel310_ia32_3dnow
.globl _nb_kernel310_ia32_3dnow
nb_kernel310_ia32_3dnow:        
_nb_kernel310_ia32_3dnow:       
.set nb310_p_nri, 8
.set nb310_iinr, 12
.set nb310_jindex, 16
.set nb310_jjnr, 20
.set nb310_shift, 24
.set nb310_shiftvec, 28
.set nb310_fshift, 32
.set nb310_gid, 36
.set nb310_pos, 40
.set nb310_faction, 44
.set nb310_charge, 48
.set nb310_p_facel, 52
.set nb310_p_krf, 56
.set nb310_p_crf, 60
.set nb310_Vc, 64
.set nb310_type, 68
.set nb310_p_ntype, 72
.set nb310_vdwparam, 76
.set nb310_Vvdw, 80
.set nb310_p_tabscale, 84
.set nb310_VFtab, 88
.set nb310_invsqrta, 92
.set nb310_dvda, 96
.set nb310_p_gbtabscale, 100
.set nb310_GBtab, 104
.set nb310_p_nthreads, 108
.set nb310_count, 112
.set nb310_mtx, 116
.set nb310_outeriter, 120
.set nb310_inneriter, 124
.set nb310_work, 128
        ## stack offsets for local variables 
.set nb310_is3, 0
.set nb310_ii3, 4
.set nb310_ix, 8
.set nb310_iy, 12
.set nb310_iz, 16
.set nb310_iq, 20
.set nb310_vctot, 28
.set nb310_Vvdwtot, 36
.set nb310_c6, 44
.set nb310_c12, 52
.set nb310_six, 60
.set nb310_twelve, 68
.set nb310_two, 76
.set nb310_n1, 84
.set nb310_tsc, 92
.set nb310_ntia, 100
.set nb310_innerjjnr, 104
.set nb310_innerk, 108
.set nb310_fix, 112
.set nb310_fiy, 116
.set nb310_fiz, 120
.set nb310_dx1, 124
.set nb310_dy1, 128
.set nb310_dz1, 132
.set nb310_dx2, 136
.set nb310_dy2, 140
.set nb310_dz2, 144
.set nb310_n, 148                           ## idx for outer loop
.set nb310_nn1, 152                         ## number of outer iterations
.set nb310_nri, 156
.set nb310_facel, 160
.set nb310_ntype, 164
.set nb310_nouter, 168
.set nb310_ninner, 172
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $176,%esp          ## local stack space 
        femms
        ## move data to local stack
        movl nb310_p_nri(%ebp),%ecx
        movl nb310_p_ntype(%ebp),%edx
        movl nb310_p_facel(%ebp),%esi
        movl nb310_p_tabscale(%ebp),%edi
        movl (%ecx),%ecx
        movl (%edx),%edx
        movl (%esi),%esi
        movl %ecx,nb310_nri(%esp)
        movl %edx,nb310_ntype(%esp)
        movl %esi,nb310_facel(%esp)

        movd  (%edi),%mm3
        punpckldq %mm3,%mm3
        movq  %mm3,nb310_tsc(%esp)
        movl $0x40000000,%eax  ## 2.0
        movl %eax,nb310_two(%esp)
        movl %eax,nb310_two+4(%esp)
        movl $0x40c00000,%ebx ## 6.0
        movl %ebx,nb310_six(%esp)
        movl %ebx,nb310_six+4(%esp)
        movl $0x41400000,%ecx  ## 12.0
        movl %ecx,nb310_twelve(%esp)
        movl %ecx,nb310_twelve+4(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb310_nouter(%esp)
        movl %eax,nb310_ninner(%esp)

_nb_kernel310_ia32_3dnow.nb310_threadloop: 
        movl  nb310_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel310_ia32_3dnow.nb310_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel310_ia32_3dnow.nb310_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb310_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb310_n(%esp)
        movl %ebx,nb310_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel310_ia32_3dnow.nb310_outerstart
        jmp _nb_kernel310_ia32_3dnow.nb310_end

_nb_kernel310_ia32_3dnow.nb310_outerstart: 
        ## ebx contains number of outer iterations
        addl nb310_nouter(%esp),%ebx
        movl %ebx,nb310_nouter(%esp)

_nb_kernel310_ia32_3dnow.nb310_outer: 
        movl  nb310_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb310_is3(%esp)      ## store is3 

        movl  nb310_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movq  (%eax,%ebx,4),%mm0        ## move shX/shY to mm0 and shZ to mm1 
        movd  8(%eax,%ebx,4),%mm1

        movl  nb310_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx=ii 

        movl  nb310_charge(%ebp),%edx
        movd  (%edx,%ebx,4),%mm2    ## mm2=charge[ii] 
        pfmul nb310_facel(%esp),%mm2
        punpckldq %mm2,%mm2         ## spread to both halves 
        movq  %mm2,nb310_iq(%esp)           ## iq =facel*charge[ii] 

        movl  nb310_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb310_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb310_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb310_pos(%ebp),%eax      ## eax = base of pos[] 

        pfadd (%eax,%ebx,4),%mm0    ## ix = shX + posX (and iy too) 
        movd  8(%eax,%ebx,4),%mm3       ## cant use direct memory add for 4 bytes (iz) 
        movl  %ebx,nb310_ii3(%esp)
        pfadd %mm3,%mm1
        movq  %mm0,nb310_ix(%esp)
        movd  %mm1,nb310_iz(%esp)

        ## clear total potential and i forces 
        pxor  %mm7,%mm7
        movq  %mm7,nb310_vctot(%esp)
        movq  %mm7,nb310_Vvdwtot(%esp)
        movq  %mm7,nb310_fix(%esp)
        movd  %mm7,nb310_fiz(%esp)

        movl  nb310_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb310_pos(%ebp),%esi
        movl  nb310_faction(%ebp),%edi
        movl  nb310_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb310_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb310_ninner(%esp),%ecx
        movl  %ecx,nb310_ninner(%esp)
        movl  %edx,nb310_innerk(%esp)      ## number of innerloop atoms 
        addl  $0,%edx
        jge   _nb_kernel310_ia32_3dnow.nb310_unroll_loop
        jmp   _nb_kernel310_ia32_3dnow.nb310_finish_inner
_nb_kernel310_ia32_3dnow.nb310_unroll_loop: 
        ## paired innerloop starts here 
        movl  nb310_innerjjnr(%esp),%ecx       ## pointer to jjnr[k] 
        movl  (%ecx),%eax
        movl  4(%ecx),%ebx           ## eax/ebx=jnr 
        addl $8,nb310_innerjjnr(%esp)             ## advance pointer (unrolled 2) 
        prefetch 16(%ecx)            ## prefetch data - trial and error says 16 is best 

        movl nb310_charge(%ebp),%ecx     ## base of charge[] 
        movq nb310_iq(%esp),%mm5
        movd (%ecx,%eax,4),%mm3      ## charge[jnr1] 
        punpckldq (%ecx,%ebx,4),%mm3     ## move charge 2 to high part of mm3 
        pfmul %mm5,%mm3              ## mm3 now has qq for both particles 

        movl nb310_type(%ebp),%ecx
        movl (%ecx,%eax,4),%edx          ## type [jnr1] 
        movl (%ecx,%ebx,4),%ecx      ## type [jnr2] 

        movl nb310_vdwparam(%ebp),%esi          ## base of vdwparam  
        shll %edx
        shll %ecx
        addl nb310_ntia(%esp),%edx           ## tja = ntia + 2*type 
        addl nb310_ntia(%esp),%ecx

        movq (%esi,%edx,4),%mm5         ## mm5 = 1st c6 / c12           
        movq (%esi,%ecx,4),%mm7         ## mm7 = 2nd c6 / c12   
        movq %mm5,%mm6
        punpckldq %mm7,%mm5             ## mm5 = 1st c6 / 2nd c6 
        punpckhdq %mm7,%mm6             ## mm6 = 1st c12 / 2nd c12 
        movq %mm5,nb310_c6(%esp)
        movq %mm6,nb310_c12(%esp)

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        movl  nb310_pos(%ebp),%esi

        movq  nb310_ix(%esp),%mm0
        movd  nb310_iz(%esp),%mm1
        movq  (%esi,%eax,4),%mm4     ## fetch first j coordinates 
        movd  8(%esi,%eax,4),%mm5
        pfsubr %mm0,%mm4             ## dr = ir - jr  
        pfsubr %mm1,%mm5
        movq  %mm4,nb310_dx1(%esp)           ## store dr 
        movd  %mm5,nb310_dz1(%esp)
        pfmul %mm4,%mm4              ## square dx,dy,dz                          
        pfmul %mm5,%mm5
        pfacc %mm5,%mm4              ## accumulate to get dx*dx+ dy*dy+ dz*dz 
        pfacc %mm5,%mm4              ## first rsq in lower mm4 

        movq  (%esi,%ebx,4),%mm6     ## fetch second j coordinates  
        movd  8(%esi,%ebx,4),%mm7

        pfsubr %mm0,%mm6             ## dr = ir - jr  
        pfsubr %mm1,%mm7
        movq  %mm6,nb310_dx2(%esp)           ## store dr 
        movd  %mm7,nb310_dz2(%esp)
        pfmul %mm6,%mm6              ## square dx,dy,dz 
        pfmul %mm7,%mm7
        pfacc %mm7,%mm6              ## accumulate to get dx*dx+ dy*dy+ dz*dz 
        pfacc %mm7,%mm6              ## second rsq in lower mm6 

        pfrsqrt %mm4,%mm0            ## lookup inverse square root seed 
        pfrsqrt %mm6,%mm1


        punpckldq %mm1,%mm0
        punpckldq %mm6,%mm4             ## now 4 has rsq and 0 the seed for both pairs. 
        movq %mm0,%mm2                  ## amd 3dnow N-R iteration to get full precision. 
        pfmul %mm0,%mm0
        pfrsqit1 %mm4,%mm0
        pfrcpit2 %mm2,%mm0
        pfmul %mm0,%mm4
        movq %mm4,%mm1
        ## mm0 is invsqrt, and mm1 r. 
        ## do potential and fscal 
        pfmul nb310_tsc(%esp),%mm1      ## mm1=rt 
        pf2iw %mm1,%mm4
        movq %mm4,nb310_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm1              ## now mm1 is eps and mm4 is n0 

        movq %mm1,%mm2
        pfmul %mm2,%mm2 ## mm1 is eps, mm2 is eps2 

        movl nb310_VFtab(%ebp),%edx
        movl nb310_n1(%esp),%ecx
        shll $2,%ecx
        ## coulomb table 
        ## load all the table values we need 
        movd (%edx,%ecx,4),%mm4
        movd 4(%edx,%ecx,4),%mm5
        movd 8(%edx,%ecx,4),%mm6
        movd 12(%edx,%ecx,4),%mm7
        movl nb310_n1+4(%esp),%ecx
        shll $2,%ecx
        punpckldq (%edx,%ecx,4),%mm4
        punpckldq 4(%edx,%ecx,4),%mm5
        punpckldq 8(%edx,%ecx,4),%mm6
        punpckldq 12(%edx,%ecx,4),%mm7

        pfmul %mm1,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 

        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 

        pfmul nb310_two(%esp),%mm7      ## two*Heps2 
        pfadd %mm6,%mm7
        pfadd %mm5,%mm7 ## mm7=FF 

        pfmul %mm1,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        pfmul %mm3,%mm5 ## vcoul=qq*VV 
        pfmul %mm7,%mm3 ## fijC=FF*qq  

        movq %mm0,%mm1
        pfmul %mm1,%mm1 ## mm1=invsq 
        movq %mm1,%mm2
        pfmul %mm1,%mm2
        pfmul %mm1,%mm2 ## mm2=rinvsix 
        movq  %mm2,%mm1
        pfmul %mm1,%mm1 ## mm1=rinvtwelve 

        pfmul nb310_tsc(%esp),%mm3

        pfmul nb310_c12(%esp),%mm1

        pfmul nb310_c6(%esp),%mm2

        movq %mm1,%mm4
        pfsub %mm2,%mm4 ## mm4 = Vvdw12-Vvdw6 

        pfmul nb310_six(%esp),%mm2
        pfmul nb310_twelve(%esp),%mm1

        pfsub %mm2,%mm1
        pfmul %mm0,%mm1 ## mm1= (12*Vvdw12-6*Vvdw6)*rinv11 

        pfsub %mm3,%mm1

        ## update vctot 
        pfadd nb310_vctot(%esp),%mm5        ## add the earlier value 
        movq %mm5,nb310_vctot(%esp)         ## store the sum       

        pfmul %mm1,%mm0   ## mm0 is total fscal now     

        prefetchw nb310_dx1(%esp)       ## prefetch i forces to cache 

        ## spread fscalar to both positions 
        movq %mm0,%mm1
        punpckldq %mm0,%mm0
        punpckhdq %mm1,%mm1

        ## calc vector force 
        prefetchw (%edi,%eax,4) ## prefetch the 1st faction to cache 
        movq nb310_dx1(%esp),%mm2       ## fetch dr 
        movd nb310_dz1(%esp),%mm3

        ## update Vvdwtot 
        pfadd nb310_Vvdwtot(%esp),%mm4        ## add the earlier value 
        movq %mm4,nb310_Vvdwtot(%esp)         ## store the sum       

        prefetchw (%edi,%ebx,4) ## prefetch the 2nd faction to cache 
        pfmul %mm0,%mm2         ## mult by fs  
        pfmul %mm0,%mm3

        movq nb310_dx2(%esp),%mm4       ## fetch dr 
        movd nb310_dz2(%esp),%mm5
        pfmul %mm1,%mm4         ## mult by fs  
        pfmul %mm1,%mm5
        ## update i forces 

        movq nb310_fix(%esp),%mm0
        movd nb310_fiz(%esp),%mm1
        pfadd %mm2,%mm0
        pfadd %mm3,%mm1

        pfadd %mm4,%mm0
        pfadd %mm5,%mm1
        movq %mm0,nb310_fix(%esp)
        movd %mm1,nb310_fiz(%esp)
        ## update j forces 

        movq (%edi,%eax,4),%mm0
        movd 8(%edi,%eax,4),%mm1
        movq (%edi,%ebx,4),%mm6
        movd 8(%edi,%ebx,4),%mm7

        pfsub %mm2,%mm0
        pfsub %mm3,%mm1
        pfsub %mm4,%mm6
        pfsub %mm5,%mm7

        movq %mm0,(%edi,%eax,4)
        movd %mm1,8(%edi,%eax,4)
        movq %mm6,(%edi,%ebx,4)
        movd %mm7,8(%edi,%ebx,4)

        ## should we do one more iteration? 
        subl $2,nb310_innerk(%esp)
        jl    _nb_kernel310_ia32_3dnow.nb310_finish_inner
        jmp   _nb_kernel310_ia32_3dnow.nb310_unroll_loop
_nb_kernel310_ia32_3dnow.nb310_finish_inner: 
        andl $1,nb310_innerk(%esp)
        jnz  _nb_kernel310_ia32_3dnow.nb310_single_inner
        jmp  _nb_kernel310_ia32_3dnow.nb310_updateouterdata
_nb_kernel310_ia32_3dnow.nb310_single_inner: 
        ## a single j particle iteration here - compare with the unrolled code for comments. 
        movl  nb310_innerjjnr(%esp),%eax
        movl  (%eax),%eax       ## eax=jnr offset 

        movl nb310_charge(%ebp),%ecx
        movd nb310_iq(%esp),%mm5
        movd (%ecx,%eax,4),%mm3
        pfmul %mm5,%mm3         ## mm3=qq 

        movl nb310_vdwparam(%ebp),%esi
        movl nb310_type(%ebp),%ecx
        movl (%ecx,%eax,4),%edx          ## type [jnr1] 
        shll %edx
        addl nb310_ntia(%esp),%edx           ## tja = ntia + 2*type 
        movd (%esi,%edx,4),%mm5         ## mm5 = 1st c6                 
        movq %mm5,nb310_c6(%esp)
        movd 4(%esi,%edx,4),%mm5        ## mm5 = 1st c12                
        movq %mm5,nb310_c12(%esp)


        movl  nb310_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        movq  nb310_ix(%esp),%mm0
        movd  nb310_iz(%esp),%mm1
        movq  (%esi,%eax,4),%mm4
        movd  8(%esi,%eax,4),%mm5
        pfsubr %mm0,%mm4
        pfsubr %mm1,%mm5
        movq  %mm4,nb310_dx1(%esp)
        pfmul %mm4,%mm4
        movd  %mm5,nb310_dz1(%esp)
        pfmul %mm5,%mm5
        pfacc %mm5,%mm4
        pfacc %mm5,%mm4         ## mm4=rsq 

        pfrsqrt %mm4,%mm0
        movq %mm0,%mm2
        pfmul %mm0,%mm0
        pfrsqit1 %mm4,%mm0
        pfrcpit2 %mm2,%mm0      ## mm1=invsqrt 
        pfmul %mm0,%mm4
        movq %mm4,%mm1
        ## mm0 is invsqrt, and mm1 r. 
        ## calculate potentials and scalar force 
        pfmul nb310_tsc(%esp),%mm1      ## mm1=rt 
        pf2iw %mm1,%mm4
        movd %mm4,nb310_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm1              ## now mm1 is eps and mm4 is n0 

        movq %mm1,%mm2
        pfmul %mm2,%mm2 ## mm1 is eps, mm2 is eps2 

        ## coulomb table 
        movl nb310_VFtab(%ebp),%edx
        movl nb310_n1(%esp),%ecx
        shll $2,%ecx
        ## load all the table values we need 
        movd (%edx,%ecx,4),%mm4
        movd 4(%edx,%ecx,4),%mm5
        movd 8(%edx,%ecx,4),%mm6
        movd 12(%edx,%ecx,4),%mm7

        pfmul %mm1,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 

        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 

        pfmul nb310_two(%esp),%mm7      ## two*Heps2 
        pfadd %mm6,%mm7
        pfadd %mm5,%mm7 ## mm7=FF 

        pfmul %mm1,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        pfmul %mm3,%mm5 ## vcoul=qq*VV 
        pfmul %mm7,%mm3 ## fijC=FF*qq  

        ## at this point mm5 contains vcoul and mm3 fijC 

        movq %mm0,%mm1
        pfmul %mm1,%mm1 ## mm1=invsq 
        movq %mm1,%mm2
        pfmul %mm1,%mm2
        pfmul %mm1,%mm2 ## mm2=rinvsix 
        movq  %mm2,%mm1
        pfmul %mm1,%mm1 ## mm1=rinvtwelve 

        pfmul nb310_tsc(%esp),%mm3

        pfmul nb310_c12(%esp),%mm1

        pfmul nb310_c6(%esp),%mm2

        movq %mm1,%mm4
        pfsub %mm2,%mm4 ## mm4 = Vvdw12-Vvdw6 

        pfmul nb310_six(%esp),%mm2
        pfmul nb310_twelve(%esp),%mm1

        pfsub %mm2,%mm1
        pfmul %mm0,%mm1 ## mm1= (12*Vvdw12-6*Vvdw6)*rinv11 

        pfsub %mm3,%mm1

        ## update vctot 
        pfadd nb310_vctot(%esp),%mm5        ## add the earlier value 
        movq %mm5,nb310_vctot(%esp)         ## store the sum       

        pfmul %mm1,%mm0   ## mm0 is total fscal now     

        ## spread fscalar to both positions 
        punpckldq %mm0,%mm0
        ## calc vectorial force 
        prefetchw (%edi,%eax,4) ## prefetch faction to cache  
        movq nb310_dx1(%esp),%mm2
        movd nb310_dz1(%esp),%mm3

        ## update Vvdwtot 
        pfadd nb310_Vvdwtot(%esp),%mm4        ## add the earlier value 
        movq %mm4,nb310_Vvdwtot(%esp)         ## store the sum       

        pfmul %mm0,%mm2
        pfmul %mm0,%mm3

        ## update i particle force 
        movq nb310_fix(%esp),%mm0
        movd nb310_fiz(%esp),%mm1
        pfadd %mm2,%mm0
        pfadd %mm3,%mm1
        movq %mm0,nb310_fix(%esp)
        movd %mm1,nb310_fiz(%esp)
        ## update j particle force 
        movq (%edi,%eax,4),%mm0
        movd 8(%edi,%eax,4),%mm1
        pfsub %mm2,%mm0
        pfsub %mm3,%mm1
        movq %mm0,(%edi,%eax,4)
        movd %mm1,8(%edi,%eax,4)
        ## done! 
_nb_kernel310_ia32_3dnow.nb310_updateouterdata: 
        movl  nb310_ii3(%esp),%ecx

        movq  (%edi,%ecx,4),%mm6       ## increment i force 
        movd  8(%edi,%ecx,4),%mm7
        pfadd nb310_fix(%esp),%mm6
        pfadd nb310_fiz(%esp),%mm7
        movq  %mm6,(%edi,%ecx,4)
        movd  %mm7,8(%edi,%ecx,4)

        movl  nb310_fshift(%ebp),%ebx      ## increment fshift force 
        movl  nb310_is3(%esp),%edx

        movq  (%ebx,%edx,4),%mm6
        movd  8(%ebx,%edx,4),%mm7
        pfadd nb310_fix(%esp),%mm6
        pfadd nb310_fiz(%esp),%mm7
        movq  %mm6,(%ebx,%edx,4)
        movd  %mm7,8(%ebx,%edx,4)

        ## get n from stack
        movl nb310_n(%esp),%esi
        ## get group index for i particle 
        movl  nb310_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        movq  nb310_vctot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb310_Vc(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment vc[gid] 

        movq  nb310_Vvdwtot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb310_Vvdw(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment Vvdw[gid] 

        ## finish if last 
        movl nb310_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel310_ia32_3dnow.nb310_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb310_n(%esp)
        jmp _nb_kernel310_ia32_3dnow.nb310_outer
_nb_kernel310_ia32_3dnow.nb310_outerend: 
        ## check if more outer neighborlists remain
        movl  nb310_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel310_ia32_3dnow.nb310_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel310_ia32_3dnow.nb310_threadloop
_nb_kernel310_ia32_3dnow.nb310_end: 
        femms

        movl nb310_nouter(%esp),%eax
        movl nb310_ninner(%esp),%ebx
        movl nb310_outeriter(%ebp),%ecx
        movl nb310_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        addl $176,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret







.globl nb_kernel310nf_ia32_3dnow
.globl _nb_kernel310nf_ia32_3dnow
nb_kernel310nf_ia32_3dnow:      
_nb_kernel310nf_ia32_3dnow:     
.set nb310nf_p_nri, 8
.set nb310nf_iinr, 12
.set nb310nf_jindex, 16
.set nb310nf_jjnr, 20
.set nb310nf_shift, 24
.set nb310nf_shiftvec, 28
.set nb310nf_fshift, 32
.set nb310nf_gid, 36
.set nb310nf_pos, 40
.set nb310nf_faction, 44
.set nb310nf_charge, 48
.set nb310nf_p_facel, 52
.set nb310nf_p_krf, 56
.set nb310nf_p_crf, 60
.set nb310nf_Vc, 64
.set nb310nf_type, 68
.set nb310nf_p_ntype, 72
.set nb310nf_vdwparam, 76
.set nb310nf_Vvdw, 80
.set nb310nf_p_tabscale, 84
.set nb310nf_VFtab, 88
.set nb310nf_invsqrta, 92
.set nb310nf_dvda, 96
.set nb310nf_p_gbtabscale, 100
.set nb310nf_GBtab, 104
.set nb310nf_p_nthreads, 108
.set nb310nf_count, 112
.set nb310nf_mtx, 116
.set nb310nf_outeriter, 120
.set nb310nf_inneriter, 124
.set nb310nf_work, 128
        ## stack offsets for local variables 
.set nb310nf_is3, 0
.set nb310nf_ii3, 4
.set nb310nf_ix, 8
.set nb310nf_iy, 12
.set nb310nf_iz, 16
.set nb310nf_iq, 20
.set nb310nf_vctot, 28
.set nb310nf_Vvdwtot, 36
.set nb310nf_c6, 44
.set nb310nf_c12, 52
.set nb310nf_n1, 60
.set nb310nf_tsc, 68
.set nb310nf_ntia, 76
.set nb310nf_innerjjnr, 80
.set nb310nf_innerk, 84
.set nb310nf_n, 88                         ## idx for outer loop
.set nb310nf_nn1, 92                       ## number of outer iterations
.set nb310nf_nri, 96
.set nb310nf_facel, 100
.set nb310nf_ntype, 104
.set nb310nf_nouter, 108
.set nb310nf_ninner, 112
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $116,%esp          ## local stack space 
        femms
        ## move data to local stack  
        movl nb310nf_p_nri(%ebp),%ecx
        movl nb310nf_p_ntype(%ebp),%edx
        movl nb310nf_p_facel(%ebp),%esi
        movl nb310nf_p_tabscale(%ebp),%edi
        movl (%ecx),%ecx
        movl (%edx),%edx
        movl (%esi),%esi
        movl %ecx,nb310nf_nri(%esp)
        movl %edx,nb310nf_ntype(%esp)
        movl %esi,nb310nf_facel(%esp)

        movd  (%edi),%mm3
        punpckldq %mm3,%mm3
        movq  %mm3,nb310nf_tsc(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb310nf_nouter(%esp)
        movl %eax,nb310nf_ninner(%esp)

_nb_kernel310nf_ia32_3dnow.nb310nf_threadloop: 
        movl  nb310nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel310nf_ia32_3dnow.nb310nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel310nf_ia32_3dnow.nb310nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb310nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb310nf_n(%esp)
        movl %ebx,nb310nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel310nf_ia32_3dnow.nb310nf_outerstart
        jmp _nb_kernel310nf_ia32_3dnow.nb310nf_end

_nb_kernel310nf_ia32_3dnow.nb310nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb310nf_nouter(%esp),%ebx
        movl %ebx,nb310nf_nouter(%esp)

_nb_kernel310nf_ia32_3dnow.nb310nf_outer: 
        movl  nb310nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 

        movl  nb310nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movq  (%eax,%ebx,4),%mm0        ## move shX/shY to mm0 and shZ to mm1 
        movd  8(%eax,%ebx,4),%mm1

        movl  nb310nf_iinr(%esp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx=ii 

        movl  nb310nf_charge(%ebp),%edx
        movd  (%edx,%ebx,4),%mm2    ## mm2=charge[ii] 
        pfmul nb310nf_facel(%esp),%mm2
        punpckldq %mm2,%mm2         ## spread to both halves 
        movq  %mm2,nb310nf_iq(%esp)         ## iq =facel*charge[ii] 

        movl  nb310nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb310nf_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb310nf_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb310nf_pos(%ebp),%eax      ## eax = base of pos[] 

        pfadd (%eax,%ebx,4),%mm0    ## ix = shX + posX (and iy too) 
        movd  8(%eax,%ebx,4),%mm3       ## cant use direct memory add for 4 bytes (iz) 
        pfadd %mm3,%mm1
        movq  %mm0,nb310nf_ix(%esp)
        movd  %mm1,nb310nf_iz(%esp)

        ## clear total potential and i forces 
        pxor  %mm7,%mm7
        movq  %mm7,nb310nf_vctot(%esp)
        movq  %mm7,nb310nf_Vvdwtot(%esp)

        movl  nb310nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb310nf_pos(%ebp),%esi
        movl  nb310nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb310nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb310nf_ninner(%esp),%ecx
        movl  %ecx,nb310nf_ninner(%esp)
        movl  %edx,nb310nf_innerk(%esp)      ## number of innerloop atoms 
        addl  $0,%edx
        jmp _nb_kernel310nf_ia32_3dnow.nb310nf_finish_inner
        jge   _nb_kernel310nf_ia32_3dnow.nb310nf_unroll_loop
        jmp   _nb_kernel310nf_ia32_3dnow.nb310nf_finish_inner
_nb_kernel310nf_ia32_3dnow.nb310nf_unroll_loop: 
        ## paired innerloop starts here 
        movl  nb310nf_innerjjnr(%esp),%ecx       ## pointer to jjnr[k] 
        movl  (%ecx),%eax
        movl  4(%ecx),%ebx           ## eax/ebx=jnr 
        addl $8,nb310nf_innerjjnr(%esp)             ## advance pointer (unrolled 2) 
        prefetch 16(%ecx)            ## prefetch data - trial and error says 16 is best 

        movl nb310nf_charge(%ebp),%ecx     ## base of charge[] 
        movq nb310nf_iq(%esp),%mm5
        movd (%ecx,%eax,4),%mm3      ## charge[jnr1] 
        punpckldq (%ecx,%ebx,4),%mm3     ## move charge 2 to high part of mm3 
        pfmul %mm5,%mm3              ## mm3 now has qq for both particles 

        movl nb310nf_type(%ebp),%ecx
        movl (%ecx,%eax,4),%edx          ## type [jnr1] 
        movl (%ecx,%ebx,4),%ecx      ## type [jnr2] 

        movl nb310nf_vdwparam(%ebp),%esi                ## base of vdwparam  
        shll %edx
        shll %ecx
        addl nb310nf_ntia(%esp),%edx         ## tja = ntia + 2*type 
        addl nb310nf_ntia(%esp),%ecx

        movq (%esi,%edx,4),%mm5         ## mm5 = 1st c6 / c12           
        movq (%esi,%ecx,4),%mm7         ## mm7 = 2nd c6 / c12   
        movq %mm5,%mm6
        punpckldq %mm7,%mm5             ## mm5 = 1st c6 / 2nd c6 
        punpckhdq %mm7,%mm6             ## mm6 = 1st c12 / 2nd c12 
        movq %mm5,nb310nf_c6(%esp)
        movq %mm6,nb310nf_c12(%esp)

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        movl  nb310nf_pos(%ebp),%esi

        movq  nb310nf_ix(%esp),%mm0
        movd  nb310nf_iz(%esp),%mm1
        movq  (%esi,%eax,4),%mm4     ## fetch first j coordinates 
        movd  8(%esi,%eax,4),%mm5
        pfsubr %mm0,%mm4             ## dr = ir - jr  
        pfsubr %mm1,%mm5
        pfmul %mm4,%mm4              ## square dx,dy,dz                          
        pfmul %mm5,%mm5
        pfacc %mm5,%mm4              ## accumulate to get dx*dx+ dy*dy+ dz*dz 
        pfacc %mm5,%mm4              ## first rsq in lower mm4 

        movq  (%esi,%ebx,4),%mm6     ## fetch second j coordinates  
        movd  8(%esi,%ebx,4),%mm7

        pfsubr %mm0,%mm6             ## dr = ir - jr  
        pfsubr %mm1,%mm7
        pfmul %mm6,%mm6              ## square dx,dy,dz 
        pfmul %mm7,%mm7
        pfacc %mm7,%mm6              ## accumulate to get dx*dx+ dy*dy+ dz*dz 
        pfacc %mm7,%mm6              ## second rsq in lower mm6 

        pfrsqrt %mm4,%mm0            ## lookup inverse square root seed 
        pfrsqrt %mm6,%mm1

        punpckldq %mm1,%mm0
        punpckldq %mm6,%mm4             ## now 4 has rsq and 0 the seed for both pairs. 
        movq %mm0,%mm2                  ## amd 3dnow N-R iteration to get full precision. 
        pfmul %mm0,%mm0
        pfrsqit1 %mm4,%mm0
        pfrcpit2 %mm2,%mm0
        pfmul %mm0,%mm4
        movq %mm4,%mm1
        ## mm0 is invsqrt, and mm1 r. 
        ## do potential and fscal 
        pfmul nb310nf_tsc(%esp),%mm1    ## mm1=rt 
        pf2iw %mm1,%mm4
        movq %mm4,nb310nf_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm1              ## now mm1 is eps and mm4 is n0 

        movq %mm1,%mm2
        pfmul %mm2,%mm2 ## mm1 is eps, mm2 is eps2 

        movl nb310nf_VFtab(%ebp),%edx
        movl nb310nf_n1(%esp),%ecx
        shll $2,%ecx
        ## coulomb table 
        ## load all the table values we need 
        movd (%edx,%ecx,4),%mm4
        movd 4(%edx,%ecx,4),%mm5
        movd 8(%edx,%ecx,4),%mm6
        movd 12(%edx,%ecx,4),%mm7
        movl nb310nf_n1+4(%esp),%ecx
        shll $2,%ecx
        punpckldq (%edx,%ecx,4),%mm4
        punpckldq 4(%edx,%ecx,4),%mm5
        punpckldq 8(%edx,%ecx,4),%mm6
        punpckldq 12(%edx,%ecx,4),%mm7

        pfmul %mm1,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 

        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 

        pfmul %mm1,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        pfmul %mm3,%mm5 ## vcoul=qq*VV 

        movq %mm0,%mm1
        pfmul %mm1,%mm1 ## mm1=invsq 
        movq %mm1,%mm2
        pfmul %mm1,%mm2
        pfmul %mm1,%mm2 ## mm2=rinvsix 
        movq  %mm2,%mm1
        pfmul %mm1,%mm1 ## mm1=rinvtwelve 

        pfmul nb310nf_tsc(%esp),%mm3

        pfmul nb310nf_c12(%esp),%mm1

        pfmul nb310nf_c6(%esp),%mm2

        movq %mm1,%mm4
        pfsub %mm2,%mm4 ## mm4 = Vvdw12-Vvdw6 
        ## update vctot 
        pfadd nb310nf_vctot(%esp),%mm5        ## add the earlier value 
        movq %mm5,nb310nf_vctot(%esp)         ## store the sum       
        ## update Vvdwtot 
        pfadd nb310nf_Vvdwtot(%esp),%mm4        ## add the earlier value 
        movq %mm4,nb310nf_Vvdwtot(%esp)         ## store the sum       

        ## should we do one more iteration? 
        subl $2,nb310nf_innerk(%esp)
        jl    _nb_kernel310nf_ia32_3dnow.nb310nf_finish_inner
        jmp   _nb_kernel310nf_ia32_3dnow.nb310nf_unroll_loop
_nb_kernel310nf_ia32_3dnow.nb310nf_finish_inner: 
        andl $1,nb310nf_innerk(%esp)
        jnz  _nb_kernel310nf_ia32_3dnow.nb310nf_single_inner
        jmp  _nb_kernel310nf_ia32_3dnow.nb310nf_updateouterdata
_nb_kernel310nf_ia32_3dnow.nb310nf_single_inner: 
        ## a single j particle iteration here - compare with the unrolled code for comments. 
        movl  nb310nf_innerjjnr(%esp),%eax
        movl  (%eax),%eax       ## eax=jnr offset 

        movl nb310nf_charge(%ebp),%ecx
        movd nb310nf_iq(%esp),%mm5
        movd (%ecx,%eax,4),%mm3
        pfmul %mm5,%mm3         ## mm3=qq 

        movl nb310nf_vdwparam(%ebp),%esi
        movl nb310nf_type(%ebp),%ecx
        movl (%ecx,%eax,4),%edx          ## type [jnr1] 
        shll %edx
        addl nb310nf_ntia(%esp),%edx         ## tja = ntia + 2*type 
        movd (%esi,%edx,4),%mm5         ## mm5 = 1st c6                 
        movq %mm5,nb310nf_c6(%esp)
        movd 4(%esi,%edx,4),%mm5        ## mm5 = 1st c12                
        movq %mm5,nb310nf_c12(%esp)


        movl  nb310nf_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        movq  nb310nf_ix(%esp),%mm0
        movd  nb310nf_iz(%esp),%mm1
        movq  (%esi,%eax,4),%mm4
        movd  8(%esi,%eax,4),%mm5
        pfsubr %mm0,%mm4
        pfsubr %mm1,%mm5
        pfmul %mm4,%mm4
        pfmul %mm5,%mm5
        pfacc %mm5,%mm4
        pfacc %mm5,%mm4         ## mm4=rsq 

        pfrsqrt %mm4,%mm0
        movq %mm0,%mm2
        pfmul %mm0,%mm0
        pfrsqit1 %mm4,%mm0
        pfrcpit2 %mm2,%mm0      ## mm1=invsqrt 
        pfmul %mm0,%mm4
        movq %mm4,%mm1
        ## mm0 is invsqrt, and mm1 r. 
        ## calculate potentials and scalar force 
        pfmul nb310nf_tsc(%esp),%mm1    ## mm1=rt 
        pf2iw %mm1,%mm4
        movd %mm4,nb310nf_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm1              ## now mm1 is eps and mm4 is n0 

        movq %mm1,%mm2
        pfmul %mm2,%mm2 ## mm1 is eps, mm2 is eps2 

        ## coulomb table 
        movl nb310nf_VFtab(%ebp),%edx
        movl nb310nf_n1(%esp),%ecx
        shll $2,%ecx
        ## load all the table values we need 
        movd (%edx,%ecx,4),%mm4
        movd 4(%edx,%ecx,4),%mm5
        movd 8(%edx,%ecx,4),%mm6
        movd 12(%edx,%ecx,4),%mm7

        pfmul %mm1,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 

        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 

        pfmul %mm1,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        pfmul %mm3,%mm5 ## vcoul=qq*VV 
        ## at this point mm5 contains vcoul 

        movq %mm0,%mm1
        pfmul %mm1,%mm1 ## mm1=invsq 
        movq %mm1,%mm2
        pfmul %mm1,%mm2
        pfmul %mm1,%mm2 ## mm2=rinvsix 
        movq  %mm2,%mm1
        pfmul %mm1,%mm1 ## mm1=rinvtwelve 

        pfmul nb310nf_tsc(%esp),%mm3

        pfmul nb310nf_c12(%esp),%mm1

        pfmul nb310nf_c6(%esp),%mm2

        movq %mm1,%mm4
        pfsub %mm2,%mm4 ## mm4 = Vvdw12-Vvdw6 
        ## update vctot 
        pfadd nb310nf_vctot(%esp),%mm5        ## add the earlier value 
        movq %mm5,nb310nf_vctot(%esp)         ## store the sum       
        ## update Vvdwtot 
        pfadd nb310nf_Vvdwtot(%esp),%mm4        ## add the earlier value 
        movq %mm4,nb310nf_Vvdwtot(%esp)         ## store the sum       

_nb_kernel310nf_ia32_3dnow.nb310nf_updateouterdata: 
        ## get n from stack
        movl nb310nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb310nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        movq  nb310nf_vctot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb310nf_Vc(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment vc[gid] 

        movq  nb310nf_Vvdwtot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb310nf_Vvdw(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment Vvdw[gid] 

        ## finish if last 
        movl nb310nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel310nf_ia32_3dnow.nb310nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb310nf_n(%esp)
        jmp _nb_kernel310nf_ia32_3dnow.nb310nf_outer
_nb_kernel310nf_ia32_3dnow.nb310nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb310nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel310nf_ia32_3dnow.nb310nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel310nf_ia32_3dnow.nb310nf_threadloop
_nb_kernel310nf_ia32_3dnow.nb310nf_end: 
        femms

        movl nb310nf_nouter(%esp),%eax
        movl nb310nf_ninner(%esp),%ebx
        movl nb310nf_outeriter(%ebp),%ecx
        movl nb310nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        addl $116,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret


