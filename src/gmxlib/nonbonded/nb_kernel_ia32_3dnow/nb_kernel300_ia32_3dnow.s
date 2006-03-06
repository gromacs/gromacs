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




.globl nb_kernel300_ia32_3dnow
.globl _nb_kernel300_ia32_3dnow
nb_kernel300_ia32_3dnow:        
_nb_kernel300_ia32_3dnow:       
.set nb300_p_nri, 8
.set nb300_iinr, 12
.set nb300_jindex, 16
.set nb300_jjnr, 20
.set nb300_shift, 24
.set nb300_shiftvec, 28
.set nb300_fshift, 32
.set nb300_gid, 36
.set nb300_pos, 40
.set nb300_faction, 44
.set nb300_charge, 48
.set nb300_p_facel, 52
.set nb300_p_krf, 56
.set nb300_p_crf, 60
.set nb300_Vc, 64
.set nb300_type, 68
.set nb300_p_ntype, 72
.set nb300_vdwparam, 76
.set nb300_Vvdw, 80
.set nb300_p_tabscale, 84
.set nb300_VFtab, 88
.set nb300_invsqrta, 92
.set nb300_dvda, 96
.set nb300_p_gbtabscale, 100
.set nb300_GBtab, 104
.set nb300_p_nthreads, 108
.set nb300_count, 112
.set nb300_mtx, 116
.set nb300_outeriter, 120
.set nb300_inneriter, 124
.set nb300_work, 128
        ## stack offsets for local variables 
.set nb300_is3, 0
.set nb300_ii3, 4
.set nb300_ix, 8
.set nb300_iy, 12
.set nb300_iz, 16
.set nb300_iq, 20
.set nb300_vctot, 28
.set nb300_two, 36
.set nb300_n1, 44
.set nb300_tsc, 52
.set nb300_ntia, 60
.set nb300_innerjjnr, 64
.set nb300_innerk, 68
.set nb300_fix, 72
.set nb300_fiy, 76
.set nb300_fiz, 80
.set nb300_dx1, 84
.set nb300_dy1, 88
.set nb300_dz1, 92
.set nb300_dx2, 96
.set nb300_dy2, 100
.set nb300_dz2, 104
.set nb300_n, 108                           ## idx for outer loop
.set nb300_nn1, 112                         ## number of outer iterations
.set nb300_nri, 116
.set nb300_facel, 120
.set nb300_nouter, 124
.set nb300_ninner, 128
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $132,%esp          ## local stack space 
        femms
        ## move data to local stack  
        movl nb300_p_nri(%ebp),%ecx
        movl nb300_p_facel(%ebp),%esi
        movl nb300_p_tabscale(%ebp),%edi
        movl (%ecx),%ecx
        movl (%esi),%esi
        movl %ecx,nb300_nri(%esp)
        movl %esi,nb300_facel(%esp)

        movd  (%edi),%mm3
        punpckldq %mm3,%mm3
        movq  %mm3,nb300_tsc(%esp)
        movl $0x40000000,%eax
        movl %eax,nb300_two(%esp)
        movl %eax,nb300_two+4(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb300_nouter(%esp)
        movl %eax,nb300_ninner(%esp)

_nb_kernel300_ia32_3dnow.nb300_threadloop: 
        movl  nb300_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel300_ia32_3dnow.nb300_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel300_ia32_3dnow.nb300_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb300_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb300_n(%esp)
        movl %ebx,nb300_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel300_ia32_3dnow.nb300_outerstart
        jmp _nb_kernel300_ia32_3dnow.nb300_end

_nb_kernel300_ia32_3dnow.nb300_outerstart: 
        ## ebx contains number of outer iterations
        addl nb300_nouter(%esp),%ebx
        movl %ebx,nb300_nouter(%esp)

_nb_kernel300_ia32_3dnow.nb300_outer: 
        movl  nb300_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb300_is3(%esp)      ## store is3 

        movl  nb300_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movq  (%eax,%ebx,4),%mm0        ## move shX/shY to mm0 and shZ to mm1 
        movd  8(%eax,%ebx,4),%mm1

        movl  nb300_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx=ii 

        movl  nb300_charge(%ebp),%edx
        movd  (%edx,%ebx,4),%mm2    ## mm2=charge[ii] 
        pfmul nb300_facel(%esp),%mm2
        punpckldq %mm2,%mm2         ## spread to both halves 
        movq  %mm2,nb300_iq(%esp)           ## iq =facel*charge[ii] 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb300_pos(%ebp),%eax      ## eax = base of pos[] 

        pfadd (%eax,%ebx,4),%mm0    ## ix = shX + posX (and iy too) 
        movd  8(%eax,%ebx,4),%mm3       ## cant use direct memory add for 4 bytes (iz) 
        movl  %ebx,nb300_ii3(%esp)
        pfadd %mm3,%mm1
        movq  %mm0,nb300_ix(%esp)
        movd  %mm1,nb300_iz(%esp)

        ## clear total potential and i forces 
        pxor  %mm7,%mm7
        movq  %mm7,nb300_vctot(%esp)
        movq  %mm7,nb300_fix(%esp)
        movd  %mm7,nb300_fiz(%esp)

        movl  nb300_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb300_pos(%ebp),%esi
        movl  nb300_faction(%ebp),%edi
        movl  nb300_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb300_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb300_ninner(%esp),%ecx
        movl  %ecx,nb300_ninner(%esp)
        movl  %edx,nb300_innerk(%esp)      ## number of innerloop atoms 
        addl  $0,%edx
        jge   _nb_kernel300_ia32_3dnow.nb300_unroll_loop
        jmp   _nb_kernel300_ia32_3dnow.nb300_finish_inner
_nb_kernel300_ia32_3dnow.nb300_unroll_loop: 
        ## paired innerloop starts here 
        movl  nb300_innerjjnr(%esp),%ecx       ## pointer to jjnr[k] 
        movl  (%ecx),%eax
        movl  4(%ecx),%ebx           ## eax/ebx=jnr 
        addl $8,nb300_innerjjnr(%esp)             ## advance pointer (unrolled 2) 
        prefetch 16(%ecx)            ## prefetch data - trial and error says 16 is best 

        movl nb300_charge(%ebp),%ecx     ## base of charge[] 
        movq nb300_iq(%esp),%mm5
        movd (%ecx,%eax,4),%mm3      ## charge[jnr1] 
    punpckldq (%ecx,%ebx,4),%mm3     ## move charge 2 to high part of mm3 
        pfmul %mm5,%mm3              ## mm3 now has qq for both particles 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        movl  nb300_pos(%ebp),%esi

        movq  nb300_ix(%esp),%mm0
        movd  nb300_iz(%esp),%mm1
        movq  (%esi,%eax,4),%mm4     ## fetch first j coordinates 
        movd  8(%esi,%eax,4),%mm5
        pfsubr %mm0,%mm4             ## dr = ir - jr  
        pfsubr %mm1,%mm5
        movq  %mm4,nb300_dx1(%esp)           ## store dr 
        movd  %mm5,nb300_dz1(%esp)
        pfmul %mm4,%mm4              ## square dx,dy,dz                          
        pfmul %mm5,%mm5
        pfacc %mm5,%mm4              ## accumulate to get dx*dx+ dy*dy+ dz*dz 
        pfacc %mm5,%mm4              ## first rsq in lower mm4 

        movq  (%esi,%ebx,4),%mm6     ## fetch second j coordinates  
        movd  8(%esi,%ebx,4),%mm7

        pfsubr %mm0,%mm6             ## dr = ir - jr  
        pfsubr %mm1,%mm7
        movq  %mm6,nb300_dx2(%esp)           ## store dr 
        movd  %mm7,nb300_dz2(%esp)
        pfmul %mm6,%mm6              ## square dx,dy,dz 
        pfmul %mm7,%mm7
        pfacc %mm7,%mm6              ## accumulate to get dx*dx+ dy*dy+ dz*dz 
        pfacc %mm7,%mm6              ## second rsq in lower mm6 

    pfrsqrt %mm4,%mm0        ## lookup inverse square root seed 
    pfrsqrt %mm6,%mm1


        punpckldq %mm1,%mm0
        punpckldq %mm6,%mm4             ## now 4 has rsq and 0 the seed for both pairs. 
    movq %mm0,%mm2                      ## amd 3dnow N-R iteration to get full precision. 
        pfmul %mm0,%mm0
    pfrsqit1 %mm4,%mm0
    pfrcpit2 %mm2,%mm0
        pfmul %mm0,%mm4
        movq %mm4,%mm1
        ## mm0 is invsqrt, and mm1 r. 
        ## do potential and fscal 
        pfmul nb300_tsc(%esp),%mm1      ## mm1=rt 
        pf2iw %mm1,%mm4
        movq %mm4,nb300_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm1              ## now mm1 is eps and mm4 is n0 

        movq %mm1,%mm2
        pfmul %mm2,%mm2 ## mm1 is eps, mm2 is eps2 

        movl nb300_VFtab(%ebp),%edx
        movl nb300_n1(%esp),%ecx
        shll $2,%ecx
        ## coulomb table 
        ## load all the table values we need 
        movd (%edx,%ecx,4),%mm4
        movd 4(%edx,%ecx,4),%mm5
        movd 8(%edx,%ecx,4),%mm6
        movd 12(%edx,%ecx,4),%mm7
        movl nb300_n1+4(%esp),%ecx
        shll $2,%ecx
        punpckldq (%edx,%ecx,4),%mm4
        punpckldq 4(%edx,%ecx,4),%mm5
        punpckldq 8(%edx,%ecx,4),%mm6
        punpckldq 12(%edx,%ecx,4),%mm7

        pfmul %mm1,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 

        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 

        pfmul nb300_two(%esp),%mm7      ## two*Heps2 
        pfadd %mm6,%mm7
        pfadd %mm5,%mm7 ## mm7=FF 

        pfmul %mm1,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        pfmul %mm3,%mm5 ## vcoul=qq*VV 
        pfmul %mm7,%mm3 ## fijC=FF*qq  

        ## at this point mm5 contains vcoul and mm3 fijC. 
        ## increment vcoul - then we can get rid of mm5. 
        ## update vctot 
        pfadd nb300_vctot(%esp),%mm5        ## add the earlier value 
        movq %mm5,nb300_vctot(%esp)         ## store the sum       

        ## change sign of mm3 
    pxor %mm1,%mm1
        pfsub %mm3,%mm1
        pfmul nb300_tsc(%esp),%mm1
        pfmul %mm1,%mm0   ## mm0 is total fscal now     

        prefetchw nb300_dx1(%esp)       ## prefetch i forces to cache 

        ## spread fscalar to both positions 
        movq %mm0,%mm1
        punpckldq %mm0,%mm0
        punpckhdq %mm1,%mm1

        ## calc vector force 
        prefetchw (%edi,%eax,4) ## prefetch the 1st faction to cache 
        movq nb300_dx1(%esp),%mm2       ## fetch dr 
        movd nb300_dz1(%esp),%mm3

        prefetchw (%edi,%ebx,4) ## prefetch the 2nd faction to cache 
        pfmul %mm0,%mm2         ## mult by fs  
        pfmul %mm0,%mm3

        movq nb300_dx2(%esp),%mm4       ## fetch dr 
        movd nb300_dz2(%esp),%mm5
        pfmul %mm1,%mm4         ## mult by fs  
        pfmul %mm1,%mm5
        ## update i forces 

        movq nb300_fix(%esp),%mm0
        movd nb300_fiz(%esp),%mm1
        pfadd %mm2,%mm0
        pfadd %mm3,%mm1

        pfadd %mm4,%mm0
        pfadd %mm5,%mm1
        movq %mm0,nb300_fix(%esp)
        movd %mm1,nb300_fiz(%esp)
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
        subl $2,nb300_innerk(%esp)
        jl    _nb_kernel300_ia32_3dnow.nb300_finish_inner
        jmp   _nb_kernel300_ia32_3dnow.nb300_unroll_loop
_nb_kernel300_ia32_3dnow.nb300_finish_inner: 
        andl $1,nb300_innerk(%esp)
        jnz  _nb_kernel300_ia32_3dnow.nb300_single_inner
        jmp  _nb_kernel300_ia32_3dnow.nb300_updateouterdata
_nb_kernel300_ia32_3dnow.nb300_single_inner: 
        ## a single j particle iteration here - compare with the unrolled code for comments. 
        movl  nb300_innerjjnr(%esp),%eax
        movl  (%eax),%eax       ## eax=jnr offset 

        movl nb300_charge(%ebp),%ecx
        movd nb300_iq(%esp),%mm5
        movd (%ecx,%eax,4),%mm3
        pfmul %mm5,%mm3         ## mm3=qq 

        movl  nb300_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        movq  nb300_ix(%esp),%mm0
        movd  nb300_iz(%esp),%mm1
        movq  (%esi,%eax,4),%mm4
        movd  8(%esi,%eax,4),%mm5
        pfsubr %mm0,%mm4
        pfsubr %mm1,%mm5
        movq  %mm4,nb300_dx1(%esp)
        pfmul %mm4,%mm4
        movd  %mm5,nb300_dz1(%esp)
        pfmul %mm5,%mm5
        pfacc %mm5,%mm4
        pfacc %mm5,%mm4         ## mm0=rsq 

    pfrsqrt %mm4,%mm0
    movq %mm0,%mm2
    pfmul %mm0,%mm0
    pfrsqit1 %mm4,%mm0
    pfrcpit2 %mm2,%mm0  ## mm1=invsqrt 
        pfmul %mm0,%mm4
        movq %mm4,%mm1
        ## mm0 is invsqrt, and mm1 r. 

        ## calculate potentials and scalar force 
        pfmul nb300_tsc(%esp),%mm1      ## mm1=rt 
        pf2iw %mm1,%mm4
        movd %mm4,nb300_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm1              ## now mm1 is eps and mm4 is n0 

        movq %mm1,%mm2
        pfmul %mm2,%mm2 ## mm1 is eps, mm2 is eps2 

        ## coulomb table 
        movl nb300_VFtab(%ebp),%edx
        movl nb300_n1(%esp),%ecx
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

        pfmul nb300_two(%esp),%mm7      ## two*Heps2 
        pfadd %mm6,%mm7
        pfadd %mm5,%mm7 ## mm7=FF 

        pfmul %mm1,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        pfmul %mm3,%mm5 ## vcoul=qq*VV 
        pfmul %mm7,%mm3 ## fijC=FF*qq  

        ## at this point mm5 contains vcoul and mm3 fijC 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        pfadd nb300_vctot(%esp),%mm5        ## add the earlier value 
        movq %mm5,nb300_vctot(%esp)         ## store the sum       

        ## change sign of mm3 
    pxor %mm1,%mm1
        pfsub %mm3,%mm1
        pfmul nb300_tsc(%esp),%mm0
        pfmul %mm1,%mm0   ## mm0 is total fscal now     

        ## spread fscalar to both positions 
        punpckldq %mm0,%mm0
        ## calc vectorial force 
        prefetchw (%edi,%eax,4) ## prefetch faction to cache  
        movq nb300_dx1(%esp),%mm2
        movd nb300_dz1(%esp),%mm3


        pfmul %mm0,%mm2
        pfmul %mm0,%mm3

        ## update i particle force 
        movq nb300_fix(%esp),%mm0
        movd nb300_fiz(%esp),%mm1
        pfadd %mm2,%mm0
        pfadd %mm3,%mm1
        movq %mm0,nb300_fix(%esp)
        movd %mm1,nb300_fiz(%esp)
        ## update j particle force 
        movq (%edi,%eax,4),%mm0
        movd 8(%edi,%eax,4),%mm1
        pfsub %mm2,%mm0
        pfsub %mm3,%mm1
        movq %mm0,(%edi,%eax,4)
        movd %mm1,8(%edi,%eax,4)
        ## done! 
_nb_kernel300_ia32_3dnow.nb300_updateouterdata: 
        movl  nb300_ii3(%esp),%ecx

        movq  (%edi,%ecx,4),%mm6       ## increment i force 
        movd  8(%edi,%ecx,4),%mm7
        pfadd nb300_fix(%esp),%mm6
        pfadd nb300_fiz(%esp),%mm7
        movq  %mm6,(%edi,%ecx,4)
        movd  %mm7,8(%edi,%ecx,4)

        movl  nb300_fshift(%ebp),%ebx      ## increment fshift force 
        movl  nb300_is3(%esp),%edx

        movq  (%ebx,%edx,4),%mm6
        movd  8(%ebx,%edx,4),%mm7
        pfadd nb300_fix(%esp),%mm6
        pfadd nb300_fiz(%esp),%mm7
        movq  %mm6,(%ebx,%edx,4)
        movd  %mm7,8(%ebx,%edx,4)

        ## get n from stack
        movl nb300_n(%esp),%esi
        ## get group index for i particle 
        movl  nb300_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        movq  nb300_vctot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb300_Vc(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment vc[gid] 

        ## finish if last 
        movl nb300_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel300_ia32_3dnow.nb300_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb300_n(%esp)
        jmp _nb_kernel300_ia32_3dnow.nb300_outer
_nb_kernel300_ia32_3dnow.nb300_outerend: 
        ## check if more outer neighborlists remain
        movl  nb300_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel300_ia32_3dnow.nb300_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel300_ia32_3dnow.nb300_threadloop
_nb_kernel300_ia32_3dnow.nb300_end: 
        femms

        movl nb300_nouter(%esp),%eax
        movl nb300_ninner(%esp),%ebx
        movl nb300_outeriter(%ebp),%ecx
        movl nb300_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        addl $132,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret



.globl nb_kernel300nf_ia32_3dnow
.globl _nb_kernel300nf_ia32_3dnow
nb_kernel300nf_ia32_3dnow:      
_nb_kernel300nf_ia32_3dnow:     
.set nb300nf_p_nri, 8
.set nb300nf_iinr, 12
.set nb300nf_jindex, 16
.set nb300nf_jjnr, 20
.set nb300nf_shift, 24
.set nb300nf_shiftvec, 28
.set nb300nf_fshift, 32
.set nb300nf_gid, 36
.set nb300nf_pos, 40
.set nb300nf_faction, 44
.set nb300nf_charge, 48
.set nb300nf_p_facel, 52
.set nb300nf_p_krf, 56
.set nb300nf_p_crf, 60
.set nb300nf_Vc, 64
.set nb300nf_type, 68
.set nb300nf_p_ntype, 72
.set nb300nf_vdwparam, 76
.set nb300nf_Vvdw, 80
.set nb300nf_p_tabscale, 84
.set nb300nf_VFtab, 88
.set nb300nf_invsqrta, 92
.set nb300nf_dvda, 96
.set nb300nf_p_gbtabscale, 100
.set nb300nf_GBtab, 104
.set nb300nf_p_nthreads, 108
.set nb300nf_count, 112
.set nb300nf_mtx, 116
.set nb300nf_outeriter, 120
.set nb300nf_inneriter, 124
.set nb300nf_work, 128
        ## stack offsets for local variables 
.set nb300nf_is3, 0
.set nb300nf_ii3, 4
.set nb300nf_ix, 8
.set nb300nf_iy, 12
.set nb300nf_iz, 16
.set nb300nf_iq, 20
.set nb300nf_vctot, 28
.set nb300nf_n1, 36
.set nb300nf_tsc, 44
.set nb300nf_ntia, 52
.set nb300nf_innerjjnr, 56
.set nb300nf_innerk, 60
.set nb300nf_n, 64                         ## idx for outer loop
.set nb300nf_nn1, 68                       ## number of outer iterations
.set nb300nf_nri, 72
.set nb300nf_facel, 76
.set nb300nf_nouter, 80
.set nb300nf_ninner, 84
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $88,%esp           ## local stack space 
        femms
        ## move data to local stack  
        movl nb300nf_p_nri(%ebp),%ecx
        movl nb300nf_p_facel(%ebp),%esi
        movl nb300nf_p_tabscale(%ebp),%edi
        movl (%ecx),%ecx
        movl (%esi),%esi
        movl %ecx,nb300nf_nri(%esp)
        movl %esi,nb300nf_facel(%esp)

        movd  (%edi),%mm3
        punpckldq %mm3,%mm3
        movq  %mm3,nb300nf_tsc(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb300nf_nouter(%esp)
        movl %eax,nb300nf_ninner(%esp)

_nb_kernel300nf_ia32_3dnow.nb300nf_threadloop: 
        movl  nb300nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel300nf_ia32_3dnow.nb300nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel300nf_ia32_3dnow.nb300nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb300nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb300nf_n(%esp)
        movl %ebx,nb300nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel300nf_ia32_3dnow.nb300nf_outerstart
        jmp _nb_kernel300nf_ia32_3dnow.nb300nf_end

_nb_kernel300nf_ia32_3dnow.nb300nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb300nf_nouter(%esp),%ebx
        movl %ebx,nb300nf_nouter(%esp)

_nb_kernel300nf_ia32_3dnow.nb300nf_outer: 
        movl  nb300nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb300nf_is3(%esp)            ## store is3 

        movl  nb300nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movq  (%eax,%ebx,4),%mm0        ## move shX/shY to mm0 and shZ to mm1 
        movd  8(%eax,%ebx,4),%mm1

        movl  nb300nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx            ## ebx=ii 

        movl  nb300nf_charge(%ebp),%edx
        movd  (%edx,%ebx,4),%mm2    ## mm2=charge[ii] 
        pfmul nb300nf_facel(%esp),%mm2
        punpckldq %mm2,%mm2         ## spread to both halves 
        movq  %mm2,nb300nf_iq(%esp)         ## iq =facel*charge[ii] 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb300nf_pos(%ebp),%eax      ## eax = base of pos[] 

        pfadd (%eax,%ebx,4),%mm0    ## ix = shX + posX (and iy too) 
        movd  8(%eax,%ebx,4),%mm3       ## cant use direct memory add for 4 bytes (iz) 
        movl  %ebx,nb300nf_ii3(%esp)
        pfadd %mm3,%mm1
        movq  %mm0,nb300nf_ix(%esp)
        movd  %mm1,nb300nf_iz(%esp)

        ## clear total potential and i forces 
        pxor  %mm7,%mm7
        movq  %mm7,nb300nf_vctot(%esp)

        movl  nb300nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb300nf_pos(%ebp),%esi
        movl  nb300nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb300nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb300nf_ninner(%esp),%ecx
        movl  %ecx,nb300nf_ninner(%esp)
        movl  %edx,nb300nf_innerk(%esp)      ## number of innerloop atoms 
        addl  $0,%edx
        jge   _nb_kernel300nf_ia32_3dnow.nb300nf_unroll_loop
        jmp   _nb_kernel300nf_ia32_3dnow.nb300nf_finish_inner
_nb_kernel300nf_ia32_3dnow.nb300nf_unroll_loop: 
        ## paired innerloop starts here 
        movl  nb300nf_innerjjnr(%esp),%ecx       ## pointer to jjnr[k] 
        movl  (%ecx),%eax
        movl  4(%ecx),%ebx           ## eax/ebx=jnr 
        addl $8,nb300nf_innerjjnr(%esp)             ## advance pointer (unrolled 2) 
        prefetch 16(%ecx)            ## prefetch data - trial and error says 16 is best 

        movl nb300nf_charge(%ebp),%ecx     ## base of charge[] 
        movq nb300nf_iq(%esp),%mm5
        movd (%ecx,%eax,4),%mm3      ## charge[jnr1] 
    punpckldq (%ecx,%ebx,4),%mm3     ## move charge 2 to high part of mm3 
        pfmul %mm5,%mm3              ## mm3 now has qq for both particles 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        movl  nb300nf_pos(%ebp),%esi

        movq  nb300nf_ix(%esp),%mm0
        movd  nb300nf_iz(%esp),%mm1
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

    pfrsqrt %mm4,%mm0        ## lookup inverse square root seed 
    pfrsqrt %mm6,%mm1


        punpckldq %mm1,%mm0
        punpckldq %mm6,%mm4             ## now 4 has rsq and 0 the seed for both pairs. 
    movq %mm0,%mm2                      ## amd 3dnow N-R iteration to get full precision. 
        pfmul %mm0,%mm0
    pfrsqit1 %mm4,%mm0
    pfrcpit2 %mm2,%mm0
        pfmul %mm0,%mm4
        movq %mm4,%mm1
        ## mm0 is invsqrt, and mm1 r. 
        ## do potential and fscal 
        pfmul nb300nf_tsc(%esp),%mm1    ## mm1=rt 
        pf2iw %mm1,%mm4
        movq %mm4,nb300nf_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm1              ## now mm1 is eps and mm4 is n0 

        movq %mm1,%mm2
        pfmul %mm2,%mm2 ## mm1 is eps, mm2 is eps2 

        movl nb300nf_VFtab(%ebp),%edx
        movl nb300nf_n1(%esp),%ecx
        shll $2,%ecx
        ## coulomb table 
        ## load all the table values we need 
        movd (%edx,%ecx,4),%mm4
        movd 4(%edx,%ecx,4),%mm5
        movd 8(%edx,%ecx,4),%mm6
        movd 12(%edx,%ecx,4),%mm7
        movl nb300nf_n1+4(%esp),%ecx
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
        ## at this point mm5 contains vcoul 
        ## increment vcoul - then we can get rid of mm5. 
        ## update vctot 
        pfadd nb300nf_vctot(%esp),%mm5        ## add the earlier value 
        movq %mm5,nb300nf_vctot(%esp)         ## store the sum       

        ## should we do one more iteration? 
        subl $2,nb300nf_innerk(%esp)
        jl    _nb_kernel300nf_ia32_3dnow.nb300nf_finish_inner
        jmp   _nb_kernel300nf_ia32_3dnow.nb300nf_unroll_loop
_nb_kernel300nf_ia32_3dnow.nb300nf_finish_inner: 
        andl $1,nb300nf_innerk(%esp)
        jnz  _nb_kernel300nf_ia32_3dnow.nb300nf_single_inner
        jmp  _nb_kernel300nf_ia32_3dnow.nb300nf_updateouterdata
_nb_kernel300nf_ia32_3dnow.nb300nf_single_inner: 
        ## a single j particle iteration here - compare with the unrolled code for comments. 
        movl  nb300nf_innerjjnr(%esp),%eax
        movl  (%eax),%eax       ## eax=jnr offset 

        movl nb300nf_charge(%ebp),%ecx
        movd nb300nf_iq(%esp),%mm5
        movd (%ecx,%eax,4),%mm3
        pfmul %mm5,%mm3         ## mm3=qq 

        movl  nb300nf_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        movq  nb300nf_ix(%esp),%mm0
        movd  nb300nf_iz(%esp),%mm1
        movq  (%esi,%eax,4),%mm4
        movd  8(%esi,%eax,4),%mm5
        pfsubr %mm0,%mm4
        pfsubr %mm1,%mm5
        pfmul %mm4,%mm4
        pfmul %mm5,%mm5
        pfacc %mm5,%mm4
        pfacc %mm5,%mm4         ## mm0=rsq 

    pfrsqrt %mm4,%mm0
    movq %mm0,%mm2
    pfmul %mm0,%mm0
    pfrsqit1 %mm4,%mm0
    pfrcpit2 %mm2,%mm0  ## mm1=invsqrt 
        pfmul %mm0,%mm4
        movq %mm4,%mm1
        ## mm0 is invsqrt, and mm1 r. 

        ## calculate potentials and scalar force 
        pfmul nb300nf_tsc(%esp),%mm1    ## mm1=rt 
        pf2iw %mm1,%mm4
        movd %mm4,nb300nf_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm1              ## now mm1 is eps and mm4 is n0 

        movq %mm1,%mm2
        pfmul %mm2,%mm2 ## mm1 is eps, mm2 is eps2 

        ## coulomb table 
        movl nb300nf_VFtab(%ebp),%edx
        movl nb300nf_n1(%esp),%ecx
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
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        pfadd nb300nf_vctot(%esp),%mm5        ## add the earlier value 
        movq %mm5,nb300nf_vctot(%esp)         ## store the sum       

_nb_kernel300nf_ia32_3dnow.nb300nf_updateouterdata: 
        ## get n from stack
        movl nb300nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb300nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        movq  nb300nf_vctot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb300nf_Vc(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment vc[gid] 

        ## finish if last 
        movl nb300nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel300nf_ia32_3dnow.nb300nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb300nf_n(%esp)
        jmp _nb_kernel300nf_ia32_3dnow.nb300nf_outer
_nb_kernel300nf_ia32_3dnow.nb300nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb300nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel300nf_ia32_3dnow.nb300nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel300nf_ia32_3dnow.nb300nf_threadloop
_nb_kernel300nf_ia32_3dnow.nb300nf_end: 
        femms

        movl nb300nf_nouter(%esp),%eax
        movl nb300nf_ninner(%esp),%ebx
        movl nb300nf_outeriter(%ebp),%ecx
        movl nb300nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        addl $88,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret




