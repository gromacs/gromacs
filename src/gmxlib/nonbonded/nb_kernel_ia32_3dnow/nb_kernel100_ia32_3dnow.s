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


.globl nb_kernel100_ia32_3dnow
.globl _nb_kernel100_ia32_3dnow
nb_kernel100_ia32_3dnow:        
_nb_kernel100_ia32_3dnow:       
.set nb100_p_nri, 8
.set nb100_iinr, 12
.set nb100_jindex, 16
.set nb100_jjnr, 20
.set nb100_shift, 24
.set nb100_shiftvec, 28
.set nb100_fshift, 32
.set nb100_gid, 36
.set nb100_pos, 40
.set nb100_faction, 44
.set nb100_charge, 48
.set nb100_p_facel, 52
.set nb100_p_krf, 56
.set nb100_p_crf, 60
.set nb100_Vc, 64
.set nb100_type, 68
.set nb100_p_ntype, 72
.set nb100_vdwparam, 76
.set nb100_Vvdw, 80
.set nb100_p_tabscale, 84
.set nb100_VFtab, 88
.set nb100_invsqrta, 92
.set nb100_dvda, 96
.set nb100_p_gbtabscale, 100
.set nb100_GBtab, 104
.set nb100_p_nthreads, 108
.set nb100_count, 112
.set nb100_mtx, 116
.set nb100_outeriter, 120
.set nb100_inneriter, 124
.set nb100_work, 128
        ## stack offsets for local variables 
.set nb100_is3, 0
.set nb100_ii3, 4
.set nb100_ix, 8
.set nb100_iy, 12
.set nb100_iz, 16
.set nb100_iq, 20
.set nb100_vctot, 28
.set nb100_innerjjnr, 36
.set nb100_innerk, 40
.set nb100_fix, 44
.set nb100_fiy, 48
.set nb100_fiz, 52
.set nb100_dx1, 56
.set nb100_dy1, 60
.set nb100_dz1, 64
.set nb100_dx2, 68
.set nb100_dy2, 72
.set nb100_dz2, 76
.set nb100_n, 80                           ## idx for outer loop
.set nb100_nn1, 84                         ## number of outer iterations
.set nb100_nri, 88
.set nb100_facel, 92
.set nb100_nouter, 96
.set nb100_ninner, 100
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $104,%esp          ## local stack space 
        femms

        movl nb100_p_nri(%ebp),%ecx
        movl nb100_p_facel(%ebp),%esi
        movl (%ecx),%ecx
        movl (%esi),%esi
        movl %ecx,nb100_nri(%esp)
        movl %esi,nb100_facel(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb100_nouter(%esp)
        movl %eax,nb100_ninner(%esp)

_nb_kernel100_ia32_3dnow.nb100_threadloop: 
        movl  nb100_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel100_ia32_3dnow.nb100_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel100_ia32_3dnow.nb100_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb100_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb100_n(%esp)
        movl %ebx,nb100_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi

        jg  _nb_kernel100_ia32_3dnow.nb100_outerstart
        jmp _nb_kernel100_ia32_3dnow.nb100_end

_nb_kernel100_ia32_3dnow.nb100_outerstart: 
        ## ebx contains number of outer iterations
        addl nb100_nouter(%esp),%ebx
        movl %ebx,nb100_nouter(%esp)

_nb_kernel100_ia32_3dnow.nb100_outer: 
        movl  nb100_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb100_is3(%esp)      ## store is3 

        movl  nb100_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movq  (%eax,%ebx,4),%mm0        ## move shX/shY to mm0 and shZ to mm1 
        movd  8(%eax,%ebx,4),%mm1

        movl  nb100_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx            ## ebx=ii 

        movl  nb100_charge(%ebp),%edx
        movd  (%edx,%ebx,4),%mm2    ## mm2=charge[ii] 
        pfmul nb100_facel(%esp),%mm2
        punpckldq %mm2,%mm2         ## spread to both halves 
        movq  %mm2,nb100_iq(%esp)           ## iq =facel*charge[ii] 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb100_pos(%ebp),%eax      ## eax = base of pos[] 

        pfadd (%eax,%ebx,4),%mm0    ## ix = shX + posX (and iy too) 
        movd  8(%eax,%ebx,4),%mm3       ## cant use direct memory add for 4 bytes (iz) 
        movl  %ebx,nb100_ii3(%esp)
        pfadd %mm3,%mm1
        movq  %mm0,nb100_ix(%esp)
        movd  %mm1,nb100_iz(%esp)

        ## clear vctot and i forces 
        pxor  %mm7,%mm7
        movq  %mm7,nb100_vctot(%esp)
        movq  %mm7,nb100_fix(%esp)
        movd  %mm7,nb100_fiz(%esp)

        movl  nb100_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb100_pos(%ebp),%esi
        movl  nb100_faction(%ebp),%edi
        movl  nb100_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb100_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb100_ninner(%esp),%ecx
        movl  %ecx,nb100_ninner(%esp)
        movl  %edx,nb100_innerk(%esp)      ## number of innerloop atoms 
        addl  $0,%edx
        jge   _nb_kernel100_ia32_3dnow.nb100_unroll_loop
        jmp   _nb_kernel100_ia32_3dnow.nb100_finish_inner
_nb_kernel100_ia32_3dnow.nb100_unroll_loop: 
        ## paired innerloop starts here 
        movl  nb100_innerjjnr(%esp),%ecx       ## pointer to jjnr[k] 
        movl  (%ecx),%eax
        movl  4(%ecx),%ebx           ## eax/ebx=jnr 
        addl $8,nb100_innerjjnr(%esp)             ## advance pointer (unrolled 2) 
        prefetch 16(%ecx)            ## prefetch data - trial and error says 16 is best 

        movl nb100_charge(%ebp),%ecx     ## base of charge[] 
        movq nb100_iq(%esp),%mm5
        movd (%ecx,%eax,4),%mm3      ## charge[jnr1] 
        movd (%ecx,%ebx,4),%mm7          ## charge[jnr2] 
        punpckldq %mm7,%mm3          ## move charge 2 to high part of mm3 
        pfmul %mm5,%mm3              ## mm3 now has qq for both particles 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        movq  nb100_ix(%esp),%mm0
        movd  nb100_iz(%esp),%mm1
        movq  (%esi,%eax,4),%mm4     ## fetch first j coordinates 
        movd  8(%esi,%eax,4),%mm5
        pfsubr %mm0,%mm4             ## dr = ir - jr  
        pfsubr %mm1,%mm5
        movq  %mm4,nb100_dx1(%esp)           ## store dr 
        movd  %mm5,nb100_dz1(%esp)
        pfmul %mm4,%mm4              ## square dx,dy,dz                          
        pfmul %mm5,%mm5
        pfacc %mm5,%mm4              ## accumulate to get dx*dx+ dy*dy+ dz*dz 
        pfacc %mm5,%mm4              ## first rsq in lower mm4 

        movq  (%esi,%ebx,4),%mm6     ## fetch second j coordinates  
        movd  8(%esi,%ebx,4),%mm7

        pfsubr %mm0,%mm6             ## dr = ir - jr  
        pfsubr %mm1,%mm7
        movq  %mm6,nb100_dx2(%esp)           ## store dr 
        movd  %mm7,nb100_dz2(%esp)
        pfmul %mm6,%mm6              ## square dx,dy,dz 
        pfmul %mm7,%mm7
        pfacc %mm7,%mm6              ## accumulate to get dx*dx+ dy*dy+ dz*dz 
        pfacc %mm7,%mm6              ## second rsq in lower mm6 

        pfrsqrt %mm4,%mm0            ## lookup inverse square root seed 
        pfrsqrt %mm6,%mm1

        punpckldq %mm1,%mm0
        punpckldq %mm6,%mm4             ## now 4 has rsq and 0 the seed for both pairs 
        movq %mm0,%mm2                  ## amd 3dnow N-R iteration to get full precision 
        pfmul %mm0,%mm0
        pfrsqit1 %mm4,%mm0
        pfrcpit2 %mm2,%mm0
        movq %mm0,%mm1
        pfmul %mm0,%mm0
        ## mm0 now contains invsq, and mm1 invsqrt
         ## do potential and fscal

        prefetchw nb100_dx1(%esp)       ## prefetch i forces to cache 

        pfmul %mm1,%mm3         ## mm3 has both vcoul 
        pfmul %mm3,%mm0         ## mm0 has both fscal 

        ## update vctot 

        pfadd nb100_vctot(%esp),%mm3        ## add the earlier value  
        movq %mm3,nb100_vctot(%esp)         ## store the sum 
        ## spread fscalar to both positions 
        movq %mm0,%mm1
        punpckldq %mm0,%mm0
        punpckhdq %mm1,%mm1
        ## calc vector force 
        prefetchw (%edi,%eax,4) ## prefetch the 1st faction to cache 
        movq nb100_dx1(%esp),%mm2       ## fetch dr 
        movd nb100_dz1(%esp),%mm3
        prefetchw (%edi,%ebx,4) ## prefetch the 2nd faction to cache 
        pfmul %mm0,%mm2         ## mult by fs  
        pfmul %mm0,%mm3

        movq nb100_dx2(%esp),%mm4       ## fetch dr 
        movd nb100_dz2(%esp),%mm5
        pfmul %mm1,%mm4         ## mult by fs  
        pfmul %mm1,%mm5
        ## update i forces 

        movq nb100_fix(%esp),%mm0
        movd nb100_fiz(%esp),%mm1
        pfadd %mm2,%mm0
        pfadd %mm3,%mm1

        pfadd %mm4,%mm0
        pfadd %mm5,%mm1
        movq %mm0,nb100_fix(%esp)
        movd %mm1,nb100_fiz(%esp)
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
        subl $2,nb100_innerk(%esp)
        jl    _nb_kernel100_ia32_3dnow.nb100_finish_inner
        jmp   _nb_kernel100_ia32_3dnow.nb100_unroll_loop
_nb_kernel100_ia32_3dnow.nb100_finish_inner: 
        andl $1,nb100_innerk(%esp)
        jnz  _nb_kernel100_ia32_3dnow.nb100_single_inner
        jmp  _nb_kernel100_ia32_3dnow.nb100_updateouterdata
_nb_kernel100_ia32_3dnow.nb100_single_inner: 
        ## a single j particle iteration here - compare with the unrolled code for comments 
        movl  nb100_innerjjnr(%esp),%eax
        movl  (%eax),%eax       ## eax=jnr offset 

        movl nb100_charge(%ebp),%ecx
        pxor %mm6,%mm6
        movd nb100_iq(%esp),%mm6
        movd (%ecx,%eax,4),%mm7
        pfmul %mm7,%mm6         ## mm6=qq 

        leal  (%eax,%eax,2),%eax

        movq  nb100_ix(%esp),%mm0
        movd  nb100_iz(%esp),%mm1
        movq  (%esi,%eax,4),%mm2
        movd  8(%esi,%eax,4),%mm3
        pfsub %mm2,%mm0
        pfsub %mm3,%mm1
        movq  %mm0,nb100_dx1(%esp)
        pfmul %mm0,%mm0
        movd  %mm1,nb100_dz1(%esp)
        pfmul %mm1,%mm1
        pfacc %mm1,%mm0
        pfacc %mm1,%mm0         ## mm0=rsq 

        pfrsqrt %mm0,%mm1
        movq %mm1,%mm2
        pfmul %mm1,%mm1
        pfrsqit1 %mm0,%mm1
        pfrcpit2 %mm2,%mm1      ## mm1=invsqrt 
        movq  %mm1,%mm4
        pfmul %mm4,%mm4         ## mm4=invsq 
        ## calculate potential and scalar force 
        pfmul %mm1,%mm6         ## mm6=vcoul 
        pfmul %mm6,%mm4         ## mm4=fscalar  

        ## update vctot 
        pfadd nb100_vctot(%esp),%mm6
        movq %mm6,nb100_vctot(%esp)
        ## spread fscalar to both positions 
        punpckldq %mm4,%mm4
        ## calc vectorial force 
        prefetchw (%edi,%eax,4) ## prefetch faction to cache  
        movq nb100_dx1(%esp),%mm0
        movd nb100_dz1(%esp),%mm1
        pfmul %mm4,%mm0
        pfmul %mm4,%mm1
        ## update i particle force 
        movq nb100_fix(%esp),%mm2
        movd nb100_fiz(%esp),%mm3
        pfadd %mm0,%mm2
        pfadd %mm1,%mm3
        movq %mm2,nb100_fix(%esp)
        movd %mm3,nb100_fiz(%esp)
        ## update j particle force 
        movq (%edi,%eax,4),%mm2
        movd 8(%edi,%eax,4),%mm3
        pfsub %mm0,%mm2
        pfsub %mm1,%mm3
        movq %mm2,(%edi,%eax,4)
        movd %mm3,8(%edi,%eax,4)
        ## done! 
_nb_kernel100_ia32_3dnow.nb100_updateouterdata: 
        movl  nb100_ii3(%esp),%ecx

        movq  (%edi,%ecx,4),%mm6       ## increment i force 
        movd  8(%edi,%ecx,4),%mm7
        pfadd nb100_fix(%esp),%mm6
        pfadd nb100_fiz(%esp),%mm7
        movq  %mm6,(%edi,%ecx,4)
        movd  %mm7,8(%edi,%ecx,4)

        movl  nb100_fshift(%ebp),%ebx      ## increment fshift force 
        movl  nb100_is3(%esp),%edx

        movq  (%ebx,%edx,4),%mm6
        movd  8(%ebx,%edx,4),%mm7
        pfadd nb100_fix(%esp),%mm6
        pfadd nb100_fiz(%esp),%mm7
        movq  %mm6,(%ebx,%edx,4)
        movd  %mm7,8(%ebx,%edx,4)

        ## get n from stack
        movl nb100_n(%esp),%esi
        ## get group index for i particle 
        movl  nb100_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        movq  nb100_vctot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb100_Vc(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment vc[gid] 
        ## finish if last 
        movl nb100_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel100_ia32_3dnow.nb100_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb100_n(%esp)
        jmp _nb_kernel100_ia32_3dnow.nb100_outer
_nb_kernel100_ia32_3dnow.nb100_outerend: 
        ## check if more outer neighborlists remain
        movl  nb100_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel100_ia32_3dnow.nb100_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel100_ia32_3dnow.nb100_threadloop
_nb_kernel100_ia32_3dnow.nb100_end: 
        femms

        movl nb100_nouter(%esp),%eax
        movl nb100_ninner(%esp),%ebx
        movl nb100_outeriter(%ebp),%ecx
        movl nb100_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        addl $104,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret










.globl nb_kernel100nf_ia32_3dnow
.globl _nb_kernel100nf_ia32_3dnow
nb_kernel100nf_ia32_3dnow:      
_nb_kernel100nf_ia32_3dnow:     
.set nb100nf_p_nri, 8
.set nb100nf_iinr, 12
.set nb100nf_jindex, 16
.set nb100nf_jjnr, 20
.set nb100nf_shift, 24
.set nb100nf_shiftvec, 28
.set nb100nf_fshift, 32
.set nb100nf_gid, 36
.set nb100nf_pos, 40
.set nb100nf_faction, 44
.set nb100nf_charge, 48
.set nb100nf_p_facel, 52
.set nb100nf_p_krf, 56
.set nb100nf_p_crf, 60
.set nb100nf_Vc, 64
.set nb100nf_type, 68
.set nb100nf_p_ntype, 72
.set nb100nf_vdwparam, 76
.set nb100nf_Vvdw, 80
.set nb100nf_p_tabscale, 84
.set nb100nf_VFtab, 88
.set nb100nf_invsqrta, 92
.set nb100nf_dvda, 96
.set nb100nf_p_gbtabscale, 100
.set nb100nf_GBtab, 104
.set nb100nf_p_nthreads, 108
.set nb100nf_count, 112
.set nb100nf_mtx, 116
.set nb100nf_outeriter, 120
.set nb100nf_inneriter, 124
.set nb100nf_work, 128
        ## stack offsets for local variables 
.set nb100nf_is3, 0
.set nb100nf_ii3, 4
.set nb100nf_ix, 8
.set nb100nf_iy, 12
.set nb100nf_iz, 16
.set nb100nf_iq, 20
.set nb100nf_vctot, 28
.set nb100nf_innerjjnr, 36
.set nb100nf_innerk, 40
.set nb100nf_n, 44                         ## idx for outer loop
.set nb100nf_nn1, 48                       ## number of outer iterations
.set nb100nf_nri, 52
.set nb100nf_facel, 56
.set nb100nf_nouter, 60
.set nb100nf_ninner, 64
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $68,%esp           ## local stack space 
        femms

        movl nb100nf_p_nri(%ebp),%ecx
        movl nb100nf_p_facel(%ebp),%esi
        movl (%ecx),%ecx
        movl (%esi),%esi
        movl %ecx,nb100nf_nri(%esp)
        movl %esi,nb100nf_facel(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb100nf_nouter(%esp)
        movl %eax,nb100nf_ninner(%esp)

_nb_kernel100nf_ia32_3dnow.nb100nf_threadloop: 
        movl  nb100nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel100nf_ia32_3dnow.nb100nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel100nf_ia32_3dnow.nb100nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb100nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb100nf_n(%esp)
        movl %ebx,nb100nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi

        jg  _nb_kernel100nf_ia32_3dnow.nb100nf_outerstart
        jmp _nb_kernel100nf_ia32_3dnow.nb100nf_end

_nb_kernel100nf_ia32_3dnow.nb100nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb100nf_nouter(%esp),%ebx
        movl %ebx,nb100nf_nouter(%esp)

_nb_kernel100nf_ia32_3dnow.nb100nf_outer: 
        movl  nb100nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb100nf_is3(%esp)            ## store is3 

        movl  nb100nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movq  (%eax,%ebx,4),%mm0        ## move shX/shY to mm0 and shZ to mm1 
        movd  8(%eax,%ebx,4),%mm1

        movl  nb100nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx            ## ebx=ii 

        movl  nb100nf_charge(%ebp),%edx
        movd  (%edx,%ebx,4),%mm2    ## mm2=charge[ii] 
        pfmul nb100nf_facel(%esp),%mm2
        punpckldq %mm2,%mm2         ## spread to both halves 
        movq  %mm2,nb100nf_iq(%esp)         ## iq =facel*charge[ii] 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb100nf_pos(%ebp),%eax      ## eax = base of pos[] 

        pfadd (%eax,%ebx,4),%mm0    ## ix = shX + posX (and iy too) 
        movd  8(%eax,%ebx,4),%mm3       ## cant use direct memory add for 4 bytes (iz) 
        movl  %ebx,nb100nf_ii3(%esp)
        pfadd %mm3,%mm1
        movq  %mm0,nb100nf_ix(%esp)
        movd  %mm1,nb100nf_iz(%esp)

        ## clear vctot
        pxor  %mm7,%mm7
        movq  %mm7,nb100nf_vctot(%esp)

        movl  nb100nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb100nf_pos(%ebp),%esi
        movl  nb100nf_faction(%ebp),%edi
        movl  nb100nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb100nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb100nf_ninner(%esp),%ecx
        movl  %ecx,nb100nf_ninner(%esp)
        movl  %edx,nb100nf_innerk(%esp)      ## number of innerloop atoms 
        addl  $0,%edx
        jge   _nb_kernel100nf_ia32_3dnow.nb100nf_unroll_loop
        jmp   _nb_kernel100nf_ia32_3dnow.nb100nf_finish_inner
_nb_kernel100nf_ia32_3dnow.nb100nf_unroll_loop: 
        ## paired innerloop starts here 
        movl  nb100nf_innerjjnr(%esp),%ecx       ## pointer to jjnr[k] 
        movl  (%ecx),%eax
        movl  4(%ecx),%ebx           ## eax/ebx=jnr 
        addl $8,nb100nf_innerjjnr(%esp)             ## advance pointer (unrolled 2) 
        prefetch 16(%ecx)            ## prefetch data - trial and error says 16 is best 

        movl nb100nf_charge(%ebp),%ecx     ## base of charge[] 
        movq nb100nf_iq(%esp),%mm5
        movd (%ecx,%eax,4),%mm3      ## charge[jnr1] 
        movd (%ecx,%ebx,4),%mm7          ## charge[jnr2] 
        punpckldq %mm7,%mm3          ## move charge 2 to high part of mm3 
        pfmul %mm5,%mm3              ## mm3 now has qq for both particles 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        movq  nb100nf_ix(%esp),%mm0
        movd  nb100nf_iz(%esp),%mm1
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
        punpckldq %mm6,%mm4             ## now 4 has rsq and 0 the seed for both pairs 
        movq %mm0,%mm2                  ## amd 3dnow N-R iteration to get full precision 
        pfmul %mm0,%mm0
        pfrsqit1 %mm4,%mm0
        pfrcpit2 %mm2,%mm0
        ## mm0=invsqrt
         ## do potential and fscal

        pfmul %mm0,%mm3         ## mm3 has both vcoul 
        ## update vctot 
        pfadd nb100nf_vctot(%esp),%mm3        ## add the earlier value  
        movq %mm3,nb100nf_vctot(%esp)         ## store the sum 

        ## should we do one more iteration? 
        subl $2,nb100nf_innerk(%esp)
        jl    _nb_kernel100nf_ia32_3dnow.nb100nf_finish_inner
        jmp   _nb_kernel100nf_ia32_3dnow.nb100nf_unroll_loop
_nb_kernel100nf_ia32_3dnow.nb100nf_finish_inner: 
        andl $1,nb100nf_innerk(%esp)
        jnz  _nb_kernel100nf_ia32_3dnow.nb100nf_single_inner
        jmp  _nb_kernel100nf_ia32_3dnow.nb100nf_updateouterdata
_nb_kernel100nf_ia32_3dnow.nb100nf_single_inner: 
        ## a single j particle iteration here - compare with the unrolled code for comments 
        movl  nb100nf_innerjjnr(%esp),%eax
        movl  (%eax),%eax       ## eax=jnr offset 

        movl nb100nf_charge(%ebp),%ecx
        pxor %mm6,%mm6
        movd nb100nf_iq(%esp),%mm6
        movd (%ecx,%eax,4),%mm7
        pfmul %mm7,%mm6         ## mm6=qq 

        leal  (%eax,%eax,2),%eax

        movq  nb100nf_ix(%esp),%mm0
        movd  nb100nf_iz(%esp),%mm1
        movq  (%esi,%eax,4),%mm2
        movd  8(%esi,%eax,4),%mm3
        pfsub %mm2,%mm0
        pfsub %mm3,%mm1
        pfmul %mm0,%mm0
        pfmul %mm1,%mm1
        pfacc %mm1,%mm0
        pfacc %mm1,%mm0         ## mm0=rsq 

        pfrsqrt %mm0,%mm1
        movq %mm1,%mm2
        pfmul %mm1,%mm1
        pfrsqit1 %mm0,%mm1
        pfrcpit2 %mm2,%mm1      ## mm1=invsqrt 
        ## calculate potential 
        pfmul %mm1,%mm6         ## mm6=vcoul 
        ## update vctot 
        pfadd nb100nf_vctot(%esp),%mm6
        movq %mm6,nb100nf_vctot(%esp)
        ## done! 
_nb_kernel100nf_ia32_3dnow.nb100nf_updateouterdata: 
        ## get n from stack
        movl nb100nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb100nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        movq  nb100nf_vctot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb100nf_Vc(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment vc[gid] 
        ## finish if last 
        movl nb100nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel100nf_ia32_3dnow.nb100nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb100nf_n(%esp)
        jmp _nb_kernel100nf_ia32_3dnow.nb100nf_outer
_nb_kernel100nf_ia32_3dnow.nb100nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb100nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel100nf_ia32_3dnow.nb100nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel100nf_ia32_3dnow.nb100nf_threadloop
_nb_kernel100nf_ia32_3dnow.nb100nf_end: 
        femms

        movl nb100nf_nouter(%esp),%eax
        movl nb100nf_ninner(%esp),%ebx
        movl nb100nf_outeriter(%ebp),%ecx
        movl nb100nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        addl $68,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret






