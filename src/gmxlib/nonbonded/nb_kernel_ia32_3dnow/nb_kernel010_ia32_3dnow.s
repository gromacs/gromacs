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




.globl nb_kernel010_ia32_3dnow
.globl _nb_kernel010_ia32_3dnow
nb_kernel010_ia32_3dnow:        
_nb_kernel010_ia32_3dnow:       
.set nb010_p_nri, 8
.set nb010_iinr, 12
.set nb010_jindex, 16
.set nb010_jjnr, 20
.set nb010_shift, 24
.set nb010_shiftvec, 28
.set nb010_fshift, 32
.set nb010_gid, 36
.set nb010_pos, 40
.set nb010_faction, 44
.set nb010_charge, 48
.set nb010_p_facel, 52
.set nb010_p_krf, 56
.set nb010_p_crf, 60
.set nb010_Vc, 64
.set nb010_type, 68
.set nb010_p_ntype, 72
.set nb010_vdwparam, 76
.set nb010_Vvdw, 80
.set nb010_p_tabscale, 84
.set nb010_VFtab, 88
.set nb010_invsqrta, 92
.set nb010_dvda, 96
.set nb010_p_gbtabscale, 100
.set nb010_GBtab, 104
.set nb010_p_nthreads, 108
.set nb010_count, 112
.set nb010_mtx, 116
.set nb010_outeriter, 120
.set nb010_inneriter, 124
.set nb010_work, 128
        ## stack offsets for local variables 
.set nb010_is3, 0
.set nb010_ii3, 4
.set nb010_ix, 8
.set nb010_iy, 12
.set nb010_iz, 16
.set nb010_Vvdwtot, 20
.set nb010_c6, 28
.set nb010_c12, 36
.set nb010_six, 44
.set nb010_twelve, 52
.set nb010_ntia, 60
.set nb010_innerjjnr, 64
.set nb010_innerk, 68
.set nb010_fix, 72
.set nb010_fiy, 76
.set nb010_fiz, 80
.set nb010_dx1, 84
.set nb010_dy1, 88
.set nb010_dz1, 92
.set nb010_dx2, 96
.set nb010_dy2, 100
.set nb010_dz2, 104
.set nb010_n, 108                           ## idx for outer loop
.set nb010_nn1, 112                         ## number of outer iterations
.set nb010_nri, 116
.set nb010_ntype, 120
.set nb010_nouter, 124
.set nb010_ninner, 128

        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $132,%esp                                          ## local stack space 
        femms

        movl $0x40c00000,%eax ## fp 6.0
        movl $0x41400000,%ebx ## fp 12.0

        movl %eax,nb010_six(%esp)
        movl %eax,nb010_six+4(%esp)
        movl %ebx,nb010_twelve(%esp)
        movl %ebx,nb010_twelve+4(%esp)

        movl nb010_p_nri(%ebp),%ecx
        movl nb010_p_ntype(%ebp),%edx
        movl (%ecx),%ecx
        movl (%edx),%edx
        movl %ecx,nb010_nri(%esp)
        movl %edx,nb010_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb010_nouter(%esp)
        movl %eax,nb010_ninner(%esp)

_nb_kernel010_ia32_3dnow.nb010_threadloop: 
        movl  nb010_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel010_ia32_3dnow.nb010_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $10,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel010_ia32_3dnow.nb010_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb010_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb010_n(%esp)
        movl %ebx,nb010_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists in ecx
        movl %eax,%esi                          ## copy n to esi

        jg  _nb_kernel010_ia32_3dnow.nb010_outerstart
        jmp _nb_kernel010_ia32_3dnow.nb010_end

_nb_kernel010_ia32_3dnow.nb010_outerstart: 
        ## ebx contains number of outer iterations
        addl nb010_nouter(%esp),%ebx
        movl %ebx,nb010_nouter(%esp)

_nb_kernel010_ia32_3dnow.nb010_outer: 

        movl  nb010_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx                        ## ebx=3*is 
        movl  %ebx,nb010_is3(%esp)              ## store is3 

        movl  nb010_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movq  (%eax,%ebx,4),%mm0                        ## move shX/shY to mm0 and shZ to mm1. 
        movd  8(%eax,%ebx,4),%mm1

        movl  nb010_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx                                        ## ebx =ii 

        movl  nb010_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb010_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb010_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx                        ## ebx = 3*ii=ii3 
        movl  nb010_pos(%ebp),%eax              ## eax = base of pos[] 

        pfadd (%eax,%ebx,4),%mm0                        ## ix = shX + posX (and iy too) 
        movd  8(%eax,%ebx,4),%mm3               ## cant use direct memory add for 4 bytes (iz) 
        movl  %ebx,nb010_ii3(%esp)
        pfadd %mm3,%mm1
        movq  %mm0,nb010_ix(%esp)
        movd  %mm1,nb010_iz(%esp)

        ## clear total potential and i forces 
        pxor  %mm7,%mm7
        movq  %mm7,nb010_Vvdwtot(%esp)
        movq  %mm7,nb010_fix(%esp)
        movd  %mm7,nb010_fiz(%esp)

        movl  nb010_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx                                        ## jindex[n] 
        movl  4(%eax,%esi,4),%edx                               ## jindex[n+1] 
        subl  %ecx,%edx                                         ## number of innerloop atoms 

        movl  nb010_pos(%ebp),%esi
        movl  nb010_faction(%ebp),%edi
        movl  nb010_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb010_innerjjnr(%esp)       ##  pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb010_ninner(%esp),%ecx
        movl  %ecx,nb010_ninner(%esp)
        movl  %edx,nb010_innerk(%esp)      ## number of innerloop atoms 
        addl  $0,%edx
        jge   _nb_kernel010_ia32_3dnow.nb010_unroll_loop
        jmp   _nb_kernel010_ia32_3dnow.nb010_finish_inner
_nb_kernel010_ia32_3dnow.nb010_unroll_loop: 
        ## paired innerloop starts here 
        movl  nb010_innerjjnr(%esp),%ecx        ## pointer to jjnr[k] 
        movl  (%ecx),%eax
        movl  4(%ecx),%ebx                              ## eax/ebx=jnr 
        addl  $8,nb010_innerjjnr(%esp)                  ## advance pointer (unrolled 2) 
        prefetch 16(%ecx)                               ## prefetch data - trial and error says 16 is best      

        movl nb010_type(%ebp),%ecx
        movl (%ecx,%eax,4),%edx                 ## type [jnr1] 
        movl (%ecx,%ebx,4),%ecx                         ## type [jnr2] 

        movl nb010_vdwparam(%ebp),%esi                  ## base of vdwparam  
        shll %edx
        shll %ecx
        addl nb010_ntia(%esp),%edx                      ## tja = ntia + 2*type 
        addl nb010_ntia(%esp),%ecx

        movq (%esi,%edx,4),%mm5                         ## mm5 = 1st c6 / c12           
        movq (%esi,%ecx,4),%mm7                         ## mm7 = 2nd c6 / c12   
        movq %mm5,%mm6
        punpckldq %mm7,%mm5                                     ## mm5 = 1st c6 / 2nd c6 
        punpckhdq %mm7,%mm6                                     ## mm6 = 1st c12 / 2nd c12 
        movq %mm5,nb010_c6(%esp)
        movq %mm6,nb010_c12(%esp)

        leal  (%eax,%eax,2),%eax                        ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        movl  nb010_pos(%ebp),%esi

        movq  nb010_ix(%esp),%mm0
        movd  nb010_iz(%esp),%mm1
        movq  (%esi,%eax,4),%mm4                        ## fetch first j coordinates 
        movd  8(%esi,%eax,4),%mm5
        pfsubr %mm0,%mm4                                        ## dr = ir - jr  
        pfsubr %mm1,%mm5
        movq  %mm4,nb010_dx1(%esp)          ## store dr 
        movd  %mm5,nb010_dz1(%esp)
        pfmul %mm4,%mm4                                         ## square dx,dy,dz                       
        pfmul %mm5,%mm5
        pfacc %mm5,%mm4                                         ## accumulate to get dx*dx+ dy*dy+ dz*dz 
        pfacc %mm5,%mm4                                         ## first rsq in lower mm4 

        movq  (%esi,%ebx,4),%mm6                        ## fetch second j coordinates  
        movd  8(%esi,%ebx,4),%mm7

        pfsubr %mm0,%mm6                                        ## dr = ir - jr  
        pfsubr %mm1,%mm7
        movq  %mm6,nb010_dx2(%esp)          ## store dr 
        movd  %mm7,nb010_dz2(%esp)
        pfmul %mm6,%mm6                                         ## square dx,dy,dz 
        pfmul %mm7,%mm7
        pfacc %mm7,%mm6                                         ## accumulate to get dx*dx+ dy*dy+ dz*dz 
        pfacc %mm7,%mm6                                         ## second rsq in lower mm6 

        pfrcp %mm4,%mm0                                         ## lookup reciprocal seed  
        pfrcp %mm6,%mm1

        punpckldq %mm1,%mm0
        punpckldq %mm6,%mm4                             ## now 4 has rsq and 0 the seed for both pairs. 
                                                        ## amd 3dnow N-R iteration to get full precision. 
        pfrcpit1 %mm0,%mm4
        pfrcpit2 %mm0,%mm4
        ## mm4 now contains invsq,
        ## do potential and fscal

        movq  %mm4,%mm0
        pfmul %mm0,%mm4
        pfmul %mm0,%mm4                                 ## mm4=rinvsix 
        movq  %mm4,%mm5
        pfmul %mm5,%mm5                                         ## mm5=rinvtwelve 

        pfmul nb010_c12(%esp),%mm5
        pfmul nb010_c6(%esp),%mm4
        movq %mm5,%mm6                                          ## mm6 is Vvdw12-Vvdw6  
        pfsub %mm4,%mm6

        pfmul nb010_six(%esp),%mm4

        pfmul nb010_twelve(%esp),%mm5
        pfsub %mm4,%mm5
        pfmul %mm5,%mm0                                         ## mm0 is total fscal now       

        prefetchw nb010_dx1(%esp)                       ## prefetch i forces to cache 

        ## spread fscalar to both positions 
        movq %mm0,%mm1
        punpckldq %mm0,%mm0
        punpckhdq %mm1,%mm1

        ## calc vector force 
        prefetchw (%edi,%eax,4)                         ## prefetch the 1st faction to cache 
        movq nb010_dx1(%esp),%mm2               ## fetch dr 
        movd nb010_dz1(%esp),%mm3

        ## update Vvdwtot  
        pfadd nb010_Vvdwtot(%esp),%mm6       ## add the earlier value 
        movq %mm6,nb010_Vvdwtot(%esp)        ## store the sum 

        prefetchw (%edi,%ebx,4)                         ## prefetch the 2nd faction to cache 
        pfmul %mm0,%mm2                                         ## mult by fs  
        pfmul %mm0,%mm3

        movq nb010_dx2(%esp),%mm4               ## fetch dr 
        movd nb010_dz2(%esp),%mm5
        pfmul %mm1,%mm4                                         ## mult by fs  
        pfmul %mm1,%mm5
        ## update i forces 

        movq nb010_fix(%esp),%mm0
        movd nb010_fiz(%esp),%mm1
        pfadd %mm2,%mm0
        pfadd %mm3,%mm1

        pfadd %mm4,%mm0
        pfadd %mm5,%mm1
        movq %mm0,nb010_fix(%esp)
        movd %mm1,nb010_fiz(%esp)
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
        subl $2,nb010_innerk(%esp)
        jl    _nb_kernel010_ia32_3dnow.nb010_finish_inner
        jmp   _nb_kernel010_ia32_3dnow.nb010_unroll_loop
_nb_kernel010_ia32_3dnow.nb010_finish_inner: 
        andl $1,nb010_innerk(%esp)
        jnz  _nb_kernel010_ia32_3dnow.nb010_single_inner
        jmp  _nb_kernel010_ia32_3dnow.nb010_updateouterdata
_nb_kernel010_ia32_3dnow.nb010_single_inner: 
        ## a single j particle iteration here - compare with the unrolled code for comments 
        movl  nb010_innerjjnr(%esp),%eax
        movl  (%eax),%eax                                       ## eax=jnr offset 

        movl nb010_vdwparam(%ebp),%esi
        movl nb010_type(%ebp),%ecx
        movl (%ecx,%eax,4),%edx                 ## type [jnr1] 
        shll %edx
        addl nb010_ntia(%esp),%edx                      ## tja = ntia + 2*type 
        movd (%esi,%edx,4),%mm5                         ## mm5 = 1st c6                 
        movq %mm5,nb010_c6(%esp)
        movd 4(%esi,%edx,4),%mm5                        ## mm5 = 1st c12                
        movq %mm5,nb010_c12(%esp)

        movl  nb010_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        movq  nb010_ix(%esp),%mm0
        movd  nb010_iz(%esp),%mm1
        movq  (%esi,%eax,4),%mm4
        movd  8(%esi,%eax,4),%mm5
        pfsubr %mm0,%mm4
        pfsubr %mm1,%mm5
        movq  %mm4,nb010_dx1(%esp)
        pfmul %mm4,%mm4
        movd  %mm5,nb010_dz1(%esp)
        pfmul %mm5,%mm5
        pfacc %mm5,%mm4
        pfacc %mm5,%mm4                                         ## mm4=rsq 

        pfrcp %mm4,%mm0
        pfrcpit1 %mm0,%mm4
        pfrcpit2 %mm0,%mm4                                      ## mm4=invsq 
        ## calculate potentials and scalar force 
        movq  %mm4,%mm0

        pfmul %mm0,%mm4
        pfmul %mm0,%mm4                 ## mm4=rinvsix 
        movq  %mm4,%mm5
        pfmul %mm5,%mm5             ## mm5=rinvtwelve 

        pfmul nb010_c12(%esp),%mm5
        pfmul nb010_c6(%esp),%mm4
        movq %mm5,%mm6  ## mm6 is Vvdw12-Vvdw6 
        pfsub %mm4,%mm6

        pfmul nb010_six(%esp),%mm4

        pfmul nb010_twelve(%esp),%mm5
        pfsub %mm4,%mm5
        pfmul %mm5,%mm0   ## mm0 is total fscal now 

        ## update Vvdwtot 
        pfadd nb010_Vvdwtot(%esp),%mm6        ## add the earlier value 
        movq %mm6,nb010_Vvdwtot(%esp)         ## store the sum   

        ## spread fscalar to both positions 
        punpckldq %mm0,%mm0
        ## calc vectorial force 
        prefetchw (%edi,%eax,4) ## prefetch faction to cache 
        movq nb010_dx1(%esp),%mm2
        movd nb010_dz1(%esp),%mm3

        pfmul %mm0,%mm2
        pfmul %mm0,%mm3

        ## update i particle force 
        movq nb010_fix(%esp),%mm0
        movd nb010_fiz(%esp),%mm1
        pfadd %mm2,%mm0
        pfadd %mm3,%mm1
        movq %mm0,nb010_fix(%esp)
        movd %mm1,nb010_fiz(%esp)
        ## update j particle force 
        movq (%edi,%eax,4),%mm0
        movd 8(%edi,%eax,4),%mm1
        pfsub %mm2,%mm0
        pfsub %mm3,%mm1
        movq %mm0,(%edi,%eax,4)
        movd %mm1,8(%edi,%eax,4)
        ## done! 
_nb_kernel010_ia32_3dnow.nb010_updateouterdata: 
        movl  nb010_ii3(%esp),%ecx

        movq  (%edi,%ecx,4),%mm6       ## increment i force 
        movd  8(%edi,%ecx,4),%mm7
        pfadd nb010_fix(%esp),%mm6
        pfadd nb010_fiz(%esp),%mm7
        movq  %mm6,(%edi,%ecx,4)
        movd  %mm7,8(%edi,%ecx,4)

        movl  nb010_fshift(%ebp),%ebx      ## increment fshift force 
        movl  nb010_is3(%esp),%edx

        movq  (%ebx,%edx,4),%mm6
        movd  8(%ebx,%edx,4),%mm7
        pfadd nb010_fix(%esp),%mm6
        pfadd nb010_fiz(%esp),%mm7
        movq  %mm6,(%ebx,%edx,4)
        movd  %mm7,8(%ebx,%edx,4)

        ## get n from stack
        movl nb010_n(%esp),%esi
        ## get group index for i particle 
        movl  nb010_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        movq  nb010_Vvdwtot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb010_Vvdw(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment Vvdw[gid] 

       ## finish if last 
        movl nb010_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel010_ia32_3dnow.nb010_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb010_n(%esp)
        jmp _nb_kernel010_ia32_3dnow.nb010_outer
_nb_kernel010_ia32_3dnow.nb010_outerend: 
        ## check if more outer neighborlists remain
        movl  nb010_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel010_ia32_3dnow.nb010_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel010_ia32_3dnow.nb010_threadloop
_nb_kernel010_ia32_3dnow.nb010_end: 
        femms

        movl nb010_nouter(%esp),%eax
        movl nb010_ninner(%esp),%ebx
        movl nb010_outeriter(%ebp),%ecx
        movl nb010_inneriter(%ebp),%edx
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





.globl nb_kernel010nf_ia32_3dnow
.globl _nb_kernel010nf_ia32_3dnow
nb_kernel010nf_ia32_3dnow:      
_nb_kernel010nf_ia32_3dnow:     
.set nb010nf_p_nri, 8
.set nb010nf_iinr, 12
.set nb010nf_jindex, 16
.set nb010nf_jjnr, 20
.set nb010nf_shift, 24
.set nb010nf_shiftvec, 28
.set nb010nf_fshift, 32
.set nb010nf_gid, 36
.set nb010nf_pos, 40
.set nb010nf_faction, 44
.set nb010nf_charge, 48
.set nb010nf_p_facel, 52
.set nb010nf_p_krf, 56
.set nb010nf_p_crf, 60
.set nb010nf_Vc, 64
.set nb010nf_type, 68
.set nb010nf_p_ntype, 72
.set nb010nf_vdwparam, 76
.set nb010nf_Vvdw, 80
.set nb010nf_p_tabscale, 84
.set nb010nf_VFtab, 88
.set nb010nf_invsqrta, 92
.set nb010nf_dvda, 96
.set nb010nf_p_gbtabscale, 100
.set nb010nf_GBtab, 104
.set nb010nf_p_nthreads, 108
.set nb010nf_count, 112
.set nb010nf_mtx, 116
.set nb010nf_outeriter, 120
.set nb010nf_inneriter, 124
.set nb010nf_work, 128
        ## stack offsets for local variables 
.set nb010nf_is3, 0
.set nb010nf_ii3, 4
.set nb010nf_ix, 8
.set nb010nf_iy, 12
.set nb010nf_iz, 16
.set nb010nf_Vvdwtot, 20
.set nb010nf_c6, 28
.set nb010nf_c12, 36
.set nb010nf_ntia, 44
.set nb010nf_innerjjnr, 48
.set nb010nf_innerk, 52
.set nb010nf_n, 56                         ## idx for outer loop
.set nb010nf_nn1, 60                       ## number of outer iterations
.set nb010nf_nri, 64
.set nb010nf_ntype, 68
.set nb010nf_nouter, 72
.set nb010nf_ninner, 76
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $80,%esp           ## local stack space 
        femms

        movl nb010nf_p_nri(%ebp),%ecx
        movl nb010nf_p_ntype(%ebp),%edx
        movl (%ecx),%ecx
        movl (%edx),%edx
        movl %ecx,nb010nf_nri(%esp)
        movl %edx,nb010nf_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb010nf_nouter(%esp)
        movl %eax,nb010nf_ninner(%esp)

_nb_kernel010nf_ia32_3dnow.nb010nf_threadloop: 
        movl  nb010nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel010nf_ia32_3dnow.nb010nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $10,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel010nf_ia32_3dnow.nb010nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb010nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb010nf_n(%esp)
        movl %ebx,nb010nf_nn1(%esp)

        subl %eax,%ebx          # # calc number of outer lists in ecx
        movl %eax,%esi  # # copy n to esi

        jg  _nb_kernel010nf_ia32_3dnow.nb010nf_outerstart
        jmp _nb_kernel010nf_ia32_3dnow.nb010nf_end

_nb_kernel010nf_ia32_3dnow.nb010nf_outerstart: 
        ## # ebx contains number of outer iterations
        addl nb010nf_nouter(%esp),%ebx
        movl %ebx,nb010nf_nouter(%esp)

_nb_kernel010nf_ia32_3dnow.nb010nf_outer: 
        movl  nb010nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb010nf_is3(%esp)            ## store is3 

        movl  nb010nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movq  (%eax,%ebx,4),%mm0        ## move shX/shY to mm0 and shZ to mm1. 
        movd  8(%eax,%ebx,4),%mm1

        movl  nb010nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movl  nb010nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb010nf_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb010nf_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb010nf_pos(%ebp),%eax      ## eax = base of pos[] 

        pfadd (%eax,%ebx,4),%mm0    ## ix = shX + posX (and iy too) 
        movd  8(%eax,%ebx,4),%mm3       ## cant use direct memory add for 4 bytes (iz) 
        movl  %ebx,nb010nf_ii3(%esp)
        pfadd %mm3,%mm1
        movq  %mm0,nb010nf_ix(%esp)
        movd  %mm1,nb010nf_iz(%esp)

        ## clear total potential 
        pxor  %mm7,%mm7
        movq  %mm7,nb010nf_Vvdwtot(%esp)

        movl  nb010nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb010nf_pos(%ebp),%esi
        movl  nb010nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb010nf_innerjjnr(%esp)       ##  pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb010nf_ninner(%esp),%ecx
        movl  %ecx,nb010nf_ninner(%esp)
        movl  %edx,nb010nf_innerk(%esp)      ## number of innerloop atoms 
        addl  $0,%edx
        jge   _nb_kernel010nf_ia32_3dnow.nb010nf_unroll_loop
        jmp   _nb_kernel010nf_ia32_3dnow.nb010nf_finish_inner
_nb_kernel010nf_ia32_3dnow.nb010nf_unroll_loop: 
        ## paired innerloop starts here 
        movl  nb010nf_innerjjnr(%esp),%ecx       ## pointer to jjnr[k] 
        movl  (%ecx),%eax
        movl  4(%ecx),%ebx           ## eax/ebx=jnr 
        addl $8,nb010nf_innerjjnr(%esp)             ## advance pointer (unrolled 2) 
        prefetch 16(%ecx)            ## prefetch data - trial and error says 16 is best         

        movl nb010nf_type(%ebp),%ecx
        movl (%ecx,%eax,4),%edx          ## type [jnr1] 
        movl (%ecx,%ebx,4),%ecx      ## type [jnr2] 

        movl nb010nf_vdwparam(%ebp),%esi                ## base of vdwparam  
        shll %edx
        shll %ecx
        addl nb010nf_ntia(%esp),%edx         ## tja = ntia + 2*type 
        addl nb010nf_ntia(%esp),%ecx

        movq (%esi,%edx,4),%mm5         ## mm5 = 1st c6 / c12           
        movq (%esi,%ecx,4),%mm7         ## mm7 = 2nd c6 / c12   
        movq %mm5,%mm6
        punpckldq %mm7,%mm5             ## mm5 = 1st c6 / 2nd c6 
        punpckhdq %mm7,%mm6             ## mm6 = 1st c12 / 2nd c12 
        movq %mm5,nb010nf_c6(%esp)
        movq %mm6,nb010nf_c12(%esp)

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        movl  nb010nf_pos(%ebp),%esi

        movq  nb010nf_ix(%esp),%mm0
        movd  nb010nf_iz(%esp),%mm1
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

        pfrcp %mm4,%mm0              ## lookup reciprocal seed  
        pfrcp %mm6,%mm1

        punpckldq %mm1,%mm0
        punpckldq %mm6,%mm4             ## now 4 has rsq and 0 the seed for both pairs. 
                                        ## amd 3dnow N-R iteration to get full precision. 
        pfrcpit1 %mm0,%mm4
        pfrcpit2 %mm0,%mm4
        ## mm4 now contains invsq,
         ## do potential and fscal

        movq  %mm4,%mm0
        pfmul %mm0,%mm4
        pfmul %mm0,%mm4                 ## mm4=rinvsix 
        movq  %mm4,%mm5
        pfmul %mm5,%mm5             ## mm5=rinvtwelve 

        pfmul nb010nf_c12(%esp),%mm5
        pfmul nb010nf_c6(%esp),%mm4
        movq %mm5,%mm6  ## mm6 is Vvdw12-Vvdw6  
        pfsub %mm4,%mm6
        ## update Vvdwtot  
        pfadd nb010nf_Vvdwtot(%esp),%mm6        ## add the earlier value 
        movq %mm6,nb010nf_Vvdwtot(%esp)         ## store the sum 

        ## should we do one more iteration? 
        subl $2,nb010nf_innerk(%esp)
        jl    _nb_kernel010nf_ia32_3dnow.nb010nf_finish_inner
        jmp   _nb_kernel010nf_ia32_3dnow.nb010nf_unroll_loop
_nb_kernel010nf_ia32_3dnow.nb010nf_finish_inner: 
        andl $1,nb010nf_innerk(%esp)
        jnz  _nb_kernel010nf_ia32_3dnow.nb010nf_single_inner
        jmp  _nb_kernel010nf_ia32_3dnow.nb010nf_updateouterdata
_nb_kernel010nf_ia32_3dnow.nb010nf_single_inner: 
        ## a single j particle iteration here - compare with the unrolled code for comments 
        movl  nb010nf_innerjjnr(%esp),%eax
        movl  (%eax),%eax       ## eax=jnr offset 

        movl nb010nf_vdwparam(%ebp),%esi
        movl nb010nf_type(%ebp),%ecx
        movl (%ecx,%eax,4),%edx         ## type [jnr1] 
        shll %edx
        addl nb010nf_ntia(%esp),%edx        ## tja = ntia + 2*type 
        movd (%esi,%edx,4),%mm5         ## mm5 = 1st c6                 
        movq %mm5,nb010nf_c6(%esp)
        movd 4(%esi,%edx,4),%mm5        ## mm5 = 1st c12                
        movq %mm5,nb010nf_c12(%esp)

        movl  nb010nf_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        movq  nb010nf_ix(%esp),%mm0
        movd  nb010nf_iz(%esp),%mm1
        movq  (%esi,%eax,4),%mm4
        movd  8(%esi,%eax,4),%mm5
        pfsubr %mm0,%mm4
        pfsubr %mm1,%mm5
        pfmul %mm4,%mm4
        pfmul %mm5,%mm5
        pfacc %mm5,%mm4
        pfacc %mm5,%mm4         ## mm4=rsq 

        pfrcp %mm4,%mm0
        pfrcpit1 %mm0,%mm4
        pfrcpit2 %mm0,%mm4      ## mm4=invsq 
        ## calculate potentials and scalar force 
        movq  %mm4,%mm0

        pfmul %mm0,%mm4
        pfmul %mm0,%mm4                 ## mm4=rinvsix 
        movq  %mm4,%mm5
        pfmul %mm5,%mm5             ## mm5=rinvtwelve 

        pfmul nb010nf_c12(%esp),%mm5
        pfmul nb010nf_c6(%esp),%mm4
        movq %mm5,%mm6  ## mm6 is Vvdw12-Vvdw6 
        pfsub %mm4,%mm6
        ## update Vvdwtot 
        pfadd nb010nf_Vvdwtot(%esp),%mm6        ## add the earlier value 
        movq %mm6,nb010nf_Vvdwtot(%esp)         ## store the sum   

_nb_kernel010nf_ia32_3dnow.nb010nf_updateouterdata: 
        ## get n from stack
        movl nb010nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb010nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        movq  nb010nf_Vvdwtot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb010nf_Vvdw(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment Vvdw[gid] 

       ## finish if last 
        movl nb010nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel010nf_ia32_3dnow.nb010nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb010nf_n(%esp)
        jmp _nb_kernel010nf_ia32_3dnow.nb010nf_outer
_nb_kernel010nf_ia32_3dnow.nb010nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb010nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel010nf_ia32_3dnow.nb010nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel010nf_ia32_3dnow.nb010nf_threadloop
_nb_kernel010nf_ia32_3dnow.nb010nf_end: 
        femms

        movl nb010nf_nouter(%esp),%eax
        movl nb010nf_ninner(%esp),%ebx
        movl nb010nf_outeriter(%ebp),%ecx
        movl nb010nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        addl $80,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret





