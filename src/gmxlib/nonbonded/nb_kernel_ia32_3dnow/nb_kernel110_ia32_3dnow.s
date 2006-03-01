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




.globl nb_kernel110_ia32_3dnow
.globl _nb_kernel110_ia32_3dnow
nb_kernel110_ia32_3dnow:        
_nb_kernel110_ia32_3dnow:       
.set nb110_p_nri, 8
.set nb110_iinr, 12
.set nb110_jindex, 16
.set nb110_jjnr, 20
.set nb110_shift, 24
.set nb110_shiftvec, 28
.set nb110_fshift, 32
.set nb110_gid, 36
.set nb110_pos, 40
.set nb110_faction, 44
.set nb110_charge, 48
.set nb110_p_facel, 52
.set nb110_p_krf, 56
.set nb110_p_crf, 60
.set nb110_Vc, 64
.set nb110_type, 68
.set nb110_p_ntype, 72
.set nb110_vdwparam, 76
.set nb110_Vvdw, 80
.set nb110_p_tabscale, 84
.set nb110_VFtab, 88
.set nb110_invsqrta, 92
.set nb110_dvda, 96
.set nb110_p_gbtabscale, 100
.set nb110_GBtab, 104
.set nb110_p_nthreads, 108
.set nb110_count, 112
.set nb110_mtx, 116
.set nb110_outeriter, 120
.set nb110_inneriter, 124
.set nb110_work, 128
        ## stack offsets for local variables 
.set nb110_is3, 0
.set nb110_ii3, 4
.set nb110_ix, 8
.set nb110_iy, 12
.set nb110_iz, 16
.set nb110_iq, 20
.set nb110_vctot, 28
.set nb110_Vvdwtot, 36
.set nb110_c6, 44
.set nb110_c12, 52
.set nb110_six, 60
.set nb110_twelve, 68
.set nb110_ntia, 76
.set nb110_innerjjnr, 80
.set nb110_innerk, 84
.set nb110_fix, 88
.set nb110_fiy, 92
.set nb110_fiz, 96
.set nb110_dx1, 100
.set nb110_dy1, 104
.set nb110_dz1, 108
.set nb110_dx2, 112
.set nb110_dy2, 116
.set nb110_dz2, 120
.set nb110_n, 124                           ## idx for outer loop
.set nb110_nn1, 128                         ## number of outer iterations
.set nb110_nri, 132
.set nb110_facel, 136
.set nb110_ntype, 140
.set nb110_nouter, 144
.set nb110_ninner, 148

        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $152,%esp          ## local stack space 
        femms
        ## move data to local stack  
        movl nb110_p_nri(%ebp),%ecx
        movl nb110_p_ntype(%ebp),%edx
        movl nb110_p_facel(%ebp),%esi
        movl (%ecx),%ecx
        movl (%edx),%edx
        movl (%esi),%esi
        movl %ecx,nb110_nri(%esp)
        movl %edx,nb110_ntype(%esp)
        movl %esi,nb110_facel(%esp)

        movl $0x40c00000,%eax ## fp 6.0
        movl $0x41400000,%ebx ## fp 12.0

        movl %eax,nb110_six(%esp)
        movl %eax,nb110_six+4(%esp)
        movl %ebx,nb110_twelve(%esp)
        movl %ebx,nb110_twelve+4(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb110_nouter(%esp)
        movl %eax,nb110_ninner(%esp)

_nb_kernel110_ia32_3dnow.nb110_threadloop: 
        movl  nb110_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel110_ia32_3dnow.nb110_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel110_ia32_3dnow.nb110_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb110_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb110_n(%esp)
        movl %ebx,nb110_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi

        jg  _nb_kernel110_ia32_3dnow.nb110_outerstart
        jmp _nb_kernel110_ia32_3dnow.nb110_end

_nb_kernel110_ia32_3dnow.nb110_outerstart: 
        ## ebx contains number of outer iterations
        addl nb110_nouter(%esp),%ebx
        movl %ebx,nb110_nouter(%esp)

_nb_kernel110_ia32_3dnow.nb110_outer: 
        movl  nb110_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb110_is3(%esp)      ## store is3 

        movl  nb110_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movq  (%eax,%ebx,4),%mm0        ## move shX/shY to mm0 and shZ to mm1 
        movd  8(%eax,%ebx,4),%mm1

        movl  nb110_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx            ## ebx=ii 

        movl  nb110_charge(%ebp),%edx
        movd  (%edx,%ebx,4),%mm2    ## mm2=charge[ii] 
        pfmul nb110_facel(%esp),%mm2
        punpckldq %mm2,%mm2         ## spread to both halves 
        movq  %mm2,nb110_iq(%esp)           ## iq =facel*charge[ii] 

        movl  nb110_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb110_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb110_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb110_pos(%ebp),%eax      ## eax = base of pos[] 

        pfadd (%eax,%ebx,4),%mm0    ## ix = shX + posX (and iy too) 
        movd  8(%eax,%ebx,4),%mm3       ## cant use direct memory add for 4 bytes (iz) 
        movl  %ebx,nb110_ii3(%esp)
        pfadd %mm3,%mm1
        movq  %mm0,nb110_ix(%esp)
        movd  %mm1,nb110_iz(%esp)

        ## clear total potential and i forces 
        pxor  %mm7,%mm7
        movq  %mm7,nb110_vctot(%esp)
        movq  %mm7,nb110_Vvdwtot(%esp)
        movq  %mm7,nb110_fix(%esp)
        movd  %mm7,nb110_fiz(%esp)

        movl  nb110_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb110_pos(%ebp),%esi
        movl  nb110_faction(%ebp),%edi
        movl  nb110_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb110_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb110_ninner(%esp),%ecx
        movl  %ecx,nb110_ninner(%esp)
        movl  %edx,nb110_innerk(%esp)      ## number of innerloop atoms 
        addl  $0,%edx
        jge   _nb_kernel110_ia32_3dnow.nb110_unroll_loop
        jmp   _nb_kernel110_ia32_3dnow.nb110_finish_inner
_nb_kernel110_ia32_3dnow.nb110_unroll_loop: 
        ## paired innerloop starts here 
        movl  nb110_innerjjnr(%esp),%ecx       ## pointer to jjnr[k] 
        movl  (%ecx),%eax
        movl  4(%ecx),%ebx           ## eax/ebx=jnr 
        addl $8,nb110_innerjjnr(%esp)             ## advance pointer (unrolled 2) 
        prefetch 16(%ecx)            ## prefetch data - trial and error says 16 is best 

        movl nb110_charge(%ebp),%ecx     ## base of charge[] 
        movq nb110_iq(%esp),%mm5
        movd (%ecx,%eax,4),%mm3      ## charge[jnr1] 
    punpckldq (%ecx,%ebx,4),%mm3     ## move charge 2 to high part of mm3 
        pfmul %mm5,%mm3              ## mm3 now has qq for both particles 

        movl nb110_type(%ebp),%ecx
        movl (%ecx,%eax,4),%edx          ## type [jnr1] 
        movl (%ecx,%ebx,4),%ecx      ## type [jnr2] 

        movl nb110_vdwparam(%ebp),%esi          ## base of vdwparam  
        shll %edx
        shll %ecx
        addl nb110_ntia(%esp),%edx           ## tja = ntia + 2*type 
        addl nb110_ntia(%esp),%ecx

        movq (%esi,%edx,4),%mm5         ## mm5 = 1st c6 / c12           
        movq (%esi,%ecx,4),%mm7         ## mm7 = 2nd c6 / c12   
        movq %mm5,%mm6
        punpckldq %mm7,%mm5             ## mm5 = 1st c6 / 2nd c6 
        punpckhdq %mm7,%mm6             ## mm6 = 1st c12 / 2nd c12 
        movq %mm5,nb110_c6(%esp)
        movq %mm6,nb110_c12(%esp)

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        movl  nb110_pos(%ebp),%esi

        movq  nb110_ix(%esp),%mm0
        movd  nb110_iz(%esp),%mm1
        movq  (%esi,%eax,4),%mm4     ## fetch first j coordinates 
        movd  8(%esi,%eax,4),%mm5
        pfsubr %mm0,%mm4             ## dr = ir - jr  
        pfsubr %mm1,%mm5
        movq  %mm4,nb110_dx1(%esp)           ## store dr 
        movd  %mm5,nb110_dz1(%esp)
        pfmul %mm4,%mm4              ## square dx,dy,dz                          
        pfmul %mm5,%mm5
        pfacc %mm5,%mm4              ## accumulate to get dx*dx+ dy*dy+ dz*dz 
        pfacc %mm5,%mm4              ## first rsq in lower mm4 

        movq  (%esi,%ebx,4),%mm6     ## fetch second j coordinates  
        movd  8(%esi,%ebx,4),%mm7

        pfsubr %mm0,%mm6             ## dr = ir - jr  
        pfsubr %mm1,%mm7
        movq  %mm6,nb110_dx2(%esp)           ## store dr 
        movd  %mm7,nb110_dz2(%esp)
        pfmul %mm6,%mm6              ## square dx,dy,dz 
        pfmul %mm7,%mm7
        pfacc %mm7,%mm6              ## accumulate to get dx*dx+ dy*dy+ dz*dz 
        pfacc %mm7,%mm6              ## second rsq in lower mm6 

    pfrsqrt %mm4,%mm0        ## lookup inverse square root seed 
    pfrsqrt %mm6,%mm1

        punpckldq %mm1,%mm0
        punpckldq %mm6,%mm4             ## now 4 has rsq and 0 the seed for both pairs 
    movq %mm0,%mm2                      ## amd 3dnow N-R iteration to get full precision 
        pfmul %mm0,%mm0
    pfrsqit1 %mm4,%mm0
    pfrcpit2 %mm2,%mm0
        movq %mm0,%mm1
        pfmul %mm0,%mm0
        ## mm0 now contains invsq, and mm1 invsqrt 
        ## do potential and fscal 
        movq %mm0,%mm4
        pfmul %mm0,%mm4
        pfmul %mm0,%mm4                 ## mm4=rinvsix 
        movq  %mm4,%mm5
        pfmul %mm5,%mm5             ## mm5=rinvtwelve 

        pfmul %mm1,%mm3         ## mm3 has vcoul for both interactions 
        movq  %mm3,%mm7     ## use mm7 for sum to make fscal  

        pfmul nb110_c12(%esp),%mm5
        pfmul nb110_c6(%esp),%mm4
        movq %mm5,%mm6  ## mm6 is Vvdw12-Vvdw6  
        pfsub %mm4,%mm6

        pfmul nb110_six(%esp),%mm4

        pfmul nb110_twelve(%esp),%mm5
        pfsub %mm4,%mm7
        pfadd %mm5,%mm7
        pfmul %mm7,%mm0   ## mm0 is total fscal now     

        prefetchw nb110_dx1(%esp)       ## prefetch i forces to cache 

        ## update vctot 
        pfadd nb110_vctot(%esp),%mm3        ## add the earlier value 
        movq %mm3,nb110_vctot(%esp)         ## store the sum       

        ## spread fscalar to both positions 
        movq %mm0,%mm1
        punpckldq %mm0,%mm0
        punpckhdq %mm1,%mm1

        ## calc vector force 
        prefetchw (%edi,%eax,4) ## prefetch the 1st faction to cache 
        movq nb110_dx1(%esp),%mm2       ## fetch dr 
        movd nb110_dz1(%esp),%mm3

        ## update Vvdwtot 
        pfadd nb110_Vvdwtot(%esp),%mm6        ## add the earlier value 
        movq %mm6,nb110_Vvdwtot(%esp)         ## store the sum       

        prefetchw (%edi,%ebx,4) ## prefetch the 2nd faction to cache 
        pfmul %mm0,%mm2         ## mult by fs  
        pfmul %mm0,%mm3

        movq nb110_dx2(%esp),%mm4       ## fetch dr 
        movd nb110_dz2(%esp),%mm5
        pfmul %mm1,%mm4         ## mult by fs  
        pfmul %mm1,%mm5
        ## update i forces 

        movq nb110_fix(%esp),%mm0
        movd nb110_fiz(%esp),%mm1
        pfadd %mm2,%mm0
        pfadd %mm3,%mm1

        pfadd %mm4,%mm0
        pfadd %mm5,%mm1
        movq %mm0,nb110_fix(%esp)
        movd %mm1,nb110_fiz(%esp)
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
        subl $2,nb110_innerk(%esp)
        jl    _nb_kernel110_ia32_3dnow.nb110_finish_inner
        jmp   _nb_kernel110_ia32_3dnow.nb110_unroll_loop
_nb_kernel110_ia32_3dnow.nb110_finish_inner: 
        andl $1,nb110_innerk(%esp)
        jnz  _nb_kernel110_ia32_3dnow.nb110_single_inner
        jmp  _nb_kernel110_ia32_3dnow.nb110_updateouterdata
_nb_kernel110_ia32_3dnow.nb110_single_inner: 
        ## a single j particle iteration here - compare with the unrolled code for comments 
        movl  nb110_innerjjnr(%esp),%eax
        movl  (%eax),%eax       ## eax=jnr offset 

        movl nb110_charge(%ebp),%ecx
        movd nb110_iq(%esp),%mm5
        movd (%ecx,%eax,4),%mm3
        pfmul %mm5,%mm3         ## mm3=qq 

        movl nb110_vdwparam(%ebp),%esi
        movl nb110_type(%ebp),%ecx
        movl (%ecx,%eax,4),%edx          ## type [jnr1] 
        shll %edx
        addl nb110_ntia(%esp),%edx           ## tja = ntia + 2*type 
        movd (%esi,%edx,4),%mm5         ## mm5 = 1st c6                 
        movq %mm5,nb110_c6(%esp)
        movd 4(%esi,%edx,4),%mm5        ## mm5 = 1st c12                
        movq %mm5,nb110_c12(%esp)


        movl  nb110_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        movq  nb110_ix(%esp),%mm0
        movd  nb110_iz(%esp),%mm1
        movq  (%esi,%eax,4),%mm4
        movd  8(%esi,%eax,4),%mm5
        pfsubr %mm0,%mm4
        pfsubr %mm1,%mm5
        movq  %mm4,nb110_dx1(%esp)
        pfmul %mm4,%mm4
        movd  %mm5,nb110_dz1(%esp)
        pfmul %mm5,%mm5
        pfacc %mm5,%mm4
        pfacc %mm5,%mm4         ## mm0=rsq 

    pfrsqrt %mm4,%mm0
    movq %mm0,%mm2
    pfmul %mm0,%mm0
    pfrsqit1 %mm4,%mm0
    pfrcpit2 %mm2,%mm0  ## mm1=invsqrt 
        movq  %mm0,%mm1
        pfmul %mm0,%mm0         ## mm0=invsq 
        ## calculate potentials and scalar force 
        movq %mm0,%mm4
        pfmul %mm0,%mm4
        pfmul %mm0,%mm4                 ## mm4=rinvsix 
        movq  %mm4,%mm5
        pfmul %mm5,%mm5             ## mm5=rinvtwelve 

        pfmul %mm1,%mm3         ## mm3 has vcoul for both interactions 
        movq  %mm3,%mm7     ## use mm7 for sum to make fscal  

        pfmul nb110_c12(%esp),%mm5
        pfmul nb110_c6(%esp),%mm4
        movq %mm5,%mm6  ## mm6 is Vvdw12-Vvdw6  
        pfsub %mm4,%mm6

        pfmul nb110_six(%esp),%mm4

        pfmul nb110_twelve(%esp),%mm5
        pfsub %mm4,%mm7
        pfadd %mm5,%mm7
        pfmul %mm7,%mm0   ## mm0 is total fscal now 

        ## update vctot 
        pfadd nb110_vctot(%esp),%mm3
        movq %mm3,nb110_vctot(%esp)

        ## update Vvdwtot 
        pfadd nb110_Vvdwtot(%esp),%mm6        ## add the earlier value 
        movq %mm6,nb110_Vvdwtot(%esp)         ## store the sum       

        ## spread fscalar to both positions 
        punpckldq %mm0,%mm0
        ## calc vectorial force 
        prefetchw (%edi,%eax,4) ## prefetch faction to cache  
        movq nb110_dx1(%esp),%mm2
        movd nb110_dz1(%esp),%mm3


        pfmul %mm0,%mm2
        pfmul %mm0,%mm3

        ## update i particle force 
        movq nb110_fix(%esp),%mm0
        movd nb110_fiz(%esp),%mm1
        pfadd %mm2,%mm0
        pfadd %mm3,%mm1
        movq %mm0,nb110_fix(%esp)
        movd %mm1,nb110_fiz(%esp)
        ## update j particle force 
        movq (%edi,%eax,4),%mm0
        movd 8(%edi,%eax,4),%mm1
        pfsub %mm2,%mm0
        pfsub %mm3,%mm1
        movq %mm0,(%edi,%eax,4)
        movd %mm1,8(%edi,%eax,4)
        ## done! 
_nb_kernel110_ia32_3dnow.nb110_updateouterdata: 
        movl  nb110_ii3(%esp),%ecx

        movq  (%edi,%ecx,4),%mm6       ## increment i force 
        movd  8(%edi,%ecx,4),%mm7
        pfadd nb110_fix(%esp),%mm6
        pfadd nb110_fiz(%esp),%mm7
        movq  %mm6,(%edi,%ecx,4)
        movd  %mm7,8(%edi,%ecx,4)

        movl  nb110_fshift(%ebp),%ebx      ## increment fshift force 
        movl  nb110_is3(%esp),%edx

        movq  (%ebx,%edx,4),%mm6
        movd  8(%ebx,%edx,4),%mm7
        pfadd nb110_fix(%esp),%mm6
        pfadd nb110_fiz(%esp),%mm7
        movq  %mm6,(%ebx,%edx,4)
        movd  %mm7,8(%ebx,%edx,4)

        ## get n from stack
        movl nb110_n(%esp),%esi
        ## get group index for i particle 
        movl  nb110_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        movq  nb110_vctot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb110_Vc(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment vc[gid] 

        movq  nb110_Vvdwtot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb110_Vvdw(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment Vvdw[gid] 

       ## finish if last 
        movl nb110_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel110_ia32_3dnow.nb110_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb110_n(%esp)
        jmp _nb_kernel110_ia32_3dnow.nb110_outer
_nb_kernel110_ia32_3dnow.nb110_outerend: 
        ## check if more outer neighborlists remain
        movl  nb110_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel110_ia32_3dnow.nb110_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel110_ia32_3dnow.nb110_threadloop
_nb_kernel110_ia32_3dnow.nb110_end: 
        femms

        movl nb110_nouter(%esp),%eax
        movl nb110_ninner(%esp),%ebx
        movl nb110_outeriter(%ebp),%ecx
        movl nb110_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        addl $152,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret




.globl nb_kernel110nf_ia32_3dnow
.globl _nb_kernel110nf_ia32_3dnow
nb_kernel110nf_ia32_3dnow:      
_nb_kernel110nf_ia32_3dnow:     
.set nb110nf_p_nri, 8
.set nb110nf_iinr, 12
.set nb110nf_jindex, 16
.set nb110nf_jjnr, 20
.set nb110nf_shift, 24
.set nb110nf_shiftvec, 28
.set nb110nf_fshift, 32
.set nb110nf_gid, 36
.set nb110nf_pos, 40
.set nb110nf_faction, 44
.set nb110nf_charge, 48
.set nb110nf_p_facel, 52
.set nb110nf_p_krf, 56
.set nb110nf_p_crf, 60
.set nb110nf_Vc, 64
.set nb110nf_type, 68
.set nb110nf_p_ntype, 72
.set nb110nf_vdwparam, 76
.set nb110nf_Vvdw, 80
.set nb110nf_p_tabscale, 84
.set nb110nf_VFtab, 88
.set nb110nf_invsqrta, 92
.set nb110nf_dvda, 96
.set nb110nf_p_gbtabscale, 100
.set nb110nf_GBtab, 104
.set nb110nf_p_nthreads, 108
.set nb110nf_count, 112
.set nb110nf_mtx, 116
.set nb110nf_outeriter, 120
.set nb110nf_inneriter, 124
.set nb110nf_work, 128
        ## stack offsets for local variables 
.set nb110nf_is3, 0
.set nb110nf_ii3, 4
.set nb110nf_ix, 8
.set nb110nf_iy, 12
.set nb110nf_iz, 16
.set nb110nf_iq, 20
.set nb110nf_vctot, 28
.set nb110nf_Vvdwtot, 36
.set nb110nf_c6, 44
.set nb110nf_c12, 52
.set nb110nf_ntia, 60
.set nb110nf_innerjjnr, 64
.set nb110nf_innerk, 68
.set nb110nf_n, 72                         ## idx for outer loop
.set nb110nf_nn1, 76                       ## number of outer iterations
.set nb110nf_nri, 80
.set nb110nf_facel, 84
.set nb110nf_ntype, 88
.set nb110nf_nouter, 92
.set nb110nf_ninner, 96
        pushl %ebp
        movl %esp,%ebp

        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $100,%esp          ## local stack space 
        femms
        movl nb110nf_p_nri(%ebp),%ecx
        movl nb110nf_p_ntype(%ebp),%edx
        movl nb110nf_p_facel(%ebp),%esi
        movl (%ecx),%ecx
        movl (%edx),%edx
        movl (%esi),%esi
        movl %ecx,nb110nf_nri(%esp)
        movl %edx,nb110nf_ntype(%esp)
        movl %esi,nb110nf_facel(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb110nf_nouter(%esp)
        movl %eax,nb110nf_ninner(%esp)

_nb_kernel110nf_ia32_3dnow.nb110nf_threadloop: 
        movl  nb110nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel110nf_ia32_3dnow.nb110nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel110nf_ia32_3dnow.nb110nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb110nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb110nf_n(%esp)
        movl %ebx,nb110nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel110nf_ia32_3dnow.nb110nf_outerstart
        jmp _nb_kernel110nf_ia32_3dnow.nb110nf_end

_nb_kernel110nf_ia32_3dnow.nb110nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb110nf_nouter(%esp),%ebx
        movl %ebx,nb110nf_nouter(%esp)

_nb_kernel110nf_ia32_3dnow.nb110nf_outer: 
        movl  nb110nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb110nf_is3(%esp)            ## store is3 

        movl  nb110nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movq  (%eax,%ebx,4),%mm0        ## move shX/shY to mm0 and shZ to mm1 
        movd  8(%eax,%ebx,4),%mm1

        movl  nb110nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx            ## ebx=ii 

        movl  nb110nf_charge(%ebp),%edx
        movd  (%edx,%ebx,4),%mm2    ## mm2=charge[ii] 
        pfmul nb110nf_facel(%esp),%mm2
        punpckldq %mm2,%mm2         ## spread to both halves 
        movq  %mm2,nb110nf_iq(%esp)         ## iq =facel*charge[ii] 

        movl  nb110nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb110nf_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb110nf_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb110nf_pos(%ebp),%eax      ## eax = base of pos[] 

        pfadd (%eax,%ebx,4),%mm0    ## ix = shX + posX (and iy too) 
        movd  8(%eax,%ebx,4),%mm3       ## cant use direct memory add for 4 bytes (iz) 
        movl  %ebx,nb110nf_ii3(%esp)
        pfadd %mm3,%mm1
        movq  %mm0,nb110nf_ix(%esp)
        movd  %mm1,nb110nf_iz(%esp)

        ## clear total potential and i forces 
        pxor  %mm7,%mm7
        movq  %mm7,nb110nf_vctot(%esp)
        movq  %mm7,nb110nf_Vvdwtot(%esp)

        movl  nb110nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb110nf_pos(%ebp),%esi
        movl  nb110nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb110nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb110nf_ninner(%esp),%ecx
        movl  %ecx,nb110nf_ninner(%esp)
        movl  %edx,nb110nf_innerk(%esp)      ## number of innerloop atoms 
        addl  $0,%edx
        jge   _nb_kernel110nf_ia32_3dnow.nb110nf_unroll_loop
        jmp   _nb_kernel110nf_ia32_3dnow.nb110nf_finish_inner
_nb_kernel110nf_ia32_3dnow.nb110nf_unroll_loop: 
        ## paired innerloop starts here 
        movl  nb110nf_innerjjnr(%esp),%ecx       ## pointer to jjnr[k] 
        movl  (%ecx),%eax
        movl  4(%ecx),%ebx           ## eax/ebx=jnr 
        addl $8,nb110nf_innerjjnr(%esp)             ## advance pointer (unrolled 2) 
        prefetch 16(%ecx)            ## prefetch data - trial and error says 16 is best 

        movl nb110nf_charge(%ebp),%ecx     ## base of charge[] 
        movq nb110nf_iq(%esp),%mm5
        movd (%ecx,%eax,4),%mm3      ## charge[jnr1] 
        punpckldq (%ecx,%ebx,4),%mm3     ## move charge 2 to high part of mm3 
        pfmul %mm5,%mm3              ## mm3 now has qq for both particles 

        movl nb110nf_type(%ebp),%ecx
        movl (%ecx,%eax,4),%edx          ## type [jnr1] 
        movl (%ecx,%ebx,4),%ecx      ## type [jnr2] 

        movl nb110nf_vdwparam(%ebp),%esi                ## base of vdwparam  
        shll %edx
        shll %ecx
        addl nb110nf_ntia(%esp),%edx         ## tja = ntia + 2*type 
        addl nb110nf_ntia(%esp),%ecx

        movq (%esi,%edx,4),%mm5         ## mm5 = 1st c6 / c12           
        movq (%esi,%ecx,4),%mm7         ## mm7 = 2nd c6 / c12   
        movq %mm5,%mm6
        punpckldq %mm7,%mm5             ## mm5 = 1st c6 / 2nd c6 
        punpckhdq %mm7,%mm6             ## mm6 = 1st c12 / 2nd c12 
        movq %mm5,nb110nf_c6(%esp)
        movq %mm6,nb110nf_c12(%esp)

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        movl  nb110nf_pos(%ebp),%esi

        movq  nb110nf_ix(%esp),%mm0
        movd  nb110nf_iz(%esp),%mm1
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
        punpckldq %mm6,%mm4             ## now 4 has rsq and 0 the seed for both pairs 
    movq %mm0,%mm2                      ## amd 3dnow N-R iteration to get full precision 
        pfmul %mm0,%mm0
    pfrsqit1 %mm4,%mm0
    pfrcpit2 %mm2,%mm0
        movq %mm0,%mm1
        pfmul %mm0,%mm0
        ## mm0 now contains invsq, and mm1 invsqrt 
        ## do potential and fscal 
        movq %mm0,%mm4
        pfmul %mm0,%mm4
        pfmul %mm0,%mm4                 ## mm4=rinvsix 
        movq  %mm4,%mm5
        pfmul %mm5,%mm5             ## mm5=rinvtwelve 

        pfmul %mm1,%mm3         ## mm3 has vcoul for both interactions 
        pfmul nb110nf_c12(%esp),%mm5
        pfmul nb110nf_c6(%esp),%mm4
        movq %mm5,%mm6  ## mm6 is Vvdw12-Vvdw6  
        pfsub %mm4,%mm6
        ## update vctot 
        pfadd nb110nf_vctot(%esp),%mm3        ## add the earlier value 
        movq %mm3,nb110nf_vctot(%esp)         ## store the sum       
        ## update Vvdwtot 
        pfadd nb110nf_Vvdwtot(%esp),%mm6        ## add the earlier value 
        movq %mm6,nb110nf_Vvdwtot(%esp)         ## store the sum       

        ## should we do one more iteration? 
        subl $2,nb110nf_innerk(%esp)
        jl    _nb_kernel110nf_ia32_3dnow.nb110nf_finish_inner
        jmp   _nb_kernel110nf_ia32_3dnow.nb110nf_unroll_loop
_nb_kernel110nf_ia32_3dnow.nb110nf_finish_inner: 
        andl $1,nb110nf_innerk(%esp)
        jnz  _nb_kernel110nf_ia32_3dnow.nb110nf_single_inner
        jmp  _nb_kernel110nf_ia32_3dnow.nb110nf_updateouterdata
_nb_kernel110nf_ia32_3dnow.nb110nf_single_inner: 
        ## a single j particle iteration here - compare with the unrolled code for comments 
        movl  nb110nf_innerjjnr(%esp),%eax
        movl  (%eax),%eax       ## eax=jnr offset 

        movl nb110nf_charge(%ebp),%ecx
        movd nb110nf_iq(%esp),%mm5
        movd (%ecx,%eax,4),%mm3
        pfmul %mm5,%mm3         ## mm3=qq 

        movl nb110nf_vdwparam(%ebp),%esi
        movl nb110nf_type(%ebp),%ecx
        movl (%ecx,%eax,4),%edx          ## type [jnr1] 
        shll %edx
        addl nb110nf_ntia(%esp),%edx         ## tja = ntia + 2*type 
        movd (%esi,%edx,4),%mm5         ## mm5 = 1st c6                 
        movq %mm5,nb110nf_c6(%esp)
        movd 4(%esi,%edx,4),%mm5        ## mm5 = 1st c12                
        movq %mm5,nb110nf_c12(%esp)


        movl  nb110nf_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        movq  nb110nf_ix(%esp),%mm0
        movd  nb110nf_iz(%esp),%mm1
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
        pfrcpit2 %mm2,%mm0      ## mm1=invsqrt 
        movq  %mm0,%mm1
        pfmul %mm0,%mm0         ## mm0=invsq 
        ## calculate potentials and scalar force 
        movq %mm0,%mm4
        pfmul %mm0,%mm4
        pfmul %mm0,%mm4                 ## mm4=rinvsix 
        movq  %mm4,%mm5
        pfmul %mm5,%mm5             ## mm5=rinvtwelve 

        pfmul %mm1,%mm3         ## mm3 has vcoul for both interactions 
        pfmul nb110nf_c12(%esp),%mm5
        pfmul nb110nf_c6(%esp),%mm4
        movq %mm5,%mm6  ## mm6 is Vvdw12-Vvdw6  
        pfsub %mm4,%mm6
        ## update vctot 
        pfadd nb110nf_vctot(%esp),%mm3
        movq %mm3,nb110nf_vctot(%esp)
        ## update Vvdwtot 
        pfadd nb110nf_Vvdwtot(%esp),%mm6        ## add the earlier value 
        movq %mm6,nb110nf_Vvdwtot(%esp)         ## store the sum       

_nb_kernel110nf_ia32_3dnow.nb110nf_updateouterdata: 
        ## get n from stack
        movl nb110nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb110nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        movq  nb110nf_vctot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb110nf_Vc(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment vc[gid] 

        movq  nb110nf_Vvdwtot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb110nf_Vvdw(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment Vvdw[gid] 

        ## finish if last 
        movl nb110nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel110nf_ia32_3dnow.nb110nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb110nf_n(%esp)
        jmp _nb_kernel110nf_ia32_3dnow.nb110nf_outer
_nb_kernel110nf_ia32_3dnow.nb110nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb110nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel110nf_ia32_3dnow.nb110nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel110nf_ia32_3dnow.nb110nf_threadloop
_nb_kernel110nf_ia32_3dnow.nb110nf_end: 
        femms

        movl nb110nf_nouter(%esp),%eax
        movl nb110nf_ninner(%esp),%ebx
        movl nb110nf_outeriter(%ebp),%ecx
        movl nb110nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        addl $100,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret

