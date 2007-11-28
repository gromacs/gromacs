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



.globl nb_kernel030_ia32_3dnow
.globl _nb_kernel030_ia32_3dnow
nb_kernel030_ia32_3dnow:        
_nb_kernel030_ia32_3dnow:       
.set nb030_p_nri, 8
.set nb030_iinr, 12
.set nb030_jindex, 16
.set nb030_jjnr, 20
.set nb030_shift, 24
.set nb030_shiftvec, 28
.set nb030_fshift, 32
.set nb030_gid, 36
.set nb030_pos, 40
.set nb030_faction, 44
.set nb030_charge, 48
.set nb030_p_facel, 52
.set nb030_p_krf, 56
.set nb030_p_crf, 60
.set nb030_Vc, 64
.set nb030_type, 68
.set nb030_p_ntype, 72
.set nb030_vdwparam, 76
.set nb030_Vvdw, 80
.set nb030_p_tabscale, 84
.set nb030_VFtab, 88
.set nb030_invsqrta, 92
.set nb030_dvda, 96
.set nb030_p_gbtabscale, 100
.set nb030_GBtab, 104
.set nb030_p_nthreads, 108
.set nb030_count, 112
.set nb030_mtx, 116
.set nb030_outeriter, 120
.set nb030_inneriter, 124
.set nb030_work, 128
        ## stack offsets for local variables 
.set nb030_is3, 0
.set nb030_ii3, 4
.set nb030_ix, 8
.set nb030_iy, 12
.set nb030_iz, 16
.set nb030_Vvdwtot, 20
.set nb030_c6, 28
.set nb030_c12, 36
.set nb030_two, 44
.set nb030_n1, 52
.set nb030_tsc, 60
.set nb030_ntia, 68
.set nb030_innerjjnr, 72
.set nb030_innerk, 76
.set nb030_fix, 80
.set nb030_fiy, 84
.set nb030_fiz, 88
.set nb030_dx1, 92
.set nb030_dy1, 96
.set nb030_dz1, 100
.set nb030_dx2, 104
.set nb030_dy2, 108
.set nb030_dz2, 112
.set nb030_n, 116                           ## idx for outer loop
.set nb030_nn1, 120                         ## number of outer iterations
.set nb030_nri, 124
.set nb030_ntype, 128
.set nb030_nouter, 132
.set nb030_ninner, 136

        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $140,%esp          ## local stack space 
        femms
        ## move data to local stack  
        movl $0x40000000,%eax
        movl %eax,nb030_two(%esp)
        movl %eax,nb030_two+4(%esp)

        movl nb030_p_nri(%ebp),%ecx
        movl nb030_p_ntype(%ebp),%edx
        movl nb030_p_tabscale(%ebp),%esi
        movl (%ecx),%ecx
        movl (%edx),%edx
        movl %ecx,nb030_nri(%esp)
        movl %edx,nb030_ntype(%esp)

        movd  (%esi),%mm3
        punpckldq %mm3,%mm3
        movq  %mm3,nb030_tsc(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb030_nouter(%esp)
        movl %eax,nb030_ninner(%esp)

_nb_kernel030_ia32_3dnow.nb030_threadloop: 
        movl  nb030_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel030_ia32_3dnow.nb030_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel030_ia32_3dnow.nb030_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb030_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb030_n(%esp)
        movl %ebx,nb030_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel030_ia32_3dnow.nb030_outerstart
        jmp _nb_kernel030_ia32_3dnow.nb030_end

_nb_kernel030_ia32_3dnow.nb030_outerstart: 
        ## ebx contains number of outer iterations
        addl nb030_nouter(%esp),%ebx
        movl %ebx,nb030_nouter(%esp)


_nb_kernel030_ia32_3dnow.nb030_outer: 
        movl  nb030_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb030_is3(%esp)      ## store is3 

        movl  nb030_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movq  (%eax,%ebx,4),%mm0        ## move shX/shY to mm0 and shZ to mm1 
        movd  8(%eax,%ebx,4),%mm1

        movl  nb030_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx=ii 

        movl  nb030_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb030_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb030_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb030_pos(%ebp),%eax      ## eax = base of pos[] 

        pfadd (%eax,%ebx,4),%mm0    ## ix = shX + posX (and iy too) 
        movd  8(%eax,%ebx,4),%mm3       ## cant use direct memory add for 4 bytes (iz) 
        movl  %ebx,nb030_ii3(%esp)
        pfadd %mm3,%mm1
        movq  %mm0,nb030_ix(%esp)
        movd  %mm1,nb030_iz(%esp)

        ## clear total potential and i forces 
        pxor  %mm7,%mm7
        movq  %mm7,nb030_Vvdwtot(%esp)
        movq  %mm7,nb030_fix(%esp)
        movd  %mm7,nb030_fiz(%esp)

        movl  nb030_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb030_pos(%ebp),%esi
        movl  nb030_faction(%ebp),%edi
        movl  nb030_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb030_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb030_ninner(%esp),%ecx
        movl  %ecx,nb030_ninner(%esp)
        movl  %edx,nb030_innerk(%esp)      ## number of innerloop atoms 
        addl  $0,%edx
        jge   _nb_kernel030_ia32_3dnow.nb030_unroll_loop
        jmp   _nb_kernel030_ia32_3dnow.nb030_finish_inner
_nb_kernel030_ia32_3dnow.nb030_unroll_loop: 
        ## paired innerloop starts here 
        movl  nb030_innerjjnr(%esp),%ecx       ## pointer to jjnr[k] 
        movl  (%ecx),%eax
        movl  4(%ecx),%ebx           ## eax/ebx=jnr 
        addl $8,nb030_innerjjnr(%esp)             ## advance pointer (unrolled 2) 
        prefetch 16(%ecx)            ## prefetch data - trial and error says 16 is best  

        movl nb030_type(%ebp),%ecx
        movl (%ecx,%eax,4),%edx          ## type [jnr1] 
        movl (%ecx,%ebx,4),%ecx      ## type [jnr2] 

        movl nb030_vdwparam(%ebp),%esi          ## base of vdwparam  
        shll %edx
        shll %ecx
        addl nb030_ntia(%esp),%edx           ## tja = ntia + 2*type  
        addl nb030_ntia(%esp),%ecx

        movq (%esi,%edx,4),%mm5         ## mm5 = 1st c6 / c12           
        movq (%esi,%ecx,4),%mm7         ## mm7 = 2nd c6 / c12   
        movq %mm5,%mm6
        punpckldq %mm7,%mm5             ## mm5 = 1st c6 / 2nd c6 
        punpckhdq %mm7,%mm6             ## mm6 = 1st c12 / 2nd c12 
        movq %mm5,nb030_c6(%esp)
        movq %mm6,nb030_c12(%esp)

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        movl  nb030_pos(%ebp),%esi

        movq  nb030_ix(%esp),%mm0
        movd  nb030_iz(%esp),%mm1
        movq  (%esi,%eax,4),%mm4     ## fetch first j coordinates 
        movd  8(%esi,%eax,4),%mm5
        pfsubr %mm0,%mm4             ## dr = ir - jr  
        pfsubr %mm1,%mm5
        movq  %mm4,nb030_dx1(%esp)           ## store dr 
        movd  %mm5,nb030_dz1(%esp)
        pfmul %mm4,%mm4              ## square dx,dy,dz                          
        pfmul %mm5,%mm5
        pfacc %mm5,%mm4              ## accumulate to get dx*dx+ dy*dy+ dz*dz 
        pfacc %mm5,%mm4              ## first rsq in lower mm4 

        movq  (%esi,%ebx,4),%mm6     ## fetch second j coordinates  
        movd  8(%esi,%ebx,4),%mm7

        pfsubr %mm0,%mm6             ## dr = ir - jr  
        pfsubr %mm1,%mm7
        movq  %mm6,nb030_dx2(%esp)           ## store dr 
        movd  %mm7,nb030_dz2(%esp)
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
        pfmul %mm0,%mm4
        movq %mm4,%mm1
        ## mm0 is invsqrt, and mm1 r 
        ## do potential and fscal 
        pfmul nb030_tsc(%esp),%mm1      ## mm1=rt 
        pf2iw %mm1,%mm4
        movq %mm4,nb030_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm1              ## now mm1 is eps and mm4 is n0 

        movq %mm1,%mm2
        pfmul %mm2,%mm2 ## mm1 is eps, mm2 is eps2 

        movl nb030_VFtab(%ebp),%edx
        ## dispersion table 
        movl nb030_n1(%esp),%ecx
        shll $3,%ecx
        ## load all the table values we need 
        movd (%edx,%ecx,4),%mm4
        movd 4(%edx,%ecx,4),%mm5
        movd 8(%edx,%ecx,4),%mm6
        movd 12(%edx,%ecx,4),%mm7
        movl nb030_n1+4(%esp),%ecx
        shll $3,%ecx
        punpckldq (%edx,%ecx,4),%mm4
        punpckldq 4(%edx,%ecx,4),%mm5
        punpckldq 8(%edx,%ecx,4),%mm6
        punpckldq 12(%edx,%ecx,4),%mm7
        pfmul %mm1,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 
        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 
        pfmul nb030_two(%esp),%mm7      ## two*Heps2 
        pfadd %mm6,%mm7
        pfadd %mm5,%mm7 ## mm7=FF 
        pfmul %mm1,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV     

        movq nb030_c6(%esp),%mm4
        pfmul %mm4,%mm7 ## fijD 
        pfmul %mm4,%mm5 ## Vvdw6            
        movq %mm7,%mm3  ## add to fscal 

        ## update Vvdwtot to release mm5! 
        pfadd nb030_Vvdwtot(%esp),%mm5        ## add the earlier value 
        movq %mm5,nb030_Vvdwtot(%esp)         ## store the sum       

        ## repulsion table 
        movl nb030_n1(%esp),%ecx
        shll $3,%ecx
        ## load all the table values we need 
        movd 16(%edx,%ecx,4),%mm4
        movd 20(%edx,%ecx,4),%mm5
        movd 24(%edx,%ecx,4),%mm6
        movd 28(%edx,%ecx,4),%mm7
        movl nb030_n1+4(%esp),%ecx
        shll $3,%ecx
        punpckldq 16(%edx,%ecx,4),%mm4
        punpckldq 20(%edx,%ecx,4),%mm5
        punpckldq 24(%edx,%ecx,4),%mm6
        punpckldq 28(%edx,%ecx,4),%mm7

        pfmul %mm1,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 
        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 
        pfmul nb030_two(%esp),%mm7      ## two*Heps2 
        pfadd %mm6,%mm7
        pfadd %mm5,%mm7 ## mm7=FF 
        pfmul %mm1,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        movq nb030_c12(%esp),%mm6
        pfmul %mm6,%mm7 ## fijR 
        pfmul %mm6,%mm5 ## Vvdw12 
        pfadd %mm7,%mm3 ## total fscal fijD+ fijR 

        ## change sign of mm3 
    pxor %mm1,%mm1
        pfsub %mm3,%mm1
        pfmul nb030_tsc(%esp),%mm1
        pfmul %mm1,%mm0   ## mm0 is total fscal now     

        prefetchw nb030_dx1(%esp)       ## prefetch i forces to cache 

        ## spread fscalar to both positions 
        movq %mm0,%mm1
        punpckldq %mm0,%mm0
        punpckhdq %mm1,%mm1

        ## calc vector force 
        prefetchw (%edi,%eax,4) ## prefetch the 1st faction to cache 
        movq nb030_dx1(%esp),%mm2       ## fetch dr 
        movd nb030_dz1(%esp),%mm3

        ## update Vvdwtot 
        pfadd nb030_Vvdwtot(%esp),%mm5        ## add the earlier value 
        movq %mm5,nb030_Vvdwtot(%esp)         ## store the sum       

        prefetchw (%edi,%ebx,4) ## prefetch the 2nd faction to cache 
        pfmul %mm0,%mm2         ## mult by fs  
        pfmul %mm0,%mm3

        movq nb030_dx2(%esp),%mm4       ## fetch dr 
        movd nb030_dz2(%esp),%mm5
        pfmul %mm1,%mm4         ## mult by fs  
        pfmul %mm1,%mm5
        ## update i forces 

        movq nb030_fix(%esp),%mm0
        movd nb030_fiz(%esp),%mm1
        pfadd %mm2,%mm0
        pfadd %mm3,%mm1

        pfadd %mm4,%mm0
        pfadd %mm5,%mm1
        movq %mm0,nb030_fix(%esp)
        movd %mm1,nb030_fiz(%esp)
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
        subl $2,nb030_innerk(%esp)
        jl    _nb_kernel030_ia32_3dnow.nb030_finish_inner
        jmp   _nb_kernel030_ia32_3dnow.nb030_unroll_loop
_nb_kernel030_ia32_3dnow.nb030_finish_inner: 
        andl $1,nb030_innerk(%esp)
        jnz  _nb_kernel030_ia32_3dnow.nb030_single_inner
        jmp  _nb_kernel030_ia32_3dnow.nb030_updateouterdata
_nb_kernel030_ia32_3dnow.nb030_single_inner: 
        ## a single j particle iteration here - compare with the unrolled code for comments 
        movl  nb030_innerjjnr(%esp),%eax
        movl  (%eax),%eax       ## eax=jnr offset 

        movl nb030_vdwparam(%ebp),%esi
        movl nb030_type(%ebp),%ecx
        movl (%ecx,%eax,4),%edx          ## type [jnr1] 
        shll %edx
        addl nb030_ntia(%esp),%edx           ## tja = ntia + 2*type 
        movd (%esi,%edx,4),%mm5         ## mm5 = 1st c6                 
        movq %mm5,nb030_c6(%esp)
        movd 4(%esi,%edx,4),%mm5        ## mm5 = 1st c12                
        movq %mm5,nb030_c12(%esp)

        movl  nb030_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        movq  nb030_ix(%esp),%mm0
        movd  nb030_iz(%esp),%mm1
        movq  (%esi,%eax,4),%mm4
        movd  8(%esi,%eax,4),%mm5
        pfsubr %mm0,%mm4
        pfsubr %mm1,%mm5
        movq  %mm4,nb030_dx1(%esp)
        pfmul %mm4,%mm4
        movd  %mm5,nb030_dz1(%esp)
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
        ## mm0 is invsqrt, and mm1 r 

        ## calculate potentials and scalar force 
        pfmul nb030_tsc(%esp),%mm1      ## mm1=rt 
        pf2iw %mm1,%mm4
        movd %mm4,nb030_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm1              ## now mm1 is eps and mm4 is n0 

        movq %mm1,%mm2
        pfmul %mm2,%mm2 ## mm1 is eps, mm2 is eps2 

        movl nb030_VFtab(%ebp),%edx
        movl nb030_n1(%esp),%ecx
        shll $3,%ecx
        ## dispersion table 
        ## load all the table values we need 
        movd (%edx,%ecx,4),%mm4
        movd 4(%edx,%ecx,4),%mm5
        movd 8(%edx,%ecx,4),%mm6
        movd 12(%edx,%ecx,4),%mm7
        pfmul %mm1,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 
        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 
        pfmul nb030_two(%esp),%mm7      ## two*Heps2 
        pfadd %mm6,%mm7
        pfadd %mm5,%mm7 ## mm7=FF 
        pfmul %mm1,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV     

        movq nb030_c6(%esp),%mm4
        pfmul %mm4,%mm7 ## fijD 
        pfmul %mm4,%mm5 ## Vvdw6            
        movq %mm7,%mm3  ## add to fscal 

        ## update Vvdwtot to release mm5! 
        pfadd nb030_Vvdwtot(%esp),%mm5        ## add the earlier value 
        movq %mm5,nb030_Vvdwtot(%esp)         ## store the sum       

        ## repulsion table 
        ## load all the table values we need 
        movd 16(%edx,%ecx,4),%mm4
        movd 20(%edx,%ecx,4),%mm5
        movd 24(%edx,%ecx,4),%mm6
        movd 28(%edx,%ecx,4),%mm7

        pfmul %mm1,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 
        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 
        pfmul nb030_two(%esp),%mm7      ## two*Heps2 
        pfadd %mm6,%mm7
        pfadd %mm5,%mm7 ## mm7=FF 
        pfmul %mm1,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        movq nb030_c12(%esp),%mm6
        pfmul %mm6,%mm7 ## fijR 
        pfmul %mm6,%mm5 ## Vvdw12 
        pfadd %mm7,%mm3 ## total fscal fijC+ fijD+ fijR 

        ## change sign of mm3 
    pxor %mm1,%mm1
        pfsub %mm3,%mm1
        pfmul nb030_tsc(%esp),%mm0
        pfmul %mm1,%mm0   ## mm0 is total fscal now     

        ## update Vvdwtot 
        pfadd nb030_Vvdwtot(%esp),%mm5        ## add the earlier value 
        movq %mm5,nb030_Vvdwtot(%esp)         ## store the sum       

        ## spread fscalar to both positions 
        punpckldq %mm0,%mm0
        ## calc vectorial force 
        prefetchw (%edi,%eax,4) ## prefetch faction to cache  
        movq nb030_dx1(%esp),%mm2
        movd nb030_dz1(%esp),%mm3

        pfmul %mm0,%mm2
        pfmul %mm0,%mm3

        ## update i particle force 
        movq nb030_fix(%esp),%mm0
        movd nb030_fiz(%esp),%mm1
        pfadd %mm2,%mm0
        pfadd %mm3,%mm1
        movq %mm0,nb030_fix(%esp)
        movd %mm1,nb030_fiz(%esp)
        ## update j particle force 
        movq (%edi,%eax,4),%mm0
        movd 8(%edi,%eax,4),%mm1
        pfsub %mm2,%mm0
        pfsub %mm3,%mm1
        movq %mm0,(%edi,%eax,4)
        movd %mm1,8(%edi,%eax,4)
        ## done! 
_nb_kernel030_ia32_3dnow.nb030_updateouterdata: 
        movl  nb030_ii3(%esp),%ecx

        movq  (%edi,%ecx,4),%mm6       ## increment i force 
        movd  8(%edi,%ecx,4),%mm7
        pfadd nb030_fix(%esp),%mm6
        pfadd nb030_fiz(%esp),%mm7
        movq  %mm6,(%edi,%ecx,4)
        movd  %mm7,8(%edi,%ecx,4)

        movl  nb030_fshift(%ebp),%ebx      ## increment fshift force 
        movl  nb030_is3(%esp),%edx

        movq  (%ebx,%edx,4),%mm6
        movd  8(%ebx,%edx,4),%mm7
        pfadd nb030_fix(%esp),%mm6
        pfadd nb030_fiz(%esp),%mm7
        movq  %mm6,(%ebx,%edx,4)
        movd  %mm7,8(%ebx,%edx,4)

        ## get n from stack
        movl nb030_n(%esp),%esi
        ## get group index for i particle 
        movl  nb030_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        movq  nb030_Vvdwtot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb030_Vvdw(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment Vvdw[gid] 

       ## finish if last 
        movl nb030_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel030_ia32_3dnow.nb030_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb030_n(%esp)
        jmp _nb_kernel030_ia32_3dnow.nb030_outer
_nb_kernel030_ia32_3dnow.nb030_outerend: 
        ## check if more outer neighborlists remain
        movl  nb030_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel030_ia32_3dnow.nb030_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel030_ia32_3dnow.nb030_threadloop
_nb_kernel030_ia32_3dnow.nb030_end: 
        femms

        movl nb030_nouter(%esp),%eax
        movl nb030_ninner(%esp),%ebx
        movl nb030_outeriter(%ebp),%ecx
        movl nb030_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        addl $140,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret






.globl nb_kernel030nf_ia32_3dnow
.globl _nb_kernel030nf_ia32_3dnow
nb_kernel030nf_ia32_3dnow:      
_nb_kernel030nf_ia32_3dnow:     
.set nb030nf_p_nri, 8
.set nb030nf_iinr, 12
.set nb030nf_jindex, 16
.set nb030nf_jjnr, 20
.set nb030nf_shift, 24
.set nb030nf_shiftvec, 28
.set nb030nf_fshift, 32
.set nb030nf_gid, 36
.set nb030nf_pos, 40
.set nb030nf_faction, 44
.set nb030nf_charge, 48
.set nb030nf_p_facel, 52
.set nb030nf_p_krf, 56
.set nb030nf_p_crf, 60
.set nb030nf_Vc, 64
.set nb030nf_type, 68
.set nb030nf_p_ntype, 72
.set nb030nf_vdwparam, 76
.set nb030nf_Vvdw, 80
.set nb030nf_p_tabscale, 84
.set nb030nf_VFtab, 88
.set nb030nf_invsqrta, 92
.set nb030nf_dvda, 96
.set nb030nf_p_gbtabscale, 100
.set nb030nf_GBtab, 104
.set nb030nf_p_nthreads, 108
.set nb030nf_count, 112
.set nb030nf_mtx, 116
.set nb030nf_outeriter, 120
.set nb030nf_inneriter, 124
.set nb030nf_work, 128
        ## stack offsets for local variables 
.set nb030nf_is3, 0
.set nb030nf_ii3, 4
.set nb030nf_ix, 8
.set nb030nf_iy, 12
.set nb030nf_iz, 16
.set nb030nf_Vvdwtot, 20
.set nb030nf_c6, 28
.set nb030nf_c12, 36
.set nb030nf_n1, 44
.set nb030nf_tsc, 52
.set nb030nf_ntia, 60
.set nb030nf_innerjjnr, 64
.set nb030nf_innerk, 68
.set nb030nf_n, 72                         ## idx for outer loop
.set nb030nf_nn1, 76                       ## number of outer iterations
.set nb030nf_nri, 80
.set nb030nf_ntype, 84
.set nb030nf_nouter, 88
.set nb030nf_ninner, 92

        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $96,%esp           ## local stack space 
        femms
        ## move data to local stack  
        movl nb030nf_p_nri(%ebp),%ecx
        movl nb030nf_p_ntype(%ebp),%edx
        movl nb030nf_p_tabscale(%ebp),%esi
        movl (%ecx),%ecx
        movl (%edx),%edx
        movl %ecx,nb030nf_nri(%esp)
        movl %edx,nb030nf_ntype(%esp)

        movd  (%esi),%mm3
        punpckldq %mm3,%mm3
        movq  %mm3,nb030nf_tsc(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb030nf_nouter(%esp)
        movl %eax,nb030nf_ninner(%esp)

_nb_kernel030nf_ia32_3dnow.nb030nf_threadloop: 
        movl  nb030nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel030nf_ia32_3dnow.nb030nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel030nf_ia32_3dnow.nb030nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb030nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb030nf_n(%esp)
        movl %ebx,nb030nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel030nf_ia32_3dnow.nb030nf_outerstart
        jmp _nb_kernel030nf_ia32_3dnow.nb030nf_end

_nb_kernel030nf_ia32_3dnow.nb030nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb030nf_nouter(%esp),%ebx
        movl %ebx,nb030nf_nouter(%esp)

_nb_kernel030nf_ia32_3dnow.nb030nf_outer: 
        movl  nb030nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb030nf_is3(%esp)            ## store is3 

        movl  nb030nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movq  (%eax,%ebx,4),%mm0        ## move shX/shY to mm0 and shZ to mm1 
        movd  8(%eax,%ebx,4),%mm1

        movl  nb030nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx            ## ebx=ii 

        movl  nb030nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb030nf_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb030nf_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb030nf_pos(%ebp),%eax      ## eax = base of pos[] 

        pfadd (%eax,%ebx,4),%mm0    ## ix = shX + posX (and iy too) 
        movd  8(%eax,%ebx,4),%mm3       ## cant use direct memory add for 4 bytes (iz) 
        movl  %ebx,nb030nf_ii3(%esp)
        pfadd %mm3,%mm1
        movq  %mm0,nb030nf_ix(%esp)
        movd  %mm1,nb030nf_iz(%esp)

        ## clear total potential 
        pxor  %mm7,%mm7
        movq  %mm7,nb030nf_Vvdwtot(%esp)

        movl  nb030nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb030nf_pos(%ebp),%esi
        movl  nb030nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb030nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb030nf_ninner(%esp),%ecx
        movl  %ecx,nb030nf_ninner(%esp)
        movl  %edx,nb030nf_innerk(%esp)      ## number of innerloop atoms 
        addl  $0,%edx
        jge   _nb_kernel030nf_ia32_3dnow.nb030nf_unroll_loop
        jmp   _nb_kernel030nf_ia32_3dnow.nb030nf_finish_inner
_nb_kernel030nf_ia32_3dnow.nb030nf_unroll_loop: 
        ## paired innerloop starts here 
        movl  nb030nf_innerjjnr(%esp),%ecx       ## pointer to jjnr[k] 
        movl  (%ecx),%eax
        movl  4(%ecx),%ebx           ## eax/ebx=jnr 
        addl $8,nb030nf_innerjjnr(%esp)             ## advance pointer (unrolled 2) 
        prefetch 16(%ecx)            ## prefetch data - trial and error says 16 is best  

        movl nb030nf_type(%ebp),%ecx
        movl (%ecx,%eax,4),%edx          ## type [jnr1] 
        movl (%ecx,%ebx,4),%ecx      ## type [jnr2] 

        movl nb030nf_vdwparam(%ebp),%esi                ## base of vdwparam  
        shll %edx
        shll %ecx
        addl nb030nf_ntia(%esp),%edx         ## tja = ntia + 2*type  
        addl nb030nf_ntia(%esp),%ecx

        movq (%esi,%edx,4),%mm5         ## mm5 = 1st c6 / c12           
        movq (%esi,%ecx,4),%mm7         ## mm7 = 2nd c6 / c12   
        movq %mm5,%mm6
        punpckldq %mm7,%mm5             ## mm5 = 1st c6 / 2nd c6 
        punpckhdq %mm7,%mm6             ## mm6 = 1st c12 / 2nd c12 
        movq %mm5,nb030nf_c6(%esp)
        movq %mm6,nb030nf_c12(%esp)

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        movl  nb030nf_pos(%ebp),%esi

        movq  nb030nf_ix(%esp),%mm0
        movd  nb030nf_iz(%esp),%mm1
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
        pfmul %mm0,%mm4
        movq %mm4,%mm1
        ## mm0 is invsqrt, and mm1 r 
        ## do potential and fscal 
        pfmul nb030nf_tsc(%esp),%mm1    ## mm1=rt 
        pf2iw %mm1,%mm4
        movq %mm4,nb030nf_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm1              ## now mm1 is eps and mm4 is n0 

        movq %mm1,%mm2
        pfmul %mm2,%mm2 ## mm1 is eps, mm2 is eps2 

        movl nb030nf_VFtab(%ebp),%edx
        ## dispersion table 
        movl nb030nf_n1(%esp),%ecx
        shll $3,%ecx
        ## load all the table values we need 
        movd (%edx,%ecx,4),%mm4
        movd 4(%edx,%ecx,4),%mm5
        movd 8(%edx,%ecx,4),%mm6
        movd 12(%edx,%ecx,4),%mm7
        movl nb030nf_n1+4(%esp),%ecx
        shll $3,%ecx
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

        movq nb030nf_c6(%esp),%mm4
        pfmul %mm4,%mm5 ## Vvdw6        
        ## update Vvdwtot to release mm5! 
        pfadd nb030nf_Vvdwtot(%esp),%mm5        ## add the earlier value 
        movq %mm5,nb030nf_Vvdwtot(%esp)         ## store the sum       

        ## repulsion table 
        movl nb030nf_n1(%esp),%ecx
        shll $3,%ecx
        ## load all the table values we need 
        movd 16(%edx,%ecx,4),%mm4
        movd 20(%edx,%ecx,4),%mm5
        movd 24(%edx,%ecx,4),%mm6
        movd 28(%edx,%ecx,4),%mm7
        movl nb030nf_n1+4(%esp),%ecx
        shll $3,%ecx
        punpckldq 16(%edx,%ecx,4),%mm4
        punpckldq 20(%edx,%ecx,4),%mm5
        punpckldq 24(%edx,%ecx,4),%mm6
        punpckldq 28(%edx,%ecx,4),%mm7

        pfmul %mm1,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 
        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 
        pfmul %mm1,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        movq nb030nf_c12(%esp),%mm6
        pfmul %mm6,%mm5 ## Vvdw12 
        ## update Vvdwtot 
        pfadd nb030nf_Vvdwtot(%esp),%mm5        ## add the earlier value 
        movq %mm5,nb030nf_Vvdwtot(%esp)         ## store the sum       

        ## should we do one more iteration? 
        subl $2,nb030nf_innerk(%esp)
        jl    _nb_kernel030nf_ia32_3dnow.nb030nf_finish_inner
        jmp   _nb_kernel030nf_ia32_3dnow.nb030nf_unroll_loop
_nb_kernel030nf_ia32_3dnow.nb030nf_finish_inner: 
        andl $1,nb030nf_innerk(%esp)
        jnz  _nb_kernel030nf_ia32_3dnow.nb030nf_single_inner
        jmp  _nb_kernel030nf_ia32_3dnow.nb030nf_updateouterdata
_nb_kernel030nf_ia32_3dnow.nb030nf_single_inner: 
        ## a single j particle iteration here - compare with the unrolled code for comments 
        movl  nb030nf_innerjjnr(%esp),%eax
        movl  (%eax),%eax       ## eax=jnr offset 

        movl nb030nf_vdwparam(%ebp),%esi
        movl nb030nf_type(%ebp),%ecx
        movl (%ecx,%eax,4),%edx          ## type [jnr1] 
        shll %edx
        addl nb030nf_ntia(%esp),%edx         ## tja = ntia + 2*type 
        movd (%esi,%edx,4),%mm5         ## mm5 = 1st c6                 
        movq %mm5,nb030nf_c6(%esp)
        movd 4(%esi,%edx,4),%mm5        ## mm5 = 1st c12                
        movq %mm5,nb030nf_c12(%esp)

        movl  nb030nf_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        movq  nb030nf_ix(%esp),%mm0
        movd  nb030nf_iz(%esp),%mm1
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
        ## mm0 is invsqrt, and mm1 r 

        ## calculate potentials and scalar force 
        pfmul nb030nf_tsc(%esp),%mm1    ## mm1=rt 
        pf2iw %mm1,%mm4
        movd %mm4,nb030nf_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm1              ## now mm1 is eps and mm4 is n0 

        movq %mm1,%mm2
        pfmul %mm2,%mm2 ## mm1 is eps, mm2 is eps2 

        movl nb030nf_VFtab(%ebp),%edx
        movl nb030nf_n1(%esp),%ecx
        shll $3,%ecx
        ## dispersion table 
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

        movq nb030nf_c6(%esp),%mm4
        pfmul %mm4,%mm5 ## Vvdw6            
        ## update Vvdwtot to release mm5! 
        pfadd nb030nf_Vvdwtot(%esp),%mm5        ## add the earlier value 
        movq %mm5,nb030nf_Vvdwtot(%esp)         ## store the sum       

        ## repulsion table 
        ## load all the table values we need 
        movd 16(%edx,%ecx,4),%mm4
        movd 20(%edx,%ecx,4),%mm5
        movd 24(%edx,%ecx,4),%mm6
        movd 28(%edx,%ecx,4),%mm7

        pfmul %mm1,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 
        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 
        pfmul %mm1,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        movq nb030nf_c12(%esp),%mm6
        pfmul %mm6,%mm5 ## Vvdw12 
        ## update Vvdwtot 
        pfadd nb030nf_Vvdwtot(%esp),%mm5        ## add the earlier value 
        movq %mm5,nb030nf_Vvdwtot(%esp)         ## store the sum       

_nb_kernel030nf_ia32_3dnow.nb030nf_updateouterdata: 
        ## get n from stack
        movl nb030nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb030nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        movq  nb030nf_Vvdwtot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb030nf_Vvdw(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment Vvdw[gid] 

       ## finish if last 
        movl nb030nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel030nf_ia32_3dnow.nb030nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb030nf_n(%esp)
        jmp _nb_kernel030nf_ia32_3dnow.nb030nf_outer
_nb_kernel030nf_ia32_3dnow.nb030nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb030nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel030nf_ia32_3dnow.nb030nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel030nf_ia32_3dnow.nb030nf_threadloop
_nb_kernel030nf_ia32_3dnow.nb030nf_end: 
        femms

        movl nb030nf_nouter(%esp),%eax
        movl nb030nf_ninner(%esp),%ebx
        movl nb030nf_outeriter(%ebp),%ecx
        movl nb030nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        addl $96,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret

