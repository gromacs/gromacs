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



.globl nb_kernel330_ia32_3dnow
.globl _nb_kernel330_ia32_3dnow
nb_kernel330_ia32_3dnow:        
_nb_kernel330_ia32_3dnow:       
.set nb330_p_nri, 8
.set nb330_iinr, 12
.set nb330_jindex, 16
.set nb330_jjnr, 20
.set nb330_shift, 24
.set nb330_shiftvec, 28
.set nb330_fshift, 32
.set nb330_gid, 36
.set nb330_pos, 40
.set nb330_faction, 44
.set nb330_charge, 48
.set nb330_p_facel, 52
.set nb330_p_krf, 56
.set nb330_p_crf, 60
.set nb330_Vc, 64
.set nb330_type, 68
.set nb330_p_ntype, 72
.set nb330_vdwparam, 76
.set nb330_Vvdw, 80
.set nb330_p_tabscale, 84
.set nb330_VFtab, 88
.set nb330_invsqrta, 92
.set nb330_dvda, 96
.set nb330_p_gbtabscale, 100
.set nb330_GBtab, 104
.set nb330_p_nthreads, 108
.set nb330_count, 112
.set nb330_mtx, 116
.set nb330_outeriter, 120
.set nb330_inneriter, 124
.set nb330_work, 128
        ## stack offsets for local variables 
.set nb330_is3, 0
.set nb330_ii3, 4
.set nb330_ix, 8
.set nb330_iy, 12
.set nb330_iz, 16
.set nb330_iq, 20
.set nb330_vctot, 28
.set nb330_Vvdwtot, 36
.set nb330_c6, 44
.set nb330_c12, 52
.set nb330_two, 60
.set nb330_n1, 68
.set nb330_tsc, 76
.set nb330_ntia, 84
.set nb330_innerjjnr, 88
.set nb330_innerk, 92
.set nb330_fix, 96
.set nb330_fiy, 100
.set nb330_fiz, 104
.set nb330_dx1, 108
.set nb330_dy1, 112
.set nb330_dz1, 116
.set nb330_dx2, 120
.set nb330_dy2, 124
.set nb330_dz2, 128
.set nb330_n, 132                           ## idx for outer loop
.set nb330_nn1, 136                         ## number of outer iterations
.set nb330_nri, 140
.set nb330_facel, 144
.set nb330_ntype, 148
.set nb330_nouter, 152
.set nb330_ninner, 156
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $160,%esp          ## local stack space 
        femms
        ## move data to local stack  
        movl nb330_p_nri(%ebp),%ecx
        movl nb330_p_ntype(%ebp),%edx
        movl nb330_p_facel(%ebp),%esi
        movl nb330_p_tabscale(%ebp),%edi
        movl (%ecx),%ecx
        movl (%edx),%edx
        movl (%esi),%esi
        movl %ecx,nb330_nri(%esp)
        movl %edx,nb330_ntype(%esp)
        movl %esi,nb330_facel(%esp)

        movl $0x40000000,%eax
        movl %eax,nb330_two(%esp)
        movl %eax,nb330_two+4(%esp)
        movd  (%edi),%mm3
        punpckldq %mm3,%mm3
        movq  %mm3,nb330_tsc(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb330_nouter(%esp)
        movl %eax,nb330_ninner(%esp)

_nb_kernel330_ia32_3dnow.nb330_threadloop: 
        movl  nb330_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel330_ia32_3dnow.nb330_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel330_ia32_3dnow.nb330_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb330_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb330_n(%esp)
        movl %ebx,nb330_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel330_ia32_3dnow.nb330_outerstart
        jmp _nb_kernel330_ia32_3dnow.nb330_end

_nb_kernel330_ia32_3dnow.nb330_outerstart: 
        ## ebx contains number of outer iterations
        addl nb330_nouter(%esp),%ebx
        movl %ebx,nb330_nouter(%esp)

_nb_kernel330_ia32_3dnow.nb330_outer: 
        movl  nb330_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb330_is3(%esp)      ## store is3 

        movl  nb330_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movq  (%eax,%ebx,4),%mm0        ## move shX/shY to mm0 and shZ to mm1 
        movd  8(%eax,%ebx,4),%mm1

        movl  nb330_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx=ii 

        movl  nb330_charge(%ebp),%edx
        movd  (%edx,%ebx,4),%mm2    ## mm2=charge[ii] 
        pfmul nb330_facel(%esp),%mm2
        punpckldq %mm2,%mm2         ## spread to both halves 
        movq  %mm2,nb330_iq(%esp)           ## iq =facel*charge[ii] 

        movl  nb330_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb330_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb330_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb330_pos(%ebp),%eax      ## eax = base of pos[] 

        pfadd (%eax,%ebx,4),%mm0    ## ix = shX + posX (and iy too) 
        movd  8(%eax,%ebx,4),%mm3       ## cant use direct memory add for 4 bytes (iz) 
        movl  %ebx,nb330_ii3(%esp)
        pfadd %mm3,%mm1
        movq  %mm0,nb330_ix(%esp)
        movd  %mm1,nb330_iz(%esp)

        ## clear total potential and i forces 
        pxor  %mm7,%mm7
        movq  %mm7,nb330_vctot(%esp)
        movq  %mm7,nb330_Vvdwtot(%esp)
        movq  %mm7,nb330_fix(%esp)
        movd  %mm7,nb330_fiz(%esp)

        movl  nb330_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb330_pos(%ebp),%esi
        movl  nb330_faction(%ebp),%edi
        movl  nb330_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb330_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb330_ninner(%esp),%ecx
        movl  %ecx,nb330_ninner(%esp)
        movl  %edx,nb330_innerk(%esp)      ## number of innerloop atoms 
        addl  $0,%edx
        jge   _nb_kernel330_ia32_3dnow.nb330_unroll_loop
        jmp   _nb_kernel330_ia32_3dnow.nb330_finish_inner
_nb_kernel330_ia32_3dnow.nb330_unroll_loop: 
        ## paired innerloop starts here 
        movl  nb330_innerjjnr(%esp),%ecx       ## pointer to jjnr[k] 
        movl  (%ecx),%eax
        movl  4(%ecx),%ebx           ## eax/ebx=jnr 
        addl $8,nb330_innerjjnr(%esp)             ## advance pointer (unrolled 2) 
        prefetch 16(%ecx)            ## prefetch data - trial and error says 16 is best 

        movl nb330_charge(%ebp),%ecx     ## base of charge[] 
        movq nb330_iq(%esp),%mm5
        movd (%ecx,%eax,4),%mm3      ## charge[jnr1] 
        punpckldq (%ecx,%ebx,4),%mm3     ## move charge 2 to high part of mm3 
        pfmul %mm5,%mm3              ## mm3 now has qq for both particles 

        movl nb330_type(%ebp),%ecx
        movl (%ecx,%eax,4),%edx          ## type [jnr1] 
        movl (%ecx,%ebx,4),%ecx      ## type [jnr2] 

        movl nb330_vdwparam(%ebp),%esi          ## base of vdwparam  
        shll %edx
        shll %ecx
        addl nb330_ntia(%esp),%edx           ## tja = ntia + 2*type 
        addl nb330_ntia(%esp),%ecx

        movq (%esi,%edx,4),%mm5         ## mm5 = 1st c6 / c12           
        movq (%esi,%ecx,4),%mm7         ## mm7 = 2nd c6 / c12   
        movq %mm5,%mm6
        punpckldq %mm7,%mm5             ## mm5 = 1st c6 / 2nd c6 
        punpckhdq %mm7,%mm6             ## mm6 = 1st c12 / 2nd c12 
        movq %mm5,nb330_c6(%esp)
        movq %mm6,nb330_c12(%esp)

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        movl  nb330_pos(%ebp),%esi

        movq  nb330_ix(%esp),%mm0
        movd  nb330_iz(%esp),%mm1
        movq  (%esi,%eax,4),%mm4     ## fetch first j coordinates 
        movd  8(%esi,%eax,4),%mm5
        pfsubr %mm0,%mm4             ## dr = ir - jr  
        pfsubr %mm1,%mm5
        movq  %mm4,nb330_dx1(%esp)           ## store dr 
        movd  %mm5,nb330_dz1(%esp)
        pfmul %mm4,%mm4              ## square dx,dy,dz                          
        pfmul %mm5,%mm5
        pfacc %mm5,%mm4              ## accumulate to get dx*dx+ dy*dy+ dz*dz 
        pfacc %mm5,%mm4              ## first rsq in lower mm4 

        movq  (%esi,%ebx,4),%mm6     ## fetch second j coordinates  
        movd  8(%esi,%ebx,4),%mm7

        pfsubr %mm0,%mm6             ## dr = ir - jr  
        pfsubr %mm1,%mm7
        movq  %mm6,nb330_dx2(%esp)           ## store dr 
        movd  %mm7,nb330_dz2(%esp)
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
        pfmul nb330_tsc(%esp),%mm1      ## mm1=rt 
        pf2iw %mm1,%mm4
        movq %mm4,nb330_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm1              ## now mm1 is eps and mm4 is n0 

        movq %mm1,%mm2
        pfmul %mm2,%mm2 ## mm1 is eps, mm2 is eps2 

        movl nb330_VFtab(%ebp),%edx
        movl nb330_n1(%esp),%ecx
        leal (%ecx,%ecx,2),%ecx
        shll $2,%ecx
        ## load all the table values we need 
        movd (%edx,%ecx,4),%mm4
        movd 4(%edx,%ecx,4),%mm5
        movd 8(%edx,%ecx,4),%mm6
        movd 12(%edx,%ecx,4),%mm7
        movl nb330_n1+4(%esp),%ecx
        leal (%ecx,%ecx,2),%ecx
        shll $2,%ecx
        punpckldq (%edx,%ecx,4),%mm4
        punpckldq 4(%edx,%ecx,4),%mm5
        punpckldq 8(%edx,%ecx,4),%mm6
        punpckldq 12(%edx,%ecx,4),%mm7

        pfmul %mm1,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 

        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 

        pfmul nb330_two(%esp),%mm7      ## two*Heps2 
        pfadd %mm6,%mm7
        pfadd %mm5,%mm7 ## mm7=FF 

        pfmul %mm1,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        pfmul %mm3,%mm5 ## vcoul=qq*VV 
        pfmul %mm7,%mm3 ## fijC=FF*qq  

        ## at this point mm5 contains vcoul and mm3 fijC 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        pfadd nb330_vctot(%esp),%mm5        ## add the earlier value 
        movq %mm5,nb330_vctot(%esp)         ## store the sum       

        ## dispersion table 
        movl nb330_n1(%esp),%ecx
        leal (%ecx,%ecx,2),%ecx
        shll $2,%ecx
        ## load all the table values we need 
        movd 16(%edx,%ecx,4),%mm4
        movd 20(%edx,%ecx,4),%mm5
        movd 24(%edx,%ecx,4),%mm6
        movd 28(%edx,%ecx,4),%mm7
        movl nb330_n1+4(%esp),%ecx
        leal (%ecx,%ecx,2),%ecx
        shll $2,%ecx
        punpckldq 16(%edx,%ecx,4),%mm4
        punpckldq 20(%edx,%ecx,4),%mm5
        punpckldq 24(%edx,%ecx,4),%mm6
        punpckldq 28(%edx,%ecx,4),%mm7
        pfmul %mm1,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 
        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 
        pfmul nb330_two(%esp),%mm7      ## two*Heps2 
        pfadd %mm6,%mm7
        pfadd %mm5,%mm7 ## mm7=FF 
        pfmul %mm1,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV     

        movq nb330_c6(%esp),%mm4
        pfmul %mm4,%mm7 ## fijD 
        pfmul %mm4,%mm5 ## Vvdw6            
        pfadd %mm7,%mm3 ## add to fscal 

        ## update Vvdwtot to release mm5! 
        pfadd nb330_Vvdwtot(%esp),%mm5        ## add the earlier value 
        movq %mm5,nb330_Vvdwtot(%esp)         ## store the sum       

        ## repulsion table 
        movl nb330_n1(%esp),%ecx
        leal (%ecx,%ecx,2),%ecx
        shll $2,%ecx
        ## load all the table values we need 
        movd 32(%edx,%ecx,4),%mm4
        movd 36(%edx,%ecx,4),%mm5
        movd 40(%edx,%ecx,4),%mm6
        movd 44(%edx,%ecx,4),%mm7
        movl nb330_n1+4(%esp),%ecx
        leal (%ecx,%ecx,2),%ecx
        shll $2,%ecx
        punpckldq 32(%edx,%ecx,4),%mm4
        punpckldq 36(%edx,%ecx,4),%mm5
        punpckldq 40(%edx,%ecx,4),%mm6
        punpckldq 44(%edx,%ecx,4),%mm7

        pfmul %mm1,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 
        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 
        pfmul nb330_two(%esp),%mm7      ## two*Heps2 
        pfadd %mm6,%mm7
        pfadd %mm5,%mm7 ## mm7=FF 
        pfmul %mm1,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        movq nb330_c12(%esp),%mm6
        pfmul %mm6,%mm7 ## fijR 
        pfmul %mm6,%mm5 ## Vvdw12 
        pfadd %mm7,%mm3 ## total fscal fijC+ fijD+ fijR 

        ## change sign of mm3 
        pxor %mm1,%mm1
        pfsub %mm3,%mm1
        pfmul nb330_tsc(%esp),%mm0
        pfmul %mm1,%mm0   ## mm0 is total fscal now     

        prefetchw nb330_dx1(%esp)       ## prefetch i forces to cache 

        ## spread fscalar to both positions 
        movq %mm0,%mm1
        punpckldq %mm0,%mm0
        punpckhdq %mm1,%mm1

        ## calc vector force 
        prefetchw (%edi,%eax,4) ## prefetch the 1st faction to cache 
        movq nb330_dx1(%esp),%mm2       ## fetch dr 
        movd nb330_dz1(%esp),%mm3

        ## update Vvdwtot 
        pfadd nb330_Vvdwtot(%esp),%mm5        ## add the earlier value 
        movq %mm5,nb330_Vvdwtot(%esp)         ## store the sum       

        prefetchw (%edi,%ebx,4) ## prefetch the 2nd faction to cache 
        pfmul %mm0,%mm2         ## mult by fs  
        pfmul %mm0,%mm3

        movq nb330_dx2(%esp),%mm4       ## fetch dr 
        movd nb330_dz2(%esp),%mm5
        pfmul %mm1,%mm4         ## mult by fs  
        pfmul %mm1,%mm5
        ## update i forces 

        movq nb330_fix(%esp),%mm0
        movd nb330_fiz(%esp),%mm1
        pfadd %mm2,%mm0
        pfadd %mm3,%mm1

        pfadd %mm4,%mm0
        pfadd %mm5,%mm1
        movq %mm0,nb330_fix(%esp)
        movd %mm1,nb330_fiz(%esp)
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
        subl $2,nb330_innerk(%esp)
        jl    _nb_kernel330_ia32_3dnow.nb330_finish_inner
        jmp   _nb_kernel330_ia32_3dnow.nb330_unroll_loop
_nb_kernel330_ia32_3dnow.nb330_finish_inner: 
        andl $1,nb330_innerk(%esp)
        jnz  _nb_kernel330_ia32_3dnow.nb330_single_inner
        jmp  _nb_kernel330_ia32_3dnow.nb330_updateouterdata
_nb_kernel330_ia32_3dnow.nb330_single_inner: 
        ## a single j particle iteration here - compare with the unrolled code for comments. 
        movl  nb330_innerjjnr(%esp),%eax
        movl  (%eax),%eax       ## eax=jnr offset 

        movl nb330_charge(%ebp),%ecx
        movd nb330_iq(%esp),%mm5
        movd (%ecx,%eax,4),%mm3
        pfmul %mm5,%mm3         ## mm3=qq 

        movl nb330_vdwparam(%ebp),%esi
        movl nb330_type(%ebp),%ecx
        movl (%ecx,%eax,4),%edx          ## type [jnr1] 
        shll %edx
        addl nb330_ntia(%esp),%edx           ## tja = ntia + 2*type 
        movd (%esi,%edx,4),%mm5         ## mm5 = 1st c6                 
        movq %mm5,nb330_c6(%esp)
        movd 4(%esi,%edx,4),%mm5        ## mm5 = 1st c12                
        movq %mm5,nb330_c12(%esp)

        movl  nb330_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        movq  nb330_ix(%esp),%mm0
        movd  nb330_iz(%esp),%mm1
        movq  (%esi,%eax,4),%mm4
        movd  8(%esi,%eax,4),%mm5
        pfsubr %mm0,%mm4
        pfsubr %mm1,%mm5
        movq  %mm4,nb330_dx1(%esp)
        pfmul %mm4,%mm4
        movd  %mm5,nb330_dz1(%esp)
        pfmul %mm5,%mm5
        pfacc %mm5,%mm4
        pfacc %mm5,%mm4         ## mm0=rsq 

        pfrsqrt %mm4,%mm0
        movq %mm0,%mm2
        pfmul %mm0,%mm0
        pfrsqit1 %mm4,%mm0
        pfrcpit2 %mm2,%mm0      ## mm0=invsqrt 
        pfmul %mm0,%mm4
        movq %mm4,%mm1
        ## mm0 is invsqrt, and mm1 r. 

        ## calculate potentials and scalar force 
        pfmul nb330_tsc(%esp),%mm1      ## mm1=rt 
        pf2iw %mm1,%mm4
        movd %mm4,nb330_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm1              ## now mm1 is eps and mm4 is n0 

        movq %mm1,%mm2
        pfmul %mm2,%mm2 ## mm1 is eps, mm2 is eps2 

        ## coulomb table 
        movl nb330_VFtab(%ebp),%edx
        movl nb330_n1(%esp),%ecx
        leal (%ecx,%ecx,2),%ecx
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

        pfmul nb330_two(%esp),%mm7      ## two*Heps2 
        pfadd %mm6,%mm7
        pfadd %mm5,%mm7 ## mm7=FF 

        pfmul %mm1,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        pfmul %mm3,%mm5 ## vcoul=qq*VV 
        pfmul %mm7,%mm3 ## fijC=FF*qq  

        ## at this point mm5 contains vcoul and mm3 fijC 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        pfadd nb330_vctot(%esp),%mm5        ## add the earlier value 
        movq %mm5,nb330_vctot(%esp)         ## store the sum       

        ## dispersion table 
        ## load all the table values we need 
        movd 16(%edx,%ecx,4),%mm4
        movd 20(%edx,%ecx,4),%mm5
        movd 24(%edx,%ecx,4),%mm6
        movd 28(%edx,%ecx,4),%mm7
        pfmul %mm1,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 
        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 
        pfmul nb330_two(%esp),%mm7      ## two*Heps2 
        pfadd %mm6,%mm7
        pfadd %mm5,%mm7 ## mm7=FF 
        pfmul %mm1,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV     

        movq nb330_c6(%esp),%mm4
        pfmul %mm4,%mm7 ## fijD 
        pfmul %mm4,%mm5 ## Vvdw6            
        pfadd %mm7,%mm3 ## add to fscal 

        ## update Vvdwtot to release mm5! 
        pfadd nb330_Vvdwtot(%esp),%mm5        ## add the earlier value 
        movq %mm5,nb330_Vvdwtot(%esp)         ## store the sum       

        ## repulsion table 
        ## load all the table values we need 
        movd 32(%edx,%ecx,4),%mm4
        movd 36(%edx,%ecx,4),%mm5
        movd 40(%edx,%ecx,4),%mm6
        movd 44(%edx,%ecx,4),%mm7

        pfmul %mm1,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 
        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 
        pfmul nb330_two(%esp),%mm7      ## two*Heps2 
        pfadd %mm6,%mm7
        pfadd %mm5,%mm7 ## mm7=FF 
        pfmul %mm1,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        movq nb330_c12(%esp),%mm6
        pfmul %mm6,%mm7 ## fijR 
        pfmul %mm6,%mm5 ## Vvdw12 
        pfadd %mm7,%mm3 ## total fscal fijC+ fijD+ fijR 

        ## change sign of mm3 
        pxor %mm1,%mm1
        pfsub %mm3,%mm1
        pfmul nb330_tsc(%esp),%mm0
        pfmul %mm1,%mm0   ## mm0 is total fscal now     

        ## update Vvdwtot 
        pfadd nb330_Vvdwtot(%esp),%mm5        ## add the earlier value 
        movq %mm5,nb330_Vvdwtot(%esp)         ## store the sum       

        ## spread fscalar to both positions 
        punpckldq %mm0,%mm0
        ## calc vectorial force 
        prefetchw (%edi,%eax,4) ## prefetch faction to cache  
        movq nb330_dx1(%esp),%mm2
        movd nb330_dz1(%esp),%mm3


        pfmul %mm0,%mm2
        pfmul %mm0,%mm3

        ## update i particle force 
        movq nb330_fix(%esp),%mm0
        movd nb330_fiz(%esp),%mm1
        pfadd %mm2,%mm0
        pfadd %mm3,%mm1
        movq %mm0,nb330_fix(%esp)
        movd %mm1,nb330_fiz(%esp)
        ## update j particle force 
        movq (%edi,%eax,4),%mm0
        movd 8(%edi,%eax,4),%mm1
        pfsub %mm2,%mm0
        pfsub %mm3,%mm1
        movq %mm0,(%edi,%eax,4)
        movd %mm1,8(%edi,%eax,4)
        ## done! 
_nb_kernel330_ia32_3dnow.nb330_updateouterdata: 
        movl  nb330_ii3(%esp),%ecx

        movq  (%edi,%ecx,4),%mm6       ## increment i force 
        movd  8(%edi,%ecx,4),%mm7
        pfadd nb330_fix(%esp),%mm6
        pfadd nb330_fiz(%esp),%mm7
        movq  %mm6,(%edi,%ecx,4)
        movd  %mm7,8(%edi,%ecx,4)

        movl  nb330_fshift(%ebp),%ebx      ## increment fshift force 
        movl  nb330_is3(%esp),%edx

        movq  (%ebx,%edx,4),%mm6
        movd  8(%ebx,%edx,4),%mm7
        pfadd nb330_fix(%esp),%mm6
        pfadd nb330_fiz(%esp),%mm7
        movq  %mm6,(%ebx,%edx,4)
        movd  %mm7,8(%ebx,%edx,4)

        ## get n from stack
        movl nb330_n(%esp),%esi
        ## get group index for i particle 
        movl  nb330_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        movq  nb330_vctot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb330_Vc(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment vc[gid] 

        movq  nb330_Vvdwtot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb330_Vvdw(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment Vvdw[gid] 

        ## finish if last 
        movl nb330_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel330_ia32_3dnow.nb330_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb330_n(%esp)
        jmp _nb_kernel330_ia32_3dnow.nb330_outer
_nb_kernel330_ia32_3dnow.nb330_outerend: 
        ## check if more outer neighborlists remain
        movl  nb330_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel330_ia32_3dnow.nb330_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel330_ia32_3dnow.nb330_threadloop
_nb_kernel330_ia32_3dnow.nb330_end: 
        femms

        movl nb330_nouter(%esp),%eax
        movl nb330_ninner(%esp),%ebx
        movl nb330_outeriter(%ebp),%ecx
        movl nb330_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        addl $160,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret



.globl nb_kernel330nf_ia32_3dnow
.globl _nb_kernel330nf_ia32_3dnow
nb_kernel330nf_ia32_3dnow:      
_nb_kernel330nf_ia32_3dnow:     
.set nb330nf_p_nri, 8
.set nb330nf_iinr, 12
.set nb330nf_jindex, 16
.set nb330nf_jjnr, 20
.set nb330nf_shift, 24
.set nb330nf_shiftvec, 28
.set nb330nf_fshift, 32
.set nb330nf_gid, 36
.set nb330nf_pos, 40
.set nb330nf_faction, 44
.set nb330nf_charge, 48
.set nb330nf_p_facel, 52
.set nb330nf_p_krf, 56
.set nb330nf_p_crf, 60
.set nb330nf_Vc, 64
.set nb330nf_type, 68
.set nb330nf_p_ntype, 72
.set nb330nf_vdwparam, 76
.set nb330nf_Vvdw, 80
.set nb330nf_p_tabscale, 84
.set nb330nf_VFtab, 88
.set nb330nf_invsqrta, 92
.set nb330nf_dvda, 96
.set nb330nf_p_gbtabscale, 100
.set nb330nf_GBtab, 104
.set nb330nf_p_nthreads, 108
.set nb330nf_count, 112
.set nb330nf_mtx, 116
.set nb330nf_outeriter, 120
.set nb330nf_inneriter, 124
.set nb330nf_work, 128
        ## stack offsets for local variables 
.set nb330nf_is3, 0
.set nb330nf_ii3, 4
.set nb330nf_ix, 8
.set nb330nf_iy, 12
.set nb330nf_iz, 16
.set nb330nf_iq, 20
.set nb330nf_vctot, 28
.set nb330nf_Vvdwtot, 36
.set nb330nf_c6, 44
.set nb330nf_c12, 52
.set nb330nf_n1, 60
.set nb330nf_tsc, 68
.set nb330nf_ntia, 76
.set nb330nf_innerjjnr, 80
.set nb330nf_innerk, 84
.set nb330nf_n, 88                         ## idx for outer loop
.set nb330nf_nn1, 92                       ## number of outer iterations
.set nb330nf_nri, 96
.set nb330nf_facel, 100
.set nb330nf_ntype, 104
.set nb330nf_nouter, 108
.set nb330nf_ninner, 112
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
        movl nb330nf_p_nri(%ebp),%ecx
        movl nb330nf_p_ntype(%ebp),%edx
        movl nb330nf_p_facel(%ebp),%esi
        movl nb330nf_p_tabscale(%ebp),%edi
        movl (%ecx),%ecx
        movl (%edx),%edx
        movl (%esi),%esi
        movl %ecx,nb330nf_nri(%esp)
        movl %edx,nb330nf_ntype(%esp)
        movl %esi,nb330nf_facel(%esp)

        movd  (%edi),%mm3
        punpckldq %mm3,%mm3
        movq  %mm3,nb330nf_tsc(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb330nf_nouter(%esp)
        movl %eax,nb330nf_ninner(%esp)

        ## assume we have at least one i particle - start directly      
_nb_kernel330nf_ia32_3dnow.nb330nf_threadloop: 
        movl  nb330nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel330nf_ia32_3dnow.nb330nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel330nf_ia32_3dnow.nb330nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb330nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb330nf_n(%esp)
        movl %ebx,nb330nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi

        jg  _nb_kernel330nf_ia32_3dnow.nb330nf_outerstart
        jmp _nb_kernel330nf_ia32_3dnow.nb330nf_end

_nb_kernel330nf_ia32_3dnow.nb330nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb330nf_nouter(%esp),%ebx
        movl %ebx,nb330nf_nouter(%esp)

_nb_kernel330nf_ia32_3dnow.nb330nf_outer: 
        movl  nb330nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb330nf_is3(%esp)            ## store is3 

        movl  nb330nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movq  (%eax,%ebx,4),%mm0        ## move shX/shY to mm0 and shZ to mm1 
        movd  8(%eax,%ebx,4),%mm1

        movl  nb330nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx=ii 

        movl  nb330nf_charge(%ebp),%edx
        movd  (%edx,%ebx,4),%mm2    ## mm2=charge[ii] 
        pfmul nb330nf_facel(%esp),%mm2
        punpckldq %mm2,%mm2         ## spread to both halves 
        movq  %mm2,nb330nf_iq(%esp)         ## iq =facel*charge[ii] 

        movl  nb330nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb330nf_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb330nf_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb330nf_pos(%ebp),%eax      ## eax = base of pos[] 

        pfadd (%eax,%ebx,4),%mm0    ## ix = shX + posX (and iy too) 
        movd  8(%eax,%ebx,4),%mm3       ## cant use direct memory add for 4 bytes (iz) 
        movl  %ebx,nb330nf_ii3(%esp)
        pfadd %mm3,%mm1
        movq  %mm0,nb330nf_ix(%esp)
        movd  %mm1,nb330nf_iz(%esp)

        ## clear total potential and i forces 
        pxor  %mm7,%mm7
        movq  %mm7,nb330nf_vctot(%esp)
        movq  %mm7,nb330nf_Vvdwtot(%esp)

        movl  nb330nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb330nf_pos(%ebp),%esi
        movl  nb330nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb330nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb330nf_ninner(%esp),%ecx
        movl  %ecx,nb330nf_ninner(%esp)
        movl  %edx,nb330nf_innerk(%esp)      ## number of innerloop atoms 
        addl  $0,%edx
        jge   _nb_kernel330nf_ia32_3dnow.nb330nf_unroll_loop
        jmp   _nb_kernel330nf_ia32_3dnow.nb330nf_finish_inner
_nb_kernel330nf_ia32_3dnow.nb330nf_unroll_loop: 
        ## paired innerloop starts here 
        movl  nb330nf_innerjjnr(%esp),%ecx       ## pointer to jjnr[k] 
        movl  (%ecx),%eax
        movl  4(%ecx),%ebx           ## eax/ebx=jnr 
        addl $8,nb330nf_innerjjnr(%esp)             ## advance pointer (unrolled 2) 
        prefetch 16(%ecx)            ## prefetch data - trial and error says 16 is best 

        movl nb330nf_charge(%ebp),%ecx     ## base of charge[] 
        movq nb330nf_iq(%esp),%mm5
        movd (%ecx,%eax,4),%mm3      ## charge[jnr1] 
        punpckldq (%ecx,%ebx,4),%mm3     ## move charge 2 to high part of mm3 
        pfmul %mm5,%mm3              ## mm3 now has qq for both particles 

        movl nb330nf_type(%ebp),%ecx
        movl (%ecx,%eax,4),%edx          ## type [jnr1] 
        movl (%ecx,%ebx,4),%ecx      ## type [jnr2] 

        movl nb330nf_vdwparam(%ebp),%esi                ## base of vdwparam  
        shll %edx
        shll %ecx
        addl nb330nf_ntia(%esp),%edx         ## tja = ntia + 2*type 
        addl nb330nf_ntia(%esp),%ecx

        movq (%esi,%edx,4),%mm5         ## mm5 = 1st c6 / c12           
        movq (%esi,%ecx,4),%mm7         ## mm7 = 2nd c6 / c12   
        movq %mm5,%mm6
        punpckldq %mm7,%mm5             ## mm5 = 1st c6 / 2nd c6 
        punpckhdq %mm7,%mm6             ## mm6 = 1st c12 / 2nd c12 
        movq %mm5,nb330nf_c6(%esp)
        movq %mm6,nb330nf_c12(%esp)

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        movl  nb330nf_pos(%ebp),%esi

        movq  nb330nf_ix(%esp),%mm0
        movd  nb330nf_iz(%esp),%mm1
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
        pfmul nb330nf_tsc(%esp),%mm1    ## mm1=rt 
        pf2iw %mm1,%mm4
        movq %mm4,nb330nf_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm1              ## now mm1 is eps and mm4 is n0 

        movq %mm1,%mm2
        pfmul %mm2,%mm2 ## mm1 is eps, mm2 is eps2 

        movl nb330nf_VFtab(%ebp),%edx
        movl nb330nf_n1(%esp),%ecx
        leal (%ecx,%ecx,2),%ecx
        shll $2,%ecx
        ## load all the table values we need 
        movd (%edx,%ecx,4),%mm4
        movd 4(%edx,%ecx,4),%mm5
        movd 8(%edx,%ecx,4),%mm6
        movd 12(%edx,%ecx,4),%mm7
        movl nb330nf_n1+4(%esp),%ecx
        leal (%ecx,%ecx,2),%ecx
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
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        pfadd nb330nf_vctot(%esp),%mm5        ## add the earlier value 
        movq %mm5,nb330nf_vctot(%esp)         ## store the sum       

        ## dispersion table 
        movl nb330nf_n1(%esp),%ecx
        leal (%ecx,%ecx,2),%ecx
        shll $2,%ecx
        ## load all the table values we need 
        movd 16(%edx,%ecx,4),%mm4
        movd 20(%edx,%ecx,4),%mm5
        movd 24(%edx,%ecx,4),%mm6
        movd 28(%edx,%ecx,4),%mm7
        movl nb330nf_n1+4(%esp),%ecx
        leal (%ecx,%ecx,2),%ecx
        shll $2,%ecx
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

        movq nb330nf_c6(%esp),%mm4
        pfmul %mm4,%mm5 ## Vvdw6    
        ## update Vvdwtot to release mm5! 
        pfadd nb330nf_Vvdwtot(%esp),%mm5        ## add the earlier value 
        movq %mm5,nb330nf_Vvdwtot(%esp)         ## store the sum       

        ## repulsion table 
        movl nb330nf_n1(%esp),%ecx
        leal (%ecx,%ecx,2),%ecx
        shll $2,%ecx
        ## load all the table values we need 
        movd 32(%edx,%ecx,4),%mm4
        movd 36(%edx,%ecx,4),%mm5
        movd 40(%edx,%ecx,4),%mm6
        movd 44(%edx,%ecx,4),%mm7
        movl nb330nf_n1+4(%esp),%ecx
        leal (%ecx,%ecx,2),%ecx
        shll $2,%ecx
        punpckldq 32(%edx,%ecx,4),%mm4
        punpckldq 36(%edx,%ecx,4),%mm5
        punpckldq 40(%edx,%ecx,4),%mm6
        punpckldq 44(%edx,%ecx,4),%mm7

        pfmul %mm1,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 
        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 
        pfmul %mm1,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        movq nb330nf_c12(%esp),%mm6
        pfmul %mm6,%mm5 ## Vvdw12 
        ## update Vvdwtot 
        pfadd nb330nf_Vvdwtot(%esp),%mm5        ## add the earlier value 
        movq %mm5,nb330nf_Vvdwtot(%esp)         ## store the sum       

        ## should we do one more iteration? 
        subl $2,nb330nf_innerk(%esp)
        jl    _nb_kernel330nf_ia32_3dnow.nb330nf_finish_inner
        jmp   _nb_kernel330nf_ia32_3dnow.nb330nf_unroll_loop
_nb_kernel330nf_ia32_3dnow.nb330nf_finish_inner: 
        andl $1,nb330nf_innerk(%esp)
        jnz  _nb_kernel330nf_ia32_3dnow.nb330nf_single_inner
        jmp  _nb_kernel330nf_ia32_3dnow.nb330nf_updateouterdata
_nb_kernel330nf_ia32_3dnow.nb330nf_single_inner: 
        ## a single j particle iteration here - compare with the unrolled code for comments. 
        movl  nb330nf_innerjjnr(%esp),%eax
        movl  (%eax),%eax       ## eax=jnr offset 

        movl nb330nf_charge(%ebp),%ecx
        movd nb330nf_iq(%esp),%mm5
        movd (%ecx,%eax,4),%mm3
        pfmul %mm5,%mm3         ## mm3=qq 

        movl nb330nf_vdwparam(%ebp),%esi
        movl nb330nf_type(%ebp),%ecx
        movl (%ecx,%eax,4),%edx          ## type [jnr1] 
        shll %edx
        addl nb330nf_ntia(%esp),%edx         ## tja = ntia + 2*type 
        movd (%esi,%edx,4),%mm5         ## mm5 = 1st c6                 
        movq %mm5,nb330nf_c6(%esp)
        movd 4(%esi,%edx,4),%mm5        ## mm5 = 1st c12                
        movq %mm5,nb330nf_c12(%esp)

        movl  nb330nf_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        movq  nb330nf_ix(%esp),%mm0
        movd  nb330nf_iz(%esp),%mm1
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
        pfmul %mm0,%mm4
        movq %mm4,%mm1
        ## mm0 is invsqrt, and mm1 r. 

        ## calculate potentials and scalar force 
        pfmul nb330nf_tsc(%esp),%mm1    ## mm1=rt 
        pf2iw %mm1,%mm4
        movd %mm4,nb330nf_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm1              ## now mm1 is eps and mm4 is n0 

        movq %mm1,%mm2
        pfmul %mm2,%mm2 ## mm1 is eps, mm2 is eps2 

        ## coulomb table 
        movl nb330nf_VFtab(%ebp),%edx
        movl nb330nf_n1(%esp),%ecx
        leal (%ecx,%ecx,2),%ecx
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
        pfadd nb330nf_vctot(%esp),%mm5        ## add the earlier value 
        movq %mm5,nb330nf_vctot(%esp)         ## store the sum       

        ## dispersion table 
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

        movq nb330nf_c6(%esp),%mm4
        pfmul %mm4,%mm5 ## Vvdw6  

        ## update Vvdwtot to release mm5! 
        pfadd nb330nf_Vvdwtot(%esp),%mm5        ## add the earlier value 
        movq %mm5,nb330nf_Vvdwtot(%esp)         ## store the sum       

        ## repulsion table 
        ## load all the table values we need 
        movd 32(%edx,%ecx,4),%mm4
        movd 36(%edx,%ecx,4),%mm5
        movd 40(%edx,%ecx,4),%mm6
        movd 44(%edx,%ecx,4),%mm7

        pfmul %mm1,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 
        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 
        pfmul %mm1,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        movq nb330nf_c12(%esp),%mm6
        pfmul %mm6,%mm5 ## Vvdw12 
        ## update Vvdwtot 
        pfadd nb330nf_Vvdwtot(%esp),%mm5        ## add the earlier value 
        movq %mm5,nb330nf_Vvdwtot(%esp)         ## store the sum       

_nb_kernel330nf_ia32_3dnow.nb330nf_updateouterdata: 
        ## get n from stack
        movl nb330nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb330nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        movq  nb330nf_vctot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb330nf_Vc(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment vc[gid] 

        movq  nb330nf_Vvdwtot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb330nf_Vvdw(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment Vvdw[gid] 

        ## finish if last 
        movl nb330nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel330nf_ia32_3dnow.nb330nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb330nf_n(%esp)
        jmp _nb_kernel330nf_ia32_3dnow.nb330nf_outer
_nb_kernel330nf_ia32_3dnow.nb330nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb330nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel330nf_ia32_3dnow.nb330nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel330nf_ia32_3dnow.nb330nf_threadloop
_nb_kernel330nf_ia32_3dnow.nb330nf_end: 
        femms

        movl nb330nf_nouter(%esp),%eax
        movl nb330nf_ninner(%esp),%ebx
        movl nb330nf_outeriter(%ebp),%ecx
        movl nb330nf_inneriter(%ebp),%edx
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



