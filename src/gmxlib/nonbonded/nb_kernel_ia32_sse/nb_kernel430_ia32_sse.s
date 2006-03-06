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




.globl nb_kernel430_ia32_sse
.globl _nb_kernel430_ia32_sse
nb_kernel430_ia32_sse:  
_nb_kernel430_ia32_sse: 
.set nb430_p_nri, 8
.set nb430_iinr, 12
.set nb430_jindex, 16
.set nb430_jjnr, 20
.set nb430_shift, 24
.set nb430_shiftvec, 28
.set nb430_fshift, 32
.set nb430_gid, 36
.set nb430_pos, 40
.set nb430_faction, 44
.set nb430_charge, 48
.set nb430_p_facel, 52
.set nb430_argkrf, 56
.set nb430_argcrf, 60
.set nb430_Vc, 64
.set nb430_type, 68
.set nb430_p_ntype, 72
.set nb430_vdwparam, 76
.set nb430_Vvdw, 80
.set nb430_p_tabscale, 84
.set nb430_VFtab, 88
.set nb430_invsqrta, 92
.set nb430_dvda, 96
.set nb430_p_gbtabscale, 100
.set nb430_GBtab, 104
.set nb430_p_nthreads, 108
.set nb430_count, 112
.set nb430_mtx, 116
.set nb430_outeriter, 120
.set nb430_inneriter, 124
.set nb430_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb430_ix, 0
.set nb430_iy, 16
.set nb430_iz, 32
.set nb430_iq, 48
.set nb430_dx, 64
.set nb430_dy, 80
.set nb430_dz, 96
.set nb430_two, 112
.set nb430_gbtsc, 128
.set nb430_tsc, 144
.set nb430_qq, 160
.set nb430_c6, 176
.set nb430_c12, 192
.set nb430_fscal, 208
.set nb430_vctot, 224
.set nb430_Vvdwtot, 240
.set nb430_fix, 256
.set nb430_fiy, 272
.set nb430_fiz, 288
.set nb430_half, 304
.set nb430_three, 320
.set nb430_r, 336
.set nb430_isai, 352
.set nb430_isaprod, 368
.set nb430_dvdasum, 384
.set nb430_gbscale, 400
.set nb430_ii, 416
.set nb430_is3, 420
.set nb430_ii3, 424
.set nb430_ntia, 428
.set nb430_innerjjnr, 432
.set nb430_innerk, 436
.set nb430_n, 440
.set nb430_nn1, 444
.set nb430_jnra, 448
.set nb430_jnrb, 452
.set nb430_jnrc, 456
.set nb430_jnrd, 460
.set nb430_nri, 464
.set nb430_facel, 468
.set nb430_ntype, 472
.set nb430_nouter, 476
.set nb430_ninner, 480
.set nb430_salign, 484
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $488,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb430_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb430_p_nri(%ebp),%ecx
        movl nb430_p_facel(%ebp),%esi
        movl nb430_p_ntype(%ebp),%edi
        movl (%ecx),%ecx
        movl (%esi),%esi
        movl (%edi),%edi
        movl %ecx,nb430_nri(%esp)
        movl %esi,nb430_facel(%esp)
        movl %edi,nb430_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb430_nouter(%esp)
        movl %eax,nb430_ninner(%esp)


        movl nb430_p_gbtabscale(%ebp),%eax
        movss (%eax),%xmm3
        movl nb430_p_tabscale(%ebp),%eax
        movss (%eax),%xmm4
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb430_gbtsc(%esp)
        movaps %xmm4,nb430_tsc(%esp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb430_half(%esp)
        movss nb430_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb430_half(%esp)
        movaps %xmm2,nb430_two(%esp)
        movaps %xmm3,nb430_three(%esp)

_nb_kernel430_ia32_sse.nb430_threadloop: 
        movl  nb430_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel430_ia32_sse.nb430_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel430_ia32_sse.nb430_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb430_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb430_n(%esp)
        movl %ebx,nb430_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel430_ia32_sse.nb430_outerstart
        jmp _nb_kernel430_ia32_sse.nb430_end

_nb_kernel430_ia32_sse.nb430_outerstart: 
        ## ebx contains number of outer iterations
        addl nb430_nouter(%esp),%ebx
        movl %ebx,nb430_nouter(%esp)

_nb_kernel430_ia32_sse.nb430_outer: 
        movl  nb430_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb430_is3(%esp)      ## store is3 

        movl  nb430_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb430_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 
        movl  %ebx,nb430_ii(%esp)

        movl  nb430_charge(%ebp),%edx
        movss (%edx,%ebx,4),%xmm3
        mulss nb430_facel(%esp),%xmm3
        shufps $0,%xmm3,%xmm3

        movl  nb430_invsqrta(%ebp),%edx         ## load invsqrta[ii]
        movss (%edx,%ebx,4),%xmm4
        shufps $0,%xmm4,%xmm4

        movl  nb430_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb430_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb430_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb430_pos(%ebp),%eax      ## eax = base of pos[]  

        addss (%eax,%ebx,4),%xmm0
        addss 4(%eax,%ebx,4),%xmm1
        addss 8(%eax,%ebx,4),%xmm2

        movaps %xmm3,nb430_iq(%esp)
        movaps %xmm4,nb430_isai(%esp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb430_ix(%esp)
        movaps %xmm1,nb430_iy(%esp)
        movaps %xmm2,nb430_iz(%esp)

        movl  %ebx,nb430_ii3(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb430_vctot(%esp)
        movaps %xmm4,nb430_Vvdwtot(%esp)
        movaps %xmm4,nb430_dvdasum(%esp)
        movaps %xmm4,nb430_fix(%esp)
        movaps %xmm4,nb430_fiy(%esp)
        movaps %xmm4,nb430_fiz(%esp)

        movl  nb430_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb430_pos(%ebp),%esi
        movl  nb430_faction(%ebp),%edi
        movl  nb430_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb430_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb430_ninner(%esp),%ecx
        movl  %ecx,nb430_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb430_innerk(%esp)      ## number of innerloop atoms

        jge   _nb_kernel430_ia32_sse.nb430_unroll_loop
        jmp   _nb_kernel430_ia32_sse.nb430_finish_inner
_nb_kernel430_ia32_sse.nb430_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb430_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 
        addl $16,nb430_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        ## load isaj
        movl nb430_invsqrta(%ebp),%esi
        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7
        movaps nb430_isai(%esp),%xmm2
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all isaj in xmm3  
        mulps  %xmm3,%xmm2

        movaps %xmm2,nb430_isaprod(%esp)
        movaps %xmm2,%xmm1
        mulps nb430_gbtsc(%esp),%xmm1
        movaps %xmm1,nb430_gbscale(%esp)

        movl nb430_charge(%ebp),%esi     ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        mulps nb430_iq(%esp),%xmm2
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        mulps  %xmm2,%xmm3
        movaps %xmm3,nb430_qq(%esp)

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movl nb430_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb430_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb430_ntia(%esp),%edi
        addl %edi,%eax
        addl %edi,%ebx
        addl %edi,%ecx
        addl %edi,%edx

        movlps (%esi,%eax,4),%xmm6
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm6
        movhps (%esi,%edx,4),%xmm7

        movaps %xmm6,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm6 ## constant 11011101

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movd  %mm2,%ecx
        movd  %mm3,%edx

        movaps %xmm4,nb430_c6(%esp)
        movaps %xmm6,nb430_c12(%esp)

        movl nb430_pos(%ebp),%esi        ## base of pos[] 

        movl %eax,nb430_jnra(%esp)
        movl %ebx,nb430_jnrb(%esp)
        movl %ecx,nb430_jnrc(%esp)
        movl %edx,nb430_jnrd(%esp)

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx
        leal  (%edx,%edx,2),%edx

        ## move four coordinates to xmm0-xmm2   

        movlps (%esi,%eax,4),%xmm4
        movlps (%esi,%ecx,4),%xmm5
        movss 8(%esi,%eax,4),%xmm2
        movss 8(%esi,%ecx,4),%xmm6

        movhps (%esi,%ebx,4),%xmm4
        movhps (%esi,%edx,4),%xmm5

        movss 8(%esi,%ebx,4),%xmm0
        movss 8(%esi,%edx,4),%xmm1

        shufps $0,%xmm0,%xmm2
        shufps $0,%xmm1,%xmm6

        movaps %xmm4,%xmm0
        movaps %xmm4,%xmm1

        shufps $136,%xmm6,%xmm2 ## constant 10001000

        shufps $136,%xmm5,%xmm0 ## constant 10001000
        shufps $221,%xmm5,%xmm1 ## constant 11011101            

        ## move ix-iz to xmm4-xmm6 
        movaps nb430_ix(%esp),%xmm4
        movaps nb430_iy(%esp),%xmm5
        movaps nb430_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb430_dx(%esp)
        movaps %xmm5,nb430_dy(%esp)
        movaps %xmm6,nb430_dz(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        ## rsq in xmm4 

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb430_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb430_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r 
        movaps %xmm4,nb430_r(%esp)
        mulps nb430_gbscale(%esp),%xmm4

        movhlps %xmm4,%xmm5
        cvttps2pi %xmm4,%mm6
        cvttps2pi %xmm5,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        cvtpi2ps %mm7,%xmm5
        movlhps %xmm5,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $2,%mm6
        pslld $2,%mm7

        movd %eax,%mm0
        movd %ebx,%mm1
        movd %ecx,%mm2
        movd %edx,%mm3

        movl nb430_GBtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        ## load coulomb table
        movaps (%esi,%eax,4),%xmm4
        movaps (%esi,%ebx,4),%xmm5
        movaps (%esi,%ecx,4),%xmm6
        movaps (%esi,%edx,4),%xmm7
        ## transpose, using xmm3 for scratch
        movaps %xmm6,%xmm3
        shufps $0xEE,%xmm7,%xmm3
        shufps $0x44,%xmm7,%xmm6
        movaps %xmm4,%xmm7
        shufps $0xEE,%xmm5,%xmm7
        shufps $0x44,%xmm5,%xmm4
        movaps %xmm4,%xmm5
        shufps $0xDD,%xmm6,%xmm5
        shufps $0x88,%xmm6,%xmm4
        movaps %xmm7,%xmm6
        shufps $0x88,%xmm3,%xmm6
        shufps $0xDD,%xmm3,%xmm7
        ## coulomb table ready, in xmm4-xmm7            

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  nb430_two(%esp),%xmm7    ## two*Heps2 
        movaps nb430_qq(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 

        ## get jnr from stack
        movl nb430_jnra(%esp),%eax
        movl nb430_jnrb(%esp),%ebx
        movl nb430_jnrc(%esp),%ecx
        movl nb430_jnrd(%esp),%edx

        movl nb430_dvda(%ebp),%esi

        ## Calculate dVda
        xorps %xmm7,%xmm7
        mulps nb430_gbscale(%esp),%xmm3
        movaps %xmm3,%xmm6
        mulps  nb430_r(%esp),%xmm6
        addps  %xmm5,%xmm6
        addps  nb430_vctot(%esp),%xmm5
        movaps %xmm5,nb430_vctot(%esp)

        ## xmm6=(vcoul+fijC*r)
        subps  %xmm6,%xmm7
        movaps %xmm7,%xmm6

        ## update dvdasum
        addps  nb430_dvdasum(%esp),%xmm7
        movaps %xmm7,nb430_dvdasum(%esp)

        ## update j atoms dvdaj
        movhlps %xmm6,%xmm7
        movaps  %xmm6,%xmm5
        movaps  %xmm7,%xmm4
        shufps $0x1,%xmm5,%xmm5
        shufps $0x1,%xmm4,%xmm4
        ## xmm6=dvdaj1 xmm5=dvdaj2 xmm7=dvdaj3 xmm4=dvdaj4
        addss  (%esi,%eax,4),%xmm6
        addss  (%esi,%ebx,4),%xmm5
        addss  (%esi,%ecx,4),%xmm7
        addss  (%esi,%edx,4),%xmm4
        movss  %xmm6,(%esi,%eax,4)
        movss  %xmm5,(%esi,%ebx,4)
        movss  %xmm7,(%esi,%ecx,4)
        movss  %xmm4,(%esi,%edx,4)

        ## put scalar force on stack temporarily 
        movaps %xmm3,nb430_fscal(%esp)

        movaps nb430_r(%esp),%xmm4
        mulps nb430_tsc(%esp),%xmm4

        movhlps %xmm4,%xmm5
        cvttps2pi %xmm4,%mm6
        cvttps2pi %xmm5,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        cvtpi2ps %mm7,%xmm5
        movlhps %xmm5,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $3,%mm6
        pslld $3,%mm7

        movl nb430_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        ## dispersion 
        movaps (%esi,%eax,4),%xmm4
        movaps (%esi,%ebx,4),%xmm5
        movaps (%esi,%ecx,4),%xmm6
        movaps (%esi,%edx,4),%xmm7
        ## transpose, using xmm3 for scratch
        movaps %xmm6,%xmm3
        shufps $0xEE,%xmm7,%xmm3
        shufps $0x44,%xmm7,%xmm6
        movaps %xmm4,%xmm7
        shufps $0xEE,%xmm5,%xmm7
        shufps $0x44,%xmm5,%xmm4
        movaps %xmm4,%xmm5
        shufps $0xDD,%xmm6,%xmm5
        shufps $0x88,%xmm6,%xmm4
        movaps %xmm7,%xmm6
        shufps $0x88,%xmm3,%xmm6
        shufps $0xDD,%xmm3,%xmm7
        ## dispersion table ready, in xmm4-xmm7         
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  nb430_two(%esp),%xmm7    ## two*Heps2 
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb430_c6(%esp),%xmm4
        mulps  %xmm4,%xmm7       ## fijD 
        mulps  %xmm4,%xmm5       ## Vvdw6
        mulps  nb430_tsc(%esp),%xmm7
        addps  nb430_fscal(%esp),%xmm7   ## add to fscal 

        ## put scalar force on stack Update Vvdwtot directly 
        addps  nb430_Vvdwtot(%esp),%xmm5
        movaps %xmm7,nb430_fscal(%esp)
        movaps %xmm5,nb430_Vvdwtot(%esp)

        ## repulsion 
        movaps 16(%esi,%eax,4),%xmm4
        movaps 16(%esi,%ebx,4),%xmm5
        movaps 16(%esi,%ecx,4),%xmm6
        movaps 16(%esi,%edx,4),%xmm7
        ## transpose, using xmm3 for scratch
        movaps %xmm6,%xmm3
        shufps $0xEE,%xmm7,%xmm3
        shufps $0x44,%xmm7,%xmm6
        movaps %xmm4,%xmm7
        shufps $0xEE,%xmm5,%xmm7
        shufps $0x44,%xmm5,%xmm4
        movaps %xmm4,%xmm5
        shufps $0xDD,%xmm6,%xmm5
        shufps $0x88,%xmm6,%xmm4
        movaps %xmm7,%xmm6
        shufps $0x88,%xmm3,%xmm6
        shufps $0xDD,%xmm3,%xmm7
        ## table ready, in xmm4-xmm7    
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  nb430_two(%esp),%xmm7    ## two*Heps2 
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb430_c12(%esp),%xmm4
        mulps  %xmm4,%xmm7 ## fijR 
        mulps  %xmm4,%xmm5 ## Vvdw12
        mulps nb430_tsc(%esp),%xmm7
        addps  nb430_fscal(%esp),%xmm7

        addps  nb430_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb430_Vvdwtot(%esp)
        xorps  %xmm4,%xmm4

        mulps %xmm0,%xmm7
        subps  %xmm7,%xmm4

        movaps nb430_dx(%esp),%xmm0
        movaps nb430_dy(%esp),%xmm1
        movaps nb430_dz(%esp),%xmm2

        movd %mm0,%eax
        movd %mm1,%ebx
        movd %mm2,%ecx
        movd %mm3,%edx

        movl   nb430_faction(%ebp),%edi
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movaps nb430_fix(%esp),%xmm3
        movaps nb430_fiy(%esp),%xmm4
        movaps nb430_fiz(%esp),%xmm5
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm5
        movaps %xmm3,nb430_fix(%esp)
        movaps %xmm4,nb430_fiy(%esp)
        movaps %xmm5,nb430_fiz(%esp)
        ## the fj's - start by accumulating x & y forces from memory 
        movlps (%edi,%eax,4),%xmm4
        movlps (%edi,%ecx,4),%xmm6
        movhps (%edi,%ebx,4),%xmm4
        movhps (%edi,%edx,4),%xmm6

        movaps %xmm4,%xmm3
        shufps $136,%xmm6,%xmm3 ## constant 10001000
        shufps $221,%xmm6,%xmm4 ## constant 11011101                          

        ## now xmm3-xmm5 contains fjx, fjy, fjz 
        subps  %xmm0,%xmm3
        subps  %xmm1,%xmm4

        ## unpack them back so we can store them - first x & y in xmm3/xmm4 

        movaps %xmm3,%xmm6
        unpcklps %xmm4,%xmm6
        unpckhps %xmm4,%xmm3
        ## xmm6(l)=x & y for j1, (h) for j2 
        ## xmm3(l)=x & y for j3, (h) for j4 
        movlps %xmm6,(%edi,%eax,4)
        movlps %xmm3,(%edi,%ecx,4)

        movhps %xmm6,(%edi,%ebx,4)
        movhps %xmm3,(%edi,%edx,4)

        ## and the z forces 
        movss  8(%edi,%eax,4),%xmm4
        movss  8(%edi,%ebx,4),%xmm5
        movss  8(%edi,%ecx,4),%xmm6
        movss  8(%edi,%edx,4),%xmm7
        subss  %xmm2,%xmm4
        shufps $229,%xmm2,%xmm2 ## constant 11100101
        subss  %xmm2,%xmm5
        shufps $234,%xmm2,%xmm2 ## constant 11101010
        subss  %xmm2,%xmm6
        shufps $255,%xmm2,%xmm2 ## constant 11111111
        subss  %xmm2,%xmm7
        movss  %xmm4,8(%edi,%eax,4)
        movss  %xmm5,8(%edi,%ebx,4)
        movss  %xmm6,8(%edi,%ecx,4)
        movss  %xmm7,8(%edi,%edx,4)

        ## should we do one more iteration? 
        subl $4,nb430_innerk(%esp)
        jl    _nb_kernel430_ia32_sse.nb430_finish_inner
        jmp   _nb_kernel430_ia32_sse.nb430_unroll_loop
_nb_kernel430_ia32_sse.nb430_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb430_innerk(%esp)
        movl  nb430_innerk(%esp),%edx
        andl  $2,%edx
        jnz   _nb_kernel430_ia32_sse.nb430_dopair
        jmp   _nb_kernel430_ia32_sse.nb430_checksingle
_nb_kernel430_ia32_sse.nb430_dopair: 

        movl  nb430_innerjjnr(%esp),%ecx

        movl  (%ecx),%eax
        movl  4(%ecx),%ebx
        addl $8,nb430_innerjjnr(%esp)

        xorps %xmm2,%xmm2
        movaps %xmm2,%xmm6

        ## load isaj
        movl nb430_invsqrta(%ebp),%esi
        movss (%esi,%eax,4),%xmm2
        movss (%esi,%ebx,4),%xmm3
        unpcklps %xmm3,%xmm2    ## isaj in xmm3(0,1)
        mulps  nb430_isai(%esp),%xmm2
        movaps %xmm2,nb430_isaprod(%esp)
        movaps %xmm2,%xmm1
        mulps nb430_gbtsc(%esp),%xmm1
        movaps %xmm1,nb430_gbscale(%esp)

        movl nb430_charge(%ebp),%esi     ## base of charge[]    
        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ebx,4),%xmm6
        unpcklps %xmm6,%xmm3 ## constant 00001000 ;# xmm3(0,1) has the charges 

        mulps  nb430_iq(%esp),%xmm2
        mulps  %xmm2,%xmm3
        movaps %xmm3,nb430_qq(%esp)

        movl nb430_type(%ebp),%esi
        movl  %eax,%ecx
        movl  %ebx,%edx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb430_vdwparam(%ebp),%esi
        shll %ecx
        shll %edx
        movl nb430_ntia(%esp),%edi
        addl %edi,%ecx
        addl %edi,%edx
        movlps (%esi,%ecx,4),%xmm6
        movhps (%esi,%edx,4),%xmm6
        movl nb430_pos(%ebp),%edi

        movaps %xmm6,%xmm4
        shufps $8,%xmm4,%xmm4 ## constant 00001000       
        shufps $13,%xmm6,%xmm6 ## constant 00001101
        movlhps %xmm7,%xmm4
        movlhps %xmm7,%xmm6

        movaps %xmm4,nb430_c6(%esp)
        movaps %xmm6,nb430_c12(%esp)

        movd  %eax,%mm0         ## copy jnr to mm0/mm1
        movd  %ebx,%mm1

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        ## move coordinates to xmm0-xmm2 
        movlps (%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2
        movhps (%edi,%ebx,4),%xmm1
        movss 8(%edi,%ebx,4),%xmm0

        movlhps %xmm7,%xmm3

        shufps $0,%xmm0,%xmm2

        movaps %xmm1,%xmm0

        shufps $136,%xmm2,%xmm2 ## constant 10001000

        shufps $136,%xmm0,%xmm0 ## constant 10001000
        shufps $221,%xmm1,%xmm1 ## constant 11011101

        movl   nb430_faction(%ebp),%edi
        ## move ix-iz to xmm4-xmm6 
        xorps   %xmm7,%xmm7

        movaps nb430_ix(%esp),%xmm4
        movaps nb430_iy(%esp),%xmm5
        movaps nb430_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb430_dx(%esp)
        movaps %xmm5,nb430_dy(%esp)
        movaps %xmm6,nb430_dz(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        ## rsq in xmm4 

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb430_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb430_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r 
        movaps %xmm4,nb430_r(%esp)
        mulps nb430_gbscale(%esp),%xmm4

        cvttps2pi %xmm4,%mm6    ## mm6 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6

        movl nb430_GBtab(%ebp),%esi
        movd %mm6,%ecx
        psrlq $32,%mm6
        movd %mm6,%edx

        ## load coulomb table
        movaps (%esi,%ecx,4),%xmm4
        movaps (%esi,%edx,4),%xmm7
        ## transpose, using xmm3 for scratch
        movaps %xmm4,%xmm6
        unpcklps %xmm7,%xmm4    ## Y1 Y2 F1 F2 
        unpckhps %xmm7,%xmm6    ## G1 G2 H1 H2
        movhlps  %xmm4,%xmm5    ## F1 F2 
        movhlps  %xmm6,%xmm7    ## H1 H2
        ## coulomb table ready, in xmm4-xmm7    

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  nb430_two(%esp),%xmm7    ## two*Heps2 
        movaps nb430_qq(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 

        ## get jnr from mm0/mm1
        movd %mm0,%ecx
        movd %mm1,%edx

        movl nb430_dvda(%ebp),%esi

        ## Calculate dVda
        xorps %xmm7,%xmm7
        mulps nb430_gbscale(%esp),%xmm3
        movaps %xmm3,%xmm6
        mulps  nb430_r(%esp),%xmm6
        addps  %xmm5,%xmm6
        addps  nb430_vctot(%esp),%xmm5
        movaps %xmm5,nb430_vctot(%esp)

        ## xmm6=(vcoul+fijC*r)
        subps  %xmm6,%xmm7
        movaps %xmm7,%xmm6

        ## update dvdasum
        addps  nb430_dvdasum(%esp),%xmm7
        movaps %xmm7,nb430_dvdasum(%esp)

        ## update j atoms dvdaj
        movaps %xmm6,%xmm7
        shufps $0x1,%xmm7,%xmm7
        addss  (%esi,%ecx,4),%xmm6
        addss  (%esi,%edx,4),%xmm7
        movss  %xmm6,(%esi,%ecx,4)
        movss  %xmm7,(%esi,%edx,4)

        ## put scalar force on stack temporarily 
        movaps %xmm3,nb430_fscal(%esp)

        movaps nb430_r(%esp),%xmm4
        mulps nb430_tsc(%esp),%xmm4

        cvttps2pi %xmm4,%mm6
        cvtpi2ps %mm6,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $3,%mm6

        movl nb430_VFtab(%ebp),%esi
        movd %mm6,%ecx
        psrlq $32,%mm6
        movd %mm6,%edx

        ## dispersion 
        movaps (%esi,%ecx,4),%xmm4
        movaps (%esi,%edx,4),%xmm7
        ## transpose, using xmm3 for scratch
        movaps %xmm4,%xmm6
        unpcklps %xmm7,%xmm4    ## Y1 Y2 F1 F2 
        unpckhps %xmm7,%xmm6    ## G1 G2 H1 H2
        movhlps  %xmm4,%xmm5    ## F1 F2 
        movhlps  %xmm6,%xmm7    ## H1 H2
        ## dispersion table ready, in xmm4-xmm7         
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  nb430_two(%esp),%xmm7    ## two*Heps2 
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb430_c6(%esp),%xmm4
        mulps  %xmm4,%xmm7       ## fijD 
        mulps  %xmm4,%xmm5       ## Vvdw6 
        mulps  nb430_tsc(%esp),%xmm7
        addps  nb430_fscal(%esp),%xmm7   ## add to fscal 

        ## put scalar force on stack Update Vvdwtot directly 
        addps  nb430_Vvdwtot(%esp),%xmm5
        movaps %xmm7,nb430_fscal(%esp)
        movaps %xmm5,nb430_Vvdwtot(%esp)

        ## repulsion 
        movaps 16(%esi,%ecx,4),%xmm4
        movaps 16(%esi,%edx,4),%xmm7
        ## transpose, using xmm3 for scratch
        movaps %xmm4,%xmm6
        unpcklps %xmm7,%xmm4    ## Y1 Y2 F1 F2 
        unpckhps %xmm7,%xmm6    ## G1 G2 H1 H2
        movhlps  %xmm4,%xmm5    ## F1 F2 
        movhlps  %xmm6,%xmm7    ## H1 H2
        ## table ready, in xmm4-xmm7    
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  nb430_two(%esp),%xmm7    ## two*Heps2 
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb430_c12(%esp),%xmm4
        mulps  %xmm4,%xmm7 ## fijR 
        mulps  %xmm4,%xmm5 ## Vvdw12 
        mulps  nb430_tsc(%esp),%xmm7
        addps  nb430_fscal(%esp),%xmm7

        addps  nb430_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb430_Vvdwtot(%esp)
        xorps  %xmm4,%xmm4

        mulps %xmm0,%xmm7
        subps  %xmm7,%xmm4

        movaps nb430_dx(%esp),%xmm0
        movaps nb430_dy(%esp),%xmm1
        movaps nb430_dz(%esp),%xmm2

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movaps nb430_fix(%esp),%xmm3
        movaps nb430_fiy(%esp),%xmm4
        movaps nb430_fiz(%esp),%xmm5
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm5
        movaps %xmm3,nb430_fix(%esp)
        movaps %xmm4,nb430_fiy(%esp)
        movaps %xmm5,nb430_fiz(%esp)
        ## update the fj's 
        movss   (%edi,%eax,4),%xmm3
        movss   4(%edi,%eax,4),%xmm4
        movss   8(%edi,%eax,4),%xmm5
        subss   %xmm0,%xmm3
        subss   %xmm1,%xmm4
        subss   %xmm2,%xmm5
        movss   %xmm3,(%edi,%eax,4)
        movss   %xmm4,4(%edi,%eax,4)
        movss   %xmm5,8(%edi,%eax,4)

        shufps $225,%xmm0,%xmm0 ## constant 11100001
        shufps $225,%xmm1,%xmm1 ## constant 11100001
        shufps $225,%xmm2,%xmm2 ## constant 11100001

        movss   (%edi,%ebx,4),%xmm3
        movss   4(%edi,%ebx,4),%xmm4
        movss   8(%edi,%ebx,4),%xmm5
        subss   %xmm0,%xmm3
        subss   %xmm1,%xmm4
        subss   %xmm2,%xmm5
        movss   %xmm3,(%edi,%ebx,4)
        movss   %xmm4,4(%edi,%ebx,4)
        movss   %xmm5,8(%edi,%ebx,4)

_nb_kernel430_ia32_sse.nb430_checksingle:       
        movl  nb430_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel430_ia32_sse.nb430_dosingle
        jmp    _nb_kernel430_ia32_sse.nb430_updateouterdata
_nb_kernel430_ia32_sse.nb430_dosingle: 
        movl nb430_charge(%ebp),%esi
        movl nb430_invsqrta(%ebp),%edx
        movl nb430_pos(%ebp),%edi
        movl  nb430_innerjjnr(%esp),%ecx
        movl  (%ecx),%eax
        xorps  %xmm2,%xmm2
        movaps %xmm2,%xmm6
        movss (%edx,%eax,4),%xmm2       ## isaj
        mulss nb430_isai(%esp),%xmm2
        movss %xmm2,nb430_isaprod(%esp)
        movss %xmm2,%xmm1
        mulss nb430_gbtsc(%esp),%xmm1
        movss %xmm1,nb430_gbscale(%esp)

        mulss  nb430_iq(%esp),%xmm2
        movss (%esi,%eax,4),%xmm6       ## xmm6(0) has the charge       
        mulss  %xmm2,%xmm6
        movss %xmm6,nb430_qq(%esp)

        movl nb430_type(%ebp),%esi
        movl %eax,%ecx
        movl (%esi,%ecx,4),%ecx
        movl nb430_vdwparam(%ebp),%esi
        shll %ecx
        addl nb430_ntia(%esp),%ecx
        movlps (%esi,%ecx,4),%xmm6
        movaps %xmm6,%xmm4
        shufps $252,%xmm4,%xmm4 ## constant 11111100    
        shufps $253,%xmm6,%xmm6 ## constant 11111101    

        movss %xmm4,nb430_c6(%esp)
        movss %xmm6,nb430_c12(%esp)

        movd  %eax,%mm0
        leal  (%eax,%eax,2),%eax

        ## move coordinates to xmm0-xmm2 
        movss (%edi,%eax,4),%xmm0
        movss 4(%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2

        movss nb430_ix(%esp),%xmm4
        movss nb430_iy(%esp),%xmm5
        movss nb430_iz(%esp),%xmm6

        ## calc dr 
        subss %xmm0,%xmm4
        subss %xmm1,%xmm5
        subss %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb430_dx(%esp)
        movaps %xmm5,nb430_dy(%esp)
        movaps %xmm6,nb430_dz(%esp)
        ## square it 
        mulss %xmm4,%xmm4
        mulss %xmm5,%xmm5
        mulss %xmm6,%xmm6
        addss %xmm5,%xmm4
        addss %xmm6,%xmm4
        ## rsq in xmm4 

        rsqrtss %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulss %xmm5,%xmm5
        movss nb430_three(%esp),%xmm1
        mulss %xmm4,%xmm5       ## rsq*lu*lu                    
        movss nb430_half(%esp),%xmm0
        subss %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulss %xmm2,%xmm1
        mulss %xmm1,%xmm0       ## xmm0=rinv 

        mulss %xmm0,%xmm4       ## xmm4=r 
        movss %xmm4,nb430_r(%esp)
        mulss nb430_gbscale(%esp),%xmm4

        cvttss2si %xmm4,%ebx    ## mm6 contain lu indices 
        cvtsi2ss %ebx,%xmm6
        subss %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%ebx

        movl nb430_GBtab(%ebp),%esi

        movaps (%esi,%ebx,4),%xmm4
        movhlps %xmm4,%xmm6
        movaps %xmm4,%xmm5
        movaps %xmm6,%xmm7
        shufps $1,%xmm5,%xmm5
        shufps $1,%xmm7,%xmm7
        ## table ready in xmm4-xmm7 

        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp      
        mulss  nb430_two(%esp),%xmm7    ## two*Heps2 
        movss nb430_qq(%esp),%xmm3
        addss  %xmm6,%xmm7
        addss  %xmm5,%xmm7 ## xmm7=FF 
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 
        mulss  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulss  %xmm7,%xmm3 ## fijC=FF*qq 

        movd %mm0,%ebx
        movl nb430_dvda(%ebp),%esi

        ## Calculate dVda
        xorps %xmm7,%xmm7
        mulss nb430_gbscale(%esp),%xmm3
        movaps %xmm3,%xmm6
        mulss  nb430_r(%esp),%xmm6
        addss  %xmm5,%xmm6
        addss  nb430_vctot(%esp),%xmm5
        movss %xmm5,nb430_vctot(%esp)


        ## xmm6=(vcoul+fijC*r)
        subss  %xmm6,%xmm7
        movaps %xmm7,%xmm6

        ## update dvdasum
        addss  nb430_dvdasum(%esp),%xmm7
        movaps %xmm7,nb430_dvdasum(%esp)

        ## update j atoms dvdaj
        addss  (%esi,%ebx,4),%xmm6
        movss  %xmm6,(%esi,%ebx,4)

        ## put scalar force on stack temporarily 
        movss %xmm3,nb430_fscal(%esp)

        movss nb430_r(%esp),%xmm4
        mulps nb430_tsc(%esp),%xmm4

        cvttss2si %xmm4,%ebx
        cvtsi2ss %ebx,%xmm6
        subss %xmm6,%xmm4
        movss %xmm4,%xmm1       ## xmm1=eps 
        movss %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=eps2 

        shll $3,%ebx
        movl nb430_VFtab(%ebp),%esi

        ## dispersion 
        movaps (%esi,%ebx,4),%xmm4
        movhlps %xmm4,%xmm6
        movaps %xmm4,%xmm5
        movaps %xmm6,%xmm7
        shufps $1,%xmm5,%xmm5
        shufps $1,%xmm7,%xmm7
        ## table ready in xmm4-xmm7 

        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp      
        mulss  nb430_two(%esp),%xmm7    ## two*Heps2 
        addss  %xmm6,%xmm7
        addss  %xmm5,%xmm7 ## xmm7=FF 
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 

        movss nb430_c6(%esp),%xmm4
        mulss  %xmm4,%xmm7       ## fijD 
        mulss  %xmm4,%xmm5       ## Vvdw6
        mulps  nb430_tsc(%esp),%xmm7
        addss  nb430_fscal(%esp),%xmm7   ## add to fscal 

        ## put scalar force on stack Update Vvdwtot directly 
        addss  nb430_Vvdwtot(%esp),%xmm5
        movss %xmm7,nb430_fscal(%esp)
        movss %xmm5,nb430_Vvdwtot(%esp)

        ## repulsion 
        movaps 16(%esi,%ebx,4),%xmm4
        movhlps %xmm4,%xmm6
        movaps %xmm4,%xmm5
        movaps %xmm6,%xmm7
        shufps $1,%xmm5,%xmm5
        shufps $1,%xmm7,%xmm7
        ## table ready in xmm4-xmm7 

        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp      
        mulss  nb430_two(%esp),%xmm7    ## two*Heps2 
        addss  %xmm6,%xmm7
        addss  %xmm5,%xmm7 ## xmm7=FF 
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 

        movss nb430_c12(%esp),%xmm4
        mulss  %xmm4,%xmm7 ## fijR 
        mulss  %xmm4,%xmm5 ## Vvdw12 
        mulps  nb430_tsc(%esp),%xmm7
        addss  nb430_fscal(%esp),%xmm7

        addss  nb430_Vvdwtot(%esp),%xmm5
        movss %xmm5,nb430_Vvdwtot(%esp)
        xorps  %xmm4,%xmm4

        mulss %xmm0,%xmm7
        subss  %xmm7,%xmm4
        movl   nb430_faction(%ebp),%edi

        movss nb430_dx(%esp),%xmm0
        movss nb430_dy(%esp),%xmm1
        movss nb430_dz(%esp),%xmm2

        mulss  %xmm4,%xmm0
        mulss  %xmm4,%xmm1
        mulss  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movss nb430_fix(%esp),%xmm3
        movss nb430_fiy(%esp),%xmm4
        movss nb430_fiz(%esp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss %xmm3,nb430_fix(%esp)
        movss %xmm4,nb430_fiy(%esp)
        movss %xmm5,nb430_fiz(%esp)
        ## update fj 

        movss   (%edi,%eax,4),%xmm3
        movss   4(%edi,%eax,4),%xmm4
        movss   8(%edi,%eax,4),%xmm5
        subss   %xmm0,%xmm3
        subss   %xmm1,%xmm4
        subss   %xmm2,%xmm5
        movss   %xmm3,(%edi,%eax,4)
        movss   %xmm4,4(%edi,%eax,4)
        movss   %xmm5,8(%edi,%eax,4)
_nb_kernel430_ia32_sse.nb430_updateouterdata: 
        movl  nb430_ii3(%esp),%ecx
        movl  nb430_faction(%ebp),%edi
        movl  nb430_fshift(%ebp),%esi
        movl  nb430_is3(%esp),%edx

        ## accumulate i forces in xmm0, xmm1, xmm2 
        movaps nb430_fix(%esp),%xmm0
        movaps nb430_fiy(%esp),%xmm1
        movaps nb430_fiz(%esp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addps  %xmm3,%xmm0
        addps  %xmm4,%xmm1
        addps  %xmm5,%xmm2 ## sum is in 1/2 in xmm0-xmm2 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        shufps $1,%xmm3,%xmm3
        shufps $1,%xmm4,%xmm4
        shufps $1,%xmm5,%xmm5
        addss  %xmm3,%xmm0
        addss  %xmm4,%xmm1
        addss  %xmm5,%xmm2      ## xmm0-xmm2 has single force in pos0 

        ## increment i force 
        movss  (%edi,%ecx,4),%xmm3
        movss  4(%edi,%ecx,4),%xmm4
        movss  8(%edi,%ecx,4),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,(%edi,%ecx,4)
        movss  %xmm4,4(%edi,%ecx,4)
        movss  %xmm5,8(%edi,%ecx,4)

        ## increment fshift force  
        movss  (%esi,%edx,4),%xmm3
        movss  4(%esi,%edx,4),%xmm4
        movss  8(%esi,%edx,4),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,(%esi,%edx,4)
        movss  %xmm4,4(%esi,%edx,4)
        movss  %xmm5,8(%esi,%edx,4)

        ## get n from stack
        movl nb430_n(%esp),%esi
        ## get group index for i particle 
        movl  nb430_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb430_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb430_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb430_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb430_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate dVda and update it 
        movaps nb430_dvdasum(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        movl nb430_ii(%esp),%edx
        movl nb430_dvda(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb430_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel430_ia32_sse.nb430_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb430_n(%esp)
        jmp _nb_kernel430_ia32_sse.nb430_outer
_nb_kernel430_ia32_sse.nb430_outerend: 
        ## check if more outer neighborlists remain
        movl  nb430_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel430_ia32_sse.nb430_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel430_ia32_sse.nb430_threadloop
_nb_kernel430_ia32_sse.nb430_end: 
        emms

        movl nb430_nouter(%esp),%eax
        movl nb430_ninner(%esp),%ebx
        movl nb430_outeriter(%ebp),%ecx
        movl nb430_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb430_salign(%esp),%eax
        addl %eax,%esp
        addl $488,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret







.globl nb_kernel430nf_ia32_sse
.globl _nb_kernel430nf_ia32_sse
nb_kernel430nf_ia32_sse:        
_nb_kernel430nf_ia32_sse:       
.set nb430nf_p_nri, 8
.set nb430nf_iinr, 12
.set nb430nf_jindex, 16
.set nb430nf_jjnr, 20
.set nb430nf_shift, 24
.set nb430nf_shiftvec, 28
.set nb430nf_fshift, 32
.set nb430nf_gid, 36
.set nb430nf_pos, 40
.set nb430nf_faction, 44
.set nb430nf_charge, 48
.set nb430nf_p_facel, 52
.set nb430nf_argkrf, 56
.set nb430nf_argcrf, 60
.set nb430nf_Vc, 64
.set nb430nf_type, 68
.set nb430nf_p_ntype, 72
.set nb430nf_vdwparam, 76
.set nb430nf_Vvdw, 80
.set nb430nf_p_tabscale, 84
.set nb430nf_VFtab, 88
.set nb430nf_invsqrta, 92
.set nb430nf_dvda, 96
.set nb430nf_p_gbtabscale, 100
.set nb430nf_GBtab, 104
.set nb430nf_p_nthreads, 108
.set nb430nf_count, 112
.set nb430nf_mtx, 116
.set nb430nf_outeriter, 120
.set nb430nf_inneriter, 124
.set nb430nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb430nf_ix, 0
.set nb430nf_iy, 16
.set nb430nf_iz, 32
.set nb430nf_iq, 48
.set nb430nf_gbtsc, 64
.set nb430nf_tsc, 80
.set nb430nf_qq, 96
.set nb430nf_c6, 112
.set nb430nf_c12, 128
.set nb430nf_vctot, 144
.set nb430nf_Vvdwtot, 160
.set nb430nf_half, 176
.set nb430nf_three, 192
.set nb430nf_isai, 208
.set nb430nf_isaprod, 224
.set nb430nf_gbscale, 240
.set nb430nf_r, 256
.set nb430nf_is3, 272
.set nb430nf_ii3, 276
.set nb430nf_ntia, 280
.set nb430nf_innerjjnr, 284
.set nb430nf_innerk, 288
.set nb430nf_n, 292
.set nb430nf_nn1, 296
.set nb430nf_nri, 300
.set nb430nf_facel, 304
.set nb430nf_ntype, 308
.set nb430nf_nouter, 312
.set nb430nf_ninner, 316
.set nb430nf_salign, 320
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $324,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb430nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb430nf_p_nri(%ebp),%ecx
        movl nb430nf_p_facel(%ebp),%esi
        movl nb430nf_p_ntype(%ebp),%edi
        movl (%ecx),%ecx
        movl (%esi),%esi
        movl (%edi),%edi
        movl %ecx,nb430nf_nri(%esp)
        movl %esi,nb430nf_facel(%esp)
        movl %edi,nb430nf_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb430nf_nouter(%esp)
        movl %eax,nb430nf_ninner(%esp)


        movl nb430nf_p_gbtabscale(%ebp),%eax
        movss (%eax),%xmm3
        movl nb430nf_p_tabscale(%ebp),%eax
        movss (%eax),%xmm4
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb430nf_gbtsc(%esp)
        movaps %xmm4,nb430nf_tsc(%esp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb430nf_half(%esp)
        movss nb430nf_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb430nf_half(%esp)
        movaps %xmm3,nb430nf_three(%esp)

_nb_kernel430nf_ia32_sse.nb430nf_threadloop: 
        movl  nb430nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel430nf_ia32_sse.nb430nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel430nf_ia32_sse.nb430nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb430nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb430nf_n(%esp)
        movl %ebx,nb430nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel430nf_ia32_sse.nb430nf_outerstart
        jmp _nb_kernel430nf_ia32_sse.nb430nf_end

_nb_kernel430nf_ia32_sse.nb430nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb430nf_nouter(%esp),%ebx
        movl %ebx,nb430nf_nouter(%esp)

_nb_kernel430nf_ia32_sse.nb430nf_outer: 
        movl  nb430nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb430nf_is3(%esp)            ## store is3 

        movl  nb430nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb430nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 

        movl  nb430nf_charge(%ebp),%edx
        movss (%edx,%ebx,4),%xmm3
        mulss nb430nf_facel(%esp),%xmm3
        shufps $0,%xmm3,%xmm3

        movl  nb430nf_invsqrta(%ebp),%edx       ## load invsqrta[ii]
        movss (%edx,%ebx,4),%xmm4
        shufps $0,%xmm4,%xmm4

        movl  nb430nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb430nf_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb430nf_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb430nf_pos(%ebp),%eax      ## eax = base of pos[]  

        addss (%eax,%ebx,4),%xmm0
        addss 4(%eax,%ebx,4),%xmm1
        addss 8(%eax,%ebx,4),%xmm2

        movaps %xmm3,nb430nf_iq(%esp)
        movaps %xmm4,nb430nf_isai(%esp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb430nf_ix(%esp)
        movaps %xmm1,nb430nf_iy(%esp)
        movaps %xmm2,nb430nf_iz(%esp)

        movl  %ebx,nb430nf_ii3(%esp)

        ## clear vctot 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb430nf_vctot(%esp)
        movaps %xmm4,nb430nf_Vvdwtot(%esp)

        movl  nb430nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb430nf_pos(%ebp),%esi
        movl  nb430nf_faction(%ebp),%edi
        movl  nb430nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb430nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb430nf_ninner(%esp),%ecx
        movl  %ecx,nb430nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb430nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel430nf_ia32_sse.nb430nf_unroll_loop
        jmp   _nb_kernel430nf_ia32_sse.nb430nf_finish_inner
_nb_kernel430nf_ia32_sse.nb430nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb430nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 
        addl $16,nb430nf_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        ## load isa2
        movl nb430nf_invsqrta(%ebp),%esi
        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7
        movaps nb430nf_isai(%esp),%xmm2
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        mulps  %xmm3,%xmm2

        movaps %xmm2,nb430nf_isaprod(%esp)
        movaps %xmm2,%xmm1
        mulps nb430nf_gbtsc(%esp),%xmm1
        movaps %xmm1,nb430nf_gbscale(%esp)

        movl nb430nf_charge(%ebp),%esi     ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        mulps nb430nf_iq(%esp),%xmm2
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        mulps  %xmm2,%xmm3
        movaps %xmm3,nb430nf_qq(%esp)

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movl nb430nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb430nf_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb430nf_ntia(%esp),%edi
        addl %edi,%eax
        addl %edi,%ebx
        addl %edi,%ecx
        addl %edi,%edx

        movlps (%esi,%eax,4),%xmm6
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm6
        movhps (%esi,%edx,4),%xmm7

        movaps %xmm6,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm6 ## constant 11011101

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movd  %mm2,%ecx
        movd  %mm3,%edx

        movaps %xmm4,nb430nf_c6(%esp)
        movaps %xmm6,nb430nf_c12(%esp)

        movl nb430nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        leal  (%ecx,%ecx,2),%ecx     ## replace jnr with j3 
        leal  (%edx,%edx,2),%edx

        ## move four coordinates to xmm0-xmm2   

        movlps (%esi,%eax,4),%xmm4
        movlps (%esi,%ecx,4),%xmm5
        movss 8(%esi,%eax,4),%xmm2
        movss 8(%esi,%ecx,4),%xmm6

        movhps (%esi,%ebx,4),%xmm4
        movhps (%esi,%edx,4),%xmm5

        movss 8(%esi,%ebx,4),%xmm0
        movss 8(%esi,%edx,4),%xmm1

        shufps $0,%xmm0,%xmm2
        shufps $0,%xmm1,%xmm6

        movaps %xmm4,%xmm0
        movaps %xmm4,%xmm1

        shufps $136,%xmm6,%xmm2 ## constant 10001000

        shufps $136,%xmm5,%xmm0 ## constant 10001000
        shufps $221,%xmm5,%xmm1 ## constant 11011101            

        ## move ix-iz to xmm4-xmm6 
        movaps nb430nf_ix(%esp),%xmm4
        movaps nb430nf_iy(%esp),%xmm5
        movaps nb430nf_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        ## rsq in xmm4 

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb430nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb430nf_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r
        movaps %xmm4,nb430nf_r(%esp)
        mulps nb430nf_gbscale(%esp),%xmm4

        movhlps %xmm4,%xmm5
        cvttps2pi %xmm4,%mm6
        cvttps2pi %xmm5,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        cvtpi2ps %mm7,%xmm5
        movlhps %xmm5,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $2,%mm6
        pslld $2,%mm7

        movd %eax,%mm0
        movd %ebx,%mm1
        movd %ecx,%mm2
        movd %edx,%mm3

        movl nb430nf_GBtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        ## load coulomb table
        movaps (%esi,%eax,4),%xmm4
        movaps (%esi,%ebx,4),%xmm5
        movaps (%esi,%ecx,4),%xmm6
        movaps (%esi,%edx,4),%xmm7
        ## transpose, using xmm3 for scratch
        movaps %xmm6,%xmm3
        shufps $0xEE,%xmm7,%xmm3
        shufps $0x44,%xmm7,%xmm6
        movaps %xmm4,%xmm7
        shufps $0xEE,%xmm5,%xmm7
        shufps $0x44,%xmm5,%xmm4
        movaps %xmm4,%xmm5
        shufps $0xDD,%xmm6,%xmm5
        shufps $0x88,%xmm6,%xmm4
        movaps %xmm7,%xmm6
        shufps $0x88,%xmm3,%xmm6
        shufps $0xDD,%xmm3,%xmm7
        ## coulomb table ready, in xmm4-xmm7            

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        movaps nb430nf_qq(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        addps  nb430nf_vctot(%esp),%xmm5
        movaps %xmm5,nb430nf_vctot(%esp)


        movaps nb430nf_r(%esp),%xmm4
        mulps nb430nf_tsc(%esp),%xmm4

        movhlps %xmm4,%xmm5
        cvttps2pi %xmm4,%mm6
        cvttps2pi %xmm5,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        cvtpi2ps %mm7,%xmm5
        movlhps %xmm5,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $3,%mm6
        pslld $3,%mm7

        movl nb430nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        ## dispersion 
        movaps (%esi,%eax,4),%xmm4
        movaps (%esi,%ebx,4),%xmm5
        movaps (%esi,%ecx,4),%xmm6
        movaps (%esi,%edx,4),%xmm7
        ## transpose, using xmm3 for scratch
        movaps %xmm6,%xmm3
        shufps $0xEE,%xmm7,%xmm3
        shufps $0x44,%xmm7,%xmm6
        movaps %xmm4,%xmm7
        shufps $0xEE,%xmm5,%xmm7
        shufps $0x44,%xmm5,%xmm4
        movaps %xmm4,%xmm5
        shufps $0xDD,%xmm6,%xmm5
        shufps $0x88,%xmm6,%xmm4
        movaps %xmm7,%xmm6
        shufps $0x88,%xmm3,%xmm6
        shufps $0xDD,%xmm3,%xmm7
        ## dispersion table ready, in xmm4-xmm7         
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  nb430nf_c6(%esp),%xmm5    ## Vvdw6
        addps  nb430nf_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb430nf_Vvdwtot(%esp)

        ## repulsion 
        movaps 16(%esi,%eax,4),%xmm4
        movaps 16(%esi,%ebx,4),%xmm5
        movaps 16(%esi,%ecx,4),%xmm6
        movaps 16(%esi,%edx,4),%xmm7
        ## transpose, using xmm3 for scratch
        movaps %xmm6,%xmm3
        shufps $0xEE,%xmm7,%xmm3
        shufps $0x44,%xmm7,%xmm6
        movaps %xmm4,%xmm7
        shufps $0xEE,%xmm5,%xmm7
        shufps $0x44,%xmm5,%xmm4
        movaps %xmm4,%xmm5
        shufps $0xDD,%xmm6,%xmm5
        shufps $0x88,%xmm6,%xmm4
        movaps %xmm7,%xmm6
        shufps $0x88,%xmm3,%xmm6
        shufps $0xDD,%xmm3,%xmm7
        ## table ready, in xmm4-xmm7    
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        mulps  nb430nf_c12(%esp),%xmm5   ## Vvdw12
        addps  nb430nf_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb430nf_Vvdwtot(%esp)

        ## should we do one more iteration? 
        subl $4,nb430nf_innerk(%esp)
        jl    _nb_kernel430nf_ia32_sse.nb430nf_finish_inner
        jmp   _nb_kernel430nf_ia32_sse.nb430nf_unroll_loop
_nb_kernel430nf_ia32_sse.nb430nf_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb430nf_innerk(%esp)
        movl  nb430nf_innerk(%esp),%edx
        andl  $2,%edx
        jnz   _nb_kernel430nf_ia32_sse.nb430nf_dopair
        jmp   _nb_kernel430nf_ia32_sse.nb430nf_checksingle
_nb_kernel430nf_ia32_sse.nb430nf_dopair: 

        movl  nb430nf_innerjjnr(%esp),%ecx

        movl  (%ecx),%eax
        movl  4(%ecx),%ebx
        addl $8,nb430nf_innerjjnr(%esp)

        xorps %xmm2,%xmm2
        movaps %xmm2,%xmm6

        ## load isa2
        movl nb430nf_invsqrta(%ebp),%esi
        movss (%esi,%eax,4),%xmm2
        movss (%esi,%ebx,4),%xmm3
        unpcklps %xmm3,%xmm2    ## isa2 in xmm3(0,1)
        mulps  nb430nf_isai(%esp),%xmm2
        movaps %xmm2,nb430nf_isaprod(%esp)
        movaps %xmm2,%xmm1
        mulps nb430nf_gbtsc(%esp),%xmm1
        movaps %xmm1,nb430nf_gbscale(%esp)

        movl nb430nf_charge(%ebp),%esi     ## base of charge[]  
        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ebx,4),%xmm6
        unpcklps %xmm6,%xmm3 ## constant 00001000 ;# xmm3(0,1) has the charges 

        mulps  nb430nf_iq(%esp),%xmm2
        mulps  %xmm2,%xmm3
        movaps %xmm3,nb430nf_qq(%esp)

        movl nb430nf_type(%ebp),%esi
        movl  %eax,%ecx
        movl  %ebx,%edx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb430nf_vdwparam(%ebp),%esi
        shll %ecx
        shll %edx
        movl nb430nf_ntia(%esp),%edi
        addl %edi,%ecx
        addl %edi,%edx
        movlps (%esi,%ecx,4),%xmm6
        movhps (%esi,%edx,4),%xmm6
        movl nb430nf_pos(%ebp),%edi

        movaps %xmm6,%xmm4
        shufps $8,%xmm4,%xmm4 ## constant 00001000       
        shufps $13,%xmm6,%xmm6 ## constant 00001101
        movlhps %xmm7,%xmm4
        movlhps %xmm7,%xmm6

        movaps %xmm4,nb430nf_c6(%esp)
        movaps %xmm6,nb430nf_c12(%esp)

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        ## move coordinates to xmm0-xmm2 
        movlps (%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2
        movhps (%edi,%ebx,4),%xmm1
        movss 8(%edi,%ebx,4),%xmm0

        movlhps %xmm7,%xmm3

        shufps $0,%xmm0,%xmm2

        movaps %xmm1,%xmm0

        shufps $136,%xmm2,%xmm2 ## constant 10001000

        shufps $136,%xmm0,%xmm0 ## constant 10001000
        shufps $221,%xmm1,%xmm1 ## constant 11011101

        movl   nb430nf_faction(%ebp),%edi
        ## move ix-iz to xmm4-xmm6 
        xorps   %xmm7,%xmm7

        movaps nb430nf_ix(%esp),%xmm4
        movaps nb430nf_iy(%esp),%xmm5
        movaps nb430nf_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        ## rsq in xmm4 

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb430nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb430nf_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r 
        movaps %xmm4,nb430nf_r(%esp)
        mulps nb430nf_gbscale(%esp),%xmm4

        cvttps2pi %xmm4,%mm6    ## mm6 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6

        movl nb430nf_GBtab(%ebp),%esi
        movd %mm6,%ecx
        psrlq $32,%mm6
        movd %mm6,%edx

        ## load coulomb table
        movaps (%esi,%ecx,4),%xmm4
        movaps (%esi,%edx,4),%xmm7
        ## transpose, using xmm3 for scratch
        movaps %xmm4,%xmm6
        unpcklps %xmm7,%xmm4    ## Y1 Y2 F1 F2 
        unpckhps %xmm7,%xmm6    ## G1 G2 H1 H2
        movhlps  %xmm4,%xmm5    ## F1 F2 
        movhlps  %xmm6,%xmm7    ## H1 H2
        ## coulomb table ready, in xmm4-xmm7    

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        movaps nb430nf_qq(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        addps  nb430nf_vctot(%esp),%xmm5
        movaps %xmm5,nb430nf_vctot(%esp)

        movaps nb430nf_r(%esp),%xmm4
        mulps nb430nf_tsc(%esp),%xmm4

        cvttps2pi %xmm4,%mm6
        cvtpi2ps %mm6,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $3,%mm6

        movl nb430nf_VFtab(%ebp),%esi
        movd %mm6,%ecx
        psrlq $32,%mm6
        movd %mm6,%edx

        ## dispersion 
        movaps (%esi,%ecx,4),%xmm4
        movaps (%esi,%edx,4),%xmm7
        ## transpose, using xmm3 for scratch
        movaps %xmm4,%xmm6
        unpcklps %xmm7,%xmm4    ## Y1 Y2 F1 F2 
        unpckhps %xmm7,%xmm6    ## G1 G2 H1 H2
        movhlps  %xmm4,%xmm5    ## F1 F2 
        movhlps  %xmm6,%xmm7    ## H1 H2
        ## dispersion table ready, in xmm4-xmm7         
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        mulps  nb430nf_c6(%esp),%xmm5    ## Vvdw6 
        addps  nb430nf_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb430nf_Vvdwtot(%esp)

        ## repulsion 
        movaps 16(%esi,%ecx,4),%xmm4
        movaps 16(%esi,%edx,4),%xmm7
        ## transpose, using xmm3 for scratch
        movaps %xmm4,%xmm6
        unpcklps %xmm7,%xmm4    ## Y1 Y2 F1 F2 
        unpckhps %xmm7,%xmm6    ## G1 G2 H1 H2
        movhlps  %xmm4,%xmm5    ## F1 F2 
        movhlps  %xmm6,%xmm7    ## H1 H2
        ## table ready, in xmm4-xmm7    
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        mulps  nb430nf_c12(%esp),%xmm5   ## Vvdw12 

        addps  nb430nf_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb430nf_Vvdwtot(%esp)
_nb_kernel430nf_ia32_sse.nb430nf_checksingle:   
        movl  nb430nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel430nf_ia32_sse.nb430nf_dosingle
        jmp    _nb_kernel430nf_ia32_sse.nb430nf_updateouterdata
_nb_kernel430nf_ia32_sse.nb430nf_dosingle: 
        movl nb430nf_charge(%ebp),%esi
        movl nb430nf_invsqrta(%ebp),%edx
        movl nb430nf_pos(%ebp),%edi
        movl  nb430nf_innerjjnr(%esp),%ecx
        movl  (%ecx),%eax
        xorps  %xmm2,%xmm2
        movaps %xmm2,%xmm6
        movss (%edx,%eax,4),%xmm2       ## isa2
        mulss nb430nf_isai(%esp),%xmm2
        movss %xmm2,nb430nf_isaprod(%esp)
        movss %xmm2,%xmm1
        mulss nb430nf_gbtsc(%esp),%xmm1
        movss %xmm1,nb430nf_gbscale(%esp)

        mulss  nb430nf_iq(%esp),%xmm2
        movss (%esi,%eax,4),%xmm6       ## xmm6(0) has the charge       
        mulss  %xmm2,%xmm6
        movss %xmm6,nb430nf_qq(%esp)

        movl nb430nf_type(%ebp),%esi
        movl %eax,%ecx
        movl (%esi,%ecx,4),%ecx
        movl nb430nf_vdwparam(%ebp),%esi
        shll %ecx
        addl nb430nf_ntia(%esp),%ecx
        movlps (%esi,%ecx,4),%xmm6
        movaps %xmm6,%xmm4
        shufps $252,%xmm4,%xmm4 ## constant 11111100    
        shufps $253,%xmm6,%xmm6 ## constant 11111101    

        movss %xmm4,nb430nf_c6(%esp)
        movss %xmm6,nb430nf_c12(%esp)

        leal  (%eax,%eax,2),%eax

        ## move coordinates to xmm0-xmm2 
        movss (%edi,%eax,4),%xmm0
        movss 4(%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2

        movss nb430nf_ix(%esp),%xmm4
        movss nb430nf_iy(%esp),%xmm5
        movss nb430nf_iz(%esp),%xmm6

        ## calc dr 
        subss %xmm0,%xmm4
        subss %xmm1,%xmm5
        subss %xmm2,%xmm6

        ## square it 
        mulss %xmm4,%xmm4
        mulss %xmm5,%xmm5
        mulss %xmm6,%xmm6
        addss %xmm5,%xmm4
        addss %xmm6,%xmm4
        ## rsq in xmm4 

        rsqrtss %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulss %xmm5,%xmm5
        movss nb430nf_three(%esp),%xmm1
        mulss %xmm4,%xmm5       ## rsq*lu*lu                    
        movss nb430nf_half(%esp),%xmm0
        subss %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulss %xmm2,%xmm1
        mulss %xmm1,%xmm0       ## xmm0=rinv 

        mulss %xmm0,%xmm4       ## xmm4=r 
        movaps %xmm4,nb430nf_r(%esp)
        mulss nb430nf_gbscale(%esp),%xmm4

        cvttss2si %xmm4,%ebx    ## mm6 contain lu indices 
        cvtsi2ss %ebx,%xmm6
        subss %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%ebx

        movl nb430nf_GBtab(%ebp),%esi

        movaps (%esi,%ebx,4),%xmm4
        movhlps %xmm4,%xmm6
        movaps %xmm4,%xmm5
        movaps %xmm6,%xmm7
        shufps $1,%xmm5,%xmm5
        shufps $1,%xmm7,%xmm7
        ## table ready in xmm4-xmm7 

        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp      
        movss nb430nf_qq(%esp),%xmm3
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 
        mulss  %xmm3,%xmm5 ## vcoul=qq*VV  
        addss  nb430nf_vctot(%esp),%xmm5
        movss %xmm5,nb430nf_vctot(%esp)

        movss nb430nf_r(%esp),%xmm4
        mulps nb430nf_tsc(%esp),%xmm4

        cvttss2si %xmm4,%ebx
        cvtsi2ss %ebx,%xmm6
        subss %xmm6,%xmm4
        movss %xmm4,%xmm1       ## xmm1=eps 
        movss %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=eps2 

        shll $3,%ebx
        movl nb430nf_VFtab(%ebp),%esi

        ## dispersion 
        movaps (%esi,%ebx,4),%xmm4
        movhlps %xmm4,%xmm6
        movaps %xmm4,%xmm5
        movaps %xmm6,%xmm7
        shufps $1,%xmm5,%xmm5
        shufps $1,%xmm7,%xmm7
        ## table ready in xmm4-xmm7 

        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp      
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 
        mulss  nb430nf_c6(%esp),%xmm5    ## Vvdw6
        addss  nb430nf_Vvdwtot(%esp),%xmm5
        movss %xmm5,nb430nf_Vvdwtot(%esp)

        ## repulsion 
        movaps 16(%esi,%ebx,4),%xmm4
        movhlps %xmm4,%xmm6
        movaps %xmm4,%xmm5
        movaps %xmm6,%xmm7
        shufps $1,%xmm5,%xmm5
        shufps $1,%xmm7,%xmm7
        ## table ready in xmm4-xmm7 

        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp      
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 

        mulss  nb430nf_c12(%esp),%xmm5   ## Vvdw12 

        addss  nb430nf_Vvdwtot(%esp),%xmm5
        movss %xmm5,nb430nf_Vvdwtot(%esp)

_nb_kernel430nf_ia32_sse.nb430nf_updateouterdata: 
        ## get n from stack
        movl nb430nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb430nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb430nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb430nf_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb430nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb430nf_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb430nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel430nf_ia32_sse.nb430nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb430nf_n(%esp)
        jmp _nb_kernel430nf_ia32_sse.nb430nf_outer
_nb_kernel430nf_ia32_sse.nb430nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb430nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel430nf_ia32_sse.nb430nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel430nf_ia32_sse.nb430nf_threadloop
_nb_kernel430nf_ia32_sse.nb430nf_end: 
        emms

        movl nb430nf_nouter(%esp),%eax
        movl nb430nf_ninner(%esp),%ebx
        movl nb430nf_outeriter(%ebp),%ecx
        movl nb430nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb430nf_salign(%esp),%eax
        addl %eax,%esp
        addl $324,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret





