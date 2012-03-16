##
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



.globl nb_kernel303_ia32_sse
.globl _nb_kernel303_ia32_sse
nb_kernel303_ia32_sse:  
_nb_kernel303_ia32_sse: 
.set nb303_p_nri, 8
.set nb303_iinr, 12
.set nb303_jindex, 16
.set nb303_jjnr, 20
.set nb303_shift, 24
.set nb303_shiftvec, 28
.set nb303_fshift, 32
.set nb303_gid, 36
.set nb303_pos, 40
.set nb303_faction, 44
.set nb303_charge, 48
.set nb303_p_facel, 52
.set nb303_argkrf, 56
.set nb303_argcrf, 60
.set nb303_Vc, 64
.set nb303_type, 68
.set nb303_p_ntype, 72
.set nb303_vdwparam, 76
.set nb303_Vvdw, 80
.set nb303_p_tabscale, 84
.set nb303_VFtab, 88
.set nb303_invsqrta, 92
.set nb303_dvda, 96
.set nb303_p_gbtabscale, 100
.set nb303_GBtab, 104
.set nb303_p_nthreads, 108
.set nb303_count, 112
.set nb303_mtx, 116
.set nb303_outeriter, 120
.set nb303_inneriter, 124
.set nb303_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb303_ixH1, 0
.set nb303_iyH1, 16
.set nb303_izH1, 32
.set nb303_ixH2, 48
.set nb303_iyH2, 64
.set nb303_izH2, 80
.set nb303_ixM, 96
.set nb303_iyM, 112
.set nb303_izM, 128
.set nb303_iqH, 144
.set nb303_iqM, 160
.set nb303_dxH1, 176
.set nb303_dyH1, 192
.set nb303_dzH1, 208
.set nb303_dxH2, 224
.set nb303_dyH2, 240
.set nb303_dzH2, 256
.set nb303_dxM, 272
.set nb303_dyM, 288
.set nb303_dzM, 304
.set nb303_qqH, 320
.set nb303_qqM, 336
.set nb303_rinvH1, 352
.set nb303_rinvH2, 368
.set nb303_rinvM, 384
.set nb303_rH1, 400
.set nb303_rH2, 416
.set nb303_rM, 432
.set nb303_tsc, 448
.set nb303_two, 464
.set nb303_vctot, 480
.set nb303_fixH1, 496
.set nb303_fiyH1, 512
.set nb303_fizH1, 528
.set nb303_fixH2, 544
.set nb303_fiyH2, 560
.set nb303_fizH2, 576
.set nb303_fixM, 592
.set nb303_fiyM, 608
.set nb303_fizM, 624
.set nb303_fjx, 640
.set nb303_fjy, 656
.set nb303_fjz, 672
.set nb303_half, 688
.set nb303_three, 704
.set nb303_is3, 720
.set nb303_ii3, 724
.set nb303_innerjjnr, 728
.set nb303_innerk, 732
.set nb303_n, 736
.set nb303_nn1, 740
.set nb303_nri, 744
.set nb303_nouter, 748
.set nb303_ninner, 752
.set nb303_salign, 756
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $760,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb303_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb303_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb303_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb303_nouter(%esp)
        movl %eax,nb303_ninner(%esp)


        movl nb303_p_tabscale(%ebp),%eax
        movss (%eax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb303_tsc(%esp)
        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb303_half(%esp)
        movss nb303_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## 2.0
        addps  %xmm2,%xmm3      ## 3.0
        movaps %xmm1,nb303_half(%esp)
        movaps %xmm2,nb303_two(%esp)
        movaps %xmm3,nb303_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb303_iinr(%ebp),%ecx             ## ecx = pointer into iinr[]    
        movl  (%ecx),%ebx               ## ebx =ii 

        movl  nb303_charge(%ebp),%edx
        movss 4(%edx,%ebx,4),%xmm3
        movss 12(%edx,%ebx,4),%xmm4
        movl nb303_p_facel(%ebp),%esi
        movss (%esi),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb303_iqH(%esp)
        movaps %xmm4,nb303_iqM(%esp)

_nb_kernel303_ia32_sse.nb303_threadloop: 
        movl  nb303_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel303_ia32_sse.nb303_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel303_ia32_sse.nb303_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb303_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb303_n(%esp)
        movl %ebx,nb303_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel303_ia32_sse.nb303_outerstart
        jmp _nb_kernel303_ia32_sse.nb303_end

_nb_kernel303_ia32_sse.nb303_outerstart: 
        ## ebx contains number of outer iterations
        addl nb303_nouter(%esp),%ebx
        movl %ebx,nb303_nouter(%esp)

_nb_kernel303_ia32_sse.nb303_outer: 
        movl  nb303_shift(%ebp),%eax            ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx        ## ebx=3*is 
        movl  %ebx,nb303_is3(%esp)      ## store is3 

        movl  nb303_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb303_iinr(%ebp),%ecx             ## ecx = pointer into iinr[]    
        movl  (%ecx,%esi,4),%ebx                ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb303_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb303_ii3(%esp)

        addss 12(%eax,%ebx,4),%xmm3
        addss 16(%eax,%ebx,4),%xmm4
        addss 20(%eax,%ebx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb303_ixH1(%esp)
        movaps %xmm4,nb303_iyH1(%esp)
        movaps %xmm5,nb303_izH1(%esp)

        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 24(%eax,%ebx,4),%xmm0
        addss 28(%eax,%ebx,4),%xmm1
        addss 32(%eax,%ebx,4),%xmm2
        addss 36(%eax,%ebx,4),%xmm3
        addss 40(%eax,%ebx,4),%xmm4
        addss 44(%eax,%ebx,4),%xmm5

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm0,nb303_ixH2(%esp)
        movaps %xmm1,nb303_iyH2(%esp)
        movaps %xmm2,nb303_izH2(%esp)
        movaps %xmm3,nb303_ixM(%esp)
        movaps %xmm4,nb303_iyM(%esp)
        movaps %xmm5,nb303_izM(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb303_vctot(%esp)
        movaps %xmm4,nb303_fixH1(%esp)
        movaps %xmm4,nb303_fiyH1(%esp)
        movaps %xmm4,nb303_fizH1(%esp)
        movaps %xmm4,nb303_fixH2(%esp)
        movaps %xmm4,nb303_fiyH2(%esp)
        movaps %xmm4,nb303_fizH2(%esp)
        movaps %xmm4,nb303_fixM(%esp)
        movaps %xmm4,nb303_fiyM(%esp)
        movaps %xmm4,nb303_fizM(%esp)

        movl  nb303_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx                ## jindex[n] 
        movl  4(%eax,%esi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movl  nb303_pos(%ebp),%esi
        movl  nb303_faction(%ebp),%edi
        movl  nb303_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb303_innerjjnr(%esp)        ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb303_ninner(%esp),%ecx
        movl  %ecx,nb303_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb303_innerk(%esp)   ## number of innerloop atoms 
        jge   _nb_kernel303_ia32_sse.nb303_unroll_loop
        jmp   _nb_kernel303_ia32_sse.nb303_odd_inner
_nb_kernel303_ia32_sse.nb303_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb303_innerjjnr(%esp),%edx        ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx             ## eax-edx=jnr1-4 

        addl $16,nb303_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb303_charge(%ebp),%esi    ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4              ## and in xmm4 
        mulps  nb303_iqH(%esp),%xmm3
        mulps  nb303_iqM(%esp),%xmm4

        movaps  %xmm3,nb303_qqH(%esp)
        movaps  %xmm4,nb303_qqM(%esp)

        movl nb303_pos(%ebp),%esi       ## base of pos[] 

        leal  (%eax,%eax,2),%eax        ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx        ## replace jnr with j3 
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

        shufps $136,%xmm6,%xmm2 ## 10001000

        shufps $136,%xmm5,%xmm0 ## 10001000
        shufps $221,%xmm5,%xmm1 ## 11011101             

        ## move ixH1-izH1 to xmm4-xmm6 
        movaps nb303_ixH1(%esp),%xmm4
        movaps nb303_iyH1(%esp),%xmm5
        movaps nb303_izH1(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb303_dxH1(%esp)
        movaps %xmm5,nb303_dyH1(%esp)
        movaps %xmm6,nb303_dzH1(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        movaps %xmm4,%xmm7
        ## rsqH1 in xmm7 

        ## move ixH2-izH2 to xmm4-xmm6 
        movaps nb303_ixH2(%esp),%xmm4
        movaps nb303_iyH2(%esp),%xmm5
        movaps nb303_izH2(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb303_dxH2(%esp)
        movaps %xmm5,nb303_dyH2(%esp)
        movaps %xmm6,nb303_dzH2(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm6
        addps %xmm4,%xmm6
        ## rsqH2 in xmm6 

        ## move ixM-izM to xmm3-xmm5  
        movaps nb303_ixM(%esp),%xmm3
        movaps nb303_iyM(%esp),%xmm4
        movaps nb303_izM(%esp),%xmm5

        ## calc dr 
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        ## store dr 
        movaps %xmm3,nb303_dxM(%esp)
        movaps %xmm4,nb303_dyM(%esp)
        movaps %xmm5,nb303_dzM(%esp)
        ## square it 
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm4,%xmm5
        addps %xmm3,%xmm5
        ## rsqM in xmm5, rsqH2 in xmm6, rsqH1 in xmm7 

        ## start with rsqH1 - seed to xmm2      
        rsqrtps %xmm7,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb303_three(%esp),%xmm4
        mulps   %xmm7,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb303_half(%esp),%xmm4
        movaps  %xmm4,nb303_rinvH1(%esp)        ## rinvH1 in xmm4 
        mulps   %xmm4,%xmm7
        movaps  %xmm7,nb303_rH1(%esp)

        ## rsqH2 - seed in xmm2 
        rsqrtps %xmm6,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb303_three(%esp),%xmm4
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb303_half(%esp),%xmm4
        movaps  %xmm4,nb303_rinvH2(%esp)        ## rinvH2 in xmm4 
        mulps   %xmm4,%xmm6
        movaps  %xmm6,nb303_rH2(%esp)

        ## rsqM - seed to xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb303_three(%esp),%xmm4
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb303_half(%esp),%xmm4
        movaps  %xmm4,nb303_rinvM(%esp)         ## rinvM in xmm4 
        mulps   %xmm4,%xmm5
        movaps  %xmm5,nb303_rM(%esp)

        ## do H1 interactions 
        ## rH1 is still in xmm7 
        mulps   nb303_tsc(%esp),%xmm7
        movhlps %xmm7,%xmm4
        cvttps2pi %xmm7,%mm6
        cvttps2pi %xmm4,%mm7    ## mm6/mm7 contain lu indices 

        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm4
        movlhps %xmm4,%xmm3

        subps %xmm3,%xmm7

        movaps %xmm7,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $2,%mm6
        pslld $2,%mm7

        movd %eax,%mm0
        movd %ebx,%mm1
        movd %ecx,%mm2
        movd %edx,%mm3

        movl nb303_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 ## 11011101
        ## coulomb table ready, in xmm4-xmm7      

                mulps  %xmm1,%xmm6      ## xmm6=Geps 
                mulps  %xmm2,%xmm7      ## xmm7=Heps2 
                addps  %xmm6,%xmm5
                addps  %xmm7,%xmm5      ## xmm5=Fp        
                mulps  nb303_two(%esp),%xmm7            ## two*Heps2 
                movaps nb303_qqH(%esp),%xmm0
                addps  %xmm6,%xmm7
                addps  %xmm5,%xmm7 ## xmm7=FF 
                mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
                addps  %xmm4,%xmm5 ## xmm5=VV 
                mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
                mulps  %xmm7,%xmm0 ## fijC=FF*qq 

                ## at this point mm5 contains vcoul and xmm0 fijC 
                ## increment vcoul - then we can get rid of mm5 
                addps  nb303_vctot(%esp),%xmm5
                movaps %xmm5,nb303_vctot(%esp)
                xorps  %xmm4,%xmm4

                mulps  nb303_tsc(%esp),%xmm0
                mulps  nb303_rinvH1(%esp),%xmm0
                subps  %xmm0,%xmm4

                movaps nb303_dxH1(%esp),%xmm0
                movaps nb303_dyH1(%esp),%xmm1
                movaps nb303_dzH1(%esp),%xmm2
                mulps  %xmm4,%xmm0
                mulps  %xmm4,%xmm1
                mulps  %xmm4,%xmm2      ## tx in xmm0-xmm2 

                ## update H1 forces 
                movaps nb303_fixH1(%esp),%xmm3
                movaps nb303_fiyH1(%esp),%xmm4
                movaps nb303_fizH1(%esp),%xmm7
                addps  %xmm0,%xmm3
                addps  %xmm1,%xmm4
                addps  %xmm2,%xmm7
                movaps %xmm3,nb303_fixH1(%esp)
                movaps %xmm4,nb303_fiyH1(%esp)
                movaps %xmm7,nb303_fizH1(%esp)
                ## update j forces with water H1 
                movaps %xmm0,nb303_fjx(%esp)
                movaps %xmm1,nb303_fjy(%esp)
                movaps %xmm2,nb303_fjz(%esp)

                ## Done with H1 interactions - now H2! 
                movaps nb303_rH2(%esp),%xmm7
                mulps   nb303_tsc(%esp),%xmm7
                movhlps %xmm7,%xmm4
                cvttps2pi %xmm7,%mm6
                cvttps2pi %xmm4,%mm7    ## mm6/mm7 contain lu indices 

                cvtpi2ps %mm6,%xmm3
                cvtpi2ps %mm7,%xmm4
                movlhps %xmm4,%xmm3

                subps %xmm3,%xmm7
                movaps %xmm7,%xmm1      ## xmm1=eps 
                movaps %xmm1,%xmm2
                mulps  %xmm2,%xmm2      ## xmm2=eps2 
                pslld $2,%mm6
                pslld $2,%mm7

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx


        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 ## 11011101
        ## coulomb table ready, in xmm4-xmm7      

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp        
        mulps  nb303_two(%esp),%xmm7            ## two*Heps2 
        movaps nb303_qqH(%esp),%xmm0
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm0,%xmm7 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm7 fijC 
        ## increment vcoul 
        xorps  %xmm4,%xmm4
        addps  nb303_vctot(%esp),%xmm5
        mulps  nb303_rinvH2(%esp),%xmm7
        movaps %xmm5,nb303_vctot(%esp)
        mulps  nb303_tsc(%esp),%xmm7
        subps %xmm7,%xmm4

        movaps nb303_dxH2(%esp),%xmm0
        movaps nb303_dyH2(%esp),%xmm1
        movaps nb303_dzH2(%esp),%xmm2
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2

        ## update H2 forces 
        movaps nb303_fixH2(%esp),%xmm3
        movaps nb303_fiyH2(%esp),%xmm4
        movaps nb303_fizH2(%esp),%xmm7
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm7
        movaps %xmm3,nb303_fixH2(%esp)
        movaps %xmm4,nb303_fiyH2(%esp)
        movaps %xmm7,nb303_fizH2(%esp)
        ## update j forces with water H2 
        addps  nb303_fjx(%esp),%xmm0
        addps  nb303_fjy(%esp),%xmm1
        addps  nb303_fjz(%esp),%xmm2
        movaps %xmm0,nb303_fjx(%esp)
        movaps %xmm1,nb303_fjy(%esp)
        movaps %xmm2,nb303_fjz(%esp)

        ## Done with H2, finally we do M interactions 
        movaps nb303_rM(%esp),%xmm7
        mulps   nb303_tsc(%esp),%xmm7
        movhlps %xmm7,%xmm4
        cvttps2pi %xmm7,%mm6
        cvttps2pi %xmm4,%mm7    ## mm6/mm7 contain lu indices 

        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm4
        movlhps %xmm4,%xmm3

        subps %xmm3,%xmm7
        movaps %xmm7,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $2,%mm6
        pslld $2,%mm7

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 ## 11011101
        ## coulomb table ready, in xmm4-xmm7      

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp        
        mulps  nb303_two(%esp),%xmm7            ## two*Heps2 
        movaps nb303_qqM(%esp),%xmm0
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm0,%xmm7 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm0 fijC 
        ## increment vcoul 
        xorps  %xmm4,%xmm4
        addps  nb303_vctot(%esp),%xmm5
        mulps  nb303_rinvM(%esp),%xmm7
        movaps %xmm5,nb303_vctot(%esp)
        mulps  nb303_tsc(%esp),%xmm7
        subps  %xmm7,%xmm4

        movaps nb303_dxM(%esp),%xmm0
        movaps nb303_dyM(%esp),%xmm1
        movaps nb303_dzM(%esp),%xmm2
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2

        movd %mm0,%eax
        movd %mm1,%ebx
        movd %mm2,%ecx
        movd %mm3,%edx

        ## update M forces 
        movaps nb303_fixM(%esp),%xmm3
        movaps nb303_fiyM(%esp),%xmm4
        movaps nb303_fizM(%esp),%xmm7
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm7
        movaps %xmm3,nb303_fixM(%esp)
        movaps %xmm4,nb303_fiyM(%esp)
        movaps %xmm7,nb303_fizM(%esp)

        movl nb303_faction(%ebp),%edi

        ## update j forces 
        addps nb303_fjx(%esp),%xmm0
        addps nb303_fjy(%esp),%xmm1
        addps nb303_fjz(%esp),%xmm2

        movlps (%edi,%eax,4),%xmm4
        movlps (%edi,%ecx,4),%xmm7
        movhps (%edi,%ebx,4),%xmm4
        movhps (%edi,%edx,4),%xmm7

        movaps %xmm4,%xmm3
        shufps $136,%xmm7,%xmm3 ## 10001000
        shufps $221,%xmm7,%xmm4 ## 11011101                           
        ## xmm3 has fjx, xmm4 has fjy 
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        ## unpack them back for storing 
        movaps %xmm3,%xmm7
        unpcklps %xmm4,%xmm7
        unpckhps %xmm4,%xmm3
        movlps %xmm7,(%edi,%eax,4)
        movlps %xmm3,(%edi,%ecx,4)
        movhps %xmm7,(%edi,%ebx,4)
        movhps %xmm3,(%edi,%edx,4)
        ## finally z forces 
        movss  8(%edi,%eax,4),%xmm0
        movss  8(%edi,%ebx,4),%xmm1
        movss  8(%edi,%ecx,4),%xmm3
        movss  8(%edi,%edx,4),%xmm4
        subss  %xmm2,%xmm0
        shufps $229,%xmm2,%xmm2 ## 11100101
        subss  %xmm2,%xmm1
        shufps $234,%xmm2,%xmm2 ## 11101010
        subss  %xmm2,%xmm3
        shufps $255,%xmm2,%xmm2 ## 11111111
        subss  %xmm2,%xmm4
        movss  %xmm0,8(%edi,%eax,4)
        movss  %xmm1,8(%edi,%ebx,4)
        movss  %xmm3,8(%edi,%ecx,4)
        movss  %xmm4,8(%edi,%edx,4)

        ## should we do one more iteration? 
        subl $4,nb303_innerk(%esp)
        jl    _nb_kernel303_ia32_sse.nb303_odd_inner
        jmp   _nb_kernel303_ia32_sse.nb303_unroll_loop
_nb_kernel303_ia32_sse.nb303_odd_inner: 
        addl $4,nb303_innerk(%esp)
        jnz   _nb_kernel303_ia32_sse.nb303_odd_loop
        jmp   _nb_kernel303_ia32_sse.nb303_updateouterdata
_nb_kernel303_ia32_sse.nb303_odd_loop: 
        movl  nb303_innerjjnr(%esp),%edx        ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb303_innerjjnr(%esp)

        xorps %xmm4,%xmm4
        movss nb303_iqM(%esp),%xmm4
        movl nb303_charge(%ebp),%esi
        movhps nb303_iqH(%esp),%xmm4
        movss (%esi,%eax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb303_qqM(%esp)    ## use dummy qq for storage 

        movl nb303_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## move j coords to xmm0-xmm2 
        movss (%esi,%eax,4),%xmm0
        movss 4(%esi,%eax,4),%xmm1
        movss 8(%esi,%eax,4),%xmm2
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movss nb303_ixM(%esp),%xmm3
        movss nb303_iyM(%esp),%xmm4
        movss nb303_izM(%esp),%xmm5

        movlps nb303_ixH1(%esp),%xmm6
        movlps nb303_ixH2(%esp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm3
        movlps nb303_iyH1(%esp),%xmm6
        movlps nb303_iyH2(%esp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm4
        movlps nb303_izH1(%esp),%xmm6
        movlps nb303_izH2(%esp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm5

        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        ## use dummy dx for storage
        movaps %xmm3,nb303_dxM(%esp)
        movaps %xmm4,nb303_dyM(%esp)
        movaps %xmm5,nb303_dzM(%esp)

        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5

        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        ## rsq in xmm4 

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb303_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb303_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        ## a little trick to avoid NaNs: 
        ## positions 0,2,and 3 are valid, but not 1. 
        ## If it contains NaN it doesnt help to mult by 0, 
        ## So we shuffle it and copy pos 0 to pos1! 
        shufps $224,%xmm0,%xmm0 ## 11100000      

        mulps %xmm0,%xmm4       ## xmm4=r 
        movaps %xmm0,nb303_rinvM(%esp)

        mulps nb303_tsc(%esp),%xmm4
        movhlps %xmm4,%xmm7
        cvttps2pi %xmm4,%mm6
        cvttps2pi %xmm7,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm7
        movlhps %xmm7,%xmm3

        subps   %xmm3,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $2,%mm6
        pslld $2,%mm7

        movd %eax,%mm0
        movd %ecx,%mm1
        movd %edx,%mm2

        movl nb303_VFtab(%ebp),%esi
        movd %mm6,%eax
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm7,%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 ## 11011101
        ## coulomb table ready, in xmm4-xmm7      
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp        
        mulps  nb303_two(%esp),%xmm7            ## two*Heps2 
        movaps nb303_qqM(%esp),%xmm0
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm0 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm0 fijC 
        ## increment vcoul - then we can get rid of mm5 
        addps  nb303_vctot(%esp),%xmm5
        movaps %xmm5,nb303_vctot(%esp)

        xorps %xmm4,%xmm4
        mulps  nb303_tsc(%esp),%xmm0
        mulps  nb303_rinvM(%esp),%xmm0
        subps  %xmm0,%xmm4

        movd %mm0,%eax
        movd %mm1,%ecx
        movd %mm2,%edx

        movaps nb303_dxM(%esp),%xmm0
        movaps nb303_dyM(%esp),%xmm1
        movaps nb303_dzM(%esp),%xmm2

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2 ## xmm0-xmm2 now contains tx-tz (partial force) 
        movss  nb303_fixM(%esp),%xmm3
        movss  nb303_fiyM(%esp),%xmm4
        movss  nb303_fizM(%esp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,nb303_fixM(%esp)
        movss  %xmm4,nb303_fiyM(%esp)
        movss  %xmm5,nb303_fizM(%esp)   ## updated the M force now do the H's 
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        shufps $230,%xmm3,%xmm3 ## 11100110      ;# shift right 
        shufps $230,%xmm4,%xmm4 ## 11100110
        shufps $230,%xmm5,%xmm5 ## 11100110
        addss  nb303_fixH1(%esp),%xmm3
        addss  nb303_fiyH1(%esp),%xmm4
        addss  nb303_fizH1(%esp),%xmm5
        movss  %xmm3,nb303_fixH1(%esp)
        movss  %xmm4,nb303_fiyH1(%esp)
        movss  %xmm5,nb303_fizH1(%esp)          ## updated the H1 force 

        movl nb303_faction(%ebp),%edi
        shufps $231,%xmm3,%xmm3 ## 11100111      ;# shift right 
        shufps $231,%xmm4,%xmm4 ## 11100111
        shufps $231,%xmm5,%xmm5 ## 11100111
        addss  nb303_fixH2(%esp),%xmm3
        addss  nb303_fiyH2(%esp),%xmm4
        addss  nb303_fizH2(%esp),%xmm5
        movss  %xmm3,nb303_fixH2(%esp)
        movss  %xmm4,nb303_fiyH2(%esp)
        movss  %xmm5,nb303_fizH2(%esp)          ## updated the H2 force 

        ## the fj's - start by accumulating the tx/ty/tz force in xmm0, xmm1 
        xorps  %xmm5,%xmm5
        movaps %xmm0,%xmm3
        movlps (%edi,%eax,4),%xmm6
        movss  8(%edi,%eax,4),%xmm7
        unpcklps %xmm1,%xmm3
        movlhps  %xmm5,%xmm3
        unpckhps %xmm1,%xmm0
        addps    %xmm3,%xmm0
        movhlps  %xmm0,%xmm3
        addps    %xmm3,%xmm0    ## x,y sum in xmm0 

        movhlps  %xmm2,%xmm1
        addss    %xmm1,%xmm2
        shufps  $1,%xmm1,%xmm1
        addss    %xmm1,%xmm2    ## z sum in xmm2 
        subps    %xmm0,%xmm6
        subss    %xmm2,%xmm7

        movlps %xmm6,(%edi,%eax,4)
        movss  %xmm7,8(%edi,%eax,4)

        decl nb303_innerk(%esp)
        jz    _nb_kernel303_ia32_sse.nb303_updateouterdata
        jmp   _nb_kernel303_ia32_sse.nb303_odd_loop
_nb_kernel303_ia32_sse.nb303_updateouterdata: 
        movl  nb303_ii3(%esp),%ecx
        movl  nb303_faction(%ebp),%edi
        movl  nb303_fshift(%ebp),%esi
        movl  nb303_is3(%esp),%edx

        ## accumulate  H1 i forces in xmm0, xmm1, xmm2 
        movaps nb303_fixH1(%esp),%xmm0
        movaps nb303_fiyH1(%esp),%xmm1
        movaps nb303_fizH1(%esp),%xmm2

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
        movss  12(%edi,%ecx,4),%xmm3
        movss  16(%edi,%ecx,4),%xmm4
        movss  20(%edi,%ecx,4),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,12(%edi,%ecx,4)
        movss  %xmm4,16(%edi,%ecx,4)
        movss  %xmm5,20(%edi,%ecx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        movaps %xmm0,%xmm6
        movss %xmm2,%xmm7
        movlhps %xmm1,%xmm6
        shufps $8,%xmm6,%xmm6 ## constant 00001000      

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movaps nb303_fixH2(%esp),%xmm0
        movaps nb303_fiyH2(%esp),%xmm1
        movaps nb303_fizH2(%esp),%xmm2

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
        movss  24(%edi,%ecx,4),%xmm3
        movss  28(%edi,%ecx,4),%xmm4
        movss  32(%edi,%ecx,4),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,24(%edi,%ecx,4)
        movss  %xmm4,28(%edi,%ecx,4)
        movss  %xmm5,32(%edi,%ecx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        addss %xmm2,%xmm7
        movlhps %xmm1,%xmm0
        shufps $8,%xmm0,%xmm0 ## constant 00001000      
        addps   %xmm0,%xmm6

        ## accumulate M i forces in xmm0, xmm1, xmm2 
        movaps nb303_fixM(%esp),%xmm0
        movaps nb303_fiyM(%esp),%xmm1
        movaps nb303_fizM(%esp),%xmm2

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
        movss  36(%edi,%ecx,4),%xmm3
        movss  40(%edi,%ecx,4),%xmm4
        movss  44(%edi,%ecx,4),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,36(%edi,%ecx,4)
        movss  %xmm4,40(%edi,%ecx,4)
        movss  %xmm5,44(%edi,%ecx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        addss %xmm2,%xmm7
        movlhps %xmm1,%xmm0
        shufps $8,%xmm0,%xmm0 ## constant 00001000      
        addps   %xmm0,%xmm6

        ## increment fshift force  
        movlps  (%esi,%edx,4),%xmm3
        movss  8(%esi,%edx,4),%xmm4
        addps  %xmm6,%xmm3
        addss  %xmm7,%xmm4
        movlps  %xmm3,(%esi,%edx,4)
        movss  %xmm4,8(%esi,%edx,4)

        ## get n from stack
        movl nb303_n(%esp),%esi
        ## get group index for i particle 
        movl  nb303_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb303_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb303_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb303_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel303_ia32_sse.nb303_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb303_n(%esp)
        jmp _nb_kernel303_ia32_sse.nb303_outer
_nb_kernel303_ia32_sse.nb303_outerend: 
        ## check if more outer neighborlists remain
        movl  nb303_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel303_ia32_sse.nb303_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel303_ia32_sse.nb303_threadloop
_nb_kernel303_ia32_sse.nb303_end: 
        emms

        movl nb303_nouter(%esp),%eax
        movl nb303_ninner(%esp),%ebx
        movl nb303_outeriter(%ebp),%ecx
        movl nb303_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb303_salign(%esp),%eax
        addl %eax,%esp
        addl $760,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret



.globl nb_kernel303nf_ia32_sse
.globl _nb_kernel303nf_ia32_sse
nb_kernel303nf_ia32_sse:        
_nb_kernel303nf_ia32_sse:       
.set nb303nf_p_nri, 8
.set nb303nf_iinr, 12
.set nb303nf_jindex, 16
.set nb303nf_jjnr, 20
.set nb303nf_shift, 24
.set nb303nf_shiftvec, 28
.set nb303nf_fshift, 32
.set nb303nf_gid, 36
.set nb303nf_pos, 40
.set nb303nf_faction, 44
.set nb303nf_charge, 48
.set nb303nf_p_facel, 52
.set nb303nf_argkrf, 56
.set nb303nf_argcrf, 60
.set nb303nf_Vc, 64
.set nb303nf_type, 68
.set nb303nf_p_ntype, 72
.set nb303nf_vdwparam, 76
.set nb303nf_Vvdw, 80
.set nb303nf_p_tabscale, 84
.set nb303nf_VFtab, 88
.set nb303nf_invsqrta, 92
.set nb303nf_dvda, 96
.set nb303nf_p_gbtabscale, 100
.set nb303nf_GBtab, 104
.set nb303nf_p_nthreads, 108
.set nb303nf_count, 112
.set nb303nf_mtx, 116
.set nb303nf_outeriter, 120
.set nb303nf_inneriter, 124
.set nb303nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb303nf_ixH1, 0
.set nb303nf_iyH1, 16
.set nb303nf_izH1, 32
.set nb303nf_ixH2, 48
.set nb303nf_iyH2, 64
.set nb303nf_izH2, 80
.set nb303nf_ixM, 96
.set nb303nf_iyM, 112
.set nb303nf_izM, 128
.set nb303nf_iqH, 144
.set nb303nf_iqM, 160
.set nb303nf_qqH, 176
.set nb303nf_qqM, 192
.set nb303nf_rinvH1, 208
.set nb303nf_rinvH2, 224
.set nb303nf_rinvM, 240
.set nb303nf_rH1, 256
.set nb303nf_rH2, 272
.set nb303nf_rM, 288
.set nb303nf_tsc, 304
.set nb303nf_vctot, 320
.set nb303nf_half, 336
.set nb303nf_three, 352
.set nb303nf_is3, 368
.set nb303nf_ii3, 372
.set nb303nf_innerjjnr, 376
.set nb303nf_innerk, 380
.set nb303nf_n, 384
.set nb303nf_nn1, 388
.set nb303nf_nri, 392
.set nb303nf_nouter, 396
.set nb303nf_ninner, 400
.set nb303nf_salign, 404
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $408,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb303nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb303nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb303nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb303nf_nouter(%esp)
        movl %eax,nb303nf_ninner(%esp)


        movl nb303nf_p_tabscale(%ebp),%eax
        movss (%eax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb303nf_tsc(%esp)
        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb303nf_half(%esp)
        movss nb303nf_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## 2.0
        addps  %xmm2,%xmm3      ## 3.0
        movaps %xmm1,nb303nf_half(%esp)
        movaps %xmm3,nb303nf_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb303nf_iinr(%ebp),%ecx           ## ecx = pointer into iinr[]    
        movl  (%ecx),%ebx               ## ebx =ii 

        movl  nb303nf_charge(%ebp),%edx
        movss 4(%edx,%ebx,4),%xmm3
        movss 12(%edx,%ebx,4),%xmm4
        movl nb303nf_p_facel(%ebp),%esi
        movss (%esi),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb303nf_iqH(%esp)
        movaps %xmm4,nb303nf_iqM(%esp)

_nb_kernel303nf_ia32_sse.nb303nf_threadloop: 
        movl  nb303nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel303nf_ia32_sse.nb303nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel303nf_ia32_sse.nb303nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb303nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb303nf_n(%esp)
        movl %ebx,nb303nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel303nf_ia32_sse.nb303nf_outerstart
        jmp _nb_kernel303nf_ia32_sse.nb303nf_end

_nb_kernel303nf_ia32_sse.nb303nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb303nf_nouter(%esp),%ebx
        movl %ebx,nb303nf_nouter(%esp)

_nb_kernel303nf_ia32_sse.nb303nf_outer: 
        movl  nb303nf_shift(%ebp),%eax          ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx        ## ebx=3*is 
        movl  %ebx,nb303nf_is3(%esp)            ## store is3 

        movl  nb303nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb303nf_iinr(%ebp),%ecx           ## ecx = pointer into iinr[]    
        movl  (%ecx,%esi,4),%ebx                ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb303nf_pos(%ebp),%eax    ## eax = base of pos[]  
        movl  %ebx,nb303nf_ii3(%esp)

        addss 12(%eax,%ebx,4),%xmm3
        addss 16(%eax,%ebx,4),%xmm4
        addss 20(%eax,%ebx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb303nf_ixH1(%esp)
        movaps %xmm4,nb303nf_iyH1(%esp)
        movaps %xmm5,nb303nf_izH1(%esp)

        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 24(%eax,%ebx,4),%xmm0
        addss 28(%eax,%ebx,4),%xmm1
        addss 32(%eax,%ebx,4),%xmm2
        addss 36(%eax,%ebx,4),%xmm3
        addss 40(%eax,%ebx,4),%xmm4
        addss 44(%eax,%ebx,4),%xmm5

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm0,nb303nf_ixH2(%esp)
        movaps %xmm1,nb303nf_iyH2(%esp)
        movaps %xmm2,nb303nf_izH2(%esp)
        movaps %xmm3,nb303nf_ixM(%esp)
        movaps %xmm4,nb303nf_iyM(%esp)
        movaps %xmm5,nb303nf_izM(%esp)

        ## clear vctot 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb303nf_vctot(%esp)

        movl  nb303nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx                ## jindex[n] 
        movl  4(%eax,%esi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movl  nb303nf_pos(%ebp),%esi
        movl  nb303nf_faction(%ebp),%edi
        movl  nb303nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb303nf_innerjjnr(%esp)      ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb303nf_ninner(%esp),%ecx
        movl  %ecx,nb303nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb303nf_innerk(%esp)         ## number of innerloop atoms 
        jge   _nb_kernel303nf_ia32_sse.nb303nf_unroll_loop
        jmp   _nb_kernel303nf_ia32_sse.nb303nf_odd_inner
_nb_kernel303nf_ia32_sse.nb303nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb303nf_innerjjnr(%esp),%edx      ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx             ## eax-edx=jnr1-4 

        addl $16,nb303nf_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb303nf_charge(%ebp),%esi  ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4              ## and in xmm4 
        mulps  nb303nf_iqH(%esp),%xmm3
        mulps  nb303nf_iqM(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps  %xmm3,nb303nf_qqH(%esp)
        movaps  %xmm4,nb303nf_qqM(%esp)

        movl nb303nf_pos(%ebp),%esi     ## base of pos[] 

        leal  (%eax,%eax,2),%eax        ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx        ## replace jnr with j3 
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

        shufps $136,%xmm6,%xmm2 ## 10001000

        shufps $136,%xmm5,%xmm0 ## 10001000
        shufps $221,%xmm5,%xmm1 ## 11011101             

        ## move ixH1-izH1 to xmm4-xmm6 
        movaps nb303nf_ixH1(%esp),%xmm4
        movaps nb303nf_iyH1(%esp),%xmm5
        movaps nb303nf_izH1(%esp),%xmm6

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
        movaps %xmm4,%xmm7
        ## rsqH1 in xmm7 

        ## move ixH2-izH2 to xmm4-xmm6 
        movaps nb303nf_ixH2(%esp),%xmm4
        movaps nb303nf_iyH2(%esp),%xmm5
        movaps nb303nf_izH2(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm6
        addps %xmm4,%xmm6
        ## rsqH2 in xmm6 

        ## move ixM-izM to xmm3-xmm5  
        movaps nb303nf_ixM(%esp),%xmm3
        movaps nb303nf_iyM(%esp),%xmm4
        movaps nb303nf_izM(%esp),%xmm5

        ## calc dr 
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        ## square it 
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm4,%xmm5
        addps %xmm3,%xmm5
        ## rsqM in xmm5, rsqH2 in xmm6, rsqH1 in xmm7 

        ## start with rsqH1 - seed to xmm2      
        rsqrtps %xmm7,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb303nf_three(%esp),%xmm4
        mulps   %xmm7,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb303nf_half(%esp),%xmm4
        movaps  %xmm4,nb303nf_rinvH1(%esp)      ## rinvH1 in xmm4 
        mulps   %xmm4,%xmm7
        movaps  %xmm7,nb303nf_rH1(%esp)

        ## rsqH2 - seed in xmm2 
        rsqrtps %xmm6,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb303nf_three(%esp),%xmm4
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb303nf_half(%esp),%xmm4
        movaps  %xmm4,nb303nf_rinvH2(%esp)      ## rinvH2 in xmm4 
        mulps   %xmm4,%xmm6
        movaps  %xmm6,nb303nf_rH2(%esp)

        ## rsqM - seed to xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb303nf_three(%esp),%xmm4
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb303nf_half(%esp),%xmm4
        movaps  %xmm4,nb303nf_rinvM(%esp)       ## rinvM in xmm4 
        mulps   %xmm4,%xmm5
        movaps  %xmm5,nb303nf_rM(%esp)

        ## do H1 interactions 
        ## rH1 is still in xmm7 
        mulps   nb303nf_tsc(%esp),%xmm7
        movhlps %xmm7,%xmm4
        cvttps2pi %xmm7,%mm6
        cvttps2pi %xmm4,%mm7    ## mm6/mm7 contain lu indices 

        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm4
        movlhps %xmm4,%xmm3

        subps %xmm3,%xmm7

        movaps %xmm7,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $2,%mm6
        pslld $2,%mm7

        movd %eax,%mm0
        movd %ebx,%mm1
        movd %ecx,%mm2
        movd %edx,%mm3

        movl nb303nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 ## 11011101
        ## coulomb table ready, in xmm4-xmm7      

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp        
        movaps nb303nf_qqH(%esp),%xmm0
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm0,%xmm5 ## vcoul=qq*VV  

        ## at this point mm5 contains vcoul 
        ## increment vcoul 
        addps  nb303nf_vctot(%esp),%xmm5
        movaps %xmm5,nb303nf_vctot(%esp)

        ## Done with H1 interactions - now H2! 
        movaps nb303nf_rH2(%esp),%xmm7
        mulps   nb303nf_tsc(%esp),%xmm7
        movhlps %xmm7,%xmm4
        cvttps2pi %xmm7,%mm6
        cvttps2pi %xmm4,%mm7    ## mm6/mm7 contain lu indices 

        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm4
        movlhps %xmm4,%xmm3

        subps %xmm3,%xmm7
        movaps %xmm7,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $2,%mm6
        pslld $2,%mm7

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 ## 11011101
        ## coulomb table ready, in xmm4-xmm7      

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp        
        movaps nb303nf_qqH(%esp),%xmm0
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        ## increment vcoul 
        addps  nb303nf_vctot(%esp),%xmm5
        movaps %xmm5,nb303nf_vctot(%esp)

        ## Done with H2, finally we do M interactions 
        movaps nb303nf_rM(%esp),%xmm7
        mulps   nb303nf_tsc(%esp),%xmm7
        movhlps %xmm7,%xmm4
        cvttps2pi %xmm7,%mm6
        cvttps2pi %xmm4,%mm7    ## mm6/mm7 contain lu indices 

        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm4
        movlhps %xmm4,%xmm3

        subps %xmm3,%xmm7
        movaps %xmm7,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $2,%mm6
        pslld $2,%mm7

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 ## 11011101
        ## coulomb table ready, in xmm4-xmm7      

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp        
        movaps nb303nf_qqM(%esp),%xmm0
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        ## increment vcoul 
        addps  nb303nf_vctot(%esp),%xmm5
        movaps %xmm5,nb303nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $4,nb303nf_innerk(%esp)
        jl    _nb_kernel303nf_ia32_sse.nb303nf_odd_inner
        jmp   _nb_kernel303nf_ia32_sse.nb303nf_unroll_loop
_nb_kernel303nf_ia32_sse.nb303nf_odd_inner: 
        addl $4,nb303nf_innerk(%esp)
        jnz   _nb_kernel303nf_ia32_sse.nb303nf_odd_loop
        jmp   _nb_kernel303nf_ia32_sse.nb303nf_updateouterdata
_nb_kernel303nf_ia32_sse.nb303nf_odd_loop: 
        movl  nb303nf_innerjjnr(%esp),%edx      ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb303nf_innerjjnr(%esp)

        xorps %xmm4,%xmm4
        movss nb303nf_iqM(%esp),%xmm4
        movl nb303nf_charge(%ebp),%esi
        movhps nb303nf_iqH(%esp),%xmm4
        movss (%esi,%eax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb303nf_qqM(%esp)          ## use dummy qq for storage 

        movl nb303nf_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## move j coords to xmm0-xmm2 
        movss (%esi,%eax,4),%xmm0
        movss 4(%esi,%eax,4),%xmm1
        movss 8(%esi,%eax,4),%xmm2
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movss nb303nf_ixM(%esp),%xmm3
        movss nb303nf_iyM(%esp),%xmm4
        movss nb303nf_izM(%esp),%xmm5

        movlps nb303nf_ixH1(%esp),%xmm6
        movlps nb303nf_ixH2(%esp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm3
        movlps nb303nf_iyH1(%esp),%xmm6
        movlps nb303nf_iyH2(%esp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm4
        movlps nb303nf_izH1(%esp),%xmm6
        movlps nb303nf_izH2(%esp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm5

        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5

        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        ## rsq in xmm4 

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb303nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb303nf_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        ## a little trick to avoid NaNs: 
        ## positions 0,2,and 3 are valid, but not 1. 
        ## If it contains NaN it doesnt help to mult by 0, 
        ## So we shuffle it and copy pos 0 to pos1! 
        shufps $224,%xmm0,%xmm0 ## 11100000      

        mulps %xmm0,%xmm4       ## xmm4=r 
        movaps %xmm0,nb303nf_rinvM(%esp)

        mulps nb303nf_tsc(%esp),%xmm4
        movhlps %xmm4,%xmm7
        cvttps2pi %xmm4,%mm6
        cvttps2pi %xmm7,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm7
        movlhps %xmm7,%xmm3

        subps   %xmm3,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $2,%mm6
        pslld $2,%mm7

        movd %eax,%mm0
        movd %ecx,%mm1
        movd %edx,%mm2

        movl nb303nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm7,%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 ## 11011101
        ## coulomb table ready, in xmm4-xmm7      
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp        
        movaps nb303nf_qqM(%esp),%xmm0
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        ## increment vcoul 
        addps  nb303nf_vctot(%esp),%xmm5
        movaps %xmm5,nb303nf_vctot(%esp)

        decl nb303nf_innerk(%esp)
        jz    _nb_kernel303nf_ia32_sse.nb303nf_updateouterdata
        jmp   _nb_kernel303nf_ia32_sse.nb303nf_odd_loop
_nb_kernel303nf_ia32_sse.nb303nf_updateouterdata: 
        ## get n from stack
        movl nb303nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb303nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb303nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb303nf_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb303nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel303nf_ia32_sse.nb303nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb303nf_n(%esp)
        jmp _nb_kernel303nf_ia32_sse.nb303nf_outer
_nb_kernel303nf_ia32_sse.nb303nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb303nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel303nf_ia32_sse.nb303nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel303nf_ia32_sse.nb303nf_threadloop
_nb_kernel303nf_ia32_sse.nb303nf_end: 
        emms

        movl nb303nf_nouter(%esp),%eax
        movl nb303nf_ninner(%esp),%ebx
        movl nb303nf_outeriter(%ebp),%ecx
        movl nb303nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb303nf_salign(%esp),%eax
        addl %eax,%esp
        addl $408,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret


