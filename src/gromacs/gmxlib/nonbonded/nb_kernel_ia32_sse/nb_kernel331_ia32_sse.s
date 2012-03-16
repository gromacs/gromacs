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



.globl nb_kernel331_ia32_sse
.globl _nb_kernel331_ia32_sse
nb_kernel331_ia32_sse:  
_nb_kernel331_ia32_sse: 
.set nb331_p_nri, 8
.set nb331_iinr, 12
.set nb331_jindex, 16
.set nb331_jjnr, 20
.set nb331_shift, 24
.set nb331_shiftvec, 28
.set nb331_fshift, 32
.set nb331_gid, 36
.set nb331_pos, 40
.set nb331_faction, 44
.set nb331_charge, 48
.set nb331_p_facel, 52
.set nb331_argkrf, 56
.set nb331_argcrf, 60
.set nb331_Vc, 64
.set nb331_type, 68
.set nb331_p_ntype, 72
.set nb331_vdwparam, 76
.set nb331_Vvdw, 80
.set nb331_p_tabscale, 84
.set nb331_VFtab, 88
.set nb331_invsqrta, 92
.set nb331_dvda, 96
.set nb331_p_gbtabscale, 100
.set nb331_GBtab, 104
.set nb331_p_nthreads, 108
.set nb331_count, 112
.set nb331_mtx, 116
.set nb331_outeriter, 120
.set nb331_inneriter, 124
.set nb331_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb331_ixO, 0
.set nb331_iyO, 16
.set nb331_izO, 32
.set nb331_ixH1, 48
.set nb331_iyH1, 64
.set nb331_izH1, 80
.set nb331_ixH2, 96
.set nb331_iyH2, 112
.set nb331_izH2, 128
.set nb331_iqO, 144
.set nb331_iqH, 160
.set nb331_dxO, 176
.set nb331_dyO, 192
.set nb331_dzO, 208
.set nb331_dxH1, 224
.set nb331_dyH1, 240
.set nb331_dzH1, 256
.set nb331_dxH2, 272
.set nb331_dyH2, 288
.set nb331_dzH2, 304
.set nb331_qqO, 320
.set nb331_qqH, 336
.set nb331_rinvO, 352
.set nb331_rinvH1, 368
.set nb331_rinvH2, 384
.set nb331_rO, 400
.set nb331_rH1, 416
.set nb331_rH2, 432
.set nb331_tsc, 448
.set nb331_two, 464
.set nb331_c6, 480
.set nb331_c12, 496
.set nb331_vctot, 512
.set nb331_Vvdwtot, 528
.set nb331_fixO, 544
.set nb331_fiyO, 560
.set nb331_fizO, 576
.set nb331_fixH1, 592
.set nb331_fiyH1, 608
.set nb331_fizH1, 624
.set nb331_fixH2, 640
.set nb331_fiyH2, 656
.set nb331_fizH2, 672
.set nb331_fjx, 688
.set nb331_fjy, 704
.set nb331_fjz, 720
.set nb331_half, 736
.set nb331_three, 752
.set nb331_is3, 768
.set nb331_ii3, 772
.set nb331_ntia, 776
.set nb331_innerjjnr, 780
.set nb331_innerk, 784
.set nb331_n, 788
.set nb331_nn1, 792
.set nb331_nri, 796
.set nb331_nouter, 800
.set nb331_ninner, 804
.set nb331_salign, 808
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $812,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb331_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb331_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb331_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb331_nouter(%esp)
        movl %eax,nb331_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb331_half(%esp)
        movss nb331_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb331_half(%esp)
        movaps %xmm2,nb331_two(%esp)
        movaps %xmm3,nb331_three(%esp)
        movl nb331_p_tabscale(%ebp),%eax
        movss (%eax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb331_tsc(%esp)

        movl nb331_p_tabscale(%ebp),%eax
        movss (%eax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb331_tsc(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb331_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb331_charge(%ebp),%edx
        movss (%edx,%ebx,4),%xmm3
        movss 4(%edx,%ebx,4),%xmm4
        movl nb331_p_facel(%ebp),%esi
        movss (%esi),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb331_iqO(%esp)
        movaps %xmm4,nb331_iqH(%esp)

        movl  nb331_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb331_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb331_ntia(%esp)

_nb_kernel331_ia32_sse.nb331_threadloop: 
        movl  nb331_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel331_ia32_sse.nb331_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel331_ia32_sse.nb331_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb331_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb331_n(%esp)
        movl %ebx,nb331_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel331_ia32_sse.nb331_outerstart
        jmp _nb_kernel331_ia32_sse.nb331_end

_nb_kernel331_ia32_sse.nb331_outerstart: 
        ## ebx contains number of outer iterations
        addl nb331_nouter(%esp),%ebx
        movl %ebx,nb331_nouter(%esp)

_nb_kernel331_ia32_sse.nb331_outer: 
        movl  nb331_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb331_is3(%esp)      ## store is3 

        movl  nb331_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb331_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb331_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb331_ii3(%esp)

        addss (%eax,%ebx,4),%xmm3
        addss 4(%eax,%ebx,4),%xmm4
        addss 8(%eax,%ebx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb331_ixO(%esp)
        movaps %xmm4,nb331_iyO(%esp)
        movaps %xmm5,nb331_izO(%esp)

        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 12(%eax,%ebx,4),%xmm0
        addss 16(%eax,%ebx,4),%xmm1
        addss 20(%eax,%ebx,4),%xmm2
        addss 24(%eax,%ebx,4),%xmm3
        addss 28(%eax,%ebx,4),%xmm4
        addss 32(%eax,%ebx,4),%xmm5

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm0,nb331_ixH1(%esp)
        movaps %xmm1,nb331_iyH1(%esp)
        movaps %xmm2,nb331_izH1(%esp)
        movaps %xmm3,nb331_ixH2(%esp)
        movaps %xmm4,nb331_iyH2(%esp)
        movaps %xmm5,nb331_izH2(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb331_vctot(%esp)
        movaps %xmm4,nb331_Vvdwtot(%esp)
        movaps %xmm4,nb331_fixO(%esp)
        movaps %xmm4,nb331_fiyO(%esp)
        movaps %xmm4,nb331_fizO(%esp)
        movaps %xmm4,nb331_fixH1(%esp)
        movaps %xmm4,nb331_fiyH1(%esp)
        movaps %xmm4,nb331_fizH1(%esp)
        movaps %xmm4,nb331_fixH2(%esp)
        movaps %xmm4,nb331_fiyH2(%esp)
        movaps %xmm4,nb331_fizH2(%esp)

        movl  nb331_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb331_pos(%ebp),%esi
        movl  nb331_faction(%ebp),%edi
        movl  nb331_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb331_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb331_ninner(%esp),%ecx
        movl  %ecx,nb331_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb331_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel331_ia32_sse.nb331_unroll_loop
        jmp   _nb_kernel331_ia32_sse.nb331_odd_inner
_nb_kernel331_ia32_sse.nb331_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb331_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 

        addl $16,nb331_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb331_charge(%ebp),%esi     ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4           ## and in xmm4 
        mulps  nb331_iqO(%esp),%xmm3
        mulps  nb331_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps  %xmm3,nb331_qqO(%esp)
        movaps  %xmm4,nb331_qqH(%esp)

        movl nb331_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb331_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb331_ntia(%esp),%edi
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

        movaps %xmm4,nb331_c6(%esp)
        movaps %xmm6,nb331_c12(%esp)

        movl nb331_pos(%ebp),%esi        ## base of pos[] 

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

        ## move ixO-izO to xmm4-xmm6 
        movaps nb331_ixO(%esp),%xmm4
        movaps nb331_iyO(%esp),%xmm5
        movaps nb331_izO(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb331_dxO(%esp)
        movaps %xmm5,nb331_dyO(%esp)
        movaps %xmm6,nb331_dzO(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        movaps %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movaps nb331_ixH1(%esp),%xmm4
        movaps nb331_iyH1(%esp),%xmm5
        movaps nb331_izH1(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb331_dxH1(%esp)
        movaps %xmm5,nb331_dyH1(%esp)
        movaps %xmm6,nb331_dzH1(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm6
        addps %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movaps nb331_ixH2(%esp),%xmm3
        movaps nb331_iyH2(%esp),%xmm4
        movaps nb331_izH2(%esp),%xmm5

        ## calc dr 
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        ## store dr 
        movaps %xmm3,nb331_dxH2(%esp)
        movaps %xmm4,nb331_dyH2(%esp)
        movaps %xmm5,nb331_dzH2(%esp)
        ## square it 
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm4,%xmm5
        addps %xmm3,%xmm5
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## start with rsqO - seed to xmm2       
        rsqrtps %xmm7,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb331_three(%esp),%xmm4
        mulps   %xmm7,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb331_half(%esp),%xmm4
        movaps  %xmm4,nb331_rinvO(%esp)         ## rinvO in xmm4 
        mulps   %xmm4,%xmm7
        movaps  %xmm7,nb331_rO(%esp)

        ## rsqH1 - seed in xmm2 
        rsqrtps %xmm6,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb331_three(%esp),%xmm4
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb331_half(%esp),%xmm4
        movaps  %xmm4,nb331_rinvH1(%esp)        ## rinvH1 in xmm4 
        mulps   %xmm4,%xmm6
        movaps  %xmm6,nb331_rH1(%esp)

        ## rsqH2 - seed to xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb331_three(%esp),%xmm4
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb331_half(%esp),%xmm4
        movaps  %xmm4,nb331_rinvH2(%esp)        ## rinvH2 in xmm4 
        mulps   %xmm4,%xmm5
        movaps  %xmm5,nb331_rH2(%esp)

        ## do O interactions 
        ## rO is still in xmm7 
        mulps   nb331_tsc(%esp),%xmm7
        movhlps %xmm7,%xmm4
        cvttps2pi %xmm7,%mm6
        cvttps2pi %xmm4,%mm7   ## mm6/mm7 contain lu indices 

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

    movl nb331_VFtab(%ebp),%esi
    movd %mm6,%eax
    psrlq $32,%mm6
    movd %mm7,%ecx
    psrlq $32,%mm7
    movd %mm6,%ebx
    movd %mm7,%edx

    leal  (%eax,%eax,2),%eax
    leal  (%ebx,%ebx,2),%ebx
    leal  (%ecx,%ecx,2),%ecx
    leal  (%edx,%edx,2),%edx

    movlps (%esi,%eax,4),%xmm5
    movlps (%esi,%ecx,4),%xmm7
    movhps (%esi,%ebx,4),%xmm5
    movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

    movaps %xmm5,%xmm4
    shufps $136,%xmm7,%xmm4 ## constant 10001000
    shufps $221,%xmm7,%xmm5 ## constant 11011101

    movlps 8(%esi,%eax,4),%xmm7
    movlps 8(%esi,%ecx,4),%xmm3
    movhps 8(%esi,%ebx,4),%xmm7
    movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
    movaps %xmm7,%xmm6
    shufps $136,%xmm3,%xmm6 ## constant 10001000
    shufps $221,%xmm3,%xmm7 ## constant 11011101
    ## coulomb table ready, in xmm4-xmm7      

    mulps  %xmm1,%xmm6      ## xmm6=Geps 
    mulps  %xmm2,%xmm7      ## xmm7=Heps2 
    addps  %xmm6,%xmm5
    addps  %xmm7,%xmm5      ## xmm5=Fp        
    mulps  nb331_two(%esp),%xmm7         ## two*Heps2 
    movaps nb331_qqO(%esp),%xmm0
    addps  %xmm6,%xmm7
    addps  %xmm5,%xmm7 ## xmm7=FF 
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
    mulps  %xmm7,%xmm0 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm0 fijC 
    ## increment vcoul - then we can get rid of mm5 
    addps  nb331_vctot(%esp),%xmm5
    movaps %xmm5,nb331_vctot(%esp)

    ## dispersion 
    movlps 16(%esi,%eax,4),%xmm5
    movlps 16(%esi,%ecx,4),%xmm7
    movhps 16(%esi,%ebx,4),%xmm5
    movhps 16(%esi,%edx,4),%xmm7    ## got half dispersion table 
    movaps %xmm5,%xmm4
    shufps $136,%xmm7,%xmm4 ## constant 10001000
    shufps $221,%xmm7,%xmm5 ## constant 11011101

    movlps 24(%esi,%eax,4),%xmm7
    movlps 24(%esi,%ecx,4),%xmm3
    movhps 24(%esi,%ebx,4),%xmm7
    movhps 24(%esi,%edx,4),%xmm3    ## other half of dispersion table 
    movaps %xmm7,%xmm6
    shufps $136,%xmm3,%xmm6 ## constant 10001000
    shufps $221,%xmm3,%xmm7 ## constant 11011101
    ## dispersion table ready, in xmm4-xmm7  
    mulps  %xmm1,%xmm6      ## xmm6=Geps 
    mulps  %xmm2,%xmm7      ## xmm7=Heps2 
    addps  %xmm6,%xmm5
    addps  %xmm7,%xmm5      ## xmm5=Fp        
    mulps  nb331_two(%esp),%xmm7         ## two*Heps2 
    addps  %xmm6,%xmm7
    addps  %xmm5,%xmm7 ## xmm7=FF 
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 

    movaps nb331_c6(%esp),%xmm4
    mulps  %xmm4,%xmm7   ## fijD 
    mulps  %xmm4,%xmm5   ## Vvdw6 
    addps  %xmm7,%xmm0 ## add to fscal 

    ## Update Vvdwtot directly 
    addps  nb331_Vvdwtot(%esp),%xmm5
    movaps %xmm5,nb331_Vvdwtot(%esp)

    ## repulsion 
    movlps 32(%esi,%eax,4),%xmm5
    movlps 32(%esi,%ecx,4),%xmm7
    movhps 32(%esi,%ebx,4),%xmm5
    movhps 32(%esi,%edx,4),%xmm7    ## got half repulsion table 
    movaps %xmm5,%xmm4
    shufps $136,%xmm7,%xmm4 ## constant 10001000
    shufps $221,%xmm7,%xmm5 ## constant 11011101

    movlps 40(%esi,%eax,4),%xmm7
    movlps 40(%esi,%ecx,4),%xmm3
    movhps 40(%esi,%ebx,4),%xmm7
    movhps 40(%esi,%edx,4),%xmm3    ## other half of repulsion table 
    movaps %xmm7,%xmm6
    shufps $136,%xmm3,%xmm6 ## constant 10001000
    shufps $221,%xmm3,%xmm7 ## constant 11011101
    ## repulsion table ready, in xmm4-xmm7      
    mulps  %xmm1,%xmm6      ## xmm6=Geps 
    mulps  %xmm2,%xmm7      ## xmm7=Heps2 
    addps  %xmm6,%xmm5
    addps  %xmm7,%xmm5      ## xmm5=Fp        
    mulps  nb331_two(%esp),%xmm7         ## two*Heps2 
    addps  %xmm6,%xmm7
    addps  %xmm5,%xmm7 ## xmm7=FF 
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 

    movaps nb331_c12(%esp),%xmm4
    mulps  %xmm4,%xmm7   ## fijD 
    mulps  %xmm4,%xmm5   ## Vvdw12 
    addps  %xmm0,%xmm7 ## add to fscal 
    addps  nb331_Vvdwtot(%esp),%xmm5   ## total nonbonded potential in xmm5 
        xorps %xmm4,%xmm4

        mulps  nb331_rinvO(%esp),%xmm7   ## total fscal now in xmm7 

        mulps  nb331_tsc(%esp),%xmm7
    movaps %xmm5,nb331_Vvdwtot(%esp)
        subps  %xmm7,%xmm4

        movaps nb331_dxO(%esp),%xmm0
        movaps nb331_dyO(%esp),%xmm1
        movaps nb331_dzO(%esp),%xmm2
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2      ## tx in xmm0-xmm2 

        ## update O forces 
        movaps nb331_fixO(%esp),%xmm3
        movaps nb331_fiyO(%esp),%xmm4
        movaps nb331_fizO(%esp),%xmm7
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm7
        movaps %xmm3,nb331_fixO(%esp)
        movaps %xmm4,nb331_fiyO(%esp)
        movaps %xmm7,nb331_fizO(%esp)
        ## update j forces with water O 
        movaps %xmm0,nb331_fjx(%esp)
        movaps %xmm1,nb331_fjy(%esp)
        movaps %xmm2,nb331_fjz(%esp)

        ## Done with O interactions - now H1! 
        movaps nb331_rH1(%esp),%xmm7
        mulps   nb331_tsc(%esp),%xmm7
        movhlps %xmm7,%xmm4
        cvttps2pi %xmm7,%mm6
        cvttps2pi %xmm4,%mm7   ## mm6/mm7 contain lu indices 

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

    leal  (%eax,%eax,2),%eax
    leal  (%ebx,%ebx,2),%ebx
    leal  (%ecx,%ecx,2),%ecx
    leal  (%edx,%edx,2),%edx

    movlps (%esi,%eax,4),%xmm5
    movlps (%esi,%ecx,4),%xmm7
    movhps (%esi,%ebx,4),%xmm5
    movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

    movaps %xmm5,%xmm4
    shufps $136,%xmm7,%xmm4 ## constant 10001000
    shufps $221,%xmm7,%xmm5 ## constant 11011101

    movlps 8(%esi,%eax,4),%xmm7
    movlps 8(%esi,%ecx,4),%xmm3
    movhps 8(%esi,%ebx,4),%xmm7
    movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
    movaps %xmm7,%xmm6
    shufps $136,%xmm3,%xmm6 ## constant 10001000
    shufps $221,%xmm3,%xmm7 ## constant 11011101
    ## coulomb table ready, in xmm4-xmm7      

    mulps  %xmm1,%xmm6      ## xmm6=Geps 
    mulps  %xmm2,%xmm7      ## xmm7=Heps2 
    addps  %xmm6,%xmm5
    addps  %xmm7,%xmm5      ## xmm5=Fp        
    mulps  nb331_two(%esp),%xmm7         ## two*Heps2 
    movaps nb331_qqH(%esp),%xmm0
    addps  %xmm6,%xmm7
    addps  %xmm5,%xmm7 ## xmm7=FF 
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
    mulps  %xmm0,%xmm7 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm7 fijC 
    ## increment vcoul 
        xorps  %xmm4,%xmm4
    addps  nb331_vctot(%esp),%xmm5
        mulps  nb331_rinvH1(%esp),%xmm7
    movaps %xmm5,nb331_vctot(%esp)
        mulps  nb331_tsc(%esp),%xmm7
        subps %xmm7,%xmm4

        movaps nb331_dxH1(%esp),%xmm0
        movaps nb331_dyH1(%esp),%xmm1
        movaps nb331_dzH1(%esp),%xmm2
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2

        ## update H1 forces 
        movaps nb331_fixH1(%esp),%xmm3
        movaps nb331_fiyH1(%esp),%xmm4
        movaps nb331_fizH1(%esp),%xmm7
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm7
        movaps %xmm3,nb331_fixH1(%esp)
        movaps %xmm4,nb331_fiyH1(%esp)
        movaps %xmm7,nb331_fizH1(%esp)
        ## update j forces with water H1 
        addps  nb331_fjx(%esp),%xmm0
        addps  nb331_fjy(%esp),%xmm1
        addps  nb331_fjz(%esp),%xmm2
        movaps %xmm0,nb331_fjx(%esp)
        movaps %xmm1,nb331_fjy(%esp)
        movaps %xmm2,nb331_fjz(%esp)

        ## Done with H1, finally we do H2 interactions 
        movaps nb331_rH2(%esp),%xmm7
        mulps   nb331_tsc(%esp),%xmm7
        movhlps %xmm7,%xmm4
        cvttps2pi %xmm7,%mm6
        cvttps2pi %xmm4,%mm7   ## mm6/mm7 contain lu indices 

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

    leal  (%eax,%eax,2),%eax
    leal  (%ebx,%ebx,2),%ebx
    leal  (%ecx,%ecx,2),%ecx
    leal  (%edx,%edx,2),%edx

    movlps (%esi,%eax,4),%xmm5
    movlps (%esi,%ecx,4),%xmm7
    movhps (%esi,%ebx,4),%xmm5
    movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

    movaps %xmm5,%xmm4
    shufps $136,%xmm7,%xmm4 ## constant 10001000
    shufps $221,%xmm7,%xmm5 ## constant 11011101

    movlps 8(%esi,%eax,4),%xmm7
    movlps 8(%esi,%ecx,4),%xmm3
    movhps 8(%esi,%ebx,4),%xmm7
    movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
    movaps %xmm7,%xmm6
    shufps $136,%xmm3,%xmm6 ## constant 10001000
    shufps $221,%xmm3,%xmm7 ## constant 11011101
    ## coulomb table ready, in xmm4-xmm7      

    mulps  %xmm1,%xmm6      ## xmm6=Geps 
    mulps  %xmm2,%xmm7      ## xmm7=Heps2 
    addps  %xmm6,%xmm5
    addps  %xmm7,%xmm5      ## xmm5=Fp        
    mulps  nb331_two(%esp),%xmm7         ## two*Heps2 
    movaps nb331_qqH(%esp),%xmm0
    addps  %xmm6,%xmm7
    addps  %xmm5,%xmm7 ## xmm7=FF 
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
    mulps  %xmm0,%xmm7 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm0 fijC 
    ## increment vcoul 
        xorps  %xmm4,%xmm4
    addps  nb331_vctot(%esp),%xmm5
        mulps  nb331_rinvH2(%esp),%xmm7
    movaps %xmm5,nb331_vctot(%esp)
        mulps  nb331_tsc(%esp),%xmm7
        subps  %xmm7,%xmm4

        movaps nb331_dxH2(%esp),%xmm0
        movaps nb331_dyH2(%esp),%xmm1
        movaps nb331_dzH2(%esp),%xmm2
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2

    movd %mm0,%eax
    movd %mm1,%ebx
    movd %mm2,%ecx
    movd %mm3,%edx

        ## update H2 forces 
        movaps nb331_fixH2(%esp),%xmm3
        movaps nb331_fiyH2(%esp),%xmm4
        movaps nb331_fizH2(%esp),%xmm7
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm7
        movaps %xmm3,nb331_fixH2(%esp)
        movaps %xmm4,nb331_fiyH2(%esp)
        movaps %xmm7,nb331_fizH2(%esp)

        movl nb331_faction(%ebp),%edi
        ## update j forces 
        addps nb331_fjx(%esp),%xmm0
        addps nb331_fjy(%esp),%xmm1
        addps nb331_fjz(%esp),%xmm2

        movlps (%edi,%eax,4),%xmm4
        movlps (%edi,%ecx,4),%xmm7
        movhps (%edi,%ebx,4),%xmm4
        movhps (%edi,%edx,4),%xmm7

        movaps %xmm4,%xmm3
        shufps $136,%xmm7,%xmm3 ## constant 10001000
        shufps $221,%xmm7,%xmm4 ## constant 11011101                          
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
        shufps $229,%xmm2,%xmm2 ## constant 11100101
        subss  %xmm2,%xmm1
        shufps $234,%xmm2,%xmm2 ## constant 11101010
        subss  %xmm2,%xmm3
        shufps $255,%xmm2,%xmm2 ## constant 11111111
        subss  %xmm2,%xmm4
        movss  %xmm0,8(%edi,%eax,4)
        movss  %xmm1,8(%edi,%ebx,4)
        movss  %xmm3,8(%edi,%ecx,4)
        movss  %xmm4,8(%edi,%edx,4)

        ## should we do one more iteration? 
        subl $4,nb331_innerk(%esp)
        jl    _nb_kernel331_ia32_sse.nb331_odd_inner
        jmp   _nb_kernel331_ia32_sse.nb331_unroll_loop
_nb_kernel331_ia32_sse.nb331_odd_inner: 
        addl $4,nb331_innerk(%esp)
        jnz   _nb_kernel331_ia32_sse.nb331_odd_loop
        jmp   _nb_kernel331_ia32_sse.nb331_updateouterdata
_nb_kernel331_ia32_sse.nb331_odd_loop: 
        movl  nb331_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb331_innerjjnr(%esp)

        xorps %xmm4,%xmm4
        movss nb331_iqO(%esp),%xmm4
        movl nb331_charge(%ebp),%esi
        movhps nb331_iqH(%esp),%xmm4
        movss (%esi,%eax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb331_qqO(%esp)    ## use oxygen qq for storage 

        xorps %xmm6,%xmm6
        movl nb331_type(%ebp),%esi
        movl (%esi,%eax,4),%ebx
        movl nb331_vdwparam(%ebp),%esi
        shll %ebx
        addl nb331_ntia(%esp),%ebx
        movlps (%esi,%ebx,4),%xmm6
        movaps %xmm6,%xmm7
        shufps $252,%xmm6,%xmm6 ## constant 11111100
        shufps $253,%xmm7,%xmm7 ## constant 11111101
        movaps %xmm6,nb331_c6(%esp)
        movaps %xmm7,nb331_c12(%esp)

        movl nb331_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## move j coords to xmm0-xmm2 
        movss (%esi,%eax,4),%xmm0
        movss 4(%esi,%eax,4),%xmm1
        movss 8(%esi,%eax,4),%xmm2
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movss nb331_ixO(%esp),%xmm3
        movss nb331_iyO(%esp),%xmm4
        movss nb331_izO(%esp),%xmm5

        movlps nb331_ixH1(%esp),%xmm6
        movlps nb331_ixH2(%esp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm3
        movlps nb331_iyH1(%esp),%xmm6
        movlps nb331_iyH2(%esp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm4
        movlps nb331_izH1(%esp),%xmm6
        movlps nb331_izH2(%esp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm5

        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        movaps %xmm3,nb331_dxO(%esp)
        movaps %xmm4,nb331_dyO(%esp)
        movaps %xmm5,nb331_dzO(%esp)

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
        movaps nb331_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb331_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 

        ## a little trick to avoid NaNs: 
        ## positions 0,2,and 3 are valid, but not 1. 
        ## If it contains NaN it doesnt help to mult by 0, 
        ## So we shuffle it and copy pos 0 to pos1! 
        shufps $224,%xmm0,%xmm0 ## constant 11100000     

        mulps %xmm0,%xmm4       ## xmm4=r 
        movaps %xmm0,nb331_rinvO(%esp)

        mulps nb331_tsc(%esp),%xmm4
        movhlps %xmm4,%xmm7
        cvttps2pi %xmm4,%mm6
        cvttps2pi %xmm7,%mm7   ## mm6/mm7 contain lu indices 
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

    movl nb331_VFtab(%ebp),%esi
    movd %mm6,%eax
    movd %mm7,%ecx
    psrlq $32,%mm7
    movd %mm7,%edx

    leal  (%eax,%eax,2),%eax
    leal  (%ecx,%ecx,2),%ecx
    leal  (%edx,%edx,2),%edx

    movlps (%esi,%eax,4),%xmm5
    movlps (%esi,%ecx,4),%xmm7
    movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

    movaps %xmm5,%xmm4
    shufps $136,%xmm7,%xmm4 ## constant 10001000
    shufps $221,%xmm7,%xmm5 ## constant 11011101

    movlps 8(%esi,%eax,4),%xmm7
    movlps 8(%esi,%ecx,4),%xmm3
    movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
    movaps %xmm7,%xmm6
    shufps $136,%xmm3,%xmm6 ## constant 10001000
    shufps $221,%xmm3,%xmm7 ## constant 11011101
    ## coulomb table ready, in xmm4-xmm7      
    mulps  %xmm1,%xmm6      ## xmm6=Geps 
    mulps  %xmm2,%xmm7      ## xmm7=Heps2 
    addps  %xmm6,%xmm5
    addps  %xmm7,%xmm5      ## xmm5=Fp        
    mulps  nb331_two(%esp),%xmm7         ## two*Heps2 
    movaps nb331_qqO(%esp),%xmm0
    addps  %xmm6,%xmm7
    addps  %xmm5,%xmm7 ## xmm7=FF 
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
    mulps  %xmm7,%xmm0 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm0 fijC 
    ## increment vcoul - then we can get rid of mm5 
    addps  nb331_vctot(%esp),%xmm5
    movaps %xmm5,nb331_vctot(%esp)

    ## dispersion 
    movlps 16(%esi,%eax,4),%xmm5        ## half table 
    movaps %xmm5,%xmm4
    shufps $252,%xmm4,%xmm4 ## constant 11111100
    shufps $253,%xmm5,%xmm5 ## constant 11111101

    movlps 24(%esi,%eax,4),%xmm7    ## other half of dispersion table 
    movaps %xmm7,%xmm6
    shufps $252,%xmm6,%xmm6 ## constant 11111100
    shufps $253,%xmm7,%xmm7 ## constant 11111101
    ## dispersion table ready, in xmm4-xmm7  
    mulss  %xmm1,%xmm6      ## xmm6=Geps 
    mulss  %xmm2,%xmm7      ## xmm7=Heps2 
    addss  %xmm6,%xmm5  ## Update Vvdwtot directly 
    addss  %xmm7,%xmm5      ## xmm5=Fp        
    mulss  nb331_two(%esp),%xmm7         ## two*Heps2 
    addss  %xmm6,%xmm7
    addss  %xmm5,%xmm7 ## xmm7=FF 
    mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addss  %xmm4,%xmm5 ## xmm5=VV 

    movaps nb331_c6(%esp),%xmm4
    mulps  %xmm4,%xmm7   ## fijD 
    mulps  %xmm4,%xmm5   ## Vvdw6 
    addps  %xmm7,%xmm0 ## add to fscal 

    ## Update Vvdwtot directly 
    addps  nb331_Vvdwtot(%esp),%xmm5
    movaps %xmm5,nb331_Vvdwtot(%esp)

    ## repulsion 
    movlps 32(%esi,%eax,4),%xmm5    ## got half repulsion table 
    movaps %xmm5,%xmm4
    shufps $136,%xmm4,%xmm4 ## constant 10001000
    shufps $221,%xmm5,%xmm5 ## constant 11011101

    movlps 40(%esi,%eax,4),%xmm7    ## other half of repulsion table 
    movaps %xmm7,%xmm6
    shufps $136,%xmm6,%xmm6 ## constant 10001000
    shufps $221,%xmm7,%xmm7 ## constant 11011101
    ## repulsion table ready, in xmm4-xmm7      
    mulss  %xmm1,%xmm6      ## xmm6=Geps 
    mulss  %xmm2,%xmm7      ## xmm7=Heps2 
    addss  %xmm6,%xmm5
    addss  %xmm7,%xmm5      ## xmm5=Fp        
    mulss  nb331_two(%esp),%xmm7         ## two*Heps2 
    addss  %xmm6,%xmm7
    addss  %xmm5,%xmm7 ## xmm7=FF 
    mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addss  %xmm4,%xmm5 ## xmm5=VV 

    movaps nb331_c12(%esp),%xmm4
    mulps  %xmm4,%xmm7   ## fijD 
    mulps  %xmm4,%xmm5   ## Vvdw12 
    addps  %xmm0,%xmm7 ## add to fscal 
    addps  nb331_Vvdwtot(%esp),%xmm5   ## total nonbonded potential in xmm5 

        xorps  %xmm4,%xmm4
    movd %mm0,%eax
    movd %mm1,%ecx
    movd %mm2,%edx

        mulps  nb331_rinvO(%esp),%xmm7   ## total fscal now in xmm7 
    movaps %xmm5,nb331_Vvdwtot(%esp)
        mulps  nb331_tsc(%esp),%xmm7
        subps %xmm7,%xmm4

        movaps nb331_dxO(%esp),%xmm0
        movaps nb331_dyO(%esp),%xmm1
        movaps nb331_dzO(%esp),%xmm2

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2 ## xmm0-xmm2 now contains tx-tz (partial force) 
        movss  nb331_fixO(%esp),%xmm3
        movss  nb331_fiyO(%esp),%xmm4
        movss  nb331_fizO(%esp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,nb331_fixO(%esp)
        movss  %xmm4,nb331_fiyO(%esp)
        movss  %xmm5,nb331_fizO(%esp)   ## updated the O force now do the H's 
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        shufps $230,%xmm3,%xmm3 ## constant 11100110     ;# shift right 
        shufps $230,%xmm4,%xmm4 ## constant 11100110
        shufps $230,%xmm5,%xmm5 ## constant 11100110
        addss  nb331_fixH1(%esp),%xmm3
        addss  nb331_fiyH1(%esp),%xmm4
        addss  nb331_fizH1(%esp),%xmm5
        movss  %xmm3,nb331_fixH1(%esp)
        movss  %xmm4,nb331_fiyH1(%esp)
        movss  %xmm5,nb331_fizH1(%esp)          ## updated the H1 force 

        movl nb331_faction(%ebp),%edi
        shufps $231,%xmm3,%xmm3 ## constant 11100111     ;# shift right 
        shufps $231,%xmm4,%xmm4 ## constant 11100111
        shufps $231,%xmm5,%xmm5 ## constant 11100111
        addss  nb331_fixH2(%esp),%xmm3
        addss  nb331_fiyH2(%esp),%xmm4
        addss  nb331_fizH2(%esp),%xmm5
        movss  %xmm3,nb331_fixH2(%esp)
        movss  %xmm4,nb331_fiyH2(%esp)
        movss  %xmm5,nb331_fizH2(%esp)          ## updated the H2 force 

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
        addss    %xmm1,%xmm2   ## z sum in xmm2 
        subps    %xmm0,%xmm6
        subss    %xmm2,%xmm7

        movlps %xmm6,(%edi,%eax,4)
        movss  %xmm7,8(%edi,%eax,4)

        decl nb331_innerk(%esp)
        jz    _nb_kernel331_ia32_sse.nb331_updateouterdata
        jmp   _nb_kernel331_ia32_sse.nb331_odd_loop
_nb_kernel331_ia32_sse.nb331_updateouterdata: 
        movl  nb331_ii3(%esp),%ecx
        movl  nb331_faction(%ebp),%edi
        movl  nb331_fshift(%ebp),%esi
        movl  nb331_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movaps nb331_fixO(%esp),%xmm0
        movaps nb331_fiyO(%esp),%xmm1
        movaps nb331_fizO(%esp),%xmm2

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

        ## accumulate force in xmm6/xmm7 for fshift 
        movaps %xmm0,%xmm6
        movss %xmm2,%xmm7
        movlhps %xmm1,%xmm6
        shufps $8,%xmm6,%xmm6 ## constant 00001000      

        ## accumulate H1i forces in xmm0, xmm1, xmm2 
        movaps nb331_fixH1(%esp),%xmm0
        movaps nb331_fiyH1(%esp),%xmm1
        movaps nb331_fizH1(%esp),%xmm2

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
        addss %xmm2,%xmm7
        movlhps %xmm1,%xmm0
        shufps $8,%xmm0,%xmm0 ## constant 00001000      
        addps   %xmm0,%xmm6

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movaps nb331_fixH2(%esp),%xmm0
        movaps nb331_fiyH2(%esp),%xmm1
        movaps nb331_fizH2(%esp),%xmm2

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

        ## increment fshift force  
        movlps  (%esi,%edx,4),%xmm3
        movss  8(%esi,%edx,4),%xmm4
        addps  %xmm6,%xmm3
        addss  %xmm7,%xmm4
        movlps  %xmm3,(%esi,%edx,4)
        movss  %xmm4,8(%esi,%edx,4)

        ## get n from stack
        movl nb331_n(%esp),%esi
        ## get group index for i particle 
        movl  nb331_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb331_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb331_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb331_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb331_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb331_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel331_ia32_sse.nb331_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb331_n(%esp)
        jmp _nb_kernel331_ia32_sse.nb331_outer
_nb_kernel331_ia32_sse.nb331_outerend: 
        ## check if more outer neighborlists remain
        movl  nb331_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel331_ia32_sse.nb331_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel331_ia32_sse.nb331_threadloop
_nb_kernel331_ia32_sse.nb331_end: 
        emms

        movl nb331_nouter(%esp),%eax
        movl nb331_ninner(%esp),%ebx
        movl nb331_outeriter(%ebp),%ecx
        movl nb331_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb331_salign(%esp),%eax
        addl %eax,%esp
        addl $812,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret


.globl nb_kernel331nf_ia32_sse
.globl _nb_kernel331nf_ia32_sse
nb_kernel331nf_ia32_sse:        
_nb_kernel331nf_ia32_sse:       
.set nb331nf_p_nri, 8
.set nb331nf_iinr, 12
.set nb331nf_jindex, 16
.set nb331nf_jjnr, 20
.set nb331nf_shift, 24
.set nb331nf_shiftvec, 28
.set nb331nf_fshift, 32
.set nb331nf_gid, 36
.set nb331nf_pos, 40
.set nb331nf_faction, 44
.set nb331nf_charge, 48
.set nb331nf_p_facel, 52
.set nb331nf_argkrf, 56
.set nb331nf_argcrf, 60
.set nb331nf_Vc, 64
.set nb331nf_type, 68
.set nb331nf_p_ntype, 72
.set nb331nf_vdwparam, 76
.set nb331nf_Vvdw, 80
.set nb331nf_p_tabscale, 84
.set nb331nf_VFtab, 88
.set nb331nf_invsqrta, 92
.set nb331nf_dvda, 96
.set nb331nf_p_gbtabscale, 100
.set nb331nf_GBtab, 104
.set nb331nf_p_nthreads, 108
.set nb331nf_count, 112
.set nb331nf_mtx, 116
.set nb331nf_outeriter, 120
.set nb331nf_inneriter, 124
.set nb331nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb331nf_ixO, 0
.set nb331nf_iyO, 16
.set nb331nf_izO, 32
.set nb331nf_ixH1, 48
.set nb331nf_iyH1, 64
.set nb331nf_izH1, 80
.set nb331nf_ixH2, 96
.set nb331nf_iyH2, 112
.set nb331nf_izH2, 128
.set nb331nf_iqO, 144
.set nb331nf_iqH, 160
.set nb331nf_qqO, 176
.set nb331nf_qqH, 192
.set nb331nf_rinvO, 208
.set nb331nf_rinvH1, 224
.set nb331nf_rinvH2, 240
.set nb331nf_rO, 256
.set nb331nf_rH1, 272
.set nb331nf_rH2, 288
.set nb331nf_tsc, 304
.set nb331nf_c6, 320
.set nb331nf_c12, 336
.set nb331nf_vctot, 352
.set nb331nf_Vvdwtot, 368
.set nb331nf_half, 384
.set nb331nf_three, 400
.set nb331nf_is3, 416
.set nb331nf_ii3, 420
.set nb331nf_ntia, 424
.set nb331nf_innerjjnr, 428
.set nb331nf_innerk, 432
.set nb331nf_n, 436
.set nb331nf_nn1, 440
.set nb331nf_nri, 444
.set nb331nf_nouter, 448
.set nb331nf_ninner, 452
.set nb331nf_salign, 456
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $460,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb331nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb331nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb331nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb331nf_nouter(%esp)
        movl %eax,nb331nf_ninner(%esp)

        movl nb331nf_p_tabscale(%ebp),%eax
        movss (%eax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb331nf_tsc(%esp)
        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb331nf_half(%esp)
        movss nb331nf_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb331nf_half(%esp)
        movaps %xmm3,nb331nf_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb331nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb331nf_charge(%ebp),%edx
        movss (%edx,%ebx,4),%xmm3
        movss 4(%edx,%ebx,4),%xmm4
        movl nb331nf_p_facel(%ebp),%esi
        movss (%esi),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb331nf_iqO(%esp)
        movaps %xmm4,nb331nf_iqH(%esp)

        movl  nb331nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb331nf_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb331nf_ntia(%esp)

_nb_kernel331nf_ia32_sse.nb331nf_threadloop: 
        movl  nb331nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel331nf_ia32_sse.nb331nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel331nf_ia32_sse.nb331nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb331nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb331nf_n(%esp)
        movl %ebx,nb331nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel331nf_ia32_sse.nb331nf_outerstart
        jmp _nb_kernel331nf_ia32_sse.nb331nf_end

_nb_kernel331nf_ia32_sse.nb331nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb331nf_nouter(%esp),%ebx
        movl %ebx,nb331nf_nouter(%esp)

_nb_kernel331nf_ia32_sse.nb331nf_outer: 
        movl  nb331nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb331nf_is3(%esp)            ## store is3 

        movl  nb331nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb331nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb331nf_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb331nf_ii3(%esp)

        addss (%eax,%ebx,4),%xmm3
        addss 4(%eax,%ebx,4),%xmm4
        addss 8(%eax,%ebx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb331nf_ixO(%esp)
        movaps %xmm4,nb331nf_iyO(%esp)
        movaps %xmm5,nb331nf_izO(%esp)

        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 12(%eax,%ebx,4),%xmm0
        addss 16(%eax,%ebx,4),%xmm1
        addss 20(%eax,%ebx,4),%xmm2
        addss 24(%eax,%ebx,4),%xmm3
        addss 28(%eax,%ebx,4),%xmm4
        addss 32(%eax,%ebx,4),%xmm5

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm0,nb331nf_ixH1(%esp)
        movaps %xmm1,nb331nf_iyH1(%esp)
        movaps %xmm2,nb331nf_izH1(%esp)
        movaps %xmm3,nb331nf_ixH2(%esp)
        movaps %xmm4,nb331nf_iyH2(%esp)
        movaps %xmm5,nb331nf_izH2(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb331nf_vctot(%esp)
        movaps %xmm4,nb331nf_Vvdwtot(%esp)

        movl  nb331nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb331nf_pos(%ebp),%esi
        movl  nb331nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb331nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb331nf_ninner(%esp),%ecx
        movl  %ecx,nb331nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb331nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel331nf_ia32_sse.nb331nf_unroll_loop
        jmp   _nb_kernel331nf_ia32_sse.nb331nf_odd_inner
_nb_kernel331nf_ia32_sse.nb331nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb331nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 

        addl $16,nb331nf_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb331nf_charge(%ebp),%esi     ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4           ## and in xmm4 
        mulps  nb331nf_iqO(%esp),%xmm3
        mulps  nb331nf_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps  %xmm3,nb331nf_qqO(%esp)
        movaps  %xmm4,nb331nf_qqH(%esp)

        movl nb331nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb331nf_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb331nf_ntia(%esp),%edi
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

        movaps %xmm4,nb331nf_c6(%esp)
        movaps %xmm6,nb331nf_c12(%esp)

        movl nb331nf_pos(%ebp),%esi        ## base of pos[] 

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

        ## move ixO-izO to xmm4-xmm6 
        movaps nb331nf_ixO(%esp),%xmm4
        movaps nb331nf_iyO(%esp),%xmm5
        movaps nb331nf_izO(%esp),%xmm6

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
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movaps nb331nf_ixH1(%esp),%xmm4
        movaps nb331nf_iyH1(%esp),%xmm5
        movaps nb331nf_izH1(%esp),%xmm6

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
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movaps nb331nf_ixH2(%esp),%xmm3
        movaps nb331nf_iyH2(%esp),%xmm4
        movaps nb331nf_izH2(%esp),%xmm5

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
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## start with rsqO - seed to xmm2       
        rsqrtps %xmm7,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb331nf_three(%esp),%xmm4
        mulps   %xmm7,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb331nf_half(%esp),%xmm4
        movaps  %xmm4,nb331nf_rinvO(%esp)       ## rinvO in xmm4 
        mulps   %xmm4,%xmm7
        movaps  %xmm7,nb331nf_rO(%esp)

        ## rsqH1 - seed in xmm2 
        rsqrtps %xmm6,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb331nf_three(%esp),%xmm4
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb331nf_half(%esp),%xmm4
        movaps  %xmm4,nb331nf_rinvH1(%esp)      ## rinvH1 in xmm4 
        mulps   %xmm4,%xmm6
        movaps  %xmm6,nb331nf_rH1(%esp)

        ## rsqH2 - seed to xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb331nf_three(%esp),%xmm4
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb331nf_half(%esp),%xmm4
        movaps  %xmm4,nb331nf_rinvH2(%esp)      ## rinvH2 in xmm4 
        mulps   %xmm4,%xmm5
        movaps  %xmm5,nb331nf_rH2(%esp)

        ## do O interactions 
        ## rO is still in xmm7 
        mulps   nb331nf_tsc(%esp),%xmm7
        movhlps %xmm7,%xmm4
        cvttps2pi %xmm7,%mm6
        cvttps2pi %xmm4,%mm7   ## mm6/mm7 contain lu indices 

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

    movl nb331nf_VFtab(%ebp),%esi
    movd %mm6,%eax
    psrlq $32,%mm6
    movd %mm7,%ecx
    psrlq $32,%mm7
    movd %mm6,%ebx
    movd %mm7,%edx

    leal  (%eax,%eax,2),%eax
    leal  (%ebx,%ebx,2),%ebx
    leal  (%ecx,%ecx,2),%ecx
    leal  (%edx,%edx,2),%edx

    movlps (%esi,%eax,4),%xmm5
    movlps (%esi,%ecx,4),%xmm7
    movhps (%esi,%ebx,4),%xmm5
    movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

    movaps %xmm5,%xmm4
    shufps $136,%xmm7,%xmm4 ## constant 10001000
    shufps $221,%xmm7,%xmm5 ## constant 11011101

    movlps 8(%esi,%eax,4),%xmm7
    movlps 8(%esi,%ecx,4),%xmm3
    movhps 8(%esi,%ebx,4),%xmm7
    movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
    movaps %xmm7,%xmm6
    shufps $136,%xmm3,%xmm6 ## constant 10001000
    shufps $221,%xmm3,%xmm7 ## constant 11011101
    ## coulomb table ready, in xmm4-xmm7      

    mulps  %xmm1,%xmm6      ## xmm6=Geps 
    mulps  %xmm2,%xmm7      ## xmm7=Heps2 
    addps  %xmm6,%xmm5
    addps  %xmm7,%xmm5      ## xmm5=Fp        
    movaps nb331nf_qqO(%esp),%xmm0
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 
    ## increment vcoul - then we can get rid of mm5 
    addps  nb331nf_vctot(%esp),%xmm5
    movaps %xmm5,nb331nf_vctot(%esp)

    ## dispersion 
    movlps 16(%esi,%eax,4),%xmm5
    movlps 16(%esi,%ecx,4),%xmm7
    movhps 16(%esi,%ebx,4),%xmm5
    movhps 16(%esi,%edx,4),%xmm7    ## got half dispersion table 
    movaps %xmm5,%xmm4
    shufps $136,%xmm7,%xmm4 ## constant 10001000
    shufps $221,%xmm7,%xmm5 ## constant 11011101

    movlps 24(%esi,%eax,4),%xmm7
    movlps 24(%esi,%ecx,4),%xmm3
    movhps 24(%esi,%ebx,4),%xmm7
    movhps 24(%esi,%edx,4),%xmm3    ## other half of dispersion table 
    movaps %xmm7,%xmm6
    shufps $136,%xmm3,%xmm6 ## constant 10001000
    shufps $221,%xmm3,%xmm7 ## constant 11011101
    ## dispersion table ready, in xmm4-xmm7  
    mulps  %xmm1,%xmm6      ## xmm6=Geps 
    mulps  %xmm2,%xmm7      ## xmm7=Heps2 
    addps  %xmm6,%xmm5
    addps  %xmm7,%xmm5      ## xmm5=Fp        
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 

    movaps nb331nf_c6(%esp),%xmm4
    mulps  %xmm4,%xmm5   ## Vvdw6 
    ## Update Vvdwtot directly 
    addps  nb331nf_Vvdwtot(%esp),%xmm5
    movaps %xmm5,nb331nf_Vvdwtot(%esp)

    ## repulsion 
    movlps 32(%esi,%eax,4),%xmm5
    movlps 32(%esi,%ecx,4),%xmm7
    movhps 32(%esi,%ebx,4),%xmm5
    movhps 32(%esi,%edx,4),%xmm7    ## got half repulsion table 
    movaps %xmm5,%xmm4
    shufps $136,%xmm7,%xmm4 ## constant 10001000
    shufps $221,%xmm7,%xmm5 ## constant 11011101

    movlps 40(%esi,%eax,4),%xmm7
    movlps 40(%esi,%ecx,4),%xmm3
    movhps 40(%esi,%ebx,4),%xmm7
    movhps 40(%esi,%edx,4),%xmm3    ## other half of repulsion table 
    movaps %xmm7,%xmm6
    shufps $136,%xmm3,%xmm6 ## constant 10001000
    shufps $221,%xmm3,%xmm7 ## constant 11011101
    ## repulsion table ready, in xmm4-xmm7      
    mulps  %xmm1,%xmm6      ## xmm6=Geps 
    mulps  %xmm2,%xmm7      ## xmm7=Heps2 
    addps  %xmm6,%xmm5
    addps  %xmm7,%xmm5      ## xmm5=Fp        
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 

    movaps nb331nf_c12(%esp),%xmm4
    mulps  %xmm4,%xmm5   ## Vvdw12 
    addps  nb331nf_Vvdwtot(%esp),%xmm5   ## total nonbonded potential in xmm5 
    movaps %xmm5,nb331nf_Vvdwtot(%esp)

        ## Done with O interactions - now H1! 
        movaps nb331nf_rH1(%esp),%xmm7
        mulps   nb331nf_tsc(%esp),%xmm7
        movhlps %xmm7,%xmm4
        cvttps2pi %xmm7,%mm6
        cvttps2pi %xmm4,%mm7   ## mm6/mm7 contain lu indices 

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

    leal  (%eax,%eax,2),%eax
    leal  (%ebx,%ebx,2),%ebx
    leal  (%ecx,%ecx,2),%ecx
    leal  (%edx,%edx,2),%edx

    movlps (%esi,%eax,4),%xmm5
    movlps (%esi,%ecx,4),%xmm7
    movhps (%esi,%ebx,4),%xmm5
    movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

    movaps %xmm5,%xmm4
    shufps $136,%xmm7,%xmm4 ## constant 10001000
    shufps $221,%xmm7,%xmm5 ## constant 11011101

    movlps 8(%esi,%eax,4),%xmm7
    movlps 8(%esi,%ecx,4),%xmm3
    movhps 8(%esi,%ebx,4),%xmm7
    movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
    movaps %xmm7,%xmm6
    shufps $136,%xmm3,%xmm6 ## constant 10001000
    shufps $221,%xmm3,%xmm7 ## constant 11011101
    ## coulomb table ready, in xmm4-xmm7      

    mulps  %xmm1,%xmm6      ## xmm6=Geps 
    mulps  %xmm2,%xmm7      ## xmm7=Heps2 
    addps  %xmm6,%xmm5
    addps  %xmm7,%xmm5      ## xmm5=Fp        
    movaps nb331nf_qqH(%esp),%xmm0
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 
    ## increment vcoul 
    addps  nb331nf_vctot(%esp),%xmm5
    movaps %xmm5,nb331nf_vctot(%esp)

        ## Done with H1, finally we do H2 interactions 
        movaps nb331nf_rH2(%esp),%xmm7
        mulps   nb331nf_tsc(%esp),%xmm7
        movhlps %xmm7,%xmm4
        cvttps2pi %xmm7,%mm6
        cvttps2pi %xmm4,%mm7   ## mm6/mm7 contain lu indices 

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

    leal  (%eax,%eax,2),%eax
    leal  (%ebx,%ebx,2),%ebx
    leal  (%ecx,%ecx,2),%ecx
    leal  (%edx,%edx,2),%edx

    movlps (%esi,%eax,4),%xmm5
    movlps (%esi,%ecx,4),%xmm7
    movhps (%esi,%ebx,4),%xmm5
    movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

    movaps %xmm5,%xmm4
    shufps $136,%xmm7,%xmm4 ## constant 10001000
    shufps $221,%xmm7,%xmm5 ## constant 11011101

    movlps 8(%esi,%eax,4),%xmm7
    movlps 8(%esi,%ecx,4),%xmm3
    movhps 8(%esi,%ebx,4),%xmm7
    movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
    movaps %xmm7,%xmm6
    shufps $136,%xmm3,%xmm6 ## constant 10001000
    shufps $221,%xmm3,%xmm7 ## constant 11011101
    ## coulomb table ready, in xmm4-xmm7      

    mulps  %xmm1,%xmm6      ## xmm6=Geps 
    mulps  %xmm2,%xmm7      ## xmm7=Heps2 
    addps  %xmm6,%xmm5
    addps  %xmm7,%xmm5      ## xmm5=Fp        
    movaps nb331nf_qqH(%esp),%xmm0
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm0,%xmm5 ## vcoul=qq*VV 
    ## at this point mm5 contains vcoul 
    ## increment vcoul 
    addps  nb331nf_vctot(%esp),%xmm5
    movaps %xmm5,nb331nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $4,nb331nf_innerk(%esp)
        jl    _nb_kernel331nf_ia32_sse.nb331nf_odd_inner
        jmp   _nb_kernel331nf_ia32_sse.nb331nf_unroll_loop
_nb_kernel331nf_ia32_sse.nb331nf_odd_inner: 
        addl $4,nb331nf_innerk(%esp)
        jnz   _nb_kernel331nf_ia32_sse.nb331nf_odd_loop
        jmp   _nb_kernel331nf_ia32_sse.nb331nf_updateouterdata
_nb_kernel331nf_ia32_sse.nb331nf_odd_loop: 
        movl  nb331nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb331nf_innerjjnr(%esp)

        xorps %xmm4,%xmm4
        movss nb331nf_iqO(%esp),%xmm4
        movl nb331nf_charge(%ebp),%esi
        movhps nb331nf_iqH(%esp),%xmm4
        movss (%esi,%eax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb331nf_qqO(%esp)          ## use oxygen qq for storage 

        xorps %xmm6,%xmm6
        movl nb331nf_type(%ebp),%esi
        movl (%esi,%eax,4),%ebx
        movl nb331nf_vdwparam(%ebp),%esi
        shll %ebx
        addl nb331nf_ntia(%esp),%ebx
        movlps (%esi,%ebx,4),%xmm6
        movaps %xmm6,%xmm7
        shufps $252,%xmm6,%xmm6 ## constant 11111100
        shufps $253,%xmm7,%xmm7 ## constant 11111101
        movaps %xmm6,nb331nf_c6(%esp)
        movaps %xmm7,nb331nf_c12(%esp)

        movl nb331nf_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## move j coords to xmm0-xmm2 
        movss (%esi,%eax,4),%xmm0
        movss 4(%esi,%eax,4),%xmm1
        movss 8(%esi,%eax,4),%xmm2
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movss nb331nf_ixO(%esp),%xmm3
        movss nb331nf_iyO(%esp),%xmm4
        movss nb331nf_izO(%esp),%xmm5

        movlps nb331nf_ixH1(%esp),%xmm6
        movlps nb331nf_ixH2(%esp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm3
        movlps nb331nf_iyH1(%esp),%xmm6
        movlps nb331nf_iyH2(%esp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm4
        movlps nb331nf_izH1(%esp),%xmm6
        movlps nb331nf_izH2(%esp),%xmm7
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
        movaps nb331nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb331nf_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        ## a little trick to avoid NaNs: 
        ## positions 0,2,and 3 are valid, but not 1. 
        ## If it contains NaN it doesnt help to mult by 0, 
        ## So we shuffle it and copy pos 0 to pos1! 
        shufps $224,%xmm0,%xmm0 ## constant 11100000     

        mulps %xmm0,%xmm4       ## xmm4=r 
        movaps %xmm0,nb331nf_rinvO(%esp)

        mulps nb331nf_tsc(%esp),%xmm4
        movhlps %xmm4,%xmm7
        cvttps2pi %xmm4,%mm6
        cvttps2pi %xmm7,%mm7   ## mm6/mm7 contain lu indices 
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

    movl nb331nf_VFtab(%ebp),%esi
    movd %mm6,%eax
    movd %mm7,%ecx
    psrlq $32,%mm7
    movd %mm7,%edx

    leal  (%eax,%eax,2),%eax
    leal  (%ecx,%ecx,2),%ecx
    leal  (%edx,%edx,2),%edx

    movlps (%esi,%eax,4),%xmm5
    movlps (%esi,%ecx,4),%xmm7
    movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

    movaps %xmm5,%xmm4
    shufps $136,%xmm7,%xmm4 ## constant 10001000
    shufps $221,%xmm7,%xmm5 ## constant 11011101

    movlps 8(%esi,%eax,4),%xmm7
    movlps 8(%esi,%ecx,4),%xmm3
    movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
    movaps %xmm7,%xmm6
    shufps $136,%xmm3,%xmm6 ## constant 10001000
    shufps $221,%xmm3,%xmm7 ## constant 11011101
    ## coulomb table ready, in xmm4-xmm7      
    mulps  %xmm1,%xmm6      ## xmm6=Geps 
    mulps  %xmm2,%xmm7      ## xmm7=Heps2 
    addps  %xmm6,%xmm5
    addps  %xmm7,%xmm5      ## xmm5=Fp        
    movaps nb331nf_qqO(%esp),%xmm0
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm0,%xmm5 ## vcoul=qq*VV 
    ## at this point mm5 contains vcoul 
    ## increment vcoul - then we can get rid of mm5 
    addps  nb331nf_vctot(%esp),%xmm5
    movaps %xmm5,nb331nf_vctot(%esp)

    ## dispersion 
    movlps 16(%esi,%eax,4),%xmm5        ## half table 
    movaps %xmm5,%xmm4
    shufps $252,%xmm4,%xmm4 ## constant 11111100
    shufps $253,%xmm5,%xmm5 ## constant 11111101

    movlps 24(%esi,%eax,4),%xmm7    ## other half of dispersion table 
    movaps %xmm7,%xmm6
    shufps $252,%xmm6,%xmm6 ## constant 11111100
    shufps $253,%xmm7,%xmm7 ## constant 11111101
    ## dispersion table ready, in xmm4-xmm7  
    mulss  %xmm1,%xmm6      ## xmm6=Geps 
    mulss  %xmm2,%xmm7      ## xmm7=Heps2 
    addss  %xmm6,%xmm5  ## Update Vvdwtot directly 
    addss  %xmm7,%xmm5      ## xmm5=Fp        
    mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addss  %xmm4,%xmm5 ## xmm5=VV 

    movaps nb331nf_c6(%esp),%xmm4
    mulps  %xmm4,%xmm5   ## Vvdw6 
    ## Update Vvdwtot directly 
    addps  nb331nf_Vvdwtot(%esp),%xmm5
    movaps %xmm5,nb331nf_Vvdwtot(%esp)

    ## repulsion 
    movlps 32(%esi,%eax,4),%xmm5    ## got half repulsion table 
    movaps %xmm5,%xmm4
    shufps $136,%xmm4,%xmm4 ## constant 10001000
    shufps $221,%xmm5,%xmm5 ## constant 11011101

    movlps 40(%esi,%eax,4),%xmm7    ## other half of repulsion table 
    movaps %xmm7,%xmm6
    shufps $136,%xmm6,%xmm6 ## constant 10001000
    shufps $221,%xmm7,%xmm7 ## constant 11011101
    ## repulsion table ready, in xmm4-xmm7      
    mulss  %xmm1,%xmm6      ## xmm6=Geps 
    mulss  %xmm2,%xmm7      ## xmm7=Heps2 
    addss  %xmm6,%xmm5
    addss  %xmm7,%xmm5      ## xmm5=Fp        
    mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addss  %xmm4,%xmm5 ## xmm5=VV 

    movaps nb331nf_c12(%esp),%xmm4
    mulps  %xmm4,%xmm5   ## Vvdw12 
    addps  nb331nf_Vvdwtot(%esp),%xmm5   ## total nonbonded potential in xmm5 
    movaps %xmm5,nb331nf_Vvdwtot(%esp)

        decl nb331nf_innerk(%esp)
        jz    _nb_kernel331nf_ia32_sse.nb331nf_updateouterdata
        jmp   _nb_kernel331nf_ia32_sse.nb331nf_odd_loop
_nb_kernel331nf_ia32_sse.nb331nf_updateouterdata: 
        ## get n from stack
        movl nb331nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb331nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb331nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb331nf_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb331nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb331nf_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb331nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel331nf_ia32_sse.nb331nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb331nf_n(%esp)
        jmp _nb_kernel331nf_ia32_sse.nb331nf_outer
_nb_kernel331nf_ia32_sse.nb331nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb331nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel331nf_ia32_sse.nb331nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel331nf_ia32_sse.nb331nf_threadloop
_nb_kernel331nf_ia32_sse.nb331nf_end: 
        emms

        movl nb331nf_nouter(%esp),%eax
        movl nb331nf_ninner(%esp),%ebx
        movl nb331nf_outeriter(%ebp),%ecx
        movl nb331nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb331nf_salign(%esp),%eax
        addl %eax,%esp
        addl $460,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret


