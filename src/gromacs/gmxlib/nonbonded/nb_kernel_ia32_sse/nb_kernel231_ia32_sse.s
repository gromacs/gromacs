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



.globl nb_kernel231_ia32_sse
.globl _nb_kernel231_ia32_sse
nb_kernel231_ia32_sse:  
_nb_kernel231_ia32_sse: 
.set nb231_p_nri, 8
.set nb231_iinr, 12
.set nb231_jindex, 16
.set nb231_jjnr, 20
.set nb231_shift, 24
.set nb231_shiftvec, 28
.set nb231_fshift, 32
.set nb231_gid, 36
.set nb231_pos, 40
.set nb231_faction, 44
.set nb231_charge, 48
.set nb231_p_facel, 52
.set nb231_argkrf, 56
.set nb231_argcrf, 60
.set nb231_Vc, 64
.set nb231_type, 68
.set nb231_p_ntype, 72
.set nb231_vdwparam, 76
.set nb231_Vvdw, 80
.set nb231_p_tabscale, 84
.set nb231_VFtab, 88
.set nb231_invsqrta, 92
.set nb231_dvda, 96
.set nb231_p_gbtabscale, 100
.set nb231_GBtab, 104
.set nb231_p_nthreads, 108
.set nb231_count, 112
.set nb231_mtx, 116
.set nb231_outeriter, 120
.set nb231_inneriter, 124
.set nb231_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb231_ixO, 0
.set nb231_iyO, 16
.set nb231_izO, 32
.set nb231_ixH1, 48
.set nb231_iyH1, 64
.set nb231_izH1, 80
.set nb231_ixH2, 96
.set nb231_iyH2, 112
.set nb231_izH2, 128
.set nb231_iqO, 144
.set nb231_iqH, 160
.set nb231_dxO, 176
.set nb231_dyO, 192
.set nb231_dzO, 208
.set nb231_dxH1, 224
.set nb231_dyH1, 240
.set nb231_dzH1, 256
.set nb231_dxH2, 272
.set nb231_dyH2, 288
.set nb231_dzH2, 304
.set nb231_qqO, 320
.set nb231_qqH, 336
.set nb231_c6, 352
.set nb231_c12, 368
.set nb231_tsc, 384
.set nb231_fstmp, 400
.set nb231_vctot, 416
.set nb231_Vvdwtot, 432
.set nb231_fixO, 448
.set nb231_fiyO, 464
.set nb231_fizO, 480
.set nb231_fixH1, 496
.set nb231_fiyH1, 512
.set nb231_fizH1, 528
.set nb231_fixH2, 544
.set nb231_fiyH2, 560
.set nb231_fizH2, 576
.set nb231_fjx, 592
.set nb231_fjy, 608
.set nb231_fjz, 624
.set nb231_half, 640
.set nb231_three, 656
.set nb231_two, 672
.set nb231_krf, 688
.set nb231_crf, 704
.set nb231_krsqO, 720
.set nb231_krsqH1, 736
.set nb231_krsqH2, 752
.set nb231_rinvO, 768
.set nb231_rinvH1, 784
.set nb231_rinvH2, 800
.set nb231_is3, 816
.set nb231_ii3, 820
.set nb231_ntia, 824
.set nb231_innerjjnr, 828
.set nb231_innerk, 832
.set nb231_n, 836
.set nb231_nn1, 840
.set nb231_nri, 844
.set nb231_nouter, 848
.set nb231_ninner, 852
.set nb231_salign, 856
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $860,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb231_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb231_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb231_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb231_nouter(%esp)
        movl %eax,nb231_ninner(%esp)

        movl nb231_p_tabscale(%ebp),%eax
        movss (%eax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb231_tsc(%esp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb231_half(%esp)
        movss nb231_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb231_half(%esp)
        movaps %xmm2,nb231_two(%esp)
        movaps %xmm3,nb231_three(%esp)

        movl nb231_argkrf(%ebp),%esi
        movl nb231_argcrf(%ebp),%edi
        movss (%esi),%xmm5
        movss (%edi),%xmm6
        shufps $0,%xmm5,%xmm5
        shufps $0,%xmm6,%xmm6
        movaps %xmm5,nb231_krf(%esp)
        movaps %xmm6,nb231_crf(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb231_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb231_charge(%ebp),%edx
        movss (%edx,%ebx,4),%xmm3
        movss 4(%edx,%ebx,4),%xmm4
        movl nb231_p_facel(%ebp),%esi
        movss (%esi),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb231_iqO(%esp)
        movaps %xmm4,nb231_iqH(%esp)

        movl  nb231_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb231_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb231_ntia(%esp)

_nb_kernel231_ia32_sse.nb231_threadloop: 
        movl  nb231_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel231_ia32_sse.nb231_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel231_ia32_sse.nb231_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb231_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb231_n(%esp)
        movl %ebx,nb231_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel231_ia32_sse.nb231_outerstart
        jmp _nb_kernel231_ia32_sse.nb231_end

_nb_kernel231_ia32_sse.nb231_outerstart: 
        ## ebx contains number of outer iterations
        addl nb231_nouter(%esp),%ebx
        movl %ebx,nb231_nouter(%esp)

_nb_kernel231_ia32_sse.nb231_outer: 
        movl  nb231_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb231_is3(%esp)      ## store is3 

        movl  nb231_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb231_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb231_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb231_ii3(%esp)

        addss (%eax,%ebx,4),%xmm3
        addss 4(%eax,%ebx,4),%xmm4
        addss 8(%eax,%ebx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb231_ixO(%esp)
        movaps %xmm4,nb231_iyO(%esp)
        movaps %xmm5,nb231_izO(%esp)

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
        movaps %xmm0,nb231_ixH1(%esp)
        movaps %xmm1,nb231_iyH1(%esp)
        movaps %xmm2,nb231_izH1(%esp)
        movaps %xmm3,nb231_ixH2(%esp)
        movaps %xmm4,nb231_iyH2(%esp)
        movaps %xmm5,nb231_izH2(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb231_vctot(%esp)
        movaps %xmm4,nb231_Vvdwtot(%esp)
        movaps %xmm4,nb231_fixO(%esp)
        movaps %xmm4,nb231_fiyO(%esp)
        movaps %xmm4,nb231_fizO(%esp)
        movaps %xmm4,nb231_fixH1(%esp)
        movaps %xmm4,nb231_fiyH1(%esp)
        movaps %xmm4,nb231_fizH1(%esp)
        movaps %xmm4,nb231_fixH2(%esp)
        movaps %xmm4,nb231_fiyH2(%esp)
        movaps %xmm4,nb231_fizH2(%esp)

        movl  nb231_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb231_pos(%ebp),%esi
        movl  nb231_faction(%ebp),%edi
        movl  nb231_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb231_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb231_ninner(%esp),%ecx
        movl  %ecx,nb231_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb231_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel231_ia32_sse.nb231_unroll_loop
        jmp   _nb_kernel231_ia32_sse.nb231_odd_inner
_nb_kernel231_ia32_sse.nb231_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb231_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 

        addl $16,nb231_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb231_charge(%ebp),%esi     ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4           ## and in xmm4 
        mulps  nb231_iqO(%esp),%xmm3
        mulps  nb231_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps  %xmm3,nb231_qqO(%esp)
        movaps  %xmm4,nb231_qqH(%esp)

        movl nb231_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb231_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb231_ntia(%esp),%edi
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

        movaps %xmm4,nb231_c6(%esp)
        movaps %xmm6,nb231_c12(%esp)

        movl nb231_pos(%ebp),%esi        ## base of pos[] 

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
        movaps nb231_ixO(%esp),%xmm4
        movaps nb231_iyO(%esp),%xmm5
        movaps nb231_izO(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb231_dxO(%esp)
        movaps %xmm5,nb231_dyO(%esp)
        movaps %xmm6,nb231_dzO(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        movaps %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movaps nb231_ixH1(%esp),%xmm4
        movaps nb231_iyH1(%esp),%xmm5
        movaps nb231_izH1(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb231_dxH1(%esp)
        movaps %xmm5,nb231_dyH1(%esp)
        movaps %xmm6,nb231_dzH1(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm6
        addps %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movaps nb231_ixH2(%esp),%xmm3
        movaps nb231_iyH2(%esp),%xmm4
        movaps nb231_izH2(%esp),%xmm5

        ## calc dr 
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        ## store dr 
        movaps %xmm3,nb231_dxH2(%esp)
        movaps %xmm4,nb231_dyH2(%esp)
        movaps %xmm5,nb231_dzH2(%esp)
        ## square it 
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm4,%xmm5
        addps %xmm3,%xmm5
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7


        movaps %xmm5,%xmm0
        movaps %xmm6,%xmm1
        movaps %xmm7,%xmm2

        mulps  nb231_krf(%esp),%xmm0
        mulps  nb231_krf(%esp),%xmm1
        mulps  nb231_krf(%esp),%xmm2

        movaps %xmm0,nb231_krsqH2(%esp)
        movaps %xmm1,nb231_krsqH1(%esp)
        movaps %xmm2,nb231_krsqO(%esp)

        ## start with rsqO (still in xmm7) - seed to xmm2       
        rsqrtps %xmm7,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb231_three(%esp),%xmm4
        mulps   %xmm7,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb231_half(%esp),%xmm4
        movaps  %xmm4,nb231_rinvO(%esp)         ## rinvO in xmm0, rsqO in xmm7

        ## rsqH1 - seed in xmm2 
        rsqrtps %xmm6,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb231_three(%esp),%xmm4
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb231_half(%esp),%xmm4
        movaps  %xmm4,nb231_rinvH1(%esp)        ## rinvH1 in xmm6 

        ## rsqH2 - seed in xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb231_three(%esp),%xmm4
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb231_half(%esp),%xmm4
        movaps  %xmm4,nb231_rinvH2(%esp)        ## rinvH2 in xmm5 


        ## do O table interactions - rsqO in xmm7.
        mulps nb231_rinvO(%esp),%xmm7
        mulps nb231_tsc(%esp),%xmm7   ## rtab

        movhlps %xmm7,%xmm5
        cvttps2pi %xmm7,%mm6
        cvttps2pi %xmm5,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        cvtpi2ps %mm7,%xmm5
        movlhps %xmm5,%xmm6
        subps  %xmm6,%xmm7
        movaps %xmm7,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $3,%mm6
        pslld $3,%mm7

        movd %eax,%mm0
        movd %ebx,%mm1
        movd %ecx,%mm2
        movd %edx,%mm3

        movl nb231_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        ## dispersion 
        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half dispersion table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of dispersion table 
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## dispersion table ready, in xmm4-xmm7         

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  nb231_two(%esp),%xmm7    ## two*Heps2 
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb231_c6(%esp),%xmm4
        mulps  %xmm4,%xmm7       ## fijD 
        mulps  %xmm4,%xmm5       ## Vvdw6 
        mulps  nb231_tsc(%esp),%xmm7
        ## put scalar force on stack Update Vvdwtot directly 
        addps  nb231_Vvdwtot(%esp),%xmm5
        movaps %xmm7,nb231_fstmp(%esp)
        movaps %xmm5,nb231_Vvdwtot(%esp)

        ## repulsion 
        movlps 16(%esi,%eax,4),%xmm5
        movlps 16(%esi,%ecx,4),%xmm7
        movhps 16(%esi,%ebx,4),%xmm5
        movhps 16(%esi,%edx,4),%xmm7    ## got half repulsion table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 24(%esi,%eax,4),%xmm7
        movlps 24(%esi,%ecx,4),%xmm3
        movhps 24(%esi,%ebx,4),%xmm7
        movhps 24(%esi,%edx,4),%xmm3    ## other half of repulsion table 
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## table ready, in xmm4-xmm7    
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  nb231_two(%esp),%xmm7    ## two*Heps2 
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb231_c12(%esp),%xmm4
        mulps  %xmm4,%xmm7 ## fijR 
        mulps  %xmm4,%xmm5 ## Vvdw12 
        mulps  nb231_tsc(%esp),%xmm7
        addps  nb231_fstmp(%esp),%xmm7
        movaps %xmm7,nb231_fstmp(%esp)

        addps  nb231_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb231_Vvdwtot(%esp)

        movaps nb231_rinvO(%esp),%xmm2
        movaps nb231_krsqO(%esp),%xmm1
        movaps %xmm2,%xmm3
        addps  %xmm1,%xmm2 ## rinv+krsq
        subps  nb231_crf(%esp),%xmm2   ## rinv+krsq-crf
        movaps %xmm3,%xmm0
        subps  %xmm1,%xmm3
        subps  %xmm1,%xmm3 ## rinv-2*krsq
        mulps  nb231_qqO(%esp),%xmm2
        mulps  nb231_qqO(%esp),%xmm3

        mulps  %xmm0,%xmm3
        subps  nb231_fstmp(%esp),%xmm3
        mulps  %xmm0,%xmm3

        addps  nb231_vctot(%esp),%xmm2
        movaps %xmm2,nb231_vctot(%esp)

        movaps nb231_dxO(%esp),%xmm0
        movaps nb231_dyO(%esp),%xmm1
        movaps nb231_dzO(%esp),%xmm2
        mulps  %xmm3,%xmm0
        mulps  %xmm3,%xmm1
        mulps  %xmm3,%xmm2

        movd %mm0,%eax
        movd %mm1,%ebx
        movd %mm2,%ecx
        movd %mm3,%edx

        ## update O forces 
        movaps nb231_fixO(%esp),%xmm3
        movaps nb231_fiyO(%esp),%xmm4
        movaps nb231_fizO(%esp),%xmm7
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm7
        movaps %xmm3,nb231_fixO(%esp)
        movaps %xmm4,nb231_fiyO(%esp)
        movaps %xmm7,nb231_fizO(%esp)
        ## update j forces with water O 
        movaps %xmm0,nb231_fjx(%esp)
        movaps %xmm1,nb231_fjy(%esp)
        movaps %xmm2,nb231_fjz(%esp)

        ## H1 interactions 
        movaps nb231_rinvH1(%esp),%xmm2
        movaps nb231_krsqH1(%esp),%xmm1
        movaps %xmm2,%xmm3
        addps  %xmm1,%xmm2 ## rinv+krsq
        subps  nb231_crf(%esp),%xmm2   ## rinv+krsq-crf
        movaps %xmm3,%xmm0
        subps  %xmm1,%xmm3
        subps  %xmm1,%xmm3 ## rinv-2*krsq
        mulps  %xmm0,%xmm0
        mulps  nb231_qqH(%esp),%xmm2
        mulps  nb231_qqH(%esp),%xmm3

        mulps  %xmm0,%xmm3
        addps  nb231_vctot(%esp),%xmm2
        movaps %xmm2,nb231_vctot(%esp)

        movaps nb231_dxH1(%esp),%xmm0
        movaps nb231_dyH1(%esp),%xmm1
        movaps nb231_dzH1(%esp),%xmm2

        mulps  %xmm3,%xmm0
        mulps  %xmm3,%xmm1
        mulps  %xmm3,%xmm2

        ## update H1 forces 
        movaps nb231_fixH1(%esp),%xmm3
        movaps nb231_fiyH1(%esp),%xmm4
        movaps nb231_fizH1(%esp),%xmm7
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm7
        movaps %xmm3,nb231_fixH1(%esp)
        movaps %xmm4,nb231_fiyH1(%esp)
        movaps %xmm7,nb231_fizH1(%esp)
        ## update j forces with water H1 
        addps  nb231_fjx(%esp),%xmm0
        addps  nb231_fjy(%esp),%xmm1
        addps  nb231_fjz(%esp),%xmm2
        movaps %xmm0,nb231_fjx(%esp)
        movaps %xmm1,nb231_fjy(%esp)
        movaps %xmm2,nb231_fjz(%esp)

        ## H2 interactions 
        movaps nb231_rinvH2(%esp),%xmm2
        movaps nb231_krsqH2(%esp),%xmm1
        movaps %xmm2,%xmm3
        addps  %xmm1,%xmm2 ## rinv+krsq
        subps  nb231_crf(%esp),%xmm2   ## rinv+krsq-crf
        movaps %xmm3,%xmm0
        subps  %xmm1,%xmm3
        subps  %xmm1,%xmm3 ## rinv-2*krsq
        mulps  %xmm0,%xmm0
        mulps  nb231_qqH(%esp),%xmm2
        mulps  nb231_qqH(%esp),%xmm3

        mulps  %xmm0,%xmm3
        addps  nb231_vctot(%esp),%xmm2
        movaps %xmm2,nb231_vctot(%esp)

        movaps nb231_dxH2(%esp),%xmm0
        movaps nb231_dyH2(%esp),%xmm1
        movaps nb231_dzH2(%esp),%xmm2

        mulps  %xmm3,%xmm0
        mulps  %xmm3,%xmm1
        mulps  %xmm3,%xmm2
        ## update H2 forces 
        movaps nb231_fixH2(%esp),%xmm3
        movaps nb231_fiyH2(%esp),%xmm4
        movaps nb231_fizH2(%esp),%xmm7
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm7
        movaps %xmm3,nb231_fixH2(%esp)
        movaps %xmm4,nb231_fiyH2(%esp)
        movaps %xmm7,nb231_fizH2(%esp)

        movl nb231_faction(%ebp),%edi
        ## update j forces 
        addps nb231_fjx(%esp),%xmm0
        addps nb231_fjy(%esp),%xmm1
        addps nb231_fjz(%esp),%xmm2

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
        subl $4,nb231_innerk(%esp)
        jl    _nb_kernel231_ia32_sse.nb231_odd_inner
        jmp   _nb_kernel231_ia32_sse.nb231_unroll_loop
_nb_kernel231_ia32_sse.nb231_odd_inner: 
        addl $4,nb231_innerk(%esp)
        jnz   _nb_kernel231_ia32_sse.nb231_odd_loop
        jmp   _nb_kernel231_ia32_sse.nb231_updateouterdata
_nb_kernel231_ia32_sse.nb231_odd_loop: 
        movl  nb231_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb231_innerjjnr(%esp)

        xorps %xmm4,%xmm4
        movss nb231_iqO(%esp),%xmm4
        movl nb231_charge(%ebp),%esi
        movhps nb231_iqH(%esp),%xmm4
        movss (%esi,%eax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb231_qqO(%esp)    ## use oxygen qq for storage 

        xorps %xmm6,%xmm6
        movl nb231_type(%ebp),%esi
        movl (%esi,%eax,4),%ebx
        movl nb231_vdwparam(%ebp),%esi
        shll %ebx
        addl nb231_ntia(%esp),%ebx
        movlps (%esi,%ebx,4),%xmm6
        movaps %xmm6,%xmm7
        shufps $252,%xmm6,%xmm6 ## constant 11111100
        shufps $253,%xmm7,%xmm7 ## constant 11111101
        movaps %xmm6,nb231_c6(%esp)
        movaps %xmm7,nb231_c12(%esp)

        movl nb231_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## move j coords to xmm0-xmm2 
        movss (%esi,%eax,4),%xmm0
        movss 4(%esi,%eax,4),%xmm1
        movss 8(%esi,%eax,4),%xmm2
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movss nb231_ixO(%esp),%xmm3
        movss nb231_iyO(%esp),%xmm4
        movss nb231_izO(%esp),%xmm5

        movlps nb231_ixH1(%esp),%xmm6
        movlps nb231_ixH2(%esp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm3
        movlps nb231_iyH1(%esp),%xmm6
        movlps nb231_iyH2(%esp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm4
        movlps nb231_izH1(%esp),%xmm6
        movlps nb231_izH2(%esp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm5

        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        movaps %xmm3,nb231_dxO(%esp)
        movaps %xmm4,nb231_dyO(%esp)
        movaps %xmm5,nb231_dzO(%esp)

        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5

        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        ## rsq in xmm4 

        movaps %xmm4,%xmm0
        mulps nb231_krf(%esp),%xmm0
        movaps %xmm0,nb231_krsqO(%esp)

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb231_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb231_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        ## a little trick to avoid NaNs: 
        ## positions 0,2,and 3 are valid, but not 1. 
        ## If it contains NaN it doesnt help to mult by 0, 
        ## So we shuffle it and copy pos 0 to pos1! 
        shufps $224,%xmm0,%xmm0 ## constant 11100000

        ## rsq still in xmm4, rinv in xmm0.
        mulps  %xmm0,%xmm4
        mulps  nb231_tsc(%esp),%xmm4   ## rtab

        cvttps2pi %xmm4,%mm6
        cvtpi2ps %mm6,%xmm6
        subss  %xmm6,%xmm4
        movss %xmm4,%xmm1       ## xmm1=eps 
        movss %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $3,%mm6

        movd %eax,%mm0

        movl nb231_VFtab(%ebp),%esi
        movd %mm6,%eax

        ## dispersion 
        movlps (%esi,%eax,4),%xmm5
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## dispersion table ready, in xmm4-xmm7         

        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp      
        mulss  nb231_two(%esp),%xmm7    ## two*Heps2 
        addss  %xmm6,%xmm7
        addss  %xmm5,%xmm7 ## xmm7=FF 
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 

        movss nb231_c6(%esp),%xmm4
        mulss  %xmm4,%xmm7       ## fijD 
        mulss  %xmm4,%xmm5       ## Vvdw6 
        mulss  nb231_tsc(%esp),%xmm7
        ## put scalar force on stack Update Vvdwtot directly 
        addss  nb231_Vvdwtot(%esp),%xmm5
        movss %xmm7,nb231_fstmp(%esp)
        movss %xmm5,nb231_Vvdwtot(%esp)

        ## repulsion 
        movlps 16(%esi,%eax,4),%xmm5
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 24(%esi,%eax,4),%xmm7
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## table ready, in xmm4-xmm7    
        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp      
        mulss  nb231_two(%esp),%xmm7    ## two*Heps2 
        addss  %xmm6,%xmm7
        addss  %xmm5,%xmm7 ## xmm7=FF 
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 

        movss nb231_c12(%esp),%xmm4
        mulss  %xmm4,%xmm7 ## fijR 
        mulss  %xmm4,%xmm5 ## Vvdw12 
        mulss  nb231_tsc(%esp),%xmm7
        addss  nb231_fstmp(%esp),%xmm7
        movss %xmm7,nb231_fstmp(%esp)

        addss  nb231_Vvdwtot(%esp),%xmm5
        movss %xmm5,nb231_Vvdwtot(%esp)

        movaps %xmm0,%xmm1      ## xmm1=rinv 
        movaps %xmm0,%xmm4
        movaps nb231_krsqO(%esp),%xmm3
        addps  %xmm3,%xmm0      ## xmm0=rinv+ krsq 
        mulps  nb231_two(%esp),%xmm3
        subps  nb231_crf(%esp),%xmm0   ## xmm0=rinv+ krsq-crf 
        subps  %xmm3,%xmm1      ## xmm1=rinv-2*krsq 
        mulps  nb231_qqO(%esp),%xmm0    ## xmm0=vcoul 
        mulps  nb231_qqO(%esp),%xmm1    ## xmm1=coul part of fs 

        mulps  %xmm4,%xmm1
        subss  nb231_fstmp(%esp),%xmm1
        mulps  %xmm1,%xmm4

        addps  nb231_vctot(%esp),%xmm0
        movaps %xmm0,nb231_vctot(%esp)

        movaps nb231_dxO(%esp),%xmm0
        movaps nb231_dyO(%esp),%xmm1
        movaps nb231_dzO(%esp),%xmm2

        movd %mm0,%eax

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        movss  nb231_fixO(%esp),%xmm3
        movss  nb231_fiyO(%esp),%xmm4
        movss  nb231_fizO(%esp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,nb231_fixO(%esp)
        movss  %xmm4,nb231_fiyO(%esp)
        movss  %xmm5,nb231_fizO(%esp)   ## updated the O force now do the H's 
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        shufps $230,%xmm3,%xmm3 ## constant 11100110     ;# shift right 
        shufps $230,%xmm4,%xmm4 ## constant 11100110
        shufps $230,%xmm5,%xmm5 ## constant 11100110
        addss  nb231_fixH1(%esp),%xmm3
        addss  nb231_fiyH1(%esp),%xmm4
        addss  nb231_fizH1(%esp),%xmm5
        movss  %xmm3,nb231_fixH1(%esp)
        movss  %xmm4,nb231_fiyH1(%esp)
        movss  %xmm5,nb231_fizH1(%esp)          ## updated the H1 force 

        movl nb231_faction(%ebp),%edi
        shufps $231,%xmm3,%xmm3 ## constant 11100111     ;# shift right 
        shufps $231,%xmm4,%xmm4 ## constant 11100111
        shufps $231,%xmm5,%xmm5 ## constant 11100111
        addss  nb231_fixH2(%esp),%xmm3
        addss  nb231_fiyH2(%esp),%xmm4
        addss  nb231_fizH2(%esp),%xmm5
        movss  %xmm3,nb231_fixH2(%esp)
        movss  %xmm4,nb231_fiyH2(%esp)
        movss  %xmm5,nb231_fizH2(%esp)          ## updated the H2 force 

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

        decl nb231_innerk(%esp)
        jz    _nb_kernel231_ia32_sse.nb231_updateouterdata
        jmp   _nb_kernel231_ia32_sse.nb231_odd_loop
_nb_kernel231_ia32_sse.nb231_updateouterdata: 
        movl  nb231_ii3(%esp),%ecx
        movl  nb231_faction(%ebp),%edi
        movl  nb231_fshift(%ebp),%esi
        movl  nb231_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movaps nb231_fixO(%esp),%xmm0
        movaps nb231_fiyO(%esp),%xmm1
        movaps nb231_fizO(%esp),%xmm2

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
        movaps nb231_fixH1(%esp),%xmm0
        movaps nb231_fiyH1(%esp),%xmm1
        movaps nb231_fizH1(%esp),%xmm2

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
        movaps nb231_fixH2(%esp),%xmm0
        movaps nb231_fiyH2(%esp),%xmm1
        movaps nb231_fizH2(%esp),%xmm2

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
        movl nb231_n(%esp),%esi
        ## get group index for i particle 
        movl  nb231_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb231_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb231_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb231_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb231_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb231_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel231_ia32_sse.nb231_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb231_n(%esp)
        jmp _nb_kernel231_ia32_sse.nb231_outer
_nb_kernel231_ia32_sse.nb231_outerend: 
        ## check if more outer neighborlists remain
        movl  nb231_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel231_ia32_sse.nb231_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel231_ia32_sse.nb231_threadloop
_nb_kernel231_ia32_sse.nb231_end: 
        emms

        movl nb231_nouter(%esp),%eax
        movl nb231_ninner(%esp),%ebx
        movl nb231_outeriter(%ebp),%ecx
        movl nb231_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb231_salign(%esp),%eax
        addl %eax,%esp
        addl $860,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret



.globl nb_kernel231nf_ia32_sse
.globl _nb_kernel231nf_ia32_sse
nb_kernel231nf_ia32_sse:        
_nb_kernel231nf_ia32_sse:       
.set nb231nf_p_nri, 8
.set nb231nf_iinr, 12
.set nb231nf_jindex, 16
.set nb231nf_jjnr, 20
.set nb231nf_shift, 24
.set nb231nf_shiftvec, 28
.set nb231nf_fshift, 32
.set nb231nf_gid, 36
.set nb231nf_pos, 40
.set nb231nf_faction, 44
.set nb231nf_charge, 48
.set nb231nf_p_facel, 52
.set nb231nf_argkrf, 56
.set nb231nf_argcrf, 60
.set nb231nf_Vc, 64
.set nb231nf_type, 68
.set nb231nf_p_ntype, 72
.set nb231nf_vdwparam, 76
.set nb231nf_Vvdw, 80
.set nb231nf_p_tabscale, 84
.set nb231nf_VFtab, 88
.set nb231nf_invsqrta, 92
.set nb231nf_dvda, 96
.set nb231nf_p_gbtabscale, 100
.set nb231nf_GBtab, 104
.set nb231nf_p_nthreads, 108
.set nb231nf_count, 112
.set nb231nf_mtx, 116
.set nb231nf_outeriter, 120
.set nb231nf_inneriter, 124
.set nb231nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb231nf_ixO, 0
.set nb231nf_iyO, 16
.set nb231nf_izO, 32
.set nb231nf_ixH1, 48
.set nb231nf_iyH1, 64
.set nb231nf_izH1, 80
.set nb231nf_ixH2, 96
.set nb231nf_iyH2, 112
.set nb231nf_izH2, 128
.set nb231nf_iqO, 144
.set nb231nf_iqH, 160
.set nb231nf_qqO, 176
.set nb231nf_qqH, 192
.set nb231nf_c6, 208
.set nb231nf_c12, 224
.set nb231nf_tsc, 240
.set nb231nf_vctot, 256
.set nb231nf_Vvdwtot, 272
.set nb231nf_half, 288
.set nb231nf_three, 304
.set nb231nf_krf, 320
.set nb231nf_crf, 336
.set nb231nf_krsqO, 352
.set nb231nf_krsqH1, 368
.set nb231nf_krsqH2, 384
.set nb231nf_rinvO, 400
.set nb231nf_rinvH1, 416
.set nb231nf_rinvH2, 432
.set nb231nf_is3, 448
.set nb231nf_ii3, 452
.set nb231nf_ntia, 456
.set nb231nf_innerjjnr, 460
.set nb231nf_innerk, 464
.set nb231nf_n, 468
.set nb231nf_nn1, 472
.set nb231nf_nri, 476
.set nb231nf_nouter, 480
.set nb231nf_ninner, 484
.set nb231nf_salign, 488
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $492,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb231nf_salign(%esp)


        emms

        ## Move args passed by reference to stack
        movl nb231nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb231nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb231nf_nouter(%esp)
        movl %eax,nb231nf_ninner(%esp)

        movl nb231nf_p_tabscale(%ebp),%eax
        movss (%eax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb231nf_tsc(%esp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb231nf_half(%esp)
        movss nb231nf_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb231nf_half(%esp)
        movaps %xmm3,nb231nf_three(%esp)

        movl nb231nf_argkrf(%ebp),%esi
        movl nb231nf_argcrf(%ebp),%edi
        movss (%esi),%xmm5
        movss (%edi),%xmm6
        shufps $0,%xmm5,%xmm5
        shufps $0,%xmm6,%xmm6
        movaps %xmm5,nb231nf_krf(%esp)
        movaps %xmm6,nb231nf_crf(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb231nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb231nf_charge(%ebp),%edx
        movss (%edx,%ebx,4),%xmm3
        movss 4(%edx,%ebx,4),%xmm4
        movl nb231nf_p_facel(%ebp),%esi
        movss (%esi),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb231nf_iqO(%esp)
        movaps %xmm4,nb231nf_iqH(%esp)

        movl  nb231nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb231nf_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb231nf_ntia(%esp)

_nb_kernel231nf_ia32_sse.nb231nf_threadloop: 
        movl  nb231nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel231nf_ia32_sse.nb231nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel231nf_ia32_sse.nb231nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb231nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb231nf_n(%esp)
        movl %ebx,nb231nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel231nf_ia32_sse.nb231nf_outerstart
        jmp _nb_kernel231nf_ia32_sse.nb231nf_end

_nb_kernel231nf_ia32_sse.nb231nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb231nf_nouter(%esp),%ebx
        movl %ebx,nb231nf_nouter(%esp)

_nb_kernel231nf_ia32_sse.nb231nf_outer: 
        movl  nb231nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb231nf_is3(%esp)            ## store is3 

        movl  nb231nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb231nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb231nf_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb231nf_ii3(%esp)

        addss (%eax,%ebx,4),%xmm3
        addss 4(%eax,%ebx,4),%xmm4
        addss 8(%eax,%ebx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb231nf_ixO(%esp)
        movaps %xmm4,nb231nf_iyO(%esp)
        movaps %xmm5,nb231nf_izO(%esp)

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
        movaps %xmm0,nb231nf_ixH1(%esp)
        movaps %xmm1,nb231nf_iyH1(%esp)
        movaps %xmm2,nb231nf_izH1(%esp)
        movaps %xmm3,nb231nf_ixH2(%esp)
        movaps %xmm4,nb231nf_iyH2(%esp)
        movaps %xmm5,nb231nf_izH2(%esp)

        ## clear vctot
        xorps %xmm4,%xmm4
        movaps %xmm4,nb231nf_vctot(%esp)
        movaps %xmm4,nb231nf_Vvdwtot(%esp)

        movl  nb231nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb231nf_pos(%ebp),%esi
        movl  nb231nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb231nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb231nf_ninner(%esp),%ecx
        movl  %ecx,nb231nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb231nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel231nf_ia32_sse.nb231nf_unroll_loop
        jmp   _nb_kernel231nf_ia32_sse.nb231nf_odd_inner
_nb_kernel231nf_ia32_sse.nb231nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb231nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 

        addl $16,nb231nf_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb231nf_charge(%ebp),%esi     ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4           ## and in xmm4 
        mulps  nb231nf_iqO(%esp),%xmm3
        mulps  nb231nf_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps  %xmm3,nb231nf_qqO(%esp)
        movaps  %xmm4,nb231nf_qqH(%esp)

        movl nb231nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb231nf_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb231nf_ntia(%esp),%edi
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

        movaps %xmm4,nb231nf_c6(%esp)
        movaps %xmm6,nb231nf_c12(%esp)

        movl nb231nf_pos(%ebp),%esi        ## base of pos[] 

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
        movaps nb231nf_ixO(%esp),%xmm4
        movaps nb231nf_iyO(%esp),%xmm5
        movaps nb231nf_izO(%esp),%xmm6

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
        movaps nb231nf_ixH1(%esp),%xmm4
        movaps nb231nf_iyH1(%esp),%xmm5
        movaps nb231nf_izH1(%esp),%xmm6

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
        movaps nb231nf_ixH2(%esp),%xmm3
        movaps nb231nf_iyH2(%esp),%xmm4
        movaps nb231nf_izH2(%esp),%xmm5

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


        movaps %xmm5,%xmm0
        movaps %xmm6,%xmm1
        movaps %xmm7,%xmm2

        mulps  nb231nf_krf(%esp),%xmm0
        mulps  nb231nf_krf(%esp),%xmm1
        mulps  nb231nf_krf(%esp),%xmm2

        movaps %xmm0,nb231nf_krsqH2(%esp)
        movaps %xmm1,nb231nf_krsqH1(%esp)
        movaps %xmm2,nb231nf_krsqO(%esp)

        ## start with rsqO (still in xmm7) - seed to xmm2       
        rsqrtps %xmm7,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb231nf_three(%esp),%xmm4
        mulps   %xmm7,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb231nf_half(%esp),%xmm4
        movaps  %xmm4,nb231nf_rinvO(%esp)       ## rinvO in xmm0, rsqO in xmm7

        ## rsqH1 - seed in xmm2 
        rsqrtps %xmm6,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb231nf_three(%esp),%xmm4
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb231nf_half(%esp),%xmm4
        movaps  %xmm4,nb231nf_rinvH1(%esp)      ## rinvH1 in xmm6 

        ## rsqH2 - seed in xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb231nf_three(%esp),%xmm4
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb231nf_half(%esp),%xmm4
        movaps  %xmm4,nb231nf_rinvH2(%esp)      ## rinvH2 in xmm5 


        ## do O table interactions - rsqO in xmm7.
        mulps nb231nf_rinvO(%esp),%xmm7
        mulps nb231nf_tsc(%esp),%xmm7   ## rtab

        movhlps %xmm7,%xmm5
        cvttps2pi %xmm7,%mm6
        cvttps2pi %xmm5,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        cvtpi2ps %mm7,%xmm5
        movlhps %xmm5,%xmm6
        subps  %xmm6,%xmm7
        movaps %xmm7,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $3,%mm6
        pslld $3,%mm7

        movl nb231nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        ## dispersion 
        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half dispersion table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of dispersion table 
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

        movaps nb231nf_c6(%esp),%xmm4
        mulps  %xmm4,%xmm5       ## Vvdw6 

        addps  nb231nf_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb231nf_Vvdwtot(%esp)

        ## repulsion 
        movlps 16(%esi,%eax,4),%xmm5
        movlps 16(%esi,%ecx,4),%xmm7
        movhps 16(%esi,%ebx,4),%xmm5
        movhps 16(%esi,%edx,4),%xmm7    ## got half repulsion table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 24(%esi,%eax,4),%xmm7
        movlps 24(%esi,%ecx,4),%xmm3
        movhps 24(%esi,%ebx,4),%xmm7
        movhps 24(%esi,%edx,4),%xmm3    ## other half of repulsion table 
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## table ready, in xmm4-xmm7    
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb231nf_c12(%esp),%xmm4
        mulps  %xmm4,%xmm5 ## Vvdw12 

        addps  nb231nf_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb231nf_Vvdwtot(%esp)

        movaps nb231nf_rinvO(%esp),%xmm2
        movaps nb231nf_krsqO(%esp),%xmm1
        addps  %xmm1,%xmm2 ## rinv+krsq
        subps  nb231nf_crf(%esp),%xmm2   ## rinv+krsq-crf
        mulps  nb231nf_qqO(%esp),%xmm2

        addps  nb231nf_vctot(%esp),%xmm2
        movaps %xmm2,nb231nf_vctot(%esp)

        ## H1 interactions 
        movaps nb231nf_rinvH1(%esp),%xmm2
        movaps nb231nf_krsqH1(%esp),%xmm1
        addps  %xmm1,%xmm2 ## rinv+krsq
        subps  nb231nf_crf(%esp),%xmm2   ## rinv+krsq-crf
        mulps  nb231nf_qqH(%esp),%xmm2

        addps  nb231nf_vctot(%esp),%xmm2
        movaps %xmm2,nb231nf_vctot(%esp)

        ## H2 interactions 
        movaps nb231nf_rinvH2(%esp),%xmm2
        movaps nb231nf_krsqH2(%esp),%xmm1
        addps  %xmm1,%xmm2 ## rinv+krsq
        subps  nb231nf_crf(%esp),%xmm2   ## rinv+krsq-crf
        mulps  nb231nf_qqH(%esp),%xmm2
        addps  nb231nf_vctot(%esp),%xmm2
        movaps %xmm2,nb231nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $4,nb231nf_innerk(%esp)
        jl    _nb_kernel231nf_ia32_sse.nb231nf_odd_inner
        jmp   _nb_kernel231nf_ia32_sse.nb231nf_unroll_loop
_nb_kernel231nf_ia32_sse.nb231nf_odd_inner: 
        addl $4,nb231nf_innerk(%esp)
        jnz   _nb_kernel231nf_ia32_sse.nb231nf_odd_loop
        jmp   _nb_kernel231nf_ia32_sse.nb231nf_updateouterdata
_nb_kernel231nf_ia32_sse.nb231nf_odd_loop: 
        movl  nb231nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb231nf_innerjjnr(%esp)

        xorps %xmm4,%xmm4
        movss nb231nf_iqO(%esp),%xmm4
        movl nb231nf_charge(%ebp),%esi
        movhps nb231nf_iqH(%esp),%xmm4
        movss (%esi,%eax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb231nf_qqO(%esp)          ## use oxygen qq for storage 

        xorps %xmm6,%xmm6
        movl nb231nf_type(%ebp),%esi
        movl (%esi,%eax,4),%ebx
        movl nb231nf_vdwparam(%ebp),%esi
        shll %ebx
        addl nb231nf_ntia(%esp),%ebx
        movlps (%esi,%ebx,4),%xmm6
        movaps %xmm6,%xmm7
        shufps $252,%xmm6,%xmm6 ## constant 11111100
        shufps $253,%xmm7,%xmm7 ## constant 11111101
        movaps %xmm6,nb231nf_c6(%esp)
        movaps %xmm7,nb231nf_c12(%esp)

        movl nb231nf_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## move j coords to xmm0-xmm2 
        movss (%esi,%eax,4),%xmm0
        movss 4(%esi,%eax,4),%xmm1
        movss 8(%esi,%eax,4),%xmm2
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movss nb231nf_ixO(%esp),%xmm3
        movss nb231nf_iyO(%esp),%xmm4
        movss nb231nf_izO(%esp),%xmm5

        movlps nb231nf_ixH1(%esp),%xmm6
        movlps nb231nf_ixH2(%esp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm3
        movlps nb231nf_iyH1(%esp),%xmm6
        movlps nb231nf_iyH2(%esp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm4
        movlps nb231nf_izH1(%esp),%xmm6
        movlps nb231nf_izH2(%esp),%xmm7
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

        movaps %xmm4,%xmm0
        mulps nb231nf_krf(%esp),%xmm0
        movaps %xmm0,nb231nf_krsqO(%esp)

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb231nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb231nf_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        ## a little trick to avoid NaNs: 
        ## positions 0,2,and 3 are valid, but not 1. 
        ## If it contains NaN it doesnt help to mult by 0, 
        ## So we shuffle it and copy pos 0 to pos1! 
        shufps $224,%xmm0,%xmm0 ## constant 11100000

        ## rsq still in xmm4, rinv in xmm0.
        mulps  %xmm0,%xmm4
        mulps  nb231nf_tsc(%esp),%xmm4   ## rtab

        cvttps2pi %xmm4,%mm6
        cvtpi2ps %mm6,%xmm6
        subss  %xmm6,%xmm4
        movss %xmm4,%xmm1       ## xmm1=eps 
        movss %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $3,%mm6

        movl nb231nf_VFtab(%ebp),%esi
        movd %mm6,%eax

        ## dispersion 
        movlps (%esi,%eax,4),%xmm5
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## dispersion table ready, in xmm4-xmm7         

        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp      
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 

        movss nb231nf_c6(%esp),%xmm4
        mulss  %xmm4,%xmm5       ## Vvdw6 

        addss  nb231nf_Vvdwtot(%esp),%xmm5
        movss %xmm5,nb231nf_Vvdwtot(%esp)

        ## repulsion 
        movlps 16(%esi,%eax,4),%xmm5
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 24(%esi,%eax,4),%xmm7
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## table ready, in xmm4-xmm7    
        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp      
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 

        movss nb231nf_c12(%esp),%xmm4
        mulss  %xmm4,%xmm5 ## Vvdw12 
        addss  nb231nf_Vvdwtot(%esp),%xmm5
        movss %xmm5,nb231nf_Vvdwtot(%esp)

        movaps %xmm0,%xmm1      ## xmm1=rinv 
        movaps %xmm0,%xmm4
        movaps nb231nf_krsqO(%esp),%xmm3
        addps  %xmm3,%xmm0      ## xmm0=rinv+ krsq 

        subps  nb231nf_crf(%esp),%xmm0   ## xmm0=rinv+ krsq-crf 

        mulps  nb231nf_qqO(%esp),%xmm0          ## xmm0=vcoul 

        addps  nb231nf_vctot(%esp),%xmm0
        movaps %xmm0,nb231nf_vctot(%esp)

        decl nb231nf_innerk(%esp)
        jz    _nb_kernel231nf_ia32_sse.nb231nf_updateouterdata
        jmp   _nb_kernel231nf_ia32_sse.nb231nf_odd_loop
_nb_kernel231nf_ia32_sse.nb231nf_updateouterdata: 
        ## get n from stack
        movl nb231nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb231nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb231nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb231nf_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb231nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb231nf_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb231nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel231nf_ia32_sse.nb231nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb231nf_n(%esp)
        jmp _nb_kernel231nf_ia32_sse.nb231nf_outer
_nb_kernel231nf_ia32_sse.nb231nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb231nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel231nf_ia32_sse.nb231nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel231nf_ia32_sse.nb231nf_threadloop
_nb_kernel231nf_ia32_sse.nb231nf_end: 
        emms

        movl nb231nf_nouter(%esp),%eax
        movl nb231nf_ninner(%esp),%ebx
        movl nb231nf_outeriter(%ebp),%ecx
        movl nb231nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb231nf_salign(%esp),%eax
        addl %eax,%esp
        addl $492,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret



