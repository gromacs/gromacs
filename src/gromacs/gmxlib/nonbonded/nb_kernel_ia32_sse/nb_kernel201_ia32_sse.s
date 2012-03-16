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




.globl nb_kernel201_ia32_sse
.globl _nb_kernel201_ia32_sse
nb_kernel201_ia32_sse:  
_nb_kernel201_ia32_sse: 
.set nb201_p_nri, 8
.set nb201_iinr, 12
.set nb201_jindex, 16
.set nb201_jjnr, 20
.set nb201_shift, 24
.set nb201_shiftvec, 28
.set nb201_fshift, 32
.set nb201_gid, 36
.set nb201_pos, 40
.set nb201_faction, 44
.set nb201_charge, 48
.set nb201_p_facel, 52
.set nb201_argkrf, 56
.set nb201_argcrf, 60
.set nb201_Vc, 64
.set nb201_type, 68
.set nb201_p_ntype, 72
.set nb201_vdwparam, 76
.set nb201_Vvdw, 80
.set nb201_p_tabscale, 84
.set nb201_VFtab, 88
.set nb201_invsqrta, 92
.set nb201_dvda, 96
.set nb201_p_gbtabscale, 100
.set nb201_GBtab, 104
.set nb201_p_nthreads, 108
.set nb201_count, 112
.set nb201_mtx, 116
.set nb201_outeriter, 120
.set nb201_inneriter, 124
.set nb201_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb201_ixO, 0
.set nb201_iyO, 16
.set nb201_izO, 32
.set nb201_ixH1, 48
.set nb201_iyH1, 64
.set nb201_izH1, 80
.set nb201_ixH2, 96
.set nb201_iyH2, 112
.set nb201_izH2, 128
.set nb201_iqO, 144
.set nb201_iqH, 160
.set nb201_dxO, 176
.set nb201_dyO, 192
.set nb201_dzO, 208
.set nb201_dxH1, 224
.set nb201_dyH1, 240
.set nb201_dzH1, 256
.set nb201_dxH2, 272
.set nb201_dyH2, 288
.set nb201_dzH2, 304
.set nb201_qqO, 320
.set nb201_qqH, 336
.set nb201_vctot, 352
.set nb201_fixO, 384
.set nb201_fiyO, 400
.set nb201_fizO, 416
.set nb201_fixH1, 432
.set nb201_fiyH1, 448
.set nb201_fizH1, 464
.set nb201_fixH2, 480
.set nb201_fiyH2, 496
.set nb201_fizH2, 512
.set nb201_fjx, 528
.set nb201_fjy, 544
.set nb201_fjz, 560
.set nb201_half, 576
.set nb201_three, 592
.set nb201_two, 608
.set nb201_krf, 624
.set nb201_crf, 640
.set nb201_krsqO, 656
.set nb201_krsqH1, 672
.set nb201_krsqH2, 688
.set nb201_is3, 704
.set nb201_ii3, 708
.set nb201_innerjjnr, 712
.set nb201_innerk, 716
.set nb201_n, 720
.set nb201_nn1, 724
.set nb201_nri, 728
.set nb201_nouter, 732
.set nb201_ninner, 736
.set nb201_salign, 740
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $744,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb201_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb201_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb201_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb201_nouter(%esp)
        movl %eax,nb201_ninner(%esp)


        movl nb201_argkrf(%ebp),%esi
        movl nb201_argcrf(%ebp),%edi
        movss (%esi),%xmm5
        movss (%edi),%xmm6
        shufps $0,%xmm5,%xmm5
        shufps $0,%xmm6,%xmm6
        movaps %xmm5,nb201_krf(%esp)
        movaps %xmm6,nb201_crf(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb201_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb201_charge(%ebp),%edx
        movss (%edx,%ebx,4),%xmm3
        movss 4(%edx,%ebx,4),%xmm4
        movl nb201_p_facel(%ebp),%esi
        movss (%esi),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb201_iqO(%esp)
        movaps %xmm4,nb201_iqH(%esp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb201_half(%esp)
        movss nb201_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb201_half(%esp)
        movaps %xmm2,nb201_two(%esp)
        movaps %xmm3,nb201_three(%esp)


_nb_kernel201_ia32_sse.nb201_threadloop: 
        movl  nb201_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel201_ia32_sse.nb201_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel201_ia32_sse.nb201_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb201_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb201_n(%esp)
        movl %ebx,nb201_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel201_ia32_sse.nb201_outerstart
        jmp _nb_kernel201_ia32_sse.nb201_end

_nb_kernel201_ia32_sse.nb201_outerstart: 
        ## ebx contains number of outer iterations
        addl nb201_nouter(%esp),%ebx
        movl %ebx,nb201_nouter(%esp)

_nb_kernel201_ia32_sse.nb201_outer: 
        movl  nb201_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb201_is3(%esp)      ## store is3 

        movl  nb201_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb201_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb201_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb201_ii3(%esp)

        addss (%eax,%ebx,4),%xmm3
        addss 4(%eax,%ebx,4),%xmm4
        addss 8(%eax,%ebx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb201_ixO(%esp)
        movaps %xmm4,nb201_iyO(%esp)
        movaps %xmm5,nb201_izO(%esp)

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
        movaps %xmm0,nb201_ixH1(%esp)
        movaps %xmm1,nb201_iyH1(%esp)
        movaps %xmm2,nb201_izH1(%esp)
        movaps %xmm3,nb201_ixH2(%esp)
        movaps %xmm4,nb201_iyH2(%esp)
        movaps %xmm5,nb201_izH2(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb201_vctot(%esp)
        movaps %xmm4,nb201_fixO(%esp)
        movaps %xmm4,nb201_fiyO(%esp)
        movaps %xmm4,nb201_fizO(%esp)
        movaps %xmm4,nb201_fixH1(%esp)
        movaps %xmm4,nb201_fiyH1(%esp)
        movaps %xmm4,nb201_fizH1(%esp)
        movaps %xmm4,nb201_fixH2(%esp)
        movaps %xmm4,nb201_fiyH2(%esp)
        movaps %xmm4,nb201_fizH2(%esp)

        movl  nb201_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb201_pos(%ebp),%esi
        movl  nb201_faction(%ebp),%edi
        movl  nb201_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb201_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb201_ninner(%esp),%ecx
        movl  %ecx,nb201_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb201_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel201_ia32_sse.nb201_unroll_loop
        jmp   _nb_kernel201_ia32_sse.nb201_odd_inner
_nb_kernel201_ia32_sse.nb201_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb201_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 

        addl $16,nb201_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb201_charge(%ebp),%esi     ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4           ## and in xmm4 
        mulps  nb201_iqO(%esp),%xmm3
        mulps  nb201_iqH(%esp),%xmm4

        movaps  %xmm3,nb201_qqO(%esp)
        movaps  %xmm4,nb201_qqH(%esp)

        movl nb201_pos(%ebp),%esi        ## base of pos[] 

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
        movaps nb201_ixO(%esp),%xmm4
        movaps nb201_iyO(%esp),%xmm5
        movaps nb201_izO(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb201_dxO(%esp)
        movaps %xmm5,nb201_dyO(%esp)
        movaps %xmm6,nb201_dzO(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        movaps %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movaps nb201_ixH1(%esp),%xmm4
        movaps nb201_iyH1(%esp),%xmm5
        movaps nb201_izH1(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb201_dxH1(%esp)
        movaps %xmm5,nb201_dyH1(%esp)
        movaps %xmm6,nb201_dzH1(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm6
        addps %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movaps nb201_ixH2(%esp),%xmm3
        movaps nb201_iyH2(%esp),%xmm4
        movaps nb201_izH2(%esp),%xmm5

        ## calc dr 
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        ## store dr 
        movaps %xmm3,nb201_dxH2(%esp)
        movaps %xmm4,nb201_dyH2(%esp)
        movaps %xmm5,nb201_dzH2(%esp)
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

        mulps  nb201_krf(%esp),%xmm0
        mulps  nb201_krf(%esp),%xmm1
        mulps  nb201_krf(%esp),%xmm2

        movaps %xmm0,nb201_krsqH2(%esp)
        movaps %xmm1,nb201_krsqH1(%esp)
        movaps %xmm2,nb201_krsqO(%esp)

        ## start with rsqO - seed in xmm2       
        rsqrtps %xmm7,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb201_three(%esp),%xmm4
        mulps   %xmm7,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb201_half(%esp),%xmm4
        movaps  %xmm4,%xmm7     ## rinvO in xmm7 
        ## rsqH1 - seed in xmm2 
        rsqrtps %xmm6,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb201_three(%esp),%xmm4
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb201_half(%esp),%xmm4
        movaps  %xmm4,%xmm6     ## rinvH1 in xmm6 
        ## rsqH2 - seed in xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb201_three(%esp),%xmm4
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb201_half(%esp),%xmm4
        movaps  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do O interactions 
        movaps  %xmm7,%xmm4
        mulps   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq 

        movaps %xmm7,%xmm0
        movaps nb201_krsqO(%esp),%xmm1
        addps  %xmm1,%xmm0
        subps  nb201_crf(%esp),%xmm0   ## xmm0=rinv+ krsq-crf 
        mulps  nb201_two(%esp),%xmm1
        subps  %xmm1,%xmm7
        mulps  nb201_qqO(%esp),%xmm0
        mulps  nb201_qqO(%esp),%xmm7

        mulps  %xmm7,%xmm4      ## total fsO in xmm4 

        addps  nb201_vctot(%esp),%xmm0
        movaps %xmm0,nb201_vctot(%esp)

        movaps nb201_dxO(%esp),%xmm0
        movaps nb201_dyO(%esp),%xmm1
        movaps nb201_dzO(%esp),%xmm2
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2

        ## update O forces 
        movaps nb201_fixO(%esp),%xmm3
        movaps nb201_fiyO(%esp),%xmm4
        movaps nb201_fizO(%esp),%xmm7
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm7
        movaps %xmm3,nb201_fixO(%esp)
        movaps %xmm4,nb201_fiyO(%esp)
        movaps %xmm7,nb201_fizO(%esp)
        ## update j forces with water O 
        movaps %xmm0,nb201_fjx(%esp)
        movaps %xmm1,nb201_fjy(%esp)
        movaps %xmm2,nb201_fjz(%esp)

        ## H1 interactions 
        movaps  %xmm6,%xmm4
        mulps   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        movaps  %xmm6,%xmm7
        movaps  nb201_krsqH1(%esp),%xmm0
        addps   %xmm0,%xmm6     ## xmm6=rinv+ krsq 
        subps   nb201_crf(%esp),%xmm6   ## xmm6=rinv+ krsq-crf 
        mulps   nb201_two(%esp),%xmm0
        subps   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulps   nb201_qqH(%esp),%xmm6   ## vcoul 
        mulps   nb201_qqH(%esp),%xmm7
        mulps  %xmm7,%xmm4              ## total fsH1 in xmm4 

        addps  nb201_vctot(%esp),%xmm6

        movaps nb201_dxH1(%esp),%xmm0
        movaps nb201_dyH1(%esp),%xmm1
        movaps nb201_dzH1(%esp),%xmm2
        movaps %xmm6,nb201_vctot(%esp)
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2

        ## update H1 forces 
        movaps nb201_fixH1(%esp),%xmm3
        movaps nb201_fiyH1(%esp),%xmm4
        movaps nb201_fizH1(%esp),%xmm7
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm7
        movaps %xmm3,nb201_fixH1(%esp)
        movaps %xmm4,nb201_fiyH1(%esp)
        movaps %xmm7,nb201_fizH1(%esp)
        ## update j forces with water H1 
        addps  nb201_fjx(%esp),%xmm0
        addps  nb201_fjy(%esp),%xmm1
        addps  nb201_fjz(%esp),%xmm2
        movaps %xmm0,nb201_fjx(%esp)
        movaps %xmm1,nb201_fjy(%esp)
        movaps %xmm2,nb201_fjz(%esp)

        ## H2 interactions 
        movaps  %xmm5,%xmm4
        mulps   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        movaps  %xmm5,%xmm7
        movaps  nb201_krsqH2(%esp),%xmm0
        addps   %xmm0,%xmm5     ## xmm6=rinv+ krsq 
        subps   nb201_crf(%esp),%xmm5   ## xmm5=rinv+ krsq-crf 
        mulps   nb201_two(%esp),%xmm0
        subps   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulps   nb201_qqH(%esp),%xmm5   ## vcoul 
        mulps   nb201_qqH(%esp),%xmm7
        mulps  %xmm7,%xmm4              ## total fsH2 in xmm4 

        addps  nb201_vctot(%esp),%xmm5

        movaps nb201_dxH2(%esp),%xmm0
        movaps nb201_dyH2(%esp),%xmm1
        movaps nb201_dzH2(%esp),%xmm2
        movaps %xmm5,nb201_vctot(%esp)
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2

        ## update H2 forces 
        movaps nb201_fixH2(%esp),%xmm3
        movaps nb201_fiyH2(%esp),%xmm4
        movaps nb201_fizH2(%esp),%xmm7
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm7
        movaps %xmm3,nb201_fixH2(%esp)
        movaps %xmm4,nb201_fiyH2(%esp)
        movaps %xmm7,nb201_fizH2(%esp)

        movl nb201_faction(%ebp),%edi
        ## update j forces 
        addps nb201_fjx(%esp),%xmm0
        addps nb201_fjy(%esp),%xmm1
        addps nb201_fjz(%esp),%xmm2

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
        subl $4,nb201_innerk(%esp)
        jl    _nb_kernel201_ia32_sse.nb201_odd_inner
        jmp   _nb_kernel201_ia32_sse.nb201_unroll_loop
_nb_kernel201_ia32_sse.nb201_odd_inner: 
        addl $4,nb201_innerk(%esp)
        jnz   _nb_kernel201_ia32_sse.nb201_odd_loop
        jmp   _nb_kernel201_ia32_sse.nb201_updateouterdata
_nb_kernel201_ia32_sse.nb201_odd_loop: 
        movl  nb201_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb201_innerjjnr(%esp)

        xorps %xmm4,%xmm4
        movss nb201_iqO(%esp),%xmm4
        movl nb201_charge(%ebp),%esi
        movhps nb201_iqH(%esp),%xmm4
        movss (%esi,%eax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb201_qqO(%esp)    ## use oxygen qq for storage 

        movl nb201_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## move j coords to xmm0-xmm2 
        movss (%esi,%eax,4),%xmm0
        movss 4(%esi,%eax,4),%xmm1
        movss 8(%esi,%eax,4),%xmm2
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movss nb201_ixO(%esp),%xmm3
        movss nb201_iyO(%esp),%xmm4
        movss nb201_izO(%esp),%xmm5

        movlps nb201_ixH1(%esp),%xmm6
        movlps nb201_ixH2(%esp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm3
        movlps nb201_iyH1(%esp),%xmm6
        movlps nb201_iyH2(%esp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm4
        movlps nb201_izH1(%esp),%xmm6
        movlps nb201_izH2(%esp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm5

        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        movaps %xmm3,nb201_dxO(%esp)
        movaps %xmm4,nb201_dyO(%esp)
        movaps %xmm5,nb201_dzO(%esp)

        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5

        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        ## rsq in xmm4 

        movaps %xmm4,%xmm0
        mulps nb201_krf(%esp),%xmm0
        movaps %xmm0,nb201_krsqO(%esp)

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb201_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb201_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        ## a little trick to avoid NaNs: 
        ## positions 0,2,and 3 are valid, but not 1. 
        ## If it contains NaN it doesnt help to mult by 0, 
        ## So we shuffle it and copy pos 0 to pos1! 
        shufps $224,%xmm0,%xmm0 ## constant 11100000

        movaps %xmm0,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvsq 

        movaps %xmm0,%xmm1      ## xmm1=rinv 
        movaps nb201_krsqO(%esp),%xmm3
        addps  %xmm3,%xmm0      ## xmm0=rinv+ krsq 
        subps  nb201_crf(%esp),%xmm0   ## xmm0=rinv+ krsq-crf 
        mulps  nb201_two(%esp),%xmm3
        subps  %xmm3,%xmm1      ## xmm1=rinv-2*krsq 
        mulps  nb201_qqO(%esp),%xmm0    ## xmm0=vcoul 
        mulps  nb201_qqO(%esp),%xmm1    ## xmm1=coul part of fs 


        mulps  %xmm1,%xmm4      ## xmm4=total fscal 
        addps  nb201_vctot(%esp),%xmm0
        movaps %xmm0,nb201_vctot(%esp)

        movaps nb201_dxO(%esp),%xmm0
        movaps nb201_dyO(%esp),%xmm1
        movaps nb201_dzO(%esp),%xmm2

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        movss  nb201_fixO(%esp),%xmm3
        movss  nb201_fiyO(%esp),%xmm4
        movss  nb201_fizO(%esp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,nb201_fixO(%esp)
        movss  %xmm4,nb201_fiyO(%esp)
        movss  %xmm5,nb201_fizO(%esp)   ## updated the O force now do the H's 
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        shufps $230,%xmm3,%xmm3 ## constant 11100110     ;# shift right 
        shufps $230,%xmm4,%xmm4 ## constant 11100110
        shufps $230,%xmm5,%xmm5 ## constant 11100110
        addss  nb201_fixH1(%esp),%xmm3
        addss  nb201_fiyH1(%esp),%xmm4
        addss  nb201_fizH1(%esp),%xmm5
        movss  %xmm3,nb201_fixH1(%esp)
        movss  %xmm4,nb201_fiyH1(%esp)
        movss  %xmm5,nb201_fizH1(%esp)          ## updated the H1 force 

        movl nb201_faction(%ebp),%edi
        shufps $231,%xmm3,%xmm3 ## constant 11100111     ;# shift right 
        shufps $231,%xmm4,%xmm4 ## constant 11100111
        shufps $231,%xmm5,%xmm5 ## constant 11100111
        addss  nb201_fixH2(%esp),%xmm3
        addss  nb201_fiyH2(%esp),%xmm4
        addss  nb201_fizH2(%esp),%xmm5
        movss  %xmm3,nb201_fixH2(%esp)
        movss  %xmm4,nb201_fiyH2(%esp)
        movss  %xmm5,nb201_fizH2(%esp)          ## updated the H2 force 

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

        decl nb201_innerk(%esp)
        jz    _nb_kernel201_ia32_sse.nb201_updateouterdata
        jmp   _nb_kernel201_ia32_sse.nb201_odd_loop
_nb_kernel201_ia32_sse.nb201_updateouterdata: 
        movl  nb201_ii3(%esp),%ecx
        movl  nb201_faction(%ebp),%edi
        movl  nb201_fshift(%ebp),%esi
        movl  nb201_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movaps nb201_fixO(%esp),%xmm0
        movaps nb201_fiyO(%esp),%xmm1
        movaps nb201_fizO(%esp),%xmm2

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
        movaps nb201_fixH1(%esp),%xmm0
        movaps nb201_fiyH1(%esp),%xmm1
        movaps nb201_fizH1(%esp),%xmm2

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
        movaps nb201_fixH2(%esp),%xmm0
        movaps nb201_fiyH2(%esp),%xmm1
        movaps nb201_fizH2(%esp),%xmm2

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
        movl nb201_n(%esp),%esi
        ## get group index for i particle 
        movl  nb201_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb201_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb201_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb201_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel201_ia32_sse.nb201_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb201_n(%esp)
        jmp _nb_kernel201_ia32_sse.nb201_outer
_nb_kernel201_ia32_sse.nb201_outerend: 
        ## check if more outer neighborlists remain
        movl  nb201_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel201_ia32_sse.nb201_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel201_ia32_sse.nb201_threadloop
_nb_kernel201_ia32_sse.nb201_end: 
        emms

        movl nb201_nouter(%esp),%eax
        movl nb201_ninner(%esp),%ebx
        movl nb201_outeriter(%ebp),%ecx
        movl nb201_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb201_salign(%esp),%eax
        addl %eax,%esp
        addl $744,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret




.globl nb_kernel201nf_ia32_sse
.globl _nb_kernel201nf_ia32_sse
nb_kernel201nf_ia32_sse:        
_nb_kernel201nf_ia32_sse:       
.set nb201nf_p_nri, 8
.set nb201nf_iinr, 12
.set nb201nf_jindex, 16
.set nb201nf_jjnr, 20
.set nb201nf_shift, 24
.set nb201nf_shiftvec, 28
.set nb201nf_fshift, 32
.set nb201nf_gid, 36
.set nb201nf_pos, 40
.set nb201nf_faction, 44
.set nb201nf_charge, 48
.set nb201nf_p_facel, 52
.set nb201nf_argkrf, 56
.set nb201nf_argcrf, 60
.set nb201nf_Vc, 64
.set nb201nf_type, 68
.set nb201nf_p_ntype, 72
.set nb201nf_vdwparam, 76
.set nb201nf_Vvdw, 80
.set nb201nf_p_tabscale, 84
.set nb201nf_VFtab, 88
.set nb201nf_invsqrta, 92
.set nb201nf_dvda, 96
.set nb201nf_p_gbtabscale, 100
.set nb201nf_GBtab, 104
.set nb201nf_p_nthreads, 108
.set nb201nf_count, 112
.set nb201nf_mtx, 116
.set nb201nf_outeriter, 120
.set nb201nf_inneriter, 124
.set nb201nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb201nf_ixO, 0
.set nb201nf_iyO, 16
.set nb201nf_izO, 32
.set nb201nf_ixH1, 48
.set nb201nf_iyH1, 64
.set nb201nf_izH1, 80
.set nb201nf_ixH2, 96
.set nb201nf_iyH2, 112
.set nb201nf_izH2, 128
.set nb201nf_iqO, 144
.set nb201nf_iqH, 160
.set nb201nf_qqO, 176
.set nb201nf_qqH, 192
.set nb201nf_vctot, 208
.set nb201nf_half, 224
.set nb201nf_three, 240
.set nb201nf_krf, 256
.set nb201nf_crf, 272
.set nb201nf_krsqO, 288
.set nb201nf_krsqH1, 304
.set nb201nf_krsqH2, 320
.set nb201nf_is3, 336
.set nb201nf_ii3, 340
.set nb201nf_innerjjnr, 344
.set nb201nf_innerk, 348
.set nb201nf_n, 352
.set nb201nf_nn1, 356
.set nb201nf_nri, 360
.set nb201nf_nouter, 364
.set nb201nf_ninner, 368
.set nb201nf_salign, 372
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $376,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb201nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb201nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb201nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb201nf_nouter(%esp)
        movl %eax,nb201nf_ninner(%esp)

        movl nb201nf_argkrf(%ebp),%esi
        movl nb201nf_argcrf(%ebp),%edi
        movss (%esi),%xmm5
        movss (%edi),%xmm6
        shufps $0,%xmm5,%xmm5
        shufps $0,%xmm6,%xmm6
        movaps %xmm5,nb201nf_krf(%esp)
        movaps %xmm6,nb201nf_crf(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb201nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb201nf_charge(%ebp),%edx
        movss (%edx,%ebx,4),%xmm3
        movss 4(%edx,%ebx,4),%xmm4
        movl nb201nf_p_facel(%ebp),%esi
        movss (%esi),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb201nf_iqO(%esp)
        movaps %xmm4,nb201nf_iqH(%esp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb201nf_half(%esp)
        movss nb201nf_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb201nf_half(%esp)
        movaps %xmm3,nb201nf_three(%esp)

_nb_kernel201nf_ia32_sse.nb201nf_threadloop: 
        movl  nb201nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel201nf_ia32_sse.nb201nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel201nf_ia32_sse.nb201nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb201nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb201nf_n(%esp)
        movl %ebx,nb201nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel201nf_ia32_sse.nb201nf_outerstart
        jmp _nb_kernel201nf_ia32_sse.nb201nf_end

_nb_kernel201nf_ia32_sse.nb201nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb201nf_nouter(%esp),%ebx
        movl %ebx,nb201nf_nouter(%esp)

_nb_kernel201nf_ia32_sse.nb201nf_outer: 
        movl  nb201nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb201nf_is3(%esp)            ## store is3 

        movl  nb201nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb201nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb201nf_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb201nf_ii3(%esp)

        addss (%eax,%ebx,4),%xmm3
        addss 4(%eax,%ebx,4),%xmm4
        addss 8(%eax,%ebx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb201nf_ixO(%esp)
        movaps %xmm4,nb201nf_iyO(%esp)
        movaps %xmm5,nb201nf_izO(%esp)

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
        movaps %xmm0,nb201nf_ixH1(%esp)
        movaps %xmm1,nb201nf_iyH1(%esp)
        movaps %xmm2,nb201nf_izH1(%esp)
        movaps %xmm3,nb201nf_ixH2(%esp)
        movaps %xmm4,nb201nf_iyH2(%esp)
        movaps %xmm5,nb201nf_izH2(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb201nf_vctot(%esp)

        movl  nb201nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb201nf_pos(%ebp),%esi
        movl  nb201nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb201nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb201nf_ninner(%esp),%ecx
        movl  %ecx,nb201nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb201nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel201nf_ia32_sse.nb201nf_unroll_loop
        jmp   _nb_kernel201nf_ia32_sse.nb201nf_odd_inner
_nb_kernel201nf_ia32_sse.nb201nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb201nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 

        addl $16,nb201nf_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb201nf_charge(%ebp),%esi     ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4           ## and in xmm4 
        mulps  nb201nf_iqO(%esp),%xmm3
        mulps  nb201nf_iqH(%esp),%xmm4

        movaps  %xmm3,nb201nf_qqO(%esp)
        movaps  %xmm4,nb201nf_qqH(%esp)

        movl nb201nf_pos(%ebp),%esi        ## base of pos[] 

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
        movaps nb201nf_ixO(%esp),%xmm4
        movaps nb201nf_iyO(%esp),%xmm5
        movaps nb201nf_izO(%esp),%xmm6

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
        movaps nb201nf_ixH1(%esp),%xmm4
        movaps nb201nf_iyH1(%esp),%xmm5
        movaps nb201nf_izH1(%esp),%xmm6

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
        movaps nb201nf_ixH2(%esp),%xmm3
        movaps nb201nf_iyH2(%esp),%xmm4
        movaps nb201nf_izH2(%esp),%xmm5

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

        mulps  nb201nf_krf(%esp),%xmm0
        mulps  nb201nf_krf(%esp),%xmm1
        mulps  nb201nf_krf(%esp),%xmm2

        movaps %xmm0,nb201nf_krsqH2(%esp)
        movaps %xmm1,nb201nf_krsqH1(%esp)
        movaps %xmm2,nb201nf_krsqO(%esp)

        ## start with rsqO - seed in xmm2       
        rsqrtps %xmm7,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb201nf_three(%esp),%xmm4
        mulps   %xmm7,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb201nf_half(%esp),%xmm4
        movaps  %xmm4,%xmm7     ## rinvO in xmm7 
        ## rsqH1 - seed in xmm2 
        rsqrtps %xmm6,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb201nf_three(%esp),%xmm4
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb201nf_half(%esp),%xmm4
        movaps  %xmm4,%xmm6     ## rinvH1 in xmm6 
        ## rsqH2 - seed in xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb201nf_three(%esp),%xmm4
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb201nf_half(%esp),%xmm4
        movaps  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do O interactions 

        movaps %xmm7,%xmm0
        movaps nb201nf_krsqO(%esp),%xmm1
        addps  %xmm1,%xmm0
        subps  nb201nf_crf(%esp),%xmm0   ## xmm0=rinv+ krsq-crf 
        mulps  nb201nf_qqO(%esp),%xmm0

        addps  nb201nf_vctot(%esp),%xmm0
        movaps %xmm0,nb201nf_vctot(%esp)

        ## H1 interactions 
        movaps  %xmm6,%xmm7
        movaps  nb201nf_krsqH1(%esp),%xmm0
        addps   %xmm0,%xmm6     ## xmm6=rinv+ krsq 
        subps   nb201nf_crf(%esp),%xmm6   ## xmm6=rinv+ krsq-crf 
        mulps   nb201nf_qqH(%esp),%xmm6   ## vcoul 
        addps   nb201nf_vctot(%esp),%xmm6

        ## H2 interactions 
        movaps  %xmm5,%xmm7
        movaps  nb201nf_krsqH2(%esp),%xmm0
        addps   %xmm0,%xmm5     ## xmm6=rinv+ krsq 
        subps   nb201nf_crf(%esp),%xmm5   ## xmm5=rinv+ krsq-crf 
        mulps   nb201nf_qqH(%esp),%xmm5   ## vcoul 
        addps  %xmm5,%xmm6
        movaps %xmm6,nb201nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $4,nb201nf_innerk(%esp)
        jl    _nb_kernel201nf_ia32_sse.nb201nf_odd_inner
        jmp   _nb_kernel201nf_ia32_sse.nb201nf_unroll_loop
_nb_kernel201nf_ia32_sse.nb201nf_odd_inner: 
        addl $4,nb201nf_innerk(%esp)
        jnz   _nb_kernel201nf_ia32_sse.nb201nf_odd_loop
        jmp   _nb_kernel201nf_ia32_sse.nb201nf_updateouterdata
_nb_kernel201nf_ia32_sse.nb201nf_odd_loop: 
        movl  nb201nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb201nf_innerjjnr(%esp)

        xorps %xmm4,%xmm4
        movss nb201nf_iqO(%esp),%xmm4
        movl nb201nf_charge(%ebp),%esi
        movhps nb201nf_iqH(%esp),%xmm4
        movss (%esi,%eax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb201nf_qqO(%esp)          ## use oxygen qq for storage 

        movl nb201nf_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## move j coords to xmm0-xmm2 
        movss (%esi,%eax,4),%xmm0
        movss 4(%esi,%eax,4),%xmm1
        movss 8(%esi,%eax,4),%xmm2
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movss nb201nf_ixO(%esp),%xmm3
        movss nb201nf_iyO(%esp),%xmm4
        movss nb201nf_izO(%esp),%xmm5

        movlps nb201nf_ixH1(%esp),%xmm6
        movlps nb201nf_ixH2(%esp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm3
        movlps nb201nf_iyH1(%esp),%xmm6
        movlps nb201nf_iyH2(%esp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm4
        movlps nb201nf_izH1(%esp),%xmm6
        movlps nb201nf_izH2(%esp),%xmm7
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
        mulps nb201nf_krf(%esp),%xmm0
        movaps %xmm0,nb201nf_krsqO(%esp)

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb201nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb201nf_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 

        ## a little trick to avoid NaNs: 
        ## positions 0,2,and 3 are valid, but not 1. 
        ## If it contains NaN it doesnt help to mult by 0, 
        ## So we shuffle it and copy pos 0 to pos1! 
        shufps $224,%xmm0,%xmm0 ## constant 11100000     

        movaps nb201nf_krsqO(%esp),%xmm3
        addps  %xmm3,%xmm0      ## xmm0=rinv+ krsq 
        subps  nb201nf_crf(%esp),%xmm0   ## xmm0=rinv+ krsq-crf 
        mulps  nb201nf_qqO(%esp),%xmm0          ## xmm0=vcoul 
        addps  nb201nf_vctot(%esp),%xmm0
        movaps %xmm0,nb201nf_vctot(%esp)

        decl nb201nf_innerk(%esp)
        jz    _nb_kernel201nf_ia32_sse.nb201nf_updateouterdata
        jmp   _nb_kernel201nf_ia32_sse.nb201nf_odd_loop
_nb_kernel201nf_ia32_sse.nb201nf_updateouterdata: 
        ## get n from stack
        movl nb201nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb201nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb201nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb201nf_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb201nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel201nf_ia32_sse.nb201nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb201nf_n(%esp)
        jmp _nb_kernel201nf_ia32_sse.nb201nf_outer
_nb_kernel201nf_ia32_sse.nb201nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb201nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel201nf_ia32_sse.nb201nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel201nf_ia32_sse.nb201nf_threadloop
_nb_kernel201nf_ia32_sse.nb201nf_end: 
        emms

        movl nb201nf_nouter(%esp),%eax
        movl nb201nf_ninner(%esp),%ebx
        movl nb201nf_outeriter(%ebp),%ecx
        movl nb201nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb201nf_salign(%esp),%eax
        addl %eax,%esp
        addl $376,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret

