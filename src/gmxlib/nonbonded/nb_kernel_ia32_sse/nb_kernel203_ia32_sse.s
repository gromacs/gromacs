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


.globl nb_kernel203_ia32_sse
.globl _nb_kernel203_ia32_sse
nb_kernel203_ia32_sse:  
_nb_kernel203_ia32_sse: 
.set nb203_p_nri, 8
.set nb203_iinr, 12
.set nb203_jindex, 16
.set nb203_jjnr, 20
.set nb203_shift, 24
.set nb203_shiftvec, 28
.set nb203_fshift, 32
.set nb203_gid, 36
.set nb203_pos, 40
.set nb203_faction, 44
.set nb203_charge, 48
.set nb203_p_facel, 52
.set nb203_argkrf, 56
.set nb203_argcrf, 60
.set nb203_Vc, 64
.set nb203_type, 68
.set nb203_p_ntype, 72
.set nb203_vdwparam, 76
.set nb203_Vvdw, 80
.set nb203_p_tabscale, 84
.set nb203_VFtab, 88
.set nb203_invsqrta, 92
.set nb203_dvda, 96
.set nb203_p_gbtabscale, 100
.set nb203_GBtab, 104
.set nb203_p_nthreads, 108
.set nb203_count, 112
.set nb203_mtx, 116
.set nb203_outeriter, 120
.set nb203_inneriter, 124
.set nb203_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb203_ixH1, 0
.set nb203_iyH1, 16
.set nb203_izH1, 32
.set nb203_ixH2, 48
.set nb203_iyH2, 64
.set nb203_izH2, 80
.set nb203_ixM, 96
.set nb203_iyM, 112
.set nb203_izM, 128
.set nb203_iqH, 144
.set nb203_iqM, 160
.set nb203_dxH1, 176
.set nb203_dyH1, 192
.set nb203_dzH1, 208
.set nb203_dxH2, 224
.set nb203_dyH2, 240
.set nb203_dzH2, 256
.set nb203_dxM, 272
.set nb203_dyM, 288
.set nb203_dzM, 304
.set nb203_qqH, 320
.set nb203_qqM, 336
.set nb203_vctot, 352
.set nb203_fixH1, 384
.set nb203_fiyH1, 400
.set nb203_fizH1, 416
.set nb203_fixH2, 432
.set nb203_fiyH2, 448
.set nb203_fizH2, 464
.set nb203_fixM, 480
.set nb203_fiyM, 496
.set nb203_fizM, 512
.set nb203_fjx, 528
.set nb203_fjy, 544
.set nb203_fjz, 560
.set nb203_half, 576
.set nb203_three, 592
.set nb203_two, 608
.set nb203_krf, 624
.set nb203_crf, 640
.set nb203_krsqH1, 656
.set nb203_krsqH2, 672
.set nb203_krsqM, 688
.set nb203_is3, 704
.set nb203_ii3, 708
.set nb203_innerjjnr, 712
.set nb203_innerk, 716
.set nb203_n, 720
.set nb203_nn1, 724
.set nb203_nri, 728
.set nb203_nouter, 732
.set nb203_ninner, 736
.set nb203_salign, 740
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
        movl %eax,nb203_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb203_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb203_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb203_nouter(%esp)
        movl %eax,nb203_ninner(%esp)


        movl nb203_argkrf(%ebp),%esi
        movl nb203_argcrf(%ebp),%edi
        movss (%esi),%xmm5
        movss (%edi),%xmm6
        shufps $0,%xmm5,%xmm5
        shufps $0,%xmm6,%xmm6
        movaps %xmm5,nb203_krf(%esp)
        movaps %xmm6,nb203_crf(%esp)
        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb203_half(%esp)
        movss nb203_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb203_half(%esp)
        movaps %xmm2,nb203_two(%esp)
        movaps %xmm3,nb203_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb203_iinr(%ebp),%ecx             ## ecx = pointer into iinr[]    
        movl  (%ecx),%ebx               ## ebx =ii 

        movl  nb203_charge(%ebp),%edx
        movss 4(%edx,%ebx,4),%xmm3
        movss 12(%edx,%ebx,4),%xmm4
        movl nb203_p_facel(%ebp),%esi
        movss (%esi),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb203_iqH(%esp)
        movaps %xmm4,nb203_iqM(%esp)

_nb_kernel203_ia32_sse.nb203_threadloop: 
        movl  nb203_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel203_ia32_sse.nb203_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel203_ia32_sse.nb203_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb203_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb203_n(%esp)
        movl %ebx,nb203_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel203_ia32_sse.nb203_outerstart
        jmp _nb_kernel203_ia32_sse.nb203_end

_nb_kernel203_ia32_sse.nb203_outerstart: 
        ## ebx contains number of outer iterations
        addl nb203_nouter(%esp),%ebx
        movl %ebx,nb203_nouter(%esp)

_nb_kernel203_ia32_sse.nb203_outer: 
        movl  nb203_shift(%ebp),%eax            ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx        ## ebx=3*is 
        movl  %ebx,nb203_is3(%esp)      ## store is3 

        movl  nb203_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb203_iinr(%ebp),%ecx             ## ecx = pointer into iinr[]    
        movl  (%ecx,%esi,4),%ebx                ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb203_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb203_ii3(%esp)

        addss 12(%eax,%ebx,4),%xmm3
        addss 16(%eax,%ebx,4),%xmm4
        addss 20(%eax,%ebx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb203_ixH1(%esp)
        movaps %xmm4,nb203_iyH1(%esp)
        movaps %xmm5,nb203_izH1(%esp)

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
        movaps %xmm0,nb203_ixH2(%esp)
        movaps %xmm1,nb203_iyH2(%esp)
        movaps %xmm2,nb203_izH2(%esp)
        movaps %xmm3,nb203_ixM(%esp)
        movaps %xmm4,nb203_iyM(%esp)
        movaps %xmm5,nb203_izM(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb203_vctot(%esp)
        movaps %xmm4,nb203_fixH1(%esp)
        movaps %xmm4,nb203_fiyH1(%esp)
        movaps %xmm4,nb203_fizH1(%esp)
        movaps %xmm4,nb203_fixH2(%esp)
        movaps %xmm4,nb203_fiyH2(%esp)
        movaps %xmm4,nb203_fizH2(%esp)
        movaps %xmm4,nb203_fixM(%esp)
        movaps %xmm4,nb203_fiyM(%esp)
        movaps %xmm4,nb203_fizM(%esp)

        movl  nb203_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx                ## jindex[n] 
        movl  4(%eax,%esi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movl  nb203_pos(%ebp),%esi
        movl  nb203_faction(%ebp),%edi
        movl  nb203_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb203_innerjjnr(%esp)        ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb203_ninner(%esp),%ecx
        movl  %ecx,nb203_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb203_innerk(%esp)   ## number of innerloop atoms 
        jge   _nb_kernel203_ia32_sse.nb203_unroll_loop
        jmp   _nb_kernel203_ia32_sse.nb203_odd_inner
_nb_kernel203_ia32_sse.nb203_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb203_innerjjnr(%esp),%edx        ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx             ## eax-edx=jnr1-4 

        addl $16,nb203_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb203_charge(%ebp),%esi    ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4              ## and in xmm4 
        mulps  nb203_iqH(%esp),%xmm3
        mulps  nb203_iqM(%esp),%xmm4

        movaps  %xmm3,nb203_qqH(%esp)
        movaps  %xmm4,nb203_qqM(%esp)

        movl nb203_pos(%ebp),%esi       ## base of pos[] 

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

        shufps $136,%xmm6,%xmm2 ## constant 10001000

        shufps $136,%xmm5,%xmm0 ## constant 10001000
        shufps $221,%xmm5,%xmm1 ## constant 11011101            

        ## move ixH1-izH1 to xmm4-xmm6 
        movaps nb203_ixH1(%esp),%xmm4
        movaps nb203_iyH1(%esp),%xmm5
        movaps nb203_izH1(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb203_dxH1(%esp)
        movaps %xmm5,nb203_dyH1(%esp)
        movaps %xmm6,nb203_dzH1(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        movaps %xmm4,%xmm7
        ## rsqH1 in xmm7 

        ## move ixH2-izH2 to xmm4-xmm6 
        movaps nb203_ixH2(%esp),%xmm4
        movaps nb203_iyH2(%esp),%xmm5
        movaps nb203_izH2(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb203_dxH2(%esp)
        movaps %xmm5,nb203_dyH2(%esp)
        movaps %xmm6,nb203_dzH2(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm6
        addps %xmm4,%xmm6
        ## rsqH2 in xmm6 

        ## move ixM-izM to xmm3-xmm5  
        movaps nb203_ixM(%esp),%xmm3
        movaps nb203_iyM(%esp),%xmm4
        movaps nb203_izM(%esp),%xmm5

        ## calc dr 
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        ## store dr 
        movaps %xmm3,nb203_dxM(%esp)
        movaps %xmm4,nb203_dyM(%esp)
        movaps %xmm5,nb203_dzM(%esp)
        ## square it 
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm4,%xmm5
        addps %xmm3,%xmm5
        ## rsqM in xmm5, rsqH2 in xmm6, rsqH1 in xmm7 

        movaps %xmm5,%xmm0
        movaps %xmm6,%xmm1
        movaps %xmm7,%xmm2

        mulps  nb203_krf(%esp),%xmm0
        mulps  nb203_krf(%esp),%xmm1
        mulps  nb203_krf(%esp),%xmm2

        movaps %xmm0,nb203_krsqM(%esp)
        movaps %xmm1,nb203_krsqH2(%esp)
        movaps %xmm2,nb203_krsqH1(%esp)

        ## start with rsqH1 - seed in xmm2      
        rsqrtps %xmm7,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb203_three(%esp),%xmm4
        mulps   %xmm7,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb203_half(%esp),%xmm4
        movaps  %xmm4,%xmm7     ## rinvH1 in xmm7 
        ## rsqH2 - seed in xmm2 
        rsqrtps %xmm6,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb203_three(%esp),%xmm4
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb203_half(%esp),%xmm4
        movaps  %xmm4,%xmm6     ## rinvH2 in xmm6 
        ## rsqM - seed in xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb203_three(%esp),%xmm4
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb203_half(%esp),%xmm4
        movaps  %xmm4,%xmm5     ## rinvM in xmm5 

        ## do H1 interactions 
        movaps  %xmm7,%xmm4
        mulps   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq 

        movaps %xmm7,%xmm0
        movaps nb203_krsqH1(%esp),%xmm1
        addps  %xmm1,%xmm0
        subps  nb203_crf(%esp),%xmm0   ## xmm0=rinv+ krsq-crf 
        mulps  nb203_two(%esp),%xmm1
        subps  %xmm1,%xmm7
        mulps  nb203_qqH(%esp),%xmm0
        mulps  nb203_qqH(%esp),%xmm7

        mulps  %xmm7,%xmm4      ## total fs H1 in xmm4 

        addps  nb203_vctot(%esp),%xmm0
        movaps %xmm0,nb203_vctot(%esp)

        movaps nb203_dxH1(%esp),%xmm0
        movaps nb203_dyH1(%esp),%xmm1
        movaps nb203_dzH1(%esp),%xmm2
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2

        ## update H1 forces 
        movaps nb203_fixH1(%esp),%xmm3
        movaps nb203_fiyH1(%esp),%xmm4
        movaps nb203_fizH1(%esp),%xmm7
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm7
        movaps %xmm3,nb203_fixH1(%esp)
        movaps %xmm4,nb203_fiyH1(%esp)
        movaps %xmm7,nb203_fizH1(%esp)
        ## update j forces with water O 
        movaps %xmm0,nb203_fjx(%esp)
        movaps %xmm1,nb203_fjy(%esp)
        movaps %xmm2,nb203_fjz(%esp)

        ## H2 interactions 
        movaps  %xmm6,%xmm4
        mulps   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        movaps  %xmm6,%xmm7
        movaps  nb203_krsqH2(%esp),%xmm0
        addps   %xmm0,%xmm6     ## xmm6=rinv+ krsq 
        subps   nb203_crf(%esp),%xmm6   ## xmm6=rinv+ krsq-crf 
        mulps   nb203_two(%esp),%xmm0
        subps   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulps   nb203_qqH(%esp),%xmm6   ## vcoul 
        mulps   nb203_qqH(%esp),%xmm7
        mulps  %xmm7,%xmm4              ## total fsH2 in xmm4 

        addps  nb203_vctot(%esp),%xmm6

        movaps nb203_dxH2(%esp),%xmm0
        movaps nb203_dyH2(%esp),%xmm1
        movaps nb203_dzH2(%esp),%xmm2
        movaps %xmm6,nb203_vctot(%esp)
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2

        ## update H2 forces 
        movaps nb203_fixH2(%esp),%xmm3
        movaps nb203_fiyH2(%esp),%xmm4
        movaps nb203_fizH2(%esp),%xmm7
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm7
        movaps %xmm3,nb203_fixH2(%esp)
        movaps %xmm4,nb203_fiyH2(%esp)
        movaps %xmm7,nb203_fizH2(%esp)
        ## update j forces with water H2
        addps  nb203_fjx(%esp),%xmm0
        addps  nb203_fjy(%esp),%xmm1
        addps  nb203_fjz(%esp),%xmm2
        movaps %xmm0,nb203_fjx(%esp)
        movaps %xmm1,nb203_fjy(%esp)
        movaps %xmm2,nb203_fjz(%esp)

        ## M interactions 
        movaps  %xmm5,%xmm4
        mulps   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        movaps  %xmm5,%xmm7
        movaps  nb203_krsqM(%esp),%xmm0
        addps   %xmm0,%xmm5     ## xmm6=rinv+ krsq 
        subps   nb203_crf(%esp),%xmm5   ## xmm5=rinv+ krsq-crf 
        mulps   nb203_two(%esp),%xmm0
        subps   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulps   nb203_qqM(%esp),%xmm5   ## vcoul 
        mulps   nb203_qqM(%esp),%xmm7
        mulps  %xmm7,%xmm4              ## total fsM in xmm4 

        addps  nb203_vctot(%esp),%xmm5

        movaps nb203_dxM(%esp),%xmm0
        movaps nb203_dyM(%esp),%xmm1
        movaps nb203_dzM(%esp),%xmm2
        movaps %xmm5,nb203_vctot(%esp)
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2

        ## update M forces 
        movaps nb203_fixM(%esp),%xmm3
        movaps nb203_fiyM(%esp),%xmm4
        movaps nb203_fizM(%esp),%xmm7
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm7
        movaps %xmm3,nb203_fixM(%esp)
        movaps %xmm4,nb203_fiyM(%esp)
        movaps %xmm7,nb203_fizM(%esp)

        movl nb203_faction(%ebp),%edi

        ## update j forces 
        addps nb203_fjx(%esp),%xmm0
        addps nb203_fjy(%esp),%xmm1
        addps nb203_fjz(%esp),%xmm2

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
        subl $4,nb203_innerk(%esp)
        jl    _nb_kernel203_ia32_sse.nb203_odd_inner
        jmp   _nb_kernel203_ia32_sse.nb203_unroll_loop
_nb_kernel203_ia32_sse.nb203_odd_inner: 
        addl $4,nb203_innerk(%esp)
        jnz   _nb_kernel203_ia32_sse.nb203_odd_loop
        jmp   _nb_kernel203_ia32_sse.nb203_updateouterdata
_nb_kernel203_ia32_sse.nb203_odd_loop: 
        movl  nb203_innerjjnr(%esp),%edx        ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb203_innerjjnr(%esp)

        xorps %xmm4,%xmm4
        movss nb203_iqM(%esp),%xmm4
        movl nb203_charge(%ebp),%esi
        movhps nb203_iqH(%esp),%xmm4
        movss (%esi,%eax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb203_qqM(%esp)    ## use dummy qq for storage 

        movl nb203_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## move j coords to xmm0-xmm2 
        movss (%esi,%eax,4),%xmm0
        movss 4(%esi,%eax,4),%xmm1
        movss 8(%esi,%eax,4),%xmm2
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movss nb203_ixM(%esp),%xmm3
        movss nb203_iyM(%esp),%xmm4
        movss nb203_izM(%esp),%xmm5

        movlps nb203_ixH1(%esp),%xmm6
        movlps nb203_ixH2(%esp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm3
        movlps nb203_iyH1(%esp),%xmm6
        movlps nb203_iyH2(%esp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm4
        movlps nb203_izH1(%esp),%xmm6
        movlps nb203_izH2(%esp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm5

        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        ## use dummy dx for storage
        movaps %xmm3,nb203_dxM(%esp)
        movaps %xmm4,nb203_dyM(%esp)
        movaps %xmm5,nb203_dzM(%esp)

        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5

        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        ## rsq in xmm4 

        movaps %xmm4,%xmm0
        mulps nb203_krf(%esp),%xmm0
        movaps %xmm0,nb203_krsqM(%esp)

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb203_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb203_half(%esp),%xmm0
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
        movaps nb203_krsqM(%esp),%xmm3
        addps  %xmm3,%xmm0      ## xmm0=rinv+ krsq 
        subps  nb203_crf(%esp),%xmm0   ## xmm0=rinv+ krsq-crf 
        mulps  nb203_two(%esp),%xmm3
        subps  %xmm3,%xmm1      ## xmm1=rinv-2*krsq 
        mulps  nb203_qqM(%esp),%xmm0    ## xmm0=vcoul 
        mulps  nb203_qqM(%esp),%xmm1    ## xmm1=coul part of fs 


        mulps  %xmm1,%xmm4      ## xmm4=total fscal 
        addps  nb203_vctot(%esp),%xmm0
        movaps %xmm0,nb203_vctot(%esp)

        movaps nb203_dxM(%esp),%xmm0
        movaps nb203_dyM(%esp),%xmm1
        movaps nb203_dzM(%esp),%xmm2

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        movss  nb203_fixM(%esp),%xmm3
        movss  nb203_fiyM(%esp),%xmm4
        movss  nb203_fizM(%esp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,nb203_fixM(%esp)
        movss  %xmm4,nb203_fiyM(%esp)
        movss  %xmm5,nb203_fizM(%esp)   ## updated the O force now do the H's 
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        shufps $230,%xmm3,%xmm3 ## constant 11100110     ;# shift right 
        shufps $230,%xmm4,%xmm4 ## constant 11100110
        shufps $230,%xmm5,%xmm5 ## constant 11100110
        addss  nb203_fixH1(%esp),%xmm3
        addss  nb203_fiyH1(%esp),%xmm4
        addss  nb203_fizH1(%esp),%xmm5
        movss  %xmm3,nb203_fixH1(%esp)
        movss  %xmm4,nb203_fiyH1(%esp)
        movss  %xmm5,nb203_fizH1(%esp)          ## updated the H1 force 

        movl nb203_faction(%ebp),%edi
        shufps $231,%xmm3,%xmm3 ## constant 11100111     ;# shift right 
        shufps $231,%xmm4,%xmm4 ## constant 11100111
        shufps $231,%xmm5,%xmm5 ## constant 11100111
        addss  nb203_fixH2(%esp),%xmm3
        addss  nb203_fiyH2(%esp),%xmm4
        addss  nb203_fizH2(%esp),%xmm5
        movss  %xmm3,nb203_fixH2(%esp)
        movss  %xmm4,nb203_fiyH2(%esp)
        movss  %xmm5,nb203_fizH2(%esp)          ## updated the H2 force 

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

        decl nb203_innerk(%esp)
        jz    _nb_kernel203_ia32_sse.nb203_updateouterdata
        jmp   _nb_kernel203_ia32_sse.nb203_odd_loop
_nb_kernel203_ia32_sse.nb203_updateouterdata: 
        movl  nb203_ii3(%esp),%ecx
        movl  nb203_faction(%ebp),%edi
        movl  nb203_fshift(%ebp),%esi
        movl  nb203_is3(%esp),%edx

        ## accumulate  H1 i forces in xmm0, xmm1, xmm2 
        movaps nb203_fixH1(%esp),%xmm0
        movaps nb203_fiyH1(%esp),%xmm1
        movaps nb203_fizH1(%esp),%xmm2

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

        ## accumulate H2 i forces in xmm0, xmm1, xmm2 
        movaps nb203_fixH2(%esp),%xmm0
        movaps nb203_fiyH2(%esp),%xmm1
        movaps nb203_fizH2(%esp),%xmm2

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

        ## accumulate m i forces in xmm0, xmm1, xmm2 
        movaps nb203_fixM(%esp),%xmm0
        movaps nb203_fiyM(%esp),%xmm1
        movaps nb203_fizM(%esp),%xmm2

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
        movl nb203_n(%esp),%esi
        ## get group index for i particle 
        movl  nb203_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb203_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb203_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb203_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel203_ia32_sse.nb203_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb203_n(%esp)
        jmp _nb_kernel203_ia32_sse.nb203_outer
_nb_kernel203_ia32_sse.nb203_outerend: 
        ## check if more outer neighborlists remain
        movl  nb203_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel203_ia32_sse.nb203_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel203_ia32_sse.nb203_threadloop
_nb_kernel203_ia32_sse.nb203_end: 
        emms

        movl nb203_nouter(%esp),%eax
        movl nb203_ninner(%esp),%ebx
        movl nb203_outeriter(%ebp),%ecx
        movl nb203_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb203_salign(%esp),%eax
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




.globl nb_kernel203nf_ia32_sse
.globl _nb_kernel203nf_ia32_sse
nb_kernel203nf_ia32_sse:        
_nb_kernel203nf_ia32_sse:       
.set nb203nf_p_nri, 8
.set nb203nf_iinr, 12
.set nb203nf_jindex, 16
.set nb203nf_jjnr, 20
.set nb203nf_shift, 24
.set nb203nf_shiftvec, 28
.set nb203nf_fshift, 32
.set nb203nf_gid, 36
.set nb203nf_pos, 40
.set nb203nf_faction, 44
.set nb203nf_charge, 48
.set nb203nf_p_facel, 52
.set nb203nf_argkrf, 56
.set nb203nf_argcrf, 60
.set nb203nf_Vc, 64
.set nb203nf_type, 68
.set nb203nf_p_ntype, 72
.set nb203nf_vdwparam, 76
.set nb203nf_Vvdw, 80
.set nb203nf_p_tabscale, 84
.set nb203nf_VFtab, 88
.set nb203nf_invsqrta, 92
.set nb203nf_dvda, 96
.set nb203nf_p_gbtabscale, 100
.set nb203nf_GBtab, 104
.set nb203nf_p_nthreads, 108
.set nb203nf_count, 112
.set nb203nf_mtx, 116
.set nb203nf_outeriter, 120
.set nb203nf_inneriter, 124
.set nb203nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb203nf_ixH1, 0
.set nb203nf_iyH1, 16
.set nb203nf_izH1, 32
.set nb203nf_ixH2, 48
.set nb203nf_iyH2, 64
.set nb203nf_izH2, 80
.set nb203nf_ixM, 96
.set nb203nf_iyM, 112
.set nb203nf_izM, 128
.set nb203nf_iqH, 144
.set nb203nf_iqM, 160
.set nb203nf_qqH, 176
.set nb203nf_qqM, 192
.set nb203nf_vctot, 208
.set nb203nf_half, 224
.set nb203nf_three, 240
.set nb203nf_krf, 256
.set nb203nf_crf, 272
.set nb203nf_krsqH1, 288
.set nb203nf_krsqH2, 304
.set nb203nf_krsqM, 320
.set nb203nf_is3, 336
.set nb203nf_ii3, 340
.set nb203nf_innerjjnr, 344
.set nb203nf_innerk, 348
.set nb203nf_n, 352
.set nb203nf_nn1, 356
.set nb203nf_nri, 360
.set nb203nf_nouter, 364
.set nb203nf_ninner, 368
.set nb203nf_salign, 372
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
        movl %eax,nb203nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb203nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb203nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb203nf_nouter(%esp)
        movl %eax,nb203nf_ninner(%esp)


        movl nb203nf_argkrf(%ebp),%esi
        movl nb203nf_argcrf(%ebp),%edi
        movss (%esi),%xmm5
        movss (%edi),%xmm6
        shufps $0,%xmm5,%xmm5
        shufps $0,%xmm6,%xmm6
        movaps %xmm5,nb203nf_krf(%esp)
        movaps %xmm6,nb203nf_crf(%esp)
        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb203nf_half(%esp)
        movss nb203nf_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb203nf_half(%esp)
        movaps %xmm3,nb203nf_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb203nf_iinr(%ebp),%ecx           ## ecx = pointer into iinr[]    
        movl  (%ecx),%ebx               ## ebx =ii 

        movl  nb203nf_charge(%ebp),%edx
        movss 4(%edx,%ebx,4),%xmm3
        movss 12(%edx,%ebx,4),%xmm4
        movl nb203nf_p_facel(%ebp),%esi
        movss (%esi),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb203nf_iqH(%esp)
        movaps %xmm4,nb203nf_iqM(%esp)

_nb_kernel203nf_ia32_sse.nb203nf_threadloop: 
        movl  nb203nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel203nf_ia32_sse.nb203nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel203nf_ia32_sse.nb203nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb203nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb203nf_n(%esp)
        movl %ebx,nb203nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel203nf_ia32_sse.nb203nf_outerstart
        jmp _nb_kernel203nf_ia32_sse.nb203nf_end

_nb_kernel203nf_ia32_sse.nb203nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb203nf_nouter(%esp),%ebx
        movl %ebx,nb203nf_nouter(%esp)

_nb_kernel203nf_ia32_sse.nb203nf_outer: 
        movl  nb203nf_shift(%ebp),%eax          ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx        ## ebx=3*is 
        movl  %ebx,nb203nf_is3(%esp)            ## store is3 

        movl  nb203nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb203nf_iinr(%ebp),%ecx           ## ecx = pointer into iinr[]    
        movl  (%ecx,%esi,4),%ebx                ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb203nf_pos(%ebp),%eax    ## eax = base of pos[]  
        movl  %ebx,nb203nf_ii3(%esp)

        addss 12(%eax,%ebx,4),%xmm3
        addss 16(%eax,%ebx,4),%xmm4
        addss 20(%eax,%ebx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb203nf_ixH1(%esp)
        movaps %xmm4,nb203nf_iyH1(%esp)
        movaps %xmm5,nb203nf_izH1(%esp)

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
        movaps %xmm0,nb203nf_ixH2(%esp)
        movaps %xmm1,nb203nf_iyH2(%esp)
        movaps %xmm2,nb203nf_izH2(%esp)
        movaps %xmm3,nb203nf_ixM(%esp)
        movaps %xmm4,nb203nf_iyM(%esp)
        movaps %xmm5,nb203nf_izM(%esp)

        ## clear vctot 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb203nf_vctot(%esp)

        movl  nb203nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx                ## jindex[n] 
        movl  4(%eax,%esi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movl  nb203nf_pos(%ebp),%esi
        movl  nb203nf_faction(%ebp),%edi
        movl  nb203nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb203nf_innerjjnr(%esp)      ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb203nf_ninner(%esp),%ecx
        movl  %ecx,nb203nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb203nf_innerk(%esp)         ## number of innerloop atoms 
        jge   _nb_kernel203nf_ia32_sse.nb203nf_unroll_loop
        jmp   _nb_kernel203nf_ia32_sse.nb203nf_odd_inner
_nb_kernel203nf_ia32_sse.nb203nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb203nf_innerjjnr(%esp),%edx      ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx             ## eax-edx=jnr1-4 

        addl $16,nb203nf_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb203nf_charge(%ebp),%esi  ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4              ## and in xmm4 
        mulps  nb203nf_iqH(%esp),%xmm3
        mulps  nb203nf_iqM(%esp),%xmm4

        movaps  %xmm3,nb203nf_qqH(%esp)
        movaps  %xmm4,nb203nf_qqM(%esp)

        movl nb203nf_pos(%ebp),%esi     ## base of pos[] 

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

        shufps $136,%xmm6,%xmm2 ## constant 10001000

        shufps $136,%xmm5,%xmm0 ## constant 10001000
        shufps $221,%xmm5,%xmm1 ## constant 11011101            

        ## move ixH1-izH1 to xmm4-xmm6 
        movaps nb203nf_ixH1(%esp),%xmm4
        movaps nb203nf_iyH1(%esp),%xmm5
        movaps nb203nf_izH1(%esp),%xmm6

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
        movaps nb203nf_ixH2(%esp),%xmm4
        movaps nb203nf_iyH2(%esp),%xmm5
        movaps nb203nf_izH2(%esp),%xmm6

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
        movaps nb203nf_ixM(%esp),%xmm3
        movaps nb203nf_iyM(%esp),%xmm4
        movaps nb203nf_izM(%esp),%xmm5

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

        movaps %xmm5,%xmm0
        movaps %xmm6,%xmm1
        movaps %xmm7,%xmm2

        mulps  nb203nf_krf(%esp),%xmm0
        mulps  nb203nf_krf(%esp),%xmm1
        mulps  nb203nf_krf(%esp),%xmm2

        movaps %xmm0,nb203nf_krsqM(%esp)
        movaps %xmm1,nb203nf_krsqH2(%esp)
        movaps %xmm2,nb203nf_krsqH1(%esp)

        ## start with rsqH1 - seed in xmm2      
        rsqrtps %xmm7,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb203nf_three(%esp),%xmm4
        mulps   %xmm7,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb203nf_half(%esp),%xmm4
        movaps  %xmm4,%xmm7     ## rinvH1 in xmm7 
        ## rsqH2 - seed in xmm2 
        rsqrtps %xmm6,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb203nf_three(%esp),%xmm4
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb203nf_half(%esp),%xmm4
        movaps  %xmm4,%xmm6     ## rinvH2 in xmm6 
        ## rsqM - seed in xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb203nf_three(%esp),%xmm4
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb203nf_half(%esp),%xmm4
        movaps  %xmm4,%xmm5     ## rinvM in xmm5 

        ## do H1 interactions - xmm7=rinv
        addps nb203nf_krsqH1(%esp),%xmm7
        subps nb203nf_crf(%esp),%xmm7   ## xmm7=rinv+ krsq-crf 
        mulps nb203nf_qqH(%esp),%xmm7
        addps nb203nf_vctot(%esp),%xmm7

        ## H2 interactions - xmm6=rinv
        addps nb203nf_krsqH2(%esp),%xmm6
        subps nb203nf_crf(%esp),%xmm6   ## xmm6=rinv+ krsq-crf 
        mulps nb203nf_qqH(%esp),%xmm6
        addps %xmm7,%xmm6

        ## M interactions - xmm5=rinv
        addps nb203nf_krsqM(%esp),%xmm5
        subps nb203nf_crf(%esp),%xmm5   ## xmm5=rinv+ krsq-crf 
        mulps nb203nf_qqM(%esp),%xmm5
        addps %xmm6,%xmm5
        movaps %xmm5,nb203nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $4,nb203nf_innerk(%esp)
        jl    _nb_kernel203nf_ia32_sse.nb203nf_odd_inner
        jmp   _nb_kernel203nf_ia32_sse.nb203nf_unroll_loop
_nb_kernel203nf_ia32_sse.nb203nf_odd_inner: 
        addl $4,nb203nf_innerk(%esp)
        jnz   _nb_kernel203nf_ia32_sse.nb203nf_odd_loop
        jmp   _nb_kernel203nf_ia32_sse.nb203nf_updateouterdata
_nb_kernel203nf_ia32_sse.nb203nf_odd_loop: 
        movl  nb203nf_innerjjnr(%esp),%edx      ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb203nf_innerjjnr(%esp)

        xorps %xmm4,%xmm4
        movss nb203nf_iqM(%esp),%xmm4
        movl nb203nf_charge(%ebp),%esi
        movhps nb203nf_iqH(%esp),%xmm4
        movss (%esi,%eax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb203nf_qqM(%esp)          ## use dummy qq for storage 

        movl nb203nf_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## move j coords to xmm0-xmm2 
        movss (%esi,%eax,4),%xmm0
        movss 4(%esi,%eax,4),%xmm1
        movss 8(%esi,%eax,4),%xmm2
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movss nb203nf_ixM(%esp),%xmm3
        movss nb203nf_iyM(%esp),%xmm4
        movss nb203nf_izM(%esp),%xmm5

        movlps nb203nf_ixH1(%esp),%xmm6
        movlps nb203nf_ixH2(%esp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm3
        movlps nb203nf_iyH1(%esp),%xmm6
        movlps nb203nf_iyH2(%esp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm4
        movlps nb203nf_izH1(%esp),%xmm6
        movlps nb203nf_izH2(%esp),%xmm7
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
        mulps nb203nf_krf(%esp),%xmm0
        movaps %xmm0,nb203nf_krsqM(%esp)

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb203nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb203nf_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        ## a little trick to avoid NaNs: 
        ## positions 0,2,and 3 are valid, but not 1. 
        ## If it contains NaN it doesnt help to mult by 0, 
        ## So we shuffle it and copy pos 0 to pos1! 
        shufps $224,%xmm0,%xmm0 ## constant 11100000

        ## xmm0=rinv 
        addps  nb203nf_krsqM(%esp),%xmm0
        subps  nb203nf_crf(%esp),%xmm0   ## xmm0=rinv+ krsq-crf 
        mulps  nb203nf_qqM(%esp),%xmm0          ## xmm0=vcoul 
        addps  nb203nf_vctot(%esp),%xmm0
        movaps %xmm0,nb203nf_vctot(%esp)

        decl nb203nf_innerk(%esp)
        jz    _nb_kernel203nf_ia32_sse.nb203nf_updateouterdata
        jmp   _nb_kernel203nf_ia32_sse.nb203nf_odd_loop
_nb_kernel203nf_ia32_sse.nb203nf_updateouterdata: 
        ## get n from stack
        movl nb203nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb203nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb203nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb203nf_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb203nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel203nf_ia32_sse.nb203nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb203nf_n(%esp)
        jmp _nb_kernel203nf_ia32_sse.nb203nf_outer
_nb_kernel203nf_ia32_sse.nb203nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb203nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel203nf_ia32_sse.nb203nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel203nf_ia32_sse.nb203nf_threadloop
_nb_kernel203nf_ia32_sse.nb203nf_end: 
        emms

        movl nb203nf_nouter(%esp),%eax
        movl nb203nf_ninner(%esp),%ebx
        movl nb203nf_outeriter(%ebp),%ecx
        movl nb203nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb203nf_salign(%esp),%eax
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

