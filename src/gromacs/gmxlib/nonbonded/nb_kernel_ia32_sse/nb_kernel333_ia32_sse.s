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



.globl nb_kernel333_ia32_sse
.globl _nb_kernel333_ia32_sse
nb_kernel333_ia32_sse:  
_nb_kernel333_ia32_sse: 
.set nb333_p_nri, 8
.set nb333_iinr, 12
.set nb333_jindex, 16
.set nb333_jjnr, 20
.set nb333_shift, 24
.set nb333_shiftvec, 28
.set nb333_fshift, 32
.set nb333_gid, 36
.set nb333_pos, 40
.set nb333_faction, 44
.set nb333_charge, 48
.set nb333_p_facel, 52
.set nb333_argkrf, 56
.set nb333_argcrf, 60
.set nb333_Vc, 64
.set nb333_type, 68
.set nb333_p_ntype, 72
.set nb333_vdwparam, 76
.set nb333_Vvdw, 80
.set nb333_p_tabscale, 84
.set nb333_VFtab, 88
.set nb333_invsqrta, 92
.set nb333_dvda, 96
.set nb333_p_gbtabscale, 100
.set nb333_GBtab, 104
.set nb333_p_nthreads, 108
.set nb333_count, 112
.set nb333_mtx, 116
.set nb333_outeriter, 120
.set nb333_inneriter, 124
.set nb333_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb333_ixO, 0
.set nb333_iyO, 16
.set nb333_izO, 32
.set nb333_ixH1, 48
.set nb333_iyH1, 64
.set nb333_izH1, 80
.set nb333_ixH2, 96
.set nb333_iyH2, 112
.set nb333_izH2, 128
.set nb333_ixM, 144
.set nb333_iyM, 160
.set nb333_izM, 176
.set nb333_iqM, 192
.set nb333_iqH, 208
.set nb333_dxO, 224
.set nb333_dyO, 240
.set nb333_dzO, 256
.set nb333_dxH1, 272
.set nb333_dyH1, 288
.set nb333_dzH1, 304
.set nb333_dxH2, 320
.set nb333_dyH2, 336
.set nb333_dzH2, 352
.set nb333_dxM, 368
.set nb333_dyM, 384
.set nb333_dzM, 400
.set nb333_qqM, 416
.set nb333_qqH, 432
.set nb333_rinvO, 448
.set nb333_rinvH1, 464
.set nb333_rinvH2, 480
.set nb333_rinvM, 496
.set nb333_rO, 512
.set nb333_rH1, 528
.set nb333_rH2, 544
.set nb333_rM, 560
.set nb333_tsc, 576
.set nb333_two, 592
.set nb333_c6, 608
.set nb333_c12, 624
.set nb333_vctot, 640
.set nb333_Vvdwtot, 656
.set nb333_fixO, 672
.set nb333_fiyO, 688
.set nb333_fizO, 704
.set nb333_fixH1, 720
.set nb333_fiyH1, 736
.set nb333_fizH1, 752
.set nb333_fixH2, 768
.set nb333_fiyH2, 784
.set nb333_fizH2, 800
.set nb333_fixM, 816
.set nb333_fiyM, 832
.set nb333_fizM, 848
.set nb333_fjx, 864
.set nb333_fjy, 880
.set nb333_fjz, 896
.set nb333_half, 912
.set nb333_three, 928
.set nb333_is3, 944
.set nb333_ii3, 948
.set nb333_ntia, 952
.set nb333_innerjjnr, 956
.set nb333_innerk, 960
.set nb333_n, 964
.set nb333_nn1, 968
.set nb333_nri, 972
.set nb333_nouter, 976
.set nb333_ninner, 980
.set nb333_salign, 984
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $988,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb333_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb333_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb333_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb333_nouter(%esp)
        movl %eax,nb333_ninner(%esp)


        movl nb333_p_tabscale(%ebp),%eax
        movss (%eax),%xmm5
        shufps $0,%xmm5,%xmm5
        movaps %xmm5,nb333_tsc(%esp)
        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb333_half(%esp)
        movss nb333_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb333_half(%esp)
        movaps %xmm2,nb333_two(%esp)
        movaps %xmm3,nb333_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb333_iinr(%ebp),%ecx             ## ecx = pointer into iinr[]    
        movl  (%ecx),%ebx               ## ebx =ii 

        movl  nb333_charge(%ebp),%edx
        movss 4(%edx,%ebx,4),%xmm4
        movss 12(%edx,%ebx,4),%xmm3
        movl nb333_p_facel(%ebp),%esi
        movss (%esi),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb333_iqM(%esp)
        movaps %xmm4,nb333_iqH(%esp)

        movl  nb333_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb333_p_ntype(%ebp),%edi
        imull (%edi),%ecx       ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb333_ntia(%esp)

_nb_kernel333_ia32_sse.nb333_threadloop: 
        movl  nb333_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel333_ia32_sse.nb333_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel333_ia32_sse.nb333_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb333_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb333_n(%esp)
        movl %ebx,nb333_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel333_ia32_sse.nb333_outerstart
        jmp _nb_kernel333_ia32_sse.nb333_end

_nb_kernel333_ia32_sse.nb333_outerstart: 
        ## ebx contains number of outer iterations
        addl nb333_nouter(%esp),%ebx
        movl %ebx,nb333_nouter(%esp)

_nb_kernel333_ia32_sse.nb333_outer: 
        movl  nb333_shift(%ebp),%eax            ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx        ## ebx=3*is 
        movl  %ebx,nb333_is3(%esp)      ## store is3 

        movl  nb333_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb333_iinr(%ebp),%ecx             ## ecx = pointer into iinr[]    
        movl  (%ecx,%esi,4),%ebx                ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        movaps %xmm0,%xmm6
        movaps %xmm1,%xmm7

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb333_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb333_ii3(%esp)

        addss (%eax,%ebx,4),%xmm3       ## ox
        addss 4(%eax,%ebx,4),%xmm4     ## oy
        addss 8(%eax,%ebx,4),%xmm5     ## oz
        addss 12(%eax,%ebx,4),%xmm6    ## h1x
        addss 16(%eax,%ebx,4),%xmm7    ## h1y
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        shufps $0,%xmm6,%xmm6
        shufps $0,%xmm7,%xmm7
        movaps %xmm3,nb333_ixO(%esp)
        movaps %xmm4,nb333_iyO(%esp)
        movaps %xmm5,nb333_izO(%esp)
        movaps %xmm6,nb333_ixH1(%esp)
        movaps %xmm7,nb333_iyH1(%esp)

        movss %xmm2,%xmm6
        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 20(%eax,%ebx,4),%xmm6    ## h1z
        addss 24(%eax,%ebx,4),%xmm0    ## h2x
        addss 28(%eax,%ebx,4),%xmm1    ## h2y
        addss 32(%eax,%ebx,4),%xmm2    ## h2z
        addss 36(%eax,%ebx,4),%xmm3    ## mx
        addss 40(%eax,%ebx,4),%xmm4    ## my
        addss 44(%eax,%ebx,4),%xmm5    ## mz

        shufps $0,%xmm6,%xmm6
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm6,nb333_izH1(%esp)
        movaps %xmm0,nb333_ixH2(%esp)
        movaps %xmm1,nb333_iyH2(%esp)
        movaps %xmm2,nb333_izH2(%esp)
        movaps %xmm3,nb333_ixM(%esp)
        movaps %xmm4,nb333_iyM(%esp)
        movaps %xmm5,nb333_izM(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb333_vctot(%esp)
        movaps %xmm4,nb333_Vvdwtot(%esp)
        movaps %xmm4,nb333_fixO(%esp)
        movaps %xmm4,nb333_fiyO(%esp)
        movaps %xmm4,nb333_fizO(%esp)
        movaps %xmm4,nb333_fixH1(%esp)
        movaps %xmm4,nb333_fiyH1(%esp)
        movaps %xmm4,nb333_fizH1(%esp)
        movaps %xmm4,nb333_fixH2(%esp)
        movaps %xmm4,nb333_fiyH2(%esp)
        movaps %xmm4,nb333_fizH2(%esp)
        movaps %xmm4,nb333_fixM(%esp)
        movaps %xmm4,nb333_fiyM(%esp)
        movaps %xmm4,nb333_fizM(%esp)

        movl  nb333_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx                ## jindex[n] 
        movl  4(%eax,%esi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movl  nb333_pos(%ebp),%esi
        movl  nb333_faction(%ebp),%edi
        movl  nb333_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb333_innerjjnr(%esp)        ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb333_ninner(%esp),%ecx
        movl  %ecx,nb333_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb333_innerk(%esp)   ## number of innerloop atoms 
        jge   _nb_kernel333_ia32_sse.nb333_unroll_loop
        jmp   _nb_kernel333_ia32_sse.nb333_odd_inner
_nb_kernel333_ia32_sse.nb333_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb333_innerjjnr(%esp),%edx        ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx             ## eax-edx=jnr1-4

        addl $16,nb333_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb333_charge(%ebp),%esi    ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4              ## and in xmm4 
        mulps  nb333_iqM(%esp),%xmm3
        mulps  nb333_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps  %xmm3,nb333_qqM(%esp)
        movaps  %xmm4,nb333_qqH(%esp)

        movl nb333_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb333_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb333_ntia(%esp),%edi
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

        movaps %xmm4,nb333_c6(%esp)
        movaps %xmm6,nb333_c12(%esp)

        movl nb333_pos(%ebp),%esi       ## base of pos[] 

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

        ## move ixO-izO to xmm4-xmm6 
        movaps nb333_ixO(%esp),%xmm4
        movaps nb333_iyO(%esp),%xmm5
        movaps nb333_izO(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb333_dxO(%esp)
        movaps %xmm5,nb333_dyO(%esp)
        movaps %xmm6,nb333_dzO(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        movaps %xmm4,%xmm7
        ## rsqO in xmm7

        ## move ixH1-izH1 to xmm4-xmm6 
        movaps nb333_ixH1(%esp),%xmm4
        movaps nb333_iyH1(%esp),%xmm5
        movaps nb333_izH1(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb333_dxH1(%esp)
        movaps %xmm5,nb333_dyH1(%esp)
        movaps %xmm6,nb333_dzH1(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm6
        addps %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movaps nb333_ixH2(%esp),%xmm3
        movaps nb333_iyH2(%esp),%xmm4
        movaps nb333_izH2(%esp),%xmm5

        ## calc dr 
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        ## store dr 
        movaps %xmm3,nb333_dxH2(%esp)
        movaps %xmm4,nb333_dyH2(%esp)
        movaps %xmm5,nb333_dzH2(%esp)
        ## square it 
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm4,%xmm5
        addps %xmm3,%xmm5

        ## move ixM-izM to xmm2-xmm4  
        movaps nb333_iyM(%esp),%xmm3
        movaps nb333_izM(%esp),%xmm4
        subps  %xmm1,%xmm3
        subps  %xmm2,%xmm4
        movaps nb333_ixM(%esp),%xmm2
        subps  %xmm0,%xmm2

        ## store dr 
        movaps %xmm2,nb333_dxM(%esp)
        movaps %xmm3,nb333_dyM(%esp)
        movaps %xmm4,nb333_dzM(%esp)
        ## square it 
        mulps %xmm2,%xmm2
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        addps %xmm3,%xmm4
        addps %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7

        ## rsqH1 - seed in xmm2 
        rsqrtps %xmm6,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb333_three(%esp),%xmm0
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb333_half(%esp),%xmm0
        movaps  %xmm0,nb333_rinvH1(%esp)        ## rinvH1  
        mulps   %xmm0,%xmm6
        movaps  %xmm6,nb333_rH1(%esp)

        ## rsqH2 - seed to xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb333_three(%esp),%xmm0
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb333_half(%esp),%xmm0
        movaps  %xmm0,nb333_rinvH2(%esp)        ## rinvH2 
        mulps   %xmm0,%xmm5
        movaps  %xmm5,nb333_rH2(%esp)

        ## rsqM - seed to xmm2 
        rsqrtps %xmm4,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb333_three(%esp),%xmm0
        mulps   %xmm4,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb333_half(%esp),%xmm0
        movaps  %xmm0,nb333_rinvM(%esp)         ## rinvM 
        mulps   %xmm0,%xmm4
        movaps  %xmm4,nb333_rM(%esp)

        ## Do the O LJ table interaction directly.
        rsqrtps %xmm7,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb333_three(%esp),%xmm0
        mulps   %xmm7,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb333_half(%esp),%xmm0   ## rinv

        movaps %xmm0,%xmm1
        mulps  %xmm7,%xmm1      ## xmm1=r
        mulps  nb333_tsc(%esp),%xmm1   ## r*tabscale

        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %eax,%mm0
        movd %ebx,%mm1
        movd %ecx,%mm2
        movd %edx,%mm3

        movl nb333_VFtab(%ebp),%esi
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

        ## load dispersion table data into xmm4-xmm7 
        movlps 16(%esi,%eax,4),%xmm5
        movlps 16(%esi,%ecx,4),%xmm7
        movhps 16(%esi,%ebx,4),%xmm5
        movhps 16(%esi,%edx,4),%xmm7    ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 24(%esi,%eax,4),%xmm7
        movlps 24(%esi,%ecx,4),%xmm3
        movhps 24(%esi,%ebx,4),%xmm7
        movhps 24(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101

        ## dispersion table YFGH ready in xmm4-xmm7
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb333_two(%esp),%xmm7            ## two*Heps2 
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb333_c6(%esp),%xmm4
        mulps  %xmm4,%xmm7      ## fijD 
        mulps  %xmm4,%xmm5      ## Vvdw6 

        ## put scalar force on stack (borrow rinvO) 
        ## Update Vvdwtot directly      
        addps  nb333_Vvdwtot(%esp),%xmm5
        movaps %xmm7,nb333_rinvO(%esp)   ## fscal 
        movaps %xmm5,nb333_Vvdwtot(%esp)

        ## load repulsion table data into xmm4-xmm7
        movlps 32(%esi,%eax,4),%xmm5
        movlps 32(%esi,%ecx,4),%xmm7
        movhps 32(%esi,%ebx,4),%xmm5
        movhps 32(%esi,%edx,4),%xmm7    ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 40(%esi,%eax,4),%xmm7
        movlps 40(%esi,%ecx,4),%xmm3
        movhps 40(%esi,%ebx,4),%xmm7
        movhps 40(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101

        ## repulsion table YFGH ready in xmm4-xmm7

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb333_two(%esp),%xmm7            ## two*Heps2 
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb333_c12(%esp),%xmm4
        mulps  %xmm4,%xmm7 ## fijR 
        mulps  %xmm4,%xmm5 ## Vvdw12 
        addps  nb333_rinvO(%esp),%xmm7   ## fscal was temp. stored in rinvO

        addps  nb333_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb333_Vvdwtot(%esp)

        xorps %xmm1,%xmm1
        mulps nb333_tsc(%esp),%xmm7
        mulps %xmm0,%xmm7
        subps  %xmm7,%xmm1      ## fscal
        movaps nb333_dxO(%esp),%xmm3
        movaps nb333_dyO(%esp),%xmm4
        movaps nb333_dzO(%esp),%xmm5
        mulps  %xmm1,%xmm3
        mulps  %xmm1,%xmm4
        mulps  %xmm1,%xmm5      ## tx in xmm3-xmm5

        ## update O forces 
        movaps nb333_fixO(%esp),%xmm0
        movaps nb333_fiyO(%esp),%xmm1
        movaps nb333_fizO(%esp),%xmm2
        addps  %xmm3,%xmm0
        addps  %xmm4,%xmm1
        addps  %xmm5,%xmm2
        movaps %xmm0,nb333_fixO(%esp)
        movaps %xmm1,nb333_fiyO(%esp)
        movaps %xmm2,nb333_fizO(%esp)
        ## update j forces with water O 
        movaps %xmm3,nb333_fjx(%esp)
        movaps %xmm4,nb333_fjy(%esp)
        movaps %xmm5,nb333_fjz(%esp)

        ## Do H1 interaction
        movl nb333_VFtab(%ebp),%esi

        movaps nb333_rH1(%esp),%xmm7
        mulps   nb333_tsc(%esp),%xmm7
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
        mulps  nb333_two(%esp),%xmm7            ## two*Heps2 
        movaps nb333_qqH(%esp),%xmm0
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm0,%xmm7 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm7 fijC 
        ## increment vcoul 
        xorps  %xmm4,%xmm4
        addps  nb333_vctot(%esp),%xmm5
        mulps  nb333_rinvH1(%esp),%xmm7
        movaps %xmm5,nb333_vctot(%esp)
        mulps  nb333_tsc(%esp),%xmm7
        subps %xmm7,%xmm4

        movaps nb333_dxH1(%esp),%xmm0
        movaps nb333_dyH1(%esp),%xmm1
        movaps nb333_dzH1(%esp),%xmm2
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2

        ## update H1 forces 
        movaps nb333_fixH1(%esp),%xmm3
        movaps nb333_fiyH1(%esp),%xmm4
        movaps nb333_fizH1(%esp),%xmm7
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm7
        movaps %xmm3,nb333_fixH1(%esp)
        movaps %xmm4,nb333_fiyH1(%esp)
        movaps %xmm7,nb333_fizH1(%esp)
        ## update j forces with water H1 
        addps  nb333_fjx(%esp),%xmm0
        addps  nb333_fjy(%esp),%xmm1
        addps  nb333_fjz(%esp),%xmm2
        movaps %xmm0,nb333_fjx(%esp)
        movaps %xmm1,nb333_fjy(%esp)
        movaps %xmm2,nb333_fjz(%esp)

        ## Done with H1, do H2 interactions 
        movaps nb333_rH2(%esp),%xmm7
        mulps   nb333_tsc(%esp),%xmm7
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
        mulps  nb333_two(%esp),%xmm7            ## two*Heps2 
        movaps nb333_qqH(%esp),%xmm0
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm0,%xmm7 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm0 fijC 
        ## increment vcoul 
        xorps  %xmm4,%xmm4
        addps  nb333_vctot(%esp),%xmm5
        mulps  nb333_rinvH2(%esp),%xmm7
        movaps %xmm5,nb333_vctot(%esp)
        mulps  nb333_tsc(%esp),%xmm7
        subps  %xmm7,%xmm4

        movaps nb333_dxH2(%esp),%xmm0
        movaps nb333_dyH2(%esp),%xmm1
        movaps nb333_dzH2(%esp),%xmm2
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2

        movd %mm0,%eax
        movd %mm1,%ebx
        movd %mm2,%ecx
        movd %mm3,%edx

        ## update H2 forces 
        movaps nb333_fixH2(%esp),%xmm3
        movaps nb333_fiyH2(%esp),%xmm4
        movaps nb333_fizH2(%esp),%xmm7
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm7
        movaps %xmm3,nb333_fixH2(%esp)
        movaps %xmm4,nb333_fiyH2(%esp)
        movaps %xmm7,nb333_fizH2(%esp)
        addps nb333_fjx(%esp),%xmm0
        addps nb333_fjy(%esp),%xmm1
        addps nb333_fjz(%esp),%xmm2
        movaps %xmm0,nb333_fjx(%esp)
        movaps %xmm1,nb333_fjy(%esp)
        movaps %xmm2,nb333_fjz(%esp)

        ## Done with H2, do M interactions 
        movaps nb333_rM(%esp),%xmm7
        mulps   nb333_tsc(%esp),%xmm7
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
        mulps  nb333_two(%esp),%xmm7            ## two*Heps2 
        movaps nb333_qqM(%esp),%xmm0
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm0,%xmm7 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm0 fijC 
        ## increment vcoul 
        xorps  %xmm4,%xmm4
        addps  nb333_vctot(%esp),%xmm5
        mulps  nb333_rinvM(%esp),%xmm7
        movaps %xmm5,nb333_vctot(%esp)
        mulps  nb333_tsc(%esp),%xmm7
        subps  %xmm7,%xmm4

        movaps nb333_dxM(%esp),%xmm0
        movaps nb333_dyM(%esp),%xmm1
        movaps nb333_dzM(%esp),%xmm2
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2

        movd %mm0,%eax
        movd %mm1,%ebx
        movd %mm2,%ecx
        movd %mm3,%edx

        ## update M forces 
        movaps nb333_fixM(%esp),%xmm3
        movaps nb333_fiyM(%esp),%xmm4
        movaps nb333_fizM(%esp),%xmm7
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm7
        movaps %xmm3,nb333_fixM(%esp)
        movaps %xmm4,nb333_fiyM(%esp)
        movaps %xmm7,nb333_fizM(%esp)

        movl nb333_faction(%ebp),%edi
        ## update j forces from stored values
        addps nb333_fjx(%esp),%xmm0
        addps nb333_fjy(%esp),%xmm1
        addps nb333_fjz(%esp),%xmm2

        movlps (%edi,%eax,4),%xmm4
        movlps (%edi,%ecx,4),%xmm7
        movhps (%edi,%ebx,4),%xmm4
        movhps (%edi,%edx,4),%xmm7

        movd %mm0,%eax
        movd %mm1,%ebx
        movd %mm2,%ecx
        movd %mm3,%edx

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
        subl $4,nb333_innerk(%esp)
        jl    _nb_kernel333_ia32_sse.nb333_odd_inner
        jmp   _nb_kernel333_ia32_sse.nb333_unroll_loop
_nb_kernel333_ia32_sse.nb333_odd_inner: 
        addl $4,nb333_innerk(%esp)
        jnz   _nb_kernel333_ia32_sse.nb333_odd_loop
        jmp   _nb_kernel333_ia32_sse.nb333_updateouterdata
_nb_kernel333_ia32_sse.nb333_odd_loop: 
        movl  nb333_innerjjnr(%esp),%edx        ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb333_innerjjnr(%esp)

        xorps %xmm4,%xmm4       ## clear reg.
        movss nb333_iqM(%esp),%xmm4
        movl nb333_charge(%ebp),%esi
        movhps nb333_iqH(%esp),%xmm4    ## [qM  0  qH  qH] 
        shufps $41,%xmm4,%xmm4 ## [0 qH qH qM]

        movss (%esi,%eax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb333_qqM(%esp)    ## use dummy qq for storage 

        xorps %xmm6,%xmm6
        movl nb333_type(%ebp),%esi
        movl (%esi,%eax,4),%ebx
        movl nb333_vdwparam(%ebp),%esi
        shll %ebx
        addl nb333_ntia(%esp),%ebx
        movlps (%esi,%ebx,4),%xmm6
        movaps %xmm6,%xmm7
        shufps $252,%xmm6,%xmm6 ## constant 11111100
        shufps $253,%xmm7,%xmm7 ## constant 11111101
        movaps %xmm6,nb333_c6(%esp)
        movaps %xmm7,nb333_c12(%esp)

        movl nb333_pos(%ebp),%esi
        leal (%eax,%eax,2),%eax

        movss nb333_ixO(%esp),%xmm3
        movss nb333_iyO(%esp),%xmm4
        movss nb333_izO(%esp),%xmm5
        movss nb333_ixH1(%esp),%xmm0
        movss nb333_iyH1(%esp),%xmm1
        movss nb333_izH1(%esp),%xmm2
        unpcklps nb333_ixH2(%esp),%xmm3         ## ixO ixH2 - -
        unpcklps nb333_iyH2(%esp),%xmm4         ## iyO iyH2 - -
        unpcklps nb333_izH2(%esp),%xmm5         ## izO izH2 - -
        unpcklps nb333_ixM(%esp),%xmm0          ## ixH1 ixM - -
        unpcklps nb333_iyM(%esp),%xmm1          ## iyH1 iyM - -
        unpcklps nb333_izM(%esp),%xmm2          ## izH1 izM - -
        unpcklps %xmm0,%xmm3    ## ixO ixH1 ixH2 ixM
        unpcklps %xmm1,%xmm4    ## same for y
        unpcklps %xmm2,%xmm5    ## same for z

        ## move j coords to xmm0-xmm2 
        movss (%esi,%eax,4),%xmm0
        movss 4(%esi,%eax,4),%xmm1
        movss 8(%esi,%eax,4),%xmm2
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        ## use O distances for storage
        movaps %xmm3,nb333_dxO(%esp)
        movaps %xmm4,nb333_dyO(%esp)
        movaps %xmm5,nb333_dzO(%esp)

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
        movaps nb333_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb333_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv

        movaps %xmm0,nb333_rinvM(%esp)
        mulps  %xmm0,%xmm4      ## r     
        mulps nb333_tsc(%esp),%xmm4

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


        movl nb333_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm7,%edx

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx
        leal  (%edx,%edx,2),%edx

        ## first do LJ table for O
        ## load dispersion table data into xmm4

        movlps 16(%esi,%eax,4),%xmm4
        movlps 24(%esi,%eax,4),%xmm6
        movaps %xmm4,%xmm5
        movaps %xmm6,%xmm7
        shufps $0x1,%xmm5,%xmm5
        shufps $0x1,%xmm7,%xmm7

        ## dispersion table YFGH ready in xmm4-xmm7
        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp 
        mulss  nb333_two(%esp),%xmm7            ## two*Heps2 
        addss  %xmm6,%xmm7
        addss  %xmm5,%xmm7 ## xmm7=FF 
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb333_c6(%esp),%xmm4
        mulss  %xmm4,%xmm7      ## fijD 
        mulss  %xmm4,%xmm5      ## Vvdw6 

        ## save scalar force in xmm3. Update Vvdwtot directly 
        addss  nb333_Vvdwtot(%esp),%xmm5
        xorps %xmm3,%xmm3
        movss %xmm7,%xmm3 ## fscal 
        movss %xmm5,nb333_Vvdwtot(%esp)

        ## load repulsion table data into xmm4
        movlps 32(%esi,%eax,4),%xmm4
        movlps 40(%esi,%eax,4),%xmm6
        movaps %xmm4,%xmm5
        movaps %xmm6,%xmm7
        shufps $0x1,%xmm5,%xmm5
        shufps $0x1,%xmm7,%xmm7
        ## repulsion table YFGH ready in xmm4-xmm7

        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp 
        mulss  nb333_two(%esp),%xmm7            ## two*Heps2 
        addss  %xmm6,%xmm7
        addss  %xmm5,%xmm7 ## xmm7=FF 
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb333_c12(%esp),%xmm4
        mulss  %xmm4,%xmm7 ## fijR 
        mulss  %xmm4,%xmm5 ## Vvdw12 
        addss  %xmm7,%xmm3

        addss  nb333_Vvdwtot(%esp),%xmm5
        movss %xmm5,nb333_Vvdwtot(%esp)

        movaps %xmm3,nb333_rinvO(%esp) ## save fscal temp. in rinvO

        ## do the Coulomb interaction for H1,H2,M
        xorps  %xmm5,%xmm5
        movlps (%esi,%ecx,4),%xmm3      ## values= Y3 F3  -  - 
        movhps (%esi,%ebx,4),%xmm5      ## values= 0  0 Y2 F2
        movhps (%esi,%edx,4),%xmm3      ## values= Y3 F3 Y4 F4 

        movaps %xmm5,%xmm4              ## values= 0  0 Y2 F2 
        shufps $0x88,%xmm3,%xmm4       ## values= 0 Y2 Y3 Y3
        shufps $0xDD,%xmm3,%xmm5       ## values= 0 F2 F3 F4 

        xorps  %xmm7,%xmm7
        movlps 8(%esi,%ecx,4),%xmm3     ## values= G3 H3  -  - 
        movhps 8(%esi,%ebx,4),%xmm7     ## values= 0  0 G2 H2
        movhps 8(%esi,%edx,4),%xmm3     ## values= G3 H3 G4 H4 

        movaps %xmm7,%xmm6              ## values= 0  0 G2 H2 
        shufps $0x88,%xmm3,%xmm6       ## values= 0 G2 G3 G3
        shufps $0xDD,%xmm3,%xmm7       ## values= 0 H2 H3 H4 

        ## xmm4 =  0  Y2 Y3 Y4
        ## xmm5 =  0  F2 F3 F4
        ## xmm6 =  0  G2 G3 G4
        ## xmm7 =  0  H2 H3 H4

        ## coulomb table ready, in xmm4-xmm7      
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp        
        mulps  nb333_two(%esp),%xmm7            ## two*Heps2 
        movaps nb333_qqM(%esp),%xmm0
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm0 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm0 fijC 
        ## increment vcoul - then we can get rid of mm5 
        addps  nb333_vctot(%esp),%xmm5
        movaps %xmm5,nb333_vctot(%esp)

        addps nb333_rinvO(%esp),%xmm0 ## total fscal (temp. storage in rinvO)

        xorps %xmm4,%xmm4
        mulps  nb333_rinvM(%esp),%xmm0
        mulps  nb333_tsc(%esp),%xmm0
        subps  %xmm0,%xmm4

        movaps nb333_dxO(%esp),%xmm0
        movaps nb333_dyO(%esp),%xmm1
        movaps nb333_dzO(%esp),%xmm2

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2 ## xmm0-xmm2 now contains tx-tz (partial force)

        movss  nb333_fixO(%esp),%xmm3
        movss  nb333_fiyO(%esp),%xmm4
        movss  nb333_fizO(%esp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,nb333_fixO(%esp)
        movss  %xmm4,nb333_fiyO(%esp)
        movss  %xmm5,nb333_fizO(%esp)   ## updated the O force now do the H's

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        shufps $0x39,%xmm3,%xmm3 ## shift right 
        shufps $0x39,%xmm4,%xmm4
        shufps $0x39,%xmm5,%xmm5
        addss  nb333_fixH1(%esp),%xmm3
        addss  nb333_fiyH1(%esp),%xmm4
        addss  nb333_fizH1(%esp),%xmm5
        movss  %xmm3,nb333_fixH1(%esp)
        movss  %xmm4,nb333_fiyH1(%esp)
        movss  %xmm5,nb333_fizH1(%esp)          ## updated the H1 force 

        shufps $0x39,%xmm3,%xmm3
        shufps $0x39,%xmm4,%xmm4
        shufps $0x39,%xmm5,%xmm5
        addss  nb333_fixH2(%esp),%xmm3
        addss  nb333_fiyH2(%esp),%xmm4
        addss  nb333_fizH2(%esp),%xmm5
        movss  %xmm3,nb333_fixH2(%esp)
        movss  %xmm4,nb333_fiyH2(%esp)
        movss  %xmm5,nb333_fizH2(%esp)          ## updated the H2 force 

        movl nb333_faction(%ebp),%edi
        shufps $0x39,%xmm3,%xmm3
        shufps $0x39,%xmm4,%xmm4
        shufps $0x39,%xmm5,%xmm5
        addss  nb333_fixM(%esp),%xmm3
        addss  nb333_fiyM(%esp),%xmm4
        addss  nb333_fizM(%esp),%xmm5
        movss  %xmm3,nb333_fixM(%esp)
        movss  %xmm4,nb333_fiyM(%esp)
        movss  %xmm5,nb333_fizM(%esp)   ## updated the M force 

        movd %mm0,%eax
        ## the fj's - move in from mem start by acc. tx/ty/tz in xmm0, xmm1
        movlps (%edi,%eax,4),%xmm6
        movss  8(%edi,%eax,4),%xmm7

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movaps  %xmm3,%xmm0
        movaps  %xmm4,%xmm1
        movaps  %xmm5,%xmm2

        shufps $0x39,%xmm3,%xmm3 ## shift right 
        shufps $0x39,%xmm4,%xmm4
        shufps $0x39,%xmm5,%xmm5
        addss  %xmm3,%xmm0
        addss  %xmm4,%xmm1
        addss  %xmm5,%xmm2
        unpcklps %xmm1,%xmm0    ## x,y sum in xmm0, z sum in xmm2

        subps    %xmm0,%xmm6
        subss    %xmm2,%xmm7

        movlps %xmm6,(%edi,%eax,4)
        movss  %xmm7,8(%edi,%eax,4)

        decl nb333_innerk(%esp)
        jz    _nb_kernel333_ia32_sse.nb333_updateouterdata
        jmp   _nb_kernel333_ia32_sse.nb333_odd_loop
_nb_kernel333_ia32_sse.nb333_updateouterdata: 
        movl  nb333_ii3(%esp),%ecx
        movl  nb333_faction(%ebp),%edi
        movl  nb333_fshift(%ebp),%esi
        movl  nb333_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movaps nb333_fixO(%esp),%xmm0
        movaps nb333_fiyO(%esp),%xmm1
        movaps nb333_fizO(%esp),%xmm2

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
        movaps nb333_fixH1(%esp),%xmm0
        movaps nb333_fiyH1(%esp),%xmm1
        movaps nb333_fizH1(%esp),%xmm2

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
        movaps nb333_fixH2(%esp),%xmm0
        movaps nb333_fiyH2(%esp),%xmm1
        movaps nb333_fizH2(%esp),%xmm2

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

        ## accumulate Mi forces in xmm0, xmm1, xmm2 
        movaps nb333_fixM(%esp),%xmm0
        movaps nb333_fiyM(%esp),%xmm1
        movaps nb333_fizM(%esp),%xmm2

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
        movl nb333_n(%esp),%esi
        ## get group index for i particle 
        movl  nb333_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb333_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb333_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb333_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb333_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb333_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel333_ia32_sse.nb333_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb333_n(%esp)
        jmp _nb_kernel333_ia32_sse.nb333_outer
_nb_kernel333_ia32_sse.nb333_outerend: 
        ## check if more outer neighborlists remain
        movl  nb333_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel333_ia32_sse.nb333_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel333_ia32_sse.nb333_threadloop
_nb_kernel333_ia32_sse.nb333_end: 
        emms

        movl nb333_nouter(%esp),%eax
        movl nb333_ninner(%esp),%ebx
        movl nb333_outeriter(%ebp),%ecx
        movl nb333_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb333_salign(%esp),%eax
        addl %eax,%esp
        addl $988,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret


.globl nb_kernel333nf_ia32_sse
.globl _nb_kernel333nf_ia32_sse
nb_kernel333nf_ia32_sse:        
_nb_kernel333nf_ia32_sse:       
.set nb333nf_p_nri, 8
.set nb333nf_iinr, 12
.set nb333nf_jindex, 16
.set nb333nf_jjnr, 20
.set nb333nf_shift, 24
.set nb333nf_shiftvec, 28
.set nb333nf_fshift, 32
.set nb333nf_gid, 36
.set nb333nf_pos, 40
.set nb333nf_faction, 44
.set nb333nf_charge, 48
.set nb333nf_p_facel, 52
.set nb333nf_argkrf, 56
.set nb333nf_argcrf, 60
.set nb333nf_Vc, 64
.set nb333nf_type, 68
.set nb333nf_p_ntype, 72
.set nb333nf_vdwparam, 76
.set nb333nf_Vvdw, 80
.set nb333nf_p_tabscale, 84
.set nb333nf_VFtab, 88
.set nb333nf_invsqrta, 92
.set nb333nf_dvda, 96
.set nb333nf_p_gbtabscale, 100
.set nb333nf_GBtab, 104
.set nb333nf_p_nthreads, 108
.set nb333nf_count, 112
.set nb333nf_mtx, 116
.set nb333nf_outeriter, 120
.set nb333nf_inneriter, 124
.set nb333nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb333nf_ixO, 0
.set nb333nf_iyO, 16
.set nb333nf_izO, 32
.set nb333nf_ixH1, 48
.set nb333nf_iyH1, 64
.set nb333nf_izH1, 80
.set nb333nf_ixH2, 96
.set nb333nf_iyH2, 112
.set nb333nf_izH2, 128
.set nb333nf_ixM, 144
.set nb333nf_iyM, 160
.set nb333nf_izM, 176
.set nb333nf_iqM, 192
.set nb333nf_iqH, 208
.set nb333nf_qqM, 224
.set nb333nf_qqH, 240
.set nb333nf_rinvO, 256
.set nb333nf_rinvH1, 272
.set nb333nf_rinvH2, 288
.set nb333nf_rinvM, 304
.set nb333nf_rO, 320
.set nb333nf_rH1, 336
.set nb333nf_rH2, 352
.set nb333nf_rM, 368
.set nb333nf_tsc, 384
.set nb333nf_c6, 400
.set nb333nf_c12, 416
.set nb333nf_vctot, 432
.set nb333nf_Vvdwtot, 448
.set nb333nf_half, 464
.set nb333nf_three, 480
.set nb333nf_is3, 496
.set nb333nf_ii3, 500
.set nb333nf_ntia, 504
.set nb333nf_innerjjnr, 508
.set nb333nf_innerk, 512
.set nb333nf_n, 516
.set nb333nf_nn1, 520
.set nb333nf_nri, 524
.set nb333nf_nouter, 528
.set nb333nf_ninner, 532
.set nb333nf_salign, 536
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $540,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb333nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb333nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb333nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb333nf_nouter(%esp)
        movl %eax,nb333nf_ninner(%esp)


        movl nb333nf_p_tabscale(%ebp),%eax
        movss (%eax),%xmm5
        shufps $0,%xmm5,%xmm5
        movaps %xmm5,nb333nf_tsc(%esp)
        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb333nf_half(%esp)
        movss nb333nf_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb333nf_half(%esp)
        movaps %xmm3,nb333nf_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb333nf_iinr(%ebp),%ecx           ## ecx = pointer into iinr[]    
        movl  (%ecx),%ebx               ## ebx =ii 

        movl  nb333nf_charge(%ebp),%edx
        movss 4(%edx,%ebx,4),%xmm4
        movss 12(%edx,%ebx,4),%xmm3
        movl nb333nf_p_facel(%ebp),%esi
        movss (%esi),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb333nf_iqM(%esp)
        movaps %xmm4,nb333nf_iqH(%esp)

        movl  nb333nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb333nf_p_ntype(%ebp),%edi
        imull (%edi),%ecx       ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb333nf_ntia(%esp)

_nb_kernel333nf_ia32_sse.nb333nf_threadloop: 
        movl  nb333nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel333nf_ia32_sse.nb333nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel333nf_ia32_sse.nb333nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb333nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb333nf_n(%esp)
        movl %ebx,nb333nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel333nf_ia32_sse.nb333nf_outerstart
        jmp _nb_kernel333nf_ia32_sse.nb333nf_end

_nb_kernel333nf_ia32_sse.nb333nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb333nf_nouter(%esp),%ebx
        movl %ebx,nb333nf_nouter(%esp)

_nb_kernel333nf_ia32_sse.nb333nf_outer: 
        movl  nb333nf_shift(%ebp),%eax          ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx        ## ebx=3*is 
        movl  %ebx,nb333nf_is3(%esp)            ## store is3 

        movl  nb333nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb333nf_iinr(%ebp),%ecx           ## ecx = pointer into iinr[]    
        movl  (%ecx,%esi,4),%ebx                ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        movaps %xmm0,%xmm6
        movaps %xmm1,%xmm7

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb333nf_pos(%ebp),%eax    ## eax = base of pos[]  
        movl  %ebx,nb333nf_ii3(%esp)

        addss (%eax,%ebx,4),%xmm3       ## ox
        addss 4(%eax,%ebx,4),%xmm4     ## oy
        addss 8(%eax,%ebx,4),%xmm5     ## oz
        addss 12(%eax,%ebx,4),%xmm6    ## h1x
        addss 16(%eax,%ebx,4),%xmm7    ## h1y
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        shufps $0,%xmm6,%xmm6
        shufps $0,%xmm7,%xmm7
        movaps %xmm3,nb333nf_ixO(%esp)
        movaps %xmm4,nb333nf_iyO(%esp)
        movaps %xmm5,nb333nf_izO(%esp)
        movaps %xmm6,nb333nf_ixH1(%esp)
        movaps %xmm7,nb333nf_iyH1(%esp)

        movss %xmm2,%xmm6
        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 20(%eax,%ebx,4),%xmm6    ## h1z
        addss 24(%eax,%ebx,4),%xmm0    ## h2x
        addss 28(%eax,%ebx,4),%xmm1    ## h2y
        addss 32(%eax,%ebx,4),%xmm2    ## h2z
        addss 36(%eax,%ebx,4),%xmm3    ## mx
        addss 40(%eax,%ebx,4),%xmm4    ## my
        addss 44(%eax,%ebx,4),%xmm5    ## mz

        shufps $0,%xmm6,%xmm6
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm6,nb333nf_izH1(%esp)
        movaps %xmm0,nb333nf_ixH2(%esp)
        movaps %xmm1,nb333nf_iyH2(%esp)
        movaps %xmm2,nb333nf_izH2(%esp)
        movaps %xmm3,nb333nf_ixM(%esp)
        movaps %xmm4,nb333nf_iyM(%esp)
        movaps %xmm5,nb333nf_izM(%esp)

        ## clear vctot 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb333nf_vctot(%esp)
        movaps %xmm4,nb333nf_Vvdwtot(%esp)

        movl  nb333nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx                ## jindex[n] 
        movl  4(%eax,%esi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movl  nb333nf_pos(%ebp),%esi
        movl  nb333nf_faction(%ebp),%edi
        movl  nb333nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb333nf_innerjjnr(%esp)      ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb333nf_ninner(%esp),%ecx
        movl  %ecx,nb333nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb333nf_innerk(%esp)         ## number of innerloop atoms 
        jge   _nb_kernel333nf_ia32_sse.nb333nf_unroll_loop
        jmp   _nb_kernel333nf_ia32_sse.nb333nf_odd_inner
_nb_kernel333nf_ia32_sse.nb333nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb333nf_innerjjnr(%esp),%edx      ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx             ## eax-edx=jnr1-4

        addl $16,nb333nf_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb333nf_charge(%ebp),%esi  ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4              ## and in xmm4 
        mulps  nb333nf_iqM(%esp),%xmm3
        mulps  nb333nf_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps  %xmm3,nb333nf_qqM(%esp)
        movaps  %xmm4,nb333nf_qqH(%esp)

        movl nb333nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb333nf_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb333nf_ntia(%esp),%edi
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

        movaps %xmm4,nb333nf_c6(%esp)
        movaps %xmm6,nb333nf_c12(%esp)

        movl nb333nf_pos(%ebp),%esi     ## base of pos[] 

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

        ## move ixO-izO to xmm4-xmm6 
        movaps nb333nf_ixO(%esp),%xmm4
        movaps nb333nf_iyO(%esp),%xmm5
        movaps nb333nf_izO(%esp),%xmm6

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
        movaps nb333nf_ixH1(%esp),%xmm4
        movaps nb333nf_iyH1(%esp),%xmm5
        movaps nb333nf_izH1(%esp),%xmm6

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
        movaps nb333nf_ixH2(%esp),%xmm3
        movaps nb333nf_iyH2(%esp),%xmm4
        movaps nb333nf_izH2(%esp),%xmm5

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

        ## move ixM-izM to xmm2-xmm4  
        movaps nb333nf_iyM(%esp),%xmm3
        movaps nb333nf_izM(%esp),%xmm4
        subps  %xmm1,%xmm3
        subps  %xmm2,%xmm4
        movaps nb333nf_ixM(%esp),%xmm2
        subps  %xmm0,%xmm2

        ## square it 
        mulps %xmm2,%xmm2
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        addps %xmm3,%xmm4
        addps %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7

        ## rsqH1 - seed in xmm2 
        rsqrtps %xmm6,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb333nf_three(%esp),%xmm0
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb333nf_half(%esp),%xmm0
        movaps  %xmm0,nb333nf_rinvH1(%esp)      ## rinvH1  
        mulps   %xmm0,%xmm6
        movaps  %xmm6,nb333nf_rH1(%esp)

        ## rsqH2 - seed to xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb333nf_three(%esp),%xmm0
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb333nf_half(%esp),%xmm0
        movaps  %xmm0,nb333nf_rinvH2(%esp)      ## rinvH2 
        mulps   %xmm0,%xmm5
        movaps  %xmm5,nb333nf_rH2(%esp)

        ## rsqM - seed to xmm2 
        rsqrtps %xmm4,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb333nf_three(%esp),%xmm0
        mulps   %xmm4,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb333nf_half(%esp),%xmm0
        movaps  %xmm0,nb333nf_rinvM(%esp)       ## rinvM 
        mulps   %xmm0,%xmm4
        movaps  %xmm4,nb333nf_rM(%esp)

        ## Do the O LJ table interaction directly.
        rsqrtps %xmm7,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb333nf_three(%esp),%xmm0
        mulps   %xmm7,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb333nf_half(%esp),%xmm0   ## rinv

        movaps %xmm0,%xmm1
        mulps  %xmm7,%xmm1      ## xmm1=r
        mulps  nb333nf_tsc(%esp),%xmm1   ## r*tabscale

        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %eax,%mm0
        movd %ebx,%mm1
        movd %ecx,%mm2
        movd %edx,%mm3

        movl nb333nf_VFtab(%ebp),%esi
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

        ## load dispersion table data into xmm4-xmm7
        movlps 16(%esi,%eax,4),%xmm5
        movlps 16(%esi,%ecx,4),%xmm7
        movhps 16(%esi,%ebx,4),%xmm5
        movhps 16(%esi,%edx,4),%xmm7    ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 24(%esi,%eax,4),%xmm7
        movlps 24(%esi,%ecx,4),%xmm3
        movhps 24(%esi,%ebx,4),%xmm7
        movhps 24(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101

        ## dispersion table YFGH ready in xmm4-xmm7
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb333nf_c6(%esp),%xmm4
        mulps  %xmm4,%xmm5      ## Vvdw6 

        ## Update Vvdwtot directly      
        addps  nb333nf_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb333nf_Vvdwtot(%esp)

        ## load repulsion table data into xmm4-xmm7
        movlps 32(%esi,%eax,4),%xmm5
        movlps 32(%esi,%ecx,4),%xmm7
        movhps 32(%esi,%ebx,4),%xmm5
        movhps 32(%esi,%edx,4),%xmm7    ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 40(%esi,%eax,4),%xmm7
        movlps 40(%esi,%ecx,4),%xmm3
        movhps 40(%esi,%ebx,4),%xmm7
        movhps 40(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101

        ## repulsion table YFGH ready in xmm4-xmm7      
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb333nf_c12(%esp),%xmm4
        mulps  %xmm4,%xmm5 ## Vvdw12 
        addps  nb333nf_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb333nf_Vvdwtot(%esp)

        ## Do H1 interaction
        movl nb333nf_VFtab(%ebp),%esi

        movaps nb333nf_rH1(%esp),%xmm7
        mulps   nb333nf_tsc(%esp),%xmm7
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
        movaps nb333nf_qqH(%esp),%xmm0
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        ## increment vcoul 
        addps  nb333nf_vctot(%esp),%xmm5
        movaps %xmm5,nb333nf_vctot(%esp)

        ## Done with H1, do H2 interactions 
        movaps nb333nf_rH2(%esp),%xmm7
        mulps   nb333nf_tsc(%esp),%xmm7
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
        movaps nb333nf_qqH(%esp),%xmm0
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        ## increment vcoul 
        addps  nb333nf_vctot(%esp),%xmm5
        movaps %xmm5,nb333nf_vctot(%esp)

        ## Done with H2, do M interactions 
        movaps nb333nf_rM(%esp),%xmm7
        mulps   nb333nf_tsc(%esp),%xmm7
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
        movaps nb333nf_qqM(%esp),%xmm0
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        ## increment vcoul 
        addps  nb333nf_vctot(%esp),%xmm5
        movaps %xmm5,nb333nf_vctot(%esp)
        ## should we do one more iteration? 
        subl $4,nb333nf_innerk(%esp)
        jl    _nb_kernel333nf_ia32_sse.nb333nf_odd_inner
        jmp   _nb_kernel333nf_ia32_sse.nb333nf_unroll_loop
_nb_kernel333nf_ia32_sse.nb333nf_odd_inner: 
        addl $4,nb333nf_innerk(%esp)
        jnz   _nb_kernel333nf_ia32_sse.nb333nf_odd_loop
        jmp   _nb_kernel333nf_ia32_sse.nb333nf_updateouterdata
_nb_kernel333nf_ia32_sse.nb333nf_odd_loop: 
        movl  nb333nf_innerjjnr(%esp),%edx      ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb333nf_innerjjnr(%esp)

        xorps %xmm4,%xmm4       ## clear reg.
        movss nb333nf_iqM(%esp),%xmm4
        movl nb333nf_charge(%ebp),%esi
        movhps nb333nf_iqH(%esp),%xmm4    ## [qM  0  qH  qH] 
        shufps $41,%xmm4,%xmm4 ## [0 qH qH qM]

        movss (%esi,%eax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb333nf_qqM(%esp)          ## use dummy qq for storage 

        xorps %xmm6,%xmm6
        movl nb333nf_type(%ebp),%esi
        movl (%esi,%eax,4),%ebx
        movl nb333nf_vdwparam(%ebp),%esi
        shll %ebx
        addl nb333nf_ntia(%esp),%ebx
        movlps (%esi,%ebx,4),%xmm6
        movaps %xmm6,%xmm7
        shufps $252,%xmm6,%xmm6 ## constant 11111100
        shufps $253,%xmm7,%xmm7 ## constant 11111101
        movaps %xmm6,nb333nf_c6(%esp)
        movaps %xmm7,nb333nf_c12(%esp)

        movl nb333nf_pos(%ebp),%esi
        leal (%eax,%eax,2),%eax

        movss nb333nf_ixO(%esp),%xmm3
        movss nb333nf_iyO(%esp),%xmm4
        movss nb333nf_izO(%esp),%xmm5
        movss nb333nf_ixH1(%esp),%xmm0
        movss nb333nf_iyH1(%esp),%xmm1
        movss nb333nf_izH1(%esp),%xmm2
        unpcklps nb333nf_ixH2(%esp),%xmm3       ## ixO ixH2 - -
        unpcklps nb333nf_iyH2(%esp),%xmm4       ## iyO iyH2 - -
        unpcklps nb333nf_izH2(%esp),%xmm5       ## izO izH2 - -
        unpcklps nb333nf_ixM(%esp),%xmm0        ## ixH1 ixM - -
        unpcklps nb333nf_iyM(%esp),%xmm1        ## iyH1 iyM - -
        unpcklps nb333nf_izM(%esp),%xmm2        ## izH1 izM - -
        unpcklps %xmm0,%xmm3    ## ixO ixH1 ixH2 ixM
        unpcklps %xmm1,%xmm4    ## same for y
        unpcklps %xmm2,%xmm5    ## same for z

        ## move j coords to xmm0-xmm2 
        movss (%esi,%eax,4),%xmm0
        movss 4(%esi,%eax,4),%xmm1
        movss 8(%esi,%eax,4),%xmm2
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

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
        movaps nb333nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb333nf_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv

        movaps %xmm0,nb333nf_rinvM(%esp)
        mulps  %xmm0,%xmm4      ## r     
        mulps nb333nf_tsc(%esp),%xmm4

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

        movl nb333nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm7,%edx

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx
        leal  (%edx,%edx,2),%edx

        ## first do LJ table for O
        ## load dispersion table data into xmm4
        movlps 16(%esi,%eax,4),%xmm4
        movlps 24(%esi,%eax,4),%xmm6
        movaps %xmm4,%xmm5
        movaps %xmm6,%xmm7
        shufps $0x1,%xmm5,%xmm5
        shufps $0x1,%xmm7,%xmm7

        ## dispersion table YFGH ready in xmm4-xmm7
        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp 
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb333nf_c6(%esp),%xmm4
        mulss  %xmm4,%xmm5      ## Vvdw6 
        addss  nb333nf_Vvdwtot(%esp),%xmm5
        movss %xmm5,nb333nf_Vvdwtot(%esp)

        ## load repulsion table data into xmm4
        movlps 32(%esi,%eax,4),%xmm4
        movlps 40(%esi,%eax,4),%xmm6
        movaps %xmm4,%xmm5
        movaps %xmm6,%xmm7
        shufps $0x1,%xmm5,%xmm5
        shufps $0x1,%xmm7,%xmm7
        ## repulsion table YFGH ready in xmm4-xmm7

        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp 
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb333nf_c12(%esp),%xmm4
        mulss  %xmm4,%xmm5 ## Vvdw12 
        addss  nb333nf_Vvdwtot(%esp),%xmm5
        movss %xmm5,nb333nf_Vvdwtot(%esp)
        ## do the Coulomb interaction for H1,H2,M
        xorps  %xmm5,%xmm5
        movlps (%esi,%ecx,4),%xmm3      ## values= Y3 F3  -  - 
        movhps (%esi,%ebx,4),%xmm5      ## values= 0  0 Y2 F2
        movhps (%esi,%edx,4),%xmm3      ## values= Y3 F3 Y4 F4 

        movaps %xmm5,%xmm4              ## values= 0  0 Y2 F2 
        shufps $0x88,%xmm3,%xmm4       ## values= 0 Y2 Y3 Y3
        shufps $0xDD,%xmm3,%xmm5       ## values= 0 F2 F3 F4 

        xorps  %xmm7,%xmm7
        movlps 8(%esi,%ecx,4),%xmm3     ## values= G3 H3  -  - 
        movhps 8(%esi,%ebx,4),%xmm7     ## values= 0  0 G2 H2
        movhps 8(%esi,%edx,4),%xmm3     ## values= G3 H3 G4 H4 

        movaps %xmm7,%xmm6              ## values= 0  0 G2 H2 
        shufps $0x88,%xmm3,%xmm6       ## values= 0 G2 G3 G3
        shufps $0xDD,%xmm3,%xmm7       ## values= 0 H2 H3 H4 

        ## xmm4 =  0  Y2 Y3 Y4
        ## xmm5 =  0  F2 F3 F4
        ## xmm6 =  0  G2 G3 G4
        ## xmm7 =  0  H2 H3 H4  
        ## coulomb table ready, in xmm4-xmm7      
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp        
        movaps nb333nf_qqM(%esp),%xmm0
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        ## increment vcoul
        addps  nb333nf_vctot(%esp),%xmm5
        movaps %xmm5,nb333nf_vctot(%esp)

        decl nb333nf_innerk(%esp)
        jz    _nb_kernel333nf_ia32_sse.nb333nf_updateouterdata
        jmp   _nb_kernel333nf_ia32_sse.nb333nf_odd_loop
_nb_kernel333nf_ia32_sse.nb333nf_updateouterdata: 
        ## get n from stack
        movl nb333nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb333nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb333nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb333nf_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb333nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb333nf_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb333nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel333nf_ia32_sse.nb333nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb333nf_n(%esp)
        jmp _nb_kernel333nf_ia32_sse.nb333nf_outer
_nb_kernel333nf_ia32_sse.nb333nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb333nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel333nf_ia32_sse.nb333nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel333nf_ia32_sse.nb333nf_threadloop
_nb_kernel333nf_ia32_sse.nb333nf_end: 
        emms

        movl nb333nf_nouter(%esp),%eax
        movl nb333nf_ninner(%esp),%ebx
        movl nb333nf_outeriter(%ebp),%ecx
        movl nb333nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb333nf_salign(%esp),%eax
        addl %eax,%esp
        addl $540,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret


