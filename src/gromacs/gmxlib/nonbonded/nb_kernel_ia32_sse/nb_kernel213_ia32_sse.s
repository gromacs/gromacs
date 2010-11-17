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


.globl nb_kernel213_ia32_sse
.globl _nb_kernel213_ia32_sse
nb_kernel213_ia32_sse:  
_nb_kernel213_ia32_sse: 
.set nb213_p_nri, 8
.set nb213_iinr, 12
.set nb213_jindex, 16
.set nb213_jjnr, 20
.set nb213_shift, 24
.set nb213_shiftvec, 28
.set nb213_fshift, 32
.set nb213_gid, 36
.set nb213_pos, 40
.set nb213_faction, 44
.set nb213_charge, 48
.set nb213_p_facel, 52
.set nb213_argkrf, 56
.set nb213_argcrf, 60
.set nb213_Vc, 64
.set nb213_type, 68
.set nb213_p_ntype, 72
.set nb213_vdwparam, 76
.set nb213_Vvdw, 80
.set nb213_p_tabscale, 84
.set nb213_VFtab, 88
.set nb213_invsqrta, 92
.set nb213_dvda, 96
.set nb213_p_gbtabscale, 100
.set nb213_GBtab, 104
.set nb213_p_nthreads, 108
.set nb213_count, 112
.set nb213_mtx, 116
.set nb213_outeriter, 120
.set nb213_inneriter, 124
.set nb213_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb213_ixO, 0
.set nb213_iyO, 16
.set nb213_izO, 32
.set nb213_ixH1, 48
.set nb213_iyH1, 64
.set nb213_izH1, 80
.set nb213_ixH2, 96
.set nb213_iyH2, 112
.set nb213_izH2, 128
.set nb213_ixM, 144
.set nb213_iyM, 160
.set nb213_izM, 176
.set nb213_iqM, 192
.set nb213_iqH, 208
.set nb213_dxO, 224
.set nb213_dyO, 240
.set nb213_dzO, 256
.set nb213_dxH1, 272
.set nb213_dyH1, 288
.set nb213_dzH1, 304
.set nb213_dxH2, 320
.set nb213_dyH2, 336
.set nb213_dzH2, 352
.set nb213_dxM, 368
.set nb213_dyM, 384
.set nb213_dzM, 400
.set nb213_qqM, 416
.set nb213_qqH, 432
.set nb213_rinvH1, 448
.set nb213_rinvH2, 464
.set nb213_rinvM, 480
.set nb213_two, 496
.set nb213_c6, 512
.set nb213_c12, 528
.set nb213_six, 544
.set nb213_twelve, 560
.set nb213_krf, 576
.set nb213_crf, 592
.set nb213_krsqH1, 608
.set nb213_krsqH2, 624
.set nb213_krsqM, 640
.set nb213_vctot, 656
.set nb213_Vvdwtot, 672
.set nb213_fixO, 688
.set nb213_fiyO, 704
.set nb213_fizO, 720
.set nb213_fixH1, 736
.set nb213_fiyH1, 752
.set nb213_fizH1, 768
.set nb213_fixH2, 784
.set nb213_fiyH2, 800
.set nb213_fizH2, 816
.set nb213_fixM, 832
.set nb213_fiyM, 848
.set nb213_fizM, 864
.set nb213_fjx, 880
.set nb213_fjy, 896
.set nb213_fjz, 912
.set nb213_half, 928
.set nb213_three, 944
.set nb213_is3, 960
.set nb213_ii3, 964
.set nb213_ntia, 968
.set nb213_innerjjnr, 972
.set nb213_innerk, 976
.set nb213_n, 980
.set nb213_nn1, 984
.set nb213_nri, 988
.set nb213_nouter, 992
.set nb213_ninner, 996
.set nb213_salign, 1000
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $1004,%esp         ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb213_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb213_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb213_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb213_nouter(%esp)
        movl %eax,nb213_ninner(%esp)

        movl nb213_argkrf(%ebp),%esi
        movl nb213_argcrf(%ebp),%edi
        movss (%esi),%xmm5
        movss (%edi),%xmm6
        shufps $0,%xmm5,%xmm5
        shufps $0,%xmm6,%xmm6
        movaps %xmm5,nb213_krf(%esp)
        movaps %xmm6,nb213_crf(%esp)
        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb213_half(%esp)
        movss nb213_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm3,%xmm4
        addps  %xmm4,%xmm4      ## 6.0
        movaps %xmm4,%xmm5
        addps  %xmm5,%xmm5      ## constant 12.0
        movaps %xmm1,nb213_half(%esp)
        movaps %xmm2,nb213_two(%esp)
        movaps %xmm3,nb213_three(%esp)
        movaps %xmm4,nb213_six(%esp)
        movaps %xmm5,nb213_twelve(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb213_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb213_charge(%ebp),%edx
        movss 4(%edx,%ebx,4),%xmm4
        movss 12(%edx,%ebx,4),%xmm3
        movl nb213_p_facel(%ebp),%esi
        movss (%esi),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb213_iqM(%esp)
        movaps %xmm4,nb213_iqH(%esp)

        movl  nb213_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb213_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb213_ntia(%esp)
_nb_kernel213_ia32_sse.nb213_threadloop: 
        movl  nb213_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel213_ia32_sse.nb213_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel213_ia32_sse.nb213_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb213_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb213_n(%esp)
        movl %ebx,nb213_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel213_ia32_sse.nb213_outerstart
        jmp _nb_kernel213_ia32_sse.nb213_end

_nb_kernel213_ia32_sse.nb213_outerstart: 
        ## ebx contains number of outer iterations
        addl nb213_nouter(%esp),%ebx
        movl %ebx,nb213_nouter(%esp)

_nb_kernel213_ia32_sse.nb213_outer: 
        movl  nb213_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx        ## ebx=3*is 
        movl  %ebx,nb213_is3(%esp)      ## store is3 

        movl  nb213_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb213_iinr(%ebp),%ecx             ## ecx = pointer into iinr[]    
        movl  (%ecx,%esi,4),%ebx                ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        movaps %xmm0,%xmm6
        movaps %xmm1,%xmm7

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb213_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb213_ii3(%esp)

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
        movaps %xmm3,nb213_ixO(%esp)
        movaps %xmm4,nb213_iyO(%esp)
        movaps %xmm5,nb213_izO(%esp)
        movaps %xmm6,nb213_ixH1(%esp)
        movaps %xmm7,nb213_iyH1(%esp)

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
        movaps %xmm6,nb213_izH1(%esp)
        movaps %xmm0,nb213_ixH2(%esp)
        movaps %xmm1,nb213_iyH2(%esp)
        movaps %xmm2,nb213_izH2(%esp)
        movaps %xmm3,nb213_ixM(%esp)
        movaps %xmm4,nb213_iyM(%esp)
        movaps %xmm5,nb213_izM(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb213_vctot(%esp)
        movaps %xmm4,nb213_Vvdwtot(%esp)
        movaps %xmm4,nb213_fixO(%esp)
        movaps %xmm4,nb213_fiyO(%esp)
        movaps %xmm4,nb213_fizO(%esp)
        movaps %xmm4,nb213_fixH1(%esp)
        movaps %xmm4,nb213_fiyH1(%esp)
        movaps %xmm4,nb213_fizH1(%esp)
        movaps %xmm4,nb213_fixH2(%esp)
        movaps %xmm4,nb213_fiyH2(%esp)
        movaps %xmm4,nb213_fizH2(%esp)
        movaps %xmm4,nb213_fixM(%esp)
        movaps %xmm4,nb213_fiyM(%esp)
        movaps %xmm4,nb213_fizM(%esp)

        movl  nb213_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx                ## jindex[n] 
        movl  4(%eax,%esi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movl  nb213_pos(%ebp),%esi
        movl  nb213_faction(%ebp),%edi
        movl  nb213_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb213_innerjjnr(%esp)        ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb213_ninner(%esp),%ecx
        movl  %ecx,nb213_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb213_innerk(%esp)   ## number of innerloop atoms 
        jge   _nb_kernel213_ia32_sse.nb213_unroll_loop
        jmp   _nb_kernel213_ia32_sse.nb213_odd_inner
_nb_kernel213_ia32_sse.nb213_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb213_innerjjnr(%esp),%edx        ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx             ## eax-edx=jnr1-4 

        addl $16,nb213_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb213_charge(%ebp),%esi    ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4              ## and in xmm4 
        mulps  nb213_iqM(%esp),%xmm3
        mulps  nb213_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps  %xmm3,nb213_qqM(%esp)
        movaps  %xmm4,nb213_qqH(%esp)

        movl nb213_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb213_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb213_ntia(%esp),%edi
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

        movaps %xmm4,nb213_c6(%esp)
        movaps %xmm6,nb213_c12(%esp)

        movl nb213_pos(%ebp),%esi       ## base of pos[] 

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
        movaps nb213_ixO(%esp),%xmm4
        movaps nb213_iyO(%esp),%xmm5
        movaps nb213_izO(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb213_dxO(%esp)
        movaps %xmm5,nb213_dyO(%esp)
        movaps %xmm6,nb213_dzO(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        movaps %xmm4,%xmm7
        ## rsqO in xmm7

        ## move ixH1-izH1 to xmm4-xmm6 
        movaps nb213_ixH1(%esp),%xmm4
        movaps nb213_iyH1(%esp),%xmm5
        movaps nb213_izH1(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb213_dxH1(%esp)
        movaps %xmm5,nb213_dyH1(%esp)
        movaps %xmm6,nb213_dzH1(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm6
        addps %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movaps nb213_ixH2(%esp),%xmm3
        movaps nb213_iyH2(%esp),%xmm4
        movaps nb213_izH2(%esp),%xmm5

        ## calc dr 
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        ## store dr 
        movaps %xmm3,nb213_dxH2(%esp)
        movaps %xmm4,nb213_dyH2(%esp)
        movaps %xmm5,nb213_dzH2(%esp)
        ## square it 
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm4,%xmm5
        addps %xmm3,%xmm5

        ## move ixM-izM to xmm2-xmm4  
        movaps nb213_iyM(%esp),%xmm3
        movaps nb213_izM(%esp),%xmm4
        subps  %xmm1,%xmm3
        subps  %xmm2,%xmm4
        movaps nb213_ixM(%esp),%xmm2
        subps  %xmm0,%xmm2

        ## store dr 
        movaps %xmm2,nb213_dxM(%esp)
        movaps %xmm3,nb213_dyM(%esp)
        movaps %xmm4,nb213_dzM(%esp)
        ## square it 
        mulps %xmm2,%xmm2
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        addps %xmm3,%xmm4
        addps %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 
        movaps %xmm4,%xmm0
        movaps %xmm5,%xmm1
        movaps %xmm6,%xmm2
        mulps  nb213_krf(%esp),%xmm0
        mulps  nb213_krf(%esp),%xmm1
        mulps  nb213_krf(%esp),%xmm2
        movaps %xmm0,nb213_krsqM(%esp)
        movaps %xmm1,nb213_krsqH2(%esp)
        movaps %xmm2,nb213_krsqH1(%esp)

        ## rsqH1 - seed in xmm2 
        rsqrtps %xmm6,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb213_three(%esp),%xmm0
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb213_half(%esp),%xmm0
        movaps  %xmm0,nb213_rinvH1(%esp)        ## rinvH1 

        ## rsqH2 - seed to xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb213_three(%esp),%xmm0
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb213_half(%esp),%xmm0
        movaps  %xmm0,nb213_rinvH2(%esp)        ## rinvH2 

        ## rsqM - seed to xmm2 
        rsqrtps %xmm4,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb213_three(%esp),%xmm0
        mulps   %xmm4,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb213_half(%esp),%xmm0
        movaps  %xmm0,nb213_rinvM(%esp)

        ## Do the O LJ-only interaction directly.       
        rcpps   %xmm7,%xmm2
        movaps  nb213_two(%esp),%xmm1
        mulps   %xmm2,%xmm7
        subps   %xmm7,%xmm1
        mulps   %xmm1,%xmm2 ## rinvsq 
        movaps  %xmm2,%xmm0
        mulps   %xmm2,%xmm0     ## r4
        mulps   %xmm2,%xmm0     ## r6
        movaps  %xmm0,%xmm1
        mulps   %xmm1,%xmm1     ## r12
        mulps   nb213_c6(%esp),%xmm0
        mulps   nb213_c12(%esp),%xmm1
        movaps  %xmm1,%xmm3
        subps   %xmm0,%xmm3     ## Vvdw12-Vvdw6
        addps   nb213_Vvdwtot(%esp),%xmm3
        movaps  %xmm3,nb213_Vvdwtot(%esp)
        mulps   nb213_six(%esp),%xmm0
        mulps   nb213_twelve(%esp),%xmm1
        subps   %xmm0,%xmm1
        mulps   %xmm2,%xmm1     ## fscal
        movaps nb213_dxO(%esp),%xmm3
        movaps nb213_dyO(%esp),%xmm4
        movaps nb213_dzO(%esp),%xmm5
        mulps  %xmm1,%xmm3
        mulps  %xmm1,%xmm4
        mulps  %xmm1,%xmm5      ## tx in xmm3-xmm5

        ## update O forces 
        movaps nb213_fixO(%esp),%xmm0
        movaps nb213_fiyO(%esp),%xmm1
        movaps nb213_fizO(%esp),%xmm2
        addps  %xmm3,%xmm0
        addps  %xmm4,%xmm1
        addps  %xmm5,%xmm2
        movaps %xmm0,nb213_fixO(%esp)
        movaps %xmm1,nb213_fiyO(%esp)
        movaps %xmm2,nb213_fizO(%esp)
        ## update j forces with water O 
        movaps %xmm3,nb213_fjx(%esp)
        movaps %xmm4,nb213_fjy(%esp)
        movaps %xmm5,nb213_fjz(%esp)

        ## Do H1 interaction
        movaps  nb213_rinvH1(%esp),%xmm7
        movaps  %xmm7,%xmm4
        mulps   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq
        movaps %xmm7,%xmm0
        movaps nb213_krsqH1(%esp),%xmm1
        addps  %xmm1,%xmm0
        subps  nb213_crf(%esp),%xmm0   ## xmm0=rinv+ krsq-crf 
        mulps  nb213_two(%esp),%xmm1
        subps  %xmm1,%xmm7
        mulps  nb213_qqH(%esp),%xmm0
        mulps  nb213_qqH(%esp),%xmm7

        mulps  %xmm7,%xmm4      ## total fs H1 in xmm4 

        addps  nb213_vctot(%esp),%xmm0
        movaps %xmm0,nb213_vctot(%esp)

        movaps nb213_dxH1(%esp),%xmm0
        movaps nb213_dyH1(%esp),%xmm1
        movaps nb213_dzH1(%esp),%xmm2
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2

        ## update H1 forces 
        movaps nb213_fixH1(%esp),%xmm3
        movaps nb213_fiyH1(%esp),%xmm4
        movaps nb213_fizH1(%esp),%xmm7
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm7
        movaps %xmm3,nb213_fixH1(%esp)
        movaps %xmm4,nb213_fiyH1(%esp)
        movaps %xmm7,nb213_fizH1(%esp)
        ## update j forces with water H1 
        addps  nb213_fjx(%esp),%xmm0
        addps  nb213_fjy(%esp),%xmm1
        addps  nb213_fjz(%esp),%xmm2
        movaps %xmm0,nb213_fjx(%esp)
        movaps %xmm1,nb213_fjy(%esp)
        movaps %xmm2,nb213_fjz(%esp)

        ## Done with H1, do H2 interactions
        movaps  nb213_rinvH2(%esp),%xmm7
        movaps  %xmm7,%xmm4
        mulps   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq
        movaps %xmm7,%xmm0
        movaps nb213_krsqH2(%esp),%xmm1
        addps  %xmm1,%xmm0
        subps  nb213_crf(%esp),%xmm0   ## xmm0=rinv+ krsq-crf 
        mulps  nb213_two(%esp),%xmm1
        subps  %xmm1,%xmm7
        mulps  nb213_qqH(%esp),%xmm0
        mulps  nb213_qqH(%esp),%xmm7

        mulps  %xmm7,%xmm4      ## total fs H2 in xmm4 

        addps  nb213_vctot(%esp),%xmm0
        movaps %xmm0,nb213_vctot(%esp)

        movaps nb213_dxH2(%esp),%xmm0
        movaps nb213_dyH2(%esp),%xmm1
        movaps nb213_dzH2(%esp),%xmm2
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2

        ## update H2 forces 
        movaps nb213_fixH2(%esp),%xmm3
        movaps nb213_fiyH2(%esp),%xmm4
        movaps nb213_fizH2(%esp),%xmm7
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm7
        movaps %xmm3,nb213_fixH2(%esp)
        movaps %xmm4,nb213_fiyH2(%esp)
        movaps %xmm7,nb213_fizH2(%esp)
        addps nb213_fjx(%esp),%xmm0
        addps nb213_fjy(%esp),%xmm1
        addps nb213_fjz(%esp),%xmm2
        movaps %xmm0,nb213_fjx(%esp)
        movaps %xmm1,nb213_fjy(%esp)
        movaps %xmm2,nb213_fjz(%esp)

        ## Done with H2, do M interactions
        movaps  nb213_rinvM(%esp),%xmm7
        movaps  %xmm7,%xmm4
        mulps   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq
        movaps %xmm7,%xmm0
        movaps nb213_krsqM(%esp),%xmm1
        addps  %xmm1,%xmm0
        subps  nb213_crf(%esp),%xmm0   ## xmm0=rinv+ krsq-crf 
        mulps  nb213_two(%esp),%xmm1
        subps  %xmm1,%xmm7
        mulps  nb213_qqM(%esp),%xmm0
        mulps  nb213_qqM(%esp),%xmm7

        mulps  %xmm7,%xmm4      ## total fs M in xmm4 

        addps  nb213_vctot(%esp),%xmm0
        movaps %xmm0,nb213_vctot(%esp)

        movaps nb213_dxM(%esp),%xmm0
        movaps nb213_dyM(%esp),%xmm1
        movaps nb213_dzM(%esp),%xmm2
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2

        ## update M forces 
        movaps nb213_fixM(%esp),%xmm3
        movaps nb213_fiyM(%esp),%xmm4
        movaps nb213_fizM(%esp),%xmm7
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm7
        movaps %xmm3,nb213_fixM(%esp)
        movaps %xmm4,nb213_fiyM(%esp)
        movaps %xmm7,nb213_fizM(%esp)

        movl nb213_faction(%ebp),%edi
        ## update j forces from stored values
        addps nb213_fjx(%esp),%xmm0
        addps nb213_fjy(%esp),%xmm1
        addps nb213_fjz(%esp),%xmm2

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
        subl $4,nb213_innerk(%esp)
        jl    _nb_kernel213_ia32_sse.nb213_odd_inner
        jmp   _nb_kernel213_ia32_sse.nb213_unroll_loop
_nb_kernel213_ia32_sse.nb213_odd_inner: 
        addl $4,nb213_innerk(%esp)
        jnz   _nb_kernel213_ia32_sse.nb213_odd_loop
        jmp   _nb_kernel213_ia32_sse.nb213_updateouterdata
_nb_kernel213_ia32_sse.nb213_odd_loop: 
        movl  nb213_innerjjnr(%esp),%edx        ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb213_innerjjnr(%esp)

        xorps %xmm4,%xmm4       ## clear reg.
        movss nb213_iqM(%esp),%xmm4
        movl nb213_charge(%ebp),%esi
        movhps nb213_iqH(%esp),%xmm4    ## [qM  0  qH  qH] 
        shufps $41,%xmm4,%xmm4 ## [0 qH qH qM]

        movss (%esi,%eax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb213_qqM(%esp)    ## use dummy qq for storage 

        xorps %xmm6,%xmm6
        movl nb213_type(%ebp),%esi
        movl (%esi,%eax,4),%ebx
        movl nb213_vdwparam(%ebp),%esi
        shll %ebx
        addl nb213_ntia(%esp),%ebx
        movlps (%esi,%ebx,4),%xmm6
        movaps %xmm6,%xmm7
        shufps $252,%xmm6,%xmm6 ## constant 11111100
        shufps $253,%xmm7,%xmm7 ## constant 11111101
        movaps %xmm6,nb213_c6(%esp)
        movaps %xmm7,nb213_c12(%esp)

        movl nb213_pos(%ebp),%esi
        leal (%eax,%eax,2),%eax

        movss nb213_ixO(%esp),%xmm3
        movss nb213_iyO(%esp),%xmm4
        movss nb213_izO(%esp),%xmm5
        movss nb213_ixH1(%esp),%xmm0
        movss nb213_iyH1(%esp),%xmm1
        movss nb213_izH1(%esp),%xmm2
        unpcklps nb213_ixH2(%esp),%xmm3         ## ixO ixH2 - -
        unpcklps nb213_iyH2(%esp),%xmm4         ## iyO iyH2 - -
        unpcklps nb213_izH2(%esp),%xmm5         ## izO izH2 - -
        unpcklps nb213_ixM(%esp),%xmm0          ## ixH1 ixM - -
        unpcklps nb213_iyM(%esp),%xmm1          ## iyH1 iyM - -
        unpcklps nb213_izM(%esp),%xmm2          ## izH1 izM - -
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
        movaps %xmm3,nb213_dxO(%esp)
        movaps %xmm4,nb213_dyO(%esp)
        movaps %xmm5,nb213_dzO(%esp)

        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5

        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        ## rsq in xmm4 
        movaps %xmm4,%xmm0
        mulps nb213_krf(%esp),%xmm0
        movaps %xmm0,nb213_krsqM(%esp)

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb213_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb213_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv


        movaps %xmm0,%xmm4
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm5
        mulps  %xmm4,%xmm4      ## xmm1=rinv, xmm4=rinvsq
        movaps nb213_krsqM(%esp),%xmm3
        addps  %xmm3,%xmm5      ## xmm0=rinv+ krsq 
        subps  nb213_crf(%esp),%xmm5   ## xmm0=rinv+ krsq-crf 
        mulps  nb213_two(%esp),%xmm3
        subps  %xmm3,%xmm1      ## xmm1=rinv-2*krsq
        movaps %xmm5,%xmm7
        mulps  nb213_qqM(%esp),%xmm7    ## xmm0=vcoul 
        mulps  nb213_qqM(%esp),%xmm1    ## xmm1=coul part of fs 
        movaps %xmm1,%xmm6

        addps  nb213_vctot(%esp),%xmm7
        movaps %xmm7,nb213_vctot(%esp)

        movaps %xmm0,%xmm1
        mulps  %xmm1,%xmm1
        movaps %xmm1,%xmm2
        mulss  %xmm1,%xmm1
        mulss  %xmm2,%xmm1      ## xmm1=rinvsix
        xorps  %xmm4,%xmm4
        movss  %xmm1,%xmm4
        mulss  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulss  nb213_c6(%esp),%xmm1
        mulss  nb213_c12(%esp),%xmm4
        movaps %xmm4,%xmm3
        subss  %xmm1,%xmm3      ## xmm3=Vvdw12-Vvdw6 
        mulss  nb213_six(%esp),%xmm1
        mulss  nb213_twelve(%esp),%xmm4
        subss  %xmm1,%xmm4
        addss  nb213_Vvdwtot(%esp),%xmm3
        movss  %xmm3,nb213_Vvdwtot(%esp)
        addps  %xmm6,%xmm4
        mulps  %xmm0,%xmm4
        mulps  %xmm0,%xmm4      ## fscal

        movaps nb213_dxO(%esp),%xmm0
        movaps nb213_dyO(%esp),%xmm1
        movaps nb213_dzO(%esp),%xmm2

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2 ## xmm0-xmm2 now contains tx-tz (partial force)

        movss  nb213_fixO(%esp),%xmm3
        movss  nb213_fiyO(%esp),%xmm4
        movss  nb213_fizO(%esp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,nb213_fixO(%esp)
        movss  %xmm4,nb213_fiyO(%esp)
        movss  %xmm5,nb213_fizO(%esp)   ## updated the O force now do the H's

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        shufps $0x39,%xmm3,%xmm3 ## shift right 
        shufps $0x39,%xmm4,%xmm4
        shufps $0x39,%xmm5,%xmm5
        addss  nb213_fixH1(%esp),%xmm3
        addss  nb213_fiyH1(%esp),%xmm4
        addss  nb213_fizH1(%esp),%xmm5
        movss  %xmm3,nb213_fixH1(%esp)
        movss  %xmm4,nb213_fiyH1(%esp)
        movss  %xmm5,nb213_fizH1(%esp)          ## updated the H1 force 

        shufps $0x39,%xmm3,%xmm3
        shufps $0x39,%xmm4,%xmm4
        shufps $0x39,%xmm5,%xmm5
        addss  nb213_fixH2(%esp),%xmm3
        addss  nb213_fiyH2(%esp),%xmm4
        addss  nb213_fizH2(%esp),%xmm5
        movss  %xmm3,nb213_fixH2(%esp)
        movss  %xmm4,nb213_fiyH2(%esp)
        movss  %xmm5,nb213_fizH2(%esp)          ## updated the H2 force 

        movl nb213_faction(%ebp),%edi
        shufps $0x39,%xmm3,%xmm3
        shufps $0x39,%xmm4,%xmm4
        shufps $0x39,%xmm5,%xmm5
        addss  nb213_fixM(%esp),%xmm3
        addss  nb213_fiyM(%esp),%xmm4
        addss  nb213_fizM(%esp),%xmm5
        movss  %xmm3,nb213_fixM(%esp)
        movss  %xmm4,nb213_fiyM(%esp)
        movss  %xmm5,nb213_fizM(%esp)   ## updated the M force 

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

        decl nb213_innerk(%esp)
        jz    _nb_kernel213_ia32_sse.nb213_updateouterdata
        jmp   _nb_kernel213_ia32_sse.nb213_odd_loop
_nb_kernel213_ia32_sse.nb213_updateouterdata: 
        movl  nb213_ii3(%esp),%ecx
        movl  nb213_faction(%ebp),%edi
        movl  nb213_fshift(%ebp),%esi
        movl  nb213_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movaps nb213_fixO(%esp),%xmm0
        movaps nb213_fiyO(%esp),%xmm1
        movaps nb213_fizO(%esp),%xmm2

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
        movaps nb213_fixH1(%esp),%xmm0
        movaps nb213_fiyH1(%esp),%xmm1
        movaps nb213_fizH1(%esp),%xmm2

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
        movaps nb213_fixH2(%esp),%xmm0
        movaps nb213_fiyH2(%esp),%xmm1
        movaps nb213_fizH2(%esp),%xmm2

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
        movaps nb213_fixM(%esp),%xmm0
        movaps nb213_fiyM(%esp),%xmm1
        movaps nb213_fizM(%esp),%xmm2

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
        movl nb213_n(%esp),%esi
        ## get group index for i particle 
        movl  nb213_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb213_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb213_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb213_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb213_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb213_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel213_ia32_sse.nb213_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb213_n(%esp)
        jmp _nb_kernel213_ia32_sse.nb213_outer
_nb_kernel213_ia32_sse.nb213_outerend: 
        ## check if more outer neighborlists remain
        movl  nb213_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel213_ia32_sse.nb213_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel213_ia32_sse.nb213_threadloop
_nb_kernel213_ia32_sse.nb213_end: 
        emms

        movl nb213_nouter(%esp),%eax
        movl nb213_ninner(%esp),%ebx
        movl nb213_outeriter(%ebp),%ecx
        movl nb213_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb213_salign(%esp),%eax
        addl %eax,%esp
        addl $1004,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret




.globl nb_kernel213nf_ia32_sse
.globl _nb_kernel213nf_ia32_sse
nb_kernel213nf_ia32_sse:        
_nb_kernel213nf_ia32_sse:       
.set nb213nf_p_nri, 8
.set nb213nf_iinr, 12
.set nb213nf_jindex, 16
.set nb213nf_jjnr, 20
.set nb213nf_shift, 24
.set nb213nf_shiftvec, 28
.set nb213nf_fshift, 32
.set nb213nf_gid, 36
.set nb213nf_pos, 40
.set nb213nf_faction, 44
.set nb213nf_charge, 48
.set nb213nf_p_facel, 52
.set nb213nf_argkrf, 56
.set nb213nf_argcrf, 60
.set nb213nf_Vc, 64
.set nb213nf_type, 68
.set nb213nf_p_ntype, 72
.set nb213nf_vdwparam, 76
.set nb213nf_Vvdw, 80
.set nb213nf_p_tabscale, 84
.set nb213nf_VFtab, 88
.set nb213nf_invsqrta, 92
.set nb213nf_dvda, 96
.set nb213nf_p_gbtabscale, 100
.set nb213nf_GBtab, 104
.set nb213nf_p_nthreads, 108
.set nb213nf_count, 112
.set nb213nf_mtx, 116
.set nb213nf_outeriter, 120
.set nb213nf_inneriter, 124
.set nb213nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb213nf_ixO, 0
.set nb213nf_iyO, 16
.set nb213nf_izO, 32
.set nb213nf_ixH1, 48
.set nb213nf_iyH1, 64
.set nb213nf_izH1, 80
.set nb213nf_ixH2, 96
.set nb213nf_iyH2, 112
.set nb213nf_izH2, 128
.set nb213nf_ixM, 144
.set nb213nf_iyM, 160
.set nb213nf_izM, 176
.set nb213nf_iqM, 192
.set nb213nf_iqH, 208
.set nb213nf_qqM, 224
.set nb213nf_qqH, 240
.set nb213nf_rinvH1, 256
.set nb213nf_rinvH2, 272
.set nb213nf_rinvM, 288
.set nb213nf_two, 304
.set nb213nf_c6, 320
.set nb213nf_c12, 336
.set nb213nf_krf, 352
.set nb213nf_crf, 368
.set nb213nf_krsqH1, 384
.set nb213nf_krsqH2, 400
.set nb213nf_krsqM, 416
.set nb213nf_vctot, 432
.set nb213nf_Vvdwtot, 448
.set nb213nf_half, 464
.set nb213nf_three, 480
.set nb213nf_is3, 496
.set nb213nf_ii3, 500
.set nb213nf_ntia, 504
.set nb213nf_innerjjnr, 508
.set nb213nf_innerk, 512
.set nb213nf_n, 516
.set nb213nf_nn1, 520
.set nb213nf_nri, 524
.set nb213nf_nouter, 528
.set nb213nf_ninner, 532
.set nb213nf_salign, 536
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
        movl %eax,nb213nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb213nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb213nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb213nf_nouter(%esp)
        movl %eax,nb213nf_ninner(%esp)

        movl nb213nf_argkrf(%ebp),%esi
        movl nb213nf_argcrf(%ebp),%edi
        movss (%esi),%xmm5
        movss (%edi),%xmm6
        shufps $0,%xmm5,%xmm5
        shufps $0,%xmm6,%xmm6
        movaps %xmm5,nb213nf_krf(%esp)
        movaps %xmm6,nb213nf_crf(%esp)
        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb213nf_half(%esp)
        movss nb213nf_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb213nf_half(%esp)
        movaps %xmm2,nb213nf_two(%esp)
        movaps %xmm3,nb213nf_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb213nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb213nf_charge(%ebp),%edx
        movss 4(%edx,%ebx,4),%xmm4
        movss 12(%edx,%ebx,4),%xmm3
        movl nb213nf_p_facel(%ebp),%esi
        movss (%esi),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb213nf_iqM(%esp)
        movaps %xmm4,nb213nf_iqH(%esp)

        movl  nb213nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb213nf_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb213nf_ntia(%esp)
_nb_kernel213nf_ia32_sse.nb213nf_threadloop: 
        movl  nb213nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel213nf_ia32_sse.nb213nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel213nf_ia32_sse.nb213nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb213nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb213nf_n(%esp)
        movl %ebx,nb213nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel213nf_ia32_sse.nb213nf_outerstart
        jmp _nb_kernel213nf_ia32_sse.nb213nf_end

_nb_kernel213nf_ia32_sse.nb213nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb213nf_nouter(%esp),%ebx
        movl %ebx,nb213nf_nouter(%esp)

_nb_kernel213nf_ia32_sse.nb213nf_outer: 
        movl  nb213nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx        ## ebx=3*is 
        movl  %ebx,nb213nf_is3(%esp)            ## store is3 

        movl  nb213nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb213nf_iinr(%ebp),%ecx           ## ecx = pointer into iinr[]    
        movl  (%ecx,%esi,4),%ebx                ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        movaps %xmm0,%xmm6
        movaps %xmm1,%xmm7

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb213nf_pos(%ebp),%eax    ## eax = base of pos[]  
        movl  %ebx,nb213nf_ii3(%esp)

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
        movaps %xmm3,nb213nf_ixO(%esp)
        movaps %xmm4,nb213nf_iyO(%esp)
        movaps %xmm5,nb213nf_izO(%esp)
        movaps %xmm6,nb213nf_ixH1(%esp)
        movaps %xmm7,nb213nf_iyH1(%esp)

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
        movaps %xmm6,nb213nf_izH1(%esp)
        movaps %xmm0,nb213nf_ixH2(%esp)
        movaps %xmm1,nb213nf_iyH2(%esp)
        movaps %xmm2,nb213nf_izH2(%esp)
        movaps %xmm3,nb213nf_ixM(%esp)
        movaps %xmm4,nb213nf_iyM(%esp)
        movaps %xmm5,nb213nf_izM(%esp)

        ## clear vctot 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb213nf_vctot(%esp)
        movaps %xmm4,nb213nf_Vvdwtot(%esp)

        movl  nb213nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx                ## jindex[n] 
        movl  4(%eax,%esi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movl  nb213nf_pos(%ebp),%esi
        movl  nb213nf_faction(%ebp),%edi
        movl  nb213nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb213nf_innerjjnr(%esp)      ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb213nf_ninner(%esp),%ecx
        movl  %ecx,nb213nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb213nf_innerk(%esp)         ## number of innerloop atoms 
        jge   _nb_kernel213nf_ia32_sse.nb213nf_unroll_loop
        jmp   _nb_kernel213nf_ia32_sse.nb213nf_odd_inner
_nb_kernel213nf_ia32_sse.nb213nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb213nf_innerjjnr(%esp),%edx      ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx             ## eax-edx=jnr1-4 

        addl $16,nb213nf_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb213nf_charge(%ebp),%esi  ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4              ## and in xmm4 
        mulps  nb213nf_iqM(%esp),%xmm3
        mulps  nb213nf_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps  %xmm3,nb213nf_qqM(%esp)
        movaps  %xmm4,nb213nf_qqH(%esp)

        movl nb213nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb213nf_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb213nf_ntia(%esp),%edi
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

        movaps %xmm4,nb213nf_c6(%esp)
        movaps %xmm6,nb213nf_c12(%esp)

        movl nb213nf_pos(%ebp),%esi     ## base of pos[] 

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
        movaps nb213nf_ixO(%esp),%xmm4
        movaps nb213nf_iyO(%esp),%xmm5
        movaps nb213nf_izO(%esp),%xmm6

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
        movaps nb213nf_ixH1(%esp),%xmm4
        movaps nb213nf_iyH1(%esp),%xmm5
        movaps nb213nf_izH1(%esp),%xmm6

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
        movaps nb213nf_ixH2(%esp),%xmm3
        movaps nb213nf_iyH2(%esp),%xmm4
        movaps nb213nf_izH2(%esp),%xmm5

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
        movaps nb213nf_iyM(%esp),%xmm3
        movaps nb213nf_izM(%esp),%xmm4
        subps  %xmm1,%xmm3
        subps  %xmm2,%xmm4
        movaps nb213nf_ixM(%esp),%xmm2
        subps  %xmm0,%xmm2

        ## square it 
        mulps %xmm2,%xmm2
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        addps %xmm3,%xmm4
        addps %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 
        movaps %xmm4,%xmm0
        movaps %xmm5,%xmm1
        movaps %xmm6,%xmm2
        mulps  nb213nf_krf(%esp),%xmm0
        mulps  nb213nf_krf(%esp),%xmm1
        mulps  nb213nf_krf(%esp),%xmm2
        movaps %xmm0,nb213nf_krsqM(%esp)
        movaps %xmm1,nb213nf_krsqH2(%esp)
        movaps %xmm2,nb213nf_krsqH1(%esp)

        ## rsqH1 - seed in xmm2 
        rsqrtps %xmm6,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb213nf_three(%esp),%xmm0
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb213nf_half(%esp),%xmm0
        movaps  %xmm0,nb213nf_rinvH1(%esp)      ## rinvH1 

        ## rsqH2 - seed to xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb213nf_three(%esp),%xmm0
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb213nf_half(%esp),%xmm0
        movaps  %xmm0,nb213nf_rinvH2(%esp)      ## rinvH2 

        ## rsqM - seed to xmm2 
        rsqrtps %xmm4,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb213nf_three(%esp),%xmm0
        mulps   %xmm4,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb213nf_half(%esp),%xmm0
        movaps  %xmm0,nb213nf_rinvM(%esp)

        ## Do the O LJ-only interaction directly.       
        rcpps   %xmm7,%xmm2
        movaps  nb213nf_two(%esp),%xmm1
        mulps   %xmm2,%xmm7
        subps   %xmm7,%xmm1
        mulps   %xmm1,%xmm2 ## rinvsq 
        movaps  %xmm2,%xmm0
        mulps   %xmm2,%xmm0     ## r4
        mulps   %xmm2,%xmm0     ## r6
        movaps  %xmm0,%xmm1
        mulps   %xmm1,%xmm1     ## r12
        mulps   nb213nf_c6(%esp),%xmm0
        mulps   nb213nf_c12(%esp),%xmm1
        movaps  %xmm1,%xmm3
        subps   %xmm0,%xmm3     ## Vvdw12-Vvdw6
        addps   nb213nf_Vvdwtot(%esp),%xmm3
        movaps  %xmm3,nb213nf_Vvdwtot(%esp)

        ## do H1 interactions
        movaps  nb213nf_rinvH1(%esp),%xmm7
        addps nb213nf_krsqH1(%esp),%xmm7
        subps nb213nf_crf(%esp),%xmm7   ## xmm7=rinv+ krsq-crf 
        mulps nb213nf_qqH(%esp),%xmm7
        addps nb213nf_vctot(%esp),%xmm7

        ## H2 interactions 
        movaps  nb213nf_rinvH2(%esp),%xmm6
        addps nb213nf_krsqH2(%esp),%xmm6
        subps nb213nf_crf(%esp),%xmm6   ## xmm6=rinv+ krsq-crf 
        mulps nb213nf_qqH(%esp),%xmm6
        addps %xmm7,%xmm6

        ## M interactions 
        movaps nb213nf_rinvM(%esp),%xmm5
        addps nb213nf_krsqM(%esp),%xmm5
        subps nb213nf_crf(%esp),%xmm5   ## xmm5=rinv+ krsq-crf 
        mulps nb213nf_qqM(%esp),%xmm5
        addps %xmm6,%xmm5
        movaps %xmm5,nb213nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $4,nb213nf_innerk(%esp)
        jl    _nb_kernel213nf_ia32_sse.nb213nf_odd_inner
        jmp   _nb_kernel213nf_ia32_sse.nb213nf_unroll_loop
_nb_kernel213nf_ia32_sse.nb213nf_odd_inner: 
        addl $4,nb213nf_innerk(%esp)
        jnz   _nb_kernel213nf_ia32_sse.nb213nf_odd_loop
        jmp   _nb_kernel213nf_ia32_sse.nb213nf_updateouterdata
_nb_kernel213nf_ia32_sse.nb213nf_odd_loop: 
        movl  nb213nf_innerjjnr(%esp),%edx      ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb213nf_innerjjnr(%esp)

        xorps %xmm4,%xmm4       ## clear reg.
        movss nb213nf_iqM(%esp),%xmm4
        movl nb213nf_charge(%ebp),%esi
        movhps nb213nf_iqH(%esp),%xmm4    ## [qM  0  qH  qH] 
        shufps $41,%xmm4,%xmm4 ## [0 qH qH qM]

        movss (%esi,%eax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb213nf_qqM(%esp)          ## use dummy qq for storage 

        xorps %xmm6,%xmm6
        movl nb213nf_type(%ebp),%esi
        movl (%esi,%eax,4),%ebx
        movl nb213nf_vdwparam(%ebp),%esi
        shll %ebx
        addl nb213nf_ntia(%esp),%ebx
        movlps (%esi,%ebx,4),%xmm6
        movaps %xmm6,%xmm7
        shufps $252,%xmm6,%xmm6 ## constant 11111100
        shufps $253,%xmm7,%xmm7 ## constant 11111101
        movaps %xmm6,nb213nf_c6(%esp)
        movaps %xmm7,nb213nf_c12(%esp)

        movl nb213nf_pos(%ebp),%esi
        leal (%eax,%eax,2),%eax

        movss nb213nf_ixO(%esp),%xmm3
        movss nb213nf_iyO(%esp),%xmm4
        movss nb213nf_izO(%esp),%xmm5
        movss nb213nf_ixH1(%esp),%xmm0
        movss nb213nf_iyH1(%esp),%xmm1
        movss nb213nf_izH1(%esp),%xmm2
        unpcklps nb213nf_ixH2(%esp),%xmm3       ## ixO ixH2 - -
        unpcklps nb213nf_iyH2(%esp),%xmm4       ## iyO iyH2 - -
        unpcklps nb213nf_izH2(%esp),%xmm5       ## izO izH2 - -
        unpcklps nb213nf_ixM(%esp),%xmm0        ## ixH1 ixM - -
        unpcklps nb213nf_iyM(%esp),%xmm1        ## iyH1 iyM - -
        unpcklps nb213nf_izM(%esp),%xmm2        ## izH1 izM - -
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
        movaps %xmm4,%xmm0
        mulps nb213nf_krf(%esp),%xmm0
        movaps %xmm0,nb213nf_krsqM(%esp)

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb213nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb213nf_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv

        movaps %xmm0,%xmm1      ## xmm1=rinv
        addps  nb213nf_krsqM(%esp),%xmm1
        subps  nb213nf_crf(%esp),%xmm1   ## xmm0=rinv+ krsq-crf 
        mulps  nb213nf_qqM(%esp),%xmm1
        addps  nb213nf_vctot(%esp),%xmm1
        movaps %xmm1,nb213nf_vctot(%esp)

        movaps %xmm0,%xmm1
        mulps  %xmm1,%xmm1
        movaps %xmm1,%xmm2
        mulss  %xmm1,%xmm1
        mulss  %xmm2,%xmm1      ## xmm1=rinvsix
        xorps  %xmm4,%xmm4
        movss  %xmm1,%xmm4
        mulss  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulss  nb213nf_c6(%esp),%xmm1
        mulss  nb213nf_c12(%esp),%xmm4
        movaps %xmm4,%xmm3
        subss  %xmm1,%xmm3      ## xmm3=Vvdw12-Vvdw6 
        addss  nb213nf_Vvdwtot(%esp),%xmm3
        movss  %xmm3,nb213nf_Vvdwtot(%esp)

        decl nb213nf_innerk(%esp)
        jz    _nb_kernel213nf_ia32_sse.nb213nf_updateouterdata
        jmp   _nb_kernel213nf_ia32_sse.nb213nf_odd_loop
_nb_kernel213nf_ia32_sse.nb213nf_updateouterdata: 
        ## get n from stack
        movl nb213nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb213nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb213nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb213nf_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb213nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb213nf_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb213nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel213nf_ia32_sse.nb213nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb213nf_n(%esp)
        jmp _nb_kernel213nf_ia32_sse.nb213nf_outer
_nb_kernel213nf_ia32_sse.nb213nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb213nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel213nf_ia32_sse.nb213nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel213nf_ia32_sse.nb213nf_threadloop
_nb_kernel213nf_ia32_sse.nb213nf_end: 
        emms

        movl nb213nf_nouter(%esp),%eax
        movl nb213nf_ninner(%esp),%ebx
        movl nb213nf_outeriter(%ebp),%ecx
        movl nb213nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb213nf_salign(%esp),%eax
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


