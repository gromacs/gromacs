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


.globl nb_kernel233_ia32_sse
.globl _nb_kernel233_ia32_sse
nb_kernel233_ia32_sse:  
_nb_kernel233_ia32_sse: 
.set nb233_p_nri, 8
.set nb233_iinr, 12
.set nb233_jindex, 16
.set nb233_jjnr, 20
.set nb233_shift, 24
.set nb233_shiftvec, 28
.set nb233_fshift, 32
.set nb233_gid, 36
.set nb233_pos, 40
.set nb233_faction, 44
.set nb233_charge, 48
.set nb233_p_facel, 52
.set nb233_argkrf, 56
.set nb233_argcrf, 60
.set nb233_Vc, 64
.set nb233_type, 68
.set nb233_p_ntype, 72
.set nb233_vdwparam, 76
.set nb233_Vvdw, 80
.set nb233_p_tabscale, 84
.set nb233_VFtab, 88
.set nb233_invsqrta, 92
.set nb233_dvda, 96
.set nb233_p_gbtabscale, 100
.set nb233_GBtab, 104
.set nb233_p_nthreads, 108
.set nb233_count, 112
.set nb233_mtx, 116
.set nb233_outeriter, 120
.set nb233_inneriter, 124
.set nb233_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb233_ixO, 0
.set nb233_iyO, 16
.set nb233_izO, 32
.set nb233_ixH1, 48
.set nb233_iyH1, 64
.set nb233_izH1, 80
.set nb233_ixH2, 96
.set nb233_iyH2, 112
.set nb233_izH2, 128
.set nb233_ixM, 144
.set nb233_iyM, 160
.set nb233_izM, 176
.set nb233_iqM, 192
.set nb233_iqH, 208
.set nb233_dxO, 224
.set nb233_dyO, 240
.set nb233_dzO, 256
.set nb233_dxH1, 272
.set nb233_dyH1, 288
.set nb233_dzH1, 304
.set nb233_dxH2, 320
.set nb233_dyH2, 336
.set nb233_dzH2, 352
.set nb233_dxM, 368
.set nb233_dyM, 384
.set nb233_dzM, 400
.set nb233_qqM, 416
.set nb233_qqH, 432
.set nb233_rinvH1, 448
.set nb233_rinvH2, 464
.set nb233_rinvM, 480
.set nb233_two, 496
.set nb233_c6, 512
.set nb233_c12, 528
.set nb233_tsc, 544
.set nb233_fstmp, 560
.set nb233_krf, 576
.set nb233_crf, 592
.set nb233_krsqH1, 608
.set nb233_krsqH2, 624
.set nb233_krsqM, 640
.set nb233_vctot, 656
.set nb233_Vvdwtot, 672
.set nb233_fixO, 688
.set nb233_fiyO, 704
.set nb233_fizO, 720
.set nb233_fixH1, 736
.set nb233_fiyH1, 752
.set nb233_fizH1, 768
.set nb233_fixH2, 784
.set nb233_fiyH2, 800
.set nb233_fizH2, 816
.set nb233_fixM, 832
.set nb233_fiyM, 848
.set nb233_fizM, 864
.set nb233_fjx, 880
.set nb233_fjy, 896
.set nb233_fjz, 912
.set nb233_half, 928
.set nb233_three, 944
.set nb233_is3, 960
.set nb233_ii3, 964
.set nb233_ntia, 968
.set nb233_innerjjnr, 972
.set nb233_innerk, 976
.set nb233_n, 980
.set nb233_nn1, 984
.set nb233_nri, 988
.set nb233_nouter, 992
.set nb233_ninner, 996
.set nb233_salign, 1000
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
        movl %eax,nb233_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb233_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb233_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb233_nouter(%esp)
        movl %eax,nb233_ninner(%esp)

        movl nb233_p_tabscale(%ebp),%eax
        movss (%eax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb233_tsc(%esp)

        movl nb233_argkrf(%ebp),%esi
        movl nb233_argcrf(%ebp),%edi
        movss (%esi),%xmm5
        movss (%edi),%xmm6
        shufps $0,%xmm5,%xmm5
        shufps $0,%xmm6,%xmm6
        movaps %xmm5,nb233_krf(%esp)
        movaps %xmm6,nb233_crf(%esp)
        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb233_half(%esp)
        movss nb233_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb233_half(%esp)
        movaps %xmm2,nb233_two(%esp)
        movaps %xmm3,nb233_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb233_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb233_charge(%ebp),%edx
        movss 4(%edx,%ebx,4),%xmm4
        movss 12(%edx,%ebx,4),%xmm3
        movl nb233_p_facel(%ebp),%esi
        movss (%esi),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb233_iqM(%esp)
        movaps %xmm4,nb233_iqH(%esp)

        movl  nb233_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb233_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb233_ntia(%esp)
_nb_kernel233_ia32_sse.nb233_threadloop: 
        movl  nb233_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel233_ia32_sse.nb233_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel233_ia32_sse.nb233_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb233_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb233_n(%esp)
        movl %ebx,nb233_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel233_ia32_sse.nb233_outerstart
        jmp _nb_kernel233_ia32_sse.nb233_end

_nb_kernel233_ia32_sse.nb233_outerstart: 
        ## ebx contains number of outer iterations
        addl nb233_nouter(%esp),%ebx
        movl %ebx,nb233_nouter(%esp)

_nb_kernel233_ia32_sse.nb233_outer: 
        movl  nb233_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx        ## ebx=3*is 
        movl  %ebx,nb233_is3(%esp)      ## store is3 

        movl  nb233_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb233_iinr(%ebp),%ecx             ## ecx = pointer into iinr[]    
        movl  (%ecx,%esi,4),%ebx                ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        movaps %xmm0,%xmm6
        movaps %xmm1,%xmm7

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb233_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb233_ii3(%esp)

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
        movaps %xmm3,nb233_ixO(%esp)
        movaps %xmm4,nb233_iyO(%esp)
        movaps %xmm5,nb233_izO(%esp)
        movaps %xmm6,nb233_ixH1(%esp)
        movaps %xmm7,nb233_iyH1(%esp)

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
        movaps %xmm6,nb233_izH1(%esp)
        movaps %xmm0,nb233_ixH2(%esp)
        movaps %xmm1,nb233_iyH2(%esp)
        movaps %xmm2,nb233_izH2(%esp)
        movaps %xmm3,nb233_ixM(%esp)
        movaps %xmm4,nb233_iyM(%esp)
        movaps %xmm5,nb233_izM(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb233_vctot(%esp)
        movaps %xmm4,nb233_Vvdwtot(%esp)
        movaps %xmm4,nb233_fixO(%esp)
        movaps %xmm4,nb233_fiyO(%esp)
        movaps %xmm4,nb233_fizO(%esp)
        movaps %xmm4,nb233_fixH1(%esp)
        movaps %xmm4,nb233_fiyH1(%esp)
        movaps %xmm4,nb233_fizH1(%esp)
        movaps %xmm4,nb233_fixH2(%esp)
        movaps %xmm4,nb233_fiyH2(%esp)
        movaps %xmm4,nb233_fizH2(%esp)
        movaps %xmm4,nb233_fixM(%esp)
        movaps %xmm4,nb233_fiyM(%esp)
        movaps %xmm4,nb233_fizM(%esp)

        movl  nb233_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx                ## jindex[n] 
        movl  4(%eax,%esi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movl  nb233_pos(%ebp),%esi
        movl  nb233_faction(%ebp),%edi
        movl  nb233_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb233_innerjjnr(%esp)        ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb233_ninner(%esp),%ecx
        movl  %ecx,nb233_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb233_innerk(%esp)   ## number of innerloop atoms 
        jge   _nb_kernel233_ia32_sse.nb233_unroll_loop
        jmp   _nb_kernel233_ia32_sse.nb233_odd_inner
_nb_kernel233_ia32_sse.nb233_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb233_innerjjnr(%esp),%edx        ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx             ## eax-edx=jnr1-4 

        addl $16,nb233_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb233_charge(%ebp),%esi    ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4              ## and in xmm4 
        mulps  nb233_iqM(%esp),%xmm3
        mulps  nb233_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps  %xmm3,nb233_qqM(%esp)
        movaps  %xmm4,nb233_qqH(%esp)

        movl nb233_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb233_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb233_ntia(%esp),%edi
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

        movaps %xmm4,nb233_c6(%esp)
        movaps %xmm6,nb233_c12(%esp)

        movl nb233_pos(%ebp),%esi       ## base of pos[] 

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
        movaps nb233_ixO(%esp),%xmm4
        movaps nb233_iyO(%esp),%xmm5
        movaps nb233_izO(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb233_dxO(%esp)
        movaps %xmm5,nb233_dyO(%esp)
        movaps %xmm6,nb233_dzO(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        movaps %xmm4,%xmm7
        ## rsqO in xmm7

        ## move ixH1-izH1 to xmm4-xmm6 
        movaps nb233_ixH1(%esp),%xmm4
        movaps nb233_iyH1(%esp),%xmm5
        movaps nb233_izH1(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb233_dxH1(%esp)
        movaps %xmm5,nb233_dyH1(%esp)
        movaps %xmm6,nb233_dzH1(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm6
        addps %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movaps nb233_ixH2(%esp),%xmm3
        movaps nb233_iyH2(%esp),%xmm4
        movaps nb233_izH2(%esp),%xmm5

        ## calc dr 
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        ## store dr 
        movaps %xmm3,nb233_dxH2(%esp)
        movaps %xmm4,nb233_dyH2(%esp)
        movaps %xmm5,nb233_dzH2(%esp)
        ## square it 
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm4,%xmm5
        addps %xmm3,%xmm5

        ## move ixM-izM to xmm2-xmm4  
        movaps nb233_iyM(%esp),%xmm3
        movaps nb233_izM(%esp),%xmm4
        subps  %xmm1,%xmm3
        subps  %xmm2,%xmm4
        movaps nb233_ixM(%esp),%xmm2
        subps  %xmm0,%xmm2

        ## store dr 
        movaps %xmm2,nb233_dxM(%esp)
        movaps %xmm3,nb233_dyM(%esp)
        movaps %xmm4,nb233_dzM(%esp)
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
        mulps  nb233_krf(%esp),%xmm0
        mulps  nb233_krf(%esp),%xmm1
        mulps  nb233_krf(%esp),%xmm2
        movaps %xmm0,nb233_krsqM(%esp)
        movaps %xmm1,nb233_krsqH2(%esp)
        movaps %xmm2,nb233_krsqH1(%esp)

        ## rsqH1 - seed in xmm2 
        rsqrtps %xmm6,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb233_three(%esp),%xmm0
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb233_half(%esp),%xmm0
        movaps  %xmm0,nb233_rinvH1(%esp)        ## rinvH1 

        ## rsqH2 - seed to xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb233_three(%esp),%xmm0
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb233_half(%esp),%xmm0
        movaps  %xmm0,nb233_rinvH2(%esp)        ## rinvH2 

        ## rsqM - seed to xmm2 
        rsqrtps %xmm4,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb233_three(%esp),%xmm0
        mulps   %xmm4,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb233_half(%esp),%xmm0
        movaps  %xmm0,nb233_rinvM(%esp)

        ## Do the O LJ-only interaction directly.       
        ## rsqO is in xmm7
        rsqrtps %xmm7,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb233_three(%esp),%xmm4
        mulps   %xmm7,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb233_half(%esp),%xmm4
        movaps  %xmm4,%xmm0
        ## xmm0=rinvO

        mulps %xmm0,%xmm7
        mulps nb233_tsc(%esp),%xmm7   ## rtab

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

        movl nb233_VFtab(%ebp),%esi
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
        mulps  nb233_two(%esp),%xmm7    ## two*Heps2 
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb233_c6(%esp),%xmm4
        mulps  %xmm4,%xmm7       ## fijD 
        mulps  %xmm4,%xmm5       ## Vvdw6 
        mulps  nb233_tsc(%esp),%xmm7
        ## put scalar force on stack Update Vvdwtot directly 
        addps  nb233_Vvdwtot(%esp),%xmm5
        movaps %xmm7,nb233_fstmp(%esp)
        movaps %xmm5,nb233_Vvdwtot(%esp)

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
        mulps  nb233_two(%esp),%xmm7    ## two*Heps2 
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb233_c12(%esp),%xmm4
        mulps  %xmm4,%xmm7 ## fijR 
        mulps  %xmm4,%xmm5 ## Vvdw12 
        mulps  nb233_tsc(%esp),%xmm7
        addps  nb233_fstmp(%esp),%xmm7

        addps  nb233_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb233_Vvdwtot(%esp)

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movd  %mm2,%ecx
        movd  %mm3,%edx

        xorps  %xmm1,%xmm1
        mulps  %xmm0,%xmm7
        subps  %xmm7,%xmm1

        movaps nb233_dxO(%esp),%xmm3
        movaps nb233_dyO(%esp),%xmm4
        movaps nb233_dzO(%esp),%xmm5
        mulps  %xmm1,%xmm3
        mulps  %xmm1,%xmm4
        mulps  %xmm1,%xmm5      ## tx in xmm3-xmm5

        ## update O forces 
        movaps nb233_fixO(%esp),%xmm0
        movaps nb233_fiyO(%esp),%xmm1
        movaps nb233_fizO(%esp),%xmm2
        addps  %xmm3,%xmm0
        addps  %xmm4,%xmm1
        addps  %xmm5,%xmm2
        movaps %xmm0,nb233_fixO(%esp)
        movaps %xmm1,nb233_fiyO(%esp)
        movaps %xmm2,nb233_fizO(%esp)
        ## update j forces with water O 
        movaps %xmm3,nb233_fjx(%esp)
        movaps %xmm4,nb233_fjy(%esp)
        movaps %xmm5,nb233_fjz(%esp)

        ## Do H1 interaction
        movaps  nb233_rinvH1(%esp),%xmm7
        movaps  %xmm7,%xmm4
        mulps   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq
        movaps %xmm7,%xmm0
        movaps nb233_krsqH1(%esp),%xmm1
        addps  %xmm1,%xmm0
        subps  nb233_crf(%esp),%xmm0   ## xmm0=rinv+ krsq-crf 
        mulps  nb233_two(%esp),%xmm1
        subps  %xmm1,%xmm7
        mulps  nb233_qqH(%esp),%xmm0
        mulps  nb233_qqH(%esp),%xmm7

        mulps  %xmm7,%xmm4      ## total fs H1 in xmm4 

        addps  nb233_vctot(%esp),%xmm0
        movaps %xmm0,nb233_vctot(%esp)

        movaps nb233_dxH1(%esp),%xmm0
        movaps nb233_dyH1(%esp),%xmm1
        movaps nb233_dzH1(%esp),%xmm2
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2

        ## update H1 forces 
        movaps nb233_fixH1(%esp),%xmm3
        movaps nb233_fiyH1(%esp),%xmm4
        movaps nb233_fizH1(%esp),%xmm7
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm7
        movaps %xmm3,nb233_fixH1(%esp)
        movaps %xmm4,nb233_fiyH1(%esp)
        movaps %xmm7,nb233_fizH1(%esp)
        ## update j forces with water H1 
        addps  nb233_fjx(%esp),%xmm0
        addps  nb233_fjy(%esp),%xmm1
        addps  nb233_fjz(%esp),%xmm2
        movaps %xmm0,nb233_fjx(%esp)
        movaps %xmm1,nb233_fjy(%esp)
        movaps %xmm2,nb233_fjz(%esp)

        ## Done with H1, do H2 interactions
        movaps  nb233_rinvH2(%esp),%xmm7
        movaps  %xmm7,%xmm4
        mulps   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq
        movaps %xmm7,%xmm0
        movaps nb233_krsqH2(%esp),%xmm1
        addps  %xmm1,%xmm0
        subps  nb233_crf(%esp),%xmm0   ## xmm0=rinv+ krsq-crf 
        mulps  nb233_two(%esp),%xmm1
        subps  %xmm1,%xmm7
        mulps  nb233_qqH(%esp),%xmm0
        mulps  nb233_qqH(%esp),%xmm7

        mulps  %xmm7,%xmm4      ## total fs H2 in xmm4 

        addps  nb233_vctot(%esp),%xmm0
        movaps %xmm0,nb233_vctot(%esp)

        movaps nb233_dxH2(%esp),%xmm0
        movaps nb233_dyH2(%esp),%xmm1
        movaps nb233_dzH2(%esp),%xmm2
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2

        ## update H2 forces 
        movaps nb233_fixH2(%esp),%xmm3
        movaps nb233_fiyH2(%esp),%xmm4
        movaps nb233_fizH2(%esp),%xmm7
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm7
        movaps %xmm3,nb233_fixH2(%esp)
        movaps %xmm4,nb233_fiyH2(%esp)
        movaps %xmm7,nb233_fizH2(%esp)
        addps nb233_fjx(%esp),%xmm0
        addps nb233_fjy(%esp),%xmm1
        addps nb233_fjz(%esp),%xmm2
        movaps %xmm0,nb233_fjx(%esp)
        movaps %xmm1,nb233_fjy(%esp)
        movaps %xmm2,nb233_fjz(%esp)

        ## Done with H2, do M interactions
        movaps  nb233_rinvM(%esp),%xmm7
        movaps  %xmm7,%xmm4
        mulps   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq
        movaps %xmm7,%xmm0
        movaps nb233_krsqM(%esp),%xmm1
        addps  %xmm1,%xmm0
        subps  nb233_crf(%esp),%xmm0   ## xmm0=rinv+ krsq-crf 
        mulps  nb233_two(%esp),%xmm1
        subps  %xmm1,%xmm7
        mulps  nb233_qqM(%esp),%xmm0
        mulps  nb233_qqM(%esp),%xmm7

        mulps  %xmm7,%xmm4      ## total fs M in xmm4 

        addps  nb233_vctot(%esp),%xmm0
        movaps %xmm0,nb233_vctot(%esp)

        movaps nb233_dxM(%esp),%xmm0
        movaps nb233_dyM(%esp),%xmm1
        movaps nb233_dzM(%esp),%xmm2
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2

        ## update M forces 
        movaps nb233_fixM(%esp),%xmm3
        movaps nb233_fiyM(%esp),%xmm4
        movaps nb233_fizM(%esp),%xmm7
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm7
        movaps %xmm3,nb233_fixM(%esp)
        movaps %xmm4,nb233_fiyM(%esp)
        movaps %xmm7,nb233_fizM(%esp)

        movl nb233_faction(%ebp),%edi
        ## update j forces from stored values
        addps nb233_fjx(%esp),%xmm0
        addps nb233_fjy(%esp),%xmm1
        addps nb233_fjz(%esp),%xmm2

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
        subl $4,nb233_innerk(%esp)
        jl    _nb_kernel233_ia32_sse.nb233_odd_inner
        jmp   _nb_kernel233_ia32_sse.nb233_unroll_loop
_nb_kernel233_ia32_sse.nb233_odd_inner: 
        addl $4,nb233_innerk(%esp)
        jnz   _nb_kernel233_ia32_sse.nb233_odd_loop
        jmp   _nb_kernel233_ia32_sse.nb233_updateouterdata
_nb_kernel233_ia32_sse.nb233_odd_loop: 
        movl  nb233_innerjjnr(%esp),%edx        ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb233_innerjjnr(%esp)

        xorps %xmm4,%xmm4       ## clear reg.
        movss nb233_iqM(%esp),%xmm4
        movl nb233_charge(%ebp),%esi
        movhps nb233_iqH(%esp),%xmm4    ## [qM  0  qH  qH] 
        shufps $41,%xmm4,%xmm4 ## [0 qH qH qM]

        movss (%esi,%eax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb233_qqM(%esp)    ## use dummy qq for storage 

        xorps %xmm6,%xmm6
        movl nb233_type(%ebp),%esi
        movl (%esi,%eax,4),%ebx
        movl nb233_vdwparam(%ebp),%esi
        shll %ebx
        addl nb233_ntia(%esp),%ebx
        movlps (%esi,%ebx,4),%xmm6
        movaps %xmm6,%xmm7
        shufps $252,%xmm6,%xmm6 ## constant 11111100
        shufps $253,%xmm7,%xmm7 ## constant 11111101
        movaps %xmm6,nb233_c6(%esp)
        movaps %xmm7,nb233_c12(%esp)

        movl nb233_pos(%ebp),%esi
        leal (%eax,%eax,2),%eax

        movss nb233_ixO(%esp),%xmm3
        movss nb233_iyO(%esp),%xmm4
        movss nb233_izO(%esp),%xmm5
        movss nb233_ixH1(%esp),%xmm0
        movss nb233_iyH1(%esp),%xmm1
        movss nb233_izH1(%esp),%xmm2
        unpcklps nb233_ixH2(%esp),%xmm3         ## ixO ixH2 - -
        unpcklps nb233_iyH2(%esp),%xmm4         ## iyO iyH2 - -
        unpcklps nb233_izH2(%esp),%xmm5         ## izO izH2 - -
        unpcklps nb233_ixM(%esp),%xmm0          ## ixH1 ixM - -
        unpcklps nb233_iyM(%esp),%xmm1          ## iyH1 iyM - -
        unpcklps nb233_izM(%esp),%xmm2          ## izH1 izM - -
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
        movaps %xmm3,nb233_dxO(%esp)
        movaps %xmm4,nb233_dyO(%esp)
        movaps %xmm5,nb233_dzO(%esp)

        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5

        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        ## rsq in xmm4 
        movaps %xmm4,%xmm0
        mulps nb233_krf(%esp),%xmm0
        movaps %xmm0,nb233_krsqM(%esp)

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb233_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb233_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv, xmm4=rsq

        mulps %xmm0,%xmm4
        mulps  nb233_tsc(%esp),%xmm4   ## rtab

        cvttps2pi %xmm4,%mm6
        cvtpi2ps %mm6,%xmm6
        subss  %xmm6,%xmm4
        movss %xmm4,%xmm1       ## xmm1=eps 
        movss %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $3,%mm6

        movd %eax,%mm0

        movl nb233_VFtab(%ebp),%esi
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
        mulss  nb233_two(%esp),%xmm7    ## two*Heps2 
        addss  %xmm6,%xmm7
        addss  %xmm5,%xmm7 ## xmm7=FF 
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 

        movss nb233_c6(%esp),%xmm4
        mulss  %xmm4,%xmm7       ## fijD 
        mulss  %xmm4,%xmm5       ## Vvdw6 
        mulss  nb233_tsc(%esp),%xmm7
        ## put scalar force on stack Update Vvdwtot directly 
        addss  nb233_Vvdwtot(%esp),%xmm5
        movss %xmm7,nb233_fstmp(%esp)
        movss %xmm5,nb233_Vvdwtot(%esp)

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
        mulss  nb233_two(%esp),%xmm7    ## two*Heps2 
        addss  %xmm6,%xmm7
        addss  %xmm5,%xmm7 ## xmm7=FF 
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 

        movss nb233_c12(%esp),%xmm4
        mulss  %xmm4,%xmm7 ## fijR 
        mulss  %xmm4,%xmm5 ## Vvdw12 
        mulss  nb233_tsc(%esp),%xmm7
        addss  nb233_fstmp(%esp),%xmm7
        movss %xmm7,nb233_fstmp(%esp)
        addss  nb233_Vvdwtot(%esp),%xmm5
        movss %xmm5,nb233_Vvdwtot(%esp)

        movd %mm0,%eax

        movaps %xmm0,%xmm4
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm5
        mulps  %xmm4,%xmm4      ## xmm1=rinv, xmm4=rinvsq
        movaps nb233_krsqM(%esp),%xmm3
        addps  %xmm3,%xmm5      ## xmm0=rinv+ krsq 
        subps  nb233_crf(%esp),%xmm5   ## xmm0=rinv+ krsq-crf 
        mulps  nb233_two(%esp),%xmm3
        subps  %xmm3,%xmm1      ## xmm1=rinv-2*krsq
        movaps %xmm5,%xmm2
        mulps  nb233_qqM(%esp),%xmm2    ## xmm2=vcoul 
        mulps  nb233_qqM(%esp),%xmm1    ## xmm1=coul part of fs 
        mulps  %xmm0,%xmm1
        subss  nb233_fstmp(%esp),%xmm1
        mulps  %xmm0,%xmm1
        movaps %xmm1,%xmm4

        addps  nb233_vctot(%esp),%xmm2
        movaps %xmm2,nb233_vctot(%esp)


        movaps nb233_dxO(%esp),%xmm0
        movaps nb233_dyO(%esp),%xmm1
        movaps nb233_dzO(%esp),%xmm2

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2 ## xmm0-xmm2 now contains tx-tz (partial force)

        movss  nb233_fixO(%esp),%xmm3
        movss  nb233_fiyO(%esp),%xmm4
        movss  nb233_fizO(%esp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,nb233_fixO(%esp)
        movss  %xmm4,nb233_fiyO(%esp)
        movss  %xmm5,nb233_fizO(%esp)   ## updated the O force now do the H's

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        shufps $0x39,%xmm3,%xmm3 ## shift right 
        shufps $0x39,%xmm4,%xmm4
        shufps $0x39,%xmm5,%xmm5
        addss  nb233_fixH1(%esp),%xmm3
        addss  nb233_fiyH1(%esp),%xmm4
        addss  nb233_fizH1(%esp),%xmm5
        movss  %xmm3,nb233_fixH1(%esp)
        movss  %xmm4,nb233_fiyH1(%esp)
        movss  %xmm5,nb233_fizH1(%esp)          ## updated the H1 force 

        shufps $0x39,%xmm3,%xmm3
        shufps $0x39,%xmm4,%xmm4
        shufps $0x39,%xmm5,%xmm5
        addss  nb233_fixH2(%esp),%xmm3
        addss  nb233_fiyH2(%esp),%xmm4
        addss  nb233_fizH2(%esp),%xmm5
        movss  %xmm3,nb233_fixH2(%esp)
        movss  %xmm4,nb233_fiyH2(%esp)
        movss  %xmm5,nb233_fizH2(%esp)          ## updated the H2 force 

        movl nb233_faction(%ebp),%edi
        shufps $0x39,%xmm3,%xmm3
        shufps $0x39,%xmm4,%xmm4
        shufps $0x39,%xmm5,%xmm5
        addss  nb233_fixM(%esp),%xmm3
        addss  nb233_fiyM(%esp),%xmm4
        addss  nb233_fizM(%esp),%xmm5
        movss  %xmm3,nb233_fixM(%esp)
        movss  %xmm4,nb233_fiyM(%esp)
        movss  %xmm5,nb233_fizM(%esp)   ## updated the M force 

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

        decl nb233_innerk(%esp)
        jz    _nb_kernel233_ia32_sse.nb233_updateouterdata
        jmp   _nb_kernel233_ia32_sse.nb233_odd_loop
_nb_kernel233_ia32_sse.nb233_updateouterdata: 
        movl  nb233_ii3(%esp),%ecx
        movl  nb233_faction(%ebp),%edi
        movl  nb233_fshift(%ebp),%esi
        movl  nb233_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movaps nb233_fixO(%esp),%xmm0
        movaps nb233_fiyO(%esp),%xmm1
        movaps nb233_fizO(%esp),%xmm2

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
        movaps nb233_fixH1(%esp),%xmm0
        movaps nb233_fiyH1(%esp),%xmm1
        movaps nb233_fizH1(%esp),%xmm2

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
        movaps nb233_fixH2(%esp),%xmm0
        movaps nb233_fiyH2(%esp),%xmm1
        movaps nb233_fizH2(%esp),%xmm2

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
        movaps nb233_fixM(%esp),%xmm0
        movaps nb233_fiyM(%esp),%xmm1
        movaps nb233_fizM(%esp),%xmm2

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
        movl nb233_n(%esp),%esi
        ## get group index for i particle 
        movl  nb233_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb233_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb233_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb233_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb233_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb233_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel233_ia32_sse.nb233_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb233_n(%esp)
        jmp _nb_kernel233_ia32_sse.nb233_outer
_nb_kernel233_ia32_sse.nb233_outerend: 
        ## check if more outer neighborlists remain
        movl  nb233_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel233_ia32_sse.nb233_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel233_ia32_sse.nb233_threadloop
_nb_kernel233_ia32_sse.nb233_end: 
        emms

        movl nb233_nouter(%esp),%eax
        movl nb233_ninner(%esp),%ebx
        movl nb233_outeriter(%ebp),%ecx
        movl nb233_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb233_salign(%esp),%eax
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


.globl nb_kernel233nf_ia32_sse
.globl _nb_kernel233nf_ia32_sse
nb_kernel233nf_ia32_sse:        
_nb_kernel233nf_ia32_sse:       
.set nb233nf_p_nri, 8
.set nb233nf_iinr, 12
.set nb233nf_jindex, 16
.set nb233nf_jjnr, 20
.set nb233nf_shift, 24
.set nb233nf_shiftvec, 28
.set nb233nf_fshift, 32
.set nb233nf_gid, 36
.set nb233nf_pos, 40
.set nb233nf_faction, 44
.set nb233nf_charge, 48
.set nb233nf_p_facel, 52
.set nb233nf_argkrf, 56
.set nb233nf_argcrf, 60
.set nb233nf_Vc, 64
.set nb233nf_type, 68
.set nb233nf_p_ntype, 72
.set nb233nf_vdwparam, 76
.set nb233nf_Vvdw, 80
.set nb233nf_p_tabscale, 84
.set nb233nf_VFtab, 88
.set nb233nf_invsqrta, 92
.set nb233nf_dvda, 96
.set nb233nf_p_gbtabscale, 100
.set nb233nf_GBtab, 104
.set nb233nf_p_nthreads, 108
.set nb233nf_count, 112
.set nb233nf_mtx, 116
.set nb233nf_outeriter, 120
.set nb233nf_inneriter, 124
.set nb233nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb233nf_ixO, 0
.set nb233nf_iyO, 16
.set nb233nf_izO, 32
.set nb233nf_ixH1, 48
.set nb233nf_iyH1, 64
.set nb233nf_izH1, 80
.set nb233nf_ixH2, 96
.set nb233nf_iyH2, 112
.set nb233nf_izH2, 128
.set nb233nf_ixM, 144
.set nb233nf_iyM, 160
.set nb233nf_izM, 176
.set nb233nf_iqM, 192
.set nb233nf_iqH, 208
.set nb233nf_qqH, 224
.set nb233nf_rinvH1, 240
.set nb233nf_rinvH2, 256
.set nb233nf_rinvM, 272
.set nb233nf_c6, 288
.set nb233nf_c12, 304
.set nb233nf_tsc, 320
.set nb233nf_krf, 336
.set nb233nf_crf, 352
.set nb233nf_krsqH1, 368
.set nb233nf_krsqH2, 384
.set nb233nf_krsqM, 400
.set nb233nf_vctot, 416
.set nb233nf_Vvdwtot, 432
.set nb233nf_half, 448
.set nb233nf_three, 464
.set nb233nf_qqM, 480
.set nb233nf_is3, 496
.set nb233nf_ii3, 500
.set nb233nf_ntia, 504
.set nb233nf_innerjjnr, 508
.set nb233nf_innerk, 512
.set nb233nf_n, 516
.set nb233nf_nn1, 520
.set nb233nf_nri, 524
.set nb233nf_nouter, 528
.set nb233nf_ninner, 532
.set nb233nf_salign, 536
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
        movl %eax,nb233nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb233nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb233nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb233nf_nouter(%esp)
        movl %eax,nb233nf_ninner(%esp)

        movl nb233nf_p_tabscale(%ebp),%eax
        movss (%eax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb233nf_tsc(%esp)

        movl nb233nf_argkrf(%ebp),%esi
        movl nb233nf_argcrf(%ebp),%edi
        movss (%esi),%xmm5
        movss (%edi),%xmm6
        shufps $0,%xmm5,%xmm5
        shufps $0,%xmm6,%xmm6
        movaps %xmm5,nb233nf_krf(%esp)
        movaps %xmm6,nb233nf_crf(%esp)
        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb233nf_half(%esp)
        movss nb233nf_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb233nf_half(%esp)
        movaps %xmm3,nb233nf_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb233nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb233nf_charge(%ebp),%edx
        movss 4(%edx,%ebx,4),%xmm4
        movss 12(%edx,%ebx,4),%xmm3
        movl nb233nf_p_facel(%ebp),%esi
        movss (%esi),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb233nf_iqM(%esp)
        movaps %xmm4,nb233nf_iqH(%esp)

        movl  nb233nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb233nf_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb233nf_ntia(%esp)
_nb_kernel233nf_ia32_sse.nb233nf_threadloop: 
        movl  nb233nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel233nf_ia32_sse.nb233nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel233nf_ia32_sse.nb233nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb233nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb233nf_n(%esp)
        movl %ebx,nb233nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel233nf_ia32_sse.nb233nf_outerstart
        jmp _nb_kernel233nf_ia32_sse.nb233nf_end

_nb_kernel233nf_ia32_sse.nb233nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb233nf_nouter(%esp),%ebx
        movl %ebx,nb233nf_nouter(%esp)

_nb_kernel233nf_ia32_sse.nb233nf_outer: 
        movl  nb233nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx        ## ebx=3*is 
        movl  %ebx,nb233nf_is3(%esp)            ## store is3 

        movl  nb233nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb233nf_iinr(%ebp),%ecx           ## ecx = pointer into iinr[]    
        movl  (%ecx,%esi,4),%ebx                ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        movaps %xmm0,%xmm6
        movaps %xmm1,%xmm7

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb233nf_pos(%ebp),%eax    ## eax = base of pos[]  
        movl  %ebx,nb233nf_ii3(%esp)

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
        movaps %xmm3,nb233nf_ixO(%esp)
        movaps %xmm4,nb233nf_iyO(%esp)
        movaps %xmm5,nb233nf_izO(%esp)
        movaps %xmm6,nb233nf_ixH1(%esp)
        movaps %xmm7,nb233nf_iyH1(%esp)

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
        movaps %xmm6,nb233nf_izH1(%esp)
        movaps %xmm0,nb233nf_ixH2(%esp)
        movaps %xmm1,nb233nf_iyH2(%esp)
        movaps %xmm2,nb233nf_izH2(%esp)
        movaps %xmm3,nb233nf_ixM(%esp)
        movaps %xmm4,nb233nf_iyM(%esp)
        movaps %xmm5,nb233nf_izM(%esp)

        ## clear vctot 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb233nf_vctot(%esp)
        movaps %xmm4,nb233nf_Vvdwtot(%esp)

        movl  nb233nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx                ## jindex[n] 
        movl  4(%eax,%esi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movl  nb233nf_pos(%ebp),%esi
        movl  nb233nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb233nf_innerjjnr(%esp)      ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb233nf_ninner(%esp),%ecx
        movl  %ecx,nb233nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb233nf_innerk(%esp)         ## number of innerloop atoms 
        jge   _nb_kernel233nf_ia32_sse.nb233nf_unroll_loop
        jmp   _nb_kernel233nf_ia32_sse.nb233nf_odd_inner
_nb_kernel233nf_ia32_sse.nb233nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb233nf_innerjjnr(%esp),%edx      ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx             ## eax-edx=jnr1-4 

        addl $16,nb233nf_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb233nf_charge(%ebp),%esi  ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4              ## and in xmm4 
        mulps  nb233nf_iqM(%esp),%xmm3
        mulps  nb233nf_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps  %xmm3,nb233nf_qqM(%esp)
        movaps  %xmm4,nb233nf_qqH(%esp)

        movl nb233nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb233nf_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb233nf_ntia(%esp),%edi
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

        movaps %xmm4,nb233nf_c6(%esp)
        movaps %xmm6,nb233nf_c12(%esp)

        movl nb233nf_pos(%ebp),%esi     ## base of pos[] 

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
        movaps nb233nf_ixO(%esp),%xmm4
        movaps nb233nf_iyO(%esp),%xmm5
        movaps nb233nf_izO(%esp),%xmm6

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
        movaps nb233nf_ixH1(%esp),%xmm4
        movaps nb233nf_iyH1(%esp),%xmm5
        movaps nb233nf_izH1(%esp),%xmm6

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
        movaps nb233nf_ixH2(%esp),%xmm3
        movaps nb233nf_iyH2(%esp),%xmm4
        movaps nb233nf_izH2(%esp),%xmm5

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
        movaps nb233nf_iyM(%esp),%xmm3
        movaps nb233nf_izM(%esp),%xmm4
        subps  %xmm1,%xmm3
        subps  %xmm2,%xmm4
        movaps nb233nf_ixM(%esp),%xmm2
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
        mulps  nb233nf_krf(%esp),%xmm0
        mulps  nb233nf_krf(%esp),%xmm1
        mulps  nb233nf_krf(%esp),%xmm2
        movaps %xmm0,nb233nf_krsqM(%esp)
        movaps %xmm1,nb233nf_krsqH2(%esp)
        movaps %xmm2,nb233nf_krsqH1(%esp)

        ## rsqH1 - seed in xmm2 
        rsqrtps %xmm6,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb233nf_three(%esp),%xmm0
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb233nf_half(%esp),%xmm0
        movaps  %xmm0,nb233nf_rinvH1(%esp)      ## rinvH1 

        ## rsqH2 - seed to xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb233nf_three(%esp),%xmm0
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb233nf_half(%esp),%xmm0
        movaps  %xmm0,nb233nf_rinvH2(%esp)      ## rinvH2 

        ## rsqM - seed to xmm2 
        rsqrtps %xmm4,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb233nf_three(%esp),%xmm0
        mulps   %xmm4,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb233nf_half(%esp),%xmm0
        movaps  %xmm0,nb233nf_rinvM(%esp)

        ## Do the O LJ-only interaction directly.       
        ## rsqO is in xmm7
        rsqrtps %xmm7,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb233nf_three(%esp),%xmm4
        mulps   %xmm7,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb233nf_half(%esp),%xmm4
        movaps  %xmm4,%xmm0
        ## xmm0=rinvO

        mulps %xmm0,%xmm7
        mulps nb233nf_tsc(%esp),%xmm7   ## rtab

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

        movl nb233nf_VFtab(%ebp),%esi
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

        movaps nb233nf_c6(%esp),%xmm4
        mulps  %xmm4,%xmm5       ## Vvdw6 

        addps  nb233nf_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb233nf_Vvdwtot(%esp)

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

        movaps nb233nf_c12(%esp),%xmm4
        mulps  %xmm4,%xmm5 ## Vvdw12 

        addps  nb233nf_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb233nf_Vvdwtot(%esp)

        ## Do H1 interaction
        movaps  nb233nf_rinvH1(%esp),%xmm7
        movaps  %xmm7,%xmm4
        mulps   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq
        movaps %xmm7,%xmm0
        movaps nb233nf_krsqH1(%esp),%xmm1
        addps  %xmm1,%xmm0
        subps  nb233nf_crf(%esp),%xmm0   ## xmm0=rinv+ krsq-crf 
        mulps  nb233nf_qqH(%esp),%xmm0

        addps  nb233nf_vctot(%esp),%xmm0
        movaps %xmm0,nb233nf_vctot(%esp)

        ## Done with H1, do H2 interactions
        movaps  nb233nf_rinvH2(%esp),%xmm7
        movaps  %xmm7,%xmm4
        mulps   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq
        movaps %xmm7,%xmm0
        movaps nb233nf_krsqH2(%esp),%xmm1
        addps  %xmm1,%xmm0
        subps  nb233nf_crf(%esp),%xmm0   ## xmm0=rinv+ krsq-crf 
        mulps  nb233nf_qqH(%esp),%xmm0

        addps  nb233nf_vctot(%esp),%xmm0
        movaps %xmm0,nb233nf_vctot(%esp)

        ## Done with H2, do M interactions
        movaps  nb233nf_rinvM(%esp),%xmm7
        movaps  %xmm7,%xmm4
        mulps   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq
        movaps %xmm7,%xmm0
        movaps nb233nf_krsqM(%esp),%xmm1
        addps  %xmm1,%xmm0
        subps  nb233nf_crf(%esp),%xmm0   ## xmm0=rinv+ krsq-crf 
        mulps  nb233nf_qqM(%esp),%xmm0

        addps  nb233nf_vctot(%esp),%xmm0
        movaps %xmm0,nb233nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $4,nb233nf_innerk(%esp)
        jl    _nb_kernel233nf_ia32_sse.nb233nf_odd_inner
        jmp   _nb_kernel233nf_ia32_sse.nb233nf_unroll_loop
_nb_kernel233nf_ia32_sse.nb233nf_odd_inner: 
        addl $4,nb233nf_innerk(%esp)
        jnz   _nb_kernel233nf_ia32_sse.nb233nf_odd_loop
        jmp   _nb_kernel233nf_ia32_sse.nb233nf_updateouterdata
_nb_kernel233nf_ia32_sse.nb233nf_odd_loop: 
        movl  nb233nf_innerjjnr(%esp),%edx      ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb233nf_innerjjnr(%esp)

        xorps %xmm4,%xmm4       ## clear reg.
        movss nb233nf_iqM(%esp),%xmm4
        movl nb233nf_charge(%ebp),%esi
        movhps nb233nf_iqH(%esp),%xmm4    ## [qM  0  qH  qH] 
        shufps $41,%xmm4,%xmm4 ## [0 qH qH qM]

        movss (%esi,%eax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb233nf_qqM(%esp)          ## use dummy qq for storage 

        xorps %xmm6,%xmm6
        movl nb233nf_type(%ebp),%esi
        movl (%esi,%eax,4),%ebx
        movl nb233nf_vdwparam(%ebp),%esi
        shll %ebx
        addl nb233nf_ntia(%esp),%ebx
        movlps (%esi,%ebx,4),%xmm6
        movaps %xmm6,%xmm7
        shufps $252,%xmm6,%xmm6 ## constant 11111100
        shufps $253,%xmm7,%xmm7 ## constant 11111101
        movaps %xmm6,nb233nf_c6(%esp)
        movaps %xmm7,nb233nf_c12(%esp)

        movl nb233nf_pos(%ebp),%esi
        leal (%eax,%eax,2),%eax

        movss nb233nf_ixO(%esp),%xmm3
        movss nb233nf_iyO(%esp),%xmm4
        movss nb233nf_izO(%esp),%xmm5
        movss nb233nf_ixH1(%esp),%xmm0
        movss nb233nf_iyH1(%esp),%xmm1
        movss nb233nf_izH1(%esp),%xmm2
        unpcklps nb233nf_ixH2(%esp),%xmm3       ## ixO ixH2 - -
        unpcklps nb233nf_iyH2(%esp),%xmm4       ## iyO iyH2 - -
        unpcklps nb233nf_izH2(%esp),%xmm5       ## izO izH2 - -
        unpcklps nb233nf_ixM(%esp),%xmm0        ## ixH1 ixM - -
        unpcklps nb233nf_iyM(%esp),%xmm1        ## iyH1 iyM - -
        unpcklps nb233nf_izM(%esp),%xmm2        ## izH1 izM - -
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
        mulps nb233nf_krf(%esp),%xmm0
        movaps %xmm0,nb233nf_krsqM(%esp)

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb233nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb233nf_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv, xmm4=rsq

        mulps %xmm0,%xmm4
        mulps  nb233nf_tsc(%esp),%xmm4   ## rtab

        cvttps2pi %xmm4,%mm6
        cvtpi2ps %mm6,%xmm6
        subss  %xmm6,%xmm4
        movss %xmm4,%xmm1       ## xmm1=eps 
        movss %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $3,%mm6

        movl nb233nf_VFtab(%ebp),%esi
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

        movss nb233nf_c6(%esp),%xmm4
        mulss  %xmm4,%xmm5       ## Vvdw6 

        ## put scalar force on stack Update Vvdwtot directly 
        addss  nb233nf_Vvdwtot(%esp),%xmm5
        movss %xmm5,nb233nf_Vvdwtot(%esp)

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

        movss nb233nf_c12(%esp),%xmm4
        mulss  %xmm4,%xmm5 ## Vvdw12 

        addss  nb233nf_Vvdwtot(%esp),%xmm5
        movss %xmm5,nb233nf_Vvdwtot(%esp)

        movaps %xmm0,%xmm4
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm5
        mulps  %xmm4,%xmm4      ## xmm1=rinv, xmm4=rinvsq
        movaps nb233nf_krsqM(%esp),%xmm3
        addps  %xmm3,%xmm5      ## xmm0=rinv+ krsq 
        subps  nb233nf_crf(%esp),%xmm5   ## xmm0=rinv+ krsq-crf 
        mulps  nb233nf_qqM(%esp),%xmm5          ## xmm2=vcoul 

        addps  nb233nf_vctot(%esp),%xmm5
        movaps %xmm5,nb233nf_vctot(%esp)

        decl nb233nf_innerk(%esp)
        jz    _nb_kernel233nf_ia32_sse.nb233nf_updateouterdata
        jmp   _nb_kernel233nf_ia32_sse.nb233nf_odd_loop
_nb_kernel233nf_ia32_sse.nb233nf_updateouterdata: 
        ## get n from stack
        movl nb233nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb233nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb233nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb233nf_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb233nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb233nf_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb233nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel233nf_ia32_sse.nb233nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb233nf_n(%esp)
        jmp _nb_kernel233nf_ia32_sse.nb233nf_outer
_nb_kernel233nf_ia32_sse.nb233nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb233nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel233nf_ia32_sse.nb233nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel233nf_ia32_sse.nb233nf_threadloop
_nb_kernel233nf_ia32_sse.nb233nf_end: 
        emms

        movl nb233nf_nouter(%esp),%eax
        movl nb233nf_ninner(%esp),%ebx
        movl nb233nf_outeriter(%ebp),%ecx
        movl nb233nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb233nf_salign(%esp),%eax
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


