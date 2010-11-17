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



.globl nb_kernel313_ia32_sse
.globl _nb_kernel313_ia32_sse
nb_kernel313_ia32_sse:  
_nb_kernel313_ia32_sse: 
.set nb313_p_nri, 8
.set nb313_iinr, 12
.set nb313_jindex, 16
.set nb313_jjnr, 20
.set nb313_shift, 24
.set nb313_shiftvec, 28
.set nb313_fshift, 32
.set nb313_gid, 36
.set nb313_pos, 40
.set nb313_faction, 44
.set nb313_charge, 48
.set nb313_p_facel, 52
.set nb313_argkrf, 56
.set nb313_argcrf, 60
.set nb313_Vc, 64
.set nb313_type, 68
.set nb313_p_ntype, 72
.set nb313_vdwparam, 76
.set nb313_Vvdw, 80
.set nb313_p_tabscale, 84
.set nb313_VFtab, 88
.set nb313_invsqrta, 92
.set nb313_dvda, 96
.set nb313_p_gbtabscale, 100
.set nb313_GBtab, 104
.set nb313_p_nthreads, 108
.set nb313_count, 112
.set nb313_mtx, 116
.set nb313_outeriter, 120
.set nb313_inneriter, 124
.set nb313_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb313_ixO, 0
.set nb313_iyO, 16
.set nb313_izO, 32
.set nb313_ixH1, 48
.set nb313_iyH1, 64
.set nb313_izH1, 80
.set nb313_ixH2, 96
.set nb313_iyH2, 112
.set nb313_izH2, 128
.set nb313_ixM, 144
.set nb313_iyM, 160
.set nb313_izM, 176
.set nb313_iqM, 192
.set nb313_iqH, 208
.set nb313_dxO, 224
.set nb313_dyO, 240
.set nb313_dzO, 256
.set nb313_dxH1, 272
.set nb313_dyH1, 288
.set nb313_dzH1, 304
.set nb313_dxH2, 320
.set nb313_dyH2, 336
.set nb313_dzH2, 352
.set nb313_dxM, 368
.set nb313_dyM, 384
.set nb313_dzM, 400
.set nb313_qqM, 416
.set nb313_qqH, 432
.set nb313_rinvH1, 448
.set nb313_rinvH2, 464
.set nb313_rinvM, 480
.set nb313_rH1, 496
.set nb313_rH2, 512
.set nb313_rM, 528
.set nb313_tsc, 544
.set nb313_two, 560
.set nb313_c6, 576
.set nb313_c12, 592
.set nb313_six, 608
.set nb313_twelve, 624
.set nb313_vctot, 640
.set nb313_Vvdwtot, 656
.set nb313_fixO, 672
.set nb313_fiyO, 688
.set nb313_fizO, 704
.set nb313_fixH1, 720
.set nb313_fiyH1, 736
.set nb313_fizH1, 752
.set nb313_fixH2, 768
.set nb313_fiyH2, 784
.set nb313_fizH2, 800
.set nb313_fixM, 816
.set nb313_fiyM, 832
.set nb313_fizM, 848
.set nb313_fjx, 864
.set nb313_fjy, 880
.set nb313_fjz, 896
.set nb313_half, 912
.set nb313_three, 928
.set nb313_is3, 944
.set nb313_ii3, 948
.set nb313_ntia, 952
.set nb313_innerjjnr, 956
.set nb313_innerk, 960
.set nb313_n, 964
.set nb313_nn1, 968
.set nb313_nri, 972
.set nb313_nouter, 976
.set nb313_ninner, 980
.set nb313_salign, 984
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
        movl %eax,nb313_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb313_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb313_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb313_nouter(%esp)
        movl %eax,nb313_ninner(%esp)


        movl nb313_p_tabscale(%ebp),%eax
        movss (%eax),%xmm5
        shufps $0,%xmm5,%xmm5
        movaps %xmm5,nb313_tsc(%esp)
        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb313_half(%esp)
        movss nb313_half(%esp),%xmm1
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
        movaps %xmm1,nb313_half(%esp)
        movaps %xmm2,nb313_two(%esp)
        movaps %xmm3,nb313_three(%esp)
        movaps %xmm4,nb313_six(%esp)
        movaps %xmm5,nb313_twelve(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb313_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb313_charge(%ebp),%edx
        movss 4(%edx,%ebx,4),%xmm4
        movss 12(%edx,%ebx,4),%xmm3
        movl nb313_p_facel(%ebp),%esi
        movss (%esi),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb313_iqM(%esp)
        movaps %xmm4,nb313_iqH(%esp)

        movl  nb313_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb313_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb313_ntia(%esp)
_nb_kernel313_ia32_sse.nb313_threadloop: 
        movl  nb313_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel313_ia32_sse.nb313_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel313_ia32_sse.nb313_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb313_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb313_n(%esp)
        movl %ebx,nb313_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel313_ia32_sse.nb313_outerstart
        jmp _nb_kernel313_ia32_sse.nb313_end

_nb_kernel313_ia32_sse.nb313_outerstart: 
        ## ebx contains number of outer iterations
        addl nb313_nouter(%esp),%ebx
        movl %ebx,nb313_nouter(%esp)

_nb_kernel313_ia32_sse.nb313_outer: 
        movl  nb313_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx        ## ebx=3*is 
        movl  %ebx,nb313_is3(%esp)      ## store is3 

        movl  nb313_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb313_iinr(%ebp),%ecx             ## ecx = pointer into iinr[]    
        movl  (%ecx,%esi,4),%ebx                ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        movaps %xmm0,%xmm6
        movaps %xmm1,%xmm7

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb313_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb313_ii3(%esp)

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
        movaps %xmm3,nb313_ixO(%esp)
        movaps %xmm4,nb313_iyO(%esp)
        movaps %xmm5,nb313_izO(%esp)
        movaps %xmm6,nb313_ixH1(%esp)
        movaps %xmm7,nb313_iyH1(%esp)

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
        movaps %xmm6,nb313_izH1(%esp)
        movaps %xmm0,nb313_ixH2(%esp)
        movaps %xmm1,nb313_iyH2(%esp)
        movaps %xmm2,nb313_izH2(%esp)
        movaps %xmm3,nb313_ixM(%esp)
        movaps %xmm4,nb313_iyM(%esp)
        movaps %xmm5,nb313_izM(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb313_vctot(%esp)
        movaps %xmm4,nb313_Vvdwtot(%esp)
        movaps %xmm4,nb313_fixO(%esp)
        movaps %xmm4,nb313_fiyO(%esp)
        movaps %xmm4,nb313_fizO(%esp)
        movaps %xmm4,nb313_fixH1(%esp)
        movaps %xmm4,nb313_fiyH1(%esp)
        movaps %xmm4,nb313_fizH1(%esp)
        movaps %xmm4,nb313_fixH2(%esp)
        movaps %xmm4,nb313_fiyH2(%esp)
        movaps %xmm4,nb313_fizH2(%esp)
        movaps %xmm4,nb313_fixM(%esp)
        movaps %xmm4,nb313_fiyM(%esp)
        movaps %xmm4,nb313_fizM(%esp)

        movl  nb313_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx                ## jindex[n] 
        movl  4(%eax,%esi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movl  nb313_pos(%ebp),%esi
        movl  nb313_faction(%ebp),%edi
        movl  nb313_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb313_innerjjnr(%esp)        ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb313_ninner(%esp),%ecx
        movl  %ecx,nb313_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb313_innerk(%esp)   ## number of innerloop atoms 
        jge   _nb_kernel313_ia32_sse.nb313_unroll_loop
        jmp   _nb_kernel313_ia32_sse.nb313_odd_inner
_nb_kernel313_ia32_sse.nb313_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb313_innerjjnr(%esp),%edx        ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx             ## eax-edx=jnr1-4 

        addl $16,nb313_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb313_charge(%ebp),%esi    ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4              ## and in xmm4 
        mulps  nb313_iqM(%esp),%xmm3
        mulps  nb313_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps  %xmm3,nb313_qqM(%esp)
        movaps  %xmm4,nb313_qqH(%esp)

        movl nb313_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb313_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb313_ntia(%esp),%edi
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

        movaps %xmm4,nb313_c6(%esp)
        movaps %xmm6,nb313_c12(%esp)

        movl nb313_pos(%ebp),%esi       ## base of pos[] 

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
        movaps nb313_ixO(%esp),%xmm4
        movaps nb313_iyO(%esp),%xmm5
        movaps nb313_izO(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb313_dxO(%esp)
        movaps %xmm5,nb313_dyO(%esp)
        movaps %xmm6,nb313_dzO(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        movaps %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movaps nb313_ixH1(%esp),%xmm4
        movaps nb313_iyH1(%esp),%xmm5
        movaps nb313_izH1(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb313_dxH1(%esp)
        movaps %xmm5,nb313_dyH1(%esp)
        movaps %xmm6,nb313_dzH1(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm6
        addps %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movaps nb313_ixH2(%esp),%xmm3
        movaps nb313_iyH2(%esp),%xmm4
        movaps nb313_izH2(%esp),%xmm5

        ## calc dr 
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        ## store dr 
        movaps %xmm3,nb313_dxH2(%esp)
        movaps %xmm4,nb313_dyH2(%esp)
        movaps %xmm5,nb313_dzH2(%esp)
        ## square it 
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm4,%xmm5
        addps %xmm3,%xmm5

        ## move ixM-izM to xmm2-xmm4  
        movaps nb313_iyM(%esp),%xmm3
        movaps nb313_izM(%esp),%xmm4
        subps  %xmm1,%xmm3
        subps  %xmm2,%xmm4
        movaps nb313_ixM(%esp),%xmm2
        subps  %xmm0,%xmm2

        ## store dr 
        movaps %xmm2,nb313_dxM(%esp)
        movaps %xmm3,nb313_dyM(%esp)
        movaps %xmm4,nb313_dzM(%esp)
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
        movaps  nb313_three(%esp),%xmm0
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb313_half(%esp),%xmm0
        movaps  %xmm0,nb313_rinvH1(%esp)        ## rinvH1 in xmm4 
        mulps   %xmm0,%xmm6
        movaps  %xmm6,nb313_rH1(%esp)

        ## rsqH2 - seed to xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb313_three(%esp),%xmm0
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb313_half(%esp),%xmm0
        movaps  %xmm0,nb313_rinvH2(%esp)        ## rinvH2 in xmm4 
        mulps   %xmm0,%xmm5
        movaps  %xmm5,nb313_rH2(%esp)

        ## rsqM - seed to xmm2 
        rsqrtps %xmm4,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb313_three(%esp),%xmm0
        mulps   %xmm4,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb313_half(%esp),%xmm0
        movaps  %xmm0,nb313_rinvM(%esp)         ## rinvM in xmm5 
        mulps   %xmm0,%xmm4
        movaps  %xmm4,nb313_rM(%esp)

        ## Do the O LJ-only interaction directly.       
        rcpps   %xmm7,%xmm2
        movaps  nb313_two(%esp),%xmm1
        mulps   %xmm2,%xmm7
        subps   %xmm7,%xmm1
        mulps   %xmm1,%xmm2 ## rinvsq 
        movaps  %xmm2,%xmm0
        mulps   %xmm2,%xmm0     ## r4
        mulps   %xmm2,%xmm0     ## r6
        movaps  %xmm0,%xmm1
        mulps   %xmm1,%xmm1     ## r12
        mulps   nb313_c6(%esp),%xmm0
        mulps   nb313_c12(%esp),%xmm1
        movaps  %xmm1,%xmm3
        subps   %xmm0,%xmm3     ## Vvdw12-Vvdw6
        addps   nb313_Vvdwtot(%esp),%xmm3
        movaps  %xmm3,nb313_Vvdwtot(%esp)
        mulps   nb313_six(%esp),%xmm0
        mulps   nb313_twelve(%esp),%xmm1
        subps   %xmm0,%xmm1
        mulps   %xmm2,%xmm1     ## fscal
        movaps nb313_dxO(%esp),%xmm3
        movaps nb313_dyO(%esp),%xmm4
        movaps nb313_dzO(%esp),%xmm5
        mulps  %xmm1,%xmm3
        mulps  %xmm1,%xmm4
        mulps  %xmm1,%xmm5      ## tx in xmm3-xmm5

        ## update O forces 
        movaps nb313_fixO(%esp),%xmm0
        movaps nb313_fiyO(%esp),%xmm1
        movaps nb313_fizO(%esp),%xmm2
        addps  %xmm3,%xmm0
        addps  %xmm4,%xmm1
        addps  %xmm5,%xmm2
        movaps %xmm0,nb313_fixO(%esp)
        movaps %xmm1,nb313_fiyO(%esp)
        movaps %xmm2,nb313_fizO(%esp)
        ## update j forces with water O 
        movaps %xmm3,nb313_fjx(%esp)
        movaps %xmm4,nb313_fjy(%esp)
        movaps %xmm5,nb313_fjz(%esp)

        ## Do H1 interaction
        movl nb313_VFtab(%ebp),%esi

        movaps nb313_rH1(%esp),%xmm7
        mulps   nb313_tsc(%esp),%xmm7
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
        mulps  nb313_two(%esp),%xmm7            ## two*Heps2 
        movaps nb313_qqH(%esp),%xmm0
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm0,%xmm7 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm7 fijC 
        ## increment vcoul 
        xorps  %xmm4,%xmm4
        addps  nb313_vctot(%esp),%xmm5
        mulps  nb313_rinvH1(%esp),%xmm7
        movaps %xmm5,nb313_vctot(%esp)
        mulps  nb313_tsc(%esp),%xmm7
        subps %xmm7,%xmm4

        movaps nb313_dxH1(%esp),%xmm0
        movaps nb313_dyH1(%esp),%xmm1
        movaps nb313_dzH1(%esp),%xmm2
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2

        ## update H1 forces 
        movaps nb313_fixH1(%esp),%xmm3
        movaps nb313_fiyH1(%esp),%xmm4
        movaps nb313_fizH1(%esp),%xmm7
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm7
        movaps %xmm3,nb313_fixH1(%esp)
        movaps %xmm4,nb313_fiyH1(%esp)
        movaps %xmm7,nb313_fizH1(%esp)
        ## update j forces with water H1 
        addps  nb313_fjx(%esp),%xmm0
        addps  nb313_fjy(%esp),%xmm1
        addps  nb313_fjz(%esp),%xmm2
        movaps %xmm0,nb313_fjx(%esp)
        movaps %xmm1,nb313_fjy(%esp)
        movaps %xmm2,nb313_fjz(%esp)

        ## Done with H1, do H2 interactions 
        movaps nb313_rH2(%esp),%xmm7
        mulps   nb313_tsc(%esp),%xmm7
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
        shufps $136,%xmm7,%xmm4 ## shuffle 10001000
        shufps $221,%xmm7,%xmm5 ## shuffle 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## shuf 10001000
        shufps $221,%xmm3,%xmm7 ## shuf 11011101
        ## coulomb table ready, in xmm4-xmm7      

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp        
        mulps  nb313_two(%esp),%xmm7            ## two*Heps2 
        movaps nb313_qqH(%esp),%xmm0
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm0,%xmm7 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm0 fijC 
        ## increment vcoul 
        xorps  %xmm4,%xmm4
        addps  nb313_vctot(%esp),%xmm5
        mulps  nb313_rinvH2(%esp),%xmm7
        movaps %xmm5,nb313_vctot(%esp)
        mulps  nb313_tsc(%esp),%xmm7
        subps  %xmm7,%xmm4

        movaps nb313_dxH2(%esp),%xmm0
        movaps nb313_dyH2(%esp),%xmm1
        movaps nb313_dzH2(%esp),%xmm2
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2

        movd %mm0,%eax
        movd %mm1,%ebx
        movd %mm2,%ecx
        movd %mm3,%edx

        ## update H2 forces 
        movaps nb313_fixH2(%esp),%xmm3
        movaps nb313_fiyH2(%esp),%xmm4
        movaps nb313_fizH2(%esp),%xmm7
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm7
        movaps %xmm3,nb313_fixH2(%esp)
        movaps %xmm4,nb313_fiyH2(%esp)
        movaps %xmm7,nb313_fizH2(%esp)
        addps nb313_fjx(%esp),%xmm0
        addps nb313_fjy(%esp),%xmm1
        addps nb313_fjz(%esp),%xmm2
        movaps %xmm0,nb313_fjx(%esp)
        movaps %xmm1,nb313_fjy(%esp)
        movaps %xmm2,nb313_fjz(%esp)

        ## Done with H2, do M interactions 
        movaps nb313_rM(%esp),%xmm7
        mulps   nb313_tsc(%esp),%xmm7
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
        mulps  nb313_two(%esp),%xmm7            ## two*Heps2 
        movaps nb313_qqM(%esp),%xmm0
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm0,%xmm7 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm0 fijC 
        ## increment vcoul 
        xorps  %xmm4,%xmm4
        addps  nb313_vctot(%esp),%xmm5
        mulps  nb313_rinvM(%esp),%xmm7
        movaps %xmm5,nb313_vctot(%esp)
        mulps  nb313_tsc(%esp),%xmm7
        subps  %xmm7,%xmm4

        movaps nb313_dxM(%esp),%xmm0
        movaps nb313_dyM(%esp),%xmm1
        movaps nb313_dzM(%esp),%xmm2
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2

        movd %mm0,%eax
        movd %mm1,%ebx
        movd %mm2,%ecx
        movd %mm3,%edx

        ## update M forces 
        movaps nb313_fixM(%esp),%xmm3
        movaps nb313_fiyM(%esp),%xmm4
        movaps nb313_fizM(%esp),%xmm7
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm7
        movaps %xmm3,nb313_fixM(%esp)
        movaps %xmm4,nb313_fiyM(%esp)
        movaps %xmm7,nb313_fizM(%esp)

        movl nb313_faction(%ebp),%edi
        ## update j forces from stored values
        addps nb313_fjx(%esp),%xmm0
        addps nb313_fjy(%esp),%xmm1
        addps nb313_fjz(%esp),%xmm2

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
        subl $4,nb313_innerk(%esp)
        jl    _nb_kernel313_ia32_sse.nb313_odd_inner
        jmp   _nb_kernel313_ia32_sse.nb313_unroll_loop
_nb_kernel313_ia32_sse.nb313_odd_inner: 
        addl $4,nb313_innerk(%esp)
        jnz   _nb_kernel313_ia32_sse.nb313_odd_loop
        jmp   _nb_kernel313_ia32_sse.nb313_updateouterdata
_nb_kernel313_ia32_sse.nb313_odd_loop: 
        movl  nb313_innerjjnr(%esp),%edx        ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb313_innerjjnr(%esp)

        xorps %xmm4,%xmm4       ## clear reg.
        movss nb313_iqM(%esp),%xmm4
        movl nb313_charge(%ebp),%esi
        movhps nb313_iqH(%esp),%xmm4    ## [qM  0  qH  qH] 
        shufps $41,%xmm4,%xmm4 ## [0 qH qH qM]

        movss (%esi,%eax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb313_qqM(%esp)    ## use dummy qq for storage 

        xorps %xmm6,%xmm6
        movl nb313_type(%ebp),%esi
        movl (%esi,%eax,4),%ebx
        movl nb313_vdwparam(%ebp),%esi
        shll %ebx
        addl nb313_ntia(%esp),%ebx
        movlps (%esi,%ebx,4),%xmm6
        movaps %xmm6,%xmm7
        shufps $252,%xmm6,%xmm6 ## constant 11111100
        shufps $253,%xmm7,%xmm7 ## constant 11111101
        movaps %xmm6,nb313_c6(%esp)
        movaps %xmm7,nb313_c12(%esp)

        movl nb313_pos(%ebp),%esi
        leal (%eax,%eax,2),%eax

        movss nb313_ixO(%esp),%xmm3
        movss nb313_iyO(%esp),%xmm4
        movss nb313_izO(%esp),%xmm5
        movss nb313_ixH1(%esp),%xmm0
        movss nb313_iyH1(%esp),%xmm1
        movss nb313_izH1(%esp),%xmm2
        unpcklps nb313_ixH2(%esp),%xmm3         ## ixO ixH2 - -
        unpcklps nb313_iyH2(%esp),%xmm4         ## iyO iyH2 - -
        unpcklps nb313_izH2(%esp),%xmm5         ## izO izH2 - -
        unpcklps nb313_ixM(%esp),%xmm0          ## ixH1 ixM - -
        unpcklps nb313_iyM(%esp),%xmm1          ## iyH1 iyM - -
        unpcklps nb313_izM(%esp),%xmm2          ## izH1 izM - -
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
        movaps %xmm3,nb313_dxO(%esp)
        movaps %xmm4,nb313_dyO(%esp)
        movaps %xmm5,nb313_dzO(%esp)

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
        movaps nb313_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb313_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv

        movaps %xmm0,nb313_rinvM(%esp)
        mulps  %xmm0,%xmm4      ## r

        mulps nb313_tsc(%esp),%xmm4
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

        movl nb313_VFtab(%ebp),%esi
        psrlq $32,%mm6
        movd %mm6,%ebx
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm7,%edx


        xorps  %xmm5,%xmm5
        movlps (%esi,%ecx,4),%xmm3      ## values= Y3 F3  -  - 
        movhps (%esi,%ebx,4),%xmm5      ## values=   0  0 Y2 F2
        movhps (%esi,%edx,4),%xmm3      ## values=  Y3 F3 Y4 F4 

        movaps %xmm5,%xmm4              ## values=   0  0 Y2 F2 
        shufps $0x88,%xmm3,%xmm4       ## values=   0 Y2 Y3 Y3
        shufps $0xDD,%xmm3,%xmm5       ## values=   0 F2 F3 F4 

        xorps  %xmm7,%xmm7
        movlps 8(%esi,%ecx,4),%xmm3     ## values=  G3 H3  -  - 
        movhps 8(%esi,%ebx,4),%xmm7     ## values=   0  0 G2 H2
        movhps 8(%esi,%edx,4),%xmm3     ## values=  G3 H3 G4 H4 

        movaps %xmm7,%xmm6              ## values=   0  0 G2 H2 
        shufps $0x88,%xmm3,%xmm6       ## values=   0 G2 G3 G3
        shufps $0xDD,%xmm3,%xmm7       ## values=   0 H2 H3 H4 

        ## xmm4 =  0  Y2 Y3 Y4
        ## xmm5 =  0  F2 F3 F4
        ## xmm6 =  0  G2 G3 G4
        ## xmm7 =  0  H2 H3 H4

        ## coulomb table ready, in xmm4-xmm7      
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp        
        mulps  nb313_two(%esp),%xmm7            ## two*Heps2 
        movaps nb313_qqM(%esp),%xmm0
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm0 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm0 fijC 
        ## increment vcoul - then we can get rid of mm5 
        addps  nb313_vctot(%esp),%xmm5
        movaps %xmm5,nb313_vctot(%esp)

        ## do nontable L-J  in first element only.
        movaps nb313_rinvM(%esp),%xmm2
        mulss  %xmm2,%xmm2
        movaps %xmm2,%xmm1
        mulss  %xmm1,%xmm1
        mulss  %xmm2,%xmm1      ## xmm1=rinvsix
        xorps  %xmm4,%xmm4
        movss  %xmm1,%xmm4
        mulss  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulss  nb313_c6(%esp),%xmm1
        mulss  nb313_c12(%esp),%xmm4
        movaps %xmm4,%xmm3
        subss  %xmm1,%xmm3      ## xmm3=Vvdw12-Vvdw6 
        mulss  nb313_six(%esp),%xmm1
        mulss  nb313_twelve(%esp),%xmm4
        subss  %xmm1,%xmm4
        addss  nb313_Vvdwtot(%esp),%xmm3
        mulss  nb313_rinvM(%esp),%xmm4
        ## add back coul stuff from memory, and work on all elements again
        mulps  nb313_tsc(%esp),%xmm0
        subps  %xmm0,%xmm4
        movss %xmm3,nb313_Vvdwtot(%esp)
        mulps  nb313_rinvM(%esp),%xmm4

        movaps nb313_dxO(%esp),%xmm0
        movaps nb313_dyO(%esp),%xmm1
        movaps nb313_dzO(%esp),%xmm2

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2 ## xmm0-xmm2 now contains tx-tz (partial force)

        movss  nb313_fixO(%esp),%xmm3
        movss  nb313_fiyO(%esp),%xmm4
        movss  nb313_fizO(%esp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,nb313_fixO(%esp)
        movss  %xmm4,nb313_fiyO(%esp)
        movss  %xmm5,nb313_fizO(%esp)   ## updated the O force now do the H's


        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        shufps $0x39,%xmm3,%xmm3 ## shift right 
        shufps $0x39,%xmm4,%xmm4
        shufps $0x39,%xmm5,%xmm5
        addss  nb313_fixH1(%esp),%xmm3
        addss  nb313_fiyH1(%esp),%xmm4
        addss  nb313_fizH1(%esp),%xmm5
        movss  %xmm3,nb313_fixH1(%esp)
        movss  %xmm4,nb313_fiyH1(%esp)
        movss  %xmm5,nb313_fizH1(%esp)          ## updated the H1 force 

        shufps $0x39,%xmm3,%xmm3
        shufps $0x39,%xmm4,%xmm4
        shufps $0x39,%xmm5,%xmm5
        addss  nb313_fixH2(%esp),%xmm3
        addss  nb313_fiyH2(%esp),%xmm4
        addss  nb313_fizH2(%esp),%xmm5
        movss  %xmm3,nb313_fixH2(%esp)
        movss  %xmm4,nb313_fiyH2(%esp)
        movss  %xmm5,nb313_fizH2(%esp)          ## updated the H2 force 

        movl nb313_faction(%ebp),%edi
        shufps $0x39,%xmm3,%xmm3
        shufps $0x39,%xmm4,%xmm4
        shufps $0x39,%xmm5,%xmm5
        addss  nb313_fixM(%esp),%xmm3
        addss  nb313_fiyM(%esp),%xmm4
        addss  nb313_fizM(%esp),%xmm5
        movss  %xmm3,nb313_fixM(%esp)
        movss  %xmm4,nb313_fiyM(%esp)
        movss  %xmm5,nb313_fizM(%esp)   ## updated the M force 

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

        decl nb313_innerk(%esp)
        jz    _nb_kernel313_ia32_sse.nb313_updateouterdata
        jmp   _nb_kernel313_ia32_sse.nb313_odd_loop
_nb_kernel313_ia32_sse.nb313_updateouterdata: 
        movl  nb313_ii3(%esp),%ecx
        movl  nb313_faction(%ebp),%edi
        movl  nb313_fshift(%ebp),%esi
        movl  nb313_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movaps nb313_fixO(%esp),%xmm0
        movaps nb313_fiyO(%esp),%xmm1
        movaps nb313_fizO(%esp),%xmm2

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
        movaps nb313_fixH1(%esp),%xmm0
        movaps nb313_fiyH1(%esp),%xmm1
        movaps nb313_fizH1(%esp),%xmm2

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
        movaps nb313_fixH2(%esp),%xmm0
        movaps nb313_fiyH2(%esp),%xmm1
        movaps nb313_fizH2(%esp),%xmm2

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
        movaps nb313_fixM(%esp),%xmm0
        movaps nb313_fiyM(%esp),%xmm1
        movaps nb313_fizM(%esp),%xmm2

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
        movl nb313_n(%esp),%esi
        ## get group index for i particle 
        movl  nb313_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb313_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb313_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb313_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb313_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb313_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel313_ia32_sse.nb313_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb313_n(%esp)
        jmp _nb_kernel313_ia32_sse.nb313_outer
_nb_kernel313_ia32_sse.nb313_outerend: 
        ## check if more outer neighborlists remain
        movl  nb313_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel313_ia32_sse.nb313_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel313_ia32_sse.nb313_threadloop
_nb_kernel313_ia32_sse.nb313_end: 
        emms

        movl nb313_nouter(%esp),%eax
        movl nb313_ninner(%esp),%ebx
        movl nb313_outeriter(%ebp),%ecx
        movl nb313_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb313_salign(%esp),%eax
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







.globl nb_kernel313nf_ia32_sse
.globl _nb_kernel313nf_ia32_sse
nb_kernel313nf_ia32_sse:        
_nb_kernel313nf_ia32_sse:       
.set nb313nf_p_nri, 8
.set nb313nf_iinr, 12
.set nb313nf_jindex, 16
.set nb313nf_jjnr, 20
.set nb313nf_shift, 24
.set nb313nf_shiftvec, 28
.set nb313nf_fshift, 32
.set nb313nf_gid, 36
.set nb313nf_pos, 40
.set nb313nf_faction, 44
.set nb313nf_charge, 48
.set nb313nf_p_facel, 52
.set nb313nf_argkrf, 56
.set nb313nf_argcrf, 60
.set nb313nf_Vc, 64
.set nb313nf_type, 68
.set nb313nf_p_ntype, 72
.set nb313nf_vdwparam, 76
.set nb313nf_Vvdw, 80
.set nb313nf_p_tabscale, 84
.set nb313nf_VFtab, 88
.set nb313nf_invsqrta, 92
.set nb313nf_dvda, 96
.set nb313nf_p_gbtabscale, 100
.set nb313nf_GBtab, 104
.set nb313nf_p_nthreads, 108
.set nb313nf_count, 112
.set nb313nf_mtx, 116
.set nb313nf_outeriter, 120
.set nb313nf_inneriter, 124
.set nb313nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb313nf_ixO, 0
.set nb313nf_iyO, 16
.set nb313nf_izO, 32
.set nb313nf_ixH1, 48
.set nb313nf_iyH1, 64
.set nb313nf_izH1, 80
.set nb313nf_ixH2, 96
.set nb313nf_iyH2, 112
.set nb313nf_izH2, 128
.set nb313nf_ixM, 144
.set nb313nf_iyM, 160
.set nb313nf_izM, 176
.set nb313nf_iqM, 192
.set nb313nf_iqH, 208
.set nb313nf_qqM, 224
.set nb313nf_qqH, 240
.set nb313nf_rinvH1, 256
.set nb313nf_rinvH2, 272
.set nb313nf_rinvM, 288
.set nb313nf_rH1, 304
.set nb313nf_rH2, 320
.set nb313nf_rM, 336
.set nb313nf_tsc, 352
.set nb313nf_two, 368
.set nb313nf_c6, 384
.set nb313nf_c12, 400
.set nb313nf_vctot, 416
.set nb313nf_Vvdwtot, 432
.set nb313nf_half, 448
.set nb313nf_three, 464
.set nb313nf_is3, 480
.set nb313nf_ii3, 484
.set nb313nf_ntia, 488
.set nb313nf_innerjjnr, 492
.set nb313nf_innerk, 496
.set nb313nf_n, 500
.set nb313nf_nn1, 504
.set nb313nf_nri, 508
.set nb313nf_nouter, 512
.set nb313nf_ninner, 516
.set nb313nf_salign, 520
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $524,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb313nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb313nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb313nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb313nf_nouter(%esp)
        movl %eax,nb313nf_ninner(%esp)


        movl nb313nf_p_tabscale(%ebp),%eax
        movss (%eax),%xmm5
        shufps $0,%xmm5,%xmm5
        movaps %xmm5,nb313nf_tsc(%esp)
        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb313nf_half(%esp)
        movss nb313nf_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb313nf_half(%esp)
        movaps %xmm2,nb313nf_two(%esp)
        movaps %xmm3,nb313nf_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb313nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb313nf_charge(%ebp),%edx
        movss 4(%edx,%ebx,4),%xmm4
        movss 12(%edx,%ebx,4),%xmm3
        movl nb313nf_p_facel(%ebp),%esi
        movss (%esi),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb313nf_iqM(%esp)
        movaps %xmm4,nb313nf_iqH(%esp)

        movl  nb313nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb313nf_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb313nf_ntia(%esp)
_nb_kernel313nf_ia32_sse.nb313nf_threadloop: 
        movl  nb313nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel313nf_ia32_sse.nb313nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel313nf_ia32_sse.nb313nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb313nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb313nf_n(%esp)
        movl %ebx,nb313nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel313nf_ia32_sse.nb313nf_outerstart
        jmp _nb_kernel313nf_ia32_sse.nb313nf_end

_nb_kernel313nf_ia32_sse.nb313nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb313nf_nouter(%esp),%ebx
        movl %ebx,nb313nf_nouter(%esp)

_nb_kernel313nf_ia32_sse.nb313nf_outer: 
        movl  nb313nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx        ## ebx=3*is 
        movl  %ebx,nb313nf_is3(%esp)            ## store is3 

        movl  nb313nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb313nf_iinr(%ebp),%ecx           ## ecx = pointer into iinr[]    
        movl  (%ecx,%esi,4),%ebx                ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        movaps %xmm0,%xmm6
        movaps %xmm1,%xmm7

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb313nf_pos(%ebp),%eax    ## eax = base of pos[]  
        movl  %ebx,nb313nf_ii3(%esp)

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
        movaps %xmm3,nb313nf_ixO(%esp)
        movaps %xmm4,nb313nf_iyO(%esp)
        movaps %xmm5,nb313nf_izO(%esp)
        movaps %xmm6,nb313nf_ixH1(%esp)
        movaps %xmm7,nb313nf_iyH1(%esp)

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
        movaps %xmm6,nb313nf_izH1(%esp)
        movaps %xmm0,nb313nf_ixH2(%esp)
        movaps %xmm1,nb313nf_iyH2(%esp)
        movaps %xmm2,nb313nf_izH2(%esp)
        movaps %xmm3,nb313nf_ixM(%esp)
        movaps %xmm4,nb313nf_iyM(%esp)
        movaps %xmm5,nb313nf_izM(%esp)

        ## clear vctot
        xorps %xmm4,%xmm4
        movaps %xmm4,nb313nf_vctot(%esp)
        movaps %xmm4,nb313nf_Vvdwtot(%esp)

        movl  nb313nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx                ## jindex[n] 
        movl  4(%eax,%esi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movl  nb313nf_pos(%ebp),%esi
        movl  nb313nf_faction(%ebp),%edi
        movl  nb313nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb313nf_innerjjnr(%esp)      ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb313nf_ninner(%esp),%ecx
        movl  %ecx,nb313nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb313nf_innerk(%esp)         ## number of innerloop atoms 
        jge   _nb_kernel313nf_ia32_sse.nb313nf_unroll_loop
        jmp   _nb_kernel313nf_ia32_sse.nb313nf_odd_inner
_nb_kernel313nf_ia32_sse.nb313nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb313nf_innerjjnr(%esp),%edx      ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx             ## eax-edx=jnr1-4 

        addl $16,nb313nf_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb313nf_charge(%ebp),%esi  ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4              ## and in xmm4 
        mulps  nb313nf_iqM(%esp),%xmm3
        mulps  nb313nf_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps  %xmm3,nb313nf_qqM(%esp)
        movaps  %xmm4,nb313nf_qqH(%esp)

        movl nb313nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb313nf_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb313nf_ntia(%esp),%edi
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

        movaps %xmm4,nb313nf_c6(%esp)
        movaps %xmm6,nb313nf_c12(%esp)

        movl nb313nf_pos(%ebp),%esi     ## base of pos[] 

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
        movaps nb313nf_ixO(%esp),%xmm4
        movaps nb313nf_iyO(%esp),%xmm5
        movaps nb313nf_izO(%esp),%xmm6

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
        movaps nb313nf_ixH1(%esp),%xmm4
        movaps nb313nf_iyH1(%esp),%xmm5
        movaps nb313nf_izH1(%esp),%xmm6

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
        movaps nb313nf_ixH2(%esp),%xmm3
        movaps nb313nf_iyH2(%esp),%xmm4
        movaps nb313nf_izH2(%esp),%xmm5

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
        movaps nb313nf_iyM(%esp),%xmm3
        movaps nb313nf_izM(%esp),%xmm4
        subps  %xmm1,%xmm3
        subps  %xmm2,%xmm4
        movaps nb313nf_ixM(%esp),%xmm2
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
        movaps  nb313nf_three(%esp),%xmm0
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb313nf_half(%esp),%xmm0
        movaps  %xmm0,nb313nf_rinvH1(%esp)      ## rinvH1 in xmm4 
        mulps   %xmm0,%xmm6
        movaps  %xmm6,nb313nf_rH1(%esp)

        ## rsqH2 - seed to xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb313nf_three(%esp),%xmm0
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb313nf_half(%esp),%xmm0
        movaps  %xmm0,nb313nf_rinvH2(%esp)      ## rinvH2 in xmm4 
        mulps   %xmm0,%xmm5
        movaps  %xmm5,nb313nf_rH2(%esp)

        ## rsqM - seed to xmm2 
        rsqrtps %xmm4,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb313nf_three(%esp),%xmm0
        mulps   %xmm4,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb313nf_half(%esp),%xmm0
        movaps  %xmm0,nb313nf_rinvM(%esp)       ## rinvM in xmm5 
        mulps   %xmm0,%xmm4
        movaps  %xmm4,nb313nf_rM(%esp)

        ## Do the O LJ-only interaction directly.       
        rcpps   %xmm7,%xmm2
        movaps  nb313nf_two(%esp),%xmm1
        mulps   %xmm2,%xmm7
        subps   %xmm7,%xmm1
        mulps   %xmm1,%xmm2 ## rinvsq 
        movaps  %xmm2,%xmm0
        mulps   %xmm2,%xmm0     ## r4
        mulps   %xmm2,%xmm0     ## r6
        movaps  %xmm0,%xmm1
        mulps   %xmm1,%xmm1     ## r12
        mulps   nb313nf_c6(%esp),%xmm0
        mulps   nb313nf_c12(%esp),%xmm1
        movaps  %xmm1,%xmm3
        subps   %xmm0,%xmm3     ## Vvdw12-Vvdw6
        addps   nb313nf_Vvdwtot(%esp),%xmm3
        movaps  %xmm3,nb313nf_Vvdwtot(%esp)

        ## Do H1 interaction
        movl nb313nf_VFtab(%ebp),%esi

        movaps nb313nf_rH1(%esp),%xmm7
        mulps   nb313nf_tsc(%esp),%xmm7
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
        movaps nb313nf_qqH(%esp),%xmm0
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        ## increment vcoul 
        addps  nb313nf_vctot(%esp),%xmm5
        movaps %xmm5,nb313nf_vctot(%esp)

        ## Done with H1, do H2 interactions 
        movaps nb313nf_rH2(%esp),%xmm7
        mulps   nb313nf_tsc(%esp),%xmm7
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
        shufps $136,%xmm7,%xmm4 ## shuffle 10001000
        shufps $221,%xmm7,%xmm5 ## shuffle 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## shuf 10001000
        shufps $221,%xmm3,%xmm7 ## shuf 11011101
        ## coulomb table ready, in xmm4-xmm7      

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp        
        movaps nb313nf_qqH(%esp),%xmm0
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        ## increment vcoul 
        addps  nb313nf_vctot(%esp),%xmm5
        movaps %xmm5,nb313nf_vctot(%esp)

        ## Done with H2, do M interactions 
        movaps nb313nf_rM(%esp),%xmm7
        mulps   nb313nf_tsc(%esp),%xmm7
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
        movaps nb313nf_qqM(%esp),%xmm0
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        ## increment vcoul 
        addps  nb313nf_vctot(%esp),%xmm5
        movaps %xmm5,nb313nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $4,nb313nf_innerk(%esp)
        jl    _nb_kernel313nf_ia32_sse.nb313nf_odd_inner
        jmp   _nb_kernel313nf_ia32_sse.nb313nf_unroll_loop
_nb_kernel313nf_ia32_sse.nb313nf_odd_inner: 
        addl $4,nb313nf_innerk(%esp)
        jnz   _nb_kernel313nf_ia32_sse.nb313nf_odd_loop
        jmp   _nb_kernel313nf_ia32_sse.nb313nf_updateouterdata
_nb_kernel313nf_ia32_sse.nb313nf_odd_loop: 
        movl  nb313nf_innerjjnr(%esp),%edx      ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb313nf_innerjjnr(%esp)

        xorps %xmm4,%xmm4       ## clear reg.
        movss nb313nf_iqM(%esp),%xmm4
        movl nb313nf_charge(%ebp),%esi
        movhps nb313nf_iqH(%esp),%xmm4    ## [qM  0  qH  qH] 
        shufps $41,%xmm4,%xmm4 ## [0 qH qH qM]

        movss (%esi,%eax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb313nf_qqM(%esp)          ## use dummy qq for storage 

        xorps %xmm6,%xmm6
        movl nb313nf_type(%ebp),%esi
        movl (%esi,%eax,4),%ebx
        movl nb313nf_vdwparam(%ebp),%esi
        shll %ebx
        addl nb313nf_ntia(%esp),%ebx
        movlps (%esi,%ebx,4),%xmm6
        movaps %xmm6,%xmm7
        shufps $252,%xmm6,%xmm6 ## constant 11111100
        shufps $253,%xmm7,%xmm7 ## constant 11111101
        movaps %xmm6,nb313nf_c6(%esp)
        movaps %xmm7,nb313nf_c12(%esp)

        movl nb313nf_pos(%ebp),%esi
        leal (%eax,%eax,2),%eax

        movss nb313nf_ixO(%esp),%xmm3
        movss nb313nf_iyO(%esp),%xmm4
        movss nb313nf_izO(%esp),%xmm5
        movss nb313nf_ixH1(%esp),%xmm0
        movss nb313nf_iyH1(%esp),%xmm1
        movss nb313nf_izH1(%esp),%xmm2
        unpcklps nb313nf_ixH2(%esp),%xmm3       ## ixO ixH2 - -
        unpcklps nb313nf_iyH2(%esp),%xmm4       ## iyO iyH2 - -
        unpcklps nb313nf_izH2(%esp),%xmm5       ## izO izH2 - -
        unpcklps nb313nf_ixM(%esp),%xmm0        ## ixH1 ixM - -
        unpcklps nb313nf_iyM(%esp),%xmm1        ## iyH1 iyM - -
        unpcklps nb313nf_izM(%esp),%xmm2        ## izH1 izM - -
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
        movaps nb313nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb313nf_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv

        movaps %xmm0,nb313nf_rinvM(%esp)
        mulps  %xmm0,%xmm4      ## r

        mulps nb313nf_tsc(%esp),%xmm4
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

        movl nb313nf_VFtab(%ebp),%esi
        psrlq $32,%mm6
        movd %mm6,%ebx
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm7,%edx

        xorps  %xmm5,%xmm5
        movlps (%esi,%ecx,4),%xmm3      ## values=  Y3 F3  -  - 
        movhps (%esi,%ebx,4),%xmm5      ## values=   0  0 Y2 F2
        movhps (%esi,%edx,4),%xmm3      ## values=  Y3 F3 Y4 F4 

        movaps %xmm5,%xmm4              ## values=   0  0 Y2 F2 
        shufps $0x88,%xmm3,%xmm4       ## values=   0 Y2 Y3 Y3
        shufps $0xDD,%xmm3,%xmm5       ## values=   0 F2 F3 F4 

        xorps  %xmm7,%xmm7
        movlps 8(%esi,%ecx,4),%xmm3     ## values=  G3 H3  -  - 
        movhps 8(%esi,%ebx,4),%xmm7     ## values=   0  0 G2 H2
        movhps 8(%esi,%edx,4),%xmm3     ## values=  G3 H3 G4 H4 

        movaps %xmm7,%xmm6              ## values=   0  0 G2 H2 
        shufps $0x88,%xmm3,%xmm6       ## values=   0 G2 G3 G3
        shufps $0xDD,%xmm3,%xmm7       ## values=   0 H2 H3 H4 

        ## xmm4 =  0  Y2 Y3 Y4
        ## xmm5 =  0  F2 F3 F4
        ## xmm6 =  0  G2 G3 G4
        ## xmm7 =  0  H2 H3 H4

        ## coulomb table ready, in xmm4-xmm7      
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp        
        movaps nb313nf_qqM(%esp),%xmm0
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        ## increment vcoul 
        addps  nb313nf_vctot(%esp),%xmm5
        movaps %xmm5,nb313nf_vctot(%esp)

        ## do nontable L-J  in first element only.
        movaps nb313nf_rinvM(%esp),%xmm2
        mulss  %xmm2,%xmm2
        movaps %xmm2,%xmm1
        mulss  %xmm1,%xmm1
        mulss  %xmm2,%xmm1      ## xmm1=rinvsix
        xorps  %xmm4,%xmm4
        movss  %xmm1,%xmm4
        mulss  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulss  nb313nf_c6(%esp),%xmm1
        mulss  nb313nf_c12(%esp),%xmm4
        movaps %xmm4,%xmm3
        subss  %xmm1,%xmm3      ## xmm3=Vvdw12-Vvdw6
        addss  nb313nf_Vvdwtot(%esp),%xmm3
        movss %xmm3,nb313nf_Vvdwtot(%esp)

        decl nb313nf_innerk(%esp)
        jz    _nb_kernel313nf_ia32_sse.nb313nf_updateouterdata
        jmp   _nb_kernel313nf_ia32_sse.nb313nf_odd_loop
_nb_kernel313nf_ia32_sse.nb313nf_updateouterdata: 
        ## get n from stack
        movl nb313nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb313nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb313nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb313nf_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb313nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb313nf_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb313nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel313nf_ia32_sse.nb313nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb313nf_n(%esp)
        jmp _nb_kernel313nf_ia32_sse.nb313nf_outer
_nb_kernel313nf_ia32_sse.nb313nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb313nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel313nf_ia32_sse.nb313nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel313nf_ia32_sse.nb313nf_threadloop
_nb_kernel313nf_ia32_sse.nb313nf_end: 
        emms

        movl nb313nf_nouter(%esp),%eax
        movl nb313nf_ninner(%esp),%ebx
        movl nb313nf_outeriter(%ebp),%ecx
        movl nb313nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb313nf_salign(%esp),%eax
        addl %eax,%esp
        addl $524,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret





