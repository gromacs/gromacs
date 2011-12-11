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



.globl nb_kernel133_ia32_sse
.globl _nb_kernel133_ia32_sse
nb_kernel133_ia32_sse:  
_nb_kernel133_ia32_sse: 
.set nb133_p_nri, 8
.set nb133_iinr, 12
.set nb133_jindex, 16
.set nb133_jjnr, 20
.set nb133_shift, 24
.set nb133_shiftvec, 28
.set nb133_fshift, 32
.set nb133_gid, 36
.set nb133_pos, 40
.set nb133_faction, 44
.set nb133_charge, 48
.set nb133_p_facel, 52
.set nb133_argkrf, 56
.set nb133_argcrf, 60
.set nb133_Vc, 64
.set nb133_type, 68
.set nb133_p_ntype, 72
.set nb133_vdwparam, 76
.set nb133_Vvdw, 80
.set nb133_p_tabscale, 84
.set nb133_VFtab, 88
.set nb133_invsqrta, 92
.set nb133_dvda, 96
.set nb133_p_gbtabscale, 100
.set nb133_GBtab, 104
.set nb133_p_nthreads, 108
.set nb133_count, 112
.set nb133_mtx, 116
.set nb133_outeriter, 120
.set nb133_inneriter, 124
.set nb133_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb133_ixO, 0
.set nb133_iyO, 16
.set nb133_izO, 32
.set nb133_ixH1, 48
.set nb133_iyH1, 64
.set nb133_izH1, 80
.set nb133_ixH2, 96
.set nb133_iyH2, 112
.set nb133_izH2, 128
.set nb133_ixM, 144
.set nb133_iyM, 160
.set nb133_izM, 176
.set nb133_iqM, 192
.set nb133_iqH, 208
.set nb133_dxO, 224
.set nb133_dyO, 240
.set nb133_dzO, 256
.set nb133_dxH1, 272
.set nb133_dyH1, 288
.set nb133_dzH1, 304
.set nb133_dxH2, 320
.set nb133_dyH2, 336
.set nb133_dzH2, 352
.set nb133_dxM, 368
.set nb133_dyM, 384
.set nb133_dzM, 400
.set nb133_qqM, 416
.set nb133_qqH, 432
.set nb133_rinvH1, 448
.set nb133_rinvH2, 464
.set nb133_rinvM, 480
.set nb133_two, 496
.set nb133_c6, 512
.set nb133_c12, 528
.set nb133_tsc, 544
.set nb133_fstmp, 560
.set nb133_vctot, 656
.set nb133_Vvdwtot, 672
.set nb133_fixO, 688
.set nb133_fiyO, 704
.set nb133_fizO, 720
.set nb133_fixH1, 736
.set nb133_fiyH1, 752
.set nb133_fizH1, 768
.set nb133_fixH2, 784
.set nb133_fiyH2, 800
.set nb133_fizH2, 816
.set nb133_fixM, 832
.set nb133_fiyM, 848
.set nb133_fizM, 864
.set nb133_fjx, 880
.set nb133_fjy, 896
.set nb133_fjz, 912
.set nb133_half, 928
.set nb133_three, 944
.set nb133_is3, 960
.set nb133_ii3, 964
.set nb133_ntia, 968
.set nb133_innerjjnr, 972
.set nb133_innerk, 976
.set nb133_n, 980
.set nb133_nn1, 984
.set nb133_nri, 988
.set nb133_nouter, 992
.set nb133_ninner, 996
.set nb133_salign, 1000
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
        movl %eax,nb133_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb133_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb133_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb133_nouter(%esp)
        movl %eax,nb133_ninner(%esp)

        movl nb133_p_tabscale(%ebp),%eax
        movss (%eax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb133_tsc(%esp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb133_half(%esp)
        movss nb133_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb133_half(%esp)
        movaps %xmm2,nb133_two(%esp)
        movaps %xmm3,nb133_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb133_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb133_charge(%ebp),%edx
        movss 4(%edx,%ebx,4),%xmm4
        movss 12(%edx,%ebx,4),%xmm3
        movl nb133_p_facel(%ebp),%esi
        movss (%esi),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb133_iqM(%esp)
        movaps %xmm4,nb133_iqH(%esp)

        movl  nb133_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb133_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb133_ntia(%esp)
_nb_kernel133_ia32_sse.nb133_threadloop: 
        movl  nb133_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel133_ia32_sse.nb133_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel133_ia32_sse.nb133_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb133_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb133_n(%esp)
        movl %ebx,nb133_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel133_ia32_sse.nb133_outerstart
        jmp _nb_kernel133_ia32_sse.nb133_end

_nb_kernel133_ia32_sse.nb133_outerstart: 
        ## ebx contains number of outer iterations
        addl nb133_nouter(%esp),%ebx
        movl %ebx,nb133_nouter(%esp)

_nb_kernel133_ia32_sse.nb133_outer: 
        movl  nb133_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx        ## ebx=3*is 
        movl  %ebx,nb133_is3(%esp)      ## store is3 

        movl  nb133_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb133_iinr(%ebp),%ecx             ## ecx = pointer into iinr[]    
        movl  (%ecx,%esi,4),%ebx                ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        movaps %xmm0,%xmm6
        movaps %xmm1,%xmm7

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb133_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb133_ii3(%esp)

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
        movaps %xmm3,nb133_ixO(%esp)
        movaps %xmm4,nb133_iyO(%esp)
        movaps %xmm5,nb133_izO(%esp)
        movaps %xmm6,nb133_ixH1(%esp)
        movaps %xmm7,nb133_iyH1(%esp)

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
        movaps %xmm6,nb133_izH1(%esp)
        movaps %xmm0,nb133_ixH2(%esp)
        movaps %xmm1,nb133_iyH2(%esp)
        movaps %xmm2,nb133_izH2(%esp)
        movaps %xmm3,nb133_ixM(%esp)
        movaps %xmm4,nb133_iyM(%esp)
        movaps %xmm5,nb133_izM(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb133_vctot(%esp)
        movaps %xmm4,nb133_Vvdwtot(%esp)
        movaps %xmm4,nb133_fixO(%esp)
        movaps %xmm4,nb133_fiyO(%esp)
        movaps %xmm4,nb133_fizO(%esp)
        movaps %xmm4,nb133_fixH1(%esp)
        movaps %xmm4,nb133_fiyH1(%esp)
        movaps %xmm4,nb133_fizH1(%esp)
        movaps %xmm4,nb133_fixH2(%esp)
        movaps %xmm4,nb133_fiyH2(%esp)
        movaps %xmm4,nb133_fizH2(%esp)
        movaps %xmm4,nb133_fixM(%esp)
        movaps %xmm4,nb133_fiyM(%esp)
        movaps %xmm4,nb133_fizM(%esp)

        movl  nb133_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx                ## jindex[n] 
        movl  4(%eax,%esi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movl  nb133_pos(%ebp),%esi
        movl  nb133_faction(%ebp),%edi
        movl  nb133_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb133_innerjjnr(%esp)        ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb133_ninner(%esp),%ecx
        movl  %ecx,nb133_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb133_innerk(%esp)   ## number of innerloop atoms 
        jge   _nb_kernel133_ia32_sse.nb133_unroll_loop
        jmp   _nb_kernel133_ia32_sse.nb133_odd_inner
_nb_kernel133_ia32_sse.nb133_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb133_innerjjnr(%esp),%edx        ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx             ## eax-edx=jnr1-4 

        addl $16,nb133_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb133_charge(%ebp),%esi    ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4              ## and in xmm4 
        mulps  nb133_iqM(%esp),%xmm3
        mulps  nb133_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps  %xmm3,nb133_qqM(%esp)
        movaps  %xmm4,nb133_qqH(%esp)

        movl nb133_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb133_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb133_ntia(%esp),%edi
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

        movaps %xmm4,nb133_c6(%esp)
        movaps %xmm6,nb133_c12(%esp)

        movl nb133_pos(%ebp),%esi       ## base of pos[] 

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
        movaps nb133_ixO(%esp),%xmm4
        movaps nb133_iyO(%esp),%xmm5
        movaps nb133_izO(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb133_dxO(%esp)
        movaps %xmm5,nb133_dyO(%esp)
        movaps %xmm6,nb133_dzO(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        movaps %xmm4,%xmm7
        ## rsqO in xmm7

        ## move ixH1-izH1 to xmm4-xmm6 
        movaps nb133_ixH1(%esp),%xmm4
        movaps nb133_iyH1(%esp),%xmm5
        movaps nb133_izH1(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb133_dxH1(%esp)
        movaps %xmm5,nb133_dyH1(%esp)
        movaps %xmm6,nb133_dzH1(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm6
        addps %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movaps nb133_ixH2(%esp),%xmm3
        movaps nb133_iyH2(%esp),%xmm4
        movaps nb133_izH2(%esp),%xmm5

        ## calc dr 
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        ## store dr 
        movaps %xmm3,nb133_dxH2(%esp)
        movaps %xmm4,nb133_dyH2(%esp)
        movaps %xmm5,nb133_dzH2(%esp)
        ## square it 
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm4,%xmm5
        addps %xmm3,%xmm5

        ## move ixM-izM to xmm2-xmm4  
        movaps nb133_iyM(%esp),%xmm3
        movaps nb133_izM(%esp),%xmm4
        subps  %xmm1,%xmm3
        subps  %xmm2,%xmm4
        movaps nb133_ixM(%esp),%xmm2
        subps  %xmm0,%xmm2

        ## store dr 
        movaps %xmm2,nb133_dxM(%esp)
        movaps %xmm3,nb133_dyM(%esp)
        movaps %xmm4,nb133_dzM(%esp)
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
        movaps  nb133_three(%esp),%xmm0
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb133_half(%esp),%xmm0
        movaps  %xmm0,nb133_rinvH1(%esp)        ## rinvH1 

        ## rsqH2 - seed to xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb133_three(%esp),%xmm0
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb133_half(%esp),%xmm0
        movaps  %xmm0,nb133_rinvH2(%esp)        ## rinvH2 

        ## rsqM - seed to xmm2 
        rsqrtps %xmm4,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb133_three(%esp),%xmm0
        mulps   %xmm4,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb133_half(%esp),%xmm0
        movaps  %xmm0,nb133_rinvM(%esp)

        ## Do the O LJ-only interaction directly.       
        ## rsqO is in xmm7
        rsqrtps %xmm7,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb133_three(%esp),%xmm4
        mulps   %xmm7,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb133_half(%esp),%xmm4
        movaps  %xmm4,%xmm0
        ## xmm0=rinvO

        mulps %xmm0,%xmm7
        mulps nb133_tsc(%esp),%xmm7   ## rtab

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

        movl nb133_VFtab(%ebp),%esi
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
        mulps  nb133_two(%esp),%xmm7    ## two*Heps2 
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb133_c6(%esp),%xmm4
        mulps  %xmm4,%xmm7       ## fijD 
        mulps  %xmm4,%xmm5       ## Vvdw6 
        mulps  nb133_tsc(%esp),%xmm7
        ## put scalar force on stack Update Vvdwtot directly 
        addps  nb133_Vvdwtot(%esp),%xmm5
        movaps %xmm7,nb133_fstmp(%esp)
        movaps %xmm5,nb133_Vvdwtot(%esp)

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
        mulps  nb133_two(%esp),%xmm7    ## two*Heps2 
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb133_c12(%esp),%xmm4
        mulps  %xmm4,%xmm7 ## fijR 
        mulps  %xmm4,%xmm5 ## Vvdw12 
        mulps  nb133_tsc(%esp),%xmm7
        addps  nb133_fstmp(%esp),%xmm7

        addps  nb133_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb133_Vvdwtot(%esp)

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movd  %mm2,%ecx
        movd  %mm3,%edx

        xorps  %xmm1,%xmm1
        mulps  %xmm0,%xmm7
        subps  %xmm7,%xmm1

        movaps nb133_dxO(%esp),%xmm3
        movaps nb133_dyO(%esp),%xmm4
        movaps nb133_dzO(%esp),%xmm5
        mulps  %xmm1,%xmm3
        mulps  %xmm1,%xmm4
        mulps  %xmm1,%xmm5      ## tx in xmm3-xmm5

        ## update O forces 
        movaps nb133_fixO(%esp),%xmm0
        movaps nb133_fiyO(%esp),%xmm1
        movaps nb133_fizO(%esp),%xmm2
        addps  %xmm3,%xmm0
        addps  %xmm4,%xmm1
        addps  %xmm5,%xmm2
        movaps %xmm0,nb133_fixO(%esp)
        movaps %xmm1,nb133_fiyO(%esp)
        movaps %xmm2,nb133_fizO(%esp)
        ## update j forces with water O 
        movaps %xmm3,nb133_fjx(%esp)
        movaps %xmm4,nb133_fjy(%esp)
        movaps %xmm5,nb133_fjz(%esp)

        ## Do H1 interaction
        movaps  nb133_rinvH1(%esp),%xmm7
        movaps  %xmm7,%xmm4
        mulps   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq
        mulps  nb133_qqH(%esp),%xmm7
        mulps  %xmm7,%xmm4      ## total fs H1 in xmm4 

        addps  nb133_vctot(%esp),%xmm7
        movaps %xmm7,nb133_vctot(%esp)

        movaps nb133_dxH1(%esp),%xmm0
        movaps nb133_dyH1(%esp),%xmm1
        movaps nb133_dzH1(%esp),%xmm2
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2

        ## update H1 forces 
        movaps nb133_fixH1(%esp),%xmm3
        movaps nb133_fiyH1(%esp),%xmm4
        movaps nb133_fizH1(%esp),%xmm7
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm7
        movaps %xmm3,nb133_fixH1(%esp)
        movaps %xmm4,nb133_fiyH1(%esp)
        movaps %xmm7,nb133_fizH1(%esp)
        ## update j forces with water H1 
        addps  nb133_fjx(%esp),%xmm0
        addps  nb133_fjy(%esp),%xmm1
        addps  nb133_fjz(%esp),%xmm2
        movaps %xmm0,nb133_fjx(%esp)
        movaps %xmm1,nb133_fjy(%esp)
        movaps %xmm2,nb133_fjz(%esp)

        ## Done with H1, do H2 interactions
        movaps  nb133_rinvH2(%esp),%xmm7
        movaps  %xmm7,%xmm4
        mulps   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq
        mulps  nb133_qqH(%esp),%xmm7

        mulps  %xmm7,%xmm4      ## total fs H2 in xmm4 

        addps  nb133_vctot(%esp),%xmm7
        movaps %xmm7,nb133_vctot(%esp)

        movaps nb133_dxH2(%esp),%xmm0
        movaps nb133_dyH2(%esp),%xmm1
        movaps nb133_dzH2(%esp),%xmm2
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2

        ## update H2 forces 
        movaps nb133_fixH2(%esp),%xmm3
        movaps nb133_fiyH2(%esp),%xmm4
        movaps nb133_fizH2(%esp),%xmm7
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm7
        movaps %xmm3,nb133_fixH2(%esp)
        movaps %xmm4,nb133_fiyH2(%esp)
        movaps %xmm7,nb133_fizH2(%esp)
        addps nb133_fjx(%esp),%xmm0
        addps nb133_fjy(%esp),%xmm1
        addps nb133_fjz(%esp),%xmm2
        movaps %xmm0,nb133_fjx(%esp)
        movaps %xmm1,nb133_fjy(%esp)
        movaps %xmm2,nb133_fjz(%esp)

        ## Done with H2, do M interactions
        movaps  nb133_rinvM(%esp),%xmm7
        movaps  %xmm7,%xmm4
        mulps   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq
        mulps  nb133_qqM(%esp),%xmm7

        mulps  %xmm7,%xmm4      ## total fs M in xmm4 

        addps  nb133_vctot(%esp),%xmm7
        movaps %xmm7,nb133_vctot(%esp)

        movaps nb133_dxM(%esp),%xmm0
        movaps nb133_dyM(%esp),%xmm1
        movaps nb133_dzM(%esp),%xmm2
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2

        ## update M forces 
        movaps nb133_fixM(%esp),%xmm3
        movaps nb133_fiyM(%esp),%xmm4
        movaps nb133_fizM(%esp),%xmm7
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm7
        movaps %xmm3,nb133_fixM(%esp)
        movaps %xmm4,nb133_fiyM(%esp)
        movaps %xmm7,nb133_fizM(%esp)

        movl nb133_faction(%ebp),%edi
        ## update j forces from stored values
        addps nb133_fjx(%esp),%xmm0
        addps nb133_fjy(%esp),%xmm1
        addps nb133_fjz(%esp),%xmm2

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
        subl $4,nb133_innerk(%esp)
        jl    _nb_kernel133_ia32_sse.nb133_odd_inner
        jmp   _nb_kernel133_ia32_sse.nb133_unroll_loop
_nb_kernel133_ia32_sse.nb133_odd_inner: 
        addl $4,nb133_innerk(%esp)
        jnz   _nb_kernel133_ia32_sse.nb133_odd_loop
        jmp   _nb_kernel133_ia32_sse.nb133_updateouterdata
_nb_kernel133_ia32_sse.nb133_odd_loop: 
        movl  nb133_innerjjnr(%esp),%edx        ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb133_innerjjnr(%esp)

        xorps %xmm4,%xmm4       ## clear reg.
        movss nb133_iqM(%esp),%xmm4
        movl nb133_charge(%ebp),%esi
        movhps nb133_iqH(%esp),%xmm4    ## [qM  0  qH  qH] 
        shufps $41,%xmm4,%xmm4 ## [0 qH qH qM]

        movss (%esi,%eax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb133_qqM(%esp)    ## use dummy qq for storage 

        xorps %xmm6,%xmm6
        movl nb133_type(%ebp),%esi
        movl (%esi,%eax,4),%ebx
        movl nb133_vdwparam(%ebp),%esi
        shll %ebx
        addl nb133_ntia(%esp),%ebx
        movlps (%esi,%ebx,4),%xmm6
        movaps %xmm6,%xmm7
        shufps $252,%xmm6,%xmm6 ## constant 11111100
        shufps $253,%xmm7,%xmm7 ## constant 11111101
        movaps %xmm6,nb133_c6(%esp)
        movaps %xmm7,nb133_c12(%esp)

        movl nb133_pos(%ebp),%esi
        leal (%eax,%eax,2),%eax

        movss nb133_ixO(%esp),%xmm3
        movss nb133_iyO(%esp),%xmm4
        movss nb133_izO(%esp),%xmm5
        movss nb133_ixH1(%esp),%xmm0
        movss nb133_iyH1(%esp),%xmm1
        movss nb133_izH1(%esp),%xmm2
        unpcklps nb133_ixH2(%esp),%xmm3         ## ixO ixH2 - -
        unpcklps nb133_iyH2(%esp),%xmm4         ## iyO iyH2 - -
        unpcklps nb133_izH2(%esp),%xmm5         ## izO izH2 - -
        unpcklps nb133_ixM(%esp),%xmm0          ## ixH1 ixM - -
        unpcklps nb133_iyM(%esp),%xmm1          ## iyH1 iyM - -
        unpcklps nb133_izM(%esp),%xmm2          ## izH1 izM - -
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
        movaps %xmm3,nb133_dxO(%esp)
        movaps %xmm4,nb133_dyO(%esp)
        movaps %xmm5,nb133_dzO(%esp)

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
        movaps nb133_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb133_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv, xmm4=rsq

        mulps %xmm0,%xmm4
        mulps  nb133_tsc(%esp),%xmm4   ## rtab

        cvttps2pi %xmm4,%mm6
        cvtpi2ps %mm6,%xmm6
        subss  %xmm6,%xmm4
        movss %xmm4,%xmm1       ## xmm1=eps 
        movss %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $3,%mm6

        movd %eax,%mm0

        movl nb133_VFtab(%ebp),%esi
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
        mulss  nb133_two(%esp),%xmm7    ## two*Heps2 
        addss  %xmm6,%xmm7
        addss  %xmm5,%xmm7 ## xmm7=FF 
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 

        movss nb133_c6(%esp),%xmm4
        mulss  %xmm4,%xmm7       ## fijD 
        mulss  %xmm4,%xmm5       ## Vvdw6 
        mulss  nb133_tsc(%esp),%xmm7
        ## put scalar force on stack Update Vvdwtot directly 
        addss  nb133_Vvdwtot(%esp),%xmm5
        movss %xmm7,nb133_fstmp(%esp)
        movss %xmm5,nb133_Vvdwtot(%esp)

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
        mulss  nb133_two(%esp),%xmm7    ## two*Heps2 
        addss  %xmm6,%xmm7
        addss  %xmm5,%xmm7 ## xmm7=FF 
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 

        movss nb133_c12(%esp),%xmm4
        mulss  %xmm4,%xmm7 ## fijR 
        mulss  %xmm4,%xmm5 ## Vvdw12 
        mulss  nb133_tsc(%esp),%xmm7
        addss  nb133_fstmp(%esp),%xmm7
        movss %xmm7,nb133_fstmp(%esp)
        addss  nb133_Vvdwtot(%esp),%xmm5
        movss %xmm5,nb133_Vvdwtot(%esp)

        movd %mm0,%eax

        movaps %xmm0,%xmm4
        mulps  nb133_qqM(%esp),%xmm4
        movaps %xmm4,%xmm2
        mulps  %xmm0,%xmm4
        subss  nb133_fstmp(%esp),%xmm4
        mulps  %xmm0,%xmm4

        addps  nb133_vctot(%esp),%xmm2
        movaps %xmm2,nb133_vctot(%esp)

        movaps nb133_dxO(%esp),%xmm0
        movaps nb133_dyO(%esp),%xmm1
        movaps nb133_dzO(%esp),%xmm2

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2 ## xmm0-xmm2 now contains tx-tz (partial force)

        movss  nb133_fixO(%esp),%xmm3
        movss  nb133_fiyO(%esp),%xmm4
        movss  nb133_fizO(%esp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,nb133_fixO(%esp)
        movss  %xmm4,nb133_fiyO(%esp)
        movss  %xmm5,nb133_fizO(%esp)   ## updated the O force now do the H's

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        shufps $0x39,%xmm3,%xmm3 ## shift right 
        shufps $0x39,%xmm4,%xmm4
        shufps $0x39,%xmm5,%xmm5
        addss  nb133_fixH1(%esp),%xmm3
        addss  nb133_fiyH1(%esp),%xmm4
        addss  nb133_fizH1(%esp),%xmm5
        movss  %xmm3,nb133_fixH1(%esp)
        movss  %xmm4,nb133_fiyH1(%esp)
        movss  %xmm5,nb133_fizH1(%esp)          ## updated the H1 force 

        shufps $0x39,%xmm3,%xmm3
        shufps $0x39,%xmm4,%xmm4
        shufps $0x39,%xmm5,%xmm5
        addss  nb133_fixH2(%esp),%xmm3
        addss  nb133_fiyH2(%esp),%xmm4
        addss  nb133_fizH2(%esp),%xmm5
        movss  %xmm3,nb133_fixH2(%esp)
        movss  %xmm4,nb133_fiyH2(%esp)
        movss  %xmm5,nb133_fizH2(%esp)          ## updated the H2 force 

        movl nb133_faction(%ebp),%edi
        shufps $0x39,%xmm3,%xmm3
        shufps $0x39,%xmm4,%xmm4
        shufps $0x39,%xmm5,%xmm5
        addss  nb133_fixM(%esp),%xmm3
        addss  nb133_fiyM(%esp),%xmm4
        addss  nb133_fizM(%esp),%xmm5
        movss  %xmm3,nb133_fixM(%esp)
        movss  %xmm4,nb133_fiyM(%esp)
        movss  %xmm5,nb133_fizM(%esp)   ## updated the M force 

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

        decl nb133_innerk(%esp)
        jz    _nb_kernel133_ia32_sse.nb133_updateouterdata
        jmp   _nb_kernel133_ia32_sse.nb133_odd_loop
_nb_kernel133_ia32_sse.nb133_updateouterdata: 
        movl  nb133_ii3(%esp),%ecx
        movl  nb133_faction(%ebp),%edi
        movl  nb133_fshift(%ebp),%esi
        movl  nb133_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movaps nb133_fixO(%esp),%xmm0
        movaps nb133_fiyO(%esp),%xmm1
        movaps nb133_fizO(%esp),%xmm2

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
        movaps nb133_fixH1(%esp),%xmm0
        movaps nb133_fiyH1(%esp),%xmm1
        movaps nb133_fizH1(%esp),%xmm2

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
        movaps nb133_fixH2(%esp),%xmm0
        movaps nb133_fiyH2(%esp),%xmm1
        movaps nb133_fizH2(%esp),%xmm2

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
        movaps nb133_fixM(%esp),%xmm0
        movaps nb133_fiyM(%esp),%xmm1
        movaps nb133_fizM(%esp),%xmm2

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
        movl nb133_n(%esp),%esi
        ## get group index for i particle 
        movl  nb133_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb133_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb133_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb133_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb133_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb133_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel133_ia32_sse.nb133_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb133_n(%esp)
        jmp _nb_kernel133_ia32_sse.nb133_outer
_nb_kernel133_ia32_sse.nb133_outerend: 
        ## check if more outer neighborlists remain
        movl  nb133_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel133_ia32_sse.nb133_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel133_ia32_sse.nb133_threadloop
_nb_kernel133_ia32_sse.nb133_end: 
        emms

        movl nb133_nouter(%esp),%eax
        movl nb133_ninner(%esp),%ebx
        movl nb133_outeriter(%ebp),%ecx
        movl nb133_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb133_salign(%esp),%eax
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


.globl nb_kernel133nf_ia32_sse
.globl _nb_kernel133nf_ia32_sse
nb_kernel133nf_ia32_sse:        
_nb_kernel133nf_ia32_sse:       
.set nb133nf_p_nri, 8
.set nb133nf_iinr, 12
.set nb133nf_jindex, 16
.set nb133nf_jjnr, 20
.set nb133nf_shift, 24
.set nb133nf_shiftvec, 28
.set nb133nf_fshift, 32
.set nb133nf_gid, 36
.set nb133nf_pos, 40
.set nb133nf_faction, 44
.set nb133nf_charge, 48
.set nb133nf_p_facel, 52
.set nb133nf_argkrf, 56
.set nb133nf_argcrf, 60
.set nb133nf_Vc, 64
.set nb133nf_type, 68
.set nb133nf_p_ntype, 72
.set nb133nf_vdwparam, 76
.set nb133nf_Vvdw, 80
.set nb133nf_p_tabscale, 84
.set nb133nf_VFtab, 88
.set nb133nf_invsqrta, 92
.set nb133nf_dvda, 96
.set nb133nf_p_gbtabscale, 100
.set nb133nf_GBtab, 104
.set nb133nf_p_nthreads, 108
.set nb133nf_count, 112
.set nb133nf_mtx, 116
.set nb133nf_outeriter, 120
.set nb133nf_inneriter, 124
.set nb133nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb133nf_ixO, 0
.set nb133nf_iyO, 16
.set nb133nf_izO, 32
.set nb133nf_ixH1, 48
.set nb133nf_iyH1, 64
.set nb133nf_izH1, 80
.set nb133nf_ixH2, 96
.set nb133nf_iyH2, 112
.set nb133nf_izH2, 128
.set nb133nf_ixM, 144
.set nb133nf_iyM, 160
.set nb133nf_izM, 176
.set nb133nf_iqM, 192
.set nb133nf_iqH, 208
.set nb133nf_qqH, 224
.set nb133nf_rinvH1, 240
.set nb133nf_rinvH2, 256
.set nb133nf_rinvM, 272
.set nb133nf_c6, 288
.set nb133nf_c12, 304
.set nb133nf_tsc, 320
.set nb133nf_vctot, 416
.set nb133nf_Vvdwtot, 432
.set nb133nf_half, 448
.set nb133nf_three, 464
.set nb133nf_qqM, 480
.set nb133nf_is3, 496
.set nb133nf_ii3, 500
.set nb133nf_ntia, 504
.set nb133nf_innerjjnr, 508
.set nb133nf_innerk, 512
.set nb133nf_n, 516
.set nb133nf_nn1, 520
.set nb133nf_nri, 524
.set nb133nf_nouter, 528
.set nb133nf_ninner, 532
.set nb133nf_salign, 536
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
        movl %eax,nb133nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb133nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb133nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb133nf_nouter(%esp)
        movl %eax,nb133nf_ninner(%esp)

        movl nb133nf_p_tabscale(%ebp),%eax
        movss (%eax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb133nf_tsc(%esp)


        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb133nf_half(%esp)
        movss nb133nf_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb133nf_half(%esp)
        movaps %xmm3,nb133nf_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb133nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb133nf_charge(%ebp),%edx
        movss 4(%edx,%ebx,4),%xmm4
        movss 12(%edx,%ebx,4),%xmm3
        movl nb133nf_p_facel(%ebp),%esi
        movss (%esi),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb133nf_iqM(%esp)
        movaps %xmm4,nb133nf_iqH(%esp)

        movl  nb133nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb133nf_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb133nf_ntia(%esp)
_nb_kernel133nf_ia32_sse.nb133nf_threadloop: 
        movl  nb133nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel133nf_ia32_sse.nb133nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel133nf_ia32_sse.nb133nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb133nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb133nf_n(%esp)
        movl %ebx,nb133nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel133nf_ia32_sse.nb133nf_outerstart
        jmp _nb_kernel133nf_ia32_sse.nb133nf_end

_nb_kernel133nf_ia32_sse.nb133nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb133nf_nouter(%esp),%ebx
        movl %ebx,nb133nf_nouter(%esp)

_nb_kernel133nf_ia32_sse.nb133nf_outer: 
        movl  nb133nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx        ## ebx=3*is 
        movl  %ebx,nb133nf_is3(%esp)            ## store is3 

        movl  nb133nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb133nf_iinr(%ebp),%ecx           ## ecx = pointer into iinr[]    
        movl  (%ecx,%esi,4),%ebx                ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        movaps %xmm0,%xmm6
        movaps %xmm1,%xmm7

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb133nf_pos(%ebp),%eax    ## eax = base of pos[]  
        movl  %ebx,nb133nf_ii3(%esp)

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
        movaps %xmm3,nb133nf_ixO(%esp)
        movaps %xmm4,nb133nf_iyO(%esp)
        movaps %xmm5,nb133nf_izO(%esp)
        movaps %xmm6,nb133nf_ixH1(%esp)
        movaps %xmm7,nb133nf_iyH1(%esp)

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
        movaps %xmm6,nb133nf_izH1(%esp)
        movaps %xmm0,nb133nf_ixH2(%esp)
        movaps %xmm1,nb133nf_iyH2(%esp)
        movaps %xmm2,nb133nf_izH2(%esp)
        movaps %xmm3,nb133nf_ixM(%esp)
        movaps %xmm4,nb133nf_iyM(%esp)
        movaps %xmm5,nb133nf_izM(%esp)

        ## clear vctot 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb133nf_vctot(%esp)
        movaps %xmm4,nb133nf_Vvdwtot(%esp)

        movl  nb133nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx                ## jindex[n] 
        movl  4(%eax,%esi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movl  nb133nf_pos(%ebp),%esi
        movl  nb133nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb133nf_innerjjnr(%esp)      ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb133nf_ninner(%esp),%ecx
        movl  %ecx,nb133nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb133nf_innerk(%esp)         ## number of innerloop atoms 
        jge   _nb_kernel133nf_ia32_sse.nb133nf_unroll_loop
        jmp   _nb_kernel133nf_ia32_sse.nb133nf_odd_inner
_nb_kernel133nf_ia32_sse.nb133nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb133nf_innerjjnr(%esp),%edx      ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx             ## eax-edx=jnr1-4 

        addl $16,nb133nf_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb133nf_charge(%ebp),%esi  ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4              ## and in xmm4 
        mulps  nb133nf_iqM(%esp),%xmm3
        mulps  nb133nf_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps  %xmm3,nb133nf_qqM(%esp)
        movaps  %xmm4,nb133nf_qqH(%esp)

        movl nb133nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb133nf_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb133nf_ntia(%esp),%edi
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

        movaps %xmm4,nb133nf_c6(%esp)
        movaps %xmm6,nb133nf_c12(%esp)

        movl nb133nf_pos(%ebp),%esi     ## base of pos[] 

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
        movaps nb133nf_ixO(%esp),%xmm4
        movaps nb133nf_iyO(%esp),%xmm5
        movaps nb133nf_izO(%esp),%xmm6

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
        movaps nb133nf_ixH1(%esp),%xmm4
        movaps nb133nf_iyH1(%esp),%xmm5
        movaps nb133nf_izH1(%esp),%xmm6

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
        movaps nb133nf_ixH2(%esp),%xmm3
        movaps nb133nf_iyH2(%esp),%xmm4
        movaps nb133nf_izH2(%esp),%xmm5

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
        movaps nb133nf_iyM(%esp),%xmm3
        movaps nb133nf_izM(%esp),%xmm4
        subps  %xmm1,%xmm3
        subps  %xmm2,%xmm4
        movaps nb133nf_ixM(%esp),%xmm2
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
        movaps  nb133nf_three(%esp),%xmm0
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb133nf_half(%esp),%xmm0
        movaps  %xmm0,nb133nf_rinvH1(%esp)      ## rinvH1 

        ## rsqH2 - seed to xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb133nf_three(%esp),%xmm0
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb133nf_half(%esp),%xmm0
        movaps  %xmm0,nb133nf_rinvH2(%esp)      ## rinvH2 

        ## rsqM - seed to xmm2 
        rsqrtps %xmm4,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb133nf_three(%esp),%xmm0
        mulps   %xmm4,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb133nf_half(%esp),%xmm0
        movaps  %xmm0,nb133nf_rinvM(%esp)

        ## Do the O LJ-only interaction directly.       
        ## rsqO is in xmm7
        rsqrtps %xmm7,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb133nf_three(%esp),%xmm4
        mulps   %xmm7,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb133nf_half(%esp),%xmm4
        movaps  %xmm4,%xmm0
        ## xmm0=rinvO

        mulps %xmm0,%xmm7
        mulps nb133nf_tsc(%esp),%xmm7   ## rtab

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

        movl nb133nf_VFtab(%ebp),%esi
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

        movaps nb133nf_c6(%esp),%xmm4
        mulps  %xmm4,%xmm5       ## Vvdw6 

        addps  nb133nf_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb133nf_Vvdwtot(%esp)

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

        movaps nb133nf_c12(%esp),%xmm4
        mulps  %xmm4,%xmm5 ## Vvdw12 

        addps  nb133nf_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb133nf_Vvdwtot(%esp)

        ## Do H1-H2-M interactions      
        movaps  nb133nf_rinvH1(%esp),%xmm7
        addps   nb133nf_rinvH2(%esp),%xmm7
        movaps  nb133nf_rinvM(%esp),%xmm6

        mulps   nb133nf_qqH(%esp),%xmm7
        mulps   nb133nf_qqM(%esp),%xmm6
        addps   %xmm6,%xmm7

        addps  nb133nf_vctot(%esp),%xmm7
        movaps %xmm7,nb133nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $4,nb133nf_innerk(%esp)
        jl    _nb_kernel133nf_ia32_sse.nb133nf_odd_inner
        jmp   _nb_kernel133nf_ia32_sse.nb133nf_unroll_loop
_nb_kernel133nf_ia32_sse.nb133nf_odd_inner: 
        addl $4,nb133nf_innerk(%esp)
        jnz   _nb_kernel133nf_ia32_sse.nb133nf_odd_loop
        jmp   _nb_kernel133nf_ia32_sse.nb133nf_updateouterdata
_nb_kernel133nf_ia32_sse.nb133nf_odd_loop: 
        movl  nb133nf_innerjjnr(%esp),%edx      ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb133nf_innerjjnr(%esp)

        xorps %xmm4,%xmm4       ## clear reg.
        movss nb133nf_iqM(%esp),%xmm4
        movl nb133nf_charge(%ebp),%esi
        movhps nb133nf_iqH(%esp),%xmm4    ## [qM  0  qH  qH] 
        shufps $41,%xmm4,%xmm4 ## [0 qH qH qM]

        movss (%esi,%eax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb133nf_qqM(%esp)          ## use dummy qq for storage 

        xorps %xmm6,%xmm6
        movl nb133nf_type(%ebp),%esi
        movl (%esi,%eax,4),%ebx
        movl nb133nf_vdwparam(%ebp),%esi
        shll %ebx
        addl nb133nf_ntia(%esp),%ebx
        movlps (%esi,%ebx,4),%xmm6
        movaps %xmm6,%xmm7
        shufps $252,%xmm6,%xmm6 ## constant 11111100
        shufps $253,%xmm7,%xmm7 ## constant 11111101
        movaps %xmm6,nb133nf_c6(%esp)
        movaps %xmm7,nb133nf_c12(%esp)

        movl nb133nf_pos(%ebp),%esi
        leal (%eax,%eax,2),%eax

        movss nb133nf_ixO(%esp),%xmm3
        movss nb133nf_iyO(%esp),%xmm4
        movss nb133nf_izO(%esp),%xmm5
        movss nb133nf_ixH1(%esp),%xmm0
        movss nb133nf_iyH1(%esp),%xmm1
        movss nb133nf_izH1(%esp),%xmm2
        unpcklps nb133nf_ixH2(%esp),%xmm3       ## ixO ixH2 - -
        unpcklps nb133nf_iyH2(%esp),%xmm4       ## iyO iyH2 - -
        unpcklps nb133nf_izH2(%esp),%xmm5       ## izO izH2 - -
        unpcklps nb133nf_ixM(%esp),%xmm0        ## ixH1 ixM - -
        unpcklps nb133nf_iyM(%esp),%xmm1        ## iyH1 iyM - -
        unpcklps nb133nf_izM(%esp),%xmm2        ## izH1 izM - -
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
        movaps nb133nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb133nf_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv, xmm4=rsq

        mulps %xmm0,%xmm4
        mulps  nb133nf_tsc(%esp),%xmm4   ## rtab

        cvttps2pi %xmm4,%mm6
        cvtpi2ps %mm6,%xmm6
        subss  %xmm6,%xmm4
        movss %xmm4,%xmm1       ## xmm1=eps 
        movss %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $3,%mm6

        movl nb133nf_VFtab(%ebp),%esi
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

        movss nb133nf_c6(%esp),%xmm4
        mulss  %xmm4,%xmm5       ## Vvdw6 

        ## put scalar force on stack Update Vvdwtot directly 
        addss  nb133nf_Vvdwtot(%esp),%xmm5
        movss %xmm5,nb133nf_Vvdwtot(%esp)

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

        movss nb133nf_c12(%esp),%xmm4
        mulss  %xmm4,%xmm5 ## Vvdw12 

        addss  nb133nf_Vvdwtot(%esp),%xmm5
        movss %xmm5,nb133nf_Vvdwtot(%esp)

        mulps  nb133nf_qqM(%esp),%xmm0          ## xmm0=vcoul 

        addps  nb133nf_vctot(%esp),%xmm0
        movaps %xmm0,nb133nf_vctot(%esp)

        decl nb133nf_innerk(%esp)
        jz    _nb_kernel133nf_ia32_sse.nb133nf_updateouterdata
        jmp   _nb_kernel133nf_ia32_sse.nb133nf_odd_loop
_nb_kernel133nf_ia32_sse.nb133nf_updateouterdata: 
        ## get n from stack
        movl nb133nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb133nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb133nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb133nf_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb133nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb133nf_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb133nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel133nf_ia32_sse.nb133nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb133nf_n(%esp)
        jmp _nb_kernel133nf_ia32_sse.nb133nf_outer
_nb_kernel133nf_ia32_sse.nb133nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb133nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel133nf_ia32_sse.nb133nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel133nf_ia32_sse.nb133nf_threadloop
_nb_kernel133nf_ia32_sse.nb133nf_end: 
        emms

        movl nb133nf_nouter(%esp),%eax
        movl nb133nf_ninner(%esp),%ebx
        movl nb133nf_outeriter(%ebp),%ecx
        movl nb133nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb133nf_salign(%esp),%eax
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


