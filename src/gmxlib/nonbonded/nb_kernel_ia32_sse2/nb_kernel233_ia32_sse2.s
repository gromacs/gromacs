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



.globl nb_kernel233_ia32_sse2
.globl _nb_kernel233_ia32_sse2
nb_kernel233_ia32_sse2: 
_nb_kernel233_ia32_sse2:        
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
        ## bottom of stack is cache-aligned for sse2 use 
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
.set nb233_iqH, 192
.set nb233_iqM, 208
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
.set nb233_qqH, 416
.set nb233_qqM, 432
.set nb233_c6, 448
.set nb233_c12, 464
.set nb233_tsc, 480
.set nb233_fstmp, 496
.set nb233_vctot, 512
.set nb233_Vvdwtot, 528
.set nb233_fixO, 544
.set nb233_fiyO, 560
.set nb233_fizO, 576
.set nb233_fixH1, 592
.set nb233_fiyH1, 608
.set nb233_fizH1, 624
.set nb233_fixH2, 640
.set nb233_fiyH2, 656
.set nb233_fizH2, 672
.set nb233_fixM, 688
.set nb233_fiyM, 704
.set nb233_fizM, 720
.set nb233_fjx, 736
.set nb233_fjy, 752
.set nb233_fjz, 768
.set nb233_half, 784
.set nb233_three, 800
.set nb233_two, 816
.set nb233_rinvH1, 832
.set nb233_rinvH2, 848
.set nb233_rinvM, 864
.set nb233_krsqH1, 880
.set nb233_krsqH2, 896
.set nb233_krsqM, 912
.set nb233_krf, 928
.set nb233_crf, 944
.set nb233_rsqO, 960
.set nb233_is3, 976
.set nb233_ii3, 980
.set nb233_ntia, 984
.set nb233_innerjjnr, 988
.set nb233_innerk, 992
.set nb233_n, 996
.set nb233_nn1, 1000
.set nb233_nri, 1004
.set nb233_nouter, 1008
.set nb233_ninner, 1012
.set nb233_salign, 1016
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $1020,%esp         ## local stack space 
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
        movsd (%eax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb233_tsc(%esp)

        movl nb233_argkrf(%ebp),%esi
        movl nb233_argcrf(%ebp),%edi
        movsd (%esi),%xmm5
        movsd (%edi),%xmm6
        shufpd $0,%xmm5,%xmm5
        shufpd $0,%xmm6,%xmm6
        movapd %xmm5,nb233_krf(%esp)
        movapd %xmm6,nb233_crf(%esp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb233_half(%esp)
        movl %ebx,nb233_half+4(%esp)
        movsd nb233_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb233_half(%esp)
        movapd %xmm2,nb233_two(%esp)
        movapd %xmm3,nb233_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb233_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb233_charge(%ebp),%edx
        movsd 8(%edx,%ebx,8),%xmm3
        movsd 24(%edx,%ebx,8),%xmm4
        movl nb233_p_facel(%ebp),%esi
        movsd (%esi),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb233_iqH(%esp)
        movapd %xmm4,nb233_iqM(%esp)

        movl  nb233_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb233_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb233_ntia(%esp)
_nb_kernel233_ia32_sse2.nb233_threadloop: 
        movl  nb233_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel233_ia32_sse2.nb233_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel233_ia32_sse2.nb233_spinlock

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
        jg  _nb_kernel233_ia32_sse2.nb233_outerstart
        jmp _nb_kernel233_ia32_sse2.nb233_end

_nb_kernel233_ia32_sse2.nb233_outerstart: 
        ## ebx contains number of outer iterations
        addl nb233_nouter(%esp),%ebx
        movl %ebx,nb233_nouter(%esp)

_nb_kernel233_ia32_sse2.nb233_outer: 
        movl  nb233_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb233_is3(%esp)      ## store is3 

        movl  nb233_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb233_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        movapd %xmm0,%xmm6
        movapd %xmm1,%xmm7

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb233_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb233_ii3(%esp)

        addsd (%eax,%ebx,8),%xmm3       ## ox
        addsd 8(%eax,%ebx,8),%xmm4      ## oy
        addsd 16(%eax,%ebx,8),%xmm5     ## oz   
        addsd 24(%eax,%ebx,8),%xmm6     ## h1x
        addsd 32(%eax,%ebx,8),%xmm7     ## h1y
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        shufpd $0,%xmm6,%xmm6
        shufpd $0,%xmm7,%xmm7
        movapd %xmm3,nb233_ixO(%esp)
        movapd %xmm4,nb233_iyO(%esp)
        movapd %xmm5,nb233_izO(%esp)
        movapd %xmm6,nb233_ixH1(%esp)
        movapd %xmm7,nb233_iyH1(%esp)

        movsd %xmm2,%xmm6
        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 40(%eax,%ebx,8),%xmm6    ## h1z
        addsd 48(%eax,%ebx,8),%xmm0    ## h2x
        addsd 56(%eax,%ebx,8),%xmm1    ## h2y
        addsd 64(%eax,%ebx,8),%xmm2    ## h2z
        addsd 72(%eax,%ebx,8),%xmm3    ## mx
        addsd 80(%eax,%ebx,8),%xmm4    ## my
        addsd 88(%eax,%ebx,8),%xmm5    ## mz

        shufpd $0,%xmm6,%xmm6
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm6,nb233_izH1(%esp)
        movapd %xmm0,nb233_ixH2(%esp)
        movapd %xmm1,nb233_iyH2(%esp)
        movapd %xmm2,nb233_izH2(%esp)
        movapd %xmm3,nb233_ixM(%esp)
        movapd %xmm4,nb233_iyM(%esp)
        movapd %xmm5,nb233_izM(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb233_vctot(%esp)
        movapd %xmm4,nb233_Vvdwtot(%esp)
        movapd %xmm4,nb233_fixO(%esp)
        movapd %xmm4,nb233_fiyO(%esp)
        movapd %xmm4,nb233_fizO(%esp)
        movapd %xmm4,nb233_fixH1(%esp)
        movapd %xmm4,nb233_fiyH1(%esp)
        movapd %xmm4,nb233_fizH1(%esp)
        movapd %xmm4,nb233_fixH2(%esp)
        movapd %xmm4,nb233_fiyH2(%esp)
        movapd %xmm4,nb233_fizH2(%esp)
        movapd %xmm4,nb233_fixM(%esp)
        movapd %xmm4,nb233_fiyM(%esp)
        movapd %xmm4,nb233_fizM(%esp)

        movl  nb233_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb233_pos(%ebp),%esi
        movl  nb233_faction(%ebp),%edi
        movl  nb233_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb233_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb233_ninner(%esp),%ecx
        movl  %ecx,nb233_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb233_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel233_ia32_sse2.nb233_unroll_loop
        jmp   _nb_kernel233_ia32_sse2.nb233_checksingle
_nb_kernel233_ia32_sse2.nb233_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb233_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb233_innerjjnr(%esp)                   ## advance pointer (unrolled 2) 

        movl nb233_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb233_iqM(%esp),%xmm3
        mulpd  nb233_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movapd  %xmm3,nb233_qqM(%esp)
        movapd  %xmm4,nb233_qqH(%esp)

        movl nb233_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb233_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb233_ntia(%esp),%edi
        addl %edi,%eax
        addl %edi,%ebx

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movlpd (%esi,%ebx,8),%xmm7      ## c6b
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 
        movhpd 8(%esi,%ebx,8),%xmm7     ## c6b c12b 
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movapd %xmm4,nb233_c6(%esp)
        movapd %xmm6,nb233_c12(%esp)

        movl nb233_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move two coordinates to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2
        movhpd (%esi,%ebx,8),%xmm0
        movhpd 8(%esi,%ebx,8),%xmm1
        movhpd 16(%esi,%ebx,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb233_ixO(%esp),%xmm4
        movapd nb233_iyO(%esp),%xmm5
        movapd nb233_izO(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb233_dxO(%esp)
        movapd %xmm5,nb233_dyO(%esp)
        movapd %xmm6,nb233_dzO(%esp)
        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb233_ixH1(%esp),%xmm4
        movapd nb233_iyH1(%esp),%xmm5
        movapd nb233_izH1(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb233_dxH1(%esp)
        movapd %xmm5,nb233_dyH1(%esp)
        movapd %xmm6,nb233_dzH1(%esp)
        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm6
        addpd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb233_ixH2(%esp),%xmm3
        movapd nb233_iyH2(%esp),%xmm4
        movapd nb233_izH2(%esp),%xmm5

        ## calc dr 
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5

        ## store dr 
        movapd %xmm3,nb233_dxH2(%esp)
        movapd %xmm4,nb233_dyH2(%esp)
        movapd %xmm5,nb233_dzH2(%esp)
        ## square it 
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        addpd %xmm4,%xmm5
        addpd %xmm3,%xmm5

        ## move ixM-izM to xmm2-xmm4  
        movapd nb233_iyM(%esp),%xmm3
        movapd nb233_izM(%esp),%xmm4
        subpd  %xmm1,%xmm3
        subpd  %xmm2,%xmm4
        movapd nb233_ixM(%esp),%xmm2
        subpd  %xmm0,%xmm2

        ## store dr 
        movapd %xmm2,nb233_dxM(%esp)
        movapd %xmm3,nb233_dyM(%esp)
        movapd %xmm4,nb233_dzM(%esp)
        ## square it 
        mulpd %xmm2,%xmm2
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        addpd %xmm3,%xmm4
        addpd %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 
        movapd %xmm7,nb233_rsqO(%esp)

        ## calculate krsq
        movapd nb233_krf(%esp),%xmm0
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2
        mulpd %xmm4,%xmm0
        mulpd %xmm5,%xmm1
        mulpd %xmm6,%xmm2
        movapd %xmm0,nb233_krsqM(%esp)
        movapd %xmm1,nb233_krsqH2(%esp)
        movapd %xmm2,nb233_krsqH1(%esp)

        ## start with rsqH1 - put seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb233_three(%esp),%xmm1
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulpd   nb233_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulpd %xmm1,%xmm1       ## lu*lu 
        mulpd %xmm1,%xmm6       ## rsq*lu*lu 
        movapd nb233_three(%esp),%xmm1
        subpd %xmm6,%xmm1       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulpd nb233_half(%esp),%xmm1   ## rinv 
        movapd  %xmm1,nb233_rinvH1(%esp)

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb233_three(%esp),%xmm1
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulpd   nb233_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulpd %xmm1,%xmm1       ## lu*lu 
        mulpd %xmm1,%xmm5       ## rsq*lu*lu 
        movapd nb233_three(%esp),%xmm1
        subpd %xmm5,%xmm1       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulpd nb233_half(%esp),%xmm1   ## rinv 
        movapd  %xmm1,nb233_rinvH2(%esp)

        ## rsqM - seed in xmm2 
        cvtpd2ps %xmm4,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb233_three(%esp),%xmm1
        mulpd   %xmm4,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulpd   nb233_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulpd %xmm1,%xmm1       ## lu*lu 
        mulpd %xmm1,%xmm4       ## rsq*lu*lu 
        movapd nb233_three(%esp),%xmm1
        subpd %xmm4,%xmm1       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulpd nb233_half(%esp),%xmm1   ## rinv 
        movapd  %xmm1,nb233_rinvM(%esp)


        ## rsqO - put seed in xmm2 
        cvtpd2ps %xmm7,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb233_three(%esp),%xmm4
        mulpd   %xmm7,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb233_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb233_three(%esp),%xmm4
        subpd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb233_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvO in xmm7 



        movapd nb233_rsqO(%esp),%xmm4
        movapd %xmm7,%xmm0
        ## LJ table interaction.
        mulpd %xmm7,%xmm4       ## xmm4=r 
        mulpd nb233_tsc(%esp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $3,%mm6           ## idx *= 8 

        movd %eax,%mm0
        movd %ebx,%mm1

        movl nb233_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx

        ## dispersion 
        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2 
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## dispersion table ready, in xmm4-xmm7         
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb233_two(%esp),%xmm7    ## two*Heps2 
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb233_c6(%esp),%xmm4
        mulpd  %xmm4,%xmm7       ## fijD 
        mulpd  %xmm4,%xmm5       ## Vvdw6 

        ## put scalar force on stack Update Vvdwtot directly 
        addpd  nb233_Vvdwtot(%esp),%xmm5
        xorpd  %xmm3,%xmm3
        mulpd  nb233_tsc(%esp),%xmm7
        subpd  %xmm7,%xmm3
        movapd %xmm3,nb233_fstmp(%esp)
        movapd %xmm5,nb233_Vvdwtot(%esp)

        ## repulsion 
        movlpd 32(%esi,%eax,8),%xmm4    ## Y1   
        movlpd 32(%esi,%ebx,8),%xmm3    ## Y2 
        movhpd 40(%esi,%eax,8),%xmm4    ## Y1 F1        
        movhpd 40(%esi,%ebx,8),%xmm3    ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 48(%esi,%eax,8),%xmm6    ## G1
        movlpd 48(%esi,%ebx,8),%xmm3    ## G2
        movhpd 56(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 56(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 

        ## table ready, in xmm4-xmm7    
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb233_two(%esp),%xmm7    ## two*Heps2 
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb233_c12(%esp),%xmm4
        mulpd  %xmm4,%xmm7
        mulpd  %xmm4,%xmm5

        addpd  nb233_Vvdwtot(%esp),%xmm5
        movapd nb233_fstmp(%esp),%xmm3
        mulpd  nb233_tsc(%esp),%xmm7
        subpd  %xmm7,%xmm3
        movapd %xmm5,nb233_Vvdwtot(%esp)

        mulpd  %xmm0,%xmm3


        movapd nb233_dxO(%esp),%xmm0
        movapd nb233_dyO(%esp),%xmm1
        movapd nb233_dzO(%esp),%xmm2

        movd %mm0,%eax
        movd %mm1,%ebx

        movl   nb233_faction(%ebp),%edi
        mulpd  %xmm3,%xmm0
        mulpd  %xmm3,%xmm1
        mulpd  %xmm3,%xmm2

        ## update O forces 
        movapd nb233_fixO(%esp),%xmm3
        movapd nb233_fiyO(%esp),%xmm4
        movapd nb233_fizO(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb233_fixO(%esp)
        movapd %xmm4,nb233_fiyO(%esp)
        movapd %xmm7,nb233_fizO(%esp)
        ## update j forces with water O 
        movapd %xmm0,nb233_fjx(%esp)
        movapd %xmm1,nb233_fjy(%esp)
        movapd %xmm2,nb233_fjz(%esp)

        ## H1 interactions 
        movapd  nb233_rinvH1(%esp),%xmm6
        movapd  %xmm6,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        movapd  %xmm6,%xmm7
        movapd  nb233_krsqH1(%esp),%xmm0
        addpd   %xmm0,%xmm6     ## xmm6=rinv+ krsq 
        mulpd   nb233_two(%esp),%xmm0
        subpd   nb233_crf(%esp),%xmm6
        subpd   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulpd   nb233_qqH(%esp),%xmm6   ## vcoul 
        mulpd   nb233_qqH(%esp),%xmm7
        mulpd  %xmm7,%xmm4              ## total fsH1 in xmm4 
        addpd  nb233_vctot(%esp),%xmm6

        movapd nb233_dxH1(%esp),%xmm0
        movapd nb233_dyH1(%esp),%xmm1
        movapd nb233_dzH1(%esp),%xmm2
        movapd %xmm6,nb233_vctot(%esp)
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb233_fixH1(%esp),%xmm3
        movapd nb233_fiyH1(%esp),%xmm4
        movapd nb233_fizH1(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb233_fixH1(%esp)
        movapd %xmm4,nb233_fiyH1(%esp)
        movapd %xmm7,nb233_fizH1(%esp)
        ## update j forces with water H1 
        addpd  nb233_fjx(%esp),%xmm0
        addpd  nb233_fjy(%esp),%xmm1
        addpd  nb233_fjz(%esp),%xmm2
        movapd %xmm0,nb233_fjx(%esp)
        movapd %xmm1,nb233_fjy(%esp)
        movapd %xmm2,nb233_fjz(%esp)

        ## H2 interactions 
        movapd  nb233_rinvH2(%esp),%xmm5
        movapd  %xmm5,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        movapd  %xmm5,%xmm7
        movapd  nb233_krsqH2(%esp),%xmm0
        addpd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        mulpd   nb233_two(%esp),%xmm0
        subpd   nb233_crf(%esp),%xmm5
        subpd   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulpd   nb233_qqH(%esp),%xmm5   ## vcoul 
        mulpd   nb233_qqH(%esp),%xmm7
        mulpd  %xmm7,%xmm4              ## total fsH2 in xmm4 

        addpd  nb233_vctot(%esp),%xmm5

        movapd nb233_dxH2(%esp),%xmm0
        movapd nb233_dyH2(%esp),%xmm1
        movapd nb233_dzH2(%esp),%xmm2
        movapd %xmm5,nb233_vctot(%esp)
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        ## update H2 forces 
        movapd nb233_fixH2(%esp),%xmm3
        movapd nb233_fiyH2(%esp),%xmm4
        movapd nb233_fizH2(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb233_fixH2(%esp)
        movapd %xmm4,nb233_fiyH2(%esp)
        movapd %xmm7,nb233_fizH2(%esp)
        ## update j forces with water H2
        addpd  nb233_fjx(%esp),%xmm0
        addpd  nb233_fjy(%esp),%xmm1
        addpd  nb233_fjz(%esp),%xmm2
        movapd %xmm0,nb233_fjx(%esp)
        movapd %xmm1,nb233_fjy(%esp)
        movapd %xmm2,nb233_fjz(%esp)

        ## M interactions 
        movapd  nb233_rinvM(%esp),%xmm5
        movapd  %xmm5,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        movapd  %xmm5,%xmm7
        movapd  nb233_krsqM(%esp),%xmm0
        addpd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        mulpd   nb233_two(%esp),%xmm0
        subpd   nb233_crf(%esp),%xmm5
        subpd   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulpd   nb233_qqM(%esp),%xmm5   ## vcoul 
        mulpd   nb233_qqM(%esp),%xmm7
        mulpd  %xmm7,%xmm4              ## total fsH2 in xmm4 

        addpd  nb233_vctot(%esp),%xmm5

        movapd nb233_dxM(%esp),%xmm0
        movapd nb233_dyM(%esp),%xmm1
        movapd nb233_dzM(%esp),%xmm2
        movapd %xmm5,nb233_vctot(%esp)
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        ## update H2 forces 
        movapd nb233_fixM(%esp),%xmm3
        movapd nb233_fiyM(%esp),%xmm4
        movapd nb233_fizM(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb233_fixM(%esp)
        movapd %xmm4,nb233_fiyM(%esp)
        movapd %xmm7,nb233_fizM(%esp)

        movl nb233_faction(%ebp),%edi
        ## update j forces 
        addpd  nb233_fjx(%esp),%xmm0
        addpd  nb233_fjy(%esp),%xmm1
        addpd  nb233_fjz(%esp),%xmm2
        movlpd (%edi,%eax,8),%xmm3
        movlpd 8(%edi,%eax,8),%xmm4
        movlpd 16(%edi,%eax,8),%xmm5
        movhpd (%edi,%ebx,8),%xmm3
        movhpd 8(%edi,%ebx,8),%xmm4
        movhpd 16(%edi,%ebx,8),%xmm5
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        movlpd %xmm3,(%edi,%eax,8)
        movlpd %xmm4,8(%edi,%eax,8)
        movlpd %xmm5,16(%edi,%eax,8)
        movhpd %xmm3,(%edi,%ebx,8)
        movhpd %xmm4,8(%edi,%ebx,8)
        movhpd %xmm5,16(%edi,%ebx,8)

        ## should we do one more iteration? 
        subl $2,nb233_innerk(%esp)
        jl   _nb_kernel233_ia32_sse2.nb233_checksingle
        jmp  _nb_kernel233_ia32_sse2.nb233_unroll_loop
_nb_kernel233_ia32_sse2.nb233_checksingle: 
        movl  nb233_innerk(%esp),%edx
        andl  $1,%edx
        jnz  _nb_kernel233_ia32_sse2.nb233_dosingle
        jmp  _nb_kernel233_ia32_sse2.nb233_updateouterdata
_nb_kernel233_ia32_sse2.nb233_dosingle: 
        movl  nb233_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb233_innerjjnr(%esp)

        movl nb233_charge(%ebp),%esi     ## base of charge[] 

        xorpd %xmm3,%xmm3
        movlpd (%esi,%eax,8),%xmm3
        movapd %xmm3,%xmm4
        mulsd  nb233_iqM(%esp),%xmm3
        mulsd  nb233_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 

        movapd  %xmm3,nb233_qqM(%esp)
        movapd  %xmm4,nb233_qqH(%esp)

        movl nb233_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb233_vdwparam(%ebp),%esi
        shll %eax
        movl nb233_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movapd %xmm4,nb233_c6(%esp)
        movapd %xmm6,nb233_c12(%esp)

        movl nb233_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb233_ixO(%esp),%xmm4
        movapd nb233_iyO(%esp),%xmm5
        movapd nb233_izO(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb233_dxO(%esp)
        movapd %xmm5,nb233_dyO(%esp)
        movapd %xmm6,nb233_dzO(%esp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 
        movapd %xmm7,nb233_rsqO(%esp)

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb233_ixH1(%esp),%xmm4
        movapd nb233_iyH1(%esp),%xmm5
        movapd nb233_izH1(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb233_dxH1(%esp)
        movapd %xmm5,nb233_dyH1(%esp)
        movapd %xmm6,nb233_dzH1(%esp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm6
        addsd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb233_ixH2(%esp),%xmm3
        movapd nb233_iyH2(%esp),%xmm4
        movapd nb233_izH2(%esp),%xmm5

        ## calc dr 
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5

        ## store dr 
        movapd %xmm3,nb233_dxH2(%esp)
        movapd %xmm4,nb233_dyH2(%esp)
        movapd %xmm5,nb233_dzH2(%esp)
        ## square it 
        mulsd %xmm3,%xmm3
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        addsd %xmm4,%xmm5
        addsd %xmm3,%xmm5
        ## move ixM-izM to xmm2-xmm4  
        movapd nb233_iyM(%esp),%xmm3
        movapd nb233_izM(%esp),%xmm4
        subpd  %xmm1,%xmm3
        subpd  %xmm2,%xmm4
        movapd nb233_ixM(%esp),%xmm2
        subpd  %xmm0,%xmm2

        ## store dr 
        movapd %xmm2,nb233_dxM(%esp)
        movapd %xmm3,nb233_dyM(%esp)
        movapd %xmm4,nb233_dzM(%esp)
        ## square it 
        mulpd %xmm2,%xmm2
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        addpd %xmm3,%xmm4
        addpd %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## calculate krsq
        movsd nb233_krf(%esp),%xmm0
        movsd %xmm0,%xmm1
        movsd %xmm0,%xmm2
        mulsd %xmm4,%xmm0
        mulsd %xmm5,%xmm1
        mulsd %xmm6,%xmm2
        movsd %xmm0,nb233_krsqM(%esp)
        movsd %xmm1,nb233_krsqH2(%esp)
        movsd %xmm2,nb233_krsqH1(%esp)

        ## start with rsqH1 - put seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb233_three(%esp),%xmm1
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb233_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm6       ## rsq*lu*lu 
        movapd nb233_three(%esp),%xmm1
        subsd %xmm6,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb233_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb233_rinvH1(%esp)

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb233_three(%esp),%xmm1
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb233_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm5       ## rsq*lu*lu 
        movapd nb233_three(%esp),%xmm1
        subsd %xmm5,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb233_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb233_rinvH2(%esp)

        ## rsqM - seed in xmm2 
        cvtsd2ss %xmm4,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb233_three(%esp),%xmm1
        mulsd   %xmm4,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb233_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm4       ## rsq*lu*lu 
        movapd nb233_three(%esp),%xmm1
        subsd %xmm4,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb233_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb233_rinvM(%esp)

        ## rsqO - put seed in xmm2 
        cvtsd2ss %xmm7,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movsd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movsd  nb233_three(%esp),%xmm4
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb233_half(%esp),%xmm4   ## iter1 ( new lu) 

        movsd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm7       ## rsq*lu*lu 
        movsd nb233_three(%esp),%xmm4
        subsd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb233_half(%esp),%xmm4   ## rinv 
        movsd  %xmm4,%xmm7      ## rinvO in xmm7 

        movsd nb233_rsqO(%esp),%xmm4
        movapd %xmm7,%xmm0
        ## LJ table interaction.
        mulsd %xmm7,%xmm4       ## xmm4=r 
        mulsd nb233_tsc(%esp),%xmm4

        cvttsd2si %xmm4,%ebx    ## mm6 = lu idx 
        cvtsi2sd %ebx,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $3,%ebx

        movl nb233_VFtab(%ebp),%esi

        ## dispersion 
        movlpd (%esi,%ebx,8),%xmm4      ## Y1   
        movhpd 8(%esi,%ebx,8),%xmm4     ## Y1 F1        
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%ebx,8),%xmm6    ## G1
        movhpd 24(%esi,%ebx,8),%xmm6    ## G1 H1        
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## dispersion table ready, in xmm4-xmm7         
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb233_two(%esp),%xmm7    ## two*Heps2 
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movsd nb233_c6(%esp),%xmm4
        mulsd  %xmm4,%xmm7       ## fijD 
        mulsd  %xmm4,%xmm5       ## Vvdw6 

        ## put scalar force on stack Update Vvdwtot directly 
        addsd  nb233_Vvdwtot(%esp),%xmm5
        xorpd  %xmm3,%xmm3
        mulsd  nb233_tsc(%esp),%xmm7
        subsd  %xmm7,%xmm3
        movsd %xmm3,nb233_fstmp(%esp)
        movsd %xmm5,nb233_Vvdwtot(%esp)

        ## repulsion 
        movlpd 32(%esi,%ebx,8),%xmm4    ## Y1   
        movhpd 40(%esi,%ebx,8),%xmm4    ## Y1 F1        

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 48(%esi,%ebx,8),%xmm6    ## G1
        movhpd 56(%esi,%ebx,8),%xmm6    ## G1 H1        

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 

        ## table ready, in xmm4-xmm7    
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb233_two(%esp),%xmm7    ## two*Heps2 
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movsd nb233_c12(%esp),%xmm4
        mulsd  %xmm4,%xmm7
        mulsd  %xmm4,%xmm5

        addsd  nb233_Vvdwtot(%esp),%xmm5
        movsd nb233_fstmp(%esp),%xmm3
        mulsd  nb233_tsc(%esp),%xmm7
        subsd  %xmm7,%xmm3
        movsd %xmm5,nb233_Vvdwtot(%esp)

        mulsd  %xmm0,%xmm3


        movsd nb233_dxO(%esp),%xmm0
        movsd nb233_dyO(%esp),%xmm1
        movsd nb233_dzO(%esp),%xmm2

        movl   nb233_faction(%ebp),%edi
        mulsd  %xmm3,%xmm0
        mulsd  %xmm3,%xmm1
        mulsd  %xmm3,%xmm2

        ## update O forces 
        movapd nb233_fixO(%esp),%xmm3
        movapd nb233_fiyO(%esp),%xmm4
        movapd nb233_fizO(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movsd %xmm3,nb233_fixO(%esp)
        movsd %xmm4,nb233_fiyO(%esp)
        movsd %xmm7,nb233_fizO(%esp)
        ## update j forces with water O 
        movsd %xmm0,nb233_fjx(%esp)
        movsd %xmm1,nb233_fjy(%esp)
        movsd %xmm2,nb233_fjz(%esp)

        ## H1 interactions
        movsd  nb233_rinvH1(%esp),%xmm6
        movsd  %xmm6,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        movsd  %xmm6,%xmm7
        movsd  nb233_krsqH1(%esp),%xmm0
        addsd   %xmm0,%xmm6     ## xmm6=rinv+ krsq 
        mulsd   nb233_two(%esp),%xmm0
        subsd   nb233_crf(%esp),%xmm6
        subsd   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulsd   nb233_qqH(%esp),%xmm6   ## vcoul 
        mulsd   nb233_qqH(%esp),%xmm7
        mulsd  %xmm7,%xmm4              ## total fsH1 in xmm4 

        addsd  nb233_vctot(%esp),%xmm6

        movapd nb233_dxH1(%esp),%xmm0
        movapd nb233_dyH1(%esp),%xmm1
        movapd nb233_dzH1(%esp),%xmm2
        movsd %xmm6,nb233_vctot(%esp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb233_fixH1(%esp),%xmm3
        movapd nb233_fiyH1(%esp),%xmm4
        movapd nb233_fizH1(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movsd %xmm3,nb233_fixH1(%esp)
        movsd %xmm4,nb233_fiyH1(%esp)
        movsd %xmm7,nb233_fizH1(%esp)
        ## update j forces with water H1 
        addsd  nb233_fjx(%esp),%xmm0
        addsd  nb233_fjy(%esp),%xmm1
        addsd  nb233_fjz(%esp),%xmm2
        movsd %xmm0,nb233_fjx(%esp)
        movsd %xmm1,nb233_fjy(%esp)
        movsd %xmm2,nb233_fjz(%esp)

        ## H2 interactions 
        movsd  nb233_rinvH2(%esp),%xmm5
        movsd  %xmm5,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        movsd  %xmm5,%xmm7
        movsd  nb233_krsqH2(%esp),%xmm0
        addsd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        mulsd   nb233_two(%esp),%xmm0
        subsd   nb233_crf(%esp),%xmm5
        subsd   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulsd   nb233_qqH(%esp),%xmm5   ## vcoul 
        mulsd   nb233_qqH(%esp),%xmm7
        mulsd  %xmm7,%xmm4              ## total fsH2 in xmm4 

        addsd  nb233_vctot(%esp),%xmm5

        movapd nb233_dxH2(%esp),%xmm0
        movapd nb233_dyH2(%esp),%xmm1
        movapd nb233_dzH2(%esp),%xmm2
        movsd %xmm5,nb233_vctot(%esp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H2 forces 
        movapd nb233_fixH2(%esp),%xmm3
        movapd nb233_fiyH2(%esp),%xmm4
        movapd nb233_fizH2(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movsd %xmm3,nb233_fixH2(%esp)
        movsd %xmm4,nb233_fiyH2(%esp)
        movsd %xmm7,nb233_fizH2(%esp)
        ## update j forces with water H2 
        addsd  nb233_fjx(%esp),%xmm0
        addsd  nb233_fjy(%esp),%xmm1
        addsd  nb233_fjz(%esp),%xmm2
        movsd %xmm0,nb233_fjx(%esp)
        movsd %xmm1,nb233_fjy(%esp)
        movsd %xmm2,nb233_fjz(%esp)

        ## M interactions 
        movsd  nb233_rinvM(%esp),%xmm5
        movsd  %xmm5,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        movsd  %xmm5,%xmm7
        movsd  nb233_krsqM(%esp),%xmm0
        addsd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        mulsd   nb233_two(%esp),%xmm0
        subsd   nb233_crf(%esp),%xmm5
        subsd   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulsd   nb233_qqM(%esp),%xmm5   ## vcoul 
        mulsd   nb233_qqM(%esp),%xmm7
        mulsd  %xmm7,%xmm4              ## total fsH2 in xmm4 

        addsd  nb233_vctot(%esp),%xmm5

        movapd nb233_dxM(%esp),%xmm0
        movapd nb233_dyM(%esp),%xmm1
        movapd nb233_dzM(%esp),%xmm2
        movsd %xmm5,nb233_vctot(%esp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update M forces 
        movapd nb233_fixM(%esp),%xmm3
        movapd nb233_fiyM(%esp),%xmm4
        movapd nb233_fizM(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movsd %xmm3,nb233_fixM(%esp)
        movsd %xmm4,nb233_fiyM(%esp)
        movsd %xmm7,nb233_fizM(%esp)

        movl nb233_faction(%ebp),%edi
        ## update j forces 
        addsd  nb233_fjx(%esp),%xmm0
        addsd  nb233_fjy(%esp),%xmm1
        addsd  nb233_fjz(%esp),%xmm2
        movlpd (%edi,%eax,8),%xmm3
        movlpd 8(%edi,%eax,8),%xmm4
        movlpd 16(%edi,%eax,8),%xmm5
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        movlpd %xmm3,(%edi,%eax,8)
        movlpd %xmm4,8(%edi,%eax,8)
        movlpd %xmm5,16(%edi,%eax,8)

_nb_kernel233_ia32_sse2.nb233_updateouterdata: 
        movl  nb233_ii3(%esp),%ecx
        movl  nb233_faction(%ebp),%edi
        movl  nb233_fshift(%ebp),%esi
        movl  nb233_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movapd nb233_fixO(%esp),%xmm0
        movapd nb233_fiyO(%esp),%xmm1
        movapd nb233_fizO(%esp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        ## increment i force 
        movsd  (%edi,%ecx,8),%xmm3
        movsd  8(%edi,%ecx,8),%xmm4
        movsd  16(%edi,%ecx,8),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movsd  %xmm3,(%edi,%ecx,8)
        movsd  %xmm4,8(%edi,%ecx,8)
        movsd  %xmm5,16(%edi,%ecx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        movapd %xmm0,%xmm6
        movsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm6

        ## accumulate H1i forces in xmm0, xmm1, xmm2 
        movapd nb233_fixH1(%esp),%xmm0
        movapd nb233_fiyH1(%esp),%xmm1
        movapd nb233_fizH1(%esp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        ## increment i force 
        movsd  24(%edi,%ecx,8),%xmm3
        movsd  32(%edi,%ecx,8),%xmm4
        movsd  40(%edi,%ecx,8),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movsd  %xmm3,24(%edi,%ecx,8)
        movsd  %xmm4,32(%edi,%ecx,8)
        movsd  %xmm5,40(%edi,%ecx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        addsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm0
        addpd %xmm0,%xmm6

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movapd nb233_fixH2(%esp),%xmm0
        movapd nb233_fiyH2(%esp),%xmm1
        movapd nb233_fizH2(%esp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        ## increment i force 
        movsd  48(%edi,%ecx,8),%xmm3
        movsd  56(%edi,%ecx,8),%xmm4
        movsd  64(%edi,%ecx,8),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movsd  %xmm3,48(%edi,%ecx,8)
        movsd  %xmm4,56(%edi,%ecx,8)
        movsd  %xmm5,64(%edi,%ecx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        addsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm0
        addpd %xmm0,%xmm6

        ## accumulate Mi forces in xmm0, xmm1, xmm2 
        movapd nb233_fixM(%esp),%xmm0
        movapd nb233_fiyM(%esp),%xmm1
        movapd nb233_fizM(%esp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        ## increment i force 
        movsd  72(%edi,%ecx,8),%xmm3
        movsd  80(%edi,%ecx,8),%xmm4
        movsd  88(%edi,%ecx,8),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movsd  %xmm3,72(%edi,%ecx,8)
        movsd  %xmm4,80(%edi,%ecx,8)
        movsd  %xmm5,88(%edi,%ecx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        addsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm0
        addpd %xmm0,%xmm6

        ## increment fshift force 
        movlpd (%esi,%edx,8),%xmm3
        movhpd 8(%esi,%edx,8),%xmm3
        movsd  16(%esi,%edx,8),%xmm4
        addpd  %xmm6,%xmm3
        addsd  %xmm7,%xmm4
        movlpd %xmm3,(%esi,%edx,8)
        movhpd %xmm3,8(%esi,%edx,8)
        movsd  %xmm4,16(%esi,%edx,8)

        ## get n from stack
        movl nb233_n(%esp),%esi
        ## get group index for i particle 
        movl  nb233_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb233_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb233_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb233_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb233_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

       ## finish if last 
        movl nb233_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel233_ia32_sse2.nb233_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb233_n(%esp)
        jmp _nb_kernel233_ia32_sse2.nb233_outer
_nb_kernel233_ia32_sse2.nb233_outerend: 
        ## check if more outer neighborlists remain
        movl  nb233_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel233_ia32_sse2.nb233_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel233_ia32_sse2.nb233_threadloop
_nb_kernel233_ia32_sse2.nb233_end: 
        emms

        movl nb233_nouter(%esp),%eax
        movl nb233_ninner(%esp),%ebx
        movl nb233_outeriter(%ebp),%ecx
        movl nb233_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb233_salign(%esp),%eax
        addl %eax,%esp
        addl $1020,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret




.globl nb_kernel233nf_ia32_sse2
.globl _nb_kernel233nf_ia32_sse2
nb_kernel233nf_ia32_sse2:       
_nb_kernel233nf_ia32_sse2:      
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
        ## bottom of stack is cache-aligned for sse2 use 
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
.set nb233nf_iqH, 192
.set nb233nf_iqM, 208
.set nb233nf_qqH, 416
.set nb233nf_qqM, 432
.set nb233nf_c6, 448
.set nb233nf_c12, 464
.set nb233nf_tsc, 480
.set nb233nf_vctot, 512
.set nb233nf_Vvdwtot, 528
.set nb233nf_half, 784
.set nb233nf_three, 800
.set nb233nf_two, 816
.set nb233nf_rinvH1, 832
.set nb233nf_rinvH2, 848
.set nb233nf_rinvM, 864
.set nb233nf_krsqH1, 880
.set nb233nf_krsqH2, 896
.set nb233nf_krsqM, 912
.set nb233nf_krf, 928
.set nb233nf_crf, 944
.set nb233nf_rsqO, 960
.set nb233nf_is3, 976
.set nb233nf_ii3, 980
.set nb233nf_ntia, 984
.set nb233nf_innerjjnr, 988
.set nb233nf_innerk, 992
.set nb233nf_n, 996
.set nb233nf_nn1, 1000
.set nb233nf_nri, 1004
.set nb233nf_nouter, 1008
.set nb233nf_ninner, 1012
.set nb233nf_salign, 1016
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $1020,%esp         ## local stack space 
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
        movsd (%eax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb233nf_tsc(%esp)

        movl nb233nf_argkrf(%ebp),%esi
        movl nb233nf_argcrf(%ebp),%edi
        movsd (%esi),%xmm5
        movsd (%edi),%xmm6
        shufpd $0,%xmm5,%xmm5
        shufpd $0,%xmm6,%xmm6
        movapd %xmm5,nb233nf_krf(%esp)
        movapd %xmm6,nb233nf_crf(%esp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb233nf_half(%esp)
        movl %ebx,nb233nf_half+4(%esp)
        movsd nb233nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb233nf_half(%esp)
        movapd %xmm2,nb233nf_two(%esp)
        movapd %xmm3,nb233nf_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb233nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb233nf_charge(%ebp),%edx
        movsd 8(%edx,%ebx,8),%xmm3
        movsd 24(%edx,%ebx,8),%xmm4
        movl nb233nf_p_facel(%ebp),%esi
        movsd (%esi),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb233nf_iqH(%esp)
        movapd %xmm4,nb233nf_iqM(%esp)

        movl  nb233nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb233nf_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb233nf_ntia(%esp)
_nb_kernel233nf_ia32_sse2.nb233nf_threadloop: 
        movl  nb233nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel233nf_ia32_sse2.nb233nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel233nf_ia32_sse2.nb233nf_spinlock

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
        jg  _nb_kernel233nf_ia32_sse2.nb233nf_outerstart
        jmp _nb_kernel233nf_ia32_sse2.nb233nf_end

_nb_kernel233nf_ia32_sse2.nb233nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb233nf_nouter(%esp),%ebx
        movl %ebx,nb233nf_nouter(%esp)

_nb_kernel233nf_ia32_sse2.nb233nf_outer: 
        movl  nb233nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb233nf_is3(%esp)            ## store is3 

        movl  nb233nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb233nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        movapd %xmm0,%xmm6
        movapd %xmm1,%xmm7

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb233nf_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb233nf_ii3(%esp)

        addsd (%eax,%ebx,8),%xmm3       ## ox
        addsd 8(%eax,%ebx,8),%xmm4      ## oy
        addsd 16(%eax,%ebx,8),%xmm5     ## oz   
        addsd 24(%eax,%ebx,8),%xmm6     ## h1x
        addsd 32(%eax,%ebx,8),%xmm7     ## h1y
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        shufpd $0,%xmm6,%xmm6
        shufpd $0,%xmm7,%xmm7
        movapd %xmm3,nb233nf_ixO(%esp)
        movapd %xmm4,nb233nf_iyO(%esp)
        movapd %xmm5,nb233nf_izO(%esp)
        movapd %xmm6,nb233nf_ixH1(%esp)
        movapd %xmm7,nb233nf_iyH1(%esp)

        movsd %xmm2,%xmm6
        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 40(%eax,%ebx,8),%xmm6    ## h1z
        addsd 48(%eax,%ebx,8),%xmm0    ## h2x
        addsd 56(%eax,%ebx,8),%xmm1    ## h2y
        addsd 64(%eax,%ebx,8),%xmm2    ## h2z
        addsd 72(%eax,%ebx,8),%xmm3    ## mx
        addsd 80(%eax,%ebx,8),%xmm4    ## my
        addsd 88(%eax,%ebx,8),%xmm5    ## mz

        shufpd $0,%xmm6,%xmm6
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm6,nb233nf_izH1(%esp)
        movapd %xmm0,nb233nf_ixH2(%esp)
        movapd %xmm1,nb233nf_iyH2(%esp)
        movapd %xmm2,nb233nf_izH2(%esp)
        movapd %xmm3,nb233nf_ixM(%esp)
        movapd %xmm4,nb233nf_iyM(%esp)
        movapd %xmm5,nb233nf_izM(%esp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb233nf_vctot(%esp)
        movapd %xmm4,nb233nf_Vvdwtot(%esp)

        movl  nb233nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb233nf_pos(%ebp),%esi
        movl  nb233nf_faction(%ebp),%edi
        movl  nb233nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb233nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb233nf_ninner(%esp),%ecx
        movl  %ecx,nb233nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb233nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel233nf_ia32_sse2.nb233nf_unroll_loop
        jmp   _nb_kernel233nf_ia32_sse2.nb233nf_checksingle
_nb_kernel233nf_ia32_sse2.nb233nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb233nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb233nf_innerjjnr(%esp)                 ## advance pointer (unrolled 2) 

        movl nb233nf_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb233nf_iqM(%esp),%xmm3
        mulpd  nb233nf_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movapd  %xmm3,nb233nf_qqM(%esp)
        movapd  %xmm4,nb233nf_qqH(%esp)

        movl nb233nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb233nf_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb233nf_ntia(%esp),%edi
        addl %edi,%eax
        addl %edi,%ebx

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movlpd (%esi,%ebx,8),%xmm7      ## c6b
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 
        movhpd 8(%esi,%ebx,8),%xmm7     ## c6b c12b 
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movapd %xmm4,nb233nf_c6(%esp)
        movapd %xmm6,nb233nf_c12(%esp)

        movl nb233nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move two coordinates to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2
        movhpd (%esi,%ebx,8),%xmm0
        movhpd 8(%esi,%ebx,8),%xmm1
        movhpd 16(%esi,%ebx,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb233nf_ixO(%esp),%xmm4
        movapd nb233nf_iyO(%esp),%xmm5
        movapd nb233nf_izO(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb233nf_ixH1(%esp),%xmm4
        movapd nb233nf_iyH1(%esp),%xmm5
        movapd nb233nf_izH1(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm6
        addpd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb233nf_ixH2(%esp),%xmm3
        movapd nb233nf_iyH2(%esp),%xmm4
        movapd nb233nf_izH2(%esp),%xmm5

        ## calc dr 
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5

        ## square it 
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        addpd %xmm4,%xmm5
        addpd %xmm3,%xmm5

        ## move ixM-izM to xmm2-xmm4  
        movapd nb233nf_iyM(%esp),%xmm3
        movapd nb233nf_izM(%esp),%xmm4
        subpd  %xmm1,%xmm3
        subpd  %xmm2,%xmm4
        movapd nb233nf_ixM(%esp),%xmm2
        subpd  %xmm0,%xmm2


        ## square it 
        mulpd %xmm2,%xmm2
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        addpd %xmm3,%xmm4
        addpd %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 
        movapd %xmm7,nb233nf_rsqO(%esp)

        ## calculate krsq
        movapd nb233nf_krf(%esp),%xmm0
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2
        mulpd %xmm4,%xmm0
        mulpd %xmm5,%xmm1
        mulpd %xmm6,%xmm2
        movapd %xmm0,nb233nf_krsqM(%esp)
        movapd %xmm1,nb233nf_krsqH2(%esp)
        movapd %xmm2,nb233nf_krsqH1(%esp)

        ## start with rsqH1 - put seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb233nf_three(%esp),%xmm1
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulpd   nb233nf_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulpd %xmm1,%xmm1       ## lu*lu 
        mulpd %xmm1,%xmm6       ## rsq*lu*lu 
        movapd nb233nf_three(%esp),%xmm1
        subpd %xmm6,%xmm1       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulpd nb233nf_half(%esp),%xmm1   ## rinv 
        movapd  %xmm1,nb233nf_rinvH1(%esp)

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb233nf_three(%esp),%xmm1
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulpd   nb233nf_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulpd %xmm1,%xmm1       ## lu*lu 
        mulpd %xmm1,%xmm5       ## rsq*lu*lu 
        movapd nb233nf_three(%esp),%xmm1
        subpd %xmm5,%xmm1       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulpd nb233nf_half(%esp),%xmm1   ## rinv 
        movapd  %xmm1,nb233nf_rinvH2(%esp)

        ## rsqM - seed in xmm2 
        cvtpd2ps %xmm4,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb233nf_three(%esp),%xmm1
        mulpd   %xmm4,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulpd   nb233nf_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulpd %xmm1,%xmm1       ## lu*lu 
        mulpd %xmm1,%xmm4       ## rsq*lu*lu 
        movapd nb233nf_three(%esp),%xmm1
        subpd %xmm4,%xmm1       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulpd nb233nf_half(%esp),%xmm1   ## rinv 
        movapd  %xmm1,nb233nf_rinvM(%esp)


        ## rsqO - put seed in xmm2 
        cvtpd2ps %xmm7,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb233nf_three(%esp),%xmm4
        mulpd   %xmm7,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb233nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb233nf_three(%esp),%xmm4
        subpd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb233nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvO in xmm7 



        movapd nb233nf_rsqO(%esp),%xmm4
        movapd %xmm7,%xmm0
        ## LJ table interaction.
        mulpd %xmm7,%xmm4       ## xmm4=r 
        mulpd nb233nf_tsc(%esp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $3,%mm6           ## idx *= 8 

        movd %eax,%mm0
        movd %ebx,%mm1

        movl nb233nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx

        ## dispersion 
        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2 
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## dispersion table ready, in xmm4-xmm7         
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb233nf_c6(%esp),%xmm4
        mulpd  %xmm4,%xmm5       ## Vvdw6 

        ## Update Vvdwtot directly 
        addpd  nb233nf_Vvdwtot(%esp),%xmm5
        movapd %xmm5,nb233nf_Vvdwtot(%esp)

        ## repulsion 
        movlpd 32(%esi,%eax,8),%xmm4    ## Y1   
        movlpd 32(%esi,%ebx,8),%xmm3    ## Y2 
        movhpd 40(%esi,%eax,8),%xmm4    ## Y1 F1        
        movhpd 40(%esi,%ebx,8),%xmm3    ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 48(%esi,%eax,8),%xmm6    ## G1
        movlpd 48(%esi,%ebx,8),%xmm3    ## G2
        movhpd 56(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 56(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 

        ## table ready, in xmm4-xmm7    
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb233nf_c12(%esp),%xmm4
        mulpd  %xmm4,%xmm5

        addpd  nb233nf_Vvdwtot(%esp),%xmm5
        movapd %xmm5,nb233nf_Vvdwtot(%esp)

        ## H1 interactions 
        movapd  nb233nf_rinvH1(%esp),%xmm6
        movapd  %xmm6,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        movapd  %xmm6,%xmm7
        movapd  nb233nf_krsqH1(%esp),%xmm0
        addpd   %xmm0,%xmm6     ## xmm6=rinv+ krsq 
        subpd   nb233nf_crf(%esp),%xmm6
        mulpd   nb233nf_qqH(%esp),%xmm6   ## vcoul 
        addpd  nb233nf_vctot(%esp),%xmm6
        movapd %xmm6,nb233nf_vctot(%esp)

        ## H2 interactions 
        movapd  nb233nf_rinvH2(%esp),%xmm5
        movapd  %xmm5,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        movapd  %xmm5,%xmm7
        movapd  nb233nf_krsqH2(%esp),%xmm0
        addpd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        subpd   nb233nf_crf(%esp),%xmm5
        mulpd   nb233nf_qqH(%esp),%xmm5   ## vcoul 
        addpd  nb233nf_vctot(%esp),%xmm5
        movapd %xmm5,nb233nf_vctot(%esp)

        ## M interactions 
        movapd  nb233nf_rinvM(%esp),%xmm5
        movapd  %xmm5,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        movapd  %xmm5,%xmm7
        movapd  nb233nf_krsqM(%esp),%xmm0
        addpd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        subpd   nb233nf_crf(%esp),%xmm5
        mulpd   nb233nf_qqM(%esp),%xmm5   ## vcoul 
        addpd  nb233nf_vctot(%esp),%xmm5
        movapd %xmm5,nb233nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $2,nb233nf_innerk(%esp)
        jl   _nb_kernel233nf_ia32_sse2.nb233nf_checksingle
        jmp  _nb_kernel233nf_ia32_sse2.nb233nf_unroll_loop
_nb_kernel233nf_ia32_sse2.nb233nf_checksingle: 
        movl  nb233nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz  _nb_kernel233nf_ia32_sse2.nb233nf_dosingle
        jmp  _nb_kernel233nf_ia32_sse2.nb233nf_updateouterdata
_nb_kernel233nf_ia32_sse2.nb233nf_dosingle: 
        movl  nb233nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb233nf_innerjjnr(%esp)

        movl nb233nf_charge(%ebp),%esi     ## base of charge[] 

        xorpd %xmm3,%xmm3
        movlpd (%esi,%eax,8),%xmm3
        movapd %xmm3,%xmm4
        mulsd  nb233nf_iqM(%esp),%xmm3
        mulsd  nb233nf_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 

        movapd  %xmm3,nb233nf_qqM(%esp)
        movapd  %xmm4,nb233nf_qqH(%esp)

        movl nb233nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb233nf_vdwparam(%ebp),%esi
        shll %eax
        movl nb233nf_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movapd %xmm4,nb233nf_c6(%esp)
        movapd %xmm6,nb233nf_c12(%esp)

        movl nb233nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb233nf_ixO(%esp),%xmm4
        movapd nb233nf_iyO(%esp),%xmm5
        movapd nb233nf_izO(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 
        movapd %xmm7,nb233nf_rsqO(%esp)

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb233nf_ixH1(%esp),%xmm4
        movapd nb233nf_iyH1(%esp),%xmm5
        movapd nb233nf_izH1(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm6
        addsd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb233nf_ixH2(%esp),%xmm3
        movapd nb233nf_iyH2(%esp),%xmm4
        movapd nb233nf_izH2(%esp),%xmm5

        ## calc dr 
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5

        ## square it 
        mulsd %xmm3,%xmm3
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        addsd %xmm4,%xmm5
        addsd %xmm3,%xmm5
        ## move ixM-izM to xmm2-xmm4  
        movapd nb233nf_iyM(%esp),%xmm3
        movapd nb233nf_izM(%esp),%xmm4
        subpd  %xmm1,%xmm3
        subpd  %xmm2,%xmm4
        movapd nb233nf_ixM(%esp),%xmm2
        subpd  %xmm0,%xmm2

        ## square it 
        mulpd %xmm2,%xmm2
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        addpd %xmm3,%xmm4
        addpd %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## calculate krsq
        movsd nb233nf_krf(%esp),%xmm0
        movsd %xmm0,%xmm1
        movsd %xmm0,%xmm2
        mulsd %xmm4,%xmm0
        mulsd %xmm5,%xmm1
        mulsd %xmm6,%xmm2
        movsd %xmm0,nb233nf_krsqM(%esp)
        movsd %xmm1,nb233nf_krsqH2(%esp)
        movsd %xmm2,nb233nf_krsqH1(%esp)

        ## start with rsqH1 - put seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb233nf_three(%esp),%xmm1
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb233nf_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm6       ## rsq*lu*lu 
        movapd nb233nf_three(%esp),%xmm1
        subsd %xmm6,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb233nf_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb233nf_rinvH1(%esp)

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb233nf_three(%esp),%xmm1
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb233nf_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm5       ## rsq*lu*lu 
        movapd nb233nf_three(%esp),%xmm1
        subsd %xmm5,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb233nf_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb233nf_rinvH2(%esp)

        ## rsqM - seed in xmm2 
        cvtsd2ss %xmm4,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb233nf_three(%esp),%xmm1
        mulsd   %xmm4,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb233nf_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm4       ## rsq*lu*lu 
        movapd nb233nf_three(%esp),%xmm1
        subsd %xmm4,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb233nf_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb233nf_rinvM(%esp)

        ## rsqO - put seed in xmm2 
        cvtsd2ss %xmm7,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movsd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movsd  nb233nf_three(%esp),%xmm4
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb233nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movsd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm7       ## rsq*lu*lu 
        movsd nb233nf_three(%esp),%xmm4
        subsd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb233nf_half(%esp),%xmm4   ## rinv 
        movsd  %xmm4,%xmm7      ## rinvO in xmm7 

        movsd nb233nf_rsqO(%esp),%xmm4
        movapd %xmm7,%xmm0
        ## LJ table interaction.
        mulsd %xmm7,%xmm4       ## xmm4=r 
        mulsd nb233nf_tsc(%esp),%xmm4

        cvttsd2si %xmm4,%ebx    ## mm6 = lu idx 
        cvtsi2sd %ebx,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $3,%ebx

        movl nb233nf_VFtab(%ebp),%esi

        ## dispersion 
        movlpd (%esi,%ebx,8),%xmm4      ## Y1   
        movhpd 8(%esi,%ebx,8),%xmm4     ## Y1 F1        
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%ebx,8),%xmm6    ## G1
        movhpd 24(%esi,%ebx,8),%xmm6    ## G1 H1        
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## dispersion table ready, in xmm4-xmm7         
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movsd nb233nf_c6(%esp),%xmm4
        mulsd  %xmm4,%xmm5       ## Vvdw6 

        ## Update Vvdwtot directly 
        addsd  nb233nf_Vvdwtot(%esp),%xmm5
        movsd %xmm5,nb233nf_Vvdwtot(%esp)

        ## repulsion 
        movlpd 32(%esi,%ebx,8),%xmm4    ## Y1   
        movhpd 40(%esi,%ebx,8),%xmm4    ## Y1 F1        

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 48(%esi,%ebx,8),%xmm6    ## G1
        movhpd 56(%esi,%ebx,8),%xmm6    ## G1 H1        

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 

        ## table ready, in xmm4-xmm7    
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movsd nb233nf_c12(%esp),%xmm4
        mulsd  %xmm4,%xmm5

        addsd  nb233nf_Vvdwtot(%esp),%xmm5
        movsd %xmm5,nb233nf_Vvdwtot(%esp)

        ## H1 interactions
        movsd  nb233nf_rinvH1(%esp),%xmm6
        movsd  %xmm6,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        movsd  %xmm6,%xmm7
        movsd  nb233nf_krsqH1(%esp),%xmm0
        addsd   %xmm0,%xmm6     ## xmm6=rinv+ krsq 
        subsd   nb233nf_crf(%esp),%xmm6
        mulsd   nb233nf_qqH(%esp),%xmm6   ## vcoul 
        addsd  nb233nf_vctot(%esp),%xmm6
        movsd %xmm6,nb233nf_vctot(%esp)

        ## H2 interactions 
        movsd  nb233nf_rinvH2(%esp),%xmm5
        movsd  %xmm5,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        movsd  %xmm5,%xmm7
        movsd  nb233nf_krsqH2(%esp),%xmm0
        addsd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        subsd   nb233nf_crf(%esp),%xmm5
        mulsd   nb233nf_qqH(%esp),%xmm5   ## vcoul 
        addsd  nb233nf_vctot(%esp),%xmm5
        movsd %xmm5,nb233nf_vctot(%esp)

        ## M interactions 
        movsd  nb233nf_rinvM(%esp),%xmm5
        movsd  %xmm5,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        movsd  %xmm5,%xmm7
        movsd  nb233nf_krsqM(%esp),%xmm0
        addsd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        subsd   nb233nf_crf(%esp),%xmm5
        mulsd   nb233nf_qqM(%esp),%xmm5   ## vcoul 
        addsd  nb233nf_vctot(%esp),%xmm5
        movsd %xmm5,nb233nf_vctot(%esp)

_nb_kernel233nf_ia32_sse2.nb233nf_updateouterdata: 
        ## get n from stack
        movl nb233nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb233nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb233nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb233nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb233nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb233nf_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

       ## finish if last 
        movl nb233nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel233nf_ia32_sse2.nb233nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb233nf_n(%esp)
        jmp _nb_kernel233nf_ia32_sse2.nb233nf_outer
_nb_kernel233nf_ia32_sse2.nb233nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb233nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel233nf_ia32_sse2.nb233nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel233nf_ia32_sse2.nb233nf_threadloop
_nb_kernel233nf_ia32_sse2.nb233nf_end: 
        emms

        movl nb233nf_nouter(%esp),%eax
        movl nb233nf_ninner(%esp),%ebx
        movl nb233nf_outeriter(%ebp),%ecx
        movl nb233nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb233nf_salign(%esp),%eax
        addl %eax,%esp
        addl $1020,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret


