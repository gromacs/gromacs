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

.globl nb_kernel133_ia32_sse2
.globl _nb_kernel133_ia32_sse2
nb_kernel133_ia32_sse2: 
_nb_kernel133_ia32_sse2:        
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
        ## bottom of stack is cache-aligned for sse2 use 
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
.set nb133_iqH, 192
.set nb133_iqM, 208
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
.set nb133_qqH, 416
.set nb133_qqM, 432
.set nb133_c6, 448
.set nb133_c12, 464
.set nb133_tsc, 480
.set nb133_fstmp, 496
.set nb133_vctot, 512
.set nb133_Vvdwtot, 528
.set nb133_fixO, 544
.set nb133_fiyO, 560
.set nb133_fizO, 576
.set nb133_fixH1, 592
.set nb133_fiyH1, 608
.set nb133_fizH1, 624
.set nb133_fixH2, 640
.set nb133_fiyH2, 656
.set nb133_fizH2, 672
.set nb133_fixM, 688
.set nb133_fiyM, 704
.set nb133_fizM, 720
.set nb133_fjx, 736
.set nb133_fjy, 752
.set nb133_fjz, 768
.set nb133_half, 784
.set nb133_three, 800
.set nb133_two, 816
.set nb133_rinvH1, 832
.set nb133_rinvH2, 848
.set nb133_rinvM, 864
.set nb133_krsqH1, 880
.set nb133_krsqH2, 896
.set nb133_krsqM, 912
.set nb133_rsqO, 960
.set nb133_is3, 976
.set nb133_ii3, 980
.set nb133_ntia, 984
.set nb133_innerjjnr, 988
.set nb133_innerk, 992
.set nb133_n, 996
.set nb133_nn1, 1000
.set nb133_nri, 1004
.set nb133_nouter, 1008
.set nb133_ninner, 1012
.set nb133_salign, 1016
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
        movsd (%eax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb133_tsc(%esp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb133_half(%esp)
        movl %ebx,nb133_half+4(%esp)
        movsd nb133_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb133_half(%esp)
        movapd %xmm2,nb133_two(%esp)
        movapd %xmm3,nb133_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb133_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb133_charge(%ebp),%edx
        movsd 8(%edx,%ebx,8),%xmm3
        movsd 24(%edx,%ebx,8),%xmm4
        movl nb133_p_facel(%ebp),%esi
        movsd (%esi),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb133_iqH(%esp)
        movapd %xmm4,nb133_iqM(%esp)

        movl  nb133_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb133_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb133_ntia(%esp)
_nb_kernel133_ia32_sse2.nb133_threadloop: 
        movl  nb133_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel133_ia32_sse2.nb133_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel133_ia32_sse2.nb133_spinlock

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
        jg  _nb_kernel133_ia32_sse2.nb133_outerstart
        jmp _nb_kernel133_ia32_sse2.nb133_end

_nb_kernel133_ia32_sse2.nb133_outerstart: 
        ## ebx contains number of outer iterations
        addl nb133_nouter(%esp),%ebx
        movl %ebx,nb133_nouter(%esp)

_nb_kernel133_ia32_sse2.nb133_outer: 
        movl  nb133_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb133_is3(%esp)      ## store is3 

        movl  nb133_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb133_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        movapd %xmm0,%xmm6
        movapd %xmm1,%xmm7

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb133_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb133_ii3(%esp)

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
        movapd %xmm3,nb133_ixO(%esp)
        movapd %xmm4,nb133_iyO(%esp)
        movapd %xmm5,nb133_izO(%esp)
        movapd %xmm6,nb133_ixH1(%esp)
        movapd %xmm7,nb133_iyH1(%esp)

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
        movapd %xmm6,nb133_izH1(%esp)
        movapd %xmm0,nb133_ixH2(%esp)
        movapd %xmm1,nb133_iyH2(%esp)
        movapd %xmm2,nb133_izH2(%esp)
        movapd %xmm3,nb133_ixM(%esp)
        movapd %xmm4,nb133_iyM(%esp)
        movapd %xmm5,nb133_izM(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb133_vctot(%esp)
        movapd %xmm4,nb133_Vvdwtot(%esp)
        movapd %xmm4,nb133_fixO(%esp)
        movapd %xmm4,nb133_fiyO(%esp)
        movapd %xmm4,nb133_fizO(%esp)
        movapd %xmm4,nb133_fixH1(%esp)
        movapd %xmm4,nb133_fiyH1(%esp)
        movapd %xmm4,nb133_fizH1(%esp)
        movapd %xmm4,nb133_fixH2(%esp)
        movapd %xmm4,nb133_fiyH2(%esp)
        movapd %xmm4,nb133_fizH2(%esp)
        movapd %xmm4,nb133_fixM(%esp)
        movapd %xmm4,nb133_fiyM(%esp)
        movapd %xmm4,nb133_fizM(%esp)

        movl  nb133_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb133_pos(%ebp),%esi
        movl  nb133_faction(%ebp),%edi
        movl  nb133_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb133_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb133_ninner(%esp),%ecx
        movl  %ecx,nb133_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb133_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel133_ia32_sse2.nb133_unroll_loop
        jmp   _nb_kernel133_ia32_sse2.nb133_checksingle
_nb_kernel133_ia32_sse2.nb133_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb133_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb133_innerjjnr(%esp)                   ## advance pointer (unrolled 2) 

        movl nb133_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb133_iqM(%esp),%xmm3
        mulpd  nb133_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movapd  %xmm3,nb133_qqM(%esp)
        movapd  %xmm4,nb133_qqH(%esp)

        movl nb133_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb133_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb133_ntia(%esp),%edi
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
        movapd %xmm4,nb133_c6(%esp)
        movapd %xmm6,nb133_c12(%esp)

        movl nb133_pos(%ebp),%esi        ## base of pos[] 

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
        movapd nb133_ixO(%esp),%xmm4
        movapd nb133_iyO(%esp),%xmm5
        movapd nb133_izO(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb133_dxO(%esp)
        movapd %xmm5,nb133_dyO(%esp)
        movapd %xmm6,nb133_dzO(%esp)
        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb133_ixH1(%esp),%xmm4
        movapd nb133_iyH1(%esp),%xmm5
        movapd nb133_izH1(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb133_dxH1(%esp)
        movapd %xmm5,nb133_dyH1(%esp)
        movapd %xmm6,nb133_dzH1(%esp)
        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm6
        addpd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb133_ixH2(%esp),%xmm3
        movapd nb133_iyH2(%esp),%xmm4
        movapd nb133_izH2(%esp),%xmm5

        ## calc dr 
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5

        ## store dr 
        movapd %xmm3,nb133_dxH2(%esp)
        movapd %xmm4,nb133_dyH2(%esp)
        movapd %xmm5,nb133_dzH2(%esp)
        ## square it 
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        addpd %xmm4,%xmm5
        addpd %xmm3,%xmm5

        ## move ixM-izM to xmm2-xmm4  
        movapd nb133_iyM(%esp),%xmm3
        movapd nb133_izM(%esp),%xmm4
        subpd  %xmm1,%xmm3
        subpd  %xmm2,%xmm4
        movapd nb133_ixM(%esp),%xmm2
        subpd  %xmm0,%xmm2

        ## store dr 
        movapd %xmm2,nb133_dxM(%esp)
        movapd %xmm3,nb133_dyM(%esp)
        movapd %xmm4,nb133_dzM(%esp)
        ## square it 
        mulpd %xmm2,%xmm2
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        addpd %xmm3,%xmm4
        addpd %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 
        movapd %xmm7,nb133_rsqO(%esp)

        ## start with rsqH1 - put seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb133_three(%esp),%xmm1
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulpd   nb133_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulpd %xmm1,%xmm1       ## lu*lu 
        mulpd %xmm1,%xmm6       ## rsq*lu*lu 
        movapd nb133_three(%esp),%xmm1
        subpd %xmm6,%xmm1       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulpd nb133_half(%esp),%xmm1   ## rinv 
        movapd  %xmm1,nb133_rinvH1(%esp)

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb133_three(%esp),%xmm1
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulpd   nb133_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulpd %xmm1,%xmm1       ## lu*lu 
        mulpd %xmm1,%xmm5       ## rsq*lu*lu 
        movapd nb133_three(%esp),%xmm1
        subpd %xmm5,%xmm1       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulpd nb133_half(%esp),%xmm1   ## rinv 
        movapd  %xmm1,nb133_rinvH2(%esp)

        ## rsqM - seed in xmm2 
        cvtpd2ps %xmm4,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb133_three(%esp),%xmm1
        mulpd   %xmm4,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulpd   nb133_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulpd %xmm1,%xmm1       ## lu*lu 
        mulpd %xmm1,%xmm4       ## rsq*lu*lu 
        movapd nb133_three(%esp),%xmm1
        subpd %xmm4,%xmm1       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulpd nb133_half(%esp),%xmm1   ## rinv 
        movapd  %xmm1,nb133_rinvM(%esp)


        ## rsqO - put seed in xmm2 
        cvtpd2ps %xmm7,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb133_three(%esp),%xmm4
        mulpd   %xmm7,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb133_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb133_three(%esp),%xmm4
        subpd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb133_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvO in xmm7 

        movapd nb133_rsqO(%esp),%xmm4
        movapd %xmm7,%xmm0
        ## LJ table interaction.
        mulpd %xmm7,%xmm4       ## xmm4=r 
        mulpd nb133_tsc(%esp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $3,%mm6           ## idx *= 8 

        movd %eax,%mm0
        movd %ebx,%mm1

        movl nb133_VFtab(%ebp),%esi
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
        mulpd  nb133_two(%esp),%xmm7    ## two*Heps2 
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb133_c6(%esp),%xmm4
        mulpd  %xmm4,%xmm7       ## fijD 
        mulpd  %xmm4,%xmm5       ## Vvdw6 

        ## put scalar force on stack Update Vvdwtot directly 
        addpd  nb133_Vvdwtot(%esp),%xmm5
        xorpd  %xmm3,%xmm3
        mulpd  nb133_tsc(%esp),%xmm7
        subpd  %xmm7,%xmm3
        movapd %xmm3,nb133_fstmp(%esp)
        movapd %xmm5,nb133_Vvdwtot(%esp)

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
        mulpd  nb133_two(%esp),%xmm7    ## two*Heps2 
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb133_c12(%esp),%xmm4
        mulpd  %xmm4,%xmm7
        mulpd  %xmm4,%xmm5

        addpd  nb133_Vvdwtot(%esp),%xmm5
        movapd nb133_fstmp(%esp),%xmm3
        mulpd  nb133_tsc(%esp),%xmm7
        subpd  %xmm7,%xmm3
        movapd %xmm5,nb133_Vvdwtot(%esp)

        mulpd  %xmm0,%xmm3


        movapd nb133_dxO(%esp),%xmm0
        movapd nb133_dyO(%esp),%xmm1
        movapd nb133_dzO(%esp),%xmm2

        movd %mm0,%eax
        movd %mm1,%ebx

        movl   nb133_faction(%ebp),%edi
        mulpd  %xmm3,%xmm0
        mulpd  %xmm3,%xmm1
        mulpd  %xmm3,%xmm2

        ## update O forces 
        movapd nb133_fixO(%esp),%xmm3
        movapd nb133_fiyO(%esp),%xmm4
        movapd nb133_fizO(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb133_fixO(%esp)
        movapd %xmm4,nb133_fiyO(%esp)
        movapd %xmm7,nb133_fizO(%esp)
        ## update j forces with water O 
        movapd %xmm0,nb133_fjx(%esp)
        movapd %xmm1,nb133_fjy(%esp)
        movapd %xmm2,nb133_fjz(%esp)

        ## H1 interactions 
        movapd  nb133_rinvH1(%esp),%xmm6
        movapd  %xmm6,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        mulpd   nb133_qqH(%esp),%xmm6   ## vcoul 
        mulpd   %xmm6,%xmm4   ## fscal
        addpd   nb133_vctot(%esp),%xmm6
        movapd  %xmm6,nb133_vctot(%esp)

        movapd nb133_dxH1(%esp),%xmm0
        movapd nb133_dyH1(%esp),%xmm1
        movapd nb133_dzH1(%esp),%xmm2
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb133_fixH1(%esp),%xmm3
        movapd nb133_fiyH1(%esp),%xmm4
        movapd nb133_fizH1(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb133_fixH1(%esp)
        movapd %xmm4,nb133_fiyH1(%esp)
        movapd %xmm7,nb133_fizH1(%esp)
        ## update j forces with water H1 
        addpd  nb133_fjx(%esp),%xmm0
        addpd  nb133_fjy(%esp),%xmm1
        addpd  nb133_fjz(%esp),%xmm2
        movapd %xmm0,nb133_fjx(%esp)
        movapd %xmm1,nb133_fjy(%esp)
        movapd %xmm2,nb133_fjz(%esp)

        ## H2 interactions 
        movapd  nb133_rinvH2(%esp),%xmm6
        movapd  %xmm6,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        mulpd   nb133_qqH(%esp),%xmm6   ## vcoul 
        mulpd   %xmm6,%xmm4   ## fscal
        addpd  nb133_vctot(%esp),%xmm6
        movapd %xmm6,nb133_vctot(%esp)

        movapd nb133_dxH2(%esp),%xmm0
        movapd nb133_dyH2(%esp),%xmm1
        movapd nb133_dzH2(%esp),%xmm2
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        ## update H2 forces 
        movapd nb133_fixH2(%esp),%xmm3
        movapd nb133_fiyH2(%esp),%xmm4
        movapd nb133_fizH2(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb133_fixH2(%esp)
        movapd %xmm4,nb133_fiyH2(%esp)
        movapd %xmm7,nb133_fizH2(%esp)
        ## update j forces with water H2
        addpd  nb133_fjx(%esp),%xmm0
        addpd  nb133_fjy(%esp),%xmm1
        addpd  nb133_fjz(%esp),%xmm2
        movapd %xmm0,nb133_fjx(%esp)
        movapd %xmm1,nb133_fjy(%esp)
        movapd %xmm2,nb133_fjz(%esp)

        ## M interactions 
        movapd  nb133_rinvM(%esp),%xmm6
        movapd  %xmm6,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        mulpd   nb133_qqM(%esp),%xmm6   ## vcoul 
        mulpd   %xmm6,%xmm4   ## fscal
        addpd  nb133_vctot(%esp),%xmm6
        movapd %xmm6,nb133_vctot(%esp)

        movapd nb133_dxM(%esp),%xmm0
        movapd nb133_dyM(%esp),%xmm1
        movapd nb133_dzM(%esp),%xmm2
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        ## update H2 forces 
        movapd nb133_fixM(%esp),%xmm3
        movapd nb133_fiyM(%esp),%xmm4
        movapd nb133_fizM(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb133_fixM(%esp)
        movapd %xmm4,nb133_fiyM(%esp)
        movapd %xmm7,nb133_fizM(%esp)

        movl nb133_faction(%ebp),%edi
        ## update j forces 
        addpd  nb133_fjx(%esp),%xmm0
        addpd  nb133_fjy(%esp),%xmm1
        addpd  nb133_fjz(%esp),%xmm2
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
        subl $2,nb133_innerk(%esp)
        jl   _nb_kernel133_ia32_sse2.nb133_checksingle
        jmp  _nb_kernel133_ia32_sse2.nb133_unroll_loop
_nb_kernel133_ia32_sse2.nb133_checksingle: 
        movl  nb133_innerk(%esp),%edx
        andl  $1,%edx
        jnz  _nb_kernel133_ia32_sse2.nb133_dosingle
        jmp  _nb_kernel133_ia32_sse2.nb133_updateouterdata
_nb_kernel133_ia32_sse2.nb133_dosingle: 
        movl  nb133_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb133_innerjjnr(%esp)

        movl nb133_charge(%ebp),%esi     ## base of charge[] 

        xorpd %xmm3,%xmm3
        movlpd (%esi,%eax,8),%xmm3
        movapd %xmm3,%xmm4
        mulsd  nb133_iqM(%esp),%xmm3
        mulsd  nb133_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 

        movapd  %xmm3,nb133_qqM(%esp)
        movapd  %xmm4,nb133_qqH(%esp)

        movl nb133_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb133_vdwparam(%ebp),%esi
        shll %eax
        movl nb133_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movapd %xmm4,nb133_c6(%esp)
        movapd %xmm6,nb133_c12(%esp)

        movl nb133_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb133_ixO(%esp),%xmm4
        movapd nb133_iyO(%esp),%xmm5
        movapd nb133_izO(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb133_dxO(%esp)
        movapd %xmm5,nb133_dyO(%esp)
        movapd %xmm6,nb133_dzO(%esp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 
        movapd %xmm7,nb133_rsqO(%esp)

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb133_ixH1(%esp),%xmm4
        movapd nb133_iyH1(%esp),%xmm5
        movapd nb133_izH1(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb133_dxH1(%esp)
        movapd %xmm5,nb133_dyH1(%esp)
        movapd %xmm6,nb133_dzH1(%esp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm6
        addsd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb133_ixH2(%esp),%xmm3
        movapd nb133_iyH2(%esp),%xmm4
        movapd nb133_izH2(%esp),%xmm5

        ## calc dr 
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5

        ## store dr 
        movapd %xmm3,nb133_dxH2(%esp)
        movapd %xmm4,nb133_dyH2(%esp)
        movapd %xmm5,nb133_dzH2(%esp)
        ## square it 
        mulsd %xmm3,%xmm3
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        addsd %xmm4,%xmm5
        addsd %xmm3,%xmm5
        ## move ixM-izM to xmm2-xmm4  
        movapd nb133_iyM(%esp),%xmm3
        movapd nb133_izM(%esp),%xmm4
        subpd  %xmm1,%xmm3
        subpd  %xmm2,%xmm4
        movapd nb133_ixM(%esp),%xmm2
        subpd  %xmm0,%xmm2

        ## store dr 
        movapd %xmm2,nb133_dxM(%esp)
        movapd %xmm3,nb133_dyM(%esp)
        movapd %xmm4,nb133_dzM(%esp)
        ## square it 
        mulpd %xmm2,%xmm2
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        addpd %xmm3,%xmm4
        addpd %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## start with rsqH1 - put seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb133_three(%esp),%xmm1
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb133_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm6       ## rsq*lu*lu 
        movapd nb133_three(%esp),%xmm1
        subsd %xmm6,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb133_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb133_rinvH1(%esp)

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb133_three(%esp),%xmm1
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb133_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm5       ## rsq*lu*lu 
        movapd nb133_three(%esp),%xmm1
        subsd %xmm5,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb133_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb133_rinvH2(%esp)

        ## rsqM - seed in xmm2 
        cvtsd2ss %xmm4,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb133_three(%esp),%xmm1
        mulsd   %xmm4,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb133_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm4       ## rsq*lu*lu 
        movapd nb133_three(%esp),%xmm1
        subsd %xmm4,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb133_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb133_rinvM(%esp)

        ## rsqO - put seed in xmm2 
        cvtsd2ss %xmm7,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movsd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movsd  nb133_three(%esp),%xmm4
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb133_half(%esp),%xmm4   ## iter1 ( new lu) 

        movsd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm7       ## rsq*lu*lu 
        movsd nb133_three(%esp),%xmm4
        subsd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb133_half(%esp),%xmm4   ## rinv 
        movsd  %xmm4,%xmm7      ## rinvO in xmm7 

        movsd nb133_rsqO(%esp),%xmm4
        movapd %xmm7,%xmm0
        ## LJ table interaction.
        mulsd %xmm7,%xmm4       ## xmm4=r 
        mulsd nb133_tsc(%esp),%xmm4

        cvttsd2si %xmm4,%ebx    ## mm6 = lu idx 
        cvtsi2sd %ebx,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $3,%ebx

        movl nb133_VFtab(%ebp),%esi

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
        mulsd  nb133_two(%esp),%xmm7    ## two*Heps2 
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movsd nb133_c6(%esp),%xmm4
        mulsd  %xmm4,%xmm7       ## fijD 
        mulsd  %xmm4,%xmm5       ## Vvdw6 

        ## put scalar force on stack Update Vvdwtot directly 
        addsd  nb133_Vvdwtot(%esp),%xmm5
        xorpd  %xmm3,%xmm3
        mulsd  nb133_tsc(%esp),%xmm7
        subsd  %xmm7,%xmm3
        movsd %xmm3,nb133_fstmp(%esp)
        movsd %xmm5,nb133_Vvdwtot(%esp)

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
        mulsd  nb133_two(%esp),%xmm7    ## two*Heps2 
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movsd nb133_c12(%esp),%xmm4
        mulsd  %xmm4,%xmm7
        mulsd  %xmm4,%xmm5

        addsd  nb133_Vvdwtot(%esp),%xmm5
        movsd nb133_fstmp(%esp),%xmm3
        mulsd  nb133_tsc(%esp),%xmm7
        subsd  %xmm7,%xmm3
        movsd %xmm5,nb133_Vvdwtot(%esp)

        mulsd  %xmm0,%xmm3


        movsd nb133_dxO(%esp),%xmm0
        movsd nb133_dyO(%esp),%xmm1
        movsd nb133_dzO(%esp),%xmm2

        movl   nb133_faction(%ebp),%edi
        mulsd  %xmm3,%xmm0
        mulsd  %xmm3,%xmm1
        mulsd  %xmm3,%xmm2

        ## update O forces 
        movapd nb133_fixO(%esp),%xmm3
        movapd nb133_fiyO(%esp),%xmm4
        movapd nb133_fizO(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movsd %xmm3,nb133_fixO(%esp)
        movsd %xmm4,nb133_fiyO(%esp)
        movsd %xmm7,nb133_fizO(%esp)
        ## update j forces with water O 
        movsd %xmm0,nb133_fjx(%esp)
        movsd %xmm1,nb133_fjy(%esp)
        movsd %xmm2,nb133_fjz(%esp)

        ## H1 interactions
        movsd  nb133_rinvH1(%esp),%xmm6
        movsd  %xmm6,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        mulsd   nb133_qqH(%esp),%xmm6   ## vcoul 
        mulsd   %xmm6,%xmm4   ## fscal
        addsd  nb133_vctot(%esp),%xmm6
        movsd %xmm6,nb133_vctot(%esp)

        movapd nb133_dxH1(%esp),%xmm0
        movapd nb133_dyH1(%esp),%xmm1
        movapd nb133_dzH1(%esp),%xmm2
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb133_fixH1(%esp),%xmm3
        movapd nb133_fiyH1(%esp),%xmm4
        movapd nb133_fizH1(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movsd %xmm3,nb133_fixH1(%esp)
        movsd %xmm4,nb133_fiyH1(%esp)
        movsd %xmm7,nb133_fizH1(%esp)
        ## update j forces with water H1 
        addsd  nb133_fjx(%esp),%xmm0
        addsd  nb133_fjy(%esp),%xmm1
        addsd  nb133_fjz(%esp),%xmm2
        movsd %xmm0,nb133_fjx(%esp)
        movsd %xmm1,nb133_fjy(%esp)
        movsd %xmm2,nb133_fjz(%esp)

        ## H2 interactions 
        movsd  nb133_rinvH2(%esp),%xmm6
        movsd  %xmm6,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        mulsd   nb133_qqH(%esp),%xmm6   ## vcoul 
        mulsd   %xmm6,%xmm4   ## fscal
        addsd  nb133_vctot(%esp),%xmm6
        movsd %xmm6,nb133_vctot(%esp)

        movapd nb133_dxH2(%esp),%xmm0
        movapd nb133_dyH2(%esp),%xmm1
        movapd nb133_dzH2(%esp),%xmm2
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H2 forces 
        movapd nb133_fixH2(%esp),%xmm3
        movapd nb133_fiyH2(%esp),%xmm4
        movapd nb133_fizH2(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movsd %xmm3,nb133_fixH2(%esp)
        movsd %xmm4,nb133_fiyH2(%esp)
        movsd %xmm7,nb133_fizH2(%esp)
        ## update j forces with water H2 
        addsd  nb133_fjx(%esp),%xmm0
        addsd  nb133_fjy(%esp),%xmm1
        addsd  nb133_fjz(%esp),%xmm2
        movsd %xmm0,nb133_fjx(%esp)
        movsd %xmm1,nb133_fjy(%esp)
        movsd %xmm2,nb133_fjz(%esp)

        ## M interactions 
        movsd  nb133_rinvM(%esp),%xmm6
        movsd  %xmm6,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        mulsd   nb133_qqM(%esp),%xmm6   ## vcoul 
        mulsd   %xmm6,%xmm4   ## fscal
        addsd  nb133_vctot(%esp),%xmm6
        movsd %xmm6,nb133_vctot(%esp)

        movapd nb133_dxM(%esp),%xmm0
        movapd nb133_dyM(%esp),%xmm1
        movapd nb133_dzM(%esp),%xmm2
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update M forces 
        movapd nb133_fixM(%esp),%xmm3
        movapd nb133_fiyM(%esp),%xmm4
        movapd nb133_fizM(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movsd %xmm3,nb133_fixM(%esp)
        movsd %xmm4,nb133_fiyM(%esp)
        movsd %xmm7,nb133_fizM(%esp)

        movl nb133_faction(%ebp),%edi
        ## update j forces 
        addsd  nb133_fjx(%esp),%xmm0
        addsd  nb133_fjy(%esp),%xmm1
        addsd  nb133_fjz(%esp),%xmm2
        movlpd (%edi,%eax,8),%xmm3
        movlpd 8(%edi,%eax,8),%xmm4
        movlpd 16(%edi,%eax,8),%xmm5
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        movlpd %xmm3,(%edi,%eax,8)
        movlpd %xmm4,8(%edi,%eax,8)
        movlpd %xmm5,16(%edi,%eax,8)

_nb_kernel133_ia32_sse2.nb133_updateouterdata: 
        movl  nb133_ii3(%esp),%ecx
        movl  nb133_faction(%ebp),%edi
        movl  nb133_fshift(%ebp),%esi
        movl  nb133_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movapd nb133_fixO(%esp),%xmm0
        movapd nb133_fiyO(%esp),%xmm1
        movapd nb133_fizO(%esp),%xmm2

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
        movapd nb133_fixH1(%esp),%xmm0
        movapd nb133_fiyH1(%esp),%xmm1
        movapd nb133_fizH1(%esp),%xmm2

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
        movapd nb133_fixH2(%esp),%xmm0
        movapd nb133_fiyH2(%esp),%xmm1
        movapd nb133_fizH2(%esp),%xmm2

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
        movapd nb133_fixM(%esp),%xmm0
        movapd nb133_fiyM(%esp),%xmm1
        movapd nb133_fizM(%esp),%xmm2

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
        movl nb133_n(%esp),%esi
        ## get group index for i particle 
        movl  nb133_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb133_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb133_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb133_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb133_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

       ## finish if last 
        movl nb133_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel133_ia32_sse2.nb133_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb133_n(%esp)
        jmp _nb_kernel133_ia32_sse2.nb133_outer
_nb_kernel133_ia32_sse2.nb133_outerend: 
        ## check if more outer neighborlists remain
        movl  nb133_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel133_ia32_sse2.nb133_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel133_ia32_sse2.nb133_threadloop
_nb_kernel133_ia32_sse2.nb133_end: 
        emms

        movl nb133_nouter(%esp),%eax
        movl nb133_ninner(%esp),%ebx
        movl nb133_outeriter(%ebp),%ecx
        movl nb133_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb133_salign(%esp),%eax
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




.globl nb_kernel133nf_ia32_sse2
.globl _nb_kernel133nf_ia32_sse2
nb_kernel133nf_ia32_sse2:       
_nb_kernel133nf_ia32_sse2:      
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
        ## bottom of stack is cache-aligned for sse2 use 
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
.set nb133nf_iqH, 192
.set nb133nf_iqM, 208
.set nb133nf_qqH, 224
.set nb133nf_qqM, 240
.set nb133nf_c6, 256
.set nb133nf_c12, 272
.set nb133nf_tsc, 288
.set nb133nf_vctot, 304
.set nb133nf_Vvdwtot, 320
.set nb133nf_half, 336
.set nb133nf_three, 352
.set nb133nf_two, 368
.set nb133nf_rinvH1, 384
.set nb133nf_rinvH2, 400
.set nb133nf_rinvM, 416
.set nb133nf_krsqH1, 432
.set nb133nf_krsqH2, 448
.set nb133nf_krsqM, 464
.set nb133nf_rsqO, 512
.set nb133nf_is3, 528
.set nb133nf_ii3, 532
.set nb133nf_ntia, 536
.set nb133nf_innerjjnr, 540
.set nb133nf_innerk, 544
.set nb133nf_n, 548
.set nb133nf_nn1, 552
.set nb133nf_nri, 556
.set nb133nf_nouter, 560
.set nb133nf_ninner, 564
.set nb133nf_salign, 568
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $572,%esp          ## local stack space 
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
        movsd (%eax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb133nf_tsc(%esp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb133nf_half(%esp)
        movl %ebx,nb133nf_half+4(%esp)
        movsd nb133nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb133nf_half(%esp)
        movapd %xmm2,nb133nf_two(%esp)
        movapd %xmm3,nb133nf_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb133nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb133nf_charge(%ebp),%edx
        movsd 8(%edx,%ebx,8),%xmm3
        movsd 24(%edx,%ebx,8),%xmm4
        movl nb133nf_p_facel(%ebp),%esi
        movsd (%esi),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb133nf_iqH(%esp)
        movapd %xmm4,nb133nf_iqM(%esp)

        movl  nb133nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb133nf_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb133nf_ntia(%esp)
_nb_kernel133nf_ia32_sse2.nb133nf_threadloop: 
        movl  nb133nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel133nf_ia32_sse2.nb133nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel133nf_ia32_sse2.nb133nf_spinlock

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
        jg  _nb_kernel133nf_ia32_sse2.nb133nf_outerstart
        jmp _nb_kernel133nf_ia32_sse2.nb133nf_end

_nb_kernel133nf_ia32_sse2.nb133nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb133nf_nouter(%esp),%ebx
        movl %ebx,nb133nf_nouter(%esp)

_nb_kernel133nf_ia32_sse2.nb133nf_outer: 
        movl  nb133nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb133nf_is3(%esp)            ## store is3 

        movl  nb133nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb133nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        movapd %xmm0,%xmm6
        movapd %xmm1,%xmm7

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb133nf_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb133nf_ii3(%esp)

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
        movapd %xmm3,nb133nf_ixO(%esp)
        movapd %xmm4,nb133nf_iyO(%esp)
        movapd %xmm5,nb133nf_izO(%esp)
        movapd %xmm6,nb133nf_ixH1(%esp)
        movapd %xmm7,nb133nf_iyH1(%esp)

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
        movapd %xmm6,nb133nf_izH1(%esp)
        movapd %xmm0,nb133nf_ixH2(%esp)
        movapd %xmm1,nb133nf_iyH2(%esp)
        movapd %xmm2,nb133nf_izH2(%esp)
        movapd %xmm3,nb133nf_ixM(%esp)
        movapd %xmm4,nb133nf_iyM(%esp)
        movapd %xmm5,nb133nf_izM(%esp)

        ## clear vctot
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb133nf_vctot(%esp)
        movapd %xmm4,nb133nf_Vvdwtot(%esp)

        movl  nb133nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb133nf_pos(%ebp),%esi
        movl  nb133nf_faction(%ebp),%edi
        movl  nb133nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb133nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb133nf_ninner(%esp),%ecx
        movl  %ecx,nb133nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb133nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel133nf_ia32_sse2.nb133nf_unroll_loop
        jmp   _nb_kernel133nf_ia32_sse2.nb133nf_checksingle
_nb_kernel133nf_ia32_sse2.nb133nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb133nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb133nf_innerjjnr(%esp)                 ## advance pointer (unrolled 2) 

        movl nb133nf_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb133nf_iqM(%esp),%xmm3
        mulpd  nb133nf_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movapd  %xmm3,nb133nf_qqM(%esp)
        movapd  %xmm4,nb133nf_qqH(%esp)

        movl nb133nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb133nf_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb133nf_ntia(%esp),%edi
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
        movapd %xmm4,nb133nf_c6(%esp)
        movapd %xmm6,nb133nf_c12(%esp)

        movl nb133nf_pos(%ebp),%esi        ## base of pos[] 

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
        movapd nb133nf_ixO(%esp),%xmm4
        movapd nb133nf_iyO(%esp),%xmm5
        movapd nb133nf_izO(%esp),%xmm6

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
        movapd nb133nf_ixH1(%esp),%xmm4
        movapd nb133nf_iyH1(%esp),%xmm5
        movapd nb133nf_izH1(%esp),%xmm6

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
        movapd nb133nf_ixH2(%esp),%xmm3
        movapd nb133nf_iyH2(%esp),%xmm4
        movapd nb133nf_izH2(%esp),%xmm5

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
        movapd nb133nf_iyM(%esp),%xmm3
        movapd nb133nf_izM(%esp),%xmm4
        subpd  %xmm1,%xmm3
        subpd  %xmm2,%xmm4
        movapd nb133nf_ixM(%esp),%xmm2
        subpd  %xmm0,%xmm2

        ## square it 
        mulpd %xmm2,%xmm2
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        addpd %xmm3,%xmm4
        addpd %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 
        movapd %xmm7,nb133nf_rsqO(%esp)

        ## start with rsqH1 - put seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb133nf_three(%esp),%xmm1
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulpd   nb133nf_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulpd %xmm1,%xmm1       ## lu*lu 
        mulpd %xmm1,%xmm6       ## rsq*lu*lu 
        movapd nb133nf_three(%esp),%xmm1
        subpd %xmm6,%xmm1       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulpd nb133nf_half(%esp),%xmm1   ## rinv 
        movapd  %xmm1,nb133nf_rinvH1(%esp)

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb133nf_three(%esp),%xmm1
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulpd   nb133nf_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulpd %xmm1,%xmm1       ## lu*lu 
        mulpd %xmm1,%xmm5       ## rsq*lu*lu 
        movapd nb133nf_three(%esp),%xmm1
        subpd %xmm5,%xmm1       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulpd nb133nf_half(%esp),%xmm1   ## rinv 
        movapd  %xmm1,nb133nf_rinvH2(%esp)

        ## rsqM - seed in xmm2 
        cvtpd2ps %xmm4,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb133nf_three(%esp),%xmm1
        mulpd   %xmm4,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulpd   nb133nf_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulpd %xmm1,%xmm1       ## lu*lu 
        mulpd %xmm1,%xmm4       ## rsq*lu*lu 
        movapd nb133nf_three(%esp),%xmm1
        subpd %xmm4,%xmm1       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulpd nb133nf_half(%esp),%xmm1   ## rinv 
        movapd  %xmm1,nb133nf_rinvM(%esp)


        ## rsqO - put seed in xmm2 
        cvtpd2ps %xmm7,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb133nf_three(%esp),%xmm4
        mulpd   %xmm7,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb133nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb133nf_three(%esp),%xmm4
        subpd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb133nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvO in xmm7 



        movapd nb133nf_rsqO(%esp),%xmm4
        movapd %xmm7,%xmm0
        ## LJ table interaction.
        mulpd %xmm7,%xmm4       ## xmm4=r 
        mulpd nb133nf_tsc(%esp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $3,%mm6           ## idx *= 8 

        movl nb133nf_VFtab(%ebp),%esi
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

        movapd nb133nf_c6(%esp),%xmm4
        mulpd  %xmm4,%xmm5       ## Vvdw6 

        ## Update Vvdwtot directly 
        addpd  nb133nf_Vvdwtot(%esp),%xmm5
        movapd %xmm5,nb133nf_Vvdwtot(%esp)

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

        movapd nb133nf_c12(%esp),%xmm4
        mulpd  %xmm4,%xmm5

        addpd  nb133nf_Vvdwtot(%esp),%xmm5
        movapd %xmm5,nb133nf_Vvdwtot(%esp)

        ## H1/H2/M interactions 
        movapd  nb133nf_rinvH1(%esp),%xmm6
        addpd   nb133nf_rinvH2(%esp),%xmm6
        movapd  nb133nf_rinvM(%esp),%xmm7
        mulpd   nb133nf_qqH(%esp),%xmm6
        mulpd   nb133nf_qqM(%esp),%xmm7
        addpd   %xmm7,%xmm6
        addpd   nb133nf_vctot(%esp),%xmm6
        movapd  %xmm6,nb133nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $2,nb133nf_innerk(%esp)
        jl   _nb_kernel133nf_ia32_sse2.nb133nf_checksingle
        jmp  _nb_kernel133nf_ia32_sse2.nb133nf_unroll_loop
_nb_kernel133nf_ia32_sse2.nb133nf_checksingle: 
        movl  nb133nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz  _nb_kernel133nf_ia32_sse2.nb133nf_dosingle
        jmp  _nb_kernel133nf_ia32_sse2.nb133nf_updateouterdata
_nb_kernel133nf_ia32_sse2.nb133nf_dosingle: 
        movl  nb133nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb133nf_innerjjnr(%esp)

        movl nb133nf_charge(%ebp),%esi     ## base of charge[] 

        xorpd %xmm3,%xmm3
        movlpd (%esi,%eax,8),%xmm3
        movapd %xmm3,%xmm4
        mulsd  nb133nf_iqM(%esp),%xmm3
        mulsd  nb133nf_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 

        movapd  %xmm3,nb133nf_qqM(%esp)
        movapd  %xmm4,nb133nf_qqH(%esp)

        movl nb133nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb133nf_vdwparam(%ebp),%esi
        shll %eax
        movl nb133nf_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movapd %xmm4,nb133nf_c6(%esp)
        movapd %xmm6,nb133nf_c12(%esp)

        movl nb133nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb133nf_ixO(%esp),%xmm4
        movapd nb133nf_iyO(%esp),%xmm5
        movapd nb133nf_izO(%esp),%xmm6

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
        movapd %xmm7,nb133nf_rsqO(%esp)

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb133nf_ixH1(%esp),%xmm4
        movapd nb133nf_iyH1(%esp),%xmm5
        movapd nb133nf_izH1(%esp),%xmm6

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
        movapd nb133nf_ixH2(%esp),%xmm3
        movapd nb133nf_iyH2(%esp),%xmm4
        movapd nb133nf_izH2(%esp),%xmm5

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
        movapd nb133nf_iyM(%esp),%xmm3
        movapd nb133nf_izM(%esp),%xmm4
        subpd  %xmm1,%xmm3
        subpd  %xmm2,%xmm4
        movapd nb133nf_ixM(%esp),%xmm2
        subpd  %xmm0,%xmm2

        ## square it 
        mulpd %xmm2,%xmm2
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        addpd %xmm3,%xmm4
        addpd %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## start with rsqH1 - put seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb133nf_three(%esp),%xmm1
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb133nf_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm6       ## rsq*lu*lu 
        movapd nb133nf_three(%esp),%xmm1
        subsd %xmm6,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb133nf_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb133nf_rinvH1(%esp)

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb133nf_three(%esp),%xmm1
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb133nf_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm5       ## rsq*lu*lu 
        movapd nb133nf_three(%esp),%xmm1
        subsd %xmm5,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb133nf_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb133nf_rinvH2(%esp)

        ## rsqM - seed in xmm2 
        cvtsd2ss %xmm4,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb133nf_three(%esp),%xmm1
        mulsd   %xmm4,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb133nf_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm4       ## rsq*lu*lu 
        movapd nb133nf_three(%esp),%xmm1
        subsd %xmm4,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb133nf_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb133nf_rinvM(%esp)

        ## rsqO - put seed in xmm2 
        cvtsd2ss %xmm7,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movsd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movsd  nb133nf_three(%esp),%xmm4
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb133nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movsd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm7       ## rsq*lu*lu 
        movsd nb133nf_three(%esp),%xmm4
        subsd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb133nf_half(%esp),%xmm4   ## rinv 
        movsd  %xmm4,%xmm7      ## rinvO in xmm7 

        movsd nb133nf_rsqO(%esp),%xmm4
        movapd %xmm7,%xmm0
        ## LJ table interaction.
        mulsd %xmm7,%xmm4       ## xmm4=r 
        mulsd nb133nf_tsc(%esp),%xmm4

        cvttsd2si %xmm4,%ebx    ## mm6 = lu idx 
        cvtsi2sd %ebx,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $3,%ebx

        movl nb133nf_VFtab(%ebp),%esi

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

        movsd nb133nf_c6(%esp),%xmm4
        mulsd  %xmm4,%xmm5       ## Vvdw6 

        ## put scalar force on stack Update Vvdwtot directly 
        addsd  nb133nf_Vvdwtot(%esp),%xmm5
        movsd %xmm5,nb133nf_Vvdwtot(%esp)

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

        movsd nb133nf_c12(%esp),%xmm4
        mulsd  %xmm4,%xmm5

        addsd  nb133nf_Vvdwtot(%esp),%xmm5
        movsd %xmm5,nb133nf_Vvdwtot(%esp)

        ## H1/H2/M interactions 
        movsd  nb133nf_rinvH1(%esp),%xmm6
        addsd  nb133nf_rinvH2(%esp),%xmm6
        movsd  nb133nf_rinvM(%esp),%xmm7
        mulsd  nb133nf_qqH(%esp),%xmm6
        mulsd  nb133nf_qqM(%esp),%xmm7
        addsd  %xmm7,%xmm6
        addsd  nb133nf_vctot(%esp),%xmm6
        movsd  %xmm6,nb133nf_vctot(%esp)

_nb_kernel133nf_ia32_sse2.nb133nf_updateouterdata: 
        ## get n from stack
        movl nb133nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb133nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb133nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb133nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb133nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb133nf_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

       ## finish if last 
        movl nb133nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel133nf_ia32_sse2.nb133nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb133nf_n(%esp)
        jmp _nb_kernel133nf_ia32_sse2.nb133nf_outer
_nb_kernel133nf_ia32_sse2.nb133nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb133nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel133nf_ia32_sse2.nb133nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel133nf_ia32_sse2.nb133nf_threadloop
_nb_kernel133nf_ia32_sse2.nb133nf_end: 
        emms

        movl nb133nf_nouter(%esp),%eax
        movl nb133nf_ninner(%esp),%ebx
        movl nb133nf_outeriter(%ebp),%ecx
        movl nb133nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb133nf_salign(%esp),%eax
        addl %eax,%esp
        addl $572,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret


