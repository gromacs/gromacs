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



.globl nb_kernel313_ia32_sse2
.globl _nb_kernel313_ia32_sse2
nb_kernel313_ia32_sse2: 
_nb_kernel313_ia32_sse2:        
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
        ## bottom of stack is cache-aligned for sse2 use 
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
.set nb313_rinvsqO, 448
.set nb313_rinvH1, 464
.set nb313_rinvH2, 480
.set nb313_rinvM, 496
.set nb313_rO, 512
.set nb313_rH1, 528
.set nb313_rH2, 544
.set nb313_rM, 560
.set nb313_tsc, 576
.set nb313_two, 592
.set nb313_c6, 608
.set nb313_c12, 624
.set nb313_six, 640
.set nb313_twelve, 656
.set nb313_vctot, 672
.set nb313_Vvdwtot, 688
.set nb313_fixO, 704
.set nb313_fiyO, 720
.set nb313_fizO, 736
.set nb313_fixH1, 752
.set nb313_fiyH1, 768
.set nb313_fizH1, 784
.set nb313_fixH2, 800
.set nb313_fiyH2, 816
.set nb313_fizH2, 832
.set nb313_fixM, 848
.set nb313_fiyM, 864
.set nb313_fizM, 880
.set nb313_fjx, 896
.set nb313_fjy, 912
.set nb313_fjz, 928
.set nb313_half, 944
.set nb313_three, 960
.set nb313_is3, 976
.set nb313_ii3, 980
.set nb313_ntia, 984
.set nb313_innerjjnr, 988
.set nb313_innerk, 992
.set nb313_n, 996
.set nb313_nn1, 1000
.set nb313_nri, 1004
.set nb313_nouter, 1008
.set nb313_ninner, 1012
.set nb313_salign, 1016
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


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb313_half(%esp)
        movl %ebx,nb313_half+4(%esp)
        movsd nb313_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm3,%xmm4
        addpd  %xmm4,%xmm4      ## 6.0
        movapd %xmm4,%xmm5
        addpd  %xmm5,%xmm5      ## 12.0
        movapd %xmm1,nb313_half(%esp)
        movapd %xmm2,nb313_two(%esp)
        movapd %xmm3,nb313_three(%esp)
        movapd %xmm4,nb313_six(%esp)
        movapd %xmm5,nb313_twelve(%esp)
        movl nb313_p_tabscale(%ebp),%eax
        movsd (%eax),%xmm5
        shufpd $0,%xmm5,%xmm5
        movapd %xmm5,nb313_tsc(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb313_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb313_charge(%ebp),%edx
        movsd 8(%edx,%ebx,8),%xmm3
        movsd 24(%edx,%ebx,8),%xmm4
        movl nb313_p_facel(%ebp),%esi
        movsd (%esi),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb313_iqH(%esp)
        movapd %xmm4,nb313_iqM(%esp)

        movl  nb313_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb313_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb313_ntia(%esp)
_nb_kernel313_ia32_sse2.nb313_threadloop: 
        movl  nb313_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel313_ia32_sse2.nb313_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel313_ia32_sse2.nb313_spinlock

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
        jg  _nb_kernel313_ia32_sse2.nb313_outerstart
        jmp _nb_kernel313_ia32_sse2.nb313_end

_nb_kernel313_ia32_sse2.nb313_outerstart: 
        ## ebx contains number of outer iterations
        addl nb313_nouter(%esp),%ebx
        movl %ebx,nb313_nouter(%esp)

_nb_kernel313_ia32_sse2.nb313_outer: 
        movl  nb313_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb313_is3(%esp)      ## store is3 

        movl  nb313_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb313_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        movapd %xmm0,%xmm6
        movapd %xmm1,%xmm7

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb313_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb313_ii3(%esp)

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
        movapd %xmm3,nb313_ixO(%esp)
        movapd %xmm4,nb313_iyO(%esp)
        movapd %xmm5,nb313_izO(%esp)
        movapd %xmm6,nb313_ixH1(%esp)
        movapd %xmm7,nb313_iyH1(%esp)

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
        movapd %xmm6,nb313_izH1(%esp)
        movapd %xmm0,nb313_ixH2(%esp)
        movapd %xmm1,nb313_iyH2(%esp)
        movapd %xmm2,nb313_izH2(%esp)
        movapd %xmm3,nb313_ixM(%esp)
        movapd %xmm4,nb313_iyM(%esp)
        movapd %xmm5,nb313_izM(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb313_vctot(%esp)
        movapd %xmm4,nb313_Vvdwtot(%esp)
        movapd %xmm4,nb313_fixO(%esp)
        movapd %xmm4,nb313_fiyO(%esp)
        movapd %xmm4,nb313_fizO(%esp)
        movapd %xmm4,nb313_fixH1(%esp)
        movapd %xmm4,nb313_fiyH1(%esp)
        movapd %xmm4,nb313_fizH1(%esp)
        movapd %xmm4,nb313_fixH2(%esp)
        movapd %xmm4,nb313_fiyH2(%esp)
        movapd %xmm4,nb313_fizH2(%esp)
        movapd %xmm4,nb313_fixM(%esp)
        movapd %xmm4,nb313_fiyM(%esp)
        movapd %xmm4,nb313_fizM(%esp)

        movl  nb313_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx     ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb313_pos(%ebp),%esi
        movl  nb313_faction(%ebp),%edi
        movl  nb313_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb313_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb313_ninner(%esp),%ecx
        movl  %ecx,nb313_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb313_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel313_ia32_sse2.nb313_unroll_loop
        jmp   _nb_kernel313_ia32_sse2.nb313_checksingle
_nb_kernel313_ia32_sse2.nb313_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb313_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb313_innerjjnr(%esp)             ## advance pointer (unrolled 2) 

        movl nb313_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb313_iqM(%esp),%xmm3
        mulpd  nb313_iqH(%esp),%xmm4
        movapd  %xmm3,nb313_qqM(%esp)
        movapd  %xmm4,nb313_qqH(%esp)

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movl nb313_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb313_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb313_ntia(%esp),%edi
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
        movapd %xmm4,nb313_c6(%esp)
        movapd %xmm6,nb313_c12(%esp)

        movl nb313_pos(%ebp),%esi        ## base of pos[] 

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
        movapd nb313_ixO(%esp),%xmm4
        movapd nb313_iyO(%esp),%xmm5
        movapd nb313_izO(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb313_dxO(%esp)
        movapd %xmm5,nb313_dyO(%esp)
        movapd %xmm6,nb313_dzO(%esp)
        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb313_ixH1(%esp),%xmm4
        movapd nb313_iyH1(%esp),%xmm5
        movapd nb313_izH1(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb313_dxH1(%esp)
        movapd %xmm5,nb313_dyH1(%esp)
        movapd %xmm6,nb313_dzH1(%esp)
        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm6
        addpd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb313_ixH2(%esp),%xmm3
        movapd nb313_iyH2(%esp),%xmm4
        movapd nb313_izH2(%esp),%xmm5

        ## calc dr 
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5

        ## store dr 
        movapd %xmm3,nb313_dxH2(%esp)
        movapd %xmm4,nb313_dyH2(%esp)
        movapd %xmm5,nb313_dzH2(%esp)
        ## square it 
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        addpd %xmm4,%xmm5
        addpd %xmm3,%xmm5
        ## move ixM-izM to xmm2-xmm4  
        movapd nb313_iyM(%esp),%xmm3
        movapd nb313_izM(%esp),%xmm4
        subpd  %xmm1,%xmm3
        subpd  %xmm2,%xmm4
        movapd nb313_ixM(%esp),%xmm2
        subpd  %xmm0,%xmm2

        ## store dr 
        movapd %xmm2,nb313_dxM(%esp)
        movapd %xmm3,nb313_dyM(%esp)
        movapd %xmm4,nb313_dzM(%esp)
        ## square it 
        mulpd %xmm2,%xmm2
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        addpd %xmm3,%xmm4
        addpd %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## 1/x for O - rsqO is in xmm7
        cvtpd2ps %xmm7,%xmm2
        movapd   %xmm7,%xmm3
        rcpps    %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2
        movapd   nb313_two(%esp),%xmm1
        movapd   %xmm1,%xmm0
        mulpd   %xmm2,%xmm7
        subpd   %xmm7,%xmm1
        mulpd   %xmm1,%xmm2 ## iter1 
        mulpd   %xmm2,%xmm3
        subpd   %xmm3,%xmm0
        mulpd   %xmm2,%xmm0 ## xmm0=rinvsq
        movapd  %xmm0,nb313_rinvsqO(%esp)

        ## rsqH1 - seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb313_three(%esp),%xmm0
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulpd   nb313_half(%esp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm6,%xmm2
        movapd %xmm0,%xmm3
        mulpd %xmm0,%xmm0       ## lu*lu 
        mulpd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb313_three(%esp),%xmm0
        subpd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulpd nb313_half(%esp),%xmm0   ## rinv 
        movapd %xmm0,nb313_rinvH1(%esp)         ## rinvH1 
        mulpd  %xmm0,%xmm6
        movapd %xmm6,nb313_rH1(%esp)    ## rH1 

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb313_three(%esp),%xmm0
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulpd   nb313_half(%esp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm5,%xmm2
        movapd %xmm0,%xmm3
        mulpd %xmm0,%xmm0       ## lu*lu 
        mulpd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb313_three(%esp),%xmm0
        subpd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulpd nb313_half(%esp),%xmm0   ## rinv 
        movapd %xmm0,nb313_rinvH2(%esp)   ## rinv 
        mulpd %xmm0,%xmm5
        movapd %xmm5,nb313_rH2(%esp)   ## r 

        ## rsqM - seed in xmm2 
        cvtpd2ps %xmm4,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb313_three(%esp),%xmm0
        mulpd   %xmm4,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulpd   nb313_half(%esp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm4,%xmm2
        movapd %xmm0,%xmm3
        mulpd %xmm0,%xmm0       ## lu*lu 
        mulpd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb313_three(%esp),%xmm0
        subpd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulpd nb313_half(%esp),%xmm0   ## rinv 
        movapd %xmm0,nb313_rinvM(%esp)   ## rinv 
        mulpd %xmm0,%xmm4
        movapd %xmm4,nb313_rM(%esp)   ## r 

        ## do O interactions 
        movapd nb313_rinvsqO(%esp),%xmm0
        movapd  %xmm0,%xmm1
        mulpd   %xmm1,%xmm1 ## rinv4
        mulpd   %xmm0,%xmm1 ##rinvsix
        movapd  %xmm1,%xmm2
        mulpd   %xmm2,%xmm2 ## rinvtwelve
        mulpd  nb313_c6(%esp),%xmm1
        mulpd  nb313_c12(%esp),%xmm2
        movapd %xmm2,%xmm3
        subpd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addpd  nb313_Vvdwtot(%esp),%xmm3
        mulpd  nb313_six(%esp),%xmm1
        mulpd  nb313_twelve(%esp),%xmm2
        subpd  %xmm1,%xmm2
        mulpd  %xmm0,%xmm2
        movapd %xmm2,%xmm4 ## total fsO 
        movapd %xmm3,nb313_Vvdwtot(%esp)

        movapd nb313_dxO(%esp),%xmm0
        movapd nb313_dyO(%esp),%xmm1
        movapd nb313_dzO(%esp),%xmm2
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2
        ## tx in xmm0-xmm2 

        ## update O forces 
        movapd nb313_fixO(%esp),%xmm3
        movapd nb313_fiyO(%esp),%xmm4
        movapd nb313_fizO(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb313_fixO(%esp)
        movapd %xmm4,nb313_fiyO(%esp)
        movapd %xmm7,nb313_fizO(%esp)
        ## update j forces with water O 
        movapd %xmm0,nb313_fjx(%esp)
        movapd %xmm1,nb313_fjy(%esp)
        movapd %xmm2,nb313_fjz(%esp)

        ## Done with O interactions - now H1! 
        movapd nb313_rH1(%esp),%xmm7
        mulpd nb313_tsc(%esp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        movd %eax,%mm0
        movd %ebx,%mm1

        pslld $2,%mm6           ## idx *= 4 
        movl nb313_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

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
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb313_two(%esp),%xmm7    ## two*Heps2 
        movapd nb313_qqH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 
        ## increment vcoul 
        xorpd  %xmm4,%xmm4
        addpd  nb313_vctot(%esp),%xmm5
        mulpd  nb313_rinvH1(%esp),%xmm3
        movapd %xmm5,nb313_vctot(%esp)
        mulpd  nb313_tsc(%esp),%xmm3
        subpd %xmm3,%xmm4

        movapd nb313_dxH1(%esp),%xmm0
        movapd nb313_dyH1(%esp),%xmm1
        movapd nb313_dzH1(%esp),%xmm2
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb313_fixH1(%esp),%xmm3
        movapd nb313_fiyH1(%esp),%xmm4
        movapd nb313_fizH1(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb313_fixH1(%esp)
        movapd %xmm4,nb313_fiyH1(%esp)
        movapd %xmm7,nb313_fizH1(%esp)
        ## update j forces with water H1 
        addpd  nb313_fjx(%esp),%xmm0
        addpd  nb313_fjy(%esp),%xmm1
        addpd  nb313_fjz(%esp),%xmm2
        movapd %xmm0,nb313_fjx(%esp)
        movapd %xmm1,nb313_fjy(%esp)
        movapd %xmm2,nb313_fjz(%esp)

        ## H2 interactions 
        movapd nb313_rH2(%esp),%xmm7
        mulpd   nb313_tsc(%esp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb313_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

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
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb313_two(%esp),%xmm7    ## two*Heps2 
        movapd nb313_qqH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 
        ## increment vcoul 
        xorpd  %xmm4,%xmm4
        addpd  nb313_vctot(%esp),%xmm5
        mulpd  nb313_rinvH2(%esp),%xmm3
        movapd %xmm5,nb313_vctot(%esp)
        mulpd  nb313_tsc(%esp),%xmm3
        subpd  %xmm3,%xmm4

        movapd nb313_dxH2(%esp),%xmm0
        movapd nb313_dyH2(%esp),%xmm1
        movapd nb313_dzH2(%esp),%xmm2
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        movd %mm0,%eax
        movd %mm1,%ebx

        ## update H2 forces 
        movapd nb313_fixH2(%esp),%xmm3
        movapd nb313_fiyH2(%esp),%xmm4
        movapd nb313_fizH2(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb313_fixH2(%esp)
        movapd %xmm4,nb313_fiyH2(%esp)
        movapd %xmm7,nb313_fizH2(%esp)
        ## update j forces with water H1 
        addpd  nb313_fjx(%esp),%xmm0
        addpd  nb313_fjy(%esp),%xmm1
        addpd  nb313_fjz(%esp),%xmm2
        movapd %xmm0,nb313_fjx(%esp)
        movapd %xmm1,nb313_fjy(%esp)
        movapd %xmm2,nb313_fjz(%esp)

        ## M interactions 
        movapd nb313_rM(%esp),%xmm7
        mulpd   nb313_tsc(%esp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        movd %eax,%mm0
        movd %ebx,%mm1

        pslld $2,%mm6           ## idx *= 4 
        movl nb313_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

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
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb313_two(%esp),%xmm7    ## two*Heps2 
        movapd nb313_qqM(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 
        ## increment vcoul 
        xorpd  %xmm4,%xmm4
        addpd  nb313_vctot(%esp),%xmm5
        mulpd  nb313_rinvM(%esp),%xmm3
        movapd %xmm5,nb313_vctot(%esp)
        mulpd  nb313_tsc(%esp),%xmm3
        subpd  %xmm3,%xmm4

        movapd nb313_dxM(%esp),%xmm0
        movapd nb313_dyM(%esp),%xmm1
        movapd nb313_dzM(%esp),%xmm2
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        movd %mm0,%eax
        movd %mm1,%ebx

        ## update H2 forces 
        movapd nb313_fixM(%esp),%xmm3
        movapd nb313_fiyM(%esp),%xmm4
        movapd nb313_fizM(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb313_fixM(%esp)
        movapd %xmm4,nb313_fiyM(%esp)
        movapd %xmm7,nb313_fizM(%esp)

        movl nb313_faction(%ebp),%edi
        ## update j forces 
        addpd  nb313_fjx(%esp),%xmm0
        addpd  nb313_fjy(%esp),%xmm1
        addpd  nb313_fjz(%esp),%xmm2

        ## the fj's - start by accumulating forces from memory 
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
        subl $2,nb313_innerk(%esp)
        jl    _nb_kernel313_ia32_sse2.nb313_checksingle
        jmp   _nb_kernel313_ia32_sse2.nb313_unroll_loop
_nb_kernel313_ia32_sse2.nb313_checksingle: 
        movl  nb313_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel313_ia32_sse2.nb313_dosingle
        jmp   _nb_kernel313_ia32_sse2.nb313_updateouterdata
_nb_kernel313_ia32_sse2.nb313_dosingle: 
        movl  nb313_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movl nb313_charge(%ebp),%esi     ## base of charge[] 
        xorpd %xmm3,%xmm3
        movlpd (%esi,%eax,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb313_iqM(%esp),%xmm3
        mulpd  nb313_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movapd  %xmm3,nb313_qqM(%esp)
        movapd  %xmm4,nb313_qqH(%esp)

        movl nb313_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb313_vdwparam(%ebp),%esi
        shll %eax
        movl nb313_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movapd %xmm4,nb313_c6(%esp)
        movapd %xmm6,nb313_c12(%esp)

        movl nb313_pos(%ebp),%esi        ## base of pos[] 
        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move coords to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb313_ixO(%esp),%xmm4
        movapd nb313_iyO(%esp),%xmm5
        movapd nb313_izO(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb313_dxO(%esp)
        movapd %xmm5,nb313_dyO(%esp)
        movapd %xmm6,nb313_dzO(%esp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb313_ixH1(%esp),%xmm4
        movapd nb313_iyH1(%esp),%xmm5
        movapd nb313_izH1(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb313_dxH1(%esp)
        movapd %xmm5,nb313_dyH1(%esp)
        movapd %xmm6,nb313_dzH1(%esp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm6
        addsd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb313_ixH2(%esp),%xmm3
        movapd nb313_iyH2(%esp),%xmm4
        movapd nb313_izH2(%esp),%xmm5

        ## calc dr 
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5

        ## store dr 
        movapd %xmm3,nb313_dxH2(%esp)
        movapd %xmm4,nb313_dyH2(%esp)
        movapd %xmm5,nb313_dzH2(%esp)
        ## square it 
        mulsd %xmm3,%xmm3
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        addsd %xmm4,%xmm5
        addsd %xmm3,%xmm5
        ## move ixM-izM to xmm2-xmm4  
        movapd nb313_iyM(%esp),%xmm3
        movapd nb313_izM(%esp),%xmm4
        subpd  %xmm1,%xmm3
        subpd  %xmm2,%xmm4
        movapd nb313_ixM(%esp),%xmm2
        subpd  %xmm0,%xmm2

        ## store dr 
        movapd %xmm2,nb313_dxM(%esp)
        movapd %xmm3,nb313_dyM(%esp)
        movapd %xmm4,nb313_dzM(%esp)
        ## square it 
        mulpd %xmm2,%xmm2
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        addpd %xmm3,%xmm4
        addpd %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## 1/x for O - rsqO is in xmm7
        cvtsd2ss %xmm7,%xmm2
        movsd   %xmm7,%xmm3
        rcpps    %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2
        movsd   nb313_two(%esp),%xmm1
        movsd   %xmm1,%xmm0
        mulsd   %xmm2,%xmm7
        subsd   %xmm7,%xmm1
        mulsd   %xmm1,%xmm2 ## iter1 
        mulsd   %xmm2,%xmm3
        subsd   %xmm3,%xmm0
        mulsd   %xmm2,%xmm0 ## xmm0=rinvsq
        movsd  %xmm0,nb313_rinvsqO(%esp)

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb313_three(%esp),%xmm0
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulsd   nb313_half(%esp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm6,%xmm2
        movapd %xmm0,%xmm3
        mulsd %xmm0,%xmm0       ## lu*lu 
        mulsd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb313_three(%esp),%xmm0
        subsd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulsd nb313_half(%esp),%xmm0   ## rinv 
        movapd %xmm0,nb313_rinvH1(%esp)         ## rinvH1 
        mulsd  %xmm0,%xmm6
        movapd %xmm6,nb313_rH1(%esp)    ## rH1 

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb313_three(%esp),%xmm0
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulsd   nb313_half(%esp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm5,%xmm2
        movapd %xmm0,%xmm3
        mulsd %xmm0,%xmm0       ## lu*lu 
        mulsd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb313_three(%esp),%xmm0
        subsd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulsd nb313_half(%esp),%xmm0   ## rinv 
        movapd %xmm0,nb313_rinvH2(%esp)   ## rinv 
        mulsd %xmm0,%xmm5
        movapd %xmm5,nb313_rH2(%esp)   ## r 

        ## rsqM - seed in xmm2 
        cvtsd2ss %xmm4,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb313_three(%esp),%xmm0
        mulsd   %xmm4,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulsd   nb313_half(%esp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm4,%xmm2
        movapd %xmm0,%xmm3
        mulsd %xmm0,%xmm0       ## lu*lu 
        mulsd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb313_three(%esp),%xmm0
        subsd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulsd nb313_half(%esp),%xmm0   ## rinv 
        movapd %xmm0,nb313_rinvM(%esp)   ## rinv 
        mulsd %xmm0,%xmm4
        movapd %xmm4,nb313_rM(%esp)   ## r 

        ## do O interactions 
        movapd  nb313_rinvsqO(%esp),%xmm0
        movapd  %xmm0,%xmm1
        mulsd   %xmm1,%xmm1 ## rinv4
        mulsd   %xmm0,%xmm1 ##rinvsix
        movapd  %xmm1,%xmm2
        mulsd   %xmm2,%xmm2 ## rinvtwelve
        mulsd  nb313_c6(%esp),%xmm1
        mulsd  nb313_c12(%esp),%xmm2
        movapd %xmm2,%xmm3
        subsd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addsd  nb313_Vvdwtot(%esp),%xmm3
        mulsd  nb313_six(%esp),%xmm1
        mulsd  nb313_twelve(%esp),%xmm2
        subsd  %xmm1,%xmm2
        mulsd  %xmm0,%xmm2
        movapd %xmm2,%xmm4 ## total fsO 
        movsd %xmm3,nb313_Vvdwtot(%esp)

        movapd nb313_dxO(%esp),%xmm0
        movapd nb313_dyO(%esp),%xmm1
        movapd nb313_dzO(%esp),%xmm2
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2
        ## tx in xmm0-xmm2 

        ## update O forces 
        movapd nb313_fixO(%esp),%xmm3
        movapd nb313_fiyO(%esp),%xmm4
        movapd nb313_fizO(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb313_fixO(%esp)
        movlpd %xmm4,nb313_fiyO(%esp)
        movlpd %xmm7,nb313_fizO(%esp)
        ## update j forces with water O 
        movlpd %xmm0,nb313_fjx(%esp)
        movlpd %xmm1,nb313_fjy(%esp)
        movlpd %xmm2,nb313_fjz(%esp)

        movd %eax,%mm0

        ## Done with O interactions - now H1! 
        movapd nb313_rH1(%esp),%xmm7
        mulpd nb313_tsc(%esp),%xmm7
        cvttsd2si %xmm7,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb313_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb313_two(%esp),%xmm7    ## two*Heps2 
        movapd nb313_qqH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 
    ## increment vcoul 
        xorpd  %xmm4,%xmm4
    addsd  nb313_vctot(%esp),%xmm5
        mulsd  nb313_rinvH1(%esp),%xmm3
    movlpd %xmm5,nb313_vctot(%esp)
        mulsd  nb313_tsc(%esp),%xmm3
        subsd %xmm3,%xmm4

        movapd nb313_dxH1(%esp),%xmm0
        movapd nb313_dyH1(%esp),%xmm1
        movapd nb313_dzH1(%esp),%xmm2
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb313_fixH1(%esp),%xmm3
        movapd nb313_fiyH1(%esp),%xmm4
        movapd nb313_fizH1(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb313_fixH1(%esp)
        movlpd %xmm4,nb313_fiyH1(%esp)
        movlpd %xmm7,nb313_fizH1(%esp)
        ## update j forces with water H1 
        addsd  nb313_fjx(%esp),%xmm0
        addsd  nb313_fjy(%esp),%xmm1
        addsd  nb313_fjz(%esp),%xmm2
        movlpd %xmm0,nb313_fjx(%esp)
        movlpd %xmm1,nb313_fjy(%esp)
        movlpd %xmm2,nb313_fjz(%esp)

        ##  H2 interactions 
        movapd nb313_rH2(%esp),%xmm7
        mulsd   nb313_tsc(%esp),%xmm7
        cvttsd2si %xmm7,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb313_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb313_two(%esp),%xmm7    ## two*Heps2 
        movapd nb313_qqH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 
        ## increment vcoul 
        xorpd  %xmm4,%xmm4
        addsd  nb313_vctot(%esp),%xmm5
        mulsd  nb313_rinvH2(%esp),%xmm3
        movlpd %xmm5,nb313_vctot(%esp)
        mulsd  nb313_tsc(%esp),%xmm3
        subsd  %xmm3,%xmm4

        movapd nb313_dxH2(%esp),%xmm0
        movapd nb313_dyH2(%esp),%xmm1
        movapd nb313_dzH2(%esp),%xmm2
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        movd %mm0,%eax

        ## update H2 forces 
        movapd nb313_fixH2(%esp),%xmm3
        movapd nb313_fiyH2(%esp),%xmm4
        movapd nb313_fizH2(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb313_fixH2(%esp)
        movlpd %xmm4,nb313_fiyH2(%esp)
        movlpd %xmm7,nb313_fizH2(%esp)
        ## update j forces with water H1 
        addsd  nb313_fjx(%esp),%xmm0
        addsd  nb313_fjy(%esp),%xmm1
        addsd  nb313_fjz(%esp),%xmm2
        movlpd %xmm0,nb313_fjx(%esp)
        movlpd %xmm1,nb313_fjy(%esp)
        movlpd %xmm2,nb313_fjz(%esp)

        ## M interactions 
        movapd nb313_rM(%esp),%xmm7
        mulsd   nb313_tsc(%esp),%xmm7
        cvttsd2si %xmm7,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb313_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb313_two(%esp),%xmm7    ## two*Heps2 
        movapd nb313_qqM(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 
        ## increment vcoul 
        xorpd  %xmm4,%xmm4
        addsd  nb313_vctot(%esp),%xmm5
        mulsd  nb313_rinvM(%esp),%xmm3
        movlpd %xmm5,nb313_vctot(%esp)
        mulsd  nb313_tsc(%esp),%xmm3
        subsd  %xmm3,%xmm4

        movapd nb313_dxM(%esp),%xmm0
        movapd nb313_dyM(%esp),%xmm1
        movapd nb313_dzM(%esp),%xmm2
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        movd %mm0,%eax

        ## update M forces 
        movapd nb313_fixM(%esp),%xmm3
        movapd nb313_fiyM(%esp),%xmm4
        movapd nb313_fizM(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb313_fixM(%esp)
        movlpd %xmm4,nb313_fiyM(%esp)
        movlpd %xmm7,nb313_fizM(%esp)

        movl nb313_faction(%ebp),%edi
        ## update j forces 
        addsd  nb313_fjx(%esp),%xmm0
        addsd  nb313_fjy(%esp),%xmm1
        addsd  nb313_fjz(%esp),%xmm2

        ## the fj's - start by accumulating forces from memory 
        movlpd (%edi,%eax,8),%xmm3
        movlpd 8(%edi,%eax,8),%xmm4
        movlpd 16(%edi,%eax,8),%xmm5
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        movlpd %xmm3,(%edi,%eax,8)
        movlpd %xmm4,8(%edi,%eax,8)
        movlpd %xmm5,16(%edi,%eax,8)

_nb_kernel313_ia32_sse2.nb313_updateouterdata: 
        movl  nb313_ii3(%esp),%ecx
        movl  nb313_faction(%ebp),%edi
        movl  nb313_fshift(%ebp),%esi
        movl  nb313_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movapd nb313_fixO(%esp),%xmm0
        movapd nb313_fiyO(%esp),%xmm1
        movapd nb313_fizO(%esp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

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
        movapd nb313_fixH1(%esp),%xmm0
        movapd nb313_fiyH1(%esp),%xmm1
        movapd nb313_fizH1(%esp),%xmm2

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
        movapd nb313_fixH2(%esp),%xmm0
        movapd nb313_fiyH2(%esp),%xmm1
        movapd nb313_fizH2(%esp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

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
        movapd nb313_fixM(%esp),%xmm0
        movapd nb313_fiyM(%esp),%xmm1
        movapd nb313_fizM(%esp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

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
        movl nb313_n(%esp),%esi
        ## get group index for i particle 
        movl  nb313_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb313_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb313_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb313_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb313_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb313_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel313_ia32_sse2.nb313_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb313_n(%esp)
        jmp _nb_kernel313_ia32_sse2.nb313_outer
_nb_kernel313_ia32_sse2.nb313_outerend: 
        ## check if more outer neighborlists remain
        movl  nb313_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel313_ia32_sse2.nb313_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel313_ia32_sse2.nb313_threadloop
_nb_kernel313_ia32_sse2.nb313_end: 
        emms

        movl nb313_nouter(%esp),%eax
        movl nb313_ninner(%esp),%ebx
        movl nb313_outeriter(%ebp),%ecx
        movl nb313_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb313_salign(%esp),%eax
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




.globl nb_kernel313nf_ia32_sse2
.globl _nb_kernel313nf_ia32_sse2
nb_kernel313nf_ia32_sse2:       
_nb_kernel313nf_ia32_sse2:      
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
        ## bottom of stack is cache-aligned for sse2 use 
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
.set nb313nf_rinvsqO, 256
.set nb313nf_rinvH1, 272
.set nb313nf_rinvH2, 288
.set nb313nf_rinvM, 304
.set nb313nf_rO, 320
.set nb313nf_rH1, 336
.set nb313nf_rH2, 352
.set nb313nf_rM, 368
.set nb313nf_tsc, 384
.set nb313nf_two, 400
.set nb313nf_c6, 416
.set nb313nf_c12, 432
.set nb313nf_vctot, 448
.set nb313nf_Vvdwtot, 464
.set nb313nf_half, 480
.set nb313nf_three, 496
.set nb313nf_is3, 512
.set nb313nf_ii3, 516
.set nb313nf_ntia, 520
.set nb313nf_innerjjnr, 524
.set nb313nf_innerk, 528
.set nb313nf_n, 532
.set nb313nf_nn1, 536
.set nb313nf_nri, 540
.set nb313nf_nouter, 544
.set nb313nf_ninner, 548
.set nb313nf_salign, 552
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $556,%esp          ## local stack space 
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


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb313nf_half(%esp)
        movl %ebx,nb313nf_half+4(%esp)
        movsd nb313nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb313nf_half(%esp)
        movapd %xmm2,nb313nf_two(%esp)
        movapd %xmm3,nb313nf_three(%esp)
        movl nb313nf_p_tabscale(%ebp),%eax
        movsd (%eax),%xmm5
        shufpd $0,%xmm5,%xmm5
        movapd %xmm5,nb313nf_tsc(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb313nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb313nf_charge(%ebp),%edx
        movsd 8(%edx,%ebx,8),%xmm3
        movsd 24(%edx,%ebx,8),%xmm4
        movl nb313nf_p_facel(%ebp),%esi
        movsd (%esi),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb313nf_iqH(%esp)
        movapd %xmm4,nb313nf_iqM(%esp)

        movl  nb313nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb313nf_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb313nf_ntia(%esp)
_nb_kernel313nf_ia32_sse2.nb313nf_threadloop: 
        movl  nb313nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel313nf_ia32_sse2.nb313nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel313nf_ia32_sse2.nb313nf_spinlock

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
        jg  _nb_kernel313nf_ia32_sse2.nb313nf_outerstart
        jmp _nb_kernel313nf_ia32_sse2.nb313nf_end

_nb_kernel313nf_ia32_sse2.nb313nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb313nf_nouter(%esp),%ebx
        movl %ebx,nb313nf_nouter(%esp)

_nb_kernel313nf_ia32_sse2.nb313nf_outer: 
        movl  nb313nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb313nf_is3(%esp)            ## store is3 

        movl  nb313nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb313nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        movapd %xmm0,%xmm6
        movapd %xmm1,%xmm7

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb313nf_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb313nf_ii3(%esp)

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
        movapd %xmm3,nb313nf_ixO(%esp)
        movapd %xmm4,nb313nf_iyO(%esp)
        movapd %xmm5,nb313nf_izO(%esp)
        movapd %xmm6,nb313nf_ixH1(%esp)
        movapd %xmm7,nb313nf_iyH1(%esp)

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
        movapd %xmm6,nb313nf_izH1(%esp)
        movapd %xmm0,nb313nf_ixH2(%esp)
        movapd %xmm1,nb313nf_iyH2(%esp)
        movapd %xmm2,nb313nf_izH2(%esp)
        movapd %xmm3,nb313nf_ixM(%esp)
        movapd %xmm4,nb313nf_iyM(%esp)
        movapd %xmm5,nb313nf_izM(%esp)

        ## clear vctot
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb313nf_vctot(%esp)
        movapd %xmm4,nb313nf_Vvdwtot(%esp)

        movl  nb313nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx     ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb313nf_pos(%ebp),%esi
        movl  nb313nf_faction(%ebp),%edi
        movl  nb313nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb313nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb313nf_ninner(%esp),%ecx
        movl  %ecx,nb313nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb313nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel313nf_ia32_sse2.nb313nf_unroll_loop
        jmp   _nb_kernel313nf_ia32_sse2.nb313nf_checksingle
_nb_kernel313nf_ia32_sse2.nb313nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb313nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb313nf_innerjjnr(%esp)             ## advance pointer (unrolled 2) 

        movl nb313nf_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb313nf_iqM(%esp),%xmm3
        mulpd  nb313nf_iqH(%esp),%xmm4
        movapd  %xmm3,nb313nf_qqM(%esp)
        movapd  %xmm4,nb313nf_qqH(%esp)

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movl nb313nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb313nf_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb313nf_ntia(%esp),%edi
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
        movapd %xmm4,nb313nf_c6(%esp)
        movapd %xmm6,nb313nf_c12(%esp)

        movl nb313nf_pos(%ebp),%esi        ## base of pos[] 

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
        movapd nb313nf_ixO(%esp),%xmm4
        movapd nb313nf_iyO(%esp),%xmm5
        movapd nb313nf_izO(%esp),%xmm6

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
        movapd nb313nf_ixH1(%esp),%xmm4
        movapd nb313nf_iyH1(%esp),%xmm5
        movapd nb313nf_izH1(%esp),%xmm6

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
        movapd nb313nf_ixH2(%esp),%xmm3
        movapd nb313nf_iyH2(%esp),%xmm4
        movapd nb313nf_izH2(%esp),%xmm5

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
        movapd nb313nf_iyM(%esp),%xmm3
        movapd nb313nf_izM(%esp),%xmm4
        subpd  %xmm1,%xmm3
        subpd  %xmm2,%xmm4
        movapd nb313nf_ixM(%esp),%xmm2
        subpd  %xmm0,%xmm2

        ## square it 
        mulpd %xmm2,%xmm2
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        addpd %xmm3,%xmm4
        addpd %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## 1/x for O - rsqO is in xmm7
        cvtpd2ps %xmm7,%xmm2
        movapd   %xmm7,%xmm3
        rcpps    %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2
        movapd   nb313nf_two(%esp),%xmm1
        movapd   %xmm1,%xmm0
        mulpd   %xmm2,%xmm7
        subpd   %xmm7,%xmm1
        mulpd   %xmm1,%xmm2 ## iter1 
        mulpd   %xmm2,%xmm3
        subpd   %xmm3,%xmm0
        mulpd   %xmm2,%xmm0 ## xmm0=rinvsq
        movapd  %xmm0,nb313nf_rinvsqO(%esp)

        ## rsqH1 - seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb313nf_three(%esp),%xmm0
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulpd   nb313nf_half(%esp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm6,%xmm2
        movapd %xmm0,%xmm3
        mulpd %xmm0,%xmm0       ## lu*lu 
        mulpd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb313nf_three(%esp),%xmm0
        subpd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulpd nb313nf_half(%esp),%xmm0   ## rinv 
        movapd %xmm0,nb313nf_rinvH1(%esp)       ## rinvH1 
        mulpd  %xmm0,%xmm6
        movapd %xmm6,nb313nf_rH1(%esp)          ## rH1 

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb313nf_three(%esp),%xmm0
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulpd   nb313nf_half(%esp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm5,%xmm2
        movapd %xmm0,%xmm3
        mulpd %xmm0,%xmm0       ## lu*lu 
        mulpd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb313nf_three(%esp),%xmm0
        subpd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulpd nb313nf_half(%esp),%xmm0   ## rinv 
        movapd %xmm0,nb313nf_rinvH2(%esp)   ## rinv 
        mulpd %xmm0,%xmm5
        movapd %xmm5,nb313nf_rH2(%esp)   ## r 

        ## rsqM - seed in xmm2 
        cvtpd2ps %xmm4,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb313nf_three(%esp),%xmm0
        mulpd   %xmm4,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulpd   nb313nf_half(%esp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm4,%xmm2
        movapd %xmm0,%xmm3
        mulpd %xmm0,%xmm0       ## lu*lu 
        mulpd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb313nf_three(%esp),%xmm0
        subpd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulpd nb313nf_half(%esp),%xmm0   ## rinv 
        movapd %xmm0,nb313nf_rinvM(%esp)   ## rinv 
        mulpd %xmm0,%xmm4
        movapd %xmm4,nb313nf_rM(%esp)   ## r 

        ## do O interactions 
        movapd nb313nf_rinvsqO(%esp),%xmm0
        movapd  %xmm0,%xmm1
        mulpd   %xmm1,%xmm1 ## rinv4
        mulpd   %xmm0,%xmm1 ##rinvsix
        movapd  %xmm1,%xmm2
        mulpd   %xmm2,%xmm2 ## rinvtwelve
        mulpd  nb313nf_c6(%esp),%xmm1
        mulpd  nb313nf_c12(%esp),%xmm2
        movapd %xmm2,%xmm3
        subpd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addpd  nb313nf_Vvdwtot(%esp),%xmm3
        movapd %xmm3,nb313nf_Vvdwtot(%esp)

        ## Done with O interactions - now H1! 
        movapd nb313nf_rH1(%esp),%xmm7
        mulpd nb313nf_tsc(%esp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        movd %eax,%mm0
        movd %ebx,%mm1

        pslld $2,%mm6           ## idx *= 4 
        movl nb313nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

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
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb313nf_qqH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## increment vcoul 
        addpd  nb313nf_vctot(%esp),%xmm5
        movapd %xmm5,nb313nf_vctot(%esp)

        ## H2 interactions 
        movapd nb313nf_rH2(%esp),%xmm7
        mulpd   nb313nf_tsc(%esp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb313nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

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
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb313nf_qqH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## increment vcoul 
        addpd  nb313nf_vctot(%esp),%xmm5
        movapd %xmm5,nb313nf_vctot(%esp)

        ## M interactions 
        movapd nb313nf_rM(%esp),%xmm7
        mulpd   nb313nf_tsc(%esp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        movd %eax,%mm0
        movd %ebx,%mm1

        pslld $2,%mm6           ## idx *= 4 
        movl nb313nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

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
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb313nf_qqM(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## increment vcoul 
        addpd  nb313nf_vctot(%esp),%xmm5
        movapd %xmm5,nb313nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $2,nb313nf_innerk(%esp)
        jl    _nb_kernel313nf_ia32_sse2.nb313nf_checksingle
        jmp   _nb_kernel313nf_ia32_sse2.nb313nf_unroll_loop
_nb_kernel313nf_ia32_sse2.nb313nf_checksingle: 
        movl  nb313nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel313nf_ia32_sse2.nb313nf_dosingle
        jmp   _nb_kernel313nf_ia32_sse2.nb313nf_updateouterdata
_nb_kernel313nf_ia32_sse2.nb313nf_dosingle: 
        movl  nb313nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movl nb313nf_charge(%ebp),%esi     ## base of charge[] 
        xorpd %xmm3,%xmm3
        movlpd (%esi,%eax,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb313nf_iqM(%esp),%xmm3
        mulpd  nb313nf_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movapd  %xmm3,nb313nf_qqM(%esp)
        movapd  %xmm4,nb313nf_qqH(%esp)

        movl nb313nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb313nf_vdwparam(%ebp),%esi
        shll %eax
        movl nb313nf_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movapd %xmm4,nb313nf_c6(%esp)
        movapd %xmm6,nb313nf_c12(%esp)

        movl nb313nf_pos(%ebp),%esi        ## base of pos[] 
        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move coords to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb313nf_ixO(%esp),%xmm4
        movapd nb313nf_iyO(%esp),%xmm5
        movapd nb313nf_izO(%esp),%xmm6

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

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb313nf_ixH1(%esp),%xmm4
        movapd nb313nf_iyH1(%esp),%xmm5
        movapd nb313nf_izH1(%esp),%xmm6

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
        movapd nb313nf_ixH2(%esp),%xmm3
        movapd nb313nf_iyH2(%esp),%xmm4
        movapd nb313nf_izH2(%esp),%xmm5

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
        movapd nb313nf_iyM(%esp),%xmm3
        movapd nb313nf_izM(%esp),%xmm4
        subpd  %xmm1,%xmm3
        subpd  %xmm2,%xmm4
        movapd nb313nf_ixM(%esp),%xmm2
        subpd  %xmm0,%xmm2

        ## square it 
        mulpd %xmm2,%xmm2
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        addpd %xmm3,%xmm4
        addpd %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## 1/x for O - rsqO is in xmm7
        cvtsd2ss %xmm7,%xmm2
        movsd   %xmm7,%xmm3
        rcpps    %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2
        movsd   nb313nf_two(%esp),%xmm1
        movsd   %xmm1,%xmm0
        mulsd   %xmm2,%xmm7
        subsd   %xmm7,%xmm1
        mulsd   %xmm1,%xmm2 ## iter1 
        mulsd   %xmm2,%xmm3
        subsd   %xmm3,%xmm0
        mulsd   %xmm2,%xmm0 ## xmm0=rinvsq
        movsd  %xmm0,nb313nf_rinvsqO(%esp)

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb313nf_three(%esp),%xmm0
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulsd   nb313nf_half(%esp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm6,%xmm2
        movapd %xmm0,%xmm3
        mulsd %xmm0,%xmm0       ## lu*lu 
        mulsd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb313nf_three(%esp),%xmm0
        subsd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulsd nb313nf_half(%esp),%xmm0   ## rinv 
        movapd %xmm0,nb313nf_rinvH1(%esp)       ## rinvH1 
        mulsd  %xmm0,%xmm6
        movapd %xmm6,nb313nf_rH1(%esp)          ## rH1 

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb313nf_three(%esp),%xmm0
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulsd   nb313nf_half(%esp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm5,%xmm2
        movapd %xmm0,%xmm3
        mulsd %xmm0,%xmm0       ## lu*lu 
        mulsd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb313nf_three(%esp),%xmm0
        subsd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulsd nb313nf_half(%esp),%xmm0   ## rinv 
        movapd %xmm0,nb313nf_rinvH2(%esp)   ## rinv 
        mulsd %xmm0,%xmm5
        movapd %xmm5,nb313nf_rH2(%esp)   ## r 

        ## rsqM - seed in xmm2 
        cvtsd2ss %xmm4,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb313nf_three(%esp),%xmm0
        mulsd   %xmm4,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulsd   nb313nf_half(%esp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm4,%xmm2
        movapd %xmm0,%xmm3
        mulsd %xmm0,%xmm0       ## lu*lu 
        mulsd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb313nf_three(%esp),%xmm0
        subsd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulsd nb313nf_half(%esp),%xmm0   ## rinv 
        movapd %xmm0,nb313nf_rinvM(%esp)   ## rinv 
        mulsd %xmm0,%xmm4
        movapd %xmm4,nb313nf_rM(%esp)   ## r 

        ## do O interactions 
        movapd  nb313nf_rinvsqO(%esp),%xmm0
        movapd  %xmm0,%xmm1
        mulsd   %xmm1,%xmm1 ## rinv4
        mulsd   %xmm0,%xmm1 ##rinvsix
        movapd  %xmm1,%xmm2
        mulsd   %xmm2,%xmm2 ## rinvtwelve
        mulsd  nb313nf_c6(%esp),%xmm1
        mulsd  nb313nf_c12(%esp),%xmm2
        movapd %xmm2,%xmm3
        subsd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addsd  nb313nf_Vvdwtot(%esp),%xmm3
        movsd %xmm3,nb313nf_Vvdwtot(%esp)

        ## Done with O interactions - now H1! 
        movapd nb313nf_rH1(%esp),%xmm7
        mulpd nb313nf_tsc(%esp),%xmm7
        cvttsd2si %xmm7,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb313nf_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb313nf_qqH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## increment vcoul 
        addsd  nb313nf_vctot(%esp),%xmm5
        movlpd %xmm5,nb313nf_vctot(%esp)

        ##  H2 interactions 
        movapd nb313nf_rH2(%esp),%xmm7
        mulsd   nb313nf_tsc(%esp),%xmm7
        cvttsd2si %xmm7,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb313nf_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb313nf_qqH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## increment vcoul 
        addsd  nb313nf_vctot(%esp),%xmm5
        movlpd %xmm5,nb313nf_vctot(%esp)

        ## M interactions 
        movapd nb313nf_rM(%esp),%xmm7
        mulsd   nb313nf_tsc(%esp),%xmm7
        cvttsd2si %xmm7,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb313nf_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   

        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb313nf_qqM(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## increment vcoul 
        addsd  nb313nf_vctot(%esp),%xmm5
        movlpd %xmm5,nb313nf_vctot(%esp)

_nb_kernel313nf_ia32_sse2.nb313nf_updateouterdata: 
        ## get n from stack
        movl nb313nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb313nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb313nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb313nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb313nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb313nf_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb313nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel313nf_ia32_sse2.nb313nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb313nf_n(%esp)
        jmp _nb_kernel313nf_ia32_sse2.nb313nf_outer
_nb_kernel313nf_ia32_sse2.nb313nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb313nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel313nf_ia32_sse2.nb313nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel313nf_ia32_sse2.nb313nf_threadloop
_nb_kernel313nf_ia32_sse2.nb313nf_end: 
        emms

        movl nb313nf_nouter(%esp),%eax
        movl nb313nf_ninner(%esp),%ebx
        movl nb313nf_outeriter(%ebp),%ecx
        movl nb313nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb313nf_salign(%esp),%eax
        addl %eax,%esp
        addl $556,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret


