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



.globl nb_kernel213_ia32_sse2
.globl _nb_kernel213_ia32_sse2
nb_kernel213_ia32_sse2: 
_nb_kernel213_ia32_sse2:        
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
        ## bottom of stack is cache-aligned for sse2 use 
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
.set nb213_iqH, 192
.set nb213_iqM, 208
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
.set nb213_qqH, 416
.set nb213_qqM, 432
.set nb213_c6, 448
.set nb213_c12, 464
.set nb213_six, 480
.set nb213_twelve, 496
.set nb213_vctot, 512
.set nb213_Vvdwtot, 528
.set nb213_fixO, 544
.set nb213_fiyO, 560
.set nb213_fizO, 576
.set nb213_fixH1, 592
.set nb213_fiyH1, 608
.set nb213_fizH1, 624
.set nb213_fixH2, 640
.set nb213_fiyH2, 656
.set nb213_fizH2, 672
.set nb213_fixM, 688
.set nb213_fiyM, 704
.set nb213_fizM, 720
.set nb213_fjx, 736
.set nb213_fjy, 752
.set nb213_fjz, 768
.set nb213_half, 784
.set nb213_three, 800
.set nb213_two, 816
.set nb213_rinvH1, 832
.set nb213_rinvH2, 848
.set nb213_rinvM, 864
.set nb213_krsqH1, 880
.set nb213_krsqH2, 896
.set nb213_krsqM, 912
.set nb213_krf, 928
.set nb213_crf, 944
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
        movsd (%esi),%xmm5
        movsd (%edi),%xmm6
        shufpd $0,%xmm5,%xmm5
        shufpd $0,%xmm6,%xmm6
        movapd %xmm5,nb213_krf(%esp)
        movapd %xmm6,nb213_crf(%esp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb213_half(%esp)
        movl %ebx,nb213_half+4(%esp)
        movsd nb213_half(%esp),%xmm1
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
        movapd %xmm1,nb213_half(%esp)
        movapd %xmm2,nb213_two(%esp)
        movapd %xmm3,nb213_three(%esp)
        movapd %xmm4,nb213_six(%esp)
        movapd %xmm5,nb213_twelve(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb213_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb213_charge(%ebp),%edx
        movsd 8(%edx,%ebx,8),%xmm3
        movsd 24(%edx,%ebx,8),%xmm4
        movl nb213_p_facel(%ebp),%esi
        movsd (%esi),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb213_iqH(%esp)
        movapd %xmm4,nb213_iqM(%esp)

        movl  nb213_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb213_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb213_ntia(%esp)
_nb_kernel213_ia32_sse2.nb213_threadloop: 
        movl  nb213_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel213_ia32_sse2.nb213_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel213_ia32_sse2.nb213_spinlock

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
        jg  _nb_kernel213_ia32_sse2.nb213_outerstart
        jmp _nb_kernel213_ia32_sse2.nb213_end

_nb_kernel213_ia32_sse2.nb213_outerstart: 
        ## ebx contains number of outer iterations
        addl nb213_nouter(%esp),%ebx
        movl %ebx,nb213_nouter(%esp)

_nb_kernel213_ia32_sse2.nb213_outer: 
        movl  nb213_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb213_is3(%esp)      ## store is3 

        movl  nb213_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb213_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        movapd %xmm0,%xmm6
        movapd %xmm1,%xmm7

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb213_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb213_ii3(%esp)

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
        movapd %xmm3,nb213_ixO(%esp)
        movapd %xmm4,nb213_iyO(%esp)
        movapd %xmm5,nb213_izO(%esp)
        movapd %xmm6,nb213_ixH1(%esp)
        movapd %xmm7,nb213_iyH1(%esp)

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
        movapd %xmm6,nb213_izH1(%esp)
        movapd %xmm0,nb213_ixH2(%esp)
        movapd %xmm1,nb213_iyH2(%esp)
        movapd %xmm2,nb213_izH2(%esp)
        movapd %xmm3,nb213_ixM(%esp)
        movapd %xmm4,nb213_iyM(%esp)
        movapd %xmm5,nb213_izM(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb213_vctot(%esp)
        movapd %xmm4,nb213_Vvdwtot(%esp)
        movapd %xmm4,nb213_fixO(%esp)
        movapd %xmm4,nb213_fiyO(%esp)
        movapd %xmm4,nb213_fizO(%esp)
        movapd %xmm4,nb213_fixH1(%esp)
        movapd %xmm4,nb213_fiyH1(%esp)
        movapd %xmm4,nb213_fizH1(%esp)
        movapd %xmm4,nb213_fixH2(%esp)
        movapd %xmm4,nb213_fiyH2(%esp)
        movapd %xmm4,nb213_fizH2(%esp)
        movapd %xmm4,nb213_fixM(%esp)
        movapd %xmm4,nb213_fiyM(%esp)
        movapd %xmm4,nb213_fizM(%esp)

        movl  nb213_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb213_pos(%ebp),%esi
        movl  nb213_faction(%ebp),%edi
        movl  nb213_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb213_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb213_ninner(%esp),%ecx
        movl  %ecx,nb213_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb213_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel213_ia32_sse2.nb213_unroll_loop
        jmp   _nb_kernel213_ia32_sse2.nb213_checksingle
_nb_kernel213_ia32_sse2.nb213_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb213_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb213_innerjjnr(%esp)                   ## advance pointer (unrolled 2) 

        movl nb213_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb213_iqM(%esp),%xmm3
        mulpd  nb213_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movapd  %xmm3,nb213_qqM(%esp)
        movapd  %xmm4,nb213_qqH(%esp)

        movl nb213_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb213_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb213_ntia(%esp),%edi
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
        movapd %xmm4,nb213_c6(%esp)
        movapd %xmm6,nb213_c12(%esp)

        movl nb213_pos(%ebp),%esi        ## base of pos[] 

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
        movapd nb213_ixO(%esp),%xmm4
        movapd nb213_iyO(%esp),%xmm5
        movapd nb213_izO(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb213_dxO(%esp)
        movapd %xmm5,nb213_dyO(%esp)
        movapd %xmm6,nb213_dzO(%esp)
        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb213_ixH1(%esp),%xmm4
        movapd nb213_iyH1(%esp),%xmm5
        movapd nb213_izH1(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb213_dxH1(%esp)
        movapd %xmm5,nb213_dyH1(%esp)
        movapd %xmm6,nb213_dzH1(%esp)
        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm6
        addpd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb213_ixH2(%esp),%xmm3
        movapd nb213_iyH2(%esp),%xmm4
        movapd nb213_izH2(%esp),%xmm5

        ## calc dr 
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5

        ## store dr 
        movapd %xmm3,nb213_dxH2(%esp)
        movapd %xmm4,nb213_dyH2(%esp)
        movapd %xmm5,nb213_dzH2(%esp)
        ## square it 
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        addpd %xmm4,%xmm5
        addpd %xmm3,%xmm5

        ## move ixM-izM to xmm2-xmm4  
        movapd nb213_iyM(%esp),%xmm3
        movapd nb213_izM(%esp),%xmm4
        subpd  %xmm1,%xmm3
        subpd  %xmm2,%xmm4
        movapd nb213_ixM(%esp),%xmm2
        subpd  %xmm0,%xmm2

        ## store dr 
        movapd %xmm2,nb213_dxM(%esp)
        movapd %xmm3,nb213_dyM(%esp)
        movapd %xmm4,nb213_dzM(%esp)
        ## square it 
        mulpd %xmm2,%xmm2
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        addpd %xmm3,%xmm4
        addpd %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## calculate krsq
        movapd nb213_krf(%esp),%xmm0
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2
        mulpd %xmm4,%xmm0
        mulpd %xmm5,%xmm1
        mulpd %xmm6,%xmm2
        movapd %xmm0,nb213_krsqM(%esp)
        movapd %xmm1,nb213_krsqH2(%esp)
        movapd %xmm2,nb213_krsqH1(%esp)

        ## start with rsqH1 - put seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb213_three(%esp),%xmm1
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulpd   nb213_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulpd %xmm1,%xmm1       ## lu*lu 
        mulpd %xmm1,%xmm6       ## rsq*lu*lu 
        movapd nb213_three(%esp),%xmm1
        subpd %xmm6,%xmm1       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulpd nb213_half(%esp),%xmm1   ## rinv 
        movapd  %xmm1,nb213_rinvH1(%esp)

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb213_three(%esp),%xmm1
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulpd   nb213_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulpd %xmm1,%xmm1       ## lu*lu 
        mulpd %xmm1,%xmm5       ## rsq*lu*lu 
        movapd nb213_three(%esp),%xmm1
        subpd %xmm5,%xmm1       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulpd nb213_half(%esp),%xmm1   ## rinv 
        movapd  %xmm1,nb213_rinvH2(%esp)

        ## rsqM - seed in xmm2 
        cvtpd2ps %xmm4,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb213_three(%esp),%xmm1
        mulpd   %xmm4,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulpd   nb213_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulpd %xmm1,%xmm1       ## lu*lu 
        mulpd %xmm1,%xmm4       ## rsq*lu*lu 
        movapd nb213_three(%esp),%xmm1
        subpd %xmm4,%xmm1       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulpd nb213_half(%esp),%xmm1   ## rinv 
        movapd  %xmm1,nb213_rinvM(%esp)

        ## do O interactions directly - rsqO is in xmm7
        cvtpd2ps %xmm7,%xmm2
        movapd   %xmm7,%xmm6
        rcpps    %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2
        movapd   nb213_two(%esp),%xmm1
        movapd   %xmm1,%xmm0
        mulpd   %xmm2,%xmm7
        subpd   %xmm7,%xmm1
        mulpd   %xmm1,%xmm2 ## iter1 
        mulpd   %xmm2,%xmm6
        subpd   %xmm6,%xmm0
        mulpd   %xmm2,%xmm0 ## xmm0=rinvsq
        movapd  %xmm0,%xmm1
        mulpd   %xmm1,%xmm1 ## rinv4
        mulpd   %xmm0,%xmm1 ##rinvsix
        movapd  %xmm1,%xmm2
        mulpd   %xmm2,%xmm2 ## rinvtwelve
        mulpd  nb213_c6(%esp),%xmm1
        mulpd  nb213_c12(%esp),%xmm2
        movapd %xmm2,%xmm3
        subpd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addpd  nb213_Vvdwtot(%esp),%xmm3
        mulpd  nb213_six(%esp),%xmm1
        mulpd  nb213_twelve(%esp),%xmm2
        subpd  %xmm1,%xmm2
        mulpd  %xmm0,%xmm2
        movapd %xmm2,%xmm4 ## total fsO 
        movapd %xmm3,nb213_Vvdwtot(%esp)

        movapd nb213_dxO(%esp),%xmm0
        movapd nb213_dyO(%esp),%xmm1
        movapd nb213_dzO(%esp),%xmm2
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        ## update O forces 
        movapd nb213_fixO(%esp),%xmm3
        movapd nb213_fiyO(%esp),%xmm4
        movapd nb213_fizO(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb213_fixO(%esp)
        movapd %xmm4,nb213_fiyO(%esp)
        movapd %xmm7,nb213_fizO(%esp)
        ## update j forces with water O 
        movapd %xmm0,nb213_fjx(%esp)
        movapd %xmm1,nb213_fjy(%esp)
        movapd %xmm2,nb213_fjz(%esp)

        ## H1 interactions 
        movapd  nb213_rinvH1(%esp),%xmm6
        movapd  %xmm6,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        movapd  %xmm6,%xmm7
        movapd  nb213_krsqH1(%esp),%xmm0
        addpd   %xmm0,%xmm6     ## xmm6=rinv+ krsq 
        mulpd   nb213_two(%esp),%xmm0
        subpd   nb213_crf(%esp),%xmm6
        subpd   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulpd   nb213_qqH(%esp),%xmm6   ## vcoul 
        mulpd   nb213_qqH(%esp),%xmm7
        mulpd  %xmm7,%xmm4              ## total fsH1 in xmm4 
        addpd  nb213_vctot(%esp),%xmm6

        movapd nb213_dxH1(%esp),%xmm0
        movapd nb213_dyH1(%esp),%xmm1
        movapd nb213_dzH1(%esp),%xmm2
        movapd %xmm6,nb213_vctot(%esp)
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb213_fixH1(%esp),%xmm3
        movapd nb213_fiyH1(%esp),%xmm4
        movapd nb213_fizH1(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb213_fixH1(%esp)
        movapd %xmm4,nb213_fiyH1(%esp)
        movapd %xmm7,nb213_fizH1(%esp)
        ## update j forces with water H1 
        addpd  nb213_fjx(%esp),%xmm0
        addpd  nb213_fjy(%esp),%xmm1
        addpd  nb213_fjz(%esp),%xmm2
        movapd %xmm0,nb213_fjx(%esp)
        movapd %xmm1,nb213_fjy(%esp)
        movapd %xmm2,nb213_fjz(%esp)

        ## H2 interactions 
        movapd  nb213_rinvH2(%esp),%xmm5
        movapd  %xmm5,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        movapd  %xmm5,%xmm7
        movapd  nb213_krsqH2(%esp),%xmm0
        addpd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        mulpd   nb213_two(%esp),%xmm0
        subpd   nb213_crf(%esp),%xmm5
        subpd   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulpd   nb213_qqH(%esp),%xmm5   ## vcoul 
        mulpd   nb213_qqH(%esp),%xmm7
        mulpd  %xmm7,%xmm4              ## total fsH2 in xmm4 

        addpd  nb213_vctot(%esp),%xmm5

        movapd nb213_dxH2(%esp),%xmm0
        movapd nb213_dyH2(%esp),%xmm1
        movapd nb213_dzH2(%esp),%xmm2
        movapd %xmm5,nb213_vctot(%esp)
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        ## update H2 forces 
        movapd nb213_fixH2(%esp),%xmm3
        movapd nb213_fiyH2(%esp),%xmm4
        movapd nb213_fizH2(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb213_fixH2(%esp)
        movapd %xmm4,nb213_fiyH2(%esp)
        movapd %xmm7,nb213_fizH2(%esp)
        ## update j forces with water H2
        addpd  nb213_fjx(%esp),%xmm0
        addpd  nb213_fjy(%esp),%xmm1
        addpd  nb213_fjz(%esp),%xmm2
        movapd %xmm0,nb213_fjx(%esp)
        movapd %xmm1,nb213_fjy(%esp)
        movapd %xmm2,nb213_fjz(%esp)

        ## M interactions 
        movapd  nb213_rinvM(%esp),%xmm5
        movapd  %xmm5,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        movapd  %xmm5,%xmm7
        movapd  nb213_krsqM(%esp),%xmm0
        addpd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        mulpd   nb213_two(%esp),%xmm0
        subpd   nb213_crf(%esp),%xmm5
        subpd   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulpd   nb213_qqM(%esp),%xmm5   ## vcoul 
        mulpd   nb213_qqM(%esp),%xmm7
        mulpd  %xmm7,%xmm4              ## total fsH2 in xmm4 

        addpd  nb213_vctot(%esp),%xmm5

        movapd nb213_dxM(%esp),%xmm0
        movapd nb213_dyM(%esp),%xmm1
        movapd nb213_dzM(%esp),%xmm2
        movapd %xmm5,nb213_vctot(%esp)
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        ## update H2 forces 
        movapd nb213_fixM(%esp),%xmm3
        movapd nb213_fiyM(%esp),%xmm4
        movapd nb213_fizM(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb213_fixM(%esp)
        movapd %xmm4,nb213_fiyM(%esp)
        movapd %xmm7,nb213_fizM(%esp)

        movl nb213_faction(%ebp),%edi
        ## update j forces 
        addpd  nb213_fjx(%esp),%xmm0
        addpd  nb213_fjy(%esp),%xmm1
        addpd  nb213_fjz(%esp),%xmm2
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
        subl $2,nb213_innerk(%esp)
        jl   _nb_kernel213_ia32_sse2.nb213_checksingle
        jmp  _nb_kernel213_ia32_sse2.nb213_unroll_loop
_nb_kernel213_ia32_sse2.nb213_checksingle: 
        movl  nb213_innerk(%esp),%edx
        andl  $1,%edx
        jnz  _nb_kernel213_ia32_sse2.nb213_dosingle
        jmp  _nb_kernel213_ia32_sse2.nb213_updateouterdata
_nb_kernel213_ia32_sse2.nb213_dosingle: 
        movl  nb213_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb213_innerjjnr(%esp)

        movl nb213_charge(%ebp),%esi     ## base of charge[] 

        xorpd %xmm3,%xmm3
        movlpd (%esi,%eax,8),%xmm3
        movapd %xmm3,%xmm4
        mulsd  nb213_iqM(%esp),%xmm3
        mulsd  nb213_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 

        movapd  %xmm3,nb213_qqM(%esp)
        movapd  %xmm4,nb213_qqH(%esp)

        movl nb213_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb213_vdwparam(%ebp),%esi
        shll %eax
        movl nb213_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movapd %xmm4,nb213_c6(%esp)
        movapd %xmm6,nb213_c12(%esp)

        movl nb213_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb213_ixO(%esp),%xmm4
        movapd nb213_iyO(%esp),%xmm5
        movapd nb213_izO(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb213_dxO(%esp)
        movapd %xmm5,nb213_dyO(%esp)
        movapd %xmm6,nb213_dzO(%esp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb213_ixH1(%esp),%xmm4
        movapd nb213_iyH1(%esp),%xmm5
        movapd nb213_izH1(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb213_dxH1(%esp)
        movapd %xmm5,nb213_dyH1(%esp)
        movapd %xmm6,nb213_dzH1(%esp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm6
        addsd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb213_ixH2(%esp),%xmm3
        movapd nb213_iyH2(%esp),%xmm4
        movapd nb213_izH2(%esp),%xmm5

        ## calc dr 
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5

        ## store dr 
        movapd %xmm3,nb213_dxH2(%esp)
        movapd %xmm4,nb213_dyH2(%esp)
        movapd %xmm5,nb213_dzH2(%esp)
        ## square it 
        mulsd %xmm3,%xmm3
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        addsd %xmm4,%xmm5
        addsd %xmm3,%xmm5
        ## move ixM-izM to xmm2-xmm4  
        movapd nb213_iyM(%esp),%xmm3
        movapd nb213_izM(%esp),%xmm4
        subpd  %xmm1,%xmm3
        subpd  %xmm2,%xmm4
        movapd nb213_ixM(%esp),%xmm2
        subpd  %xmm0,%xmm2

        ## store dr 
        movapd %xmm2,nb213_dxM(%esp)
        movapd %xmm3,nb213_dyM(%esp)
        movapd %xmm4,nb213_dzM(%esp)
        ## square it 
        mulpd %xmm2,%xmm2
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        addpd %xmm3,%xmm4
        addpd %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## calculate krsq
        movsd nb213_krf(%esp),%xmm0
        movsd %xmm0,%xmm1
        movsd %xmm0,%xmm2
        mulsd %xmm4,%xmm0
        mulsd %xmm5,%xmm1
        mulsd %xmm6,%xmm2
        movsd %xmm0,nb213_krsqM(%esp)
        movsd %xmm1,nb213_krsqH2(%esp)
        movsd %xmm2,nb213_krsqH1(%esp)

        ## start with rsqH1 - put seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb213_three(%esp),%xmm1
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb213_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm6       ## rsq*lu*lu 
        movapd nb213_three(%esp),%xmm1
        subsd %xmm6,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb213_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb213_rinvH1(%esp)

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb213_three(%esp),%xmm1
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb213_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm5       ## rsq*lu*lu 
        movapd nb213_three(%esp),%xmm1
        subsd %xmm5,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb213_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb213_rinvH2(%esp)

        ## rsqM - seed in xmm2 
        cvtsd2ss %xmm4,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb213_three(%esp),%xmm1
        mulsd   %xmm4,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb213_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm4       ## rsq*lu*lu 
        movapd nb213_three(%esp),%xmm1
        subsd %xmm4,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb213_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb213_rinvM(%esp)

        ## do O interactions directly. xmm7=rsq
        cvtsd2ss %xmm7,%xmm2
        movapd   %xmm7,%xmm6
        rcpps    %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2
        movapd   nb213_two(%esp),%xmm1
        movapd   %xmm1,%xmm0
        mulsd   %xmm2,%xmm7
        subsd   %xmm7,%xmm1
        mulsd   %xmm1,%xmm2 ## iter1 
        mulsd   %xmm2,%xmm6
        subsd   %xmm6,%xmm0
        mulsd   %xmm2,%xmm0 ## xmm0=rinvsq
        movapd  %xmm0,%xmm1
        mulsd   %xmm1,%xmm1 ## rinv4
        mulsd   %xmm0,%xmm1 ##rinvsix
        movapd  %xmm1,%xmm2
        mulsd   %xmm2,%xmm2 ## rinvtwelve
        mulsd  nb213_c6(%esp),%xmm1
        mulsd  nb213_c12(%esp),%xmm2
        movapd %xmm2,%xmm3
        subsd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addsd  nb213_Vvdwtot(%esp),%xmm3
        mulsd  nb213_six(%esp),%xmm1
        mulsd  nb213_twelve(%esp),%xmm2
        subsd  %xmm1,%xmm2
        mulsd  %xmm0,%xmm2
        movapd %xmm2,%xmm4 ## total fsO 
        movsd %xmm3,nb213_Vvdwtot(%esp)

        movapd nb213_dxO(%esp),%xmm0
        movapd nb213_dyO(%esp),%xmm1
        movapd nb213_dzO(%esp),%xmm2
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update O forces 
        movapd nb213_fixO(%esp),%xmm3
        movapd nb213_fiyO(%esp),%xmm4
        movapd nb213_fizO(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movsd %xmm3,nb213_fixO(%esp)
        movsd %xmm4,nb213_fiyO(%esp)
        movsd %xmm7,nb213_fizO(%esp)
        ## update j forces with water O 
        movsd %xmm0,nb213_fjx(%esp)
        movsd %xmm1,nb213_fjy(%esp)
        movsd %xmm2,nb213_fjz(%esp)

        ## H1 interactions
        movsd  nb213_rinvH1(%esp),%xmm6
        movsd  %xmm6,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        movsd  %xmm6,%xmm7
        movsd  nb213_krsqH1(%esp),%xmm0
        addsd   %xmm0,%xmm6     ## xmm6=rinv+ krsq 
        mulsd   nb213_two(%esp),%xmm0
        subsd   nb213_crf(%esp),%xmm6
        subsd   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulsd   nb213_qqH(%esp),%xmm6   ## vcoul 
        mulsd   nb213_qqH(%esp),%xmm7
        mulsd  %xmm7,%xmm4              ## total fsH1 in xmm4 

        addsd  nb213_vctot(%esp),%xmm6

        movapd nb213_dxH1(%esp),%xmm0
        movapd nb213_dyH1(%esp),%xmm1
        movapd nb213_dzH1(%esp),%xmm2
        movsd %xmm6,nb213_vctot(%esp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb213_fixH1(%esp),%xmm3
        movapd nb213_fiyH1(%esp),%xmm4
        movapd nb213_fizH1(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movsd %xmm3,nb213_fixH1(%esp)
        movsd %xmm4,nb213_fiyH1(%esp)
        movsd %xmm7,nb213_fizH1(%esp)
        ## update j forces with water H1 
        addsd  nb213_fjx(%esp),%xmm0
        addsd  nb213_fjy(%esp),%xmm1
        addsd  nb213_fjz(%esp),%xmm2
        movsd %xmm0,nb213_fjx(%esp)
        movsd %xmm1,nb213_fjy(%esp)
        movsd %xmm2,nb213_fjz(%esp)

        ## H2 interactions 
        movsd  nb213_rinvH2(%esp),%xmm5
        movsd  %xmm5,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        movsd  %xmm5,%xmm7
        movsd  nb213_krsqH2(%esp),%xmm0
        addsd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        mulsd   nb213_two(%esp),%xmm0
        subsd   nb213_crf(%esp),%xmm5
        subsd   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulsd   nb213_qqH(%esp),%xmm5   ## vcoul 
        mulsd   nb213_qqH(%esp),%xmm7
        mulsd  %xmm7,%xmm4              ## total fsH2 in xmm4 

        addsd  nb213_vctot(%esp),%xmm5

        movapd nb213_dxH2(%esp),%xmm0
        movapd nb213_dyH2(%esp),%xmm1
        movapd nb213_dzH2(%esp),%xmm2
        movsd %xmm5,nb213_vctot(%esp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H2 forces 
        movapd nb213_fixH2(%esp),%xmm3
        movapd nb213_fiyH2(%esp),%xmm4
        movapd nb213_fizH2(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movsd %xmm3,nb213_fixH2(%esp)
        movsd %xmm4,nb213_fiyH2(%esp)
        movsd %xmm7,nb213_fizH2(%esp)
        ## update j forces with water H2 
        addsd  nb213_fjx(%esp),%xmm0
        addsd  nb213_fjy(%esp),%xmm1
        addsd  nb213_fjz(%esp),%xmm2
        movsd %xmm0,nb213_fjx(%esp)
        movsd %xmm1,nb213_fjy(%esp)
        movsd %xmm2,nb213_fjz(%esp)

        ## M interactions 
        movsd  nb213_rinvM(%esp),%xmm5
        movsd  %xmm5,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        movsd  %xmm5,%xmm7
        movsd  nb213_krsqM(%esp),%xmm0
        addsd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        mulsd   nb213_two(%esp),%xmm0
        subsd   nb213_crf(%esp),%xmm5
        subsd   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulsd   nb213_qqM(%esp),%xmm5   ## vcoul 
        mulsd   nb213_qqM(%esp),%xmm7
        mulsd  %xmm7,%xmm4              ## total fsH2 in xmm4 

        addsd  nb213_vctot(%esp),%xmm5

        movapd nb213_dxM(%esp),%xmm0
        movapd nb213_dyM(%esp),%xmm1
        movapd nb213_dzM(%esp),%xmm2
        movsd %xmm5,nb213_vctot(%esp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update M forces 
        movapd nb213_fixM(%esp),%xmm3
        movapd nb213_fiyM(%esp),%xmm4
        movapd nb213_fizM(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movsd %xmm3,nb213_fixM(%esp)
        movsd %xmm4,nb213_fiyM(%esp)
        movsd %xmm7,nb213_fizM(%esp)

        movl nb213_faction(%ebp),%edi
        ## update j forces 
        addsd  nb213_fjx(%esp),%xmm0
        addsd  nb213_fjy(%esp),%xmm1
        addsd  nb213_fjz(%esp),%xmm2
        movlpd (%edi,%eax,8),%xmm3
        movlpd 8(%edi,%eax,8),%xmm4
        movlpd 16(%edi,%eax,8),%xmm5
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        movlpd %xmm3,(%edi,%eax,8)
        movlpd %xmm4,8(%edi,%eax,8)
        movlpd %xmm5,16(%edi,%eax,8)

_nb_kernel213_ia32_sse2.nb213_updateouterdata: 
        movl  nb213_ii3(%esp),%ecx
        movl  nb213_faction(%ebp),%edi
        movl  nb213_fshift(%ebp),%esi
        movl  nb213_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movapd nb213_fixO(%esp),%xmm0
        movapd nb213_fiyO(%esp),%xmm1
        movapd nb213_fizO(%esp),%xmm2

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
        movapd nb213_fixH1(%esp),%xmm0
        movapd nb213_fiyH1(%esp),%xmm1
        movapd nb213_fizH1(%esp),%xmm2

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
        movapd nb213_fixH2(%esp),%xmm0
        movapd nb213_fiyH2(%esp),%xmm1
        movapd nb213_fizH2(%esp),%xmm2

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
        movapd nb213_fixM(%esp),%xmm0
        movapd nb213_fiyM(%esp),%xmm1
        movapd nb213_fizM(%esp),%xmm2

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
        movl nb213_n(%esp),%esi
        ## get group index for i particle 
        movl  nb213_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb213_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb213_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb213_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb213_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

       ## finish if last 
        movl nb213_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel213_ia32_sse2.nb213_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb213_n(%esp)
        jmp _nb_kernel213_ia32_sse2.nb213_outer
_nb_kernel213_ia32_sse2.nb213_outerend: 
        ## check if more outer neighborlists remain
        movl  nb213_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel213_ia32_sse2.nb213_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel213_ia32_sse2.nb213_threadloop
_nb_kernel213_ia32_sse2.nb213_end: 
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




.globl nb_kernel213nf_ia32_sse2
.globl _nb_kernel213nf_ia32_sse2
nb_kernel213nf_ia32_sse2:       
_nb_kernel213nf_ia32_sse2:      
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
        ## bottom of stack is cache-aligned for sse2 use 
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
.set nb213nf_iqH, 192
.set nb213nf_iqM, 208
.set nb213nf_qqH, 224
.set nb213nf_qqM, 240
.set nb213nf_c6, 256
.set nb213nf_c12, 272
.set nb213nf_vctot, 288
.set nb213nf_Vvdwtot, 304
.set nb213nf_half, 320
.set nb213nf_three, 336
.set nb213nf_two, 352
.set nb213nf_rinvH1, 368
.set nb213nf_rinvH2, 384
.set nb213nf_rinvM, 400
.set nb213nf_krsqH1, 416
.set nb213nf_krsqH2, 432
.set nb213nf_krsqM, 448
.set nb213nf_krf, 464
.set nb213nf_crf, 480
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
        movsd (%esi),%xmm5
        movsd (%edi),%xmm6
        shufpd $0,%xmm5,%xmm5
        shufpd $0,%xmm6,%xmm6
        movapd %xmm5,nb213nf_krf(%esp)
        movapd %xmm6,nb213nf_crf(%esp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb213nf_half(%esp)
        movl %ebx,nb213nf_half+4(%esp)
        movsd nb213nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb213nf_half(%esp)
        movapd %xmm2,nb213nf_two(%esp)
        movapd %xmm3,nb213nf_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb213nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb213nf_charge(%ebp),%edx
        movsd 8(%edx,%ebx,8),%xmm3
        movsd 24(%edx,%ebx,8),%xmm4
        movl nb213nf_p_facel(%ebp),%esi
        movsd (%esi),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb213nf_iqH(%esp)
        movapd %xmm4,nb213nf_iqM(%esp)

        movl  nb213nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb213nf_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb213nf_ntia(%esp)
_nb_kernel213nf_ia32_sse2.nb213nf_threadloop: 
        movl  nb213nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel213nf_ia32_sse2.nb213nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel213nf_ia32_sse2.nb213nf_spinlock

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
        jg  _nb_kernel213nf_ia32_sse2.nb213nf_outerstart
        jmp _nb_kernel213nf_ia32_sse2.nb213nf_end

_nb_kernel213nf_ia32_sse2.nb213nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb213nf_nouter(%esp),%ebx
        movl %ebx,nb213nf_nouter(%esp)

_nb_kernel213nf_ia32_sse2.nb213nf_outer: 
        movl  nb213nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb213nf_is3(%esp)            ## store is3 

        movl  nb213nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb213nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        movapd %xmm0,%xmm6
        movapd %xmm1,%xmm7

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb213nf_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb213nf_ii3(%esp)

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
        movapd %xmm3,nb213nf_ixO(%esp)
        movapd %xmm4,nb213nf_iyO(%esp)
        movapd %xmm5,nb213nf_izO(%esp)
        movapd %xmm6,nb213nf_ixH1(%esp)
        movapd %xmm7,nb213nf_iyH1(%esp)

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
        movapd %xmm6,nb213nf_izH1(%esp)
        movapd %xmm0,nb213nf_ixH2(%esp)
        movapd %xmm1,nb213nf_iyH2(%esp)
        movapd %xmm2,nb213nf_izH2(%esp)
        movapd %xmm3,nb213nf_ixM(%esp)
        movapd %xmm4,nb213nf_iyM(%esp)
        movapd %xmm5,nb213nf_izM(%esp)

        ## clear vctot
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb213nf_vctot(%esp)
        movapd %xmm4,nb213nf_Vvdwtot(%esp)

        movl  nb213nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb213nf_pos(%ebp),%esi
        movl  nb213nf_faction(%ebp),%edi
        movl  nb213nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb213nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb213nf_ninner(%esp),%ecx
        movl  %ecx,nb213nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb213nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel213nf_ia32_sse2.nb213nf_unroll_loop
        jmp   _nb_kernel213nf_ia32_sse2.nb213nf_checksingle
_nb_kernel213nf_ia32_sse2.nb213nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb213nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb213nf_innerjjnr(%esp)
        ## advance pointer (unrolled 2) 

        movl nb213nf_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb213nf_iqM(%esp),%xmm3
        mulpd  nb213nf_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movapd  %xmm3,nb213nf_qqM(%esp)
        movapd  %xmm4,nb213nf_qqH(%esp)

        movl nb213nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb213nf_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb213nf_ntia(%esp),%edi
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
        movapd %xmm4,nb213nf_c6(%esp)
        movapd %xmm6,nb213nf_c12(%esp)

        movl nb213nf_pos(%ebp),%esi        ## base of pos[] 

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
        movapd nb213nf_ixO(%esp),%xmm4
        movapd nb213nf_iyO(%esp),%xmm5
        movapd nb213nf_izO(%esp),%xmm6

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
        movapd nb213nf_ixH1(%esp),%xmm4
        movapd nb213nf_iyH1(%esp),%xmm5
        movapd nb213nf_izH1(%esp),%xmm6

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
        movapd nb213nf_ixH2(%esp),%xmm3
        movapd nb213nf_iyH2(%esp),%xmm4
        movapd nb213nf_izH2(%esp),%xmm5

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
        movapd nb213nf_iyM(%esp),%xmm3
        movapd nb213nf_izM(%esp),%xmm4
        subpd  %xmm1,%xmm3
        subpd  %xmm2,%xmm4
        movapd nb213nf_ixM(%esp),%xmm2
        subpd  %xmm0,%xmm2

        ## square it 
        mulpd %xmm2,%xmm2
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        addpd %xmm3,%xmm4
        addpd %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## calculate krsq
        movapd nb213nf_krf(%esp),%xmm0
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2
        mulpd %xmm4,%xmm0
        mulpd %xmm5,%xmm1
        mulpd %xmm6,%xmm2
        movapd %xmm0,nb213nf_krsqM(%esp)
        movapd %xmm1,nb213nf_krsqH2(%esp)
        movapd %xmm2,nb213nf_krsqH1(%esp)

        ## start with rsqH1 - put seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb213nf_three(%esp),%xmm1
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulpd   nb213nf_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulpd %xmm1,%xmm1       ## lu*lu 
        mulpd %xmm1,%xmm6       ## rsq*lu*lu 
        movapd nb213nf_three(%esp),%xmm1
        subpd %xmm6,%xmm1       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulpd nb213nf_half(%esp),%xmm1   ## rinv 
        movapd  %xmm1,nb213nf_rinvH1(%esp)

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb213nf_three(%esp),%xmm1
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulpd   nb213nf_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulpd %xmm1,%xmm1       ## lu*lu 
        mulpd %xmm1,%xmm5       ## rsq*lu*lu 
        movapd nb213nf_three(%esp),%xmm1
        subpd %xmm5,%xmm1       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulpd nb213nf_half(%esp),%xmm1   ## rinv 
        movapd  %xmm1,nb213nf_rinvH2(%esp)

        ## rsqM - seed in xmm2 
        cvtpd2ps %xmm4,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb213nf_three(%esp),%xmm1
        mulpd   %xmm4,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulpd   nb213nf_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulpd %xmm1,%xmm1       ## lu*lu 
        mulpd %xmm1,%xmm4       ## rsq*lu*lu 
        movapd nb213nf_three(%esp),%xmm1
        subpd %xmm4,%xmm1       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulpd nb213nf_half(%esp),%xmm1   ## rinv 
        movapd  %xmm1,nb213nf_rinvM(%esp)

        ## do O interactions directly - rsqO is in xmm7
        cvtpd2ps %xmm7,%xmm2
        movapd   %xmm7,%xmm6
        rcpps    %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2
        movapd   nb213nf_two(%esp),%xmm1
        movapd   %xmm1,%xmm0
        mulpd   %xmm2,%xmm7
        subpd   %xmm7,%xmm1
        mulpd   %xmm1,%xmm2 ## iter1 
        mulpd   %xmm2,%xmm6
        subpd   %xmm6,%xmm0
        mulpd   %xmm2,%xmm0 ## xmm0=rinvsq
        movapd  %xmm0,%xmm1
        mulpd   %xmm1,%xmm1 ## rinv4
        mulpd   %xmm0,%xmm1 ##rinvsix
        movapd  %xmm1,%xmm2
        mulpd   %xmm2,%xmm2 ## rinvtwelve
        mulpd  nb213nf_c6(%esp),%xmm1
        mulpd  nb213nf_c12(%esp),%xmm2
        movapd %xmm2,%xmm3
        subpd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addpd  nb213nf_Vvdwtot(%esp),%xmm3
        movapd %xmm3,nb213nf_Vvdwtot(%esp)

        ## H1 interactions 
        movapd  nb213nf_rinvH1(%esp),%xmm6
        addpd   nb213nf_krsqH1(%esp),%xmm6
        subpd   nb213nf_crf(%esp),%xmm6
        mulpd   nb213nf_qqH(%esp),%xmm6   ## vcoul 
        addpd   nb213nf_vctot(%esp),%xmm6
        movapd %xmm6,nb213nf_vctot(%esp)

        ## H2 interactions 
        movapd  nb213nf_rinvH2(%esp),%xmm6
        addpd   nb213nf_krsqH2(%esp),%xmm6
        subpd   nb213nf_crf(%esp),%xmm6
        mulpd   nb213nf_qqH(%esp),%xmm6   ## vcoul 
        addpd   nb213nf_vctot(%esp),%xmm6
        movapd %xmm6,nb213nf_vctot(%esp)

        ## M interactions 
        movapd  nb213nf_rinvM(%esp),%xmm6
        addpd   nb213nf_krsqM(%esp),%xmm6
        subpd   nb213nf_crf(%esp),%xmm6
        mulpd   nb213nf_qqM(%esp),%xmm6   ## vcoul 
        addpd   nb213nf_vctot(%esp),%xmm6
        movapd %xmm6,nb213nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $2,nb213nf_innerk(%esp)
        jl   _nb_kernel213nf_ia32_sse2.nb213nf_checksingle
        jmp  _nb_kernel213nf_ia32_sse2.nb213nf_unroll_loop
_nb_kernel213nf_ia32_sse2.nb213nf_checksingle: 
        movl  nb213nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz  _nb_kernel213nf_ia32_sse2.nb213nf_dosingle
        jmp  _nb_kernel213nf_ia32_sse2.nb213nf_updateouterdata
_nb_kernel213nf_ia32_sse2.nb213nf_dosingle: 
        movl  nb213nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb213nf_innerjjnr(%esp)

        movl nb213nf_charge(%ebp),%esi     ## base of charge[] 

        xorpd %xmm3,%xmm3
        movlpd (%esi,%eax,8),%xmm3
        movapd %xmm3,%xmm4
        mulsd  nb213nf_iqM(%esp),%xmm3
        mulsd  nb213nf_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 

        movapd  %xmm3,nb213nf_qqM(%esp)
        movapd  %xmm4,nb213nf_qqH(%esp)

        movl nb213nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb213nf_vdwparam(%ebp),%esi
        shll %eax
        movl nb213nf_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movapd %xmm4,nb213nf_c6(%esp)
        movapd %xmm6,nb213nf_c12(%esp)

        movl nb213nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb213nf_ixO(%esp),%xmm4
        movapd nb213nf_iyO(%esp),%xmm5
        movapd nb213nf_izO(%esp),%xmm6

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
        movapd nb213nf_ixH1(%esp),%xmm4
        movapd nb213nf_iyH1(%esp),%xmm5
        movapd nb213nf_izH1(%esp),%xmm6

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
        movapd nb213nf_ixH2(%esp),%xmm3
        movapd nb213nf_iyH2(%esp),%xmm4
        movapd nb213nf_izH2(%esp),%xmm5

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
        movapd nb213nf_iyM(%esp),%xmm3
        movapd nb213nf_izM(%esp),%xmm4
        subpd  %xmm1,%xmm3
        subpd  %xmm2,%xmm4
        movapd nb213nf_ixM(%esp),%xmm2
        subpd  %xmm0,%xmm2

        ## square it 
        mulpd %xmm2,%xmm2
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        addpd %xmm3,%xmm4
        addpd %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## calculate krsq
        movsd nb213nf_krf(%esp),%xmm0
        movsd %xmm0,%xmm1
        movsd %xmm0,%xmm2
        mulsd %xmm4,%xmm0
        mulsd %xmm5,%xmm1
        mulsd %xmm6,%xmm2
        movsd %xmm0,nb213nf_krsqM(%esp)
        movsd %xmm1,nb213nf_krsqH2(%esp)
        movsd %xmm2,nb213nf_krsqH1(%esp)

        ## start with rsqH1 - put seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb213nf_three(%esp),%xmm1
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb213nf_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm6       ## rsq*lu*lu 
        movapd nb213nf_three(%esp),%xmm1
        subsd %xmm6,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb213nf_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb213nf_rinvH1(%esp)

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb213nf_three(%esp),%xmm1
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb213nf_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm5       ## rsq*lu*lu 
        movapd nb213nf_three(%esp),%xmm1
        subsd %xmm5,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb213nf_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb213nf_rinvH2(%esp)

        ## rsqM - seed in xmm2 
        cvtsd2ss %xmm4,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb213nf_three(%esp),%xmm1
        mulsd   %xmm4,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb213nf_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm4       ## rsq*lu*lu 
        movapd nb213nf_three(%esp),%xmm1
        subsd %xmm4,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb213nf_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb213nf_rinvM(%esp)

        ## do O interactions directly. xmm7=rsq
        cvtsd2ss %xmm7,%xmm2
        movapd   %xmm7,%xmm6
        rcpps    %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2
        movapd   nb213nf_two(%esp),%xmm1
        movapd   %xmm1,%xmm0
        mulsd   %xmm2,%xmm7
        subsd   %xmm7,%xmm1
        mulsd   %xmm1,%xmm2 ## iter1 
        mulsd   %xmm2,%xmm6
        subsd   %xmm6,%xmm0
        mulsd   %xmm2,%xmm0 ## xmm0=rinvsq
        movapd  %xmm0,%xmm1
        mulsd   %xmm1,%xmm1 ## rinv4
        mulsd   %xmm0,%xmm1 ##rinvsix
        movapd  %xmm1,%xmm2
        mulsd   %xmm2,%xmm2 ## rinvtwelve
        mulsd  nb213nf_c6(%esp),%xmm1
        mulsd  nb213nf_c12(%esp),%xmm2
        movapd %xmm2,%xmm3
        subsd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addsd  nb213nf_Vvdwtot(%esp),%xmm3
        movsd %xmm3,nb213nf_Vvdwtot(%esp)

        ## H1 interactions 
        movsd  nb213nf_rinvH1(%esp),%xmm6
        addsd   nb213nf_krsqH1(%esp),%xmm6
        subsd   nb213nf_crf(%esp),%xmm6
        mulsd   nb213nf_qqH(%esp),%xmm6   ## vcoul 
        addsd   nb213nf_vctot(%esp),%xmm6
        movsd %xmm6,nb213nf_vctot(%esp)

        ## H2 interactions 
        movsd  nb213nf_rinvH2(%esp),%xmm6
        addsd   nb213nf_krsqH2(%esp),%xmm6
        subsd   nb213nf_crf(%esp),%xmm6
        mulsd   nb213nf_qqH(%esp),%xmm6   ## vcoul 
        addsd   nb213nf_vctot(%esp),%xmm6
        movsd %xmm6,nb213nf_vctot(%esp)

        ## M interactions 
        movsd  nb213nf_rinvM(%esp),%xmm6
        addsd   nb213nf_krsqM(%esp),%xmm6
        subsd   nb213nf_crf(%esp),%xmm6
        mulsd   nb213nf_qqM(%esp),%xmm6   ## vcoul 
        addsd   nb213nf_vctot(%esp),%xmm6
        movsd %xmm6,nb213nf_vctot(%esp)

_nb_kernel213nf_ia32_sse2.nb213nf_updateouterdata: 
        ## get n from stack
        movl nb213nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb213nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb213nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb213nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb213nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb213nf_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

       ## finish if last 
        movl nb213nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel213nf_ia32_sse2.nb213nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb213nf_n(%esp)
        jmp _nb_kernel213nf_ia32_sse2.nb213nf_outer
_nb_kernel213nf_ia32_sse2.nb213nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb213nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel213nf_ia32_sse2.nb213nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel213nf_ia32_sse2.nb213nf_threadloop
_nb_kernel213nf_ia32_sse2.nb213nf_end: 
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


