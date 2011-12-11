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



.globl nb_kernel231_ia32_sse2
.globl _nb_kernel231_ia32_sse2
nb_kernel231_ia32_sse2: 
_nb_kernel231_ia32_sse2:        
.set nb231_p_nri, 8
.set nb231_iinr, 12
.set nb231_jindex, 16
.set nb231_jjnr, 20
.set nb231_shift, 24
.set nb231_shiftvec, 28
.set nb231_fshift, 32
.set nb231_gid, 36
.set nb231_pos, 40
.set nb231_faction, 44
.set nb231_charge, 48
.set nb231_p_facel, 52
.set nb231_argkrf, 56
.set nb231_argcrf, 60
.set nb231_Vc, 64
.set nb231_type, 68
.set nb231_p_ntype, 72
.set nb231_vdwparam, 76
.set nb231_Vvdw, 80
.set nb231_p_tabscale, 84
.set nb231_VFtab, 88
.set nb231_invsqrta, 92
.set nb231_dvda, 96
.set nb231_p_gbtabscale, 100
.set nb231_GBtab, 104
.set nb231_p_nthreads, 108
.set nb231_count, 112
.set nb231_mtx, 116
.set nb231_outeriter, 120
.set nb231_inneriter, 124
.set nb231_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb231_ixO, 0
.set nb231_iyO, 16
.set nb231_izO, 32
.set nb231_ixH1, 48
.set nb231_iyH1, 64
.set nb231_izH1, 80
.set nb231_ixH2, 96
.set nb231_iyH2, 112
.set nb231_izH2, 128
.set nb231_iqO, 144
.set nb231_iqH, 160
.set nb231_dxO, 176
.set nb231_dyO, 192
.set nb231_dzO, 208
.set nb231_dxH1, 224
.set nb231_dyH1, 240
.set nb231_dzH1, 256
.set nb231_dxH2, 272
.set nb231_dyH2, 288
.set nb231_dzH2, 304
.set nb231_qqO, 320
.set nb231_qqH, 336
.set nb231_c6, 352
.set nb231_c12, 368
.set nb231_tsc, 384
.set nb231_fstmp, 400
.set nb231_vctot, 416
.set nb231_Vvdwtot, 432
.set nb231_fixO, 448
.set nb231_fiyO, 464
.set nb231_fizO, 480
.set nb231_fixH1, 496
.set nb231_fiyH1, 512
.set nb231_fizH1, 528
.set nb231_fixH2, 544
.set nb231_fiyH2, 560
.set nb231_fizH2, 576
.set nb231_fjx, 592
.set nb231_fjy, 608
.set nb231_fjz, 624
.set nb231_half, 640
.set nb231_three, 656
.set nb231_two, 672
.set nb231_krf, 688
.set nb231_crf, 704
.set nb231_krsqO, 720
.set nb231_krsqH1, 736
.set nb231_krsqH2, 752
.set nb231_rsqO, 768
.set nb231_rinvH1, 784
.set nb231_rinvH2, 800
.set nb231_is3, 816
.set nb231_ii3, 820
.set nb231_ntia, 824
.set nb231_innerjjnr, 828
.set nb231_innerk, 832
.set nb231_n, 836
.set nb231_nn1, 840
.set nb231_nri, 844
.set nb231_nouter, 848
.set nb231_ninner, 852
.set nb231_salign, 856
        pushl %ebp
        movl %esp,%ebp
    pushl %eax
    pushl %ebx
    pushl %ecx
    pushl %edx
        pushl %esi
        pushl %edi
        subl $860,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb231_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb231_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb231_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb231_nouter(%esp)
        movl %eax,nb231_ninner(%esp)

        movl nb231_p_tabscale(%ebp),%eax
        movsd (%eax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb231_tsc(%esp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb231_half(%esp)
        movl %ebx,nb231_half+4(%esp)
        movsd nb231_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb231_half(%esp)
        movapd %xmm2,nb231_two(%esp)
        movapd %xmm3,nb231_three(%esp)

        movl nb231_argkrf(%ebp),%esi
        movl nb231_argcrf(%ebp),%edi
        movsd (%esi),%xmm5
        movsd (%edi),%xmm6
        shufpd $0,%xmm5,%xmm5
        shufpd $0,%xmm6,%xmm6
        movapd %xmm5,nb231_krf(%esp)
        movapd %xmm6,nb231_crf(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb231_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb231_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        movsd 8(%edx,%ebx,8),%xmm4
        movl nb231_p_facel(%ebp),%esi
        movsd (%esi),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb231_iqO(%esp)
        movapd %xmm4,nb231_iqH(%esp)

        movl  nb231_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb231_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb231_ntia(%esp)
_nb_kernel231_ia32_sse2.nb231_threadloop: 
        movl  nb231_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel231_ia32_sse2.nb231_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel231_ia32_sse2.nb231_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb231_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb231_n(%esp)
        movl %ebx,nb231_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel231_ia32_sse2.nb231_outerstart
        jmp _nb_kernel231_ia32_sse2.nb231_end

_nb_kernel231_ia32_sse2.nb231_outerstart: 
        ## ebx contains number of outer iterations
        addl nb231_nouter(%esp),%ebx
        movl %ebx,nb231_nouter(%esp)

_nb_kernel231_ia32_sse2.nb231_outer: 
        movl  nb231_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb231_is3(%esp)      ## store is3 

        movl  nb231_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb231_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb231_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb231_ii3(%esp)

        addsd (%eax,%ebx,8),%xmm3
        addsd 8(%eax,%ebx,8),%xmm4
        addsd 16(%eax,%ebx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb231_ixO(%esp)
        movapd %xmm4,nb231_iyO(%esp)
        movapd %xmm5,nb231_izO(%esp)

        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 24(%eax,%ebx,8),%xmm0
        addsd 32(%eax,%ebx,8),%xmm1
        addsd 40(%eax,%ebx,8),%xmm2
        addsd 48(%eax,%ebx,8),%xmm3
        addsd 56(%eax,%ebx,8),%xmm4
        addsd 64(%eax,%ebx,8),%xmm5

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm0,nb231_ixH1(%esp)
        movapd %xmm1,nb231_iyH1(%esp)
        movapd %xmm2,nb231_izH1(%esp)
        movapd %xmm3,nb231_ixH2(%esp)
        movapd %xmm4,nb231_iyH2(%esp)
        movapd %xmm5,nb231_izH2(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb231_vctot(%esp)
        movapd %xmm4,nb231_Vvdwtot(%esp)
        movapd %xmm4,nb231_fixO(%esp)
        movapd %xmm4,nb231_fiyO(%esp)
        movapd %xmm4,nb231_fizO(%esp)
        movapd %xmm4,nb231_fixH1(%esp)
        movapd %xmm4,nb231_fiyH1(%esp)
        movapd %xmm4,nb231_fizH1(%esp)
        movapd %xmm4,nb231_fixH2(%esp)
        movapd %xmm4,nb231_fiyH2(%esp)
        movapd %xmm4,nb231_fizH2(%esp)

        movl  nb231_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx     ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb231_pos(%ebp),%esi
        movl  nb231_faction(%ebp),%edi
        movl  nb231_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb231_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb231_ninner(%esp),%ecx
        movl  %ecx,nb231_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb231_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel231_ia32_sse2.nb231_unroll_loop
        jmp   _nb_kernel231_ia32_sse2.nb231_checksingle
_nb_kernel231_ia32_sse2.nb231_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb231_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb231_innerjjnr(%esp)                   ## advance pointer (unrolled 2) 

        movl nb231_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb231_iqO(%esp),%xmm3
        mulpd  nb231_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movapd  %xmm3,nb231_qqO(%esp)
        movapd  %xmm4,nb231_qqH(%esp)

        movl nb231_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb231_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb231_ntia(%esp),%edi
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
        movapd %xmm4,nb231_c6(%esp)
        movapd %xmm6,nb231_c12(%esp)

        movl nb231_pos(%ebp),%esi        ## base of pos[] 

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
        movapd nb231_ixO(%esp),%xmm4
        movapd nb231_iyO(%esp),%xmm5
        movapd nb231_izO(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb231_dxO(%esp)
        movapd %xmm5,nb231_dyO(%esp)
        movapd %xmm6,nb231_dzO(%esp)
        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 
        movapd %xmm7,nb231_rsqO(%esp)

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb231_ixH1(%esp),%xmm4
        movapd nb231_iyH1(%esp),%xmm5
        movapd nb231_izH1(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb231_dxH1(%esp)
        movapd %xmm5,nb231_dyH1(%esp)
        movapd %xmm6,nb231_dzH1(%esp)
        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm6
        addpd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb231_ixH2(%esp),%xmm3
        movapd nb231_iyH2(%esp),%xmm4
        movapd nb231_izH2(%esp),%xmm5

        ## calc dr 
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5

        ## store dr 
        movapd %xmm3,nb231_dxH2(%esp)
        movapd %xmm4,nb231_dyH2(%esp)
        movapd %xmm5,nb231_dzH2(%esp)
        ## square it 
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        addpd %xmm4,%xmm5
        addpd %xmm3,%xmm5
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        movapd %xmm5,%xmm0
        movapd %xmm6,%xmm1
        movapd %xmm7,%xmm2

        mulpd  nb231_krf(%esp),%xmm0
        mulpd  nb231_krf(%esp),%xmm1
        mulpd  nb231_krf(%esp),%xmm2

        movapd %xmm0,nb231_krsqH2(%esp)
        movapd %xmm1,nb231_krsqH1(%esp)
        movapd %xmm2,nb231_krsqO(%esp)

        ## start with rsqO - put seed in xmm2 
        cvtpd2ps %xmm7,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb231_three(%esp),%xmm4
        mulpd   %xmm7,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb231_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb231_three(%esp),%xmm4
        subpd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb231_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvO in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb231_three(%esp),%xmm4
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb231_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb231_three(%esp),%xmm4
        subpd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb231_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 
        movapd  %xmm6,nb231_rinvH1(%esp)

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb231_three(%esp),%xmm4
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb231_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb231_three(%esp),%xmm4
        subpd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb231_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 
        movapd  %xmm5,nb231_rinvH2(%esp)

        ## do O interactions 
        movapd %xmm7,%xmm0
        movapd %xmm7,%xmm6
        movapd nb231_krsqO(%esp),%xmm1
        addpd  %xmm1,%xmm0
        mulpd  nb231_two(%esp),%xmm1
        subpd  nb231_crf(%esp),%xmm0   ## xmm0=rinv+ krsq-crf 
        subpd  %xmm1,%xmm6
        mulpd  nb231_qqO(%esp),%xmm0
        mulpd  nb231_qqO(%esp),%xmm6

        mulpd  %xmm7,%xmm6
        movapd %xmm6,nb231_fstmp(%esp)   ## save to temp. storage

        addpd  nb231_vctot(%esp),%xmm0
        movapd %xmm0,nb231_vctot(%esp)

        movapd %xmm7,%xmm0
        movapd nb231_rsqO(%esp),%xmm4
        ## LJ table interaction. xmm0=rinv, xmm4=rsq

        mulpd %xmm0,%xmm4       ## xmm4=r 
        mulpd nb231_tsc(%esp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $3,%mm6           ## idx *= 8 

        movd %eax,%mm0
        movd %ebx,%mm1

        movl nb231_VFtab(%ebp),%esi
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
        mulpd  nb231_two(%esp),%xmm7    ## two*Heps2 
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb231_c6(%esp),%xmm4
        mulpd  %xmm4,%xmm7       ## fijD 
        mulpd  %xmm4,%xmm5       ## Vvdw6 

        ## put scalar force on stack Update Vvdwtot directly 
        addpd  nb231_Vvdwtot(%esp),%xmm5
        movapd nb231_fstmp(%esp),%xmm3
        mulpd  nb231_tsc(%esp),%xmm7
        subpd  %xmm7,%xmm3
        movapd %xmm3,nb231_fstmp(%esp)
        movapd %xmm5,nb231_Vvdwtot(%esp)

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
        mulpd  nb231_two(%esp),%xmm7    ## two*Heps2 
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb231_c12(%esp),%xmm4
        mulpd  %xmm4,%xmm7
        mulpd  %xmm4,%xmm5

        addpd  nb231_Vvdwtot(%esp),%xmm5
        movapd nb231_fstmp(%esp),%xmm3
        mulpd  nb231_tsc(%esp),%xmm7
        subpd  %xmm7,%xmm3
        movapd %xmm5,nb231_Vvdwtot(%esp)

        mulpd  %xmm0,%xmm3

        movapd nb231_dxO(%esp),%xmm0
        movapd nb231_dyO(%esp),%xmm1
        movapd nb231_dzO(%esp),%xmm2

        movd %mm0,%eax
        movd %mm1,%ebx

        movl   nb231_faction(%ebp),%edi
        mulpd  %xmm3,%xmm0
        mulpd  %xmm3,%xmm1
        mulpd  %xmm3,%xmm2

        ## update O forces 
        movapd nb231_fixO(%esp),%xmm3
        movapd nb231_fiyO(%esp),%xmm4
        movapd nb231_fizO(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb231_fixO(%esp)
        movapd %xmm4,nb231_fiyO(%esp)
        movapd %xmm7,nb231_fizO(%esp)
        ## update j forces with water O 
        movapd %xmm0,nb231_fjx(%esp)
        movapd %xmm1,nb231_fjy(%esp)
        movapd %xmm2,nb231_fjz(%esp)

        ## H1 interactions 
        movapd  nb231_rinvH1(%esp),%xmm6
        movapd  %xmm6,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        movapd  %xmm6,%xmm7
        movapd  nb231_krsqH1(%esp),%xmm0
        addpd   %xmm0,%xmm6     ## xmm6=rinv+ krsq 
        mulpd   nb231_two(%esp),%xmm0
        subpd   nb231_crf(%esp),%xmm6
        subpd   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulpd   nb231_qqH(%esp),%xmm6   ## vcoul 
        mulpd   nb231_qqH(%esp),%xmm7
        mulpd  %xmm7,%xmm4              ## total fsH1 in xmm4 

        addpd  nb231_vctot(%esp),%xmm6

        movapd nb231_dxH1(%esp),%xmm0
        movapd nb231_dyH1(%esp),%xmm1
        movapd nb231_dzH1(%esp),%xmm2
        movapd %xmm6,nb231_vctot(%esp)
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb231_fixH1(%esp),%xmm3
        movapd nb231_fiyH1(%esp),%xmm4
        movapd nb231_fizH1(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb231_fixH1(%esp)
        movapd %xmm4,nb231_fiyH1(%esp)
        movapd %xmm7,nb231_fizH1(%esp)
        ## update j forces with water H1 
        addpd  nb231_fjx(%esp),%xmm0
        addpd  nb231_fjy(%esp),%xmm1
        addpd  nb231_fjz(%esp),%xmm2
        movapd %xmm0,nb231_fjx(%esp)
        movapd %xmm1,nb231_fjy(%esp)
        movapd %xmm2,nb231_fjz(%esp)

        ## H2 interactions 
        movapd  nb231_rinvH2(%esp),%xmm5
        movapd  %xmm5,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        movapd  %xmm5,%xmm7
        movapd  nb231_krsqH2(%esp),%xmm0
        addpd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        mulpd   nb231_two(%esp),%xmm0
        subpd   nb231_crf(%esp),%xmm5
        subpd   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulpd   nb231_qqH(%esp),%xmm5   ## vcoul 
        mulpd   nb231_qqH(%esp),%xmm7
        mulpd  %xmm7,%xmm4              ## total fsH2 in xmm4 

        addpd  nb231_vctot(%esp),%xmm5

        movapd nb231_dxH2(%esp),%xmm0
        movapd nb231_dyH2(%esp),%xmm1
        movapd nb231_dzH2(%esp),%xmm2
        movapd %xmm5,nb231_vctot(%esp)
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        ## update H2 forces 
        movapd nb231_fixH2(%esp),%xmm3
        movapd nb231_fiyH2(%esp),%xmm4
        movapd nb231_fizH2(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb231_fixH2(%esp)
        movapd %xmm4,nb231_fiyH2(%esp)
        movapd %xmm7,nb231_fizH2(%esp)

        movl nb231_faction(%ebp),%edi
        ## update j forces 
        addpd  nb231_fjx(%esp),%xmm0
        addpd  nb231_fjy(%esp),%xmm1
        addpd  nb231_fjz(%esp),%xmm2
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
        subl $2,nb231_innerk(%esp)
        jl    _nb_kernel231_ia32_sse2.nb231_checksingle
        jmp   _nb_kernel231_ia32_sse2.nb231_unroll_loop
_nb_kernel231_ia32_sse2.nb231_checksingle: 
        movl  nb231_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel231_ia32_sse2.nb231_dosingle
        jmp   _nb_kernel231_ia32_sse2.nb231_updateouterdata
_nb_kernel231_ia32_sse2.nb231_dosingle: 
        movl  nb231_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb231_innerjjnr(%esp)

        movl nb231_charge(%ebp),%esi     ## base of charge[] 
        xorpd %xmm3,%xmm3
        movlpd (%esi,%eax,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb231_iqO(%esp),%xmm3
        mulpd  nb231_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 

        movapd  %xmm3,nb231_qqO(%esp)
        movapd  %xmm4,nb231_qqH(%esp)

        movl nb231_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb231_vdwparam(%ebp),%esi
        shll %eax
        movl nb231_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 
        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movapd %xmm4,nb231_c6(%esp)
        movapd %xmm6,nb231_c12(%esp)

        movl nb231_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb231_ixO(%esp),%xmm4
        movapd nb231_iyO(%esp),%xmm5
        movapd nb231_izO(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb231_dxO(%esp)
        movapd %xmm5,nb231_dyO(%esp)
        movapd %xmm6,nb231_dzO(%esp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 
        movapd %xmm7,nb231_rsqO(%esp)

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb231_ixH1(%esp),%xmm4
        movapd nb231_iyH1(%esp),%xmm5
        movapd nb231_izH1(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb231_dxH1(%esp)
        movapd %xmm5,nb231_dyH1(%esp)
        movapd %xmm6,nb231_dzH1(%esp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm6
        addsd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb231_ixH2(%esp),%xmm3
        movapd nb231_iyH2(%esp),%xmm4
        movapd nb231_izH2(%esp),%xmm5

        ## calc dr 
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5

        ## store dr 
        movapd %xmm3,nb231_dxH2(%esp)
        movapd %xmm4,nb231_dyH2(%esp)
        movapd %xmm5,nb231_dzH2(%esp)
        ## square it 
        mulsd %xmm3,%xmm3
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        addsd %xmm4,%xmm5
        addsd %xmm3,%xmm5
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        movapd %xmm5,%xmm0
        movapd %xmm6,%xmm1
        movapd %xmm7,%xmm2

        mulsd  nb231_krf(%esp),%xmm0
        mulsd  nb231_krf(%esp),%xmm1
        mulsd  nb231_krf(%esp),%xmm2

        movapd %xmm0,nb231_krsqH2(%esp)
        movapd %xmm1,nb231_krsqH1(%esp)
        movapd %xmm2,nb231_krsqO(%esp)

        ## start with rsqO - put seed in xmm2 
        cvtsd2ss %xmm7,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb231_three(%esp),%xmm4
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb231_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb231_three(%esp),%xmm4
        subsd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb231_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvO in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb231_three(%esp),%xmm4
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb231_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb231_three(%esp),%xmm4
        subsd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb231_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 
        movsd  %xmm6,nb231_rinvH1(%esp)

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb231_three(%esp),%xmm4
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb231_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb231_three(%esp),%xmm4
        subsd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb231_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 
        movsd %xmm5,nb231_rinvH2(%esp)

        ## do O interactions 
        movsd %xmm7,%xmm0
        movsd %xmm7,%xmm6
        movsd nb231_krsqO(%esp),%xmm1
        addsd  %xmm1,%xmm0
        mulsd  nb231_two(%esp),%xmm1
        subsd  nb231_crf(%esp),%xmm0   ## xmm0=rinv+ krsq-crf 
        subsd  %xmm1,%xmm6
        mulsd  nb231_qqO(%esp),%xmm0
        mulsd  nb231_qqO(%esp),%xmm6

        mulsd  %xmm7,%xmm6
        movsd %xmm6,nb231_fstmp(%esp)   ## save to temp. storage

        addsd  nb231_vctot(%esp),%xmm0
        movsd %xmm0,nb231_vctot(%esp)

        movsd %xmm7,%xmm0
        movsd nb231_rsqO(%esp),%xmm4
        ## LJ table interaction. xmm0=rinv, xmm4=rsq

        mulsd %xmm0,%xmm4       ## xmm4=r 
        mulsd nb231_tsc(%esp),%xmm4

        cvttsd2si %xmm4,%ebx    ## mm6 = lu idx 
        cvtsi2sd %ebx,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $3,%ebx

        movl nb231_VFtab(%ebp),%esi

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
        mulsd  nb231_two(%esp),%xmm7    ## two*Heps2 
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movsd nb231_c6(%esp),%xmm4
        mulsd  %xmm4,%xmm7       ## fijD 
        mulsd  %xmm4,%xmm5       ## Vvdw6 

        ## put scalar force on stack Update Vvdwtot directly 
        addsd  nb231_Vvdwtot(%esp),%xmm5
        movsd nb231_fstmp(%esp),%xmm3
        mulsd  nb231_tsc(%esp),%xmm7
        subsd  %xmm7,%xmm3
        movsd %xmm3,nb231_fstmp(%esp)
        movsd %xmm5,nb231_Vvdwtot(%esp)

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
        mulsd  nb231_two(%esp),%xmm7    ## two*Heps2 
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movsd nb231_c12(%esp),%xmm4
        mulsd  %xmm4,%xmm7
        mulsd  %xmm4,%xmm5

        addsd  nb231_Vvdwtot(%esp),%xmm5
        movsd nb231_fstmp(%esp),%xmm3
        mulsd  nb231_tsc(%esp),%xmm7
        subsd  %xmm7,%xmm3
        movsd %xmm5,nb231_Vvdwtot(%esp)

        mulsd  %xmm0,%xmm3

        movsd nb231_dxO(%esp),%xmm0
        movsd nb231_dyO(%esp),%xmm1
        movsd nb231_dzO(%esp),%xmm2

        movl   nb231_faction(%ebp),%edi
        mulsd  %xmm3,%xmm0
        mulsd  %xmm3,%xmm1
        mulsd  %xmm3,%xmm2

        ## update O forces 
        movapd nb231_fixO(%esp),%xmm3
        movapd nb231_fiyO(%esp),%xmm4
        movapd nb231_fizO(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb231_fixO(%esp)
        movlpd %xmm4,nb231_fiyO(%esp)
        movlpd %xmm7,nb231_fizO(%esp)
        ## update j forces with water O 
        movlpd %xmm0,nb231_fjx(%esp)
        movlpd %xmm1,nb231_fjy(%esp)
        movlpd %xmm2,nb231_fjz(%esp)

        ## H1 interactions 
        movsd   nb231_rinvH1(%esp),%xmm6
        movsd   %xmm6,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        movsd   %xmm6,%xmm7
        movsd   nb231_krsqH1(%esp),%xmm0
        addsd   %xmm0,%xmm6     ## xmm6=rinv+ krsq 
        mulsd   nb231_two(%esp),%xmm0
        subsd   nb231_crf(%esp),%xmm6
        subsd   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulsd   nb231_qqH(%esp),%xmm6   ## vcoul 
        mulsd   nb231_qqH(%esp),%xmm7
        mulsd  %xmm7,%xmm4              ## total fsH1 in xmm4 

        addsd  nb231_vctot(%esp),%xmm6

        movapd nb231_dxH1(%esp),%xmm0
        movapd nb231_dyH1(%esp),%xmm1
        movapd nb231_dzH1(%esp),%xmm2
        movlpd %xmm6,nb231_vctot(%esp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb231_fixH1(%esp),%xmm3
        movapd nb231_fiyH1(%esp),%xmm4
        movapd nb231_fizH1(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb231_fixH1(%esp)
        movlpd %xmm4,nb231_fiyH1(%esp)
        movlpd %xmm7,nb231_fizH1(%esp)
        ## update j forces with water H1 
        addsd  nb231_fjx(%esp),%xmm0
        addsd  nb231_fjy(%esp),%xmm1
        addsd  nb231_fjz(%esp),%xmm2
        movlpd %xmm0,nb231_fjx(%esp)
        movlpd %xmm1,nb231_fjy(%esp)
        movlpd %xmm2,nb231_fjz(%esp)

        ## H2 interactions 
        movapd  nb231_rinvH2(%esp),%xmm5
        movapd  %xmm5,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        movapd  %xmm5,%xmm7
        movsd   nb231_krsqH2(%esp),%xmm0
        addsd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        mulsd   nb231_two(%esp),%xmm0
        subsd   nb231_crf(%esp),%xmm5
        subsd   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulsd   nb231_qqH(%esp),%xmm5   ## vcoul 
        mulsd   nb231_qqH(%esp),%xmm7
        mulsd  %xmm7,%xmm4              ## total fsH2 in xmm4 

        addsd  nb231_vctot(%esp),%xmm5

        movapd nb231_dxH2(%esp),%xmm0
        movapd nb231_dyH2(%esp),%xmm1
        movapd nb231_dzH2(%esp),%xmm2
        movlpd %xmm5,nb231_vctot(%esp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H2 forces 
        movapd nb231_fixH2(%esp),%xmm3
        movapd nb231_fiyH2(%esp),%xmm4
        movapd nb231_fizH2(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb231_fixH2(%esp)
        movlpd %xmm4,nb231_fiyH2(%esp)
        movlpd %xmm7,nb231_fizH2(%esp)

        movl nb231_faction(%ebp),%edi
        ## update j forces 
        addsd  nb231_fjx(%esp),%xmm0
        addsd  nb231_fjy(%esp),%xmm1
        addsd  nb231_fjz(%esp),%xmm2
        movlpd (%edi,%eax,8),%xmm3
        movlpd 8(%edi,%eax,8),%xmm4
        movlpd 16(%edi,%eax,8),%xmm5
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        movlpd %xmm3,(%edi,%eax,8)
        movlpd %xmm4,8(%edi,%eax,8)
        movlpd %xmm5,16(%edi,%eax,8)

_nb_kernel231_ia32_sse2.nb231_updateouterdata: 
        movl  nb231_ii3(%esp),%ecx
        movl  nb231_faction(%ebp),%edi
        movl  nb231_fshift(%ebp),%esi
        movl  nb231_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movapd nb231_fixO(%esp),%xmm0
        movapd nb231_fiyO(%esp),%xmm1
        movapd nb231_fizO(%esp),%xmm2

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
        movapd nb231_fixH1(%esp),%xmm0
        movapd nb231_fiyH1(%esp),%xmm1
        movapd nb231_fizH1(%esp),%xmm2

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
        movapd nb231_fixH2(%esp),%xmm0
        movapd nb231_fiyH2(%esp),%xmm1
        movapd nb231_fizH2(%esp),%xmm2

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
        movl nb231_n(%esp),%esi
        ## get group index for i particle 
        movl  nb231_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb231_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb231_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb231_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb231_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

       ## finish if last 
        movl nb231_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel231_ia32_sse2.nb231_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb231_n(%esp)
        jmp _nb_kernel231_ia32_sse2.nb231_outer
_nb_kernel231_ia32_sse2.nb231_outerend: 
        ## check if more outer neighborlists remain
        movl  nb231_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel231_ia32_sse2.nb231_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel231_ia32_sse2.nb231_threadloop
_nb_kernel231_ia32_sse2.nb231_end: 
        emms

        movl nb231_nouter(%esp),%eax
        movl nb231_ninner(%esp),%ebx
        movl nb231_outeriter(%ebp),%ecx
        movl nb231_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb231_salign(%esp),%eax
        addl %eax,%esp
        addl $860,%esp
        popl %edi
        popl %esi
    popl %edx
    popl %ecx
    popl %ebx
    popl %eax
        leave
        ret





.globl nb_kernel231nf_ia32_sse2
.globl _nb_kernel231nf_ia32_sse2
nb_kernel231nf_ia32_sse2:       
_nb_kernel231nf_ia32_sse2:      
.set nb231nf_p_nri, 8
.set nb231nf_iinr, 12
.set nb231nf_jindex, 16
.set nb231nf_jjnr, 20
.set nb231nf_shift, 24
.set nb231nf_shiftvec, 28
.set nb231nf_fshift, 32
.set nb231nf_gid, 36
.set nb231nf_pos, 40
.set nb231nf_faction, 44
.set nb231nf_charge, 48
.set nb231nf_p_facel, 52
.set nb231nf_argkrf, 56
.set nb231nf_argcrf, 60
.set nb231nf_Vc, 64
.set nb231nf_type, 68
.set nb231nf_p_ntype, 72
.set nb231nf_vdwparam, 76
.set nb231nf_Vvdw, 80
.set nb231nf_p_tabscale, 84
.set nb231nf_VFtab, 88
.set nb231nf_invsqrta, 92
.set nb231nf_dvda, 96
.set nb231nf_p_gbtabscale, 100
.set nb231nf_GBtab, 104
.set nb231nf_p_nthreads, 108
.set nb231nf_count, 112
.set nb231nf_mtx, 116
.set nb231nf_outeriter, 120
.set nb231nf_inneriter, 124
.set nb231nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb231nf_ixO, 0
.set nb231nf_iyO, 16
.set nb231nf_izO, 32
.set nb231nf_ixH1, 48
.set nb231nf_iyH1, 64
.set nb231nf_izH1, 80
.set nb231nf_ixH2, 96
.set nb231nf_iyH2, 112
.set nb231nf_izH2, 128
.set nb231nf_iqO, 144
.set nb231nf_iqH, 160
.set nb231nf_qqO, 176
.set nb231nf_qqH, 192
.set nb231nf_c6, 208
.set nb231nf_c12, 224
.set nb231nf_tsc, 240
.set nb231nf_vctot, 256
.set nb231nf_Vvdwtot, 272
.set nb231nf_half, 288
.set nb231nf_three, 304
.set nb231nf_two, 320
.set nb231nf_krf, 336
.set nb231nf_crf, 352
.set nb231nf_krsqO, 368
.set nb231nf_krsqH1, 384
.set nb231nf_krsqH2, 400
.set nb231nf_rsqO, 416
.set nb231nf_rinvH1, 432
.set nb231nf_rinvH2, 448
.set nb231nf_is3, 464
.set nb231nf_ii3, 468
.set nb231nf_ntia, 472
.set nb231nf_innerjjnr, 476
.set nb231nf_innerk, 480
.set nb231nf_n, 484
.set nb231nf_nn1, 488
.set nb231nf_nri, 492
.set nb231nf_nouter, 496
.set nb231nf_ninner, 500
.set nb231nf_salign, 504
        pushl %ebp
        movl %esp,%ebp
    pushl %eax
    pushl %ebx
    pushl %ecx
    pushl %edx
        pushl %esi
        pushl %edi
        subl $508,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb231nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb231nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb231nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb231nf_nouter(%esp)
        movl %eax,nb231nf_ninner(%esp)

        movl nb231nf_p_tabscale(%ebp),%eax
        movsd (%eax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb231nf_tsc(%esp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb231nf_half(%esp)
        movl %ebx,nb231nf_half+4(%esp)
        movsd nb231nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb231nf_half(%esp)
        movapd %xmm2,nb231nf_two(%esp)
        movapd %xmm3,nb231nf_three(%esp)

        movl nb231nf_argkrf(%ebp),%esi
        movl nb231nf_argcrf(%ebp),%edi
        movsd (%esi),%xmm5
        movsd (%edi),%xmm6
        shufpd $0,%xmm5,%xmm5
        shufpd $0,%xmm6,%xmm6
        movapd %xmm5,nb231nf_krf(%esp)
        movapd %xmm6,nb231nf_crf(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb231nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb231nf_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        movsd 8(%edx,%ebx,8),%xmm4
        movl nb231nf_p_facel(%ebp),%esi
        movsd (%esi),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb231nf_iqO(%esp)
        movapd %xmm4,nb231nf_iqH(%esp)

        movl  nb231nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb231nf_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb231nf_ntia(%esp)
_nb_kernel231nf_ia32_sse2.nb231nf_threadloop: 
        movl  nb231nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel231nf_ia32_sse2.nb231nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel231nf_ia32_sse2.nb231nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb231nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb231nf_n(%esp)
        movl %ebx,nb231nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel231nf_ia32_sse2.nb231nf_outerstart
        jmp _nb_kernel231nf_ia32_sse2.nb231nf_end

_nb_kernel231nf_ia32_sse2.nb231nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb231nf_nouter(%esp),%ebx
        movl %ebx,nb231nf_nouter(%esp)

_nb_kernel231nf_ia32_sse2.nb231nf_outer: 
        movl  nb231nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb231nf_is3(%esp)            ## store is3 

        movl  nb231nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb231nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb231nf_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb231nf_ii3(%esp)

        addsd (%eax,%ebx,8),%xmm3
        addsd 8(%eax,%ebx,8),%xmm4
        addsd 16(%eax,%ebx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb231nf_ixO(%esp)
        movapd %xmm4,nb231nf_iyO(%esp)
        movapd %xmm5,nb231nf_izO(%esp)

        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 24(%eax,%ebx,8),%xmm0
        addsd 32(%eax,%ebx,8),%xmm1
        addsd 40(%eax,%ebx,8),%xmm2
        addsd 48(%eax,%ebx,8),%xmm3
        addsd 56(%eax,%ebx,8),%xmm4
        addsd 64(%eax,%ebx,8),%xmm5

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm0,nb231nf_ixH1(%esp)
        movapd %xmm1,nb231nf_iyH1(%esp)
        movapd %xmm2,nb231nf_izH1(%esp)
        movapd %xmm3,nb231nf_ixH2(%esp)
        movapd %xmm4,nb231nf_iyH2(%esp)
        movapd %xmm5,nb231nf_izH2(%esp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb231nf_vctot(%esp)
        movapd %xmm4,nb231nf_Vvdwtot(%esp)

        movl  nb231nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx     ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb231nf_pos(%ebp),%esi
        movl  nb231nf_faction(%ebp),%edi
        movl  nb231nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb231nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb231nf_ninner(%esp),%ecx
        movl  %ecx,nb231nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb231nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel231nf_ia32_sse2.nb231nf_unroll_loop
        jmp   _nb_kernel231nf_ia32_sse2.nb231nf_checksingle
_nb_kernel231nf_ia32_sse2.nb231nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb231nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb231nf_innerjjnr(%esp)                 ## advance pointer (unrolled 2) 

        movl nb231nf_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb231nf_iqO(%esp),%xmm3
        mulpd  nb231nf_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movapd  %xmm3,nb231nf_qqO(%esp)
        movapd  %xmm4,nb231nf_qqH(%esp)

        movl nb231nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb231nf_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb231nf_ntia(%esp),%edi
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
        movapd %xmm4,nb231nf_c6(%esp)
        movapd %xmm6,nb231nf_c12(%esp)

        movl nb231nf_pos(%ebp),%esi        ## base of pos[] 

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
        movapd nb231nf_ixO(%esp),%xmm4
        movapd nb231nf_iyO(%esp),%xmm5
        movapd nb231nf_izO(%esp),%xmm6

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
        movapd %xmm7,nb231nf_rsqO(%esp)

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb231nf_ixH1(%esp),%xmm4
        movapd nb231nf_iyH1(%esp),%xmm5
        movapd nb231nf_izH1(%esp),%xmm6

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
        movapd nb231nf_ixH2(%esp),%xmm3
        movapd nb231nf_iyH2(%esp),%xmm4
        movapd nb231nf_izH2(%esp),%xmm5

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
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        movapd %xmm5,%xmm0
        movapd %xmm6,%xmm1
        movapd %xmm7,%xmm2

        mulpd  nb231nf_krf(%esp),%xmm0
        mulpd  nb231nf_krf(%esp),%xmm1
        mulpd  nb231nf_krf(%esp),%xmm2

        movapd %xmm0,nb231nf_krsqH2(%esp)
        movapd %xmm1,nb231nf_krsqH1(%esp)
        movapd %xmm2,nb231nf_krsqO(%esp)

        ## start with rsqO - put seed in xmm2 
        cvtpd2ps %xmm7,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb231nf_three(%esp),%xmm4
        mulpd   %xmm7,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb231nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb231nf_three(%esp),%xmm4
        subpd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb231nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvO in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb231nf_three(%esp),%xmm4
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb231nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb231nf_three(%esp),%xmm4
        subpd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb231nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 
        movapd  %xmm6,nb231nf_rinvH1(%esp)

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb231nf_three(%esp),%xmm4
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb231nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb231nf_three(%esp),%xmm4
        subpd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb231nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 
        movapd  %xmm5,nb231nf_rinvH2(%esp)

        ## do O interactions 
        movapd %xmm7,%xmm0
        movapd nb231nf_krsqO(%esp),%xmm1
        addpd  %xmm1,%xmm0
        subpd  nb231nf_crf(%esp),%xmm0   ## xmm0=rinv+ krsq-crf 
        mulpd  nb231nf_qqO(%esp),%xmm0

        mulpd  %xmm7,%xmm6

        addpd  nb231nf_vctot(%esp),%xmm0
        movapd %xmm0,nb231nf_vctot(%esp)

        movapd %xmm7,%xmm0
        movapd nb231nf_rsqO(%esp),%xmm4
        ## LJ table interaction. xmm0=rinv, xmm4=rsq

        mulpd %xmm0,%xmm4       ## xmm4=r 
        mulpd nb231nf_tsc(%esp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $3,%mm6           ## idx *= 8 

        movl nb231nf_VFtab(%ebp),%esi
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

        movapd nb231nf_c6(%esp),%xmm4
        mulpd  %xmm4,%xmm5       ## Vvdw6 

        ##  Update Vvdwtot directly 
        addpd  nb231nf_Vvdwtot(%esp),%xmm5
        movapd %xmm5,nb231nf_Vvdwtot(%esp)

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

        movapd nb231nf_c12(%esp),%xmm4
        mulpd  %xmm4,%xmm5

        addpd  nb231nf_Vvdwtot(%esp),%xmm5
        movapd %xmm5,nb231nf_Vvdwtot(%esp)

        ## H1 interactions 
        movapd  nb231nf_rinvH1(%esp),%xmm6
        movapd  %xmm6,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        movapd  %xmm6,%xmm7
        movapd  nb231nf_krsqH1(%esp),%xmm0
        addpd   %xmm0,%xmm6     ## xmm6=rinv+ krsq 
        subpd   nb231nf_crf(%esp),%xmm6
        mulpd   nb231nf_qqH(%esp),%xmm6   ## vcoul 

        addpd  nb231nf_vctot(%esp),%xmm6
        movapd %xmm6,nb231nf_vctot(%esp)

        ## H2 interactions 
        movapd  nb231nf_rinvH2(%esp),%xmm5
        movapd  %xmm5,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        movapd  %xmm5,%xmm7
        movapd  nb231nf_krsqH2(%esp),%xmm0
        addpd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        subpd   nb231nf_crf(%esp),%xmm5
        mulpd   nb231nf_qqH(%esp),%xmm5   ## vcoul 

        addpd  nb231nf_vctot(%esp),%xmm5
        movapd %xmm5,nb231nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $2,nb231nf_innerk(%esp)
        jl    _nb_kernel231nf_ia32_sse2.nb231nf_checksingle
        jmp   _nb_kernel231nf_ia32_sse2.nb231nf_unroll_loop
_nb_kernel231nf_ia32_sse2.nb231nf_checksingle: 
        movl  nb231nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel231nf_ia32_sse2.nb231nf_dosingle
        jmp   _nb_kernel231nf_ia32_sse2.nb231nf_updateouterdata
_nb_kernel231nf_ia32_sse2.nb231nf_dosingle: 
        movl  nb231nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb231nf_innerjjnr(%esp)

        movl nb231nf_charge(%ebp),%esi     ## base of charge[] 
        xorpd %xmm3,%xmm3
        movlpd (%esi,%eax,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb231nf_iqO(%esp),%xmm3
        mulpd  nb231nf_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 

        movapd  %xmm3,nb231nf_qqO(%esp)
        movapd  %xmm4,nb231nf_qqH(%esp)

        movl nb231nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb231nf_vdwparam(%ebp),%esi
        shll %eax
        movl nb231nf_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 
        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movapd %xmm4,nb231nf_c6(%esp)
        movapd %xmm6,nb231nf_c12(%esp)

        movl nb231nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb231nf_ixO(%esp),%xmm4
        movapd nb231nf_iyO(%esp),%xmm5
        movapd nb231nf_izO(%esp),%xmm6

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
        movapd %xmm7,nb231nf_rsqO(%esp)

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb231nf_ixH1(%esp),%xmm4
        movapd nb231nf_iyH1(%esp),%xmm5
        movapd nb231nf_izH1(%esp),%xmm6

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
        movapd nb231nf_ixH2(%esp),%xmm3
        movapd nb231nf_iyH2(%esp),%xmm4
        movapd nb231nf_izH2(%esp),%xmm5

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
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        movapd %xmm5,%xmm0
        movapd %xmm6,%xmm1
        movapd %xmm7,%xmm2

        mulsd  nb231nf_krf(%esp),%xmm0
        mulsd  nb231nf_krf(%esp),%xmm1
        mulsd  nb231nf_krf(%esp),%xmm2

        movapd %xmm0,nb231nf_krsqH2(%esp)
        movapd %xmm1,nb231nf_krsqH1(%esp)
        movapd %xmm2,nb231nf_krsqO(%esp)

        ## start with rsqO - put seed in xmm2 
        cvtsd2ss %xmm7,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb231nf_three(%esp),%xmm4
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb231nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb231nf_three(%esp),%xmm4
        subsd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb231nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvO in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb231nf_three(%esp),%xmm4
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb231nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb231nf_three(%esp),%xmm4
        subsd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb231nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 
        movsd  %xmm6,nb231nf_rinvH1(%esp)

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb231nf_three(%esp),%xmm4
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb231nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb231nf_three(%esp),%xmm4
        subsd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb231nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 
        movsd %xmm5,nb231nf_rinvH2(%esp)

        ## do O interactions 
        movsd %xmm7,%xmm0
        movsd nb231nf_krsqO(%esp),%xmm1
        addsd  %xmm1,%xmm0
        subsd  nb231nf_crf(%esp),%xmm0   ## xmm0=rinv+ krsq-crf 
        mulsd  nb231nf_qqO(%esp),%xmm0

        addsd  nb231nf_vctot(%esp),%xmm0
        movsd %xmm0,nb231nf_vctot(%esp)

        movsd %xmm7,%xmm0
        movsd nb231nf_rsqO(%esp),%xmm4
        ## LJ table interaction. xmm0=rinv, xmm4=rsq

        mulsd %xmm0,%xmm4       ## xmm4=r 
        mulsd nb231nf_tsc(%esp),%xmm4

        cvttsd2si %xmm4,%ebx    ## mm6 = lu idx 
        cvtsi2sd %ebx,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $3,%ebx

        movl nb231nf_VFtab(%ebp),%esi

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

        movsd nb231nf_c6(%esp),%xmm4
        mulsd  %xmm4,%xmm5       ## Vvdw6 

        ## Update Vvdwtot directly 
        addsd  nb231nf_Vvdwtot(%esp),%xmm5
        movsd %xmm5,nb231nf_Vvdwtot(%esp)

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

        movsd nb231nf_c12(%esp),%xmm4
        mulsd  %xmm4,%xmm5

        addsd  nb231nf_Vvdwtot(%esp),%xmm5
        movsd %xmm5,nb231nf_Vvdwtot(%esp)

        ## H1 interactions 
        movsd   nb231nf_rinvH1(%esp),%xmm6
        movsd   %xmm6,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        movsd   nb231nf_krsqH1(%esp),%xmm0
        addsd   %xmm0,%xmm6     ## xmm6=rinv+ krsq 
        subsd   nb231nf_crf(%esp),%xmm6
        mulsd   nb231nf_qqH(%esp),%xmm6   ## vcoul 

        addsd  nb231nf_vctot(%esp),%xmm6
        movlpd %xmm6,nb231nf_vctot(%esp)

        ## H2 interactions 
        movapd  nb231nf_rinvH2(%esp),%xmm5
        movapd  %xmm5,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        movapd  %xmm5,%xmm7
        movsd   nb231nf_krsqH2(%esp),%xmm0
        addsd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        subsd   nb231nf_crf(%esp),%xmm5
        mulsd   nb231nf_qqH(%esp),%xmm5   ## vcoul 

        addsd  nb231nf_vctot(%esp),%xmm5
        movlpd %xmm5,nb231nf_vctot(%esp)

_nb_kernel231nf_ia32_sse2.nb231nf_updateouterdata: 
        ## get n from stack
        movl nb231nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb231nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb231nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb231nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb231nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb231nf_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

       ## finish if last 
        movl nb231nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel231nf_ia32_sse2.nb231nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb231nf_n(%esp)
        jmp _nb_kernel231nf_ia32_sse2.nb231nf_outer
_nb_kernel231nf_ia32_sse2.nb231nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb231nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel231nf_ia32_sse2.nb231nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel231nf_ia32_sse2.nb231nf_threadloop
_nb_kernel231nf_ia32_sse2.nb231nf_end: 
        emms

        movl nb231nf_nouter(%esp),%eax
        movl nb231nf_ninner(%esp),%ebx
        movl nb231nf_outeriter(%ebp),%ecx
        movl nb231nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb231nf_salign(%esp),%eax
        addl %eax,%esp
        addl $508,%esp
        popl %edi
        popl %esi
    popl %edx
    popl %ecx
    popl %ebx
    popl %eax
        leave
        ret



