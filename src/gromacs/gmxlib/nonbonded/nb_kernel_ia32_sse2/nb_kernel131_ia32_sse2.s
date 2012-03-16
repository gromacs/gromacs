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



.globl nb_kernel131_ia32_sse2
.globl _nb_kernel131_ia32_sse2
nb_kernel131_ia32_sse2: 
_nb_kernel131_ia32_sse2:        
.set nb131_p_nri, 8
.set nb131_iinr, 12
.set nb131_jindex, 16
.set nb131_jjnr, 20
.set nb131_shift, 24
.set nb131_shiftvec, 28
.set nb131_fshift, 32
.set nb131_gid, 36
.set nb131_pos, 40
.set nb131_faction, 44
.set nb131_charge, 48
.set nb131_p_facel, 52
.set nb131_argkrf, 56
.set nb131_argcrf, 60
.set nb131_Vc, 64
.set nb131_type, 68
.set nb131_p_ntype, 72
.set nb131_vdwparam, 76
.set nb131_Vvdw, 80
.set nb131_p_tabscale, 84
.set nb131_VFtab, 88
.set nb131_invsqrta, 92
.set nb131_dvda, 96
.set nb131_p_gbtabscale, 100
.set nb131_GBtab, 104
.set nb131_p_nthreads, 108
.set nb131_count, 112
.set nb131_mtx, 116
.set nb131_outeriter, 120
.set nb131_inneriter, 124
.set nb131_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb131_ixO, 0
.set nb131_iyO, 16
.set nb131_izO, 32
.set nb131_ixH1, 48
.set nb131_iyH1, 64
.set nb131_izH1, 80
.set nb131_ixH2, 96
.set nb131_iyH2, 112
.set nb131_izH2, 128
.set nb131_iqO, 144
.set nb131_iqH, 160
.set nb131_dxO, 176
.set nb131_dyO, 192
.set nb131_dzO, 208
.set nb131_dxH1, 224
.set nb131_dyH1, 240
.set nb131_dzH1, 256
.set nb131_dxH2, 272
.set nb131_dyH2, 288
.set nb131_dzH2, 304
.set nb131_qqO, 320
.set nb131_qqH, 336
.set nb131_c6, 352
.set nb131_c12, 368
.set nb131_tsc, 384
.set nb131_fstmp, 400
.set nb131_vctot, 416
.set nb131_Vvdwtot, 432
.set nb131_fixO, 448
.set nb131_fiyO, 464
.set nb131_fizO, 480
.set nb131_fixH1, 496
.set nb131_fiyH1, 512
.set nb131_fizH1, 528
.set nb131_fixH2, 544
.set nb131_fiyH2, 560
.set nb131_fizH2, 576
.set nb131_fjx, 592
.set nb131_fjy, 608
.set nb131_fjz, 624
.set nb131_half, 640
.set nb131_three, 656
.set nb131_two, 672
.set nb131_krsqO, 720
.set nb131_krsqH1, 736
.set nb131_krsqH2, 752
.set nb131_rsqO, 768
.set nb131_rinvH1, 784
.set nb131_rinvH2, 800
.set nb131_is3, 816
.set nb131_ii3, 820
.set nb131_ntia, 824
.set nb131_innerjjnr, 828
.set nb131_innerk, 832
.set nb131_n, 836
.set nb131_nn1, 840
.set nb131_nri, 844
.set nb131_nouter, 848
.set nb131_ninner, 852
.set nb131_salign, 856
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
        movl %eax,nb131_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb131_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb131_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb131_nouter(%esp)
        movl %eax,nb131_ninner(%esp)

        movl nb131_p_tabscale(%ebp),%eax
        movsd (%eax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb131_tsc(%esp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb131_half(%esp)
        movl %ebx,nb131_half+4(%esp)
        movsd nb131_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb131_half(%esp)
        movapd %xmm2,nb131_two(%esp)
        movapd %xmm3,nb131_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb131_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb131_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        movsd 8(%edx,%ebx,8),%xmm4
        movl nb131_p_facel(%ebp),%esi
        movsd (%esi),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb131_iqO(%esp)
        movapd %xmm4,nb131_iqH(%esp)

        movl  nb131_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb131_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb131_ntia(%esp)
_nb_kernel131_ia32_sse2.nb131_threadloop: 
        movl  nb131_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel131_ia32_sse2.nb131_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel131_ia32_sse2.nb131_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb131_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb131_n(%esp)
        movl %ebx,nb131_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel131_ia32_sse2.nb131_outerstart
        jmp _nb_kernel131_ia32_sse2.nb131_end

_nb_kernel131_ia32_sse2.nb131_outerstart: 
        ## ebx contains number of outer iterations
        addl nb131_nouter(%esp),%ebx
        movl %ebx,nb131_nouter(%esp)

_nb_kernel131_ia32_sse2.nb131_outer: 
        movl  nb131_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb131_is3(%esp)      ## store is3 

        movl  nb131_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb131_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb131_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb131_ii3(%esp)

        addsd (%eax,%ebx,8),%xmm3
        addsd 8(%eax,%ebx,8),%xmm4
        addsd 16(%eax,%ebx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb131_ixO(%esp)
        movapd %xmm4,nb131_iyO(%esp)
        movapd %xmm5,nb131_izO(%esp)

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
        movapd %xmm0,nb131_ixH1(%esp)
        movapd %xmm1,nb131_iyH1(%esp)
        movapd %xmm2,nb131_izH1(%esp)
        movapd %xmm3,nb131_ixH2(%esp)
        movapd %xmm4,nb131_iyH2(%esp)
        movapd %xmm5,nb131_izH2(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb131_vctot(%esp)
        movapd %xmm4,nb131_Vvdwtot(%esp)
        movapd %xmm4,nb131_fixO(%esp)
        movapd %xmm4,nb131_fiyO(%esp)
        movapd %xmm4,nb131_fizO(%esp)
        movapd %xmm4,nb131_fixH1(%esp)
        movapd %xmm4,nb131_fiyH1(%esp)
        movapd %xmm4,nb131_fizH1(%esp)
        movapd %xmm4,nb131_fixH2(%esp)
        movapd %xmm4,nb131_fiyH2(%esp)
        movapd %xmm4,nb131_fizH2(%esp)

        movl  nb131_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx     ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb131_pos(%ebp),%esi
        movl  nb131_faction(%ebp),%edi
        movl  nb131_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb131_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb131_ninner(%esp),%ecx
        movl  %ecx,nb131_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb131_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel131_ia32_sse2.nb131_unroll_loop
        jmp   _nb_kernel131_ia32_sse2.nb131_checksingle
_nb_kernel131_ia32_sse2.nb131_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb131_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb131_innerjjnr(%esp)                   ## advance pointer (unrolled 2) 

        movl nb131_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb131_iqO(%esp),%xmm3
        mulpd  nb131_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movapd  %xmm3,nb131_qqO(%esp)
        movapd  %xmm4,nb131_qqH(%esp)

        movl nb131_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb131_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb131_ntia(%esp),%edi
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
        movapd %xmm4,nb131_c6(%esp)
        movapd %xmm6,nb131_c12(%esp)

        movl nb131_pos(%ebp),%esi        ## base of pos[] 

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
        movapd nb131_ixO(%esp),%xmm4
        movapd nb131_iyO(%esp),%xmm5
        movapd nb131_izO(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb131_dxO(%esp)
        movapd %xmm5,nb131_dyO(%esp)
        movapd %xmm6,nb131_dzO(%esp)
        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 
        movapd %xmm7,nb131_rsqO(%esp)

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb131_ixH1(%esp),%xmm4
        movapd nb131_iyH1(%esp),%xmm5
        movapd nb131_izH1(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb131_dxH1(%esp)
        movapd %xmm5,nb131_dyH1(%esp)
        movapd %xmm6,nb131_dzH1(%esp)
        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm6
        addpd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb131_ixH2(%esp),%xmm3
        movapd nb131_iyH2(%esp),%xmm4
        movapd nb131_izH2(%esp),%xmm5

        ## calc dr 
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5

        ## store dr 
        movapd %xmm3,nb131_dxH2(%esp)
        movapd %xmm4,nb131_dyH2(%esp)
        movapd %xmm5,nb131_dzH2(%esp)
        ## square it 
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        addpd %xmm4,%xmm5
        addpd %xmm3,%xmm5
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## start with rsqO - put seed in xmm2 
        cvtpd2ps %xmm7,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb131_three(%esp),%xmm4
        mulpd   %xmm7,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb131_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb131_three(%esp),%xmm4
        subpd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb131_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvO in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb131_three(%esp),%xmm4
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb131_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb131_three(%esp),%xmm4
        subpd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb131_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 
        movapd  %xmm6,nb131_rinvH1(%esp)

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb131_three(%esp),%xmm4
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb131_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb131_three(%esp),%xmm4
        subpd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb131_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 
        movapd  %xmm5,nb131_rinvH2(%esp)

        ## do O interactions 
        movapd %xmm7,%xmm0
        mulpd  nb131_qqO(%esp),%xmm7   ## vcoul
        movapd %xmm0,%xmm6
        mulpd  %xmm7,%xmm6 ## vcoul*rinv

        movapd %xmm6,nb131_fstmp(%esp)   ## save to temp. storage

        addpd  nb131_vctot(%esp),%xmm7
        movapd %xmm7,nb131_vctot(%esp)

        movapd nb131_rsqO(%esp),%xmm4
        ## LJ table interaction. xmm0=rinv, xmm4=rsq
        mulpd %xmm0,%xmm4       ## xmm4=r 
        mulpd nb131_tsc(%esp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $3,%mm6           ## idx *= 8 

        movd %eax,%mm0
        movd %ebx,%mm1

        movl nb131_VFtab(%ebp),%esi
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
        mulpd  nb131_two(%esp),%xmm7    ## two*Heps2 
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb131_c6(%esp),%xmm4
        mulpd  %xmm4,%xmm7       ## fijD 
        mulpd  %xmm4,%xmm5       ## Vvdw6 

        ## put scalar force on stack Update Vvdwtot directly 
        addpd  nb131_Vvdwtot(%esp),%xmm5
        movapd nb131_fstmp(%esp),%xmm3
        mulpd  nb131_tsc(%esp),%xmm7
        subpd  %xmm7,%xmm3
        movapd %xmm3,nb131_fstmp(%esp)
        movapd %xmm5,nb131_Vvdwtot(%esp)

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
        mulpd  nb131_two(%esp),%xmm7    ## two*Heps2 
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb131_c12(%esp),%xmm4
        mulpd  %xmm4,%xmm7
        mulpd  %xmm4,%xmm5

        addpd  nb131_Vvdwtot(%esp),%xmm5
        movapd nb131_fstmp(%esp),%xmm3
        mulpd  nb131_tsc(%esp),%xmm7
        subpd  %xmm7,%xmm3
        movapd %xmm5,nb131_Vvdwtot(%esp)

        mulpd  %xmm0,%xmm3

        movapd nb131_dxO(%esp),%xmm0
        movapd nb131_dyO(%esp),%xmm1
        movapd nb131_dzO(%esp),%xmm2

        movd %mm0,%eax
        movd %mm1,%ebx

        movl   nb131_faction(%ebp),%edi
        mulpd  %xmm3,%xmm0
        mulpd  %xmm3,%xmm1
        mulpd  %xmm3,%xmm2

        ## update O forces 
        movapd nb131_fixO(%esp),%xmm3
        movapd nb131_fiyO(%esp),%xmm4
        movapd nb131_fizO(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb131_fixO(%esp)
        movapd %xmm4,nb131_fiyO(%esp)
        movapd %xmm7,nb131_fizO(%esp)
        ## update j forces with water O 
        movapd %xmm0,nb131_fjx(%esp)
        movapd %xmm1,nb131_fjy(%esp)
        movapd %xmm2,nb131_fjz(%esp)

        ## H1 interactions 
        movapd  nb131_rinvH1(%esp),%xmm6
        movapd  %xmm6,%xmm4
        mulpd   nb131_qqH(%esp),%xmm6   ## vcoul 
        mulpd   %xmm4,%xmm4 ## rinvsq
        mulpd   %xmm6,%xmm4 ## vcoul*rinvsq

        addpd   nb131_vctot(%esp),%xmm6
        movapd %xmm6,nb131_vctot(%esp)

        movapd nb131_dxH1(%esp),%xmm0
        movapd nb131_dyH1(%esp),%xmm1
        movapd nb131_dzH1(%esp),%xmm2

        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb131_fixH1(%esp),%xmm3
        movapd nb131_fiyH1(%esp),%xmm4
        movapd nb131_fizH1(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb131_fixH1(%esp)
        movapd %xmm4,nb131_fiyH1(%esp)
        movapd %xmm7,nb131_fizH1(%esp)
        ## update j forces with water H1 
        addpd  nb131_fjx(%esp),%xmm0
        addpd  nb131_fjy(%esp),%xmm1
        addpd  nb131_fjz(%esp),%xmm2
        movapd %xmm0,nb131_fjx(%esp)
        movapd %xmm1,nb131_fjy(%esp)
        movapd %xmm2,nb131_fjz(%esp)

        ## H2 interactions 
        movapd  nb131_rinvH2(%esp),%xmm6
        movapd  %xmm6,%xmm4
        mulpd   nb131_qqH(%esp),%xmm6   ## vcoul 
        mulpd   %xmm4,%xmm4 ## rinvsq
        mulpd   %xmm6,%xmm4 ## vcoul*rinvsq

        addpd  nb131_vctot(%esp),%xmm6

        movapd nb131_dxH2(%esp),%xmm0
        movapd nb131_dyH2(%esp),%xmm1
        movapd nb131_dzH2(%esp),%xmm2
        movapd %xmm6,nb131_vctot(%esp)
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        ## update H2 forces 
        movapd nb131_fixH2(%esp),%xmm3
        movapd nb131_fiyH2(%esp),%xmm4
        movapd nb131_fizH2(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb131_fixH2(%esp)
        movapd %xmm4,nb131_fiyH2(%esp)
        movapd %xmm7,nb131_fizH2(%esp)

        movl nb131_faction(%ebp),%edi
        ## update j forces 
        addpd  nb131_fjx(%esp),%xmm0
        addpd  nb131_fjy(%esp),%xmm1
        addpd  nb131_fjz(%esp),%xmm2
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
        subl $2,nb131_innerk(%esp)
        jl    _nb_kernel131_ia32_sse2.nb131_checksingle
        jmp   _nb_kernel131_ia32_sse2.nb131_unroll_loop
_nb_kernel131_ia32_sse2.nb131_checksingle: 
        movl  nb131_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel131_ia32_sse2.nb131_dosingle
        jmp   _nb_kernel131_ia32_sse2.nb131_updateouterdata
_nb_kernel131_ia32_sse2.nb131_dosingle: 
        movl  nb131_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb131_innerjjnr(%esp)

        movl nb131_charge(%ebp),%esi     ## base of charge[] 
        xorpd %xmm3,%xmm3
        movlpd (%esi,%eax,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb131_iqO(%esp),%xmm3
        mulpd  nb131_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 

        movapd  %xmm3,nb131_qqO(%esp)
        movapd  %xmm4,nb131_qqH(%esp)

        movl nb131_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb131_vdwparam(%ebp),%esi
        shll %eax
        movl nb131_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 
        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movapd %xmm4,nb131_c6(%esp)
        movapd %xmm6,nb131_c12(%esp)

        movl nb131_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb131_ixO(%esp),%xmm4
        movapd nb131_iyO(%esp),%xmm5
        movapd nb131_izO(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb131_dxO(%esp)
        movapd %xmm5,nb131_dyO(%esp)
        movapd %xmm6,nb131_dzO(%esp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 
        movapd %xmm7,nb131_rsqO(%esp)

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb131_ixH1(%esp),%xmm4
        movapd nb131_iyH1(%esp),%xmm5
        movapd nb131_izH1(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb131_dxH1(%esp)
        movapd %xmm5,nb131_dyH1(%esp)
        movapd %xmm6,nb131_dzH1(%esp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm6
        addsd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb131_ixH2(%esp),%xmm3
        movapd nb131_iyH2(%esp),%xmm4
        movapd nb131_izH2(%esp),%xmm5

        ## calc dr 
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5

        ## store dr 
        movapd %xmm3,nb131_dxH2(%esp)
        movapd %xmm4,nb131_dyH2(%esp)
        movapd %xmm5,nb131_dzH2(%esp)
        ## square it 
        mulsd %xmm3,%xmm3
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        addsd %xmm4,%xmm5
        addsd %xmm3,%xmm5
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## start with rsqO - put seed in xmm2 
        cvtsd2ss %xmm7,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb131_three(%esp),%xmm4
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb131_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb131_three(%esp),%xmm4
        subsd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb131_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvO in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb131_three(%esp),%xmm4
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb131_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb131_three(%esp),%xmm4
        subsd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb131_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 
        movsd  %xmm6,nb131_rinvH1(%esp)

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb131_three(%esp),%xmm4
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb131_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb131_three(%esp),%xmm4
        subsd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb131_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 
        movsd %xmm5,nb131_rinvH2(%esp)

        ## do O interactions 
        movsd %xmm7,%xmm0
        mulsd  nb131_qqO(%esp),%xmm7   ## vcoul
        movsd %xmm0,%xmm6
        mulsd  %xmm7,%xmm6 ## vcoul*rinv

        movsd %xmm6,nb131_fstmp(%esp)   ## save to temp. storage

        addsd  nb131_vctot(%esp),%xmm7
        movsd %xmm7,nb131_vctot(%esp)

        movsd nb131_rsqO(%esp),%xmm4

        ## LJ table interaction. xmm0=rinv, xmm4=rsq    
        mulsd %xmm0,%xmm4       ## xmm4=r 
        mulsd nb131_tsc(%esp),%xmm4

        cvttsd2si %xmm4,%ebx    ## mm6 = lu idx 
        cvtsi2sd %ebx,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $3,%ebx

        movl nb131_VFtab(%ebp),%esi

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
        mulsd  nb131_two(%esp),%xmm7    ## two*Heps2 
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movsd nb131_c6(%esp),%xmm4
        mulsd  %xmm4,%xmm7       ## fijD 
        mulsd  %xmm4,%xmm5       ## Vvdw6 

        ## put scalar force on stack Update Vvdwtot directly 
        addsd  nb131_Vvdwtot(%esp),%xmm5
        movsd nb131_fstmp(%esp),%xmm3
        mulsd  nb131_tsc(%esp),%xmm7
        subsd  %xmm7,%xmm3
        movsd %xmm3,nb131_fstmp(%esp)
        movsd %xmm5,nb131_Vvdwtot(%esp)

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
        mulsd  nb131_two(%esp),%xmm7    ## two*Heps2 
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movsd nb131_c12(%esp),%xmm4
        mulsd  %xmm4,%xmm7
        mulsd  %xmm4,%xmm5

        addsd  nb131_Vvdwtot(%esp),%xmm5
        movsd nb131_fstmp(%esp),%xmm3
        mulsd  nb131_tsc(%esp),%xmm7
        subsd  %xmm7,%xmm3
        movsd %xmm5,nb131_Vvdwtot(%esp)

        mulsd  %xmm0,%xmm3

        movsd nb131_dxO(%esp),%xmm0
        movsd nb131_dyO(%esp),%xmm1
        movsd nb131_dzO(%esp),%xmm2

        movl   nb131_faction(%ebp),%edi
        mulsd  %xmm3,%xmm0
        mulsd  %xmm3,%xmm1
        mulsd  %xmm3,%xmm2

        ## update O forces 
        movapd nb131_fixO(%esp),%xmm3
        movapd nb131_fiyO(%esp),%xmm4
        movapd nb131_fizO(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb131_fixO(%esp)
        movlpd %xmm4,nb131_fiyO(%esp)
        movlpd %xmm7,nb131_fizO(%esp)
        ## update j forces with water O 
        movlpd %xmm0,nb131_fjx(%esp)
        movlpd %xmm1,nb131_fjy(%esp)
        movlpd %xmm2,nb131_fjz(%esp)

        ## H1 interactions 
        movsd  nb131_rinvH1(%esp),%xmm6
        movsd  %xmm6,%xmm4
        mulsd   nb131_qqH(%esp),%xmm6   ## vcoul 
        mulsd   %xmm4,%xmm4 ## rinvsq
        mulsd   %xmm6,%xmm4 ## vcoul*rinvsq

        addsd  nb131_vctot(%esp),%xmm6

        movapd nb131_dxH1(%esp),%xmm0
        movapd nb131_dyH1(%esp),%xmm1
        movapd nb131_dzH1(%esp),%xmm2
        movlpd %xmm6,nb131_vctot(%esp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb131_fixH1(%esp),%xmm3
        movapd nb131_fiyH1(%esp),%xmm4
        movapd nb131_fizH1(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb131_fixH1(%esp)
        movlpd %xmm4,nb131_fiyH1(%esp)
        movlpd %xmm7,nb131_fizH1(%esp)
        ## update j forces with water H1 
        addsd  nb131_fjx(%esp),%xmm0
        addsd  nb131_fjy(%esp),%xmm1
        addsd  nb131_fjz(%esp),%xmm2
        movlpd %xmm0,nb131_fjx(%esp)
        movlpd %xmm1,nb131_fjy(%esp)
        movlpd %xmm2,nb131_fjz(%esp)

        ## H2 interactions 
        movsd  nb131_rinvH2(%esp),%xmm6
        movsd  %xmm6,%xmm4
        mulsd   nb131_qqH(%esp),%xmm6   ## vcoul 
        mulsd   %xmm4,%xmm4 ## rinvsq
        mulsd   %xmm6,%xmm4 ## vcoul*rinvsq

        addsd  nb131_vctot(%esp),%xmm6

        movapd nb131_dxH2(%esp),%xmm0
        movapd nb131_dyH2(%esp),%xmm1
        movapd nb131_dzH2(%esp),%xmm2
        movlpd %xmm6,nb131_vctot(%esp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H2 forces 
        movapd nb131_fixH2(%esp),%xmm3
        movapd nb131_fiyH2(%esp),%xmm4
        movapd nb131_fizH2(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb131_fixH2(%esp)
        movlpd %xmm4,nb131_fiyH2(%esp)
        movlpd %xmm7,nb131_fizH2(%esp)

        movl nb131_faction(%ebp),%edi
        ## update j forces 
        addsd  nb131_fjx(%esp),%xmm0
        addsd  nb131_fjy(%esp),%xmm1
        addsd  nb131_fjz(%esp),%xmm2
        movlpd (%edi,%eax,8),%xmm3
        movlpd 8(%edi,%eax,8),%xmm4
        movlpd 16(%edi,%eax,8),%xmm5
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        movlpd %xmm3,(%edi,%eax,8)
        movlpd %xmm4,8(%edi,%eax,8)
        movlpd %xmm5,16(%edi,%eax,8)

_nb_kernel131_ia32_sse2.nb131_updateouterdata: 
        movl  nb131_ii3(%esp),%ecx
        movl  nb131_faction(%ebp),%edi
        movl  nb131_fshift(%ebp),%esi
        movl  nb131_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movapd nb131_fixO(%esp),%xmm0
        movapd nb131_fiyO(%esp),%xmm1
        movapd nb131_fizO(%esp),%xmm2

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
        movapd nb131_fixH1(%esp),%xmm0
        movapd nb131_fiyH1(%esp),%xmm1
        movapd nb131_fizH1(%esp),%xmm2

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
        movapd nb131_fixH2(%esp),%xmm0
        movapd nb131_fiyH2(%esp),%xmm1
        movapd nb131_fizH2(%esp),%xmm2

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
        movl nb131_n(%esp),%esi
        ## get group index for i particle 
        movl  nb131_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb131_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb131_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb131_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb131_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

       ## finish if last 
        movl nb131_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel131_ia32_sse2.nb131_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb131_n(%esp)
        jmp _nb_kernel131_ia32_sse2.nb131_outer
_nb_kernel131_ia32_sse2.nb131_outerend: 
        ## check if more outer neighborlists remain
        movl  nb131_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel131_ia32_sse2.nb131_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel131_ia32_sse2.nb131_threadloop
_nb_kernel131_ia32_sse2.nb131_end: 
        emms

        movl nb131_nouter(%esp),%eax
        movl nb131_ninner(%esp),%ebx
        movl nb131_outeriter(%ebp),%ecx
        movl nb131_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb131_salign(%esp),%eax
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





.globl nb_kernel131nf_ia32_sse2
.globl _nb_kernel131nf_ia32_sse2
nb_kernel131nf_ia32_sse2:       
_nb_kernel131nf_ia32_sse2:      
.set nb131nf_p_nri, 8
.set nb131nf_iinr, 12
.set nb131nf_jindex, 16
.set nb131nf_jjnr, 20
.set nb131nf_shift, 24
.set nb131nf_shiftvec, 28
.set nb131nf_fshift, 32
.set nb131nf_gid, 36
.set nb131nf_pos, 40
.set nb131nf_faction, 44
.set nb131nf_charge, 48
.set nb131nf_p_facel, 52
.set nb131nf_argkrf, 56
.set nb131nf_argcrf, 60
.set nb131nf_Vc, 64
.set nb131nf_type, 68
.set nb131nf_p_ntype, 72
.set nb131nf_vdwparam, 76
.set nb131nf_Vvdw, 80
.set nb131nf_p_tabscale, 84
.set nb131nf_VFtab, 88
.set nb131nf_invsqrta, 92
.set nb131nf_dvda, 96
.set nb131nf_p_gbtabscale, 100
.set nb131nf_GBtab, 104
.set nb131nf_p_nthreads, 108
.set nb131nf_count, 112
.set nb131nf_mtx, 116
.set nb131nf_outeriter, 120
.set nb131nf_inneriter, 124
.set nb131nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb131nf_ixO, 0
.set nb131nf_iyO, 16
.set nb131nf_izO, 32
.set nb131nf_ixH1, 48
.set nb131nf_iyH1, 64
.set nb131nf_izH1, 80
.set nb131nf_ixH2, 96
.set nb131nf_iyH2, 112
.set nb131nf_izH2, 128
.set nb131nf_iqO, 144
.set nb131nf_iqH, 160
.set nb131nf_qqO, 176
.set nb131nf_qqH, 192
.set nb131nf_c6, 208
.set nb131nf_c12, 224
.set nb131nf_tsc, 240
.set nb131nf_vctot, 256
.set nb131nf_Vvdwtot, 272
.set nb131nf_half, 288
.set nb131nf_three, 304
.set nb131nf_two, 320
.set nb131nf_krsqO, 368
.set nb131nf_krsqH1, 384
.set nb131nf_krsqH2, 400
.set nb131nf_rsqO, 416
.set nb131nf_rinvH1, 432
.set nb131nf_rinvH2, 448
.set nb131nf_is3, 464
.set nb131nf_ii3, 468
.set nb131nf_ntia, 472
.set nb131nf_innerjjnr, 476
.set nb131nf_innerk, 480
.set nb131nf_n, 484
.set nb131nf_nn1, 488
.set nb131nf_nri, 492
.set nb131nf_nouter, 496
.set nb131nf_ninner, 500
.set nb131nf_salign, 504
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
        movl %eax,nb131nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb131nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb131nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb131nf_nouter(%esp)
        movl %eax,nb131nf_ninner(%esp)

        movl nb131nf_p_tabscale(%ebp),%eax
        movsd (%eax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb131nf_tsc(%esp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb131nf_half(%esp)
        movl %ebx,nb131nf_half+4(%esp)
        movsd nb131nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb131nf_half(%esp)
        movapd %xmm2,nb131nf_two(%esp)
        movapd %xmm3,nb131nf_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb131nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb131nf_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        movsd 8(%edx,%ebx,8),%xmm4
        movl nb131nf_p_facel(%ebp),%esi
        movsd (%esi),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb131nf_iqO(%esp)
        movapd %xmm4,nb131nf_iqH(%esp)

        movl  nb131nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb131nf_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb131nf_ntia(%esp)
_nb_kernel131nf_ia32_sse2.nb131nf_threadloop: 
        movl  nb131nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel131nf_ia32_sse2.nb131nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel131nf_ia32_sse2.nb131nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb131nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb131nf_n(%esp)
        movl %ebx,nb131nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel131nf_ia32_sse2.nb131nf_outerstart
        jmp _nb_kernel131nf_ia32_sse2.nb131nf_end

_nb_kernel131nf_ia32_sse2.nb131nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb131nf_nouter(%esp),%ebx
        movl %ebx,nb131nf_nouter(%esp)

_nb_kernel131nf_ia32_sse2.nb131nf_outer: 
        movl  nb131nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb131nf_is3(%esp)            ## store is3 

        movl  nb131nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb131nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb131nf_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb131nf_ii3(%esp)

        addsd (%eax,%ebx,8),%xmm3
        addsd 8(%eax,%ebx,8),%xmm4
        addsd 16(%eax,%ebx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb131nf_ixO(%esp)
        movapd %xmm4,nb131nf_iyO(%esp)
        movapd %xmm5,nb131nf_izO(%esp)

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
        movapd %xmm0,nb131nf_ixH1(%esp)
        movapd %xmm1,nb131nf_iyH1(%esp)
        movapd %xmm2,nb131nf_izH1(%esp)
        movapd %xmm3,nb131nf_ixH2(%esp)
        movapd %xmm4,nb131nf_iyH2(%esp)
        movapd %xmm5,nb131nf_izH2(%esp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb131nf_vctot(%esp)
        movapd %xmm4,nb131nf_Vvdwtot(%esp)

        movl  nb131nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx     ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb131nf_pos(%ebp),%esi
        movl  nb131nf_faction(%ebp),%edi
        movl  nb131nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb131nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb131nf_ninner(%esp),%ecx
        movl  %ecx,nb131nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb131nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel131nf_ia32_sse2.nb131nf_unroll_loop
        jmp   _nb_kernel131nf_ia32_sse2.nb131nf_checksingle
_nb_kernel131nf_ia32_sse2.nb131nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb131nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb131nf_innerjjnr(%esp)                 ## advance pointer (unrolled 2) 

        movl nb131nf_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb131nf_iqO(%esp),%xmm3
        mulpd  nb131nf_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movapd  %xmm3,nb131nf_qqO(%esp)
        movapd  %xmm4,nb131nf_qqH(%esp)

        movl nb131nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb131nf_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb131nf_ntia(%esp),%edi
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
        movapd %xmm4,nb131nf_c6(%esp)
        movapd %xmm6,nb131nf_c12(%esp)

        movl nb131nf_pos(%ebp),%esi        ## base of pos[] 

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
        movapd nb131nf_ixO(%esp),%xmm4
        movapd nb131nf_iyO(%esp),%xmm5
        movapd nb131nf_izO(%esp),%xmm6

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
        movapd %xmm7,nb131nf_rsqO(%esp)

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb131nf_ixH1(%esp),%xmm4
        movapd nb131nf_iyH1(%esp),%xmm5
        movapd nb131nf_izH1(%esp),%xmm6

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
        movapd nb131nf_ixH2(%esp),%xmm3
        movapd nb131nf_iyH2(%esp),%xmm4
        movapd nb131nf_izH2(%esp),%xmm5

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

        ## start with rsqO - put seed in xmm2 
        cvtpd2ps %xmm7,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb131nf_three(%esp),%xmm4
        mulpd   %xmm7,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb131nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb131nf_three(%esp),%xmm4
        subpd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb131nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvO in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb131nf_three(%esp),%xmm4
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb131nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb131nf_three(%esp),%xmm4
        subpd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb131nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 
        movapd  %xmm6,nb131nf_rinvH1(%esp)

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb131nf_three(%esp),%xmm4
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb131nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb131nf_three(%esp),%xmm4
        subpd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb131nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 
        movapd  %xmm5,nb131nf_rinvH2(%esp)

        ## do O interactions 
        movapd %xmm7,%xmm0
        mulpd  nb131nf_qqO(%esp),%xmm7

        addpd  nb131nf_vctot(%esp),%xmm7
        movapd %xmm7,nb131nf_vctot(%esp)

        movapd nb131nf_rsqO(%esp),%xmm4
        ## LJ table interaction. xmm0=rinv, xmm4=rsq

        mulpd %xmm0,%xmm4       ## xmm4=r 
        mulpd nb131nf_tsc(%esp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $3,%mm6           ## idx *= 8 

        movl nb131nf_VFtab(%ebp),%esi
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

        movapd nb131nf_c6(%esp),%xmm4
        mulpd  %xmm4,%xmm5       ## Vvdw6 

        ##  Update Vvdwtot directly 
        addpd  nb131nf_Vvdwtot(%esp),%xmm5
        movapd %xmm5,nb131nf_Vvdwtot(%esp)

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

        movapd nb131nf_c12(%esp),%xmm4
        mulpd  %xmm4,%xmm5

        addpd  nb131nf_Vvdwtot(%esp),%xmm5
        movapd %xmm5,nb131nf_Vvdwtot(%esp)

        ## H1 & H2 interactions 
        movapd  nb131nf_rinvH1(%esp),%xmm6
        addpd   nb131nf_rinvH2(%esp),%xmm6
        mulpd   nb131nf_qqH(%esp),%xmm6   ## vcoul 

        addpd  nb131nf_vctot(%esp),%xmm6
        movapd %xmm6,nb131nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $2,nb131nf_innerk(%esp)
        jl    _nb_kernel131nf_ia32_sse2.nb131nf_checksingle
        jmp   _nb_kernel131nf_ia32_sse2.nb131nf_unroll_loop
_nb_kernel131nf_ia32_sse2.nb131nf_checksingle: 
        movl  nb131nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel131nf_ia32_sse2.nb131nf_dosingle
        jmp   _nb_kernel131nf_ia32_sse2.nb131nf_updateouterdata
_nb_kernel131nf_ia32_sse2.nb131nf_dosingle: 
        movl  nb131nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb131nf_innerjjnr(%esp)

        movl nb131nf_charge(%ebp),%esi     ## base of charge[] 
        xorpd %xmm3,%xmm3
        movlpd (%esi,%eax,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb131nf_iqO(%esp),%xmm3
        mulpd  nb131nf_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 

        movapd  %xmm3,nb131nf_qqO(%esp)
        movapd  %xmm4,nb131nf_qqH(%esp)

        movl nb131nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb131nf_vdwparam(%ebp),%esi
        shll %eax
        movl nb131nf_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 
        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movapd %xmm4,nb131nf_c6(%esp)
        movapd %xmm6,nb131nf_c12(%esp)

        movl nb131nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb131nf_ixO(%esp),%xmm4
        movapd nb131nf_iyO(%esp),%xmm5
        movapd nb131nf_izO(%esp),%xmm6

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
        movapd %xmm7,nb131nf_rsqO(%esp)

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb131nf_ixH1(%esp),%xmm4
        movapd nb131nf_iyH1(%esp),%xmm5
        movapd nb131nf_izH1(%esp),%xmm6

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
        movapd nb131nf_ixH2(%esp),%xmm3
        movapd nb131nf_iyH2(%esp),%xmm4
        movapd nb131nf_izH2(%esp),%xmm5

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

        ## start with rsqO - put seed in xmm2 
        cvtsd2ss %xmm7,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb131nf_three(%esp),%xmm4
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb131nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb131nf_three(%esp),%xmm4
        subsd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb131nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvO in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb131nf_three(%esp),%xmm4
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb131nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb131nf_three(%esp),%xmm4
        subsd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb131nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 
        movsd  %xmm6,nb131nf_rinvH1(%esp)

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb131nf_three(%esp),%xmm4
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb131nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb131nf_three(%esp),%xmm4
        subsd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb131nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 
        movsd %xmm5,nb131nf_rinvH2(%esp)

        ## do O interactions 
        movsd %xmm7,%xmm0
        mulsd  nb131nf_qqO(%esp),%xmm7

        addsd  nb131nf_vctot(%esp),%xmm7
        movsd %xmm7,nb131nf_vctot(%esp)

        movsd nb131nf_rsqO(%esp),%xmm4
        ## LJ table interaction. xmm0=rinv, xmm4=rsq

        mulsd %xmm0,%xmm4       ## xmm4=r 
        mulsd nb131nf_tsc(%esp),%xmm4

        cvttsd2si %xmm4,%ebx    ## mm6 = lu idx 
        cvtsi2sd %ebx,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $3,%ebx

        movl nb131nf_VFtab(%ebp),%esi

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

        movsd nb131nf_c6(%esp),%xmm4
        mulsd  %xmm4,%xmm5       ## Vvdw6 

        ## Update Vvdwtot directly 
        addsd  nb131nf_Vvdwtot(%esp),%xmm5
        movsd %xmm5,nb131nf_Vvdwtot(%esp)

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

        movsd nb131nf_c12(%esp),%xmm4
        mulsd  %xmm4,%xmm5

        addsd  nb131nf_Vvdwtot(%esp),%xmm5
        movsd %xmm5,nb131nf_Vvdwtot(%esp)


        ## H1 & H2 interactions 
        movsd  nb131nf_rinvH1(%esp),%xmm6
        addsd   nb131nf_rinvH2(%esp),%xmm6
        mulsd   nb131nf_qqH(%esp),%xmm6   ## vcoul 
        addsd   nb131nf_vctot(%esp),%xmm6
        movsd  %xmm6,nb131nf_vctot(%esp)

_nb_kernel131nf_ia32_sse2.nb131nf_updateouterdata: 
        ## get n from stack
        movl nb131nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb131nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb131nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb131nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb131nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb131nf_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

       ## finish if last 
        movl nb131nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel131nf_ia32_sse2.nb131nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb131nf_n(%esp)
        jmp _nb_kernel131nf_ia32_sse2.nb131nf_outer
_nb_kernel131nf_ia32_sse2.nb131nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb131nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel131nf_ia32_sse2.nb131nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel131nf_ia32_sse2.nb131nf_threadloop
_nb_kernel131nf_ia32_sse2.nb131nf_end: 
        emms

        movl nb131nf_nouter(%esp),%eax
        movl nb131nf_ninner(%esp),%ebx
        movl nb131nf_outeriter(%ebp),%ecx
        movl nb131nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb131nf_salign(%esp),%eax
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



