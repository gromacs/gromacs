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




.globl nb_kernel211_ia32_sse2
.globl _nb_kernel211_ia32_sse2
nb_kernel211_ia32_sse2: 
_nb_kernel211_ia32_sse2:        
.set nb211_p_nri, 8
.set nb211_iinr, 12
.set nb211_jindex, 16
.set nb211_jjnr, 20
.set nb211_shift, 24
.set nb211_shiftvec, 28
.set nb211_fshift, 32
.set nb211_gid, 36
.set nb211_pos, 40
.set nb211_faction, 44
.set nb211_charge, 48
.set nb211_p_facel, 52
.set nb211_argkrf, 56
.set nb211_argcrf, 60
.set nb211_Vc, 64
.set nb211_type, 68
.set nb211_p_ntype, 72
.set nb211_vdwparam, 76
.set nb211_Vvdw, 80
.set nb211_p_tabscale, 84
.set nb211_VFtab, 88
.set nb211_invsqrta, 92
.set nb211_dvda, 96
.set nb211_p_gbtabscale, 100
.set nb211_GBtab, 104
.set nb211_p_nthreads, 108
.set nb211_count, 112
.set nb211_mtx, 116
.set nb211_outeriter, 120
.set nb211_inneriter, 124
.set nb211_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb211_ixO, 0
.set nb211_iyO, 16
.set nb211_izO, 32
.set nb211_ixH1, 48
.set nb211_iyH1, 64
.set nb211_izH1, 80
.set nb211_ixH2, 96
.set nb211_iyH2, 112
.set nb211_izH2, 128
.set nb211_iqO, 144
.set nb211_iqH, 160
.set nb211_dxO, 176
.set nb211_dyO, 192
.set nb211_dzO, 208
.set nb211_dxH1, 224
.set nb211_dyH1, 240
.set nb211_dzH1, 256
.set nb211_dxH2, 272
.set nb211_dyH2, 288
.set nb211_dzH2, 304
.set nb211_qqO, 320
.set nb211_qqH, 336
.set nb211_c6, 352
.set nb211_c12, 368
.set nb211_six, 384
.set nb211_twelve, 400
.set nb211_vctot, 416
.set nb211_Vvdwtot, 432
.set nb211_fixO, 448
.set nb211_fiyO, 464
.set nb211_fizO, 480
.set nb211_fixH1, 496
.set nb211_fiyH1, 512
.set nb211_fizH1, 528
.set nb211_fixH2, 544
.set nb211_fiyH2, 560
.set nb211_fizH2, 576
.set nb211_fjx, 592
.set nb211_fjy, 608
.set nb211_fjz, 624
.set nb211_half, 640
.set nb211_three, 656
.set nb211_two, 672
.set nb211_krf, 688
.set nb211_crf, 704
.set nb211_krsqO, 720
.set nb211_krsqH1, 736
.set nb211_krsqH2, 752
.set nb211_is3, 768
.set nb211_ii3, 772
.set nb211_ntia, 776
.set nb211_innerjjnr, 780
.set nb211_innerk, 784
.set nb211_n, 788
.set nb211_nn1, 792
.set nb211_nri, 796
.set nb211_nouter, 800
.set nb211_ninner, 804
.set nb211_salign, 808
        pushl %ebp
        movl %esp,%ebp
    pushl %eax
    pushl %ebx
    pushl %ecx
    pushl %edx
        pushl %esi
        pushl %edi
        subl $812,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb211_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb211_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb211_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb211_nouter(%esp)
        movl %eax,nb211_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb211_half(%esp)
        movl %ebx,nb211_half+4(%esp)
        movsd nb211_half(%esp),%xmm1
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
        movapd %xmm1,nb211_half(%esp)
        movapd %xmm2,nb211_two(%esp)
        movapd %xmm3,nb211_three(%esp)
        movapd %xmm4,nb211_six(%esp)
        movapd %xmm5,nb211_twelve(%esp)

        movl nb211_argkrf(%ebp),%esi
        movl nb211_argcrf(%ebp),%edi
        movsd (%esi),%xmm5
        movsd (%edi),%xmm6
        shufpd $0,%xmm5,%xmm5
        shufpd $0,%xmm6,%xmm6
        movapd %xmm5,nb211_krf(%esp)
        movapd %xmm6,nb211_crf(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb211_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb211_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        movsd 8(%edx,%ebx,8),%xmm4
        movl nb211_p_facel(%ebp),%esi
        movsd (%esi),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb211_iqO(%esp)
        movapd %xmm4,nb211_iqH(%esp)

        movl  nb211_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb211_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb211_ntia(%esp)
_nb_kernel211_ia32_sse2.nb211_threadloop: 
        movl  nb211_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel211_ia32_sse2.nb211_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel211_ia32_sse2.nb211_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb211_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb211_n(%esp)
        movl %ebx,nb211_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel211_ia32_sse2.nb211_outerstart
        jmp _nb_kernel211_ia32_sse2.nb211_end

_nb_kernel211_ia32_sse2.nb211_outerstart: 
        ## ebx contains number of outer iterations
        addl nb211_nouter(%esp),%ebx
        movl %ebx,nb211_nouter(%esp)

_nb_kernel211_ia32_sse2.nb211_outer: 
        movl  nb211_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb211_is3(%esp)      ## store is3 

        movl  nb211_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb211_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb211_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb211_ii3(%esp)

        addsd (%eax,%ebx,8),%xmm3
        addsd 8(%eax,%ebx,8),%xmm4
        addsd 16(%eax,%ebx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb211_ixO(%esp)
        movapd %xmm4,nb211_iyO(%esp)
        movapd %xmm5,nb211_izO(%esp)

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
        movapd %xmm0,nb211_ixH1(%esp)
        movapd %xmm1,nb211_iyH1(%esp)
        movapd %xmm2,nb211_izH1(%esp)
        movapd %xmm3,nb211_ixH2(%esp)
        movapd %xmm4,nb211_iyH2(%esp)
        movapd %xmm5,nb211_izH2(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb211_vctot(%esp)
        movapd %xmm4,nb211_Vvdwtot(%esp)
        movapd %xmm4,nb211_fixO(%esp)
        movapd %xmm4,nb211_fiyO(%esp)
        movapd %xmm4,nb211_fizO(%esp)
        movapd %xmm4,nb211_fixH1(%esp)
        movapd %xmm4,nb211_fiyH1(%esp)
        movapd %xmm4,nb211_fizH1(%esp)
        movapd %xmm4,nb211_fixH2(%esp)
        movapd %xmm4,nb211_fiyH2(%esp)
        movapd %xmm4,nb211_fizH2(%esp)

        movl  nb211_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx     ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb211_pos(%ebp),%esi
        movl  nb211_faction(%ebp),%edi
        movl  nb211_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb211_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb211_ninner(%esp),%ecx
        movl  %ecx,nb211_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb211_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel211_ia32_sse2.nb211_unroll_loop
        jmp   _nb_kernel211_ia32_sse2.nb211_checksingle
_nb_kernel211_ia32_sse2.nb211_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb211_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb211_innerjjnr(%esp)                   ## advance pointer (unrolled 2) 

        movl nb211_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb211_iqO(%esp),%xmm3
        mulpd  nb211_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movapd  %xmm3,nb211_qqO(%esp)
        movapd  %xmm4,nb211_qqH(%esp)

        movl nb211_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb211_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb211_ntia(%esp),%edi
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
        movapd %xmm4,nb211_c6(%esp)
        movapd %xmm6,nb211_c12(%esp)

        movl nb211_pos(%ebp),%esi        ## base of pos[] 

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
        movapd nb211_ixO(%esp),%xmm4
        movapd nb211_iyO(%esp),%xmm5
        movapd nb211_izO(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb211_dxO(%esp)
        movapd %xmm5,nb211_dyO(%esp)
        movapd %xmm6,nb211_dzO(%esp)
        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb211_ixH1(%esp),%xmm4
        movapd nb211_iyH1(%esp),%xmm5
        movapd nb211_izH1(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb211_dxH1(%esp)
        movapd %xmm5,nb211_dyH1(%esp)
        movapd %xmm6,nb211_dzH1(%esp)
        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm6
        addpd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb211_ixH2(%esp),%xmm3
        movapd nb211_iyH2(%esp),%xmm4
        movapd nb211_izH2(%esp),%xmm5

        ## calc dr 
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5

        ## store dr 
        movapd %xmm3,nb211_dxH2(%esp)
        movapd %xmm4,nb211_dyH2(%esp)
        movapd %xmm5,nb211_dzH2(%esp)
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

        mulpd  nb211_krf(%esp),%xmm0
        mulpd  nb211_krf(%esp),%xmm1
        mulpd  nb211_krf(%esp),%xmm2

        movapd %xmm0,nb211_krsqH2(%esp)
        movapd %xmm1,nb211_krsqH1(%esp)
        movapd %xmm2,nb211_krsqO(%esp)

        ## start with rsqO - put seed in xmm2 
        cvtpd2ps %xmm7,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb211_three(%esp),%xmm4
        mulpd   %xmm7,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb211_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb211_three(%esp),%xmm4
        subpd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb211_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvO in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb211_three(%esp),%xmm4
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb211_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb211_three(%esp),%xmm4
        subpd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb211_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb211_three(%esp),%xmm4
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb211_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb211_three(%esp),%xmm4
        subpd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb211_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do O interactions 
        movapd  %xmm7,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq 
        movapd %xmm4,%xmm1
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulpd  nb211_c6(%esp),%xmm1
        mulpd  nb211_c12(%esp),%xmm2
        movapd %xmm2,%xmm3
        subpd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addpd  nb211_Vvdwtot(%esp),%xmm3
        mulpd  nb211_six(%esp),%xmm1
        mulpd  nb211_twelve(%esp),%xmm2
        subpd  %xmm1,%xmm2      ## nb part of fs  

        movapd %xmm7,%xmm0
        movapd nb211_krsqO(%esp),%xmm1
        addpd  %xmm1,%xmm0
        mulpd  nb211_two(%esp),%xmm1
        subpd  nb211_crf(%esp),%xmm0   ## xmm0=rinv+ krsq-crf 
        subpd  %xmm1,%xmm7
        mulpd  nb211_qqO(%esp),%xmm0
        mulpd  nb211_qqO(%esp),%xmm7
        addpd  %xmm7,%xmm2

        mulpd  %xmm2,%xmm4      ## total fsO in xmm4 

        addpd  nb211_vctot(%esp),%xmm0
        movapd %xmm3,nb211_Vvdwtot(%esp)
        movapd %xmm0,nb211_vctot(%esp)

        movapd nb211_dxO(%esp),%xmm0
        movapd nb211_dyO(%esp),%xmm1
        movapd nb211_dzO(%esp),%xmm2
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        ## update O forces 
        movapd nb211_fixO(%esp),%xmm3
        movapd nb211_fiyO(%esp),%xmm4
        movapd nb211_fizO(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb211_fixO(%esp)
        movapd %xmm4,nb211_fiyO(%esp)
        movapd %xmm7,nb211_fizO(%esp)
        ## update j forces with water O 
        movapd %xmm0,nb211_fjx(%esp)
        movapd %xmm1,nb211_fjy(%esp)
        movapd %xmm2,nb211_fjz(%esp)

        ## H1 interactions 
        movapd  %xmm6,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        movapd  %xmm6,%xmm7
        movapd  nb211_krsqH1(%esp),%xmm0
        addpd   %xmm0,%xmm6     ## xmm6=rinv+ krsq 
        mulpd   nb211_two(%esp),%xmm0
        subpd   nb211_crf(%esp),%xmm6
        subpd   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulpd   nb211_qqH(%esp),%xmm6   ## vcoul 
        mulpd   nb211_qqH(%esp),%xmm7
        mulpd  %xmm7,%xmm4              ## total fsH1 in xmm4 

        addpd  nb211_vctot(%esp),%xmm6

        movapd nb211_dxH1(%esp),%xmm0
        movapd nb211_dyH1(%esp),%xmm1
        movapd nb211_dzH1(%esp),%xmm2
        movapd %xmm6,nb211_vctot(%esp)
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb211_fixH1(%esp),%xmm3
        movapd nb211_fiyH1(%esp),%xmm4
        movapd nb211_fizH1(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb211_fixH1(%esp)
        movapd %xmm4,nb211_fiyH1(%esp)
        movapd %xmm7,nb211_fizH1(%esp)
        ## update j forces with water H1 
        addpd  nb211_fjx(%esp),%xmm0
        addpd  nb211_fjy(%esp),%xmm1
        addpd  nb211_fjz(%esp),%xmm2
        movapd %xmm0,nb211_fjx(%esp)
        movapd %xmm1,nb211_fjy(%esp)
        movapd %xmm2,nb211_fjz(%esp)

        ## H2 interactions 
        movapd  %xmm5,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        movapd  %xmm5,%xmm7
        movapd  nb211_krsqH2(%esp),%xmm0
        addpd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        mulpd   nb211_two(%esp),%xmm0
        subpd   nb211_crf(%esp),%xmm5
        subpd   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulpd   nb211_qqH(%esp),%xmm5   ## vcoul 
        mulpd   nb211_qqH(%esp),%xmm7
        mulpd  %xmm7,%xmm4              ## total fsH2 in xmm4 

        addpd  nb211_vctot(%esp),%xmm5

        movapd nb211_dxH2(%esp),%xmm0
        movapd nb211_dyH2(%esp),%xmm1
        movapd nb211_dzH2(%esp),%xmm2
        movapd %xmm5,nb211_vctot(%esp)
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        ## update H2 forces 
        movapd nb211_fixH2(%esp),%xmm3
        movapd nb211_fiyH2(%esp),%xmm4
        movapd nb211_fizH2(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb211_fixH2(%esp)
        movapd %xmm4,nb211_fiyH2(%esp)
        movapd %xmm7,nb211_fizH2(%esp)

        movl nb211_faction(%ebp),%edi
        ## update j forces 
        addpd  nb211_fjx(%esp),%xmm0
        addpd  nb211_fjy(%esp),%xmm1
        addpd  nb211_fjz(%esp),%xmm2
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
        subl $2,nb211_innerk(%esp)
        jl    _nb_kernel211_ia32_sse2.nb211_checksingle
        jmp   _nb_kernel211_ia32_sse2.nb211_unroll_loop
_nb_kernel211_ia32_sse2.nb211_checksingle: 
        movl  nb211_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel211_ia32_sse2.nb211_dosingle
        jmp   _nb_kernel211_ia32_sse2.nb211_updateouterdata
_nb_kernel211_ia32_sse2.nb211_dosingle: 
        movl  nb211_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb211_innerjjnr(%esp)

        movl nb211_charge(%ebp),%esi     ## base of charge[] 
        xorpd %xmm3,%xmm3
        movlpd (%esi,%eax,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb211_iqO(%esp),%xmm3
        mulpd  nb211_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 

        movapd  %xmm3,nb211_qqO(%esp)
        movapd  %xmm4,nb211_qqH(%esp)

        movl nb211_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb211_vdwparam(%ebp),%esi
        shll %eax
        movl nb211_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 
        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movapd %xmm4,nb211_c6(%esp)
        movapd %xmm6,nb211_c12(%esp)

        movl nb211_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb211_ixO(%esp),%xmm4
        movapd nb211_iyO(%esp),%xmm5
        movapd nb211_izO(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb211_dxO(%esp)
        movapd %xmm5,nb211_dyO(%esp)
        movapd %xmm6,nb211_dzO(%esp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb211_ixH1(%esp),%xmm4
        movapd nb211_iyH1(%esp),%xmm5
        movapd nb211_izH1(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb211_dxH1(%esp)
        movapd %xmm5,nb211_dyH1(%esp)
        movapd %xmm6,nb211_dzH1(%esp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm6
        addsd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb211_ixH2(%esp),%xmm3
        movapd nb211_iyH2(%esp),%xmm4
        movapd nb211_izH2(%esp),%xmm5

        ## calc dr 
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5

        ## store dr 
        movapd %xmm3,nb211_dxH2(%esp)
        movapd %xmm4,nb211_dyH2(%esp)
        movapd %xmm5,nb211_dzH2(%esp)
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

        mulsd  nb211_krf(%esp),%xmm0
        mulsd  nb211_krf(%esp),%xmm1
        mulsd  nb211_krf(%esp),%xmm2

        movapd %xmm0,nb211_krsqH2(%esp)
        movapd %xmm1,nb211_krsqH1(%esp)
        movapd %xmm2,nb211_krsqO(%esp)

        ## start with rsqO - put seed in xmm2 
        cvtsd2ss %xmm7,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb211_three(%esp),%xmm4
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb211_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb211_three(%esp),%xmm4
        subsd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb211_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvO in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb211_three(%esp),%xmm4
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb211_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb211_three(%esp),%xmm4
        subsd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb211_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb211_three(%esp),%xmm4
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb211_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb211_three(%esp),%xmm4
        subsd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb211_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do O interactions 
        movapd  %xmm7,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq 
        movapd %xmm4,%xmm1
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulsd  nb211_c6(%esp),%xmm1
        mulsd  nb211_c12(%esp),%xmm2
        movapd %xmm2,%xmm3
        subsd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addsd  nb211_Vvdwtot(%esp),%xmm3
        mulsd  nb211_six(%esp),%xmm1
        mulsd  nb211_twelve(%esp),%xmm2
        subsd  %xmm1,%xmm2      ## nb part of fs  

        movapd %xmm7,%xmm0
        movapd nb211_krsqO(%esp),%xmm1
        addsd  %xmm1,%xmm0
        mulsd  nb211_two(%esp),%xmm1
        subsd  nb211_crf(%esp),%xmm0   ## xmm0=rinv+ krsq-crf 
        subsd  %xmm1,%xmm7
        mulsd  nb211_qqO(%esp),%xmm0
        mulsd  nb211_qqO(%esp),%xmm7
        addsd  %xmm7,%xmm2

        mulsd  %xmm2,%xmm4      ## total fsO in xmm4 

        addsd  nb211_vctot(%esp),%xmm0
        movlpd %xmm3,nb211_Vvdwtot(%esp)
        movlpd %xmm0,nb211_vctot(%esp)

        movapd nb211_dxO(%esp),%xmm0
        movapd nb211_dyO(%esp),%xmm1
        movapd nb211_dzO(%esp),%xmm2
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update O forces 
        movapd nb211_fixO(%esp),%xmm3
        movapd nb211_fiyO(%esp),%xmm4
        movapd nb211_fizO(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb211_fixO(%esp)
        movlpd %xmm4,nb211_fiyO(%esp)
        movlpd %xmm7,nb211_fizO(%esp)
        ## update j forces with water O 
        movlpd %xmm0,nb211_fjx(%esp)
        movlpd %xmm1,nb211_fjy(%esp)
        movlpd %xmm2,nb211_fjz(%esp)

        ## H1 interactions 
        movapd  %xmm6,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        movapd  %xmm6,%xmm7
        movapd  nb211_krsqH1(%esp),%xmm0
        addsd   %xmm0,%xmm6     ## xmm6=rinv+ krsq 
        mulsd   nb211_two(%esp),%xmm0
        subsd   nb211_crf(%esp),%xmm6
        subsd   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulsd   nb211_qqH(%esp),%xmm6   ## vcoul 
        mulsd   nb211_qqH(%esp),%xmm7
        mulsd  %xmm7,%xmm4              ## total fsH1 in xmm4 

        addsd  nb211_vctot(%esp),%xmm6

        movapd nb211_dxH1(%esp),%xmm0
        movapd nb211_dyH1(%esp),%xmm1
        movapd nb211_dzH1(%esp),%xmm2
        movlpd %xmm6,nb211_vctot(%esp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb211_fixH1(%esp),%xmm3
        movapd nb211_fiyH1(%esp),%xmm4
        movapd nb211_fizH1(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb211_fixH1(%esp)
        movlpd %xmm4,nb211_fiyH1(%esp)
        movlpd %xmm7,nb211_fizH1(%esp)
        ## update j forces with water H1 
        addsd  nb211_fjx(%esp),%xmm0
        addsd  nb211_fjy(%esp),%xmm1
        addsd  nb211_fjz(%esp),%xmm2
        movlpd %xmm0,nb211_fjx(%esp)
        movlpd %xmm1,nb211_fjy(%esp)
        movlpd %xmm2,nb211_fjz(%esp)

        ## H2 interactions 
        movapd  %xmm5,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        movapd  %xmm5,%xmm7
        movapd  nb211_krsqH2(%esp),%xmm0
        addsd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        mulsd   nb211_two(%esp),%xmm0
        subsd   nb211_crf(%esp),%xmm5
        subsd   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulsd   nb211_qqH(%esp),%xmm5   ## vcoul 
        mulsd   nb211_qqH(%esp),%xmm7
        mulsd  %xmm7,%xmm4              ## total fsH2 in xmm4 

        addsd  nb211_vctot(%esp),%xmm5

        movapd nb211_dxH2(%esp),%xmm0
        movapd nb211_dyH2(%esp),%xmm1
        movapd nb211_dzH2(%esp),%xmm2
        movlpd %xmm5,nb211_vctot(%esp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H2 forces 
        movapd nb211_fixH2(%esp),%xmm3
        movapd nb211_fiyH2(%esp),%xmm4
        movapd nb211_fizH2(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb211_fixH2(%esp)
        movlpd %xmm4,nb211_fiyH2(%esp)
        movlpd %xmm7,nb211_fizH2(%esp)

        movl nb211_faction(%ebp),%edi
        ## update j forces 
        addsd  nb211_fjx(%esp),%xmm0
        addsd  nb211_fjy(%esp),%xmm1
        addsd  nb211_fjz(%esp),%xmm2
        movlpd (%edi,%eax,8),%xmm3
        movlpd 8(%edi,%eax,8),%xmm4
        movlpd 16(%edi,%eax,8),%xmm5
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        movlpd %xmm3,(%edi,%eax,8)
        movlpd %xmm4,8(%edi,%eax,8)
        movlpd %xmm5,16(%edi,%eax,8)

_nb_kernel211_ia32_sse2.nb211_updateouterdata: 
        movl  nb211_ii3(%esp),%ecx
        movl  nb211_faction(%ebp),%edi
        movl  nb211_fshift(%ebp),%esi
        movl  nb211_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movapd nb211_fixO(%esp),%xmm0
        movapd nb211_fiyO(%esp),%xmm1
        movapd nb211_fizO(%esp),%xmm2

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
        movapd nb211_fixH1(%esp),%xmm0
        movapd nb211_fiyH1(%esp),%xmm1
        movapd nb211_fizH1(%esp),%xmm2

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
        movapd nb211_fixH2(%esp),%xmm0
        movapd nb211_fiyH2(%esp),%xmm1
        movapd nb211_fizH2(%esp),%xmm2

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
        movl nb211_n(%esp),%esi
        ## get group index for i particle 
        movl  nb211_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb211_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb211_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb211_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb211_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

       ## finish if last 
        movl nb211_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel211_ia32_sse2.nb211_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb211_n(%esp)
        jmp _nb_kernel211_ia32_sse2.nb211_outer
_nb_kernel211_ia32_sse2.nb211_outerend: 
        ## check if more outer neighborlists remain
        movl  nb211_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel211_ia32_sse2.nb211_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel211_ia32_sse2.nb211_threadloop
_nb_kernel211_ia32_sse2.nb211_end: 
        emms

        movl nb211_nouter(%esp),%eax
        movl nb211_ninner(%esp),%ebx
        movl nb211_outeriter(%ebp),%ecx
        movl nb211_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb211_salign(%esp),%eax
        addl %eax,%esp
        addl $812,%esp
        popl %edi
        popl %esi
    popl %edx
    popl %ecx
    popl %ebx
    popl %eax
        leave
        ret




.globl nb_kernel211nf_ia32_sse2
.globl _nb_kernel211nf_ia32_sse2
nb_kernel211nf_ia32_sse2:       
_nb_kernel211nf_ia32_sse2:      
.set nb211nf_p_nri, 8
.set nb211nf_iinr, 12
.set nb211nf_jindex, 16
.set nb211nf_jjnr, 20
.set nb211nf_shift, 24
.set nb211nf_shiftvec, 28
.set nb211nf_fshift, 32
.set nb211nf_gid, 36
.set nb211nf_pos, 40
.set nb211nf_faction, 44
.set nb211nf_charge, 48
.set nb211nf_p_facel, 52
.set nb211nf_argkrf, 56
.set nb211nf_argcrf, 60
.set nb211nf_Vc, 64
.set nb211nf_type, 68
.set nb211nf_p_ntype, 72
.set nb211nf_vdwparam, 76
.set nb211nf_Vvdw, 80
.set nb211nf_p_tabscale, 84
.set nb211nf_VFtab, 88
.set nb211nf_invsqrta, 92
.set nb211nf_dvda, 96
.set nb211nf_p_gbtabscale, 100
.set nb211nf_GBtab, 104
.set nb211nf_p_nthreads, 108
.set nb211nf_count, 112
.set nb211nf_mtx, 116
.set nb211nf_outeriter, 120
.set nb211nf_inneriter, 124
.set nb211nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb211nf_ixO, 0
.set nb211nf_iyO, 16
.set nb211nf_izO, 32
.set nb211nf_ixH1, 48
.set nb211nf_iyH1, 64
.set nb211nf_izH1, 80
.set nb211nf_ixH2, 96
.set nb211nf_iyH2, 112
.set nb211nf_izH2, 128
.set nb211nf_iqO, 144
.set nb211nf_iqH, 160
.set nb211nf_qqO, 176
.set nb211nf_qqH, 192
.set nb211nf_c6, 208
.set nb211nf_c12, 224
.set nb211nf_vctot, 240
.set nb211nf_Vvdwtot, 256
.set nb211nf_half, 272
.set nb211nf_three, 288
.set nb211nf_krf, 304
.set nb211nf_crf, 320
.set nb211nf_krsqO, 336
.set nb211nf_krsqH1, 352
.set nb211nf_krsqH2, 368
.set nb211nf_is3, 384
.set nb211nf_ii3, 388
.set nb211nf_ntia, 392
.set nb211nf_innerjjnr, 396
.set nb211nf_innerk, 400
.set nb211nf_n, 404
.set nb211nf_nn1, 408
.set nb211nf_nri, 412
.set nb211nf_nouter, 416
.set nb211nf_ninner, 420
.set nb211nf_salign, 424
        pushl %ebp
        movl %esp,%ebp
    pushl %eax
    pushl %ebx
    pushl %ecx
    pushl %edx
        pushl %esi
        pushl %edi
        subl $428,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb211nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb211nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb211nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb211nf_nouter(%esp)
        movl %eax,nb211nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb211nf_half(%esp)
        movl %ebx,nb211nf_half+4(%esp)
        movsd nb211nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb211nf_half(%esp)
        movapd %xmm3,nb211nf_three(%esp)

        movl nb211nf_argkrf(%ebp),%esi
        movl nb211nf_argcrf(%ebp),%edi
        movsd (%esi),%xmm5
        movsd (%edi),%xmm6
        shufpd $0,%xmm5,%xmm5
        shufpd $0,%xmm6,%xmm6
        movapd %xmm5,nb211nf_krf(%esp)
        movapd %xmm6,nb211nf_crf(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb211nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb211nf_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        movsd 8(%edx,%ebx,8),%xmm4
        movl nb211nf_p_facel(%ebp),%esi
        movsd (%esi),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb211nf_iqO(%esp)
        movapd %xmm4,nb211nf_iqH(%esp)

        movl  nb211nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb211nf_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb211nf_ntia(%esp)
_nb_kernel211nf_ia32_sse2.nb211nf_threadloop: 
        movl  nb211nf_count(%ebp),%esi          ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel211nf_ia32_sse2.nb211nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel211nf_ia32_sse2.nb211nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb211nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb211nf_n(%esp)
        movl %ebx,nb211nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel211nf_ia32_sse2.nb211nf_outerstart
        jmp _nb_kernel211nf_ia32_sse2.nb211nf_end

_nb_kernel211nf_ia32_sse2.nb211nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb211nf_nouter(%esp),%ebx
        movl %ebx,nb211nf_nouter(%esp)

_nb_kernel211nf_ia32_sse2.nb211nf_outer: 
        movl  nb211nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb211nf_is3(%esp)            ## store is3 

        movl  nb211nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb211nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb211nf_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb211nf_ii3(%esp)

        addsd (%eax,%ebx,8),%xmm3
        addsd 8(%eax,%ebx,8),%xmm4
        addsd 16(%eax,%ebx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb211nf_ixO(%esp)
        movapd %xmm4,nb211nf_iyO(%esp)
        movapd %xmm5,nb211nf_izO(%esp)

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
        movapd %xmm0,nb211nf_ixH1(%esp)
        movapd %xmm1,nb211nf_iyH1(%esp)
        movapd %xmm2,nb211nf_izH1(%esp)
        movapd %xmm3,nb211nf_ixH2(%esp)
        movapd %xmm4,nb211nf_iyH2(%esp)
        movapd %xmm5,nb211nf_izH2(%esp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb211nf_vctot(%esp)
        movapd %xmm4,nb211nf_Vvdwtot(%esp)

        movl  nb211nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb211nf_pos(%ebp),%esi
        movl  nb211nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb211nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb211nf_ninner(%esp),%ecx
        movl  %ecx,nb211nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb211nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel211nf_ia32_sse2.nb211nf_unroll_loop
        jmp   _nb_kernel211nf_ia32_sse2.nb211nf_checksingle
_nb_kernel211nf_ia32_sse2.nb211nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb211nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb211nf_innerjjnr(%esp)                 ## advance pointer (unrolled 2) 

        movl nb211nf_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb211nf_iqO(%esp),%xmm3
        mulpd  nb211nf_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movapd  %xmm3,nb211nf_qqO(%esp)
        movapd  %xmm4,nb211nf_qqH(%esp)

        movl nb211nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb211nf_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb211nf_ntia(%esp),%edi
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
        movapd %xmm4,nb211nf_c6(%esp)
        movapd %xmm6,nb211nf_c12(%esp)

        movl nb211nf_pos(%ebp),%esi        ## base of pos[] 

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
        movapd nb211nf_ixO(%esp),%xmm4
        movapd nb211nf_iyO(%esp),%xmm5
        movapd nb211nf_izO(%esp),%xmm6

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
        movapd nb211nf_ixH1(%esp),%xmm4
        movapd nb211nf_iyH1(%esp),%xmm5
        movapd nb211nf_izH1(%esp),%xmm6

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
        movapd nb211nf_ixH2(%esp),%xmm3
        movapd nb211nf_iyH2(%esp),%xmm4
        movapd nb211nf_izH2(%esp),%xmm5

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

        mulpd  nb211nf_krf(%esp),%xmm0
        mulpd  nb211nf_krf(%esp),%xmm1
        mulpd  nb211nf_krf(%esp),%xmm2

        movapd %xmm0,nb211nf_krsqH2(%esp)
        movapd %xmm1,nb211nf_krsqH1(%esp)
        movapd %xmm2,nb211nf_krsqO(%esp)

        ## start with rsqO - put seed in xmm2 
        cvtpd2ps %xmm7,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb211nf_three(%esp),%xmm4
        mulpd   %xmm7,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb211nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb211nf_three(%esp),%xmm4
        subpd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb211nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvO in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb211nf_three(%esp),%xmm4
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb211nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb211nf_three(%esp),%xmm4
        subpd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb211nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb211nf_three(%esp),%xmm4
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb211nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb211nf_three(%esp),%xmm4
        subpd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb211nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do O interactions 
        movapd  %xmm7,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq 
        movapd %xmm4,%xmm1
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulpd  nb211nf_c6(%esp),%xmm1
        mulpd  nb211nf_c12(%esp),%xmm2
        movapd %xmm2,%xmm3
        subpd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addpd  nb211nf_Vvdwtot(%esp),%xmm3

        movapd %xmm7,%xmm0
        movapd nb211nf_krsqO(%esp),%xmm1
        addpd  %xmm1,%xmm0
        subpd  nb211nf_crf(%esp),%xmm0   ## xmm0=rinv+ krsq-crf 
        mulpd  nb211nf_qqO(%esp),%xmm0

        addpd  nb211nf_vctot(%esp),%xmm0
        movapd %xmm3,nb211nf_Vvdwtot(%esp)
        movapd %xmm0,nb211nf_vctot(%esp)

        ## H1 interactions 
        movapd  nb211nf_krsqH1(%esp),%xmm0
        addpd   %xmm0,%xmm6     ## xmm6=rinv+ krsq 
        subpd   nb211nf_crf(%esp),%xmm6
        mulpd   nb211nf_qqH(%esp),%xmm6   ## vcoul      
        addpd  nb211nf_vctot(%esp),%xmm6
        movapd %xmm6,nb211nf_vctot(%esp)

        ## H2 interactions 
        movapd  nb211nf_krsqH2(%esp),%xmm0
        addpd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        subpd   nb211nf_crf(%esp),%xmm5
        mulpd   nb211nf_qqH(%esp),%xmm5   ## vcoul 
        addpd  nb211nf_vctot(%esp),%xmm5
        movapd %xmm5,nb211nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $2,nb211nf_innerk(%esp)
        jl    _nb_kernel211nf_ia32_sse2.nb211nf_checksingle
        jmp   _nb_kernel211nf_ia32_sse2.nb211nf_unroll_loop
_nb_kernel211nf_ia32_sse2.nb211nf_checksingle: 
        movl  nb211nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel211nf_ia32_sse2.nb211nf_dosingle
        jmp   _nb_kernel211nf_ia32_sse2.nb211nf_updateouterdata
_nb_kernel211nf_ia32_sse2.nb211nf_dosingle: 
        movl  nb211nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb211nf_innerjjnr(%esp)

        movl nb211nf_charge(%ebp),%esi     ## base of charge[] 
        xorpd %xmm3,%xmm3
        movlpd (%esi,%eax,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb211nf_iqO(%esp),%xmm3
        mulpd  nb211nf_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 

        movapd  %xmm3,nb211nf_qqO(%esp)
        movapd  %xmm4,nb211nf_qqH(%esp)

        movl nb211nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb211nf_vdwparam(%ebp),%esi
        shll %eax
        movl nb211nf_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 
        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movapd %xmm4,nb211nf_c6(%esp)
        movapd %xmm6,nb211nf_c12(%esp)

        movl nb211nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb211nf_ixO(%esp),%xmm4
        movapd nb211nf_iyO(%esp),%xmm5
        movapd nb211nf_izO(%esp),%xmm6

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
        movapd nb211nf_ixH1(%esp),%xmm4
        movapd nb211nf_iyH1(%esp),%xmm5
        movapd nb211nf_izH1(%esp),%xmm6

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
        movapd nb211nf_ixH2(%esp),%xmm3
        movapd nb211nf_iyH2(%esp),%xmm4
        movapd nb211nf_izH2(%esp),%xmm5

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

        mulsd  nb211nf_krf(%esp),%xmm0
        mulsd  nb211nf_krf(%esp),%xmm1
        mulsd  nb211nf_krf(%esp),%xmm2

        movapd %xmm0,nb211nf_krsqH2(%esp)
        movapd %xmm1,nb211nf_krsqH1(%esp)
        movapd %xmm2,nb211nf_krsqO(%esp)

        ## start with rsqO - put seed in xmm2 
        cvtsd2ss %xmm7,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb211nf_three(%esp),%xmm4
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb211nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb211nf_three(%esp),%xmm4
        subsd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb211nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvO in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb211nf_three(%esp),%xmm4
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb211nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb211nf_three(%esp),%xmm4
        subsd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb211nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb211nf_three(%esp),%xmm4
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb211nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb211nf_three(%esp),%xmm4
        subsd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb211nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do O interactions 
        movapd  %xmm7,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq 
        movapd %xmm4,%xmm1
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulsd  nb211nf_c6(%esp),%xmm1
        mulsd  nb211nf_c12(%esp),%xmm2
        movapd %xmm2,%xmm3
        subsd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addsd  nb211nf_Vvdwtot(%esp),%xmm3

        movapd %xmm7,%xmm0
        movapd nb211nf_krsqO(%esp),%xmm1
        addsd  %xmm1,%xmm0
        subsd  nb211nf_crf(%esp),%xmm0   ## xmm0=rinv+ krsq-crf 
        mulsd  nb211nf_qqO(%esp),%xmm0
        addsd  nb211nf_vctot(%esp),%xmm0
        movlpd %xmm3,nb211nf_Vvdwtot(%esp)
        movlpd %xmm0,nb211nf_vctot(%esp)

        ## H1 interactions 
        movapd  nb211nf_krsqH1(%esp),%xmm0
        addsd   %xmm0,%xmm6     ## xmm6=rinv+ krsq 
        subsd   nb211nf_crf(%esp),%xmm6
        mulsd   nb211nf_qqH(%esp),%xmm6   ## vcoul 
        addsd  nb211nf_vctot(%esp),%xmm6
        movlpd %xmm6,nb211nf_vctot(%esp)

        ## H2 interactions 
        movapd  nb211nf_krsqH2(%esp),%xmm0
        addsd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        subsd   nb211nf_crf(%esp),%xmm5
        mulsd   nb211nf_qqH(%esp),%xmm5   ## vcoul 
        addsd  nb211nf_vctot(%esp),%xmm5
        movlpd %xmm5,nb211nf_vctot(%esp)

_nb_kernel211nf_ia32_sse2.nb211nf_updateouterdata: 
        ## get n from stack
        movl nb211nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb211nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb211nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb211nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb211nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb211nf_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb211nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel211nf_ia32_sse2.nb211nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb211nf_n(%esp)
        jmp _nb_kernel211nf_ia32_sse2.nb211nf_outer
_nb_kernel211nf_ia32_sse2.nb211nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb211nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel211nf_ia32_sse2.nb211nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel211nf_ia32_sse2.nb211nf_threadloop
_nb_kernel211nf_ia32_sse2.nb211nf_end: 
        emms

        movl nb211nf_nouter(%esp),%eax
        movl nb211nf_ninner(%esp),%ebx
        movl nb211nf_outeriter(%ebp),%ecx
        movl nb211nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb211nf_salign(%esp),%eax
        addl %eax,%esp
        addl $428,%esp
        popl %edi
        popl %esi
    popl %edx
    popl %ecx
    popl %ebx
    popl %eax
        leave
        ret

