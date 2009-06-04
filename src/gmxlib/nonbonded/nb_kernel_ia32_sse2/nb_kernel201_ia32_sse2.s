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



.globl nb_kernel201_ia32_sse2
.globl _nb_kernel201_ia32_sse2
nb_kernel201_ia32_sse2: 
_nb_kernel201_ia32_sse2:        
.set nb201_p_nri, 8
.set nb201_iinr, 12
.set nb201_jindex, 16
.set nb201_jjnr, 20
.set nb201_shift, 24
.set nb201_shiftvec, 28
.set nb201_fshift, 32
.set nb201_gid, 36
.set nb201_pos, 40
.set nb201_faction, 44
.set nb201_charge, 48
.set nb201_p_facel, 52
.set nb201_argkrf, 56
.set nb201_argcrf, 60
.set nb201_Vc, 64
.set nb201_type, 68
.set nb201_p_ntype, 72
.set nb201_vdwparam, 76
.set nb201_Vvdw, 80
.set nb201_p_tabscale, 84
.set nb201_VFtab, 88
.set nb201_invsqrta, 92
.set nb201_dvda, 96
.set nb201_p_gbtabscale, 100
.set nb201_GBtab, 104
.set nb201_p_nthreads, 108
.set nb201_count, 112
.set nb201_mtx, 116
.set nb201_outeriter, 120
.set nb201_inneriter, 124
.set nb201_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb201_ixO, 0
.set nb201_iyO, 16
.set nb201_izO, 32
.set nb201_ixH1, 48
.set nb201_iyH1, 64
.set nb201_izH1, 80
.set nb201_ixH2, 96
.set nb201_iyH2, 112
.set nb201_izH2, 128
.set nb201_iqO, 144
.set nb201_iqH, 160
.set nb201_dxO, 176
.set nb201_dyO, 192
.set nb201_dzO, 208
.set nb201_dxH1, 224
.set nb201_dyH1, 240
.set nb201_dzH1, 256
.set nb201_dxH2, 272
.set nb201_dyH2, 288
.set nb201_dzH2, 304
.set nb201_qqO, 320
.set nb201_qqH, 336
.set nb201_vctot, 352
.set nb201_fixO, 384
.set nb201_fiyO, 400
.set nb201_fizO, 416
.set nb201_fixH1, 432
.set nb201_fiyH1, 448
.set nb201_fizH1, 464
.set nb201_fixH2, 480
.set nb201_fiyH2, 496
.set nb201_fizH2, 512
.set nb201_fjx, 528
.set nb201_fjy, 544
.set nb201_fjz, 560
.set nb201_half, 576
.set nb201_three, 592
.set nb201_two, 608
.set nb201_krf, 624
.set nb201_crf, 640
.set nb201_krsqO, 656
.set nb201_krsqH1, 672
.set nb201_krsqH2, 688
.set nb201_is3, 704
.set nb201_ii3, 708
.set nb201_innerjjnr, 712
.set nb201_innerk, 716
.set nb201_n, 720
.set nb201_nn1, 724
.set nb201_nri, 728
.set nb201_nouter, 732
.set nb201_ninner, 736
.set nb201_salign, 740
        pushl %ebp
        movl %esp,%ebp
    pushl %eax
    pushl %ebx
    pushl %ecx
    pushl %edx
        pushl %esi
        pushl %edi
        subl $744,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb201_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb201_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb201_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb201_nouter(%esp)
        movl %eax,nb201_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb201_half(%esp)
        movl %ebx,nb201_half+4(%esp)
        movsd nb201_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb201_half(%esp)
        movapd %xmm2,nb201_two(%esp)
        movapd %xmm3,nb201_three(%esp)

        movl nb201_argkrf(%ebp),%esi
        movl nb201_argcrf(%ebp),%edi
        movsd (%esi),%xmm5
        movsd (%edi),%xmm6
        shufpd $0,%xmm5,%xmm5
        shufpd $0,%xmm6,%xmm6
        movapd %xmm5,nb201_krf(%esp)
        movapd %xmm6,nb201_crf(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb201_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb201_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        movsd 8(%edx,%ebx,8),%xmm4
        movl nb201_p_facel(%ebp),%esi
        movsd (%esi),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb201_iqO(%esp)
        movapd %xmm4,nb201_iqH(%esp)

_nb_kernel201_ia32_sse2.nb201_threadloop: 
        movl  nb201_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel201_ia32_sse2.nb201_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel201_ia32_sse2.nb201_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb201_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb201_n(%esp)
        movl %ebx,nb201_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel201_ia32_sse2.nb201_outerstart
        jmp _nb_kernel201_ia32_sse2.nb201_end

_nb_kernel201_ia32_sse2.nb201_outerstart: 
        ## ebx contains number of outer iterations
        addl nb201_nouter(%esp),%ebx
        movl %ebx,nb201_nouter(%esp)

_nb_kernel201_ia32_sse2.nb201_outer: 
        movl  nb201_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb201_is3(%esp)      ## store is3 

        movl  nb201_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb201_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb201_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb201_ii3(%esp)

        addsd (%eax,%ebx,8),%xmm3
        addsd 8(%eax,%ebx,8),%xmm4
        addsd 16(%eax,%ebx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb201_ixO(%esp)
        movapd %xmm4,nb201_iyO(%esp)
        movapd %xmm5,nb201_izO(%esp)

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
        movapd %xmm0,nb201_ixH1(%esp)
        movapd %xmm1,nb201_iyH1(%esp)
        movapd %xmm2,nb201_izH1(%esp)
        movapd %xmm3,nb201_ixH2(%esp)
        movapd %xmm4,nb201_iyH2(%esp)
        movapd %xmm5,nb201_izH2(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb201_vctot(%esp)
        movapd %xmm4,nb201_fixO(%esp)
        movapd %xmm4,nb201_fiyO(%esp)
        movapd %xmm4,nb201_fizO(%esp)
        movapd %xmm4,nb201_fixH1(%esp)
        movapd %xmm4,nb201_fiyH1(%esp)
        movapd %xmm4,nb201_fizH1(%esp)
        movapd %xmm4,nb201_fixH2(%esp)
        movapd %xmm4,nb201_fiyH2(%esp)
        movapd %xmm4,nb201_fizH2(%esp)

        movl  nb201_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb201_pos(%ebp),%esi
        movl  nb201_faction(%ebp),%edi
        movl  nb201_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb201_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb201_ninner(%esp),%ecx
        movl  %ecx,nb201_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb201_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel201_ia32_sse2.nb201_unroll_loop
        jmp   _nb_kernel201_ia32_sse2.nb201_checksingle
_nb_kernel201_ia32_sse2.nb201_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb201_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb201_innerjjnr(%esp)                   ## advance pointer (unrolled 2) 

        movl nb201_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb201_iqO(%esp),%xmm3
        mulpd  nb201_iqH(%esp),%xmm4
        movapd  %xmm3,nb201_qqO(%esp)
        movapd  %xmm4,nb201_qqH(%esp)

        movl nb201_pos(%ebp),%esi        ## base of pos[] 

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
        movapd nb201_ixO(%esp),%xmm4
        movapd nb201_iyO(%esp),%xmm5
        movapd nb201_izO(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb201_dxO(%esp)
        movapd %xmm5,nb201_dyO(%esp)
        movapd %xmm6,nb201_dzO(%esp)

        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb201_ixH1(%esp),%xmm4
        movapd nb201_iyH1(%esp),%xmm5
        movapd nb201_izH1(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb201_dxH1(%esp)
        movapd %xmm5,nb201_dyH1(%esp)
        movapd %xmm6,nb201_dzH1(%esp)
        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm6
        addpd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb201_ixH2(%esp),%xmm3
        movapd nb201_iyH2(%esp),%xmm4
        movapd nb201_izH2(%esp),%xmm5

        ## calc dr 
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5

        ## store dr 
        movapd %xmm3,nb201_dxH2(%esp)
        movapd %xmm4,nb201_dyH2(%esp)
        movapd %xmm5,nb201_dzH2(%esp)
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

        mulpd  nb201_krf(%esp),%xmm0
        mulpd  nb201_krf(%esp),%xmm1
        mulpd  nb201_krf(%esp),%xmm2

        movapd %xmm0,nb201_krsqH2(%esp)
        movapd %xmm1,nb201_krsqH1(%esp)
        movapd %xmm2,nb201_krsqO(%esp)

        ## start with rsqO - put seed in xmm2 
        cvtpd2ps %xmm7,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb201_three(%esp),%xmm4
        mulpd   %xmm7,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb201_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb201_three(%esp),%xmm4
        subpd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb201_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvO in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb201_three(%esp),%xmm4
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb201_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb201_three(%esp),%xmm4
        subpd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb201_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb201_three(%esp),%xmm4
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb201_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb201_three(%esp),%xmm4
        subpd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb201_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do O interactions 
        movapd  %xmm7,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        movapd  %xmm7,%xmm3
        movapd  nb201_krsqO(%esp),%xmm0
        addpd   %xmm0,%xmm7     ## xmm6=rinv+ krsq 
        mulpd   nb201_two(%esp),%xmm0
        subpd   nb201_crf(%esp),%xmm7
        subpd   %xmm0,%xmm3     ## xmm7=rinv-2*krsq 
        mulpd   nb201_qqO(%esp),%xmm7   ## vcoul 
        mulpd   nb201_qqO(%esp),%xmm3
        mulpd  %xmm3,%xmm4      ## total fsH1 in xmm4 

        addpd  nb201_vctot(%esp),%xmm7

        movapd nb201_dxO(%esp),%xmm0
        movapd nb201_dyO(%esp),%xmm1
        movapd nb201_dzO(%esp),%xmm2
        movapd %xmm7,nb201_vctot(%esp)
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        ## update O forces 
        movapd nb201_fixO(%esp),%xmm3
        movapd nb201_fiyO(%esp),%xmm4
        movapd nb201_fizO(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb201_fixO(%esp)
        movapd %xmm4,nb201_fiyO(%esp)
        movapd %xmm7,nb201_fizO(%esp)
        ## update j forces with water O 
        movapd %xmm0,nb201_fjx(%esp)
        movapd %xmm1,nb201_fjy(%esp)
        movapd %xmm2,nb201_fjz(%esp)

        ## H1 interactions 
        movapd  %xmm6,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        movapd  %xmm6,%xmm7
        movapd  nb201_krsqH1(%esp),%xmm0
        addpd   %xmm0,%xmm6     ## xmm6=rinv+ krsq 
        mulpd   nb201_two(%esp),%xmm0
        subpd   nb201_crf(%esp),%xmm6
        subpd   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulpd   nb201_qqH(%esp),%xmm6   ## vcoul 
        mulpd   nb201_qqH(%esp),%xmm7
        mulpd  %xmm7,%xmm4              ## total fsH1 in xmm4 

        addpd  nb201_vctot(%esp),%xmm6

        movapd nb201_dxH1(%esp),%xmm0
        movapd nb201_dyH1(%esp),%xmm1
        movapd nb201_dzH1(%esp),%xmm2
        movapd %xmm6,nb201_vctot(%esp)
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb201_fixH1(%esp),%xmm3
        movapd nb201_fiyH1(%esp),%xmm4
        movapd nb201_fizH1(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb201_fixH1(%esp)
        movapd %xmm4,nb201_fiyH1(%esp)
        movapd %xmm7,nb201_fizH1(%esp)
        ## update j forces with water H1 
        addpd  nb201_fjx(%esp),%xmm0
        addpd  nb201_fjy(%esp),%xmm1
        addpd  nb201_fjz(%esp),%xmm2
        movapd %xmm0,nb201_fjx(%esp)
        movapd %xmm1,nb201_fjy(%esp)
        movapd %xmm2,nb201_fjz(%esp)

        ## H2 interactions 
        movapd  %xmm5,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        movapd  %xmm5,%xmm7
        movapd  nb201_krsqH2(%esp),%xmm0
        addpd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        mulpd   nb201_two(%esp),%xmm0
        subpd   nb201_crf(%esp),%xmm5
        subpd   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulpd   nb201_qqH(%esp),%xmm5   ## vcoul 
        mulpd   nb201_qqH(%esp),%xmm7
        mulpd  %xmm7,%xmm4              ## total fsH2 in xmm4 

        addpd  nb201_vctot(%esp),%xmm5

        movapd nb201_dxH2(%esp),%xmm0
        movapd nb201_dyH2(%esp),%xmm1
        movapd nb201_dzH2(%esp),%xmm2
        movapd %xmm5,nb201_vctot(%esp)
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        ## update H2 forces 
        movapd nb201_fixH2(%esp),%xmm3
        movapd nb201_fiyH2(%esp),%xmm4
        movapd nb201_fizH2(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb201_fixH2(%esp)
        movapd %xmm4,nb201_fiyH2(%esp)
        movapd %xmm7,nb201_fizH2(%esp)

        movl nb201_faction(%ebp),%edi
        ## update j forces 
        addpd  nb201_fjx(%esp),%xmm0
        addpd  nb201_fjy(%esp),%xmm1
        addpd  nb201_fjz(%esp),%xmm2
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
        subl $2,nb201_innerk(%esp)
        jl    _nb_kernel201_ia32_sse2.nb201_checksingle
        jmp   _nb_kernel201_ia32_sse2.nb201_unroll_loop
_nb_kernel201_ia32_sse2.nb201_checksingle: 
        movl  nb201_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel201_ia32_sse2.nb201_dosingle
        jmp   _nb_kernel201_ia32_sse2.nb201_updateouterdata
_nb_kernel201_ia32_sse2.nb201_dosingle: 
        movl  nb201_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb201_innerjjnr(%esp)

        movl nb201_charge(%ebp),%esi     ## base of charge[] 
        xorpd %xmm3,%xmm3
        movlpd (%esi,%eax,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb201_iqO(%esp),%xmm3
        mulpd  nb201_iqH(%esp),%xmm4
        movapd  %xmm3,nb201_qqO(%esp)
        movapd  %xmm4,nb201_qqH(%esp)

        movl nb201_pos(%ebp),%esi        ## base of pos[] 
        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb201_ixO(%esp),%xmm4
        movapd nb201_iyO(%esp),%xmm5
        movapd nb201_izO(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb201_dxO(%esp)
        movapd %xmm5,nb201_dyO(%esp)
        movapd %xmm6,nb201_dzO(%esp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb201_ixH1(%esp),%xmm4
        movapd nb201_iyH1(%esp),%xmm5
        movapd nb201_izH1(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb201_dxH1(%esp)
        movapd %xmm5,nb201_dyH1(%esp)
        movapd %xmm6,nb201_dzH1(%esp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm6
        addsd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb201_ixH2(%esp),%xmm3
        movapd nb201_iyH2(%esp),%xmm4
        movapd nb201_izH2(%esp),%xmm5

        ## calc dr 
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5

        ## store dr 
        movapd %xmm3,nb201_dxH2(%esp)
        movapd %xmm4,nb201_dyH2(%esp)
        movapd %xmm5,nb201_dzH2(%esp)
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

        mulsd  nb201_krf(%esp),%xmm0
        mulsd  nb201_krf(%esp),%xmm1
        mulsd  nb201_krf(%esp),%xmm2

        movapd %xmm0,nb201_krsqH2(%esp)
        movapd %xmm1,nb201_krsqH1(%esp)
        movapd %xmm2,nb201_krsqO(%esp)

        ## start with rsqO - put seed in xmm2 
        cvtsd2ss %xmm7,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb201_three(%esp),%xmm4
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb201_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb201_three(%esp),%xmm4
        subsd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb201_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvO in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb201_three(%esp),%xmm4
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb201_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb201_three(%esp),%xmm4
        subsd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb201_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb201_three(%esp),%xmm4
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb201_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb201_three(%esp),%xmm4
        subsd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb201_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do O interactions 
        movapd  %xmm7,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        movapd  %xmm7,%xmm3
        movapd  nb201_krsqO(%esp),%xmm0
        addsd   %xmm0,%xmm7     ## xmm6=rinv+ krsq 
        mulsd   nb201_two(%esp),%xmm0
        subsd   nb201_crf(%esp),%xmm7
        subsd   %xmm0,%xmm3     ## xmm7=rinv-2*krsq 
        mulsd   nb201_qqO(%esp),%xmm7   ## vcoul 
        mulsd   nb201_qqO(%esp),%xmm3
        mulsd  %xmm3,%xmm4      ## total fsH1 in xmm4 

        addsd  nb201_vctot(%esp),%xmm7

        movapd nb201_dxO(%esp),%xmm0
        movapd nb201_dyO(%esp),%xmm1
        movapd nb201_dzO(%esp),%xmm2
        movlpd %xmm7,nb201_vctot(%esp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update O forces 
        movapd nb201_fixO(%esp),%xmm3
        movapd nb201_fiyO(%esp),%xmm4
        movapd nb201_fizO(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb201_fixO(%esp)
        movlpd %xmm4,nb201_fiyO(%esp)
        movlpd %xmm7,nb201_fizO(%esp)
        ## update j forces with water O 
        movlpd %xmm0,nb201_fjx(%esp)
        movlpd %xmm1,nb201_fjy(%esp)
        movlpd %xmm2,nb201_fjz(%esp)

        ## H1 interactions 
        movapd  %xmm6,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        movapd  %xmm6,%xmm7
        movapd  nb201_krsqH1(%esp),%xmm0
        addsd   %xmm0,%xmm6     ## xmm6=rinv+ krsq 
        mulsd   nb201_two(%esp),%xmm0
        subsd   nb201_crf(%esp),%xmm6
        subsd   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulsd   nb201_qqH(%esp),%xmm6   ## vcoul 
        mulsd   nb201_qqH(%esp),%xmm7
        mulsd  %xmm7,%xmm4              ## total fsH1 in xmm4 

        addsd  nb201_vctot(%esp),%xmm6

        movapd nb201_dxH1(%esp),%xmm0
        movapd nb201_dyH1(%esp),%xmm1
        movapd nb201_dzH1(%esp),%xmm2
        movlpd %xmm6,nb201_vctot(%esp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb201_fixH1(%esp),%xmm3
        movapd nb201_fiyH1(%esp),%xmm4
        movapd nb201_fizH1(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb201_fixH1(%esp)
        movlpd %xmm4,nb201_fiyH1(%esp)
        movlpd %xmm7,nb201_fizH1(%esp)
        ## update j forces with water H1 
        addsd  nb201_fjx(%esp),%xmm0
        addsd  nb201_fjy(%esp),%xmm1
        addsd  nb201_fjz(%esp),%xmm2
        movlpd %xmm0,nb201_fjx(%esp)
        movlpd %xmm1,nb201_fjy(%esp)
        movlpd %xmm2,nb201_fjz(%esp)

        ## H2 interactions 
        movapd  %xmm5,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        movapd  %xmm5,%xmm7
        movapd  nb201_krsqH2(%esp),%xmm0
        addsd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        mulsd   nb201_two(%esp),%xmm0
        subsd   nb201_crf(%esp),%xmm5
        subsd   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulsd   nb201_qqH(%esp),%xmm5   ## vcoul 
        mulsd   nb201_qqH(%esp),%xmm7
        mulsd  %xmm7,%xmm4              ## total fsH2 in xmm4 

        addsd  nb201_vctot(%esp),%xmm5

        movapd nb201_dxH2(%esp),%xmm0
        movapd nb201_dyH2(%esp),%xmm1
        movapd nb201_dzH2(%esp),%xmm2
        movlpd %xmm5,nb201_vctot(%esp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H2 forces 
        movapd nb201_fixH2(%esp),%xmm3
        movapd nb201_fiyH2(%esp),%xmm4
        movapd nb201_fizH2(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb201_fixH2(%esp)
        movlpd %xmm4,nb201_fiyH2(%esp)
        movlpd %xmm7,nb201_fizH2(%esp)

        movl nb201_faction(%ebp),%edi
        ## update j forces 
        addsd  nb201_fjx(%esp),%xmm0
        addsd  nb201_fjy(%esp),%xmm1
        addsd  nb201_fjz(%esp),%xmm2
        movlpd (%edi,%eax,8),%xmm3
        movlpd 8(%edi,%eax,8),%xmm4
        movlpd 16(%edi,%eax,8),%xmm5
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        movlpd %xmm3,(%edi,%eax,8)
        movlpd %xmm4,8(%edi,%eax,8)
        movlpd %xmm5,16(%edi,%eax,8)

_nb_kernel201_ia32_sse2.nb201_updateouterdata: 
        movl  nb201_ii3(%esp),%ecx
        movl  nb201_faction(%ebp),%edi
        movl  nb201_fshift(%ebp),%esi
        movl  nb201_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movapd nb201_fixO(%esp),%xmm0
        movapd nb201_fiyO(%esp),%xmm1
        movapd nb201_fizO(%esp),%xmm2

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
        movapd nb201_fixH1(%esp),%xmm0
        movapd nb201_fiyH1(%esp),%xmm1
        movapd nb201_fizH1(%esp),%xmm2

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
        movapd nb201_fixH2(%esp),%xmm0
        movapd nb201_fiyH2(%esp),%xmm1
        movapd nb201_fizH2(%esp),%xmm2

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
        movl nb201_n(%esp),%esi
        ## get group index for i particle 
        movl  nb201_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb201_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb201_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb201_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel201_ia32_sse2.nb201_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb201_n(%esp)
        jmp _nb_kernel201_ia32_sse2.nb201_outer
_nb_kernel201_ia32_sse2.nb201_outerend: 
        ## check if more outer neighborlists remain
        movl  nb201_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel201_ia32_sse2.nb201_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel201_ia32_sse2.nb201_threadloop
_nb_kernel201_ia32_sse2.nb201_end: 
        emms

        movl nb201_nouter(%esp),%eax
        movl nb201_ninner(%esp),%ebx
        movl nb201_outeriter(%ebp),%ecx
        movl nb201_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb201_salign(%esp),%eax
        addl %eax,%esp
        addl $744,%esp
        popl %edi
        popl %esi
    popl %edx
    popl %ecx
    popl %ebx
    popl %eax
        leave
        ret




.globl nb_kernel201nf_ia32_sse2
.globl _nb_kernel201nf_ia32_sse2
nb_kernel201nf_ia32_sse2:       
_nb_kernel201nf_ia32_sse2:      
.set nb201nf_p_nri, 8
.set nb201nf_iinr, 12
.set nb201nf_jindex, 16
.set nb201nf_jjnr, 20
.set nb201nf_shift, 24
.set nb201nf_shiftvec, 28
.set nb201nf_fshift, 32
.set nb201nf_gid, 36
.set nb201nf_pos, 40
.set nb201nf_faction, 44
.set nb201nf_charge, 48
.set nb201nf_p_facel, 52
.set nb201nf_argkrf, 56
.set nb201nf_argcrf, 60
.set nb201nf_Vc, 64
.set nb201nf_type, 68
.set nb201nf_p_ntype, 72
.set nb201nf_vdwparam, 76
.set nb201nf_Vvdw, 80
.set nb201nf_p_tabscale, 84
.set nb201nf_VFtab, 88
.set nb201nf_invsqrta, 92
.set nb201nf_dvda, 96
.set nb201nf_p_gbtabscale, 100
.set nb201nf_GBtab, 104
.set nb201nf_p_nthreads, 108
.set nb201nf_count, 112
.set nb201nf_mtx, 116
.set nb201nf_outeriter, 120
.set nb201nf_inneriter, 124
.set nb201nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb201nf_ixO, 0
.set nb201nf_iyO, 16
.set nb201nf_izO, 32
.set nb201nf_ixH1, 48
.set nb201nf_iyH1, 64
.set nb201nf_izH1, 80
.set nb201nf_ixH2, 96
.set nb201nf_iyH2, 112
.set nb201nf_izH2, 128
.set nb201nf_iqO, 144
.set nb201nf_iqH, 160
.set nb201nf_qqO, 176
.set nb201nf_qqH, 192
.set nb201nf_vctot, 208
.set nb201nf_half, 224
.set nb201nf_three, 240
.set nb201nf_krf, 256
.set nb201nf_crf, 272
.set nb201nf_krsqO, 288
.set nb201nf_krsqH1, 304
.set nb201nf_krsqH2, 320
.set nb201nf_is3, 336
.set nb201nf_ii3, 340
.set nb201nf_innerjjnr, 344
.set nb201nf_innerk, 348
.set nb201nf_n, 352
.set nb201nf_nn1, 356
.set nb201nf_nri, 360
.set nb201nf_nouter, 364
.set nb201nf_ninner, 368
.set nb201nf_salign, 372
        pushl %ebp
        movl %esp,%ebp
    pushl %eax
    pushl %ebx
    pushl %ecx
    pushl %edx
        pushl %esi
        pushl %edi
        subl $376,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb201nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb201nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb201nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb201nf_nouter(%esp)
        movl %eax,nb201nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb201nf_half(%esp)
        movl %ebx,nb201nf_half+4(%esp)
        movsd nb201nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb201nf_half(%esp)
        movapd %xmm3,nb201nf_three(%esp)

        movl nb201nf_argkrf(%ebp),%esi
        movl nb201nf_argcrf(%ebp),%edi
        movsd (%esi),%xmm5
        movsd (%edi),%xmm6
        shufpd $0,%xmm5,%xmm5
        shufpd $0,%xmm6,%xmm6
        movapd %xmm5,nb201nf_krf(%esp)
        movapd %xmm6,nb201nf_crf(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb201nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb201nf_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        movsd 8(%edx,%ebx,8),%xmm4
        movl nb201nf_p_facel(%ebp),%esi
        movsd (%esi),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb201nf_iqO(%esp)
        movapd %xmm4,nb201nf_iqH(%esp)

_nb_kernel201nf_ia32_sse2.nb201nf_threadloop: 
        movl  nb201nf_count(%ebp),%esi          ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel201nf_ia32_sse2.nb201nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel201nf_ia32_sse2.nb201nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb201nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb201nf_n(%esp)
        movl %ebx,nb201nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel201nf_ia32_sse2.nb201nf_outerstart
        jmp _nb_kernel201nf_ia32_sse2.nb201nf_end

_nb_kernel201nf_ia32_sse2.nb201nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb201nf_nouter(%esp),%ebx
        movl %ebx,nb201nf_nouter(%esp)

_nb_kernel201nf_ia32_sse2.nb201nf_outer: 
        movl  nb201nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 

        movl  nb201nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb201nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb201nf_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb201nf_ii3(%esp)

        addsd (%eax,%ebx,8),%xmm3
        addsd 8(%eax,%ebx,8),%xmm4
        addsd 16(%eax,%ebx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb201nf_ixO(%esp)
        movapd %xmm4,nb201nf_iyO(%esp)
        movapd %xmm5,nb201nf_izO(%esp)

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
        movapd %xmm0,nb201nf_ixH1(%esp)
        movapd %xmm1,nb201nf_iyH1(%esp)
        movapd %xmm2,nb201nf_izH1(%esp)
        movapd %xmm3,nb201nf_ixH2(%esp)
        movapd %xmm4,nb201nf_iyH2(%esp)
        movapd %xmm5,nb201nf_izH2(%esp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb201nf_vctot(%esp)

        movl  nb201nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx     ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb201nf_pos(%ebp),%esi
        movl  nb201nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb201nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb201nf_ninner(%esp),%ecx
        movl  %ecx,nb201nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb201nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel201nf_ia32_sse2.nb201nf_unroll_loop
        jmp   _nb_kernel201nf_ia32_sse2.nb201nf_checksingle
_nb_kernel201nf_ia32_sse2.nb201nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb201nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb201nf_innerjjnr(%esp)                 ## advance pointer (unrolled 2) 

        movl nb201nf_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb201nf_iqO(%esp),%xmm3
        mulpd  nb201nf_iqH(%esp),%xmm4
        movapd  %xmm3,nb201nf_qqO(%esp)
        movapd  %xmm4,nb201nf_qqH(%esp)

        movl nb201nf_pos(%ebp),%esi        ## base of pos[] 

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
        movapd nb201nf_ixO(%esp),%xmm4
        movapd nb201nf_iyO(%esp),%xmm5
        movapd nb201nf_izO(%esp),%xmm6

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
        movapd nb201nf_ixH1(%esp),%xmm4
        movapd nb201nf_iyH1(%esp),%xmm5
        movapd nb201nf_izH1(%esp),%xmm6

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
        movapd nb201nf_ixH2(%esp),%xmm3
        movapd nb201nf_iyH2(%esp),%xmm4
        movapd nb201nf_izH2(%esp),%xmm5

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

        mulpd  nb201nf_krf(%esp),%xmm0
        mulpd  nb201nf_krf(%esp),%xmm1
        mulpd  nb201nf_krf(%esp),%xmm2

        movapd %xmm0,nb201nf_krsqH2(%esp)
        movapd %xmm1,nb201nf_krsqH1(%esp)
        movapd %xmm2,nb201nf_krsqO(%esp)

        ## start with rsqO - put seed in xmm2 
        cvtpd2ps %xmm7,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb201nf_three(%esp),%xmm4
        mulpd   %xmm7,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb201nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb201nf_three(%esp),%xmm4
        subpd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb201nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvO in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb201nf_three(%esp),%xmm4
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb201nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb201nf_three(%esp),%xmm4
        subpd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb201nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb201nf_three(%esp),%xmm4
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb201nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb201nf_three(%esp),%xmm4
        subpd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb201nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do O interactions 
        movapd  nb201nf_krsqO(%esp),%xmm0
        addpd   %xmm0,%xmm7     ## xmm7=rinv+ krsq 
        subpd   nb201nf_crf(%esp),%xmm7
        mulpd   nb201nf_qqO(%esp),%xmm7   ## vcoul      
        addpd  nb201nf_vctot(%esp),%xmm7

        ## H1 interactions 
        movapd  nb201nf_krsqH1(%esp),%xmm0
        addpd   %xmm0,%xmm6     ## xmm6=rinv+ krsq 
        subpd   nb201nf_crf(%esp),%xmm6
        mulpd   nb201nf_qqH(%esp),%xmm6   ## vcoul 
        addpd  %xmm7,%xmm6

        ## H2 interactions 
        movapd  nb201nf_krsqH2(%esp),%xmm0
        addpd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        subpd   nb201nf_crf(%esp),%xmm5
        mulpd   nb201nf_qqH(%esp),%xmm5   ## vcoul 
        addpd  %xmm6,%xmm5
        movapd %xmm5,nb201nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $2,nb201nf_innerk(%esp)
        jl    _nb_kernel201nf_ia32_sse2.nb201nf_checksingle
        jmp   _nb_kernel201nf_ia32_sse2.nb201nf_unroll_loop
_nb_kernel201nf_ia32_sse2.nb201nf_checksingle: 
        movl  nb201nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel201nf_ia32_sse2.nb201nf_dosingle
        jmp   _nb_kernel201nf_ia32_sse2.nb201nf_updateouterdata
_nb_kernel201nf_ia32_sse2.nb201nf_dosingle: 
        movl  nb201nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb201nf_innerjjnr(%esp)

        movl nb201nf_charge(%ebp),%esi     ## base of charge[] 
        xorpd %xmm3,%xmm3
        movlpd (%esi,%eax,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb201nf_iqO(%esp),%xmm3
        mulpd  nb201nf_iqH(%esp),%xmm4
        movapd  %xmm3,nb201nf_qqO(%esp)
        movapd  %xmm4,nb201nf_qqH(%esp)

        movl nb201nf_pos(%ebp),%esi        ## base of pos[] 
        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb201nf_ixO(%esp),%xmm4
        movapd nb201nf_iyO(%esp),%xmm5
        movapd nb201nf_izO(%esp),%xmm6

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
        movapd nb201nf_ixH1(%esp),%xmm4
        movapd nb201nf_iyH1(%esp),%xmm5
        movapd nb201nf_izH1(%esp),%xmm6

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
        movapd nb201nf_ixH2(%esp),%xmm3
        movapd nb201nf_iyH2(%esp),%xmm4
        movapd nb201nf_izH2(%esp),%xmm5

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

        mulsd  nb201nf_krf(%esp),%xmm0
        mulsd  nb201nf_krf(%esp),%xmm1
        mulsd  nb201nf_krf(%esp),%xmm2

        movapd %xmm0,nb201nf_krsqH2(%esp)
        movapd %xmm1,nb201nf_krsqH1(%esp)
        movapd %xmm2,nb201nf_krsqO(%esp)

        ## start with rsqO - put seed in xmm2 
        cvtsd2ss %xmm7,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb201nf_three(%esp),%xmm4
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb201nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb201nf_three(%esp),%xmm4
        subsd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb201nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvO in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb201nf_three(%esp),%xmm4
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb201nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb201nf_three(%esp),%xmm4
        subsd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb201nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb201nf_three(%esp),%xmm4
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb201nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb201nf_three(%esp),%xmm4
        subsd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb201nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do O interactions 
        movapd  nb201nf_krsqO(%esp),%xmm0
        addsd   %xmm0,%xmm7     ## xmm7=rinv+ krsq 
        subsd   nb201nf_crf(%esp),%xmm7
        mulsd   nb201nf_qqO(%esp),%xmm7   ## vcoul      
        addsd  nb201nf_vctot(%esp),%xmm7

        ## H1 interactions 
        movapd  nb201nf_krsqH1(%esp),%xmm0
        addsd   %xmm0,%xmm6     ## xmm6=rinv+ krsq 
        subsd   nb201nf_crf(%esp),%xmm6
        mulsd   nb201nf_qqH(%esp),%xmm6   ## vcoul 
        addsd  %xmm7,%xmm6

        ## H2 interactions 
        movapd  nb201nf_krsqH2(%esp),%xmm0
        addsd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        subsd   nb201nf_crf(%esp),%xmm5
        mulsd   nb201nf_qqH(%esp),%xmm5   ## vcoul 
        addsd  %xmm6,%xmm5
        movlpd %xmm5,nb201nf_vctot(%esp)

_nb_kernel201nf_ia32_sse2.nb201nf_updateouterdata: 
        ## get n from stack
        movl nb201nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb201nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb201nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb201nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb201nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel201nf_ia32_sse2.nb201nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb201nf_n(%esp)
        jmp _nb_kernel201nf_ia32_sse2.nb201nf_outer
_nb_kernel201nf_ia32_sse2.nb201nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb201nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel201nf_ia32_sse2.nb201nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel201nf_ia32_sse2.nb201nf_threadloop
_nb_kernel201nf_ia32_sse2.nb201nf_end: 
        emms

        movl nb201nf_nouter(%esp),%eax
        movl nb201nf_ninner(%esp),%ebx
        movl nb201nf_outeriter(%ebp),%ecx
        movl nb201nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb201nf_salign(%esp),%eax
        addl %eax,%esp
        addl $376,%esp
        popl %edi
        popl %esi
    popl %edx
    popl %ecx
    popl %ebx
    popl %eax
        leave
        ret

