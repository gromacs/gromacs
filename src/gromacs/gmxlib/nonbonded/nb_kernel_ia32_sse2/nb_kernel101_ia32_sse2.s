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




.globl nb_kernel101_ia32_sse2
.globl _nb_kernel101_ia32_sse2
nb_kernel101_ia32_sse2: 
_nb_kernel101_ia32_sse2:        
.set nb101_p_nri, 8
.set nb101_iinr, 12
.set nb101_jindex, 16
.set nb101_jjnr, 20
.set nb101_shift, 24
.set nb101_shiftvec, 28
.set nb101_fshift, 32
.set nb101_gid, 36
.set nb101_pos, 40
.set nb101_faction, 44
.set nb101_charge, 48
.set nb101_p_facel, 52
.set nb101_argkrf, 56
.set nb101_argcrf, 60
.set nb101_Vc, 64
.set nb101_type, 68
.set nb101_p_ntype, 72
.set nb101_vdwparam, 76
.set nb101_Vvdw, 80
.set nb101_p_tabscale, 84
.set nb101_VFtab, 88
.set nb101_invsqrta, 92
.set nb101_dvda, 96
.set nb101_p_gbtabscale, 100
.set nb101_GBtab, 104
.set nb101_p_nthreads, 108
.set nb101_count, 112
.set nb101_mtx, 116
.set nb101_outeriter, 120
.set nb101_inneriter, 124
.set nb101_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb101_ixO, 0
.set nb101_iyO, 16
.set nb101_izO, 32
.set nb101_ixH1, 48
.set nb101_iyH1, 64
.set nb101_izH1, 80
.set nb101_ixH2, 96
.set nb101_iyH2, 112
.set nb101_izH2, 128
.set nb101_iqO, 144
.set nb101_iqH, 160
.set nb101_dxO, 176
.set nb101_dyO, 192
.set nb101_dzO, 208
.set nb101_dxH1, 224
.set nb101_dyH1, 240
.set nb101_dzH1, 256
.set nb101_dxH2, 272
.set nb101_dyH2, 288
.set nb101_dzH2, 304
.set nb101_qqO, 320
.set nb101_qqH, 336
.set nb101_vctot, 352
.set nb101_fixO, 368
.set nb101_fiyO, 384
.set nb101_fizO, 400
.set nb101_fixH1, 416
.set nb101_fiyH1, 432
.set nb101_fizH1, 448
.set nb101_fixH2, 464
.set nb101_fiyH2, 480
.set nb101_fizH2, 496
.set nb101_fjx, 512
.set nb101_fjy, 528
.set nb101_fjz, 544
.set nb101_half, 560
.set nb101_three, 576
.set nb101_is3, 592
.set nb101_ii3, 596
.set nb101_innerjjnr, 600
.set nb101_innerk, 604
.set nb101_n, 608
.set nb101_nn1, 612
.set nb101_nri, 616
.set nb101_nouter, 620
.set nb101_ninner, 624
.set nb101_salign, 628
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $632,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb101_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb101_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb101_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb101_nouter(%esp)
        movl %eax,nb101_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb101_half(%esp)
        movl %ebx,nb101_half+4(%esp)
        movsd nb101_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb101_half(%esp)
        movapd %xmm3,nb101_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb101_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb101_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        movsd 8(%edx,%ebx,8),%xmm4
        movl nb101_p_facel(%ebp),%esi
        movsd (%esi),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb101_iqO(%esp)
        movapd %xmm4,nb101_iqH(%esp)

_nb_kernel101_ia32_sse2.nb101_threadloop: 
        movl  nb101_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel101_ia32_sse2.nb101_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel101_ia32_sse2.nb101_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb101_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb101_n(%esp)
        movl %ebx,nb101_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel101_ia32_sse2.nb101_outerstart
        jmp _nb_kernel101_ia32_sse2.nb101_end

_nb_kernel101_ia32_sse2.nb101_outerstart: 
        ## ebx contains number of outer iterations
        addl nb101_nouter(%esp),%ebx
        movl %ebx,nb101_nouter(%esp)

_nb_kernel101_ia32_sse2.nb101_outer: 
        movl  nb101_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb101_is3(%esp)      ## store is3 

        movl  nb101_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb101_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb101_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb101_ii3(%esp)

        addsd (%eax,%ebx,8),%xmm3
        addsd 8(%eax,%ebx,8),%xmm4
        addsd 16(%eax,%ebx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb101_ixO(%esp)
        movapd %xmm4,nb101_iyO(%esp)
        movapd %xmm5,nb101_izO(%esp)

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
        movapd %xmm0,nb101_ixH1(%esp)
        movapd %xmm1,nb101_iyH1(%esp)
        movapd %xmm2,nb101_izH1(%esp)
        movapd %xmm3,nb101_ixH2(%esp)
        movapd %xmm4,nb101_iyH2(%esp)
        movapd %xmm5,nb101_izH2(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb101_vctot(%esp)
        movapd %xmm4,nb101_fixO(%esp)
        movapd %xmm4,nb101_fiyO(%esp)
        movapd %xmm4,nb101_fizO(%esp)
        movapd %xmm4,nb101_fixH1(%esp)
        movapd %xmm4,nb101_fiyH1(%esp)
        movapd %xmm4,nb101_fizH1(%esp)
        movapd %xmm4,nb101_fixH2(%esp)
        movapd %xmm4,nb101_fiyH2(%esp)
        movapd %xmm4,nb101_fizH2(%esp)

        movl  nb101_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx     ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb101_pos(%ebp),%esi
        movl  nb101_faction(%ebp),%edi
        movl  nb101_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb101_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb101_ninner(%esp),%ecx
        movl  %ecx,nb101_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb101_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel101_ia32_sse2.nb101_unroll_loop
        jmp   _nb_kernel101_ia32_sse2.nb101_checksingle
_nb_kernel101_ia32_sse2.nb101_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb101_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb101_innerjjnr(%esp)             ## advance pointer (unrolled 2) 

        movl nb101_charge(%ebp),%esi     ## base of charge[] 


        movlpd (%esi,%eax,8),%xmm6      ## jq A 
        movhpd (%esi,%ebx,8),%xmm6      ## jq B 
        movapd nb101_iqO(%esp),%xmm3
        movapd nb101_iqH(%esp),%xmm4
        mulpd %xmm6,%xmm3               ## qqO 
        mulpd %xmm6,%xmm4               ## qqH 

        movapd  %xmm3,nb101_qqO(%esp)
        movapd  %xmm4,nb101_qqH(%esp)

        movl nb101_pos(%ebp),%esi        ## base of pos[] 

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
        movapd nb101_ixO(%esp),%xmm4
        movapd nb101_iyO(%esp),%xmm5
        movapd nb101_izO(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb101_dxO(%esp)
        movapd %xmm5,nb101_dyO(%esp)
        movapd %xmm6,nb101_dzO(%esp)
        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb101_ixH1(%esp),%xmm4
        movapd nb101_iyH1(%esp),%xmm5
        movapd nb101_izH1(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb101_dxH1(%esp)
        movapd %xmm5,nb101_dyH1(%esp)
        movapd %xmm6,nb101_dzH1(%esp)
        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm6
        addpd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb101_ixH2(%esp),%xmm3
        movapd nb101_iyH2(%esp),%xmm4
        movapd nb101_izH2(%esp),%xmm5

        ## calc dr 
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5

        ## store dr 
        movapd %xmm3,nb101_dxH2(%esp)
        movapd %xmm4,nb101_dyH2(%esp)
        movapd %xmm5,nb101_dzH2(%esp)
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
        movapd  nb101_three(%esp),%xmm4
        mulpd   %xmm7,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb101_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb101_three(%esp),%xmm4
        subpd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb101_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvO in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb101_three(%esp),%xmm4
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb101_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb101_three(%esp),%xmm4
        subpd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb101_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb101_three(%esp),%xmm4
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb101_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb101_three(%esp),%xmm4
        subpd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb101_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do O interactions 
        movapd  %xmm7,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq 
        mulpd  nb101_qqO(%esp),%xmm7    ## xmm7=vcoul 

        mulpd  %xmm7,%xmm4      ## total fsO in xmm4 

        addpd  nb101_vctot(%esp),%xmm7

        movapd %xmm7,nb101_vctot(%esp)

        movapd nb101_dxO(%esp),%xmm0
        movapd nb101_dyO(%esp),%xmm1
        movapd nb101_dzO(%esp),%xmm2
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        ## update O forces 
        movapd nb101_fixO(%esp),%xmm3
        movapd nb101_fiyO(%esp),%xmm4
        movapd nb101_fizO(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb101_fixO(%esp)
        movapd %xmm4,nb101_fiyO(%esp)
        movapd %xmm7,nb101_fizO(%esp)
        ## update j forces with water O 
        movapd %xmm0,nb101_fjx(%esp)
        movapd %xmm1,nb101_fjy(%esp)
        movapd %xmm2,nb101_fjz(%esp)

        ## H1 interactions 
        movapd  %xmm6,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        mulpd  nb101_qqH(%esp),%xmm6    ## xmm6=vcoul 
        mulpd  %xmm6,%xmm4              ## total fsH1 in xmm4 

        addpd  nb101_vctot(%esp),%xmm6

        movapd nb101_dxH1(%esp),%xmm0
        movapd nb101_dyH1(%esp),%xmm1
        movapd nb101_dzH1(%esp),%xmm2
        movapd %xmm6,nb101_vctot(%esp)
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb101_fixH1(%esp),%xmm3
        movapd nb101_fiyH1(%esp),%xmm4
        movapd nb101_fizH1(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb101_fixH1(%esp)
        movapd %xmm4,nb101_fiyH1(%esp)
        movapd %xmm7,nb101_fizH1(%esp)
        ## update j forces with water H1 
        addpd  nb101_fjx(%esp),%xmm0
        addpd  nb101_fjy(%esp),%xmm1
        addpd  nb101_fjz(%esp),%xmm2
        movapd %xmm0,nb101_fjx(%esp)
        movapd %xmm1,nb101_fjy(%esp)
        movapd %xmm2,nb101_fjz(%esp)

        ## H2 interactions 
        movapd  %xmm5,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        mulpd  nb101_qqH(%esp),%xmm5    ## xmm5=vcoul 
        mulpd  %xmm5,%xmm4              ## total fsH1 in xmm4 

        addpd  nb101_vctot(%esp),%xmm5

        movapd nb101_dxH2(%esp),%xmm0
        movapd nb101_dyH2(%esp),%xmm1
        movapd nb101_dzH2(%esp),%xmm2
        movapd %xmm5,nb101_vctot(%esp)
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        ## update H2 forces 
        movapd nb101_fixH2(%esp),%xmm3
        movapd nb101_fiyH2(%esp),%xmm4
        movapd nb101_fizH2(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb101_fixH2(%esp)
        movapd %xmm4,nb101_fiyH2(%esp)
        movapd %xmm7,nb101_fizH2(%esp)

        movl nb101_faction(%ebp),%edi
        ## update j forces 
        addpd  nb101_fjx(%esp),%xmm0
        addpd  nb101_fjy(%esp),%xmm1
        addpd  nb101_fjz(%esp),%xmm2

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
        subl $2,nb101_innerk(%esp)
        jl    _nb_kernel101_ia32_sse2.nb101_checksingle
        jmp   _nb_kernel101_ia32_sse2.nb101_unroll_loop
_nb_kernel101_ia32_sse2.nb101_checksingle:      
        movl  nb101_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel101_ia32_sse2.nb101_dosingle
        jmp    _nb_kernel101_ia32_sse2.nb101_updateouterdata
_nb_kernel101_ia32_sse2.nb101_dosingle: 
        movl  nb101_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movl nb101_charge(%ebp),%esi     ## base of charge[] 
        xorpd %xmm6,%xmm6
        movlpd (%esi,%eax,8),%xmm6      ## jq A 

        movapd nb101_iqO(%esp),%xmm3
        movapd nb101_iqH(%esp),%xmm4
        mulsd %xmm6,%xmm3               ## qqO 
        mulsd %xmm6,%xmm4               ## qqH 

        movapd  %xmm3,nb101_qqO(%esp)
        movapd  %xmm4,nb101_qqH(%esp)

        movl nb101_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2        
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb101_ixO(%esp),%xmm4
        movapd nb101_iyO(%esp),%xmm5
        movapd nb101_izO(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb101_dxO(%esp)
        movapd %xmm5,nb101_dyO(%esp)
        movapd %xmm6,nb101_dzO(%esp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb101_ixH1(%esp),%xmm4
        movapd nb101_iyH1(%esp),%xmm5
        movapd nb101_izH1(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb101_dxH1(%esp)
        movapd %xmm5,nb101_dyH1(%esp)
        movapd %xmm6,nb101_dzH1(%esp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm6
        addsd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb101_ixH2(%esp),%xmm3
        movapd nb101_iyH2(%esp),%xmm4
        movapd nb101_izH2(%esp),%xmm5

        ## calc dr 
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5

        ## store dr 
        movapd %xmm3,nb101_dxH2(%esp)
        movapd %xmm4,nb101_dyH2(%esp)
        movapd %xmm5,nb101_dzH2(%esp)
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
        movapd  nb101_three(%esp),%xmm4
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb101_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb101_three(%esp),%xmm4
        subsd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb101_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvO in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb101_three(%esp),%xmm4
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb101_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb101_three(%esp),%xmm4
        subsd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb101_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb101_three(%esp),%xmm4
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb101_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb101_three(%esp),%xmm4
        subsd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb101_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do O interactions 
        movapd  %xmm7,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq 
        mulsd  nb101_qqO(%esp),%xmm7    ## xmm7=vcoul 

        mulsd  %xmm7,%xmm4      ## total fsO in xmm4 

        addsd  nb101_vctot(%esp),%xmm7

        movlpd %xmm7,nb101_vctot(%esp)

        movapd nb101_dxO(%esp),%xmm0
        movapd nb101_dyO(%esp),%xmm1
        movapd nb101_dzO(%esp),%xmm2
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update O forces 
        movapd nb101_fixO(%esp),%xmm3
        movapd nb101_fiyO(%esp),%xmm4
        movapd nb101_fizO(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb101_fixO(%esp)
        movlpd %xmm4,nb101_fiyO(%esp)
        movlpd %xmm7,nb101_fizO(%esp)
        ## update j forces with water O 
        movlpd %xmm0,nb101_fjx(%esp)
        movlpd %xmm1,nb101_fjy(%esp)
        movlpd %xmm2,nb101_fjz(%esp)

        ## H1 interactions 
        movapd  %xmm6,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        mulsd  nb101_qqH(%esp),%xmm6    ## xmm6=vcoul 
        mulsd  %xmm6,%xmm4              ## total fsH1 in xmm4 

        addsd  nb101_vctot(%esp),%xmm6

        movapd nb101_dxH1(%esp),%xmm0
        movapd nb101_dyH1(%esp),%xmm1
        movapd nb101_dzH1(%esp),%xmm2
        movlpd %xmm6,nb101_vctot(%esp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb101_fixH1(%esp),%xmm3
        movapd nb101_fiyH1(%esp),%xmm4
        movapd nb101_fizH1(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb101_fixH1(%esp)
        movlpd %xmm4,nb101_fiyH1(%esp)
        movlpd %xmm7,nb101_fizH1(%esp)
        ## update j forces with water H1 
        addsd  nb101_fjx(%esp),%xmm0
        addsd  nb101_fjy(%esp),%xmm1
        addsd  nb101_fjz(%esp),%xmm2
        movsd %xmm0,nb101_fjx(%esp)
        movsd %xmm1,nb101_fjy(%esp)
        movsd %xmm2,nb101_fjz(%esp)

        ## H2 interactions 
        movapd  %xmm5,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        mulsd  nb101_qqH(%esp),%xmm5    ## xmm5=vcoul 
        mulsd  %xmm5,%xmm4              ## total fsH1 in xmm4 

        addsd  nb101_vctot(%esp),%xmm5

        movapd nb101_dxH2(%esp),%xmm0
        movapd nb101_dyH2(%esp),%xmm1
        movapd nb101_dzH2(%esp),%xmm2
        movlpd %xmm5,nb101_vctot(%esp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H2 forces 
        movapd nb101_fixH2(%esp),%xmm3
        movapd nb101_fiyH2(%esp),%xmm4
        movapd nb101_fizH2(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb101_fixH2(%esp)
        movlpd %xmm4,nb101_fiyH2(%esp)
        movlpd %xmm7,nb101_fizH2(%esp)

        movl nb101_faction(%ebp),%edi
        ## update j forces 
        addsd  nb101_fjx(%esp),%xmm0
        addsd  nb101_fjy(%esp),%xmm1
        addsd  nb101_fjz(%esp),%xmm2

        movlpd (%edi,%eax,8),%xmm3
        movlpd 8(%edi,%eax,8),%xmm4
        movlpd 16(%edi,%eax,8),%xmm5
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        movlpd %xmm3,(%edi,%eax,8)
        movlpd %xmm4,8(%edi,%eax,8)
        movlpd %xmm5,16(%edi,%eax,8)

_nb_kernel101_ia32_sse2.nb101_updateouterdata: 
        movl  nb101_ii3(%esp),%ecx
        movl  nb101_faction(%ebp),%edi
        movl  nb101_fshift(%ebp),%esi
        movl  nb101_is3(%esp),%edx

        ## accumulate Oi forces in xmm0, xmm1, xmm2 
        movapd nb101_fixO(%esp),%xmm0
        movapd nb101_fiyO(%esp),%xmm1
        movapd nb101_fizO(%esp),%xmm2

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
        movapd nb101_fixH1(%esp),%xmm0
        movapd nb101_fiyH1(%esp),%xmm1
        movapd nb101_fizH1(%esp),%xmm2

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
        movapd nb101_fixH2(%esp),%xmm0
        movapd nb101_fiyH2(%esp),%xmm1
        movapd nb101_fizH2(%esp),%xmm2

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
        movl nb101_n(%esp),%esi
        ## get group index for i particle 
        movl  nb101_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb101_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 

        ## add earlier value from mem 
        movl  nb101_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb101_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel101_ia32_sse2.nb101_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb101_n(%esp)
        jmp _nb_kernel101_ia32_sse2.nb101_outer
_nb_kernel101_ia32_sse2.nb101_outerend: 
        ## check if more outer neighborlists remain
        movl  nb101_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel101_ia32_sse2.nb101_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel101_ia32_sse2.nb101_threadloop
_nb_kernel101_ia32_sse2.nb101_end: 
        emms

        movl nb101_nouter(%esp),%eax
        movl nb101_ninner(%esp),%ebx
        movl nb101_outeriter(%ebp),%ecx
        movl nb101_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb101_salign(%esp),%eax
        addl %eax,%esp
        addl $632,%esp
        popl %edi
        popl %esi
    popl %edx
    popl %ecx
    popl %ebx
    popl %eax
        leave
        ret



.globl nb_kernel101nf_ia32_sse2
.globl _nb_kernel101nf_ia32_sse2
nb_kernel101nf_ia32_sse2:       
_nb_kernel101nf_ia32_sse2:      
.set nb101nf_p_nri, 8
.set nb101nf_iinr, 12
.set nb101nf_jindex, 16
.set nb101nf_jjnr, 20
.set nb101nf_shift, 24
.set nb101nf_shiftvec, 28
.set nb101nf_fshift, 32
.set nb101nf_gid, 36
.set nb101nf_pos, 40
.set nb101nf_faction, 44
.set nb101nf_charge, 48
.set nb101nf_p_facel, 52
.set nb101nf_argkrf, 56
.set nb101nf_argcrf, 60
.set nb101nf_Vc, 64
.set nb101nf_type, 68
.set nb101nf_p_ntype, 72
.set nb101nf_vdwparam, 76
.set nb101nf_Vvdw, 80
.set nb101nf_p_tabscale, 84
.set nb101nf_VFtab, 88
.set nb101nf_invsqrta, 92
.set nb101nf_dvda, 96
.set nb101nf_p_gbtabscale, 100
.set nb101nf_GBtab, 104
.set nb101nf_p_nthreads, 108
.set nb101nf_count, 112
.set nb101nf_mtx, 116
.set nb101nf_outeriter, 120
.set nb101nf_inneriter, 124
.set nb101nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb101nf_ixO, 0
.set nb101nf_iyO, 16
.set nb101nf_izO, 32
.set nb101nf_ixH1, 48
.set nb101nf_iyH1, 64
.set nb101nf_izH1, 80
.set nb101nf_ixH2, 96
.set nb101nf_iyH2, 112
.set nb101nf_izH2, 128
.set nb101nf_iqO, 144
.set nb101nf_iqH, 160
.set nb101nf_qqO, 176
.set nb101nf_qqH, 192
.set nb101nf_vctot, 208
.set nb101nf_half, 224
.set nb101nf_three, 240
.set nb101nf_is3, 256
.set nb101nf_ii3, 260
.set nb101nf_innerjjnr, 264
.set nb101nf_innerk, 268
.set nb101nf_n, 272
.set nb101nf_nn1, 276
.set nb101nf_nri, 280
.set nb101nf_nouter, 284
.set nb101nf_ninner, 288
.set nb101nf_salign, 292
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $296,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb101nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb101nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb101nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb101nf_nouter(%esp)
        movl %eax,nb101nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb101nf_half(%esp)
        movl %ebx,nb101nf_half+4(%esp)
        movsd nb101nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb101nf_half(%esp)
        movapd %xmm3,nb101nf_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb101nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb101nf_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        movsd 8(%edx,%ebx,8),%xmm4
        movl nb101nf_p_facel(%ebp),%esi
        movsd (%esi),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb101nf_iqO(%esp)
        movapd %xmm4,nb101nf_iqH(%esp)

_nb_kernel101nf_ia32_sse2.nb101nf_threadloop: 
        movl  nb101nf_count(%ebp),%esi          ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel101nf_ia32_sse2.nb101nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel101nf_ia32_sse2.nb101nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb101nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb101nf_n(%esp)
        movl %ebx,nb101nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel101nf_ia32_sse2.nb101nf_outerstart
        jmp _nb_kernel101nf_ia32_sse2.nb101nf_end

_nb_kernel101nf_ia32_sse2.nb101nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb101nf_nouter(%esp),%ebx
        movl %ebx,nb101nf_nouter(%esp)

_nb_kernel101nf_ia32_sse2.nb101nf_outer: 
        movl  nb101nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 

        movl  nb101nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb101nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb101nf_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb101nf_ii3(%esp)

        addsd (%eax,%ebx,8),%xmm3
        addsd 8(%eax,%ebx,8),%xmm4
        addsd 16(%eax,%ebx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb101nf_ixO(%esp)
        movapd %xmm4,nb101nf_iyO(%esp)
        movapd %xmm5,nb101nf_izO(%esp)

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
        movapd %xmm0,nb101nf_ixH1(%esp)
        movapd %xmm1,nb101nf_iyH1(%esp)
        movapd %xmm2,nb101nf_izH1(%esp)
        movapd %xmm3,nb101nf_ixH2(%esp)
        movapd %xmm4,nb101nf_iyH2(%esp)
        movapd %xmm5,nb101nf_izH2(%esp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb101nf_vctot(%esp)

        movl  nb101nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb101nf_pos(%ebp),%esi
        movl  nb101nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb101nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb101nf_ninner(%esp),%ecx
        movl  %ecx,nb101nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb101nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel101nf_ia32_sse2.nb101nf_unroll_loop
        jmp   _nb_kernel101nf_ia32_sse2.nb101nf_checksingle
_nb_kernel101nf_ia32_sse2.nb101nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb101nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb101nf_innerjjnr(%esp)             ## advance pointer (unrolled 2) 

        movl nb101nf_charge(%ebp),%esi     ## base of charge[] 


        movlpd (%esi,%eax,8),%xmm6      ## jq A 
        movhpd (%esi,%ebx,8),%xmm6      ## jq B 
        movapd nb101nf_iqO(%esp),%xmm3
        movapd nb101nf_iqH(%esp),%xmm4
        mulpd %xmm6,%xmm3               ## qqO 
        mulpd %xmm6,%xmm4               ## qqH 

        movapd  %xmm3,nb101nf_qqO(%esp)
        movapd  %xmm4,nb101nf_qqH(%esp)

        movl nb101nf_pos(%ebp),%esi        ## base of pos[] 

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
        movapd nb101nf_ixO(%esp),%xmm4
        movapd nb101nf_iyO(%esp),%xmm5
        movapd nb101nf_izO(%esp),%xmm6

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
        movapd nb101nf_ixH1(%esp),%xmm4
        movapd nb101nf_iyH1(%esp),%xmm5
        movapd nb101nf_izH1(%esp),%xmm6

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
        movapd nb101nf_ixH2(%esp),%xmm3
        movapd nb101nf_iyH2(%esp),%xmm4
        movapd nb101nf_izH2(%esp),%xmm5

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
        movapd  nb101nf_three(%esp),%xmm4
        mulpd   %xmm7,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb101nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb101nf_three(%esp),%xmm4
        subpd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb101nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvO in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb101nf_three(%esp),%xmm4
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb101nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb101nf_three(%esp),%xmm4
        subpd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb101nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb101nf_three(%esp),%xmm4
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb101nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb101nf_three(%esp),%xmm4
        subpd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb101nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do O interactions 
        mulpd  nb101nf_qqO(%esp),%xmm7          ## xmm7=vcoul 
        addpd  nb101nf_vctot(%esp),%xmm7
        movapd %xmm7,nb101nf_vctot(%esp)

        ## H1 interactions 
        mulpd  nb101nf_qqH(%esp),%xmm6          ## xmm6=vcoul 
        addpd  nb101nf_vctot(%esp),%xmm6
        movapd %xmm6,nb101nf_vctot(%esp)

        ## H2 interactions 
        mulpd  nb101nf_qqH(%esp),%xmm5          ## xmm5=vcoul 
        addpd  nb101nf_vctot(%esp),%xmm5
        movapd %xmm5,nb101nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $2,nb101nf_innerk(%esp)
        jl    _nb_kernel101nf_ia32_sse2.nb101nf_checksingle
        jmp   _nb_kernel101nf_ia32_sse2.nb101nf_unroll_loop
_nb_kernel101nf_ia32_sse2.nb101nf_checksingle:  
        movl  nb101nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel101nf_ia32_sse2.nb101nf_dosingle
        jmp   _nb_kernel101nf_ia32_sse2.nb101nf_updateouterdata
_nb_kernel101nf_ia32_sse2.nb101nf_dosingle: 
        movl  nb101nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movl nb101nf_charge(%ebp),%esi     ## base of charge[] 
        xorpd %xmm6,%xmm6
        movlpd (%esi,%eax,8),%xmm6      ## jq A 

        movapd nb101nf_iqO(%esp),%xmm3
        movapd nb101nf_iqH(%esp),%xmm4
        mulsd %xmm6,%xmm3               ## qqO 
        mulsd %xmm6,%xmm4               ## qqH 

        movapd  %xmm3,nb101nf_qqO(%esp)
        movapd  %xmm4,nb101nf_qqH(%esp)

        movl nb101nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2        
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb101nf_ixO(%esp),%xmm4
        movapd nb101nf_iyO(%esp),%xmm5
        movapd nb101nf_izO(%esp),%xmm6

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
        movapd nb101nf_ixH1(%esp),%xmm4
        movapd nb101nf_iyH1(%esp),%xmm5
        movapd nb101nf_izH1(%esp),%xmm6

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
        movapd nb101nf_ixH2(%esp),%xmm3
        movapd nb101nf_iyH2(%esp),%xmm4
        movapd nb101nf_izH2(%esp),%xmm5

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
        movapd  nb101nf_three(%esp),%xmm4
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb101nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb101nf_three(%esp),%xmm4
        subsd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb101nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvO in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb101nf_three(%esp),%xmm4
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb101nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb101nf_three(%esp),%xmm4
        subsd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb101nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb101nf_three(%esp),%xmm4
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb101nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb101nf_three(%esp),%xmm4
        subsd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb101nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do O interactions 
        mulsd  nb101nf_qqO(%esp),%xmm7          ## xmm7=vcoul 
        addsd  nb101nf_vctot(%esp),%xmm7
        movlpd %xmm7,nb101nf_vctot(%esp)

        ## H1 interactions 
        mulsd  nb101nf_qqH(%esp),%xmm6          ## xmm6=vcoul 
        addsd  nb101nf_vctot(%esp),%xmm6
        movlpd %xmm6,nb101nf_vctot(%esp)

        ## H2 interactions 
        mulsd  nb101nf_qqH(%esp),%xmm5          ## xmm5=vcoul 
        addsd  nb101nf_vctot(%esp),%xmm5
        movlpd %xmm5,nb101nf_vctot(%esp)

_nb_kernel101nf_ia32_sse2.nb101nf_updateouterdata: 
        ## get n from stack
        movl nb101nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb101nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        movapd nb101nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 

        ## add earlier value from mem 
        movl  nb101nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb101nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel101nf_ia32_sse2.nb101nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb101nf_n(%esp)
        jmp _nb_kernel101nf_ia32_sse2.nb101nf_outer
_nb_kernel101nf_ia32_sse2.nb101nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb101nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel101nf_ia32_sse2.nb101nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel101nf_ia32_sse2.nb101nf_threadloop
_nb_kernel101nf_ia32_sse2.nb101nf_end: 
        emms

        movl nb101nf_nouter(%esp),%eax
        movl nb101nf_ninner(%esp),%ebx
        movl nb101nf_outeriter(%ebp),%ecx
        movl nb101nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb101nf_salign(%esp),%eax
        addl %eax,%esp
        addl $296,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret


