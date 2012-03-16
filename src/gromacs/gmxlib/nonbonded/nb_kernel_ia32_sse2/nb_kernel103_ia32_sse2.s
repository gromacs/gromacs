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




.globl nb_kernel103_ia32_sse2
.globl _nb_kernel103_ia32_sse2
nb_kernel103_ia32_sse2: 
_nb_kernel103_ia32_sse2:        
.set nb103_p_nri, 8
.set nb103_iinr, 12
.set nb103_jindex, 16
.set nb103_jjnr, 20
.set nb103_shift, 24
.set nb103_shiftvec, 28
.set nb103_fshift, 32
.set nb103_gid, 36
.set nb103_pos, 40
.set nb103_faction, 44
.set nb103_charge, 48
.set nb103_p_facel, 52
.set nb103_argkrf, 56
.set nb103_argcrf, 60
.set nb103_Vc, 64
.set nb103_type, 68
.set nb103_p_ntype, 72
.set nb103_vdwparam, 76
.set nb103_Vvdw, 80
.set nb103_p_tabscale, 84
.set nb103_VFtab, 88
.set nb103_invsqrta, 92
.set nb103_dvda, 96
.set nb103_p_gbtabscale, 100
.set nb103_GBtab, 104
.set nb103_p_nthreads, 108
.set nb103_count, 112
.set nb103_mtx, 116
.set nb103_outeriter, 120
.set nb103_inneriter, 124
.set nb103_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb103_ixH1, 0
.set nb103_iyH1, 16
.set nb103_izH1, 32
.set nb103_ixH2, 48
.set nb103_iyH2, 64
.set nb103_izH2, 80
.set nb103_ixM, 96
.set nb103_iyM, 112
.set nb103_izM, 128
.set nb103_iqM, 144
.set nb103_iqH, 160
.set nb103_dxH1, 176
.set nb103_dyH1, 192
.set nb103_dzH1, 208
.set nb103_dxH2, 224
.set nb103_dyH2, 240
.set nb103_dzH2, 256
.set nb103_dxM, 272
.set nb103_dyM, 288
.set nb103_dzM, 304
.set nb103_qqM, 320
.set nb103_qqH, 336
.set nb103_vctot, 352
.set nb103_fixM, 368
.set nb103_fiyM, 384
.set nb103_fizM, 400
.set nb103_fixH1, 416
.set nb103_fiyH1, 432
.set nb103_fizH1, 448
.set nb103_fixH2, 464
.set nb103_fiyH2, 480
.set nb103_fizH2, 496
.set nb103_fjx, 512
.set nb103_fjy, 528
.set nb103_fjz, 544
.set nb103_half, 560
.set nb103_three, 576
.set nb103_is3, 592
.set nb103_ii3, 596
.set nb103_innerjjnr, 600
.set nb103_innerk, 604
.set nb103_n, 608
.set nb103_nn1, 612
.set nb103_nri, 616
.set nb103_nouter, 620
.set nb103_ninner, 624
.set nb103_salign, 628
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
        movl %eax,nb103_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb103_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb103_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb103_nouter(%esp)
        movl %eax,nb103_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb103_half(%esp)
        movl %ebx,nb103_half+4(%esp)
        movsd nb103_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb103_half(%esp)
        movapd %xmm3,nb103_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb103_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb103_charge(%ebp),%edx
        movsd 8(%edx,%ebx,8),%xmm3
        movsd 24(%edx,%ebx,8),%xmm4
        movl nb103_p_facel(%ebp),%esi
        movsd (%esi),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb103_iqH(%esp)
        movapd %xmm4,nb103_iqM(%esp)

_nb_kernel103_ia32_sse2.nb103_threadloop: 
        movl  nb103_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel103_ia32_sse2.nb103_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel103_ia32_sse2.nb103_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb103_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb103_n(%esp)
        movl %ebx,nb103_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel103_ia32_sse2.nb103_outerstart
        jmp _nb_kernel103_ia32_sse2.nb103_end

_nb_kernel103_ia32_sse2.nb103_outerstart: 
        ## ebx contains number of outer iterations
        addl nb103_nouter(%esp),%ebx
        movl %ebx,nb103_nouter(%esp)

_nb_kernel103_ia32_sse2.nb103_outer: 
        movl  nb103_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb103_is3(%esp)      ## store is3 

        movl  nb103_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb103_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb103_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb103_ii3(%esp)

        addsd 24(%eax,%ebx,8),%xmm3
        addsd 32(%eax,%ebx,8),%xmm4
        addsd 40(%eax,%ebx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb103_ixH1(%esp)
        movapd %xmm4,nb103_iyH1(%esp)
        movapd %xmm5,nb103_izH1(%esp)

        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 48(%eax,%ebx,8),%xmm0
        addsd 56(%eax,%ebx,8),%xmm1
        addsd 64(%eax,%ebx,8),%xmm2
        addsd 72(%eax,%ebx,8),%xmm3
        addsd 80(%eax,%ebx,8),%xmm4
        addsd 88(%eax,%ebx,8),%xmm5

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm0,nb103_ixH2(%esp)
        movapd %xmm1,nb103_iyH2(%esp)
        movapd %xmm2,nb103_izH2(%esp)
        movapd %xmm3,nb103_ixM(%esp)
        movapd %xmm4,nb103_iyM(%esp)
        movapd %xmm5,nb103_izM(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb103_vctot(%esp)
        movapd %xmm4,nb103_fixM(%esp)
        movapd %xmm4,nb103_fiyM(%esp)
        movapd %xmm4,nb103_fizM(%esp)
        movapd %xmm4,nb103_fixH1(%esp)
        movapd %xmm4,nb103_fiyH1(%esp)
        movapd %xmm4,nb103_fizH1(%esp)
        movapd %xmm4,nb103_fixH2(%esp)
        movapd %xmm4,nb103_fiyH2(%esp)
        movapd %xmm4,nb103_fizH2(%esp)

        movl  nb103_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx     ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb103_pos(%ebp),%esi
        movl  nb103_faction(%ebp),%edi
        movl  nb103_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb103_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb103_ninner(%esp),%ecx
        movl  %ecx,nb103_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb103_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel103_ia32_sse2.nb103_unroll_loop
        jmp   _nb_kernel103_ia32_sse2.nb103_checksingle
_nb_kernel103_ia32_sse2.nb103_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb103_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb103_innerjjnr(%esp)             ## advance pointer (unrolled 2) 

        movl nb103_charge(%ebp),%esi     ## base of charge[] 


        movlpd (%esi,%eax,8),%xmm6      ## jq A 
        movhpd (%esi,%ebx,8),%xmm6      ## jq B 
        movapd nb103_iqM(%esp),%xmm3
        movapd nb103_iqH(%esp),%xmm4
        mulpd %xmm6,%xmm3               ## qqM 
        mulpd %xmm6,%xmm4               ## qqH 

        movapd  %xmm3,nb103_qqM(%esp)
        movapd  %xmm4,nb103_qqH(%esp)

        movl nb103_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move two coordinates to xmm0-xmm2    
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2
        movhpd (%esi,%ebx,8),%xmm0
        movhpd 8(%esi,%ebx,8),%xmm1
        movhpd 16(%esi,%ebx,8),%xmm2

        ## move ixM-izM to xmm4-xmm6 
        movapd nb103_ixM(%esp),%xmm4
        movapd nb103_iyM(%esp),%xmm5
        movapd nb103_izM(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb103_dxM(%esp)
        movapd %xmm5,nb103_dyM(%esp)
        movapd %xmm6,nb103_dzM(%esp)
        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb103_ixH1(%esp),%xmm4
        movapd nb103_iyH1(%esp),%xmm5
        movapd nb103_izH1(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb103_dxH1(%esp)
        movapd %xmm5,nb103_dyH1(%esp)
        movapd %xmm6,nb103_dzH1(%esp)
        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm6
        addpd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb103_ixH2(%esp),%xmm3
        movapd nb103_iyH2(%esp),%xmm4
        movapd nb103_izH2(%esp),%xmm5

        ## calc dr 
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5

        ## store dr 
        movapd %xmm3,nb103_dxH2(%esp)
        movapd %xmm4,nb103_dyH2(%esp)
        movapd %xmm5,nb103_dzH2(%esp)
        ## square it 
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        addpd %xmm4,%xmm5
        addpd %xmm3,%xmm5
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## start with rsqM - put seed in xmm2 
        cvtpd2ps %xmm7,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb103_three(%esp),%xmm4
        mulpd   %xmm7,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb103_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb103_three(%esp),%xmm4
        subpd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb103_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvM in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb103_three(%esp),%xmm4
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb103_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb103_three(%esp),%xmm4
        subpd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb103_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb103_three(%esp),%xmm4
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb103_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb103_three(%esp),%xmm4
        subpd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb103_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do M interactions 
        movapd  %xmm7,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq 
        mulpd  nb103_qqM(%esp),%xmm7    ## xmm7=vcoul 

        mulpd  %xmm7,%xmm4      ## total fsO in xmm4 

        addpd  nb103_vctot(%esp),%xmm7

        movapd %xmm7,nb103_vctot(%esp)

        movapd nb103_dxM(%esp),%xmm0
        movapd nb103_dyM(%esp),%xmm1
        movapd nb103_dzM(%esp),%xmm2
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        ## update M forces 
        movapd nb103_fixM(%esp),%xmm3
        movapd nb103_fiyM(%esp),%xmm4
        movapd nb103_fizM(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb103_fixM(%esp)
        movapd %xmm4,nb103_fiyM(%esp)
        movapd %xmm7,nb103_fizM(%esp)
        ## update j forces with water M
        movapd %xmm0,nb103_fjx(%esp)
        movapd %xmm1,nb103_fjy(%esp)
        movapd %xmm2,nb103_fjz(%esp)

        ## H1 interactions 
        movapd  %xmm6,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        mulpd  nb103_qqH(%esp),%xmm6    ## xmm6=vcoul 
        mulpd  %xmm6,%xmm4              ## total fsH1 in xmm4 

        addpd  nb103_vctot(%esp),%xmm6

        movapd nb103_dxH1(%esp),%xmm0
        movapd nb103_dyH1(%esp),%xmm1
        movapd nb103_dzH1(%esp),%xmm2
        movapd %xmm6,nb103_vctot(%esp)
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb103_fixH1(%esp),%xmm3
        movapd nb103_fiyH1(%esp),%xmm4
        movapd nb103_fizH1(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb103_fixH1(%esp)
        movapd %xmm4,nb103_fiyH1(%esp)
        movapd %xmm7,nb103_fizH1(%esp)
        ## update j forces with water H1 
        addpd  nb103_fjx(%esp),%xmm0
        addpd  nb103_fjy(%esp),%xmm1
        addpd  nb103_fjz(%esp),%xmm2
        movapd %xmm0,nb103_fjx(%esp)
        movapd %xmm1,nb103_fjy(%esp)
        movapd %xmm2,nb103_fjz(%esp)

        ## H2 interactions 
        movapd  %xmm5,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        mulpd  nb103_qqH(%esp),%xmm5    ## xmm5=vcoul 
        mulpd  %xmm5,%xmm4              ## total fsH1 in xmm4 

        addpd  nb103_vctot(%esp),%xmm5

        movapd nb103_dxH2(%esp),%xmm0
        movapd nb103_dyH2(%esp),%xmm1
        movapd nb103_dzH2(%esp),%xmm2
        movapd %xmm5,nb103_vctot(%esp)
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        ## update H2 forces 
        movapd nb103_fixH2(%esp),%xmm3
        movapd nb103_fiyH2(%esp),%xmm4
        movapd nb103_fizH2(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb103_fixH2(%esp)
        movapd %xmm4,nb103_fiyH2(%esp)
        movapd %xmm7,nb103_fizH2(%esp)

        movl nb103_faction(%ebp),%edi
        ## update j forces 
        addpd  nb103_fjx(%esp),%xmm0
        addpd  nb103_fjy(%esp),%xmm1
        addpd  nb103_fjz(%esp),%xmm2

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
        subl $2,nb103_innerk(%esp)
        jl    _nb_kernel103_ia32_sse2.nb103_checksingle
        jmp   _nb_kernel103_ia32_sse2.nb103_unroll_loop
_nb_kernel103_ia32_sse2.nb103_checksingle:      
        movl  nb103_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel103_ia32_sse2.nb103_dosingle
        jmp    _nb_kernel103_ia32_sse2.nb103_updateouterdata
_nb_kernel103_ia32_sse2.nb103_dosingle: 
        movl  nb103_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movl nb103_charge(%ebp),%esi     ## base of charge[] 
        xorpd %xmm6,%xmm6
        movlpd (%esi,%eax,8),%xmm6      ## jq A 

        movapd nb103_iqM(%esp),%xmm3
        movapd nb103_iqH(%esp),%xmm4
        mulsd %xmm6,%xmm3               ## qqM
        mulsd %xmm6,%xmm4               ## qqH 

        movapd  %xmm3,nb103_qqM(%esp)
        movapd  %xmm4,nb103_qqH(%esp)

        movl nb103_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2        
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ixM-izM to xmm4-xmm6 
        movapd nb103_ixM(%esp),%xmm4
        movapd nb103_iyM(%esp),%xmm5
        movapd nb103_izM(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb103_dxM(%esp)
        movapd %xmm5,nb103_dyM(%esp)
        movapd %xmm6,nb103_dzM(%esp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqM in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb103_ixH1(%esp),%xmm4
        movapd nb103_iyH1(%esp),%xmm5
        movapd nb103_izH1(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb103_dxH1(%esp)
        movapd %xmm5,nb103_dyH1(%esp)
        movapd %xmm6,nb103_dzH1(%esp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm6
        addsd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb103_ixH2(%esp),%xmm3
        movapd nb103_iyH2(%esp),%xmm4
        movapd nb103_izH2(%esp),%xmm5

        ## calc dr 
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5

        ## store dr 
        movapd %xmm3,nb103_dxH2(%esp)
        movapd %xmm4,nb103_dyH2(%esp)
        movapd %xmm5,nb103_dzH2(%esp)
        ## square it 
        mulsd %xmm3,%xmm3
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        addsd %xmm4,%xmm5
        addsd %xmm3,%xmm5
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## start with rsqM - put seed in xmm2 
        cvtsd2ss %xmm7,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb103_three(%esp),%xmm4
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb103_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb103_three(%esp),%xmm4
        subsd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb103_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvM in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb103_three(%esp),%xmm4
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb103_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb103_three(%esp),%xmm4
        subsd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb103_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb103_three(%esp),%xmm4
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb103_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb103_three(%esp),%xmm4
        subsd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb103_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do M interactions 
        movapd  %xmm7,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq 
        mulsd  nb103_qqM(%esp),%xmm7    ## xmm7=vcoul 

        mulsd  %xmm7,%xmm4      ## total fsM in xmm4 

        addsd  nb103_vctot(%esp),%xmm7

        movlpd %xmm7,nb103_vctot(%esp)

        movapd nb103_dxM(%esp),%xmm0
        movapd nb103_dyM(%esp),%xmm1
        movapd nb103_dzM(%esp),%xmm2
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update M forces 
        movapd nb103_fixM(%esp),%xmm3
        movapd nb103_fiyM(%esp),%xmm4
        movapd nb103_fizM(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb103_fixM(%esp)
        movlpd %xmm4,nb103_fiyM(%esp)
        movlpd %xmm7,nb103_fizM(%esp)
        ## update j forces with water M 
        movlpd %xmm0,nb103_fjx(%esp)
        movlpd %xmm1,nb103_fjy(%esp)
        movlpd %xmm2,nb103_fjz(%esp)

        ## H1 interactions 
        movapd  %xmm6,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        mulsd  nb103_qqH(%esp),%xmm6    ## xmm6=vcoul 
        mulsd  %xmm6,%xmm4              ## total fsH1 in xmm4 

        addsd  nb103_vctot(%esp),%xmm6

        movapd nb103_dxH1(%esp),%xmm0
        movapd nb103_dyH1(%esp),%xmm1
        movapd nb103_dzH1(%esp),%xmm2
        movlpd %xmm6,nb103_vctot(%esp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb103_fixH1(%esp),%xmm3
        movapd nb103_fiyH1(%esp),%xmm4
        movapd nb103_fizH1(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb103_fixH1(%esp)
        movlpd %xmm4,nb103_fiyH1(%esp)
        movlpd %xmm7,nb103_fizH1(%esp)
        ## update j forces with water H1 
        addsd  nb103_fjx(%esp),%xmm0
        addsd  nb103_fjy(%esp),%xmm1
        addsd  nb103_fjz(%esp),%xmm2
        movsd %xmm0,nb103_fjx(%esp)
        movsd %xmm1,nb103_fjy(%esp)
        movsd %xmm2,nb103_fjz(%esp)

        ## H2 interactions 
        movapd  %xmm5,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        mulsd  nb103_qqH(%esp),%xmm5    ## xmm5=vcoul 
        mulsd  %xmm5,%xmm4              ## total fsH1 in xmm4 

        addsd  nb103_vctot(%esp),%xmm5

        movapd nb103_dxH2(%esp),%xmm0
        movapd nb103_dyH2(%esp),%xmm1
        movapd nb103_dzH2(%esp),%xmm2
        movlpd %xmm5,nb103_vctot(%esp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H2 forces 
        movapd nb103_fixH2(%esp),%xmm3
        movapd nb103_fiyH2(%esp),%xmm4
        movapd nb103_fizH2(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb103_fixH2(%esp)
        movlpd %xmm4,nb103_fiyH2(%esp)
        movlpd %xmm7,nb103_fizH2(%esp)

        movl nb103_faction(%ebp),%edi
        ## update j forces 
        addsd  nb103_fjx(%esp),%xmm0
        addsd  nb103_fjy(%esp),%xmm1
        addsd  nb103_fjz(%esp),%xmm2

        movlpd (%edi,%eax,8),%xmm3
        movlpd 8(%edi,%eax,8),%xmm4
        movlpd 16(%edi,%eax,8),%xmm5
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        movlpd %xmm3,(%edi,%eax,8)
        movlpd %xmm4,8(%edi,%eax,8)
        movlpd %xmm5,16(%edi,%eax,8)

_nb_kernel103_ia32_sse2.nb103_updateouterdata: 
        movl  nb103_ii3(%esp),%ecx
        movl  nb103_faction(%ebp),%edi
        movl  nb103_fshift(%ebp),%esi
        movl  nb103_is3(%esp),%edx

        ## accumulate H1i forces in xmm0, xmm1, xmm2 
        movapd nb103_fixH1(%esp),%xmm0
        movapd nb103_fiyH1(%esp),%xmm1
        movapd nb103_fizH1(%esp),%xmm2

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
        movapd %xmm0,%xmm6
        movsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm6

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movapd nb103_fixH2(%esp),%xmm0
        movapd nb103_fiyH2(%esp),%xmm1
        movapd nb103_fizH2(%esp),%xmm2

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

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movapd nb103_fixM(%esp),%xmm0
        movapd nb103_fiyM(%esp),%xmm1
        movapd nb103_fizM(%esp),%xmm2

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
        movl nb103_n(%esp),%esi
        ## get group index for i particle 
        movl  nb103_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb103_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 

        ## add earlier value from mem 
        movl  nb103_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb103_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel103_ia32_sse2.nb103_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb103_n(%esp)
        jmp _nb_kernel103_ia32_sse2.nb103_outer
_nb_kernel103_ia32_sse2.nb103_outerend: 
        ## check if more outer neighborlists remain
        movl  nb103_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel103_ia32_sse2.nb103_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel103_ia32_sse2.nb103_threadloop
_nb_kernel103_ia32_sse2.nb103_end: 
        emms

        movl nb103_nouter(%esp),%eax
        movl nb103_ninner(%esp),%ebx
        movl nb103_outeriter(%ebp),%ecx
        movl nb103_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb103_salign(%esp),%eax
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



.globl nb_kernel103nf_ia32_sse2
.globl _nb_kernel103nf_ia32_sse2
nb_kernel103nf_ia32_sse2:       
_nb_kernel103nf_ia32_sse2:      
.set nb103nf_p_nri, 8
.set nb103nf_iinr, 12
.set nb103nf_jindex, 16
.set nb103nf_jjnr, 20
.set nb103nf_shift, 24
.set nb103nf_shiftvec, 28
.set nb103nf_fshift, 32
.set nb103nf_gid, 36
.set nb103nf_pos, 40
.set nb103nf_faction, 44
.set nb103nf_charge, 48
.set nb103nf_p_facel, 52
.set nb103nf_argkrf, 56
.set nb103nf_argcrf, 60
.set nb103nf_Vc, 64
.set nb103nf_type, 68
.set nb103nf_p_ntype, 72
.set nb103nf_vdwparam, 76
.set nb103nf_Vvdw, 80
.set nb103nf_p_tabscale, 84
.set nb103nf_VFtab, 88
.set nb103nf_invsqrta, 92
.set nb103nf_dvda, 96
.set nb103nf_p_gbtabscale, 100
.set nb103nf_GBtab, 104
.set nb103nf_p_nthreads, 108
.set nb103nf_count, 112
.set nb103nf_mtx, 116
.set nb103nf_outeriter, 120
.set nb103nf_inneriter, 124
.set nb103nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb103nf_ixM, 0
.set nb103nf_iyM, 16
.set nb103nf_izM, 32
.set nb103nf_ixH1, 48
.set nb103nf_iyH1, 64
.set nb103nf_izH1, 80
.set nb103nf_ixH2, 96
.set nb103nf_iyH2, 112
.set nb103nf_izH2, 128
.set nb103nf_iqM, 144
.set nb103nf_iqH, 160
.set nb103nf_qqM, 176
.set nb103nf_qqH, 192
.set nb103nf_vctot, 208
.set nb103nf_half, 224
.set nb103nf_three, 240
.set nb103nf_is3, 256
.set nb103nf_ii3, 260
.set nb103nf_innerjjnr, 264
.set nb103nf_innerk, 268
.set nb103nf_n, 272
.set nb103nf_nn1, 276
.set nb103nf_nri, 280
.set nb103nf_nouter, 284
.set nb103nf_ninner, 288
.set nb103nf_salign, 292
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
        movl %eax,nb103nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb103nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb103nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb103nf_nouter(%esp)
        movl %eax,nb103nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb103nf_half(%esp)
        movl %ebx,nb103nf_half+4(%esp)
        movsd nb103nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb103nf_half(%esp)
        movapd %xmm3,nb103nf_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb103nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb103nf_charge(%ebp),%edx
        movsd 8(%edx,%ebx,8),%xmm3
        movsd 24(%edx,%ebx,8),%xmm4
        movl nb103nf_p_facel(%ebp),%esi
        movsd (%esi),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb103nf_iqH(%esp)
        movapd %xmm4,nb103nf_iqM(%esp)

_nb_kernel103nf_ia32_sse2.nb103nf_threadloop: 
        movl  nb103nf_count(%ebp),%esi          ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel103nf_ia32_sse2.nb103nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel103nf_ia32_sse2.nb103nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb103nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb103nf_n(%esp)
        movl %ebx,nb103nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel103nf_ia32_sse2.nb103nf_outerstart
        jmp _nb_kernel103nf_ia32_sse2.nb103nf_end

_nb_kernel103nf_ia32_sse2.nb103nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb103nf_nouter(%esp),%ebx
        movl %ebx,nb103nf_nouter(%esp)

_nb_kernel103nf_ia32_sse2.nb103nf_outer: 
        movl  nb103nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 

        movl  nb103nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb103nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb103nf_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb103nf_ii3(%esp)

        addsd 24(%eax,%ebx,8),%xmm3
        addsd 32(%eax,%ebx,8),%xmm4
        addsd 40(%eax,%ebx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb103nf_ixH1(%esp)
        movapd %xmm4,nb103nf_iyH1(%esp)
        movapd %xmm5,nb103nf_izH1(%esp)

        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 48(%eax,%ebx,8),%xmm0
        addsd 56(%eax,%ebx,8),%xmm1
        addsd 64(%eax,%ebx,8),%xmm2
        addsd 72(%eax,%ebx,8),%xmm3
        addsd 80(%eax,%ebx,8),%xmm4
        addsd 88(%eax,%ebx,8),%xmm5

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm0,nb103nf_ixH2(%esp)
        movapd %xmm1,nb103nf_iyH2(%esp)
        movapd %xmm2,nb103nf_izH2(%esp)
        movapd %xmm3,nb103nf_ixM(%esp)
        movapd %xmm4,nb103nf_iyM(%esp)
        movapd %xmm5,nb103nf_izM(%esp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb103nf_vctot(%esp)

        movl  nb103nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb103nf_pos(%ebp),%esi
        movl  nb103nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb103nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb103nf_ninner(%esp),%ecx
        movl  %ecx,nb103nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb103nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel103nf_ia32_sse2.nb103nf_unroll_loop
        jmp   _nb_kernel103nf_ia32_sse2.nb103nf_checksingle
_nb_kernel103nf_ia32_sse2.nb103nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb103nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb103nf_innerjjnr(%esp)             ## advance pointer (unrolled 2) 

        movl nb103nf_charge(%ebp),%esi     ## base of charge[] 


        movlpd (%esi,%eax,8),%xmm6      ## jq A 
        movhpd (%esi,%ebx,8),%xmm6      ## jq B 
        movapd nb103nf_iqM(%esp),%xmm3
        movapd nb103nf_iqH(%esp),%xmm4
        mulpd %xmm6,%xmm3               ## qqM 
        mulpd %xmm6,%xmm4               ## qqH 

        movapd  %xmm3,nb103nf_qqM(%esp)
        movapd  %xmm4,nb103nf_qqH(%esp)

        movl nb103nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move two coordinates to xmm0-xmm2    
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2
        movhpd (%esi,%ebx,8),%xmm0
        movhpd 8(%esi,%ebx,8),%xmm1
        movhpd 16(%esi,%ebx,8),%xmm2

        ## move ixM-izM to xmm4-xmm6 
        movapd nb103nf_ixM(%esp),%xmm4
        movapd nb103nf_iyM(%esp),%xmm5
        movapd nb103nf_izM(%esp),%xmm6

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
        ## rsqM in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb103nf_ixH1(%esp),%xmm4
        movapd nb103nf_iyH1(%esp),%xmm5
        movapd nb103nf_izH1(%esp),%xmm6

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
        movapd nb103nf_ixH2(%esp),%xmm3
        movapd nb103nf_iyH2(%esp),%xmm4
        movapd nb103nf_izH2(%esp),%xmm5

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
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqM in xmm7 

        ## start with rsqM - put seed in xmm2 
        cvtpd2ps %xmm7,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb103nf_three(%esp),%xmm4
        mulpd   %xmm7,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb103nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb103nf_three(%esp),%xmm4
        subpd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb103nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvM in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb103nf_three(%esp),%xmm4
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb103nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb103nf_three(%esp),%xmm4
        subpd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb103nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb103nf_three(%esp),%xmm4
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb103nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb103nf_three(%esp),%xmm4
        subpd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb103nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do M interactions 
        mulpd  nb103nf_qqM(%esp),%xmm7          ## xmm7=vcoul 
        addpd  nb103nf_vctot(%esp),%xmm7
        movapd %xmm7,nb103nf_vctot(%esp)

        ## H1 interactions 
        mulpd  nb103nf_qqH(%esp),%xmm6          ## xmm6=vcoul 
        addpd  nb103nf_vctot(%esp),%xmm6
        movapd %xmm6,nb103nf_vctot(%esp)

        ## H2 interactions 
        mulpd  nb103nf_qqH(%esp),%xmm5          ## xmm5=vcoul 
        addpd  nb103nf_vctot(%esp),%xmm5
        movapd %xmm5,nb103nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $2,nb103nf_innerk(%esp)
        jl    _nb_kernel103nf_ia32_sse2.nb103nf_checksingle
        jmp   _nb_kernel103nf_ia32_sse2.nb103nf_unroll_loop
_nb_kernel103nf_ia32_sse2.nb103nf_checksingle:  
        movl  nb103nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel103nf_ia32_sse2.nb103nf_dosingle
        jmp   _nb_kernel103nf_ia32_sse2.nb103nf_updateouterdata
_nb_kernel103nf_ia32_sse2.nb103nf_dosingle: 
        movl  nb103nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movl nb103nf_charge(%ebp),%esi     ## base of charge[] 
        xorpd %xmm6,%xmm6
        movlpd (%esi,%eax,8),%xmm6      ## jq A 

        movapd nb103nf_iqM(%esp),%xmm3
        movapd nb103nf_iqH(%esp),%xmm4
        mulsd %xmm6,%xmm3               ## qqM 
        mulsd %xmm6,%xmm4               ## qqH 

        movapd  %xmm3,nb103nf_qqM(%esp)
        movapd  %xmm4,nb103nf_qqH(%esp)

        movl nb103nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2        
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ixM-izM to xmm4-xmm6 
        movapd nb103nf_ixM(%esp),%xmm4
        movapd nb103nf_iyM(%esp),%xmm5
        movapd nb103nf_izM(%esp),%xmm6

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
        ## rsqM in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb103nf_ixH1(%esp),%xmm4
        movapd nb103nf_iyH1(%esp),%xmm5
        movapd nb103nf_izH1(%esp),%xmm6

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
        movapd nb103nf_ixH2(%esp),%xmm3
        movapd nb103nf_iyH2(%esp),%xmm4
        movapd nb103nf_izH2(%esp),%xmm5

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
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqM in xmm7 

        ## start with rsqM - put seed in xmm2 
        cvtsd2ss %xmm7,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb103nf_three(%esp),%xmm4
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb103nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb103nf_three(%esp),%xmm4
        subsd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb103nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvM in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb103nf_three(%esp),%xmm4
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb103nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb103nf_three(%esp),%xmm4
        subsd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb103nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb103nf_three(%esp),%xmm4
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb103nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb103nf_three(%esp),%xmm4
        subsd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb103nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do M interactions 
        mulsd  nb103nf_qqM(%esp),%xmm7          ## xmm7=vcoul 
        addsd  nb103nf_vctot(%esp),%xmm7
        movlpd %xmm7,nb103nf_vctot(%esp)

        ## H1 interactions 
        mulsd  nb103nf_qqH(%esp),%xmm6          ## xmm6=vcoul 
        addsd  nb103nf_vctot(%esp),%xmm6
        movlpd %xmm6,nb103nf_vctot(%esp)

        ## H2 interactions 
        mulsd  nb103nf_qqH(%esp),%xmm5          ## xmm5=vcoul 
        addsd  nb103nf_vctot(%esp),%xmm5
        movlpd %xmm5,nb103nf_vctot(%esp)

_nb_kernel103nf_ia32_sse2.nb103nf_updateouterdata: 
        ## get n from stack
        movl nb103nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb103nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        movapd nb103nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 

        ## add earlier value from mem 
        movl  nb103nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb103nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel103nf_ia32_sse2.nb103nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb103nf_n(%esp)
        jmp _nb_kernel103nf_ia32_sse2.nb103nf_outer
_nb_kernel103nf_ia32_sse2.nb103nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb103nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel103nf_ia32_sse2.nb103nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel103nf_ia32_sse2.nb103nf_threadloop
_nb_kernel103nf_ia32_sse2.nb103nf_end: 
        emms

        movl nb103nf_nouter(%esp),%eax
        movl nb103nf_ninner(%esp),%ebx
        movl nb103nf_outeriter(%ebp),%ecx
        movl nb103nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb103nf_salign(%esp),%eax
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


