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



.globl nb_kernel311_ia32_sse2
.globl _nb_kernel311_ia32_sse2
nb_kernel311_ia32_sse2: 
_nb_kernel311_ia32_sse2:        
.set nb311_p_nri, 8
.set nb311_iinr, 12
.set nb311_jindex, 16
.set nb311_jjnr, 20
.set nb311_shift, 24
.set nb311_shiftvec, 28
.set nb311_fshift, 32
.set nb311_gid, 36
.set nb311_pos, 40
.set nb311_faction, 44
.set nb311_charge, 48
.set nb311_p_facel, 52
.set nb311_argkrf, 56
.set nb311_argcrf, 60
.set nb311_Vc, 64
.set nb311_type, 68
.set nb311_p_ntype, 72
.set nb311_vdwparam, 76
.set nb311_Vvdw, 80
.set nb311_p_tabscale, 84
.set nb311_VFtab, 88
.set nb311_invsqrta, 92
.set nb311_dvda, 96
.set nb311_p_gbtabscale, 100
.set nb311_GBtab, 104
.set nb311_p_nthreads, 108
.set nb311_count, 112
.set nb311_mtx, 116
.set nb311_outeriter, 120
.set nb311_inneriter, 124
.set nb311_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb311_ixO, 0
.set nb311_iyO, 16
.set nb311_izO, 32
.set nb311_ixH1, 48
.set nb311_iyH1, 64
.set nb311_izH1, 80
.set nb311_ixH2, 96
.set nb311_iyH2, 112
.set nb311_izH2, 128
.set nb311_iqO, 144
.set nb311_iqH, 160
.set nb311_dxO, 176
.set nb311_dyO, 192
.set nb311_dzO, 208
.set nb311_dxH1, 224
.set nb311_dyH1, 240
.set nb311_dzH1, 256
.set nb311_dxH2, 272
.set nb311_dyH2, 288
.set nb311_dzH2, 304
.set nb311_qqO, 320
.set nb311_qqH, 336
.set nb311_rinvO, 352
.set nb311_rinvH1, 368
.set nb311_rinvH2, 384
.set nb311_rO, 400
.set nb311_rH1, 416
.set nb311_rH2, 432
.set nb311_tsc, 448
.set nb311_two, 464
.set nb311_c6, 480
.set nb311_c12, 496
.set nb311_six, 512
.set nb311_twelve, 528
.set nb311_vctot, 544
.set nb311_Vvdwtot, 560
.set nb311_fixO, 576
.set nb311_fiyO, 592
.set nb311_fizO, 608
.set nb311_fixH1, 624
.set nb311_fiyH1, 640
.set nb311_fizH1, 656
.set nb311_fixH2, 672
.set nb311_fiyH2, 688
.set nb311_fizH2, 704
.set nb311_fjx, 720
.set nb311_fjy, 736
.set nb311_fjz, 752
.set nb311_half, 768
.set nb311_three, 784
.set nb311_is3, 800
.set nb311_ii3, 804
.set nb311_ntia, 808
.set nb311_innerjjnr, 812
.set nb311_innerk, 816
.set nb311_n, 820
.set nb311_nn1, 824
.set nb311_nri, 828
.set nb311_nouter, 832
.set nb311_ninner, 836
.set nb311_salign, 840
        pushl %ebp
        movl %esp,%ebp
    pushl %eax
    pushl %ebx
    pushl %ecx
    pushl %edx
        pushl %esi
        pushl %edi
        subl $844,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb311_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb311_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb311_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb311_nouter(%esp)
        movl %eax,nb311_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb311_half(%esp)
        movl %ebx,nb311_half+4(%esp)
        movsd nb311_half(%esp),%xmm1
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
        movapd %xmm1,nb311_half(%esp)
        movapd %xmm2,nb311_two(%esp)
        movapd %xmm3,nb311_three(%esp)
        movapd %xmm4,nb311_six(%esp)
        movapd %xmm5,nb311_twelve(%esp)

        movl nb311_p_tabscale(%ebp),%eax
        movsd (%eax),%xmm5
        shufpd $0,%xmm5,%xmm5
        movapd %xmm5,nb311_tsc(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb311_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb311_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        movsd 8(%edx,%ebx,8),%xmm4
        movl nb311_p_facel(%ebp),%esi
        movsd (%esi),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb311_iqO(%esp)
        movapd %xmm4,nb311_iqH(%esp)

        movl  nb311_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb311_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb311_ntia(%esp)
_nb_kernel311_ia32_sse2.nb311_threadloop: 
        movl  nb311_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel311_ia32_sse2.nb311_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel311_ia32_sse2.nb311_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb311_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb311_n(%esp)
        movl %ebx,nb311_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel311_ia32_sse2.nb311_outerstart
        jmp _nb_kernel311_ia32_sse2.nb311_end

_nb_kernel311_ia32_sse2.nb311_outerstart: 
        ## ebx contains number of outer iterations
        addl nb311_nouter(%esp),%ebx
        movl %ebx,nb311_nouter(%esp)

_nb_kernel311_ia32_sse2.nb311_outer: 
        movl  nb311_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb311_is3(%esp)      ## store is3 

        movl  nb311_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb311_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb311_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb311_ii3(%esp)

        addsd (%eax,%ebx,8),%xmm3
        addsd 8(%eax,%ebx,8),%xmm4
        addsd 16(%eax,%ebx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb311_ixO(%esp)
        movapd %xmm4,nb311_iyO(%esp)
        movapd %xmm5,nb311_izO(%esp)

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
        movapd %xmm0,nb311_ixH1(%esp)
        movapd %xmm1,nb311_iyH1(%esp)
        movapd %xmm2,nb311_izH1(%esp)
        movapd %xmm3,nb311_ixH2(%esp)
        movapd %xmm4,nb311_iyH2(%esp)
        movapd %xmm5,nb311_izH2(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb311_vctot(%esp)
        movapd %xmm4,nb311_Vvdwtot(%esp)
        movapd %xmm4,nb311_fixO(%esp)
        movapd %xmm4,nb311_fiyO(%esp)
        movapd %xmm4,nb311_fizO(%esp)
        movapd %xmm4,nb311_fixH1(%esp)
        movapd %xmm4,nb311_fiyH1(%esp)
        movapd %xmm4,nb311_fizH1(%esp)
        movapd %xmm4,nb311_fixH2(%esp)
        movapd %xmm4,nb311_fiyH2(%esp)
        movapd %xmm4,nb311_fizH2(%esp)

        movl  nb311_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb311_pos(%ebp),%esi
        movl  nb311_faction(%ebp),%edi
        movl  nb311_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb311_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb311_ninner(%esp),%ecx
        movl  %ecx,nb311_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb311_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel311_ia32_sse2.nb311_unroll_loop
        jmp   _nb_kernel311_ia32_sse2.nb311_checksingle
_nb_kernel311_ia32_sse2.nb311_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb311_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb311_innerjjnr(%esp)             ## advance pointer (unrolled 2) 

        movl nb311_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb311_iqO(%esp),%xmm3
        mulpd  nb311_iqH(%esp),%xmm4
        movapd  %xmm3,nb311_qqO(%esp)
        movapd  %xmm4,nb311_qqH(%esp)

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movl nb311_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb311_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb311_ntia(%esp),%edi
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
        movapd %xmm4,nb311_c6(%esp)
        movapd %xmm6,nb311_c12(%esp)

        movl nb311_pos(%ebp),%esi        ## base of pos[] 

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
        movapd nb311_ixO(%esp),%xmm4
        movapd nb311_iyO(%esp),%xmm5
        movapd nb311_izO(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb311_dxO(%esp)
        movapd %xmm5,nb311_dyO(%esp)
        movapd %xmm6,nb311_dzO(%esp)
        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb311_ixH1(%esp),%xmm4
        movapd nb311_iyH1(%esp),%xmm5
        movapd nb311_izH1(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb311_dxH1(%esp)
        movapd %xmm5,nb311_dyH1(%esp)
        movapd %xmm6,nb311_dzH1(%esp)
        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm6
        addpd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb311_ixH2(%esp),%xmm3
        movapd nb311_iyH2(%esp),%xmm4
        movapd nb311_izH2(%esp),%xmm5

        ## calc dr 
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5

        ## store dr 
        movapd %xmm3,nb311_dxH2(%esp)
        movapd %xmm4,nb311_dyH2(%esp)
        movapd %xmm5,nb311_dzH2(%esp)
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
        movapd  nb311_three(%esp),%xmm4
        mulpd   %xmm7,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb311_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm7,%xmm2
        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb311_three(%esp),%xmm4
        subpd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb311_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,nb311_rinvO(%esp)         ## rinvO in xmm4 
        mulpd   %xmm4,%xmm7
        movapd  %xmm7,nb311_rO(%esp)    ## r in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb311_three(%esp),%xmm4
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb311_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm6,%xmm2
        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb311_three(%esp),%xmm4
        subpd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb311_half(%esp),%xmm4   ## rinv 
        movapd %xmm4,nb311_rinvH1(%esp)         ## rinvH1 
        mulpd  %xmm4,%xmm6
        movapd %xmm6,nb311_rH1(%esp)    ## rH1 

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb311_three(%esp),%xmm4
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb311_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm5,%xmm2
        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb311_three(%esp),%xmm4
        subpd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb311_half(%esp),%xmm4   ## rinv 
        movapd %xmm4,nb311_rinvH2(%esp)   ## rinv 
        mulpd %xmm4,%xmm5
        movapd %xmm5,nb311_rH2(%esp)   ## r 

        ## do O interactions 
        ## rO is still in xmm7 
        movapd nb311_rinvO(%esp),%xmm0
        mulpd   nb311_tsc(%esp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movd %eax,%mm0
        movd %ebx,%mm1
        movl nb311_VFtab(%ebp),%esi
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
    mulpd  nb311_two(%esp),%xmm7         ## two*Heps2 
    movapd nb311_qqO(%esp),%xmm0
    addpd  %xmm6,%xmm7
    addpd  %xmm5,%xmm7 ## xmm7=FF 
    mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addpd  %xmm4,%xmm5 ## xmm5=VV 
    mulpd  %xmm0,%xmm5 ## vcoul=qq*VV  
    mulpd  %xmm7,%xmm0 ## fijC=FF*qq 

        ## do nontable L-J 
        movapd nb311_rinvO(%esp),%xmm2
        mulpd  %xmm2,%xmm2

    ## at this point mm5 contains vcoul and xmm0 fijC 
    ## increment vcoul - then we can get rid of mm5 
    addpd  nb311_vctot(%esp),%xmm5
    movapd %xmm5,nb311_vctot(%esp)

        movapd %xmm2,%xmm1
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm4
        mulpd  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulpd  nb311_c6(%esp),%xmm1
        mulpd  nb311_c12(%esp),%xmm4
        movapd %xmm4,%xmm3
        subpd  %xmm1,%xmm3      ## xmm3=Vvdw12-Vvdw6 
        mulpd  nb311_six(%esp),%xmm1
        mulpd  nb311_twelve(%esp),%xmm4
        subpd  %xmm1,%xmm4
        addpd  nb311_Vvdwtot(%esp),%xmm3
        mulpd  nb311_rinvO(%esp),%xmm4
        mulpd  nb311_tsc(%esp),%xmm0
        subpd  %xmm0,%xmm4
        movapd %xmm3,nb311_Vvdwtot(%esp)
        mulpd  nb311_rinvO(%esp),%xmm4

        movapd nb311_dxO(%esp),%xmm0
        movapd nb311_dyO(%esp),%xmm1
        movapd nb311_dzO(%esp),%xmm2
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2      ## tx in xmm0-xmm2 

        ## update O forces 
        movapd nb311_fixO(%esp),%xmm3
        movapd nb311_fiyO(%esp),%xmm4
        movapd nb311_fizO(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb311_fixO(%esp)
        movapd %xmm4,nb311_fiyO(%esp)
        movapd %xmm7,nb311_fizO(%esp)
        ## update j forces with water O 
        movapd %xmm0,nb311_fjx(%esp)
        movapd %xmm1,nb311_fjy(%esp)
        movapd %xmm2,nb311_fjz(%esp)

        ## Done with O interactions - now H1! 
        movapd nb311_rH1(%esp),%xmm7
        mulpd nb311_tsc(%esp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb311_VFtab(%ebp),%esi
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
        mulpd  nb311_two(%esp),%xmm7    ## two*Heps2 
        movapd nb311_qqH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 
    ## increment vcoul 
        xorpd  %xmm4,%xmm4
    addpd  nb311_vctot(%esp),%xmm5
        mulpd  nb311_rinvH1(%esp),%xmm3
    movapd %xmm5,nb311_vctot(%esp)
        mulpd  nb311_tsc(%esp),%xmm3
        subpd %xmm3,%xmm4

        movapd nb311_dxH1(%esp),%xmm0
        movapd nb311_dyH1(%esp),%xmm1
        movapd nb311_dzH1(%esp),%xmm2
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb311_fixH1(%esp),%xmm3
        movapd nb311_fiyH1(%esp),%xmm4
        movapd nb311_fizH1(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb311_fixH1(%esp)
        movapd %xmm4,nb311_fiyH1(%esp)
        movapd %xmm7,nb311_fizH1(%esp)
        ## update j forces with water H1 
        addpd  nb311_fjx(%esp),%xmm0
        addpd  nb311_fjy(%esp),%xmm1
        addpd  nb311_fjz(%esp),%xmm2
        movapd %xmm0,nb311_fjx(%esp)
        movapd %xmm1,nb311_fjy(%esp)
        movapd %xmm2,nb311_fjz(%esp)

        ## Done with H1, finally we do H2 interactions 
        movapd nb311_rH2(%esp),%xmm7
        mulpd   nb311_tsc(%esp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb311_VFtab(%ebp),%esi
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
        mulpd  nb311_two(%esp),%xmm7    ## two*Heps2 
        movapd nb311_qqH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 
    ## increment vcoul 
        xorpd  %xmm4,%xmm4
    addpd  nb311_vctot(%esp),%xmm5
        mulpd  nb311_rinvH2(%esp),%xmm3
    movapd %xmm5,nb311_vctot(%esp)
        mulpd  nb311_tsc(%esp),%xmm3
        subpd  %xmm3,%xmm4

        movapd nb311_dxH2(%esp),%xmm0
        movapd nb311_dyH2(%esp),%xmm1
        movapd nb311_dzH2(%esp),%xmm2
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

    movd %mm0,%eax
    movd %mm1,%ebx

        ## update H2 forces 
        movapd nb311_fixH2(%esp),%xmm3
        movapd nb311_fiyH2(%esp),%xmm4
        movapd nb311_fizH2(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb311_fixH2(%esp)
        movapd %xmm4,nb311_fiyH2(%esp)
        movapd %xmm7,nb311_fizH2(%esp)

        movl nb311_faction(%ebp),%edi
        ## update j forces 
        addpd  nb311_fjx(%esp),%xmm0
        addpd  nb311_fjy(%esp),%xmm1
        addpd  nb311_fjz(%esp),%xmm2

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
        subl $2,nb311_innerk(%esp)
        jl    _nb_kernel311_ia32_sse2.nb311_checksingle
        jmp   _nb_kernel311_ia32_sse2.nb311_unroll_loop
_nb_kernel311_ia32_sse2.nb311_checksingle: 
        movl  nb311_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel311_ia32_sse2.nb311_dosingle
        jmp   _nb_kernel311_ia32_sse2.nb311_updateouterdata
_nb_kernel311_ia32_sse2.nb311_dosingle: 
        movl  nb311_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movl nb311_charge(%ebp),%esi     ## base of charge[] 
        xorpd %xmm3,%xmm3
        movlpd (%esi,%eax,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb311_iqO(%esp),%xmm3
        mulpd  nb311_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movapd  %xmm3,nb311_qqO(%esp)
        movapd  %xmm4,nb311_qqH(%esp)

        movl nb311_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb311_vdwparam(%ebp),%esi
        shll %eax
        movl nb311_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movapd %xmm4,nb311_c6(%esp)
        movapd %xmm6,nb311_c12(%esp)

        movl nb311_pos(%ebp),%esi        ## base of pos[] 
        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move coords to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb311_ixO(%esp),%xmm4
        movapd nb311_iyO(%esp),%xmm5
        movapd nb311_izO(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb311_dxO(%esp)
        movapd %xmm5,nb311_dyO(%esp)
        movapd %xmm6,nb311_dzO(%esp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb311_ixH1(%esp),%xmm4
        movapd nb311_iyH1(%esp),%xmm5
        movapd nb311_izH1(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb311_dxH1(%esp)
        movapd %xmm5,nb311_dyH1(%esp)
        movapd %xmm6,nb311_dzH1(%esp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm6
        addsd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb311_ixH2(%esp),%xmm3
        movapd nb311_iyH2(%esp),%xmm4
        movapd nb311_izH2(%esp),%xmm5

        ## calc dr 
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5

        ## store dr 
        movapd %xmm3,nb311_dxH2(%esp)
        movapd %xmm4,nb311_dyH2(%esp)
        movapd %xmm5,nb311_dzH2(%esp)
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
        movapd  nb311_three(%esp),%xmm4
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb311_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm7,%xmm2
        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb311_three(%esp),%xmm4
        subsd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb311_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,nb311_rinvO(%esp)         ## rinvO in xmm4 
        mulsd   %xmm4,%xmm7
        movapd  %xmm7,nb311_rO(%esp)    ## r in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb311_three(%esp),%xmm4
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb311_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm6,%xmm2
        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb311_three(%esp),%xmm4
        subsd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb311_half(%esp),%xmm4   ## rinv 
        movapd %xmm4,nb311_rinvH1(%esp)         ## rinvH1 
        mulsd  %xmm4,%xmm6
        movapd %xmm6,nb311_rH1(%esp)    ## rH1 

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb311_three(%esp),%xmm4
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb311_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm5,%xmm2
        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb311_three(%esp),%xmm4
        subsd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb311_half(%esp),%xmm4   ## rinv 
        movapd %xmm4,nb311_rinvH2(%esp)   ## rinv 
        mulsd %xmm4,%xmm5
        movapd %xmm5,nb311_rH2(%esp)   ## r 

        ## do O interactions 
        movd %eax,%mm0
        ## rO is still in xmm7 
        movapd nb311_rinvO(%esp),%xmm0
        mulsd   nb311_tsc(%esp),%xmm7
        cvttsd2si %xmm7,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb311_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1  

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7  
    mulsd  %xmm1,%xmm6      ## xmm6=Geps 
    mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
    addsd  %xmm6,%xmm5
    addsd  %xmm7,%xmm5      ## xmm5=Fp        
    mulsd  nb311_two(%esp),%xmm7         ## two*Heps2 
    movapd nb311_qqO(%esp),%xmm0
    addsd  %xmm6,%xmm7
    addsd  %xmm5,%xmm7 ## xmm7=FF 
    mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addsd  %xmm4,%xmm5 ## xmm5=VV 
    mulsd  %xmm0,%xmm5 ## vcoul=qq*VV  
    mulsd  %xmm7,%xmm0 ## fijC=FF*qq 

        ## do nontable L-J 
        movapd nb311_rinvO(%esp),%xmm2
        mulsd  %xmm2,%xmm2

    ## at this point mm5 contains vcoul and xmm0 fijC 
    ## increment vcoul - then we can get rid of mm5 
    addsd  nb311_vctot(%esp),%xmm5
    movlpd %xmm5,nb311_vctot(%esp)

        movapd %xmm2,%xmm1
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm4
        mulsd  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulsd  nb311_c6(%esp),%xmm1
        mulsd  nb311_c12(%esp),%xmm4
        movapd %xmm4,%xmm3
        subsd  %xmm1,%xmm3      ## xmm3=Vvdw12-Vvdw6 
        mulsd  nb311_six(%esp),%xmm1
        mulsd  nb311_twelve(%esp),%xmm4
        subsd  %xmm1,%xmm4
        addsd  nb311_Vvdwtot(%esp),%xmm3
        mulsd  nb311_rinvO(%esp),%xmm4
        mulsd  nb311_tsc(%esp),%xmm0
        subsd  %xmm0,%xmm4
        movlpd %xmm3,nb311_Vvdwtot(%esp)
        mulsd  nb311_rinvO(%esp),%xmm4

        movapd nb311_dxO(%esp),%xmm0
        movapd nb311_dyO(%esp),%xmm1
        movapd nb311_dzO(%esp),%xmm2
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2      ## tx in xmm0-xmm2 

        ## update O forces 
        movapd nb311_fixO(%esp),%xmm3
        movapd nb311_fiyO(%esp),%xmm4
        movapd nb311_fizO(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb311_fixO(%esp)
        movlpd %xmm4,nb311_fiyO(%esp)
        movlpd %xmm7,nb311_fizO(%esp)
        ## update j forces with water O 
        movlpd %xmm0,nb311_fjx(%esp)
        movlpd %xmm1,nb311_fjy(%esp)
        movlpd %xmm2,nb311_fjz(%esp)

        ## Done with O interactions - now H1! 
        movapd nb311_rH1(%esp),%xmm7
        mulpd nb311_tsc(%esp),%xmm7
        cvttsd2si %xmm7,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb311_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1  

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb311_two(%esp),%xmm7    ## two*Heps2 
        movapd nb311_qqH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 
    ## increment vcoul 
        xorpd  %xmm4,%xmm4
    addsd  nb311_vctot(%esp),%xmm5
        mulsd  nb311_rinvH1(%esp),%xmm3
    movlpd %xmm5,nb311_vctot(%esp)
        mulsd  nb311_tsc(%esp),%xmm3
        subsd %xmm3,%xmm4

        movapd nb311_dxH1(%esp),%xmm0
        movapd nb311_dyH1(%esp),%xmm1
        movapd nb311_dzH1(%esp),%xmm2
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb311_fixH1(%esp),%xmm3
        movapd nb311_fiyH1(%esp),%xmm4
        movapd nb311_fizH1(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb311_fixH1(%esp)
        movlpd %xmm4,nb311_fiyH1(%esp)
        movlpd %xmm7,nb311_fizH1(%esp)
        ## update j forces with water H1 
        addsd  nb311_fjx(%esp),%xmm0
        addsd  nb311_fjy(%esp),%xmm1
        addsd  nb311_fjz(%esp),%xmm2
        movlpd %xmm0,nb311_fjx(%esp)
        movlpd %xmm1,nb311_fjy(%esp)
        movlpd %xmm2,nb311_fjz(%esp)

        ## Done with H1, finally we do H2 interactions 
        movapd nb311_rH2(%esp),%xmm7
        mulsd   nb311_tsc(%esp),%xmm7
        cvttsd2si %xmm7,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb311_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1  

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb311_two(%esp),%xmm7    ## two*Heps2 
        movapd nb311_qqH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 
    ## increment vcoul 
        xorpd  %xmm4,%xmm4
    addsd  nb311_vctot(%esp),%xmm5
        mulsd  nb311_rinvH2(%esp),%xmm3
    movlpd %xmm5,nb311_vctot(%esp)
        mulsd  nb311_tsc(%esp),%xmm3
        subsd  %xmm3,%xmm4

        movapd nb311_dxH2(%esp),%xmm0
        movapd nb311_dyH2(%esp),%xmm1
        movapd nb311_dzH2(%esp),%xmm2
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

    movd %mm0,%eax

        ## update H2 forces 
        movapd nb311_fixH2(%esp),%xmm3
        movapd nb311_fiyH2(%esp),%xmm4
        movapd nb311_fizH2(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb311_fixH2(%esp)
        movlpd %xmm4,nb311_fiyH2(%esp)
        movlpd %xmm7,nb311_fizH2(%esp)

        movl nb311_faction(%ebp),%edi
        ## update j forces 
        addsd  nb311_fjx(%esp),%xmm0
        addsd  nb311_fjy(%esp),%xmm1
        addsd  nb311_fjz(%esp),%xmm2

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

_nb_kernel311_ia32_sse2.nb311_updateouterdata: 
        movl  nb311_ii3(%esp),%ecx
        movl  nb311_faction(%ebp),%edi
        movl  nb311_fshift(%ebp),%esi
        movl  nb311_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movapd nb311_fixO(%esp),%xmm0
        movapd nb311_fiyO(%esp),%xmm1
        movapd nb311_fizO(%esp),%xmm2

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
        movapd nb311_fixH1(%esp),%xmm0
        movapd nb311_fiyH1(%esp),%xmm1
        movapd nb311_fizH1(%esp),%xmm2

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
        movapd nb311_fixH2(%esp),%xmm0
        movapd nb311_fiyH2(%esp),%xmm1
        movapd nb311_fizH2(%esp),%xmm2

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
        movl nb311_n(%esp),%esi
        ## get group index for i particle 
        movl  nb311_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb311_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb311_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb311_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb311_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb311_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel311_ia32_sse2.nb311_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb311_n(%esp)
        jmp _nb_kernel311_ia32_sse2.nb311_outer
_nb_kernel311_ia32_sse2.nb311_outerend: 
        ## check if more outer neighborlists remain
        movl  nb311_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel311_ia32_sse2.nb311_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel311_ia32_sse2.nb311_threadloop
_nb_kernel311_ia32_sse2.nb311_end: 
        emms

        movl nb311_nouter(%esp),%eax
        movl nb311_ninner(%esp),%ebx
        movl nb311_outeriter(%ebp),%ecx
        movl nb311_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb311_salign(%esp),%eax
        addl %eax,%esp
        addl $844,%esp
        popl %edi
        popl %esi
    popl %edx
    popl %ecx
    popl %ebx
    popl %eax
        leave
        ret





.globl nb_kernel311nf_ia32_sse2
.globl _nb_kernel311nf_ia32_sse2
nb_kernel311nf_ia32_sse2:       
_nb_kernel311nf_ia32_sse2:      
.set nb311nf_p_nri, 8
.set nb311nf_iinr, 12
.set nb311nf_jindex, 16
.set nb311nf_jjnr, 20
.set nb311nf_shift, 24
.set nb311nf_shiftvec, 28
.set nb311nf_fshift, 32
.set nb311nf_gid, 36
.set nb311nf_pos, 40
.set nb311nf_faction, 44
.set nb311nf_charge, 48
.set nb311nf_p_facel, 52
.set nb311nf_argkrf, 56
.set nb311nf_argcrf, 60
.set nb311nf_Vc, 64
.set nb311nf_type, 68
.set nb311nf_p_ntype, 72
.set nb311nf_vdwparam, 76
.set nb311nf_Vvdw, 80
.set nb311nf_p_tabscale, 84
.set nb311nf_VFtab, 88
.set nb311nf_invsqrta, 92
.set nb311nf_dvda, 96
.set nb311nf_p_gbtabscale, 100
.set nb311nf_GBtab, 104
.set nb311nf_p_nthreads, 108
.set nb311nf_count, 112
.set nb311nf_mtx, 116
.set nb311nf_outeriter, 120
.set nb311nf_inneriter, 124
.set nb311nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb311nf_ixO, 0
.set nb311nf_iyO, 16
.set nb311nf_izO, 32
.set nb311nf_ixH1, 48
.set nb311nf_iyH1, 64
.set nb311nf_izH1, 80
.set nb311nf_ixH2, 96
.set nb311nf_iyH2, 112
.set nb311nf_izH2, 128
.set nb311nf_iqO, 144
.set nb311nf_iqH, 160
.set nb311nf_qqO, 176
.set nb311nf_qqH, 192
.set nb311nf_rinvO, 208
.set nb311nf_rinvH1, 224
.set nb311nf_rinvH2, 240
.set nb311nf_rO, 256
.set nb311nf_rH1, 272
.set nb311nf_rH2, 288
.set nb311nf_tsc, 304
.set nb311nf_c6, 320
.set nb311nf_c12, 336
.set nb311nf_vctot, 352
.set nb311nf_Vvdwtot, 368
.set nb311nf_half, 384
.set nb311nf_three, 400
.set nb311nf_is3, 416
.set nb311nf_ii3, 420
.set nb311nf_ntia, 424
.set nb311nf_innerjjnr, 428
.set nb311nf_innerk, 432
.set nb311nf_n, 436
.set nb311nf_nn1, 440
.set nb311nf_nri, 444
.set nb311nf_nouter, 448
.set nb311nf_ninner, 452
.set nb311nf_salign, 456
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $460,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb311nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb311nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb311nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb311nf_nouter(%esp)
        movl %eax,nb311nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb311nf_half(%esp)
        movl %ebx,nb311nf_half+4(%esp)
        movsd nb311nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb311nf_half(%esp)
        movapd %xmm3,nb311nf_three(%esp)

        movl nb311nf_p_tabscale(%ebp),%eax
        movsd (%eax),%xmm5
        shufpd $0,%xmm5,%xmm5
        movapd %xmm5,nb311nf_tsc(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb311nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb311nf_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        movsd 8(%edx,%ebx,8),%xmm4
        movl nb311nf_p_facel(%ebp),%esi
        movsd (%esi),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb311nf_iqO(%esp)
        movapd %xmm4,nb311nf_iqH(%esp)

        movl  nb311nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb311nf_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb311nf_ntia(%esp)
_nb_kernel311nf_ia32_sse2.nb311nf_threadloop: 
        movl  nb311nf_count(%ebp),%esi          ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel311nf_ia32_sse2.nb311nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel311nf_ia32_sse2.nb311nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb311nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb311nf_n(%esp)
        movl %ebx,nb311nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel311nf_ia32_sse2.nb311nf_outerstart
        jmp _nb_kernel311nf_ia32_sse2.nb311nf_end

_nb_kernel311nf_ia32_sse2.nb311nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb311nf_nouter(%esp),%ebx
        movl %ebx,nb311nf_nouter(%esp)

_nb_kernel311nf_ia32_sse2.nb311nf_outer: 
        movl  nb311nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 

        movl  nb311nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb311nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb311nf_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb311nf_ii3(%esp)

        addsd (%eax,%ebx,8),%xmm3
        addsd 8(%eax,%ebx,8),%xmm4
        addsd 16(%eax,%ebx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb311nf_ixO(%esp)
        movapd %xmm4,nb311nf_iyO(%esp)
        movapd %xmm5,nb311nf_izO(%esp)

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
        movapd %xmm0,nb311nf_ixH1(%esp)
        movapd %xmm1,nb311nf_iyH1(%esp)
        movapd %xmm2,nb311nf_izH1(%esp)
        movapd %xmm3,nb311nf_ixH2(%esp)
        movapd %xmm4,nb311nf_iyH2(%esp)
        movapd %xmm5,nb311nf_izH2(%esp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb311nf_vctot(%esp)
        movapd %xmm4,nb311nf_Vvdwtot(%esp)

        movl  nb311nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb311nf_pos(%ebp),%esi
        movl  nb311nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb311nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb311nf_ninner(%esp),%ecx
        movl  %ecx,nb311nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb311nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel311nf_ia32_sse2.nb311nf_unroll_loop
        jmp   _nb_kernel311nf_ia32_sse2.nb311nf_checksingle
_nb_kernel311nf_ia32_sse2.nb311nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb311nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb311nf_innerjjnr(%esp)             ## advance pointer (unrolled 2) 

        movl nb311nf_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb311nf_iqO(%esp),%xmm3
        mulpd  nb311nf_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movapd  %xmm3,nb311nf_qqO(%esp)
        movapd  %xmm4,nb311nf_qqH(%esp)

        movl nb311nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb311nf_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb311nf_ntia(%esp),%edi
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
        movapd %xmm4,nb311nf_c6(%esp)
        movapd %xmm6,nb311nf_c12(%esp)

        movl nb311nf_pos(%ebp),%esi        ## base of pos[] 

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
        movapd nb311nf_ixO(%esp),%xmm4
        movapd nb311nf_iyO(%esp),%xmm5
        movapd nb311nf_izO(%esp),%xmm6

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
        movapd nb311nf_ixH1(%esp),%xmm4
        movapd nb311nf_iyH1(%esp),%xmm5
        movapd nb311nf_izH1(%esp),%xmm6

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
        movapd nb311nf_ixH2(%esp),%xmm3
        movapd nb311nf_iyH2(%esp),%xmm4
        movapd nb311nf_izH2(%esp),%xmm5

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
        movapd  nb311nf_three(%esp),%xmm4
        mulpd   %xmm7,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb311nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm7,%xmm2
        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb311nf_three(%esp),%xmm4
        subpd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb311nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,nb311nf_rinvO(%esp)       ## rinvO in xmm4 
        mulpd   %xmm4,%xmm7
        movapd  %xmm7,nb311nf_rO(%esp)          ## r in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb311nf_three(%esp),%xmm4
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb311nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm6,%xmm2
        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb311nf_three(%esp),%xmm4
        subpd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb311nf_half(%esp),%xmm4   ## rinv 
        movapd %xmm4,nb311nf_rinvH1(%esp)       ## rinvH1 
        mulpd  %xmm4,%xmm6
        movapd %xmm6,nb311nf_rH1(%esp)          ## rH1 

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb311nf_three(%esp),%xmm4
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb311nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm5,%xmm2
        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb311nf_three(%esp),%xmm4
        subpd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb311nf_half(%esp),%xmm4   ## rinv 
        movapd %xmm4,nb311nf_rinvH2(%esp)   ## rinv 
        mulpd %xmm4,%xmm5
        movapd %xmm5,nb311nf_rH2(%esp)   ## r 

        ## do O interactions 
        ## rO is still in xmm7 
        movapd nb311nf_rinvO(%esp),%xmm0
        mulpd   nb311nf_tsc(%esp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movd %eax,%mm0
        movd %ebx,%mm1
        movl nb311nf_VFtab(%ebp),%esi
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
    movapd nb311nf_qqO(%esp),%xmm0
    mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addpd  %xmm4,%xmm5 ## xmm5=VV 
    mulpd  %xmm0,%xmm5 ## vcoul=qq*VV  

        ## do nontable L-J 
        movapd nb311nf_rinvO(%esp),%xmm2
        mulpd  %xmm2,%xmm2

    ## at this point mm5 contains vcoul and xmm0 fijC 
    ## increment vcoul - then we can get rid of mm5 
    addpd  nb311nf_vctot(%esp),%xmm5
    movapd %xmm5,nb311nf_vctot(%esp)

        movapd %xmm2,%xmm1
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm4
        mulpd  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulpd  nb311nf_c6(%esp),%xmm1
        mulpd  nb311nf_c12(%esp),%xmm4
        movapd %xmm4,%xmm3
        subpd  %xmm1,%xmm3      ## xmm3=Vvdw12-Vvdw6 
        addpd  nb311nf_Vvdwtot(%esp),%xmm3
        movapd %xmm3,nb311nf_Vvdwtot(%esp)


        ## Done with O interactions - now H1! 
        movapd nb311nf_rH1(%esp),%xmm7
        mulpd nb311nf_tsc(%esp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb311nf_VFtab(%ebp),%esi
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
        movapd nb311nf_qqH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV 
    ## at this point mm5 contains vcoul 
    ## increment vcoul 
    addpd  nb311nf_vctot(%esp),%xmm5
    movapd %xmm5,nb311nf_vctot(%esp)

        ## Done with H1, finally we do H2 interactions 
        movapd nb311nf_rH2(%esp),%xmm7
        mulpd   nb311nf_tsc(%esp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb311nf_VFtab(%ebp),%esi
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
        movapd nb311nf_qqH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV 
    ## at this point mm5 contains vcoul 
    ## increment vcoul 
    addpd  nb311nf_vctot(%esp),%xmm5
    movapd %xmm5,nb311nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $2,nb311nf_innerk(%esp)
        jl    _nb_kernel311nf_ia32_sse2.nb311nf_checksingle
        jmp   _nb_kernel311nf_ia32_sse2.nb311nf_unroll_loop
_nb_kernel311nf_ia32_sse2.nb311nf_checksingle: 
        movl  nb311nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel311nf_ia32_sse2.nb311nf_dosingle
        jmp   _nb_kernel311nf_ia32_sse2.nb311nf_updateouterdata
_nb_kernel311nf_ia32_sse2.nb311nf_dosingle: 
        movl  nb311nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movl nb311nf_charge(%ebp),%esi     ## base of charge[] 
        xorpd %xmm3,%xmm3
        movlpd (%esi,%eax,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb311nf_iqO(%esp),%xmm3
        mulpd  nb311nf_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movapd  %xmm3,nb311nf_qqO(%esp)
        movapd  %xmm4,nb311nf_qqH(%esp)

        movl nb311nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb311nf_vdwparam(%ebp),%esi
        shll %eax
        movl nb311nf_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movapd %xmm4,nb311nf_c6(%esp)
        movapd %xmm6,nb311nf_c12(%esp)

        movl nb311nf_pos(%ebp),%esi        ## base of pos[] 
        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move coords to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb311nf_ixO(%esp),%xmm4
        movapd nb311nf_iyO(%esp),%xmm5
        movapd nb311nf_izO(%esp),%xmm6

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
        movapd nb311nf_ixH1(%esp),%xmm4
        movapd nb311nf_iyH1(%esp),%xmm5
        movapd nb311nf_izH1(%esp),%xmm6

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
        movapd nb311nf_ixH2(%esp),%xmm3
        movapd nb311nf_iyH2(%esp),%xmm4
        movapd nb311nf_izH2(%esp),%xmm5

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
        movapd  nb311nf_three(%esp),%xmm4
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb311nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm7,%xmm2
        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb311nf_three(%esp),%xmm4
        subsd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb311nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,nb311nf_rinvO(%esp)       ## rinvO in xmm4 
        mulsd   %xmm4,%xmm7
        movapd  %xmm7,nb311nf_rO(%esp)          ## r in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb311nf_three(%esp),%xmm4
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb311nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm6,%xmm2
        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb311nf_three(%esp),%xmm4
        subsd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb311nf_half(%esp),%xmm4   ## rinv 
        movapd %xmm4,nb311nf_rinvH1(%esp)       ## rinvH1 
        mulsd  %xmm4,%xmm6
        movapd %xmm6,nb311nf_rH1(%esp)          ## rH1 

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb311nf_three(%esp),%xmm4
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb311nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm5,%xmm2
        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb311nf_three(%esp),%xmm4
        subsd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb311nf_half(%esp),%xmm4   ## rinv 
        movapd %xmm4,nb311nf_rinvH2(%esp)   ## rinv 
        mulsd %xmm4,%xmm5
        movapd %xmm5,nb311nf_rH2(%esp)   ## r 

        ## do O interactions 
        movd %eax,%mm0
        ## rO is still in xmm7 
        movapd nb311nf_rinvO(%esp),%xmm0
        mulsd   nb311nf_tsc(%esp),%xmm7
        cvttsd2si %xmm7,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb311nf_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1  

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7  
    mulsd  %xmm1,%xmm6      ## xmm6=Geps 
    mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
    addsd  %xmm6,%xmm5
    addsd  %xmm7,%xmm5      ## xmm5=Fp        
    movapd nb311nf_qqO(%esp),%xmm0
    mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addsd  %xmm4,%xmm5 ## xmm5=VV 
    mulsd  %xmm0,%xmm5 ## vcoul=qq*VV  

        ## do nontable L-J 
        movapd nb311nf_rinvO(%esp),%xmm2
        mulsd  %xmm2,%xmm2

    ## at this point mm5 contains vcoul 
    ## increment vcoul - then we can get rid of mm5 
    addsd  nb311nf_vctot(%esp),%xmm5
    movlpd %xmm5,nb311nf_vctot(%esp)

        movapd %xmm2,%xmm1
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm4
        mulsd  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulsd  nb311nf_c6(%esp),%xmm1
        mulsd  nb311nf_c12(%esp),%xmm4
        movapd %xmm4,%xmm3
        subsd  %xmm1,%xmm3      ## xmm3=Vvdw12-Vvdw6 
        addsd  nb311nf_Vvdwtot(%esp),%xmm3
        movlpd %xmm3,nb311nf_Vvdwtot(%esp)

        ## Done with O interactions - now H1! 
        movapd nb311nf_rH1(%esp),%xmm7
        mulpd nb311nf_tsc(%esp),%xmm7
        cvttsd2si %xmm7,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb311nf_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1  

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb311nf_qqH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 
    ## increment vcoul 
    addsd  nb311nf_vctot(%esp),%xmm5
    movlpd %xmm5,nb311nf_vctot(%esp)

        ## Done with H1, finally we do H2 interactions 
        movapd nb311nf_rH2(%esp),%xmm7
        mulsd   nb311nf_tsc(%esp),%xmm7
        cvttsd2si %xmm7,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb311nf_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1  

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb311nf_qqH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 
    ## increment vcoul 
    addsd  nb311nf_vctot(%esp),%xmm5
    movlpd %xmm5,nb311nf_vctot(%esp)

_nb_kernel311nf_ia32_sse2.nb311nf_updateouterdata: 
        ## get n from stack
        movl nb311nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb311nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb311nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb311nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb311nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb311nf_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb311nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel311nf_ia32_sse2.nb311nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb311nf_n(%esp)
        jmp _nb_kernel311nf_ia32_sse2.nb311nf_outer
_nb_kernel311nf_ia32_sse2.nb311nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb311nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel311nf_ia32_sse2.nb311nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel311nf_ia32_sse2.nb311nf_threadloop
_nb_kernel311nf_ia32_sse2.nb311nf_end: 
        emms

        movl nb311nf_nouter(%esp),%eax
        movl nb311nf_ninner(%esp),%ebx
        movl nb311nf_outeriter(%ebp),%ecx
        movl nb311nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb311nf_salign(%esp),%eax
        addl %eax,%esp
        addl $460,%esp
        popl %edi
        popl %esi
    popl %edx
    popl %ecx
    popl %ebx
    popl %eax
        leave
        ret


