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




.globl nb_kernel333_ia32_sse2
.globl _nb_kernel333_ia32_sse2
nb_kernel333_ia32_sse2: 
_nb_kernel333_ia32_sse2:        
.set nb333_p_nri, 8
.set nb333_iinr, 12
.set nb333_jindex, 16
.set nb333_jjnr, 20
.set nb333_shift, 24
.set nb333_shiftvec, 28
.set nb333_fshift, 32
.set nb333_gid, 36
.set nb333_pos, 40
.set nb333_faction, 44
.set nb333_charge, 48
.set nb333_p_facel, 52
.set nb333_argkrf, 56
.set nb333_argcrf, 60
.set nb333_Vc, 64
.set nb333_type, 68
.set nb333_p_ntype, 72
.set nb333_vdwparam, 76
.set nb333_Vvdw, 80
.set nb333_p_tabscale, 84
.set nb333_VFtab, 88
.set nb333_invsqrta, 92
.set nb333_dvda, 96
.set nb333_p_gbtabscale, 100
.set nb333_GBtab, 104
.set nb333_p_nthreads, 108
.set nb333_count, 112
.set nb333_mtx, 116
.set nb333_outeriter, 120
.set nb333_inneriter, 124
.set nb333_work, 128
        ## stack offsets for local variables 
        ## bottom of stack is cache-aligned for sse2 use 
.set nb333_ixO, 0
.set nb333_iyO, 16
.set nb333_izO, 32
.set nb333_ixH1, 48
.set nb333_iyH1, 64
.set nb333_izH1, 80
.set nb333_ixH2, 96
.set nb333_iyH2, 112
.set nb333_izH2, 128
.set nb333_ixM, 144
.set nb333_iyM, 160
.set nb333_izM, 176
.set nb333_iqM, 192
.set nb333_iqH, 208
.set nb333_dxO, 224
.set nb333_dyO, 240
.set nb333_dzO, 256
.set nb333_dxH1, 272
.set nb333_dyH1, 288
.set nb333_dzH1, 304
.set nb333_dxH2, 320
.set nb333_dyH2, 336
.set nb333_dzH2, 352
.set nb333_dxM, 368
.set nb333_dyM, 384
.set nb333_dzM, 400
.set nb333_qqM, 416
.set nb333_qqH, 432
.set nb333_rinvO, 448
.set nb333_rinvH1, 464
.set nb333_rinvH2, 480
.set nb333_rinvM, 496
.set nb333_rO, 512
.set nb333_rH1, 528
.set nb333_rH2, 544
.set nb333_rM, 560
.set nb333_tsc, 576
.set nb333_two, 592
.set nb333_c6, 608
.set nb333_c12, 624
.set nb333_vctot, 640
.set nb333_Vvdwtot, 656
.set nb333_fixO, 672
.set nb333_fiyO, 688
.set nb333_fizO, 704
.set nb333_fixH1, 720
.set nb333_fiyH1, 736
.set nb333_fizH1, 752
.set nb333_fixH2, 768
.set nb333_fiyH2, 784
.set nb333_fizH2, 800
.set nb333_fixM, 816
.set nb333_fiyM, 832
.set nb333_fizM, 848
.set nb333_fjx, 864
.set nb333_fjy, 880
.set nb333_fjz, 896
.set nb333_half, 912
.set nb333_three, 928
.set nb333_is3, 944
.set nb333_ii3, 948
.set nb333_ntia, 952
.set nb333_innerjjnr, 956
.set nb333_innerk, 960
.set nb333_n, 964
.set nb333_nn1, 968
.set nb333_nri, 972
.set nb333_nouter, 976
.set nb333_ninner, 980
.set nb333_salign, 984
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $988,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb333_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb333_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb333_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb333_nouter(%esp)
        movl %eax,nb333_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb333_half(%esp)
        movl %ebx,nb333_half+4(%esp)
        movsd nb333_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb333_half(%esp)
        movapd %xmm2,nb333_two(%esp)
        movapd %xmm3,nb333_three(%esp)
        movl nb333_p_tabscale(%ebp),%eax
        movsd (%eax),%xmm5
        shufpd $0,%xmm5,%xmm5
        movapd %xmm5,nb333_tsc(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb333_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb333_charge(%ebp),%edx
        movsd 8(%edx,%ebx,8),%xmm3
        movsd 24(%edx,%ebx,8),%xmm4
        movl nb333_p_facel(%ebp),%esi
        movsd (%esi),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb333_iqH(%esp)
        movapd %xmm4,nb333_iqM(%esp)

        movl  nb333_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb333_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb333_ntia(%esp)
_nb_kernel333_ia32_sse2.nb333_threadloop: 
        movl  nb333_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel333_ia32_sse2.nb333_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel333_ia32_sse2.nb333_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb333_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb333_n(%esp)
        movl %ebx,nb333_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel333_ia32_sse2.nb333_outerstart
        jmp _nb_kernel333_ia32_sse2.nb333_end

_nb_kernel333_ia32_sse2.nb333_outerstart: 
        ## ebx contains number of outer iterations
        addl nb333_nouter(%esp),%ebx
        movl %ebx,nb333_nouter(%esp)

_nb_kernel333_ia32_sse2.nb333_outer: 
        movl  nb333_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb333_is3(%esp)      ## store is3 

        movl  nb333_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb333_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        movapd %xmm0,%xmm6
        movapd %xmm1,%xmm7

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb333_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb333_ii3(%esp)

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
        movapd %xmm3,nb333_ixO(%esp)
        movapd %xmm4,nb333_iyO(%esp)
        movapd %xmm5,nb333_izO(%esp)
        movapd %xmm6,nb333_ixH1(%esp)
        movapd %xmm7,nb333_iyH1(%esp)

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
        movapd %xmm6,nb333_izH1(%esp)
        movapd %xmm0,nb333_ixH2(%esp)
        movapd %xmm1,nb333_iyH2(%esp)
        movapd %xmm2,nb333_izH2(%esp)
        movapd %xmm3,nb333_ixM(%esp)
        movapd %xmm4,nb333_iyM(%esp)
        movapd %xmm5,nb333_izM(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb333_vctot(%esp)
        movapd %xmm4,nb333_Vvdwtot(%esp)
        movapd %xmm4,nb333_fixO(%esp)
        movapd %xmm4,nb333_fiyO(%esp)
        movapd %xmm4,nb333_fizO(%esp)
        movapd %xmm4,nb333_fixH1(%esp)
        movapd %xmm4,nb333_fiyH1(%esp)
        movapd %xmm4,nb333_fizH1(%esp)
        movapd %xmm4,nb333_fixH2(%esp)
        movapd %xmm4,nb333_fiyH2(%esp)
        movapd %xmm4,nb333_fizH2(%esp)
        movapd %xmm4,nb333_fixM(%esp)
        movapd %xmm4,nb333_fiyM(%esp)
        movapd %xmm4,nb333_fizM(%esp)

        movl  nb333_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx     ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb333_pos(%ebp),%esi
        movl  nb333_faction(%ebp),%edi
        movl  nb333_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb333_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb333_ninner(%esp),%ecx
        movl  %ecx,nb333_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb333_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel333_ia32_sse2.nb333_unroll_loop
        jmp   _nb_kernel333_ia32_sse2.nb333_checksingle
_nb_kernel333_ia32_sse2.nb333_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb333_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb333_innerjjnr(%esp)             ## advance pointer (unrolled 2) 

        movl nb333_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb333_iqM(%esp),%xmm3
        mulpd  nb333_iqH(%esp),%xmm4
        movapd  %xmm3,nb333_qqM(%esp)
        movapd  %xmm4,nb333_qqH(%esp)

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movl nb333_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb333_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb333_ntia(%esp),%edi
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
        movapd %xmm4,nb333_c6(%esp)
        movapd %xmm6,nb333_c12(%esp)

        movl nb333_pos(%ebp),%esi        ## base of pos[] 

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
        movapd nb333_ixO(%esp),%xmm4
        movapd nb333_iyO(%esp),%xmm5
        movapd nb333_izO(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb333_dxO(%esp)
        movapd %xmm5,nb333_dyO(%esp)
        movapd %xmm6,nb333_dzO(%esp)
        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb333_ixH1(%esp),%xmm4
        movapd nb333_iyH1(%esp),%xmm5
        movapd nb333_izH1(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb333_dxH1(%esp)
        movapd %xmm5,nb333_dyH1(%esp)
        movapd %xmm6,nb333_dzH1(%esp)
        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm6
        addpd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb333_ixH2(%esp),%xmm3
        movapd nb333_iyH2(%esp),%xmm4
        movapd nb333_izH2(%esp),%xmm5

        ## calc dr 
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5

        ## store dr 
        movapd %xmm3,nb333_dxH2(%esp)
        movapd %xmm4,nb333_dyH2(%esp)
        movapd %xmm5,nb333_dzH2(%esp)
        ## square it 
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        addpd %xmm4,%xmm5
        addpd %xmm3,%xmm5
        ## move ixM-izM to xmm2-xmm4  
        movapd nb333_iyM(%esp),%xmm3
        movapd nb333_izM(%esp),%xmm4
        subpd  %xmm1,%xmm3
        subpd  %xmm2,%xmm4
        movapd nb333_ixM(%esp),%xmm2
        subpd  %xmm0,%xmm2

        ## store dr 
        movapd %xmm2,nb333_dxM(%esp)
        movapd %xmm3,nb333_dyM(%esp)
        movapd %xmm4,nb333_dzM(%esp)
        ## square it 
        mulpd %xmm2,%xmm2
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        addpd %xmm3,%xmm4
        addpd %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## rsqO - put seed in xmm2 
        cvtpd2ps %xmm7,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb333_three(%esp),%xmm0
        mulpd   %xmm7,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulpd   nb333_half(%esp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm7,%xmm2
        movapd %xmm0,%xmm3
        mulpd %xmm0,%xmm0       ## lu*lu 
        mulpd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb333_three(%esp),%xmm0
        subpd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulpd nb333_half(%esp),%xmm0   ## rinv 
        movapd  %xmm0,nb333_rinvO(%esp)         ## rinvO in xmm0
        mulpd   %xmm0,%xmm7
        movapd  %xmm7,nb333_rO(%esp)    ## r in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb333_three(%esp),%xmm0
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulpd   nb333_half(%esp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm6,%xmm2
        movapd %xmm0,%xmm3
        mulpd %xmm0,%xmm0       ## lu*lu 
        mulpd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb333_three(%esp),%xmm0
        subpd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulpd nb333_half(%esp),%xmm0   ## rinv 
        movapd %xmm0,nb333_rinvH1(%esp)         ## rinvH1 
        mulpd  %xmm0,%xmm6
        movapd %xmm6,nb333_rH1(%esp)    ## rH1 

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb333_three(%esp),%xmm0
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulpd   nb333_half(%esp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm5,%xmm2
        movapd %xmm0,%xmm3
        mulpd %xmm0,%xmm0       ## lu*lu 
        mulpd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb333_three(%esp),%xmm0
        subpd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulpd nb333_half(%esp),%xmm0   ## rinv 
        movapd %xmm0,nb333_rinvH2(%esp)   ## rinv 
        mulpd %xmm0,%xmm5
        movapd %xmm5,nb333_rH2(%esp)   ## r 

        ## rsqM - seed in xmm2 
        cvtpd2ps %xmm4,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb333_three(%esp),%xmm0
        mulpd   %xmm4,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulpd   nb333_half(%esp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm4,%xmm2
        movapd %xmm0,%xmm3
        mulpd %xmm0,%xmm0       ## lu*lu 
        mulpd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb333_three(%esp),%xmm0
        subpd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulpd nb333_half(%esp),%xmm0   ## rinv 
        movapd %xmm0,nb333_rinvM(%esp)   ## rinv 
        mulpd %xmm0,%xmm4
        movapd %xmm4,nb333_rM(%esp)   ## r 

        ## do O interactions 
        ## rO is still in xmm7 
        movapd nb333_rinvO(%esp),%xmm0
        mulpd   nb333_tsc(%esp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movd %eax,%mm0
        movd %ebx,%mm1
        movl nb333_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        leal  (%eax,%eax,2),%eax ## idx *= 3 (total *=12 now) 
        leal  (%ebx,%ebx,2),%ebx

        ## Dispersion 
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
        ## Dispersion table ready, in xmm4-xmm7                 
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb333_two(%esp),%xmm7    ## two*Heps2 
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb333_c6(%esp),%xmm4
        mulpd  %xmm4,%xmm7       ## fijD 
        mulpd  %xmm4,%xmm5       ## Vvdw6 
        movapd  %xmm7,%xmm0     ## fscal summation register 

        addpd  nb333_Vvdwtot(%esp),%xmm5
        movapd %xmm5,nb333_Vvdwtot(%esp)

        ## Repulsion 
        movlpd 64(%esi,%eax,8),%xmm4    ## Y1
        movlpd 64(%esi,%ebx,8),%xmm3    ## Y2
        movhpd 72(%esi,%eax,8),%xmm4    ## Y1 F1        
        movhpd 72(%esi,%ebx,8),%xmm3    ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 80(%esi,%eax,8),%xmm6    ## G1
        movlpd 80(%esi,%ebx,8),%xmm3    ## G2
        movhpd 88(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 88(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## Dispersion table ready, in xmm4-xmm7                 
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb333_two(%esp),%xmm7    ## two*Heps2 
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb333_c12(%esp),%xmm4
        mulpd  %xmm4,%xmm7 ## fijR 
        mulpd  %xmm4,%xmm5 ## Vvdw12 
        addpd  %xmm0,%xmm7      ## fscal sum
        mulpd nb333_tsc(%esp),%xmm7
        mulpd nb333_rinvO(%esp),%xmm7

        addpd  nb333_Vvdwtot(%esp),%xmm5
        movapd %xmm5,nb333_Vvdwtot(%esp)

        xorpd %xmm4,%xmm4
        subpd %xmm7,%xmm4
        movapd nb333_dxO(%esp),%xmm0
        movapd nb333_dyO(%esp),%xmm1
        movapd nb333_dzO(%esp),%xmm2

        movd %mm0,%eax
        movd %mm1,%ebx

        movapd nb333_dxO(%esp),%xmm0
        movapd nb333_dyO(%esp),%xmm1
        movapd nb333_dzO(%esp),%xmm2
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2      ## tx in xmm0-xmm2 

        ## update O forces 
        movapd nb333_fixO(%esp),%xmm3
        movapd nb333_fiyO(%esp),%xmm4
        movapd nb333_fizO(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb333_fixO(%esp)
        movapd %xmm4,nb333_fiyO(%esp)
        movapd %xmm7,nb333_fizO(%esp)
        ## update j forces with water O 
        movapd %xmm0,nb333_fjx(%esp)
        movapd %xmm1,nb333_fjy(%esp)
        movapd %xmm2,nb333_fjz(%esp)

        ## Done with O interactions - now H1! 
        movapd nb333_rH1(%esp),%xmm7
        mulpd nb333_tsc(%esp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb333_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        leal  (%eax,%eax,2),%eax ## idx *= 3 (total *=12 now)   
        leal  (%ebx,%ebx,2),%ebx

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
        mulpd  nb333_two(%esp),%xmm7    ## two*Heps2 
        movapd nb333_qqH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 
        ## increment vcoul 
        xorpd  %xmm4,%xmm4
        addpd  nb333_vctot(%esp),%xmm5
        mulpd  nb333_rinvH1(%esp),%xmm3
        movapd %xmm5,nb333_vctot(%esp)
        mulpd  nb333_tsc(%esp),%xmm3
        subpd %xmm3,%xmm4

        movapd nb333_dxH1(%esp),%xmm0
        movapd nb333_dyH1(%esp),%xmm1
        movapd nb333_dzH1(%esp),%xmm2
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb333_fixH1(%esp),%xmm3
        movapd nb333_fiyH1(%esp),%xmm4
        movapd nb333_fizH1(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb333_fixH1(%esp)
        movapd %xmm4,nb333_fiyH1(%esp)
        movapd %xmm7,nb333_fizH1(%esp)
        ## update j forces with water H1 
        addpd  nb333_fjx(%esp),%xmm0
        addpd  nb333_fjy(%esp),%xmm1
        addpd  nb333_fjz(%esp),%xmm2
        movapd %xmm0,nb333_fjx(%esp)
        movapd %xmm1,nb333_fjy(%esp)
        movapd %xmm2,nb333_fjz(%esp)

        ## H2 interactions 
        movapd nb333_rH2(%esp),%xmm7
        mulpd   nb333_tsc(%esp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb333_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        leal  (%eax,%eax,2),%eax ## idx *= 3 (total *=12 now)
        leal  (%ebx,%ebx,2),%ebx

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
        mulpd  nb333_two(%esp),%xmm7    ## two*Heps2 
        movapd nb333_qqH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 
        ## increment vcoul 
        xorpd  %xmm4,%xmm4
        addpd  nb333_vctot(%esp),%xmm5
        mulpd  nb333_rinvH2(%esp),%xmm3
        movapd %xmm5,nb333_vctot(%esp)
        mulpd  nb333_tsc(%esp),%xmm3
        subpd  %xmm3,%xmm4

        movapd nb333_dxH2(%esp),%xmm0
        movapd nb333_dyH2(%esp),%xmm1
        movapd nb333_dzH2(%esp),%xmm2
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

    movd %mm0,%eax
    movd %mm1,%ebx

        ## update H2 forces 
        movapd nb333_fixH2(%esp),%xmm3
        movapd nb333_fiyH2(%esp),%xmm4
        movapd nb333_fizH2(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb333_fixH2(%esp)
        movapd %xmm4,nb333_fiyH2(%esp)
        movapd %xmm7,nb333_fizH2(%esp)
        ## update j forces with water H1 
        addpd  nb333_fjx(%esp),%xmm0
        addpd  nb333_fjy(%esp),%xmm1
        addpd  nb333_fjz(%esp),%xmm2
        movapd %xmm0,nb333_fjx(%esp)
        movapd %xmm1,nb333_fjy(%esp)
        movapd %xmm2,nb333_fjz(%esp)

        ## M interactions 
        movapd nb333_rM(%esp),%xmm7
        mulpd   nb333_tsc(%esp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb333_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        leal  (%eax,%eax,2),%eax ## idx *= 3 (total *=12 now)
        leal  (%ebx,%ebx,2),%ebx

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
        mulpd  nb333_two(%esp),%xmm7    ## two*Heps2 
        movapd nb333_qqM(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 
        ## increment vcoul 
        xorpd  %xmm4,%xmm4
        addpd  nb333_vctot(%esp),%xmm5
        mulpd  nb333_rinvM(%esp),%xmm3
        movapd %xmm5,nb333_vctot(%esp)
        mulpd  nb333_tsc(%esp),%xmm3
        subpd  %xmm3,%xmm4

        movapd nb333_dxM(%esp),%xmm0
        movapd nb333_dyM(%esp),%xmm1
        movapd nb333_dzM(%esp),%xmm2
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        movd %mm0,%eax
        movd %mm1,%ebx

        ## update H2 forces 
        movapd nb333_fixM(%esp),%xmm3
        movapd nb333_fiyM(%esp),%xmm4
        movapd nb333_fizM(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb333_fixM(%esp)
        movapd %xmm4,nb333_fiyM(%esp)
        movapd %xmm7,nb333_fizM(%esp)

        movl nb333_faction(%ebp),%edi
        ## update j forces 
        addpd  nb333_fjx(%esp),%xmm0
        addpd  nb333_fjy(%esp),%xmm1
        addpd  nb333_fjz(%esp),%xmm2

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
        subl $2,nb333_innerk(%esp)
        jl    _nb_kernel333_ia32_sse2.nb333_checksingle
        jmp   _nb_kernel333_ia32_sse2.nb333_unroll_loop
_nb_kernel333_ia32_sse2.nb333_checksingle: 
        movl  nb333_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel333_ia32_sse2.nb333_dosingle
        jmp   _nb_kernel333_ia32_sse2.nb333_updateouterdata
_nb_kernel333_ia32_sse2.nb333_dosingle: 
        movl  nb333_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movl nb333_charge(%ebp),%esi     ## base of charge[] 
        xorpd %xmm3,%xmm3
        movlpd (%esi,%eax,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb333_iqM(%esp),%xmm3
        mulpd  nb333_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movapd  %xmm3,nb333_qqM(%esp)
        movapd  %xmm4,nb333_qqH(%esp)

        movl nb333_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb333_vdwparam(%ebp),%esi
        shll %eax
        movl nb333_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movapd %xmm4,nb333_c6(%esp)
        movapd %xmm6,nb333_c12(%esp)

        movl nb333_pos(%ebp),%esi        ## base of pos[] 
        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move coords to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb333_ixO(%esp),%xmm4
        movapd nb333_iyO(%esp),%xmm5
        movapd nb333_izO(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb333_dxO(%esp)
        movapd %xmm5,nb333_dyO(%esp)
        movapd %xmm6,nb333_dzO(%esp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb333_ixH1(%esp),%xmm4
        movapd nb333_iyH1(%esp),%xmm5
        movapd nb333_izH1(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb333_dxH1(%esp)
        movapd %xmm5,nb333_dyH1(%esp)
        movapd %xmm6,nb333_dzH1(%esp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm6
        addsd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb333_ixH2(%esp),%xmm3
        movapd nb333_iyH2(%esp),%xmm4
        movapd nb333_izH2(%esp),%xmm5

        ## calc dr 
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5

        ## store dr 
        movapd %xmm3,nb333_dxH2(%esp)
        movapd %xmm4,nb333_dyH2(%esp)
        movapd %xmm5,nb333_dzH2(%esp)
        ## square it 
        mulsd %xmm3,%xmm3
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        addsd %xmm4,%xmm5
        addsd %xmm3,%xmm5
        ## move ixM-izM to xmm2-xmm4  
        movapd nb333_iyM(%esp),%xmm3
        movapd nb333_izM(%esp),%xmm4
        subpd  %xmm1,%xmm3
        subpd  %xmm2,%xmm4
        movapd nb333_ixM(%esp),%xmm2
        subpd  %xmm0,%xmm2

        ## store dr 
        movapd %xmm2,nb333_dxM(%esp)
        movapd %xmm3,nb333_dyM(%esp)
        movapd %xmm4,nb333_dzM(%esp)
        ## square it 
        mulpd %xmm2,%xmm2
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        addpd %xmm3,%xmm4
        addpd %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## start with rsqO - put seed in xmm2 
        cvtsd2ss %xmm7,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb333_three(%esp),%xmm0
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulsd   nb333_half(%esp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm7,%xmm2
        movapd %xmm0,%xmm3
        mulsd %xmm0,%xmm0       ## lu*lu 
        mulsd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb333_three(%esp),%xmm0
        subsd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulsd nb333_half(%esp),%xmm0   ## rinv 
        movapd  %xmm0,nb333_rinvO(%esp)         ## rinvO in xmm0 
        mulsd   %xmm0,%xmm7
        movapd  %xmm7,nb333_rO(%esp)    ## r in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb333_three(%esp),%xmm0
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulsd   nb333_half(%esp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm6,%xmm2
        movapd %xmm0,%xmm3
        mulsd %xmm0,%xmm0       ## lu*lu 
        mulsd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb333_three(%esp),%xmm0
        subsd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulsd nb333_half(%esp),%xmm0   ## rinv 
        movapd %xmm0,nb333_rinvH1(%esp)         ## rinvH1 
        mulsd  %xmm0,%xmm6
        movapd %xmm6,nb333_rH1(%esp)    ## rH1 

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb333_three(%esp),%xmm0
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulsd   nb333_half(%esp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm5,%xmm2
        movapd %xmm0,%xmm3
        mulsd %xmm0,%xmm0       ## lu*lu 
        mulsd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb333_three(%esp),%xmm0
        subsd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulsd nb333_half(%esp),%xmm0   ## rinv 
        movapd %xmm0,nb333_rinvH2(%esp)   ## rinv 
        mulsd %xmm0,%xmm5
        movapd %xmm5,nb333_rH2(%esp)   ## r 

        ## rsqM - seed in xmm2 
        cvtsd2ss %xmm4,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb333_three(%esp),%xmm0
        mulsd   %xmm4,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulsd   nb333_half(%esp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm4,%xmm2
        movapd %xmm0,%xmm3
        mulsd %xmm0,%xmm0       ## lu*lu 
        mulsd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb333_three(%esp),%xmm0
        subsd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulsd nb333_half(%esp),%xmm0   ## rinv 
        movapd %xmm0,nb333_rinvM(%esp)   ## rinv 
        mulsd %xmm0,%xmm4
        movapd %xmm4,nb333_rM(%esp)   ## r 

        ## do O interactions 
        movd %eax,%mm0
        ## rO is still in xmm7 
        mulsd   nb333_tsc(%esp),%xmm7
        cvttsd2si %xmm7,%eax    ## lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm7,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax
        movl nb333_VFtab(%ebp),%esi
        leal  (%eax,%eax,2),%eax ## idx *= 3 (total *=12 now)   

        ## Dispersion 
        movsd 32(%esi,%eax,8),%xmm4     ## Y1   
        movsd 40(%esi,%eax,8),%xmm5     ## F1   
        movsd 48(%esi,%eax,8),%xmm6     ## G1   
        movsd 56(%esi,%eax,8),%xmm7     ## H1   
        ## Dispersion table ready, in xmm4-xmm7                 
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb333_two(%esp),%xmm7    ## two*Heps2 
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb333_c6(%esp),%xmm4
        mulsd  %xmm4,%xmm7       ## fijD 
        mulsd  %xmm4,%xmm5       ## Vvdw6 
        movsd  %xmm7,%xmm0      ## fscal 

        addsd  nb333_Vvdwtot(%esp),%xmm5
        movsd %xmm5,nb333_Vvdwtot(%esp)

        ## Repulsion 
        movsd 64(%esi,%eax,8),%xmm4     ## Y1   
        movsd 72(%esi,%eax,8),%xmm5     ## F1   
        movsd 80(%esi,%eax,8),%xmm6     ## G1   
        movsd 88(%esi,%eax,8),%xmm7     ## H1   

        unpcklpd %xmm3,%xmm6
        unpckhpd %xmm3,%xmm7
        ## Dispersion table ready, in xmm4-xmm7                 
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb333_two(%esp),%xmm7    ## two*Heps2 
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb333_c12(%esp),%xmm4
        mulsd  %xmm4,%xmm7 ## fijR 
        mulsd  %xmm4,%xmm5 ## Vvdw12 
        addsd  %xmm0,%xmm7      ## fscal summation 
        mulsd nb333_tsc(%esp),%xmm7
        mulsd nb333_rinvO(%esp),%xmm7

        addsd  nb333_Vvdwtot(%esp),%xmm5
        movsd %xmm5,nb333_Vvdwtot(%esp)

        xorpd %xmm4,%xmm4
        subsd %xmm7,%xmm4

        movapd nb333_dxO(%esp),%xmm0
        movapd nb333_dyO(%esp),%xmm1
        movapd nb333_dzO(%esp),%xmm2
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2      ## tx in xmm0-xmm2 

        ## update O forces 
        movapd nb333_fixO(%esp),%xmm3
        movapd nb333_fiyO(%esp),%xmm4
        movapd nb333_fizO(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb333_fixO(%esp)
        movlpd %xmm4,nb333_fiyO(%esp)
        movlpd %xmm7,nb333_fizO(%esp)
        ## update j forces with water O 
        movlpd %xmm0,nb333_fjx(%esp)
        movlpd %xmm1,nb333_fjy(%esp)
        movlpd %xmm2,nb333_fjz(%esp)

        ## Done with O interactions - now H1! 
        movapd nb333_rH1(%esp),%xmm7
        mulpd nb333_tsc(%esp),%xmm7
        cvttsd2si %xmm7,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb333_VFtab(%ebp),%esi
        leal  (%eax,%eax,2),%eax ## idx *= 3 (total *=12 now)   

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb333_two(%esp),%xmm7    ## two*Heps2 
        movapd nb333_qqH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 
    ## increment vcoul 
        xorpd  %xmm4,%xmm4
    addsd  nb333_vctot(%esp),%xmm5
        mulsd  nb333_rinvH1(%esp),%xmm3
    movlpd %xmm5,nb333_vctot(%esp)
        mulsd  nb333_tsc(%esp),%xmm3
        subsd %xmm3,%xmm4

        movapd nb333_dxH1(%esp),%xmm0
        movapd nb333_dyH1(%esp),%xmm1
        movapd nb333_dzH1(%esp),%xmm2
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb333_fixH1(%esp),%xmm3
        movapd nb333_fiyH1(%esp),%xmm4
        movapd nb333_fizH1(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb333_fixH1(%esp)
        movlpd %xmm4,nb333_fiyH1(%esp)
        movlpd %xmm7,nb333_fizH1(%esp)
        ## update j forces with water H1 
        addsd  nb333_fjx(%esp),%xmm0
        addsd  nb333_fjy(%esp),%xmm1
        addsd  nb333_fjz(%esp),%xmm2
        movlpd %xmm0,nb333_fjx(%esp)
        movlpd %xmm1,nb333_fjy(%esp)
        movlpd %xmm2,nb333_fjz(%esp)

        ##  H2 interactions 
        movapd nb333_rH2(%esp),%xmm7
        mulsd   nb333_tsc(%esp),%xmm7
        cvttsd2si %xmm7,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb333_VFtab(%ebp),%esi
        leal  (%eax,%eax,2),%eax ## idx *= 3 (total *=12 now)   

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   

        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb333_two(%esp),%xmm7    ## two*Heps2 
        movapd nb333_qqH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 
        ## increment vcoul 
        xorpd  %xmm4,%xmm4
        addsd  nb333_vctot(%esp),%xmm5
        mulsd  nb333_rinvH2(%esp),%xmm3
        movlpd %xmm5,nb333_vctot(%esp)
        mulsd  nb333_tsc(%esp),%xmm3
        subsd  %xmm3,%xmm4

        movapd nb333_dxH2(%esp),%xmm0
        movapd nb333_dyH2(%esp),%xmm1
        movapd nb333_dzH2(%esp),%xmm2
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        movd %mm0,%eax

        ## update H2 forces 
        movapd nb333_fixH2(%esp),%xmm3
        movapd nb333_fiyH2(%esp),%xmm4
        movapd nb333_fizH2(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb333_fixH2(%esp)
        movlpd %xmm4,nb333_fiyH2(%esp)
        movlpd %xmm7,nb333_fizH2(%esp)
        ## update j forces with water H1 
        addsd  nb333_fjx(%esp),%xmm0
        addsd  nb333_fjy(%esp),%xmm1
        addsd  nb333_fjz(%esp),%xmm2
        movlpd %xmm0,nb333_fjx(%esp)
        movlpd %xmm1,nb333_fjy(%esp)
        movlpd %xmm2,nb333_fjz(%esp)

        ## M interactions 
        movapd nb333_rM(%esp),%xmm7
        mulsd   nb333_tsc(%esp),%xmm7
        cvttsd2si %xmm7,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb333_VFtab(%ebp),%esi
        leal  (%eax,%eax,2),%eax ## idx *= 3 (total *=12 now)   

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   

        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb333_two(%esp),%xmm7    ## two*Heps2 
        movapd nb333_qqM(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 
        ## increment vcoul 
        xorpd  %xmm4,%xmm4
        addsd  nb333_vctot(%esp),%xmm5
        mulsd  nb333_rinvM(%esp),%xmm3
        movlpd %xmm5,nb333_vctot(%esp)
        mulsd  nb333_tsc(%esp),%xmm3
        subsd  %xmm3,%xmm4

        movapd nb333_dxM(%esp),%xmm0
        movapd nb333_dyM(%esp),%xmm1
        movapd nb333_dzM(%esp),%xmm2
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        movd %mm0,%eax

        ## update M forces 
        movapd nb333_fixM(%esp),%xmm3
        movapd nb333_fiyM(%esp),%xmm4
        movapd nb333_fizM(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb333_fixM(%esp)
        movlpd %xmm4,nb333_fiyM(%esp)
        movlpd %xmm7,nb333_fizM(%esp)

        movl nb333_faction(%ebp),%edi
        ## update j forces 
        addsd  nb333_fjx(%esp),%xmm0
        addsd  nb333_fjy(%esp),%xmm1
        addsd  nb333_fjz(%esp),%xmm2

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

_nb_kernel333_ia32_sse2.nb333_updateouterdata: 
        movl  nb333_ii3(%esp),%ecx
        movl  nb333_faction(%ebp),%edi
        movl  nb333_fshift(%ebp),%esi
        movl  nb333_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movapd nb333_fixO(%esp),%xmm0
        movapd nb333_fiyO(%esp),%xmm1
        movapd nb333_fizO(%esp),%xmm2

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
        movapd nb333_fixH1(%esp),%xmm0
        movapd nb333_fiyH1(%esp),%xmm1
        movapd nb333_fizH1(%esp),%xmm2

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
        movapd nb333_fixH2(%esp),%xmm0
        movapd nb333_fiyH2(%esp),%xmm1
        movapd nb333_fizH2(%esp),%xmm2

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
        movapd nb333_fixM(%esp),%xmm0
        movapd nb333_fiyM(%esp),%xmm1
        movapd nb333_fizM(%esp),%xmm2

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
        movl nb333_n(%esp),%esi
        ## get group index for i particle 
        movl  nb333_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb333_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb333_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb333_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb333_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb333_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel333_ia32_sse2.nb333_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb333_n(%esp)
        jmp _nb_kernel333_ia32_sse2.nb333_outer
_nb_kernel333_ia32_sse2.nb333_outerend: 
        ## check if more outer neighborlists remain
        movl  nb333_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel333_ia32_sse2.nb333_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel333_ia32_sse2.nb333_threadloop
_nb_kernel333_ia32_sse2.nb333_end: 
        emms

        movl nb333_nouter(%esp),%eax
        movl nb333_ninner(%esp),%ebx
        movl nb333_outeriter(%ebp),%ecx
        movl nb333_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb333_salign(%esp),%eax
        addl %eax,%esp
        addl $988,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret






.globl nb_kernel333nf_ia32_sse2
.globl _nb_kernel333nf_ia32_sse2
nb_kernel333nf_ia32_sse2:       
_nb_kernel333nf_ia32_sse2:      
.set nb333nf_p_nri, 8
.set nb333nf_iinr, 12
.set nb333nf_jindex, 16
.set nb333nf_jjnr, 20
.set nb333nf_shift, 24
.set nb333nf_shiftvec, 28
.set nb333nf_fshift, 32
.set nb333nf_gid, 36
.set nb333nf_pos, 40
.set nb333nf_faction, 44
.set nb333nf_charge, 48
.set nb333nf_p_facel, 52
.set nb333nf_argkrf, 56
.set nb333nf_argcrf, 60
.set nb333nf_Vc, 64
.set nb333nf_type, 68
.set nb333nf_p_ntype, 72
.set nb333nf_vdwparam, 76
.set nb333nf_Vvdw, 80
.set nb333nf_p_tabscale, 84
.set nb333nf_VFtab, 88
.set nb333nf_invsqrta, 92
.set nb333nf_dvda, 96
.set nb333nf_p_gbtabscale, 100
.set nb333nf_GBtab, 104
.set nb333nf_p_nthreads, 108
.set nb333nf_count, 112
.set nb333nf_mtx, 116
.set nb333nf_outeriter, 120
.set nb333nf_inneriter, 124
.set nb333nf_work, 128
        ## stack offsets for local variables 
        ## bottom of stack is cache-aligned for sse2 use 
.set nb333nf_ixO, 0
.set nb333nf_iyO, 16
.set nb333nf_izO, 32
.set nb333nf_ixH1, 48
.set nb333nf_iyH1, 64
.set nb333nf_izH1, 80
.set nb333nf_ixH2, 96
.set nb333nf_iyH2, 112
.set nb333nf_izH2, 128
.set nb333nf_ixM, 144
.set nb333nf_iyM, 160
.set nb333nf_izM, 176
.set nb333nf_iqM, 192
.set nb333nf_iqH, 208
.set nb333nf_qqM, 224
.set nb333nf_qqH, 240
.set nb333nf_rO, 256
.set nb333nf_rH1, 272
.set nb333nf_rH2, 288
.set nb333nf_rM, 304
.set nb333nf_tsc, 320
.set nb333nf_c6, 336
.set nb333nf_c12, 352
.set nb333nf_vctot, 368
.set nb333nf_Vvdwtot, 384
.set nb333nf_half, 400
.set nb333nf_three, 416
.set nb333nf_is3, 432
.set nb333nf_ii3, 436
.set nb333nf_ntia, 440
.set nb333nf_innerjjnr, 444
.set nb333nf_innerk, 448
.set nb333nf_n, 452
.set nb333nf_nn1, 456
.set nb333nf_nri, 460
.set nb333nf_nouter, 464
.set nb333nf_ninner, 468
.set nb333nf_salign, 472
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $476,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb333nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb333nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb333nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb333nf_nouter(%esp)
        movl %eax,nb333nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb333nf_half(%esp)
        movl %ebx,nb333nf_half+4(%esp)
        movsd nb333nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb333nf_half(%esp)
        movapd %xmm3,nb333nf_three(%esp)
        movl nb333nf_p_tabscale(%ebp),%eax
        movsd (%eax),%xmm5
        shufpd $0,%xmm5,%xmm5
        movapd %xmm5,nb333nf_tsc(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb333nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb333nf_charge(%ebp),%edx
        movsd 8(%edx,%ebx,8),%xmm3
        movsd 24(%edx,%ebx,8),%xmm4
        movl nb333nf_p_facel(%ebp),%esi
        movsd (%esi),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb333nf_iqH(%esp)
        movapd %xmm4,nb333nf_iqM(%esp)

        movl  nb333nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb333nf_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb333nf_ntia(%esp)
_nb_kernel333nf_ia32_sse2.nb333nf_threadloop: 
        movl  nb333nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel333nf_ia32_sse2.nb333nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel333nf_ia32_sse2.nb333nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb333nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb333nf_n(%esp)
        movl %ebx,nb333nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel333nf_ia32_sse2.nb333nf_outerstart
        jmp _nb_kernel333nf_ia32_sse2.nb333nf_end

_nb_kernel333nf_ia32_sse2.nb333nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb333nf_nouter(%esp),%ebx
        movl %ebx,nb333nf_nouter(%esp)

_nb_kernel333nf_ia32_sse2.nb333nf_outer: 
        movl  nb333nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb333nf_is3(%esp)            ## store is3 

        movl  nb333nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb333nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        movapd %xmm0,%xmm6
        movapd %xmm1,%xmm7

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb333nf_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb333nf_ii3(%esp)

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
        movapd %xmm3,nb333nf_ixO(%esp)
        movapd %xmm4,nb333nf_iyO(%esp)
        movapd %xmm5,nb333nf_izO(%esp)
        movapd %xmm6,nb333nf_ixH1(%esp)
        movapd %xmm7,nb333nf_iyH1(%esp)

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
        movapd %xmm6,nb333nf_izH1(%esp)
        movapd %xmm0,nb333nf_ixH2(%esp)
        movapd %xmm1,nb333nf_iyH2(%esp)
        movapd %xmm2,nb333nf_izH2(%esp)
        movapd %xmm3,nb333nf_ixM(%esp)
        movapd %xmm4,nb333nf_iyM(%esp)
        movapd %xmm5,nb333nf_izM(%esp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb333nf_vctot(%esp)
        movapd %xmm4,nb333nf_Vvdwtot(%esp)

        movl  nb333nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx     ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb333nf_pos(%ebp),%esi
        movl  nb333nf_faction(%ebp),%edi
        movl  nb333nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb333nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb333nf_ninner(%esp),%ecx
        movl  %ecx,nb333nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb333nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel333nf_ia32_sse2.nb333nf_unroll_loop
        jmp   _nb_kernel333nf_ia32_sse2.nb333nf_checksingle
_nb_kernel333nf_ia32_sse2.nb333nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb333nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb333nf_innerjjnr(%esp)             ## advance pointer (unrolled 2) 

        movl nb333nf_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb333nf_iqM(%esp),%xmm3
        mulpd  nb333nf_iqH(%esp),%xmm4
        movapd  %xmm3,nb333nf_qqM(%esp)
        movapd  %xmm4,nb333nf_qqH(%esp)

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movl nb333nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb333nf_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb333nf_ntia(%esp),%edi
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
        movapd %xmm4,nb333nf_c6(%esp)
        movapd %xmm6,nb333nf_c12(%esp)

        movl nb333nf_pos(%ebp),%esi        ## base of pos[] 

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
        movapd nb333nf_ixO(%esp),%xmm4
        movapd nb333nf_iyO(%esp),%xmm5
        movapd nb333nf_izO(%esp),%xmm6

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
        movapd nb333nf_ixH1(%esp),%xmm4
        movapd nb333nf_iyH1(%esp),%xmm5
        movapd nb333nf_izH1(%esp),%xmm6

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
        movapd nb333nf_ixH2(%esp),%xmm3
        movapd nb333nf_iyH2(%esp),%xmm4
        movapd nb333nf_izH2(%esp),%xmm5

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
        movapd nb333nf_iyM(%esp),%xmm3
        movapd nb333nf_izM(%esp),%xmm4
        subpd  %xmm1,%xmm3
        subpd  %xmm2,%xmm4
        movapd nb333nf_ixM(%esp),%xmm2
        subpd  %xmm0,%xmm2

        ## square it 
        mulpd %xmm2,%xmm2
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        addpd %xmm3,%xmm4
        addpd %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## rsqO - put seed in xmm2 
        cvtpd2ps %xmm7,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb333nf_three(%esp),%xmm0
        mulpd   %xmm7,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulpd   nb333nf_half(%esp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm7,%xmm2
        movapd %xmm0,%xmm3
        mulpd %xmm0,%xmm0       ## lu*lu 
        mulpd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb333nf_three(%esp),%xmm0
        subpd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulpd nb333nf_half(%esp),%xmm0   ## rinv 
        mulpd   %xmm0,%xmm7
        movapd  %xmm7,nb333nf_rO(%esp)          ## r in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb333nf_three(%esp),%xmm0
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulpd   nb333nf_half(%esp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm6,%xmm2
        movapd %xmm0,%xmm3
        mulpd %xmm0,%xmm0       ## lu*lu 
        mulpd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb333nf_three(%esp),%xmm0
        subpd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulpd nb333nf_half(%esp),%xmm0   ## rinv 
        mulpd  %xmm0,%xmm6
        movapd %xmm6,nb333nf_rH1(%esp)          ## rH1 

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb333nf_three(%esp),%xmm0
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulpd   nb333nf_half(%esp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm5,%xmm2
        movapd %xmm0,%xmm3
        mulpd %xmm0,%xmm0       ## lu*lu 
        mulpd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb333nf_three(%esp),%xmm0
        subpd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulpd nb333nf_half(%esp),%xmm0   ## rinv 
        mulpd %xmm0,%xmm5
        movapd %xmm5,nb333nf_rH2(%esp)   ## r 

        ## rsqM - seed in xmm2 
        cvtpd2ps %xmm4,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb333nf_three(%esp),%xmm0
        mulpd   %xmm4,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulpd   nb333nf_half(%esp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm4,%xmm2
        movapd %xmm0,%xmm3
        mulpd %xmm0,%xmm0       ## lu*lu 
        mulpd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb333nf_three(%esp),%xmm0
        subpd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulpd nb333nf_half(%esp),%xmm0   ## rinv 
        mulpd %xmm0,%xmm4
        movapd %xmm4,nb333nf_rM(%esp)   ## r 

        ## do O interactions 
        ## rO is still in xmm7 
        mulpd   nb333nf_tsc(%esp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movd %eax,%mm0
        movd %ebx,%mm1
        movl nb333nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        leal  (%eax,%eax,2),%eax ## idx *= 3 (total *=12 now) 
        leal  (%ebx,%ebx,2),%ebx

        ## Dispersion 
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
        ## Dispersion table ready, in xmm4-xmm7                 
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb333nf_c6(%esp),%xmm4
        mulpd  %xmm4,%xmm5       ## Vvdw6               
        addpd  nb333nf_Vvdwtot(%esp),%xmm5
        movapd %xmm5,nb333nf_Vvdwtot(%esp)

        ## Repulsion 
        movlpd 64(%esi,%eax,8),%xmm4    ## Y1
        movlpd 64(%esi,%ebx,8),%xmm3    ## Y2
        movhpd 72(%esi,%eax,8),%xmm4    ## Y1 F1        
        movhpd 72(%esi,%ebx,8),%xmm3    ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 80(%esi,%eax,8),%xmm6    ## G1
        movlpd 80(%esi,%ebx,8),%xmm3    ## G2
        movhpd 88(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 88(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## Dispersion table ready, in xmm4-xmm7                 
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb333nf_c12(%esp),%xmm4
        mulpd  %xmm4,%xmm5 ## Vvdw12 
        addpd  nb333nf_Vvdwtot(%esp),%xmm5
        movapd %xmm5,nb333nf_Vvdwtot(%esp)

        ## Done with O interactions - now H1! 
        movapd nb333nf_rH1(%esp),%xmm7
        mulpd nb333nf_tsc(%esp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb333nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        leal  (%eax,%eax,2),%eax ## idx *= 3 (total *=12 now)   
        leal  (%ebx,%ebx,2),%ebx

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
        movapd nb333nf_qqH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## increment vcoul 
        xorpd  %xmm4,%xmm4
        addpd  nb333nf_vctot(%esp),%xmm5
        movapd %xmm5,nb333nf_vctot(%esp)

        ## H2 interactions 
        movapd nb333nf_rH2(%esp),%xmm7
        mulpd   nb333nf_tsc(%esp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb333nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        leal  (%eax,%eax,2),%eax ## idx *= 3 (total *=12 now)
        leal  (%ebx,%ebx,2),%ebx

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
        movapd nb333nf_qqH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## increment vcoul 
        xorpd  %xmm4,%xmm4
        addpd  nb333nf_vctot(%esp),%xmm5
        movapd %xmm5,nb333nf_vctot(%esp)

        ## M interactions 
        movapd nb333nf_rM(%esp),%xmm7
        mulpd   nb333nf_tsc(%esp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb333nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        leal  (%eax,%eax,2),%eax ## idx *= 3 (total *=12 now)
        leal  (%ebx,%ebx,2),%ebx

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
        movapd nb333nf_qqM(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## increment vcoul 
        xorpd  %xmm4,%xmm4
        addpd  nb333nf_vctot(%esp),%xmm5
        movapd %xmm5,nb333nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $2,nb333nf_innerk(%esp)
        jl    _nb_kernel333nf_ia32_sse2.nb333nf_checksingle
        jmp   _nb_kernel333nf_ia32_sse2.nb333nf_unroll_loop
_nb_kernel333nf_ia32_sse2.nb333nf_checksingle: 
        movl  nb333nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel333nf_ia32_sse2.nb333nf_dosingle
        jmp   _nb_kernel333nf_ia32_sse2.nb333nf_updateouterdata
_nb_kernel333nf_ia32_sse2.nb333nf_dosingle: 
        movl  nb333nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movl nb333nf_charge(%ebp),%esi     ## base of charge[] 
        xorpd %xmm3,%xmm3
        movlpd (%esi,%eax,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb333nf_iqM(%esp),%xmm3
        mulpd  nb333nf_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movapd  %xmm3,nb333nf_qqM(%esp)
        movapd  %xmm4,nb333nf_qqH(%esp)

        movl nb333nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb333nf_vdwparam(%ebp),%esi
        shll %eax
        movl nb333nf_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movapd %xmm4,nb333nf_c6(%esp)
        movapd %xmm6,nb333nf_c12(%esp)

        movl nb333nf_pos(%ebp),%esi        ## base of pos[] 
        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move coords to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb333nf_ixO(%esp),%xmm4
        movapd nb333nf_iyO(%esp),%xmm5
        movapd nb333nf_izO(%esp),%xmm6

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
        movapd nb333nf_ixH1(%esp),%xmm4
        movapd nb333nf_iyH1(%esp),%xmm5
        movapd nb333nf_izH1(%esp),%xmm6

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
        movapd nb333nf_ixH2(%esp),%xmm3
        movapd nb333nf_iyH2(%esp),%xmm4
        movapd nb333nf_izH2(%esp),%xmm5

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
        movapd nb333nf_iyM(%esp),%xmm3
        movapd nb333nf_izM(%esp),%xmm4
        subpd  %xmm1,%xmm3
        subpd  %xmm2,%xmm4
        movapd nb333nf_ixM(%esp),%xmm2
        subpd  %xmm0,%xmm2

        ## square it 
        mulpd %xmm2,%xmm2
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        addpd %xmm3,%xmm4
        addpd %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## start with rsqO - put seed in xmm2 
        cvtsd2ss %xmm7,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb333nf_three(%esp),%xmm0
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulsd   nb333nf_half(%esp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm7,%xmm2
        movapd %xmm0,%xmm3
        mulsd %xmm0,%xmm0       ## lu*lu 
        mulsd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb333nf_three(%esp),%xmm0
        subsd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulsd nb333nf_half(%esp),%xmm0   ## rinv 
        mulsd   %xmm0,%xmm7
        movapd  %xmm7,nb333nf_rO(%esp)          ## r in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb333nf_three(%esp),%xmm0
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulsd   nb333nf_half(%esp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm6,%xmm2
        movapd %xmm0,%xmm3
        mulsd %xmm0,%xmm0       ## lu*lu 
        mulsd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb333nf_three(%esp),%xmm0
        subsd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulsd nb333nf_half(%esp),%xmm0   ## rinv 
        mulsd  %xmm0,%xmm6
        movapd %xmm6,nb333nf_rH1(%esp)          ## rH1 

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb333nf_three(%esp),%xmm0
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulsd   nb333nf_half(%esp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm5,%xmm2
        movapd %xmm0,%xmm3
        mulsd %xmm0,%xmm0       ## lu*lu 
        mulsd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb333nf_three(%esp),%xmm0
        subsd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulsd nb333nf_half(%esp),%xmm0   ## rinv 
        mulsd %xmm0,%xmm5
        movapd %xmm5,nb333nf_rH2(%esp)   ## r 

        ## rsqM - seed in xmm2 
        cvtsd2ss %xmm4,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb333nf_three(%esp),%xmm0
        mulsd   %xmm4,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulsd   nb333nf_half(%esp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm4,%xmm2
        movapd %xmm0,%xmm3
        mulsd %xmm0,%xmm0       ## lu*lu 
        mulsd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb333nf_three(%esp),%xmm0
        subsd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulsd nb333nf_half(%esp),%xmm0   ## rinv 
        mulsd %xmm0,%xmm4
        movapd %xmm4,nb333nf_rM(%esp)   ## r 

        ## do O interactions 
        movd %eax,%mm0
        ## rO is still in xmm7 
        mulsd   nb333nf_tsc(%esp),%xmm7
        cvttsd2si %xmm7,%eax    ## lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm7,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax
        movl nb333nf_VFtab(%ebp),%esi
        leal  (%eax,%eax,2),%eax ## idx *= 3 (total *=12 now)   

        ## Dispersion 
        movsd 32(%esi,%eax,8),%xmm4     ## Y1   
        movsd 40(%esi,%eax,8),%xmm5     ## F1   
        movsd 48(%esi,%eax,8),%xmm6     ## G1   
        movsd 56(%esi,%eax,8),%xmm7     ## H1   
        ## Dispersion table ready, in xmm4-xmm7                 
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb333nf_c6(%esp),%xmm4
        mulsd  %xmm4,%xmm5       ## Vvdw6 

        addsd  nb333nf_Vvdwtot(%esp),%xmm5
        movsd %xmm5,nb333nf_Vvdwtot(%esp)

        ## Repulsion 
        movsd 64(%esi,%eax,8),%xmm4     ## Y1   
        movsd 72(%esi,%eax,8),%xmm5     ## F1   
        movsd 80(%esi,%eax,8),%xmm6     ## G1   
        movsd 88(%esi,%eax,8),%xmm7     ## H1   
        ## Dispersion table ready, in xmm4-xmm7                 
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb333nf_c12(%esp),%xmm4
        mulsd  %xmm4,%xmm5 ## Vvdw12 
        addsd  nb333nf_Vvdwtot(%esp),%xmm5
        movsd %xmm5,nb333nf_Vvdwtot(%esp)

        ## Done with O interactions - now H1! 
        movapd nb333nf_rH1(%esp),%xmm7
        mulpd nb333nf_tsc(%esp),%xmm7
        cvttsd2si %xmm7,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb333nf_VFtab(%ebp),%esi
        leal  (%eax,%eax,2),%eax ## idx *= 3 (total *=12 now)   

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb333nf_qqH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## increment vcoul 
        xorpd  %xmm4,%xmm4
        addsd  nb333nf_vctot(%esp),%xmm5
        movlpd %xmm5,nb333nf_vctot(%esp)

        ##  H2 interactions 
        movapd nb333nf_rH2(%esp),%xmm7
        mulsd   nb333nf_tsc(%esp),%xmm7
        cvttsd2si %xmm7,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb333nf_VFtab(%ebp),%esi
        leal  (%eax,%eax,2),%eax ## idx *= 3 (total *=12 now)   

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   

        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb333nf_qqH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## increment vcoul 
        xorpd  %xmm4,%xmm4
        addsd  nb333nf_vctot(%esp),%xmm5
        movlpd %xmm5,nb333nf_vctot(%esp)

        ## M interactions 
        movapd nb333nf_rM(%esp),%xmm7
        mulsd   nb333nf_tsc(%esp),%xmm7
        cvttsd2si %xmm7,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb333nf_VFtab(%ebp),%esi
        leal  (%eax,%eax,2),%eax ## idx *= 3 (total *=12 now)   

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb333nf_qqM(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## increment vcoul 
        xorpd  %xmm4,%xmm4
        addsd  nb333nf_vctot(%esp),%xmm5
        movlpd %xmm5,nb333nf_vctot(%esp)

_nb_kernel333nf_ia32_sse2.nb333nf_updateouterdata: 
        ## get n from stack
        movl nb333nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb333nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb333nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb333nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb333nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb333nf_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb333nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel333nf_ia32_sse2.nb333nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb333nf_n(%esp)
        jmp _nb_kernel333nf_ia32_sse2.nb333nf_outer
_nb_kernel333nf_ia32_sse2.nb333nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb333nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel333nf_ia32_sse2.nb333nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel333nf_ia32_sse2.nb333nf_threadloop
_nb_kernel333nf_ia32_sse2.nb333nf_end: 
        emms

        movl nb333nf_nouter(%esp),%eax
        movl nb333nf_ninner(%esp),%ebx
        movl nb333nf_outeriter(%ebp),%ecx
        movl nb333nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb333nf_salign(%esp),%eax
        addl %eax,%esp
        addl $476,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret


