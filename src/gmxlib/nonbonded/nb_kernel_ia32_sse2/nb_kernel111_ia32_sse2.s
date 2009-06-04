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


.globl nb_kernel111_ia32_sse2
.globl _nb_kernel111_ia32_sse2
nb_kernel111_ia32_sse2: 
_nb_kernel111_ia32_sse2:        
.set nb111_p_nri, 8
.set nb111_iinr, 12
.set nb111_jindex, 16
.set nb111_jjnr, 20
.set nb111_shift, 24
.set nb111_shiftvec, 28
.set nb111_fshift, 32
.set nb111_gid, 36
.set nb111_pos, 40
.set nb111_faction, 44
.set nb111_charge, 48
.set nb111_p_facel, 52
.set nb111_argkrf, 56
.set nb111_argcrf, 60
.set nb111_Vc, 64
.set nb111_type, 68
.set nb111_p_ntype, 72
.set nb111_vdwparam, 76
.set nb111_Vvdw, 80
.set nb111_p_tabscale, 84
.set nb111_VFtab, 88
.set nb111_invsqrta, 92
.set nb111_dvda, 96
.set nb111_p_gbtabscale, 100
.set nb111_GBtab, 104
.set nb111_p_nthreads, 108
.set nb111_count, 112
.set nb111_mtx, 116
.set nb111_outeriter, 120
.set nb111_inneriter, 124
.set nb111_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb111_ixO, 0
.set nb111_iyO, 16
.set nb111_izO, 32
.set nb111_ixH1, 48
.set nb111_iyH1, 64
.set nb111_izH1, 80
.set nb111_ixH2, 96
.set nb111_iyH2, 112
.set nb111_izH2, 128
.set nb111_iqO, 144
.set nb111_iqH, 160
.set nb111_dxO, 176
.set nb111_dyO, 192
.set nb111_dzO, 208
.set nb111_dxH1, 224
.set nb111_dyH1, 240
.set nb111_dzH1, 256
.set nb111_dxH2, 272
.set nb111_dyH2, 288
.set nb111_dzH2, 304
.set nb111_qqO, 320
.set nb111_qqH, 336
.set nb111_c6, 352
.set nb111_c12, 368
.set nb111_six, 384
.set nb111_twelve, 400
.set nb111_vctot, 416
.set nb111_Vvdwtot, 432
.set nb111_fixO, 448
.set nb111_fiyO, 464
.set nb111_fizO, 480
.set nb111_fixH1, 496
.set nb111_fiyH1, 512
.set nb111_fizH1, 528
.set nb111_fixH2, 544
.set nb111_fiyH2, 560
.set nb111_fizH2, 576
.set nb111_fjx, 592
.set nb111_fjy, 608
.set nb111_fjz, 624
.set nb111_half, 640
.set nb111_three, 656
.set nb111_is3, 672
.set nb111_ii3, 676
.set nb111_ntia, 680
.set nb111_innerjjnr, 684
.set nb111_innerk, 688
.set nb111_n, 692
.set nb111_nn1, 696
.set nb111_nri, 700
.set nb111_nouter, 704
.set nb111_ninner, 708
.set nb111_salign, 712
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $716,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb111_salign(%esp)
        emms

        ## Move args passed by reference to stack
        movl nb111_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb111_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb111_nouter(%esp)
        movl %eax,nb111_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb111_half(%esp)
        movl %ebx,nb111_half+4(%esp)
        movsd nb111_half(%esp),%xmm1
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
        movapd %xmm1,nb111_half(%esp)
        movapd %xmm3,nb111_three(%esp)
        movapd %xmm4,nb111_six(%esp)
        movapd %xmm5,nb111_twelve(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb111_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb111_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        movsd 8(%edx,%ebx,8),%xmm4
        movl nb111_p_facel(%ebp),%esi
        movsd (%esi),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb111_iqO(%esp)
        movapd %xmm4,nb111_iqH(%esp)

        movl  nb111_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb111_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb111_ntia(%esp)
_nb_kernel111_ia32_sse2.nb111_threadloop: 
        movl  nb111_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel111_ia32_sse2.nb111_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel111_ia32_sse2.nb111_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb111_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb111_n(%esp)
        movl %ebx,nb111_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel111_ia32_sse2.nb111_outerstart
        jmp _nb_kernel111_ia32_sse2.nb111_end

_nb_kernel111_ia32_sse2.nb111_outerstart: 
        ## ebx contains number of outer iterations
        addl nb111_nouter(%esp),%ebx
        movl %ebx,nb111_nouter(%esp)

_nb_kernel111_ia32_sse2.nb111_outer: 
        movl  nb111_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb111_is3(%esp)      ## store is3 

        movl  nb111_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb111_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb111_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb111_ii3(%esp)

        addsd (%eax,%ebx,8),%xmm3
        addsd 8(%eax,%ebx,8),%xmm4
        addsd 16(%eax,%ebx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb111_ixO(%esp)
        movapd %xmm4,nb111_iyO(%esp)
        movapd %xmm5,nb111_izO(%esp)

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
        movapd %xmm0,nb111_ixH1(%esp)
        movapd %xmm1,nb111_iyH1(%esp)
        movapd %xmm2,nb111_izH1(%esp)
        movapd %xmm3,nb111_ixH2(%esp)
        movapd %xmm4,nb111_iyH2(%esp)
        movapd %xmm5,nb111_izH2(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb111_vctot(%esp)
        movapd %xmm4,nb111_Vvdwtot(%esp)
        movapd %xmm4,nb111_fixO(%esp)
        movapd %xmm4,nb111_fiyO(%esp)
        movapd %xmm4,nb111_fizO(%esp)
        movapd %xmm4,nb111_fixH1(%esp)
        movapd %xmm4,nb111_fiyH1(%esp)
        movapd %xmm4,nb111_fizH1(%esp)
        movapd %xmm4,nb111_fixH2(%esp)
        movapd %xmm4,nb111_fiyH2(%esp)
        movapd %xmm4,nb111_fizH2(%esp)

        movl  nb111_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb111_pos(%ebp),%esi
        movl  nb111_faction(%ebp),%edi
        movl  nb111_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb111_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb111_ninner(%esp),%ecx
        movl  %ecx,nb111_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb111_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel111_ia32_sse2.nb111_unroll_loop
        jmp   _nb_kernel111_ia32_sse2.nb111_checksingle
_nb_kernel111_ia32_sse2.nb111_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb111_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb111_innerjjnr(%esp)                   ## advance pointer (unrolled 2) 

        movl nb111_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb111_iqO(%esp),%xmm3
        mulpd  nb111_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movapd  %xmm3,nb111_qqO(%esp)
        movapd  %xmm4,nb111_qqH(%esp)

        movl nb111_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb111_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb111_ntia(%esp),%edi
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
        movapd %xmm4,nb111_c6(%esp)
        movapd %xmm6,nb111_c12(%esp)

        movl nb111_pos(%ebp),%esi        ## base of pos[] 

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
        movapd nb111_ixO(%esp),%xmm4
        movapd nb111_iyO(%esp),%xmm5
        movapd nb111_izO(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb111_dxO(%esp)
        movapd %xmm5,nb111_dyO(%esp)
        movapd %xmm6,nb111_dzO(%esp)
        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb111_ixH1(%esp),%xmm4
        movapd nb111_iyH1(%esp),%xmm5
        movapd nb111_izH1(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb111_dxH1(%esp)
        movapd %xmm5,nb111_dyH1(%esp)
        movapd %xmm6,nb111_dzH1(%esp)
        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm6
        addpd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb111_ixH2(%esp),%xmm3
        movapd nb111_iyH2(%esp),%xmm4
        movapd nb111_izH2(%esp),%xmm5

        ## calc dr 
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5

        ## store dr 
        movapd %xmm3,nb111_dxH2(%esp)
        movapd %xmm4,nb111_dyH2(%esp)
        movapd %xmm5,nb111_dzH2(%esp)
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
        movapd  nb111_three(%esp),%xmm4
        mulpd   %xmm7,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb111_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb111_three(%esp),%xmm4
        subpd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb111_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvO in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb111_three(%esp),%xmm4
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb111_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb111_three(%esp),%xmm4
        subpd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb111_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb111_three(%esp),%xmm4
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb111_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb111_three(%esp),%xmm4
        subpd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb111_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do O interactions 
        movapd  %xmm7,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq 
        movapd %xmm4,%xmm1
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulpd  nb111_qqO(%esp),%xmm7    ## xmm7=vcoul 

        mulpd  nb111_c6(%esp),%xmm1
        mulpd  nb111_c12(%esp),%xmm2
        movapd %xmm2,%xmm3
        subpd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addpd  nb111_Vvdwtot(%esp),%xmm3
        mulpd  nb111_six(%esp),%xmm1
        mulpd  nb111_twelve(%esp),%xmm2
        subpd  %xmm1,%xmm2
        addpd  %xmm7,%xmm2
        mulpd  %xmm2,%xmm4      ## total fsO in xmm4 

        addpd  nb111_vctot(%esp),%xmm7

        movapd %xmm3,nb111_Vvdwtot(%esp)
        movapd %xmm7,nb111_vctot(%esp)

        movapd nb111_dxO(%esp),%xmm0
        movapd nb111_dyO(%esp),%xmm1
        movapd nb111_dzO(%esp),%xmm2
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        ## update O forces 
        movapd nb111_fixO(%esp),%xmm3
        movapd nb111_fiyO(%esp),%xmm4
        movapd nb111_fizO(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb111_fixO(%esp)
        movapd %xmm4,nb111_fiyO(%esp)
        movapd %xmm7,nb111_fizO(%esp)
        ## update j forces with water O 
        movapd %xmm0,nb111_fjx(%esp)
        movapd %xmm1,nb111_fjy(%esp)
        movapd %xmm2,nb111_fjz(%esp)

        ## H1 interactions 
        movapd  %xmm6,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        mulpd  nb111_qqH(%esp),%xmm6    ## xmm6=vcoul 
        mulpd  %xmm6,%xmm4              ## total fsH1 in xmm4 

        addpd  nb111_vctot(%esp),%xmm6

        movapd nb111_dxH1(%esp),%xmm0
        movapd nb111_dyH1(%esp),%xmm1
        movapd nb111_dzH1(%esp),%xmm2
        movapd %xmm6,nb111_vctot(%esp)
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb111_fixH1(%esp),%xmm3
        movapd nb111_fiyH1(%esp),%xmm4
        movapd nb111_fizH1(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb111_fixH1(%esp)
        movapd %xmm4,nb111_fiyH1(%esp)
        movapd %xmm7,nb111_fizH1(%esp)
        ## update j forces with water H1 
        addpd  nb111_fjx(%esp),%xmm0
        addpd  nb111_fjy(%esp),%xmm1
        addpd  nb111_fjz(%esp),%xmm2
        movapd %xmm0,nb111_fjx(%esp)
        movapd %xmm1,nb111_fjy(%esp)
        movapd %xmm2,nb111_fjz(%esp)

        ## H2 interactions 
        movapd  %xmm5,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        mulpd  nb111_qqH(%esp),%xmm5    ## xmm5=vcoul 
        mulpd  %xmm5,%xmm4              ## total fsH1 in xmm4 

        addpd  nb111_vctot(%esp),%xmm5

        movapd nb111_dxH2(%esp),%xmm0
        movapd nb111_dyH2(%esp),%xmm1
        movapd nb111_dzH2(%esp),%xmm2
        movapd %xmm5,nb111_vctot(%esp)
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        ## update H2 forces 
        movapd nb111_fixH2(%esp),%xmm3
        movapd nb111_fiyH2(%esp),%xmm4
        movapd nb111_fizH2(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb111_fixH2(%esp)
        movapd %xmm4,nb111_fiyH2(%esp)
        movapd %xmm7,nb111_fizH2(%esp)

        movl nb111_faction(%ebp),%edi
        ## update j forces 
        addpd  nb111_fjx(%esp),%xmm0
        addpd  nb111_fjy(%esp),%xmm1
        addpd  nb111_fjz(%esp),%xmm2
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
        subl $2,nb111_innerk(%esp)
        jl   _nb_kernel111_ia32_sse2.nb111_checksingle
        jmp  _nb_kernel111_ia32_sse2.nb111_unroll_loop
_nb_kernel111_ia32_sse2.nb111_checksingle: 
        movl  nb111_innerk(%esp),%edx
        andl  $1,%edx
        jnz  _nb_kernel111_ia32_sse2.nb111_dosingle
        jmp  _nb_kernel111_ia32_sse2.nb111_updateouterdata
_nb_kernel111_ia32_sse2.nb111_dosingle: 
        movl  nb111_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb111_innerjjnr(%esp)

        movl nb111_charge(%ebp),%esi     ## base of charge[] 

        xorpd %xmm3,%xmm3
        movlpd (%esi,%eax,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb111_iqO(%esp),%xmm3
        mulpd  nb111_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 

        movapd  %xmm3,nb111_qqO(%esp)
        movapd  %xmm4,nb111_qqH(%esp)

        movl nb111_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb111_vdwparam(%ebp),%esi
        shll %eax
        movl nb111_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movapd %xmm4,nb111_c6(%esp)
        movapd %xmm6,nb111_c12(%esp)

        movl nb111_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb111_ixO(%esp),%xmm4
        movapd nb111_iyO(%esp),%xmm5
        movapd nb111_izO(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb111_dxO(%esp)
        movapd %xmm5,nb111_dyO(%esp)
        movapd %xmm6,nb111_dzO(%esp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb111_ixH1(%esp),%xmm4
        movapd nb111_iyH1(%esp),%xmm5
        movapd nb111_izH1(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb111_dxH1(%esp)
        movapd %xmm5,nb111_dyH1(%esp)
        movapd %xmm6,nb111_dzH1(%esp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm6
        addsd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb111_ixH2(%esp),%xmm3
        movapd nb111_iyH2(%esp),%xmm4
        movapd nb111_izH2(%esp),%xmm5

        ## calc dr 
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5

        ## store dr 
        movapd %xmm3,nb111_dxH2(%esp)
        movapd %xmm4,nb111_dyH2(%esp)
        movapd %xmm5,nb111_dzH2(%esp)
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
        movapd  nb111_three(%esp),%xmm4
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb111_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb111_three(%esp),%xmm4
        subsd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb111_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvO in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb111_three(%esp),%xmm4
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb111_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb111_three(%esp),%xmm4
        subsd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb111_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb111_three(%esp),%xmm4
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb111_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb111_three(%esp),%xmm4
        subsd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb111_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do O interactions 
        movapd  %xmm7,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq 
        movapd %xmm4,%xmm1
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulsd  nb111_qqO(%esp),%xmm7    ## xmm7=vcoul 

        mulsd  nb111_c6(%esp),%xmm1
        mulsd  nb111_c12(%esp),%xmm2
        movapd %xmm2,%xmm3
        subsd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addsd  nb111_Vvdwtot(%esp),%xmm3
        mulsd  nb111_six(%esp),%xmm1
        mulsd  nb111_twelve(%esp),%xmm2
        subsd  %xmm1,%xmm2
        addsd  %xmm7,%xmm2
        mulsd  %xmm2,%xmm4      ## total fsO in xmm4 

        addsd  nb111_vctot(%esp),%xmm7

        movsd %xmm3,nb111_Vvdwtot(%esp)
        movsd %xmm7,nb111_vctot(%esp)

        movapd nb111_dxO(%esp),%xmm0
        movapd nb111_dyO(%esp),%xmm1
        movapd nb111_dzO(%esp),%xmm2
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update O forces 
        movapd nb111_fixO(%esp),%xmm3
        movapd nb111_fiyO(%esp),%xmm4
        movapd nb111_fizO(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movsd %xmm3,nb111_fixO(%esp)
        movsd %xmm4,nb111_fiyO(%esp)
        movsd %xmm7,nb111_fizO(%esp)
        ## update j forces with water O 
        movsd %xmm0,nb111_fjx(%esp)
        movsd %xmm1,nb111_fjy(%esp)
        movsd %xmm2,nb111_fjz(%esp)

        ## H1 interactions 
        movapd  %xmm6,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        mulsd  nb111_qqH(%esp),%xmm6    ## xmm6=vcoul 
        mulsd  %xmm6,%xmm4              ## total fsH1 in xmm4 

        addsd  nb111_vctot(%esp),%xmm6

        movapd nb111_dxH1(%esp),%xmm0
        movapd nb111_dyH1(%esp),%xmm1
        movapd nb111_dzH1(%esp),%xmm2
        movsd %xmm6,nb111_vctot(%esp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb111_fixH1(%esp),%xmm3
        movapd nb111_fiyH1(%esp),%xmm4
        movapd nb111_fizH1(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movsd %xmm3,nb111_fixH1(%esp)
        movsd %xmm4,nb111_fiyH1(%esp)
        movsd %xmm7,nb111_fizH1(%esp)
        ## update j forces with water H1 
        addsd  nb111_fjx(%esp),%xmm0
        addsd  nb111_fjy(%esp),%xmm1
        addsd  nb111_fjz(%esp),%xmm2
        movsd %xmm0,nb111_fjx(%esp)
        movsd %xmm1,nb111_fjy(%esp)
        movsd %xmm2,nb111_fjz(%esp)

        ## H2 interactions 
        movapd  %xmm5,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        mulsd  nb111_qqH(%esp),%xmm5    ## xmm5=vcoul 
        mulsd  %xmm5,%xmm4              ## total fsH1 in xmm4 

        addsd  nb111_vctot(%esp),%xmm5

        movapd nb111_dxH2(%esp),%xmm0
        movapd nb111_dyH2(%esp),%xmm1
        movapd nb111_dzH2(%esp),%xmm2
        movsd %xmm5,nb111_vctot(%esp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H2 forces 
        movapd nb111_fixH2(%esp),%xmm3
        movapd nb111_fiyH2(%esp),%xmm4
        movapd nb111_fizH2(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movsd %xmm3,nb111_fixH2(%esp)
        movsd %xmm4,nb111_fiyH2(%esp)
        movsd %xmm7,nb111_fizH2(%esp)

        movl nb111_faction(%ebp),%edi
        ## update j forces 
        addsd  nb111_fjx(%esp),%xmm0
        addsd  nb111_fjy(%esp),%xmm1
        addsd  nb111_fjz(%esp),%xmm2
        movlpd (%edi,%eax,8),%xmm3
        movlpd 8(%edi,%eax,8),%xmm4
        movlpd 16(%edi,%eax,8),%xmm5
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        movlpd %xmm3,(%edi,%eax,8)
        movlpd %xmm4,8(%edi,%eax,8)
        movlpd %xmm5,16(%edi,%eax,8)

_nb_kernel111_ia32_sse2.nb111_updateouterdata: 
        movl  nb111_ii3(%esp),%ecx
        movl  nb111_faction(%ebp),%edi
        movl  nb111_fshift(%ebp),%esi
        movl  nb111_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movapd nb111_fixO(%esp),%xmm0
        movapd nb111_fiyO(%esp),%xmm1
        movapd nb111_fizO(%esp),%xmm2

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
        movapd nb111_fixH1(%esp),%xmm0
        movapd nb111_fiyH1(%esp),%xmm1
        movapd nb111_fizH1(%esp),%xmm2

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
        movapd nb111_fixH2(%esp),%xmm0
        movapd nb111_fiyH2(%esp),%xmm1
        movapd nb111_fizH2(%esp),%xmm2

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
        movl nb111_n(%esp),%esi
        ## get group index for i particle 
        movl  nb111_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb111_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb111_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb111_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb111_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

       ## finish if last 
        movl nb111_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel111_ia32_sse2.nb111_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb111_n(%esp)
        jmp _nb_kernel111_ia32_sse2.nb111_outer
_nb_kernel111_ia32_sse2.nb111_outerend: 
        ## check if more outer neighborlists remain
        movl  nb111_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel111_ia32_sse2.nb111_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel111_ia32_sse2.nb111_threadloop
_nb_kernel111_ia32_sse2.nb111_end: 
        emms

        movl nb111_nouter(%esp),%eax
        movl nb111_ninner(%esp),%ebx
        movl nb111_outeriter(%ebp),%ecx
        movl nb111_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb111_salign(%esp),%eax
        addl %eax,%esp
        addl $716,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret


.globl nb_kernel111nf_ia32_sse2
.globl _nb_kernel111nf_ia32_sse2
nb_kernel111nf_ia32_sse2:       
_nb_kernel111nf_ia32_sse2:      
.set nb111nf_p_nri, 8
.set nb111nf_iinr, 12
.set nb111nf_jindex, 16
.set nb111nf_jjnr, 20
.set nb111nf_shift, 24
.set nb111nf_shiftvec, 28
.set nb111nf_fshift, 32
.set nb111nf_gid, 36
.set nb111nf_pos, 40
.set nb111nf_faction, 44
.set nb111nf_charge, 48
.set nb111nf_p_facel, 52
.set nb111nf_argkrf, 56
.set nb111nf_argcrf, 60
.set nb111nf_Vc, 64
.set nb111nf_type, 68
.set nb111nf_p_ntype, 72
.set nb111nf_vdwparam, 76
.set nb111nf_Vvdw, 80
.set nb111nf_p_tabscale, 84
.set nb111nf_VFtab, 88
.set nb111nf_invsqrta, 92
.set nb111nf_dvda, 96
.set nb111nf_p_gbtabscale, 100
.set nb111nf_GBtab, 104
.set nb111nf_p_nthreads, 108
.set nb111nf_count, 112
.set nb111nf_mtx, 116
.set nb111nf_outeriter, 120
.set nb111nf_inneriter, 124
.set nb111nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb111nf_ixO, 0
.set nb111nf_iyO, 16
.set nb111nf_izO, 32
.set nb111nf_ixH1, 48
.set nb111nf_iyH1, 64
.set nb111nf_izH1, 80
.set nb111nf_ixH2, 96
.set nb111nf_iyH2, 112
.set nb111nf_izH2, 128
.set nb111nf_iqO, 144
.set nb111nf_iqH, 160
.set nb111nf_qqO, 176
.set nb111nf_qqH, 192
.set nb111nf_c6, 208
.set nb111nf_c12, 224
.set nb111nf_vctot, 240
.set nb111nf_Vvdwtot, 256
.set nb111nf_half, 272
.set nb111nf_three, 288
.set nb111nf_is3, 304
.set nb111nf_ii3, 308
.set nb111nf_ntia, 312
.set nb111nf_innerjjnr, 316
.set nb111nf_innerk, 320
.set nb111nf_n, 324
.set nb111nf_nn1, 328
.set nb111nf_nri, 332
.set nb111nf_nouter, 336
.set nb111nf_ninner, 340
.set nb111nf_salign, 344
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $348,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb111nf_salign(%esp)
        emms

        ## Move args passed by reference to stack
        movl nb111nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb111nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb111nf_nouter(%esp)
        movl %eax,nb111nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb111nf_half(%esp)
        movl %ebx,nb111nf_half+4(%esp)
        movsd nb111nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb111nf_half(%esp)
        movapd %xmm3,nb111nf_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb111nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb111nf_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        movsd 8(%edx,%ebx,8),%xmm4
        movl nb111nf_p_facel(%ebp),%esi
        movsd (%esi),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb111nf_iqO(%esp)
        movapd %xmm4,nb111nf_iqH(%esp)

        movl  nb111nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb111nf_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb111nf_ntia(%esp)
_nb_kernel111nf_ia32_sse2.nb111nf_threadloop: 
        movl  nb111nf_count(%ebp),%esi          ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel111nf_ia32_sse2.nb111nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel111nf_ia32_sse2.nb111nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb111nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb111nf_n(%esp)
        movl %ebx,nb111nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel111nf_ia32_sse2.nb111nf_outerstart
        jmp _nb_kernel111nf_ia32_sse2.nb111nf_end

_nb_kernel111nf_ia32_sse2.nb111nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb111nf_nouter(%esp),%ebx
        movl %ebx,nb111nf_nouter(%esp)

_nb_kernel111nf_ia32_sse2.nb111nf_outer: 
        movl  nb111nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 

        movl  nb111nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb111nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb111nf_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb111nf_ii3(%esp)

        addsd (%eax,%ebx,8),%xmm3
        addsd 8(%eax,%ebx,8),%xmm4
        addsd 16(%eax,%ebx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb111nf_ixO(%esp)
        movapd %xmm4,nb111nf_iyO(%esp)
        movapd %xmm5,nb111nf_izO(%esp)

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
        movapd %xmm0,nb111nf_ixH1(%esp)
        movapd %xmm1,nb111nf_iyH1(%esp)
        movapd %xmm2,nb111nf_izH1(%esp)
        movapd %xmm3,nb111nf_ixH2(%esp)
        movapd %xmm4,nb111nf_iyH2(%esp)
        movapd %xmm5,nb111nf_izH2(%esp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb111nf_vctot(%esp)
        movapd %xmm4,nb111nf_Vvdwtot(%esp)

        movl  nb111nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb111nf_pos(%ebp),%esi
        movl  nb111nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb111nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb111nf_ninner(%esp),%ecx
        movl  %ecx,nb111nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb111nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel111nf_ia32_sse2.nb111nf_unroll_loop
        jmp   _nb_kernel111nf_ia32_sse2.nb111nf_checksingle
_nb_kernel111nf_ia32_sse2.nb111nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb111nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb111nf_innerjjnr(%esp)                 ## advance pointer (unrolled 2) 

        movl nb111nf_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb111nf_iqO(%esp),%xmm3
        mulpd  nb111nf_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movapd  %xmm3,nb111nf_qqO(%esp)
        movapd  %xmm4,nb111nf_qqH(%esp)

        movl nb111nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb111nf_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb111nf_ntia(%esp),%edi
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
        movapd %xmm4,nb111nf_c6(%esp)
        movapd %xmm6,nb111nf_c12(%esp)

        movl nb111nf_pos(%ebp),%esi        ## base of pos[] 

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
        movapd nb111nf_ixO(%esp),%xmm4
        movapd nb111nf_iyO(%esp),%xmm5
        movapd nb111nf_izO(%esp),%xmm6

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
        movapd nb111nf_ixH1(%esp),%xmm4
        movapd nb111nf_iyH1(%esp),%xmm5
        movapd nb111nf_izH1(%esp),%xmm6

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
        movapd nb111nf_ixH2(%esp),%xmm3
        movapd nb111nf_iyH2(%esp),%xmm4
        movapd nb111nf_izH2(%esp),%xmm5

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
        movapd  nb111nf_three(%esp),%xmm4
        mulpd   %xmm7,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb111nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb111nf_three(%esp),%xmm4
        subpd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb111nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvO in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb111nf_three(%esp),%xmm4
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb111nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb111nf_three(%esp),%xmm4
        subpd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb111nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb111nf_three(%esp),%xmm4
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb111nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb111nf_three(%esp),%xmm4
        subpd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb111nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do O interactions 
        movapd  %xmm7,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq 
        movapd %xmm4,%xmm1
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulpd  nb111nf_qqO(%esp),%xmm7          ## xmm7=vcoul 

        mulpd  nb111nf_c6(%esp),%xmm1
        mulpd  nb111nf_c12(%esp),%xmm2
        movapd %xmm2,%xmm3
        subpd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addpd  nb111nf_Vvdwtot(%esp),%xmm3
        addpd  nb111nf_vctot(%esp),%xmm7
        movapd %xmm3,nb111nf_Vvdwtot(%esp)
        movapd %xmm7,nb111nf_vctot(%esp)

        ## H1 interactions 
        mulpd  nb111nf_qqH(%esp),%xmm6          ## xmm6=vcoul 
        addpd  nb111nf_vctot(%esp),%xmm6
        movapd %xmm6,nb111nf_vctot(%esp)

        ## H2 interactions 
        mulpd  nb111nf_qqH(%esp),%xmm5          ## xmm5=vcoul 
        addpd  nb111nf_vctot(%esp),%xmm5
        movapd %xmm5,nb111nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $2,nb111nf_innerk(%esp)
        jl    _nb_kernel111nf_ia32_sse2.nb111nf_checksingle
        jmp   _nb_kernel111nf_ia32_sse2.nb111nf_unroll_loop
_nb_kernel111nf_ia32_sse2.nb111nf_checksingle: 
        movl  nb111nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel111nf_ia32_sse2.nb111nf_dosingle
        jmp   _nb_kernel111nf_ia32_sse2.nb111nf_updateouterdata
_nb_kernel111nf_ia32_sse2.nb111nf_dosingle: 
        movl  nb111nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb111nf_innerjjnr(%esp)

        movl nb111nf_charge(%ebp),%esi     ## base of charge[] 

        xorpd %xmm3,%xmm3
        movlpd (%esi,%eax,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb111nf_iqO(%esp),%xmm3
        mulpd  nb111nf_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 

        movapd  %xmm3,nb111nf_qqO(%esp)
        movapd  %xmm4,nb111nf_qqH(%esp)

        movl nb111nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb111nf_vdwparam(%ebp),%esi
        shll %eax
        movl nb111nf_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movapd %xmm4,nb111nf_c6(%esp)
        movapd %xmm6,nb111nf_c12(%esp)

        movl nb111nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb111nf_ixO(%esp),%xmm4
        movapd nb111nf_iyO(%esp),%xmm5
        movapd nb111nf_izO(%esp),%xmm6

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
        movapd nb111nf_ixH1(%esp),%xmm4
        movapd nb111nf_iyH1(%esp),%xmm5
        movapd nb111nf_izH1(%esp),%xmm6

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
        movapd nb111nf_ixH2(%esp),%xmm3
        movapd nb111nf_iyH2(%esp),%xmm4
        movapd nb111nf_izH2(%esp),%xmm5

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
        movapd  nb111nf_three(%esp),%xmm4
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb111nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb111nf_three(%esp),%xmm4
        subsd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb111nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvO in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb111nf_three(%esp),%xmm4
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb111nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb111nf_three(%esp),%xmm4
        subsd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb111nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb111nf_three(%esp),%xmm4
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb111nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb111nf_three(%esp),%xmm4
        subsd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb111nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do O interactions 
        movapd  %xmm7,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq 
        movapd %xmm4,%xmm1
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulsd  nb111nf_qqO(%esp),%xmm7          ## xmm7=vcoul 

        mulsd  nb111nf_c6(%esp),%xmm1
        mulsd  nb111nf_c12(%esp),%xmm2
        movapd %xmm2,%xmm3
        subsd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addsd  nb111nf_Vvdwtot(%esp),%xmm3
        addsd  nb111nf_vctot(%esp),%xmm7
        movsd %xmm3,nb111nf_Vvdwtot(%esp)
        movsd %xmm7,nb111nf_vctot(%esp)

        ## H1 interactions 
        mulsd  nb111nf_qqH(%esp),%xmm6          ## xmm6=vcoul 
        addsd  nb111nf_vctot(%esp),%xmm6
        movsd %xmm6,nb111nf_vctot(%esp)

        ## H2 interactions 
        mulsd  nb111nf_qqH(%esp),%xmm5          ## xmm5=vcoul 
        addsd  nb111nf_vctot(%esp),%xmm5
        movsd %xmm5,nb111nf_vctot(%esp)

_nb_kernel111nf_ia32_sse2.nb111nf_updateouterdata: 
        ## get n from stack
        movl nb111nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb111nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb111nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb111nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb111nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb111nf_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb111nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel111nf_ia32_sse2.nb111nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb111nf_n(%esp)
        jmp _nb_kernel111nf_ia32_sse2.nb111nf_outer
_nb_kernel111nf_ia32_sse2.nb111nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb111nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel111nf_ia32_sse2.nb111nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel111nf_ia32_sse2.nb111nf_threadloop
_nb_kernel111nf_ia32_sse2.nb111nf_end: 
        emms

        movl nb111nf_nouter(%esp),%eax
        movl nb111nf_ninner(%esp),%ebx
        movl nb111nf_outeriter(%ebp),%ecx
        movl nb111nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb111nf_salign(%esp),%eax
        addl %eax,%esp
        addl $348,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret


