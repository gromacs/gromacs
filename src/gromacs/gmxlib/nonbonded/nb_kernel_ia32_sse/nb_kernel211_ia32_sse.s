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


.globl nb_kernel211_ia32_sse
.globl _nb_kernel211_ia32_sse
nb_kernel211_ia32_sse:  
_nb_kernel211_ia32_sse: 
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
        ## bottom of stack is cache-aligned for sse use 
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
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb211_half(%esp)
        movss nb211_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm3,%xmm4
        addps  %xmm4,%xmm4      ## 6.0
        movaps %xmm4,%xmm5
        addps  %xmm5,%xmm5      ## constant 12.0
        movaps %xmm1,nb211_half(%esp)
        movaps %xmm2,nb211_two(%esp)
        movaps %xmm3,nb211_three(%esp)
        movaps %xmm4,nb211_six(%esp)
        movaps %xmm5,nb211_twelve(%esp)

        movl nb211_argkrf(%ebp),%esi
        movl nb211_argcrf(%ebp),%edi
        movss (%esi),%xmm5
        movss (%edi),%xmm6
        shufps $0,%xmm5,%xmm5
        shufps $0,%xmm6,%xmm6
        movaps %xmm5,nb211_krf(%esp)
        movaps %xmm6,nb211_crf(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb211_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb211_charge(%ebp),%edx
        movss (%edx,%ebx,4),%xmm3
        movss 4(%edx,%ebx,4),%xmm4
        movl nb211_p_facel(%ebp),%esi
        movss (%esi),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb211_iqO(%esp)
        movaps %xmm4,nb211_iqH(%esp)

        movl  nb211_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb211_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb211_ntia(%esp)

_nb_kernel211_ia32_sse.nb211_threadloop: 
        movl  nb211_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel211_ia32_sse.nb211_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel211_ia32_sse.nb211_spinlock

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
        jg  _nb_kernel211_ia32_sse.nb211_outerstart
        jmp _nb_kernel211_ia32_sse.nb211_end

_nb_kernel211_ia32_sse.nb211_outerstart: 
        ## ebx contains number of outer iterations
        addl nb211_nouter(%esp),%ebx
        movl %ebx,nb211_nouter(%esp)

_nb_kernel211_ia32_sse.nb211_outer: 
        movl  nb211_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb211_is3(%esp)      ## store is3 

        movl  nb211_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb211_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb211_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb211_ii3(%esp)

        addss (%eax,%ebx,4),%xmm3
        addss 4(%eax,%ebx,4),%xmm4
        addss 8(%eax,%ebx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb211_ixO(%esp)
        movaps %xmm4,nb211_iyO(%esp)
        movaps %xmm5,nb211_izO(%esp)

        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 12(%eax,%ebx,4),%xmm0
        addss 16(%eax,%ebx,4),%xmm1
        addss 20(%eax,%ebx,4),%xmm2
        addss 24(%eax,%ebx,4),%xmm3
        addss 28(%eax,%ebx,4),%xmm4
        addss 32(%eax,%ebx,4),%xmm5

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm0,nb211_ixH1(%esp)
        movaps %xmm1,nb211_iyH1(%esp)
        movaps %xmm2,nb211_izH1(%esp)
        movaps %xmm3,nb211_ixH2(%esp)
        movaps %xmm4,nb211_iyH2(%esp)
        movaps %xmm5,nb211_izH2(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb211_vctot(%esp)
        movaps %xmm4,nb211_Vvdwtot(%esp)
        movaps %xmm4,nb211_fixO(%esp)
        movaps %xmm4,nb211_fiyO(%esp)
        movaps %xmm4,nb211_fizO(%esp)
        movaps %xmm4,nb211_fixH1(%esp)
        movaps %xmm4,nb211_fiyH1(%esp)
        movaps %xmm4,nb211_fizH1(%esp)
        movaps %xmm4,nb211_fixH2(%esp)
        movaps %xmm4,nb211_fiyH2(%esp)
        movaps %xmm4,nb211_fizH2(%esp)

        movl  nb211_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb211_pos(%ebp),%esi
        movl  nb211_faction(%ebp),%edi
        movl  nb211_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb211_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb211_ninner(%esp),%ecx
        movl  %ecx,nb211_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb211_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel211_ia32_sse.nb211_unroll_loop
        jmp   _nb_kernel211_ia32_sse.nb211_odd_inner
_nb_kernel211_ia32_sse.nb211_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb211_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 

        addl $16,nb211_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb211_charge(%ebp),%esi     ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4           ## and in xmm4 
        mulps  nb211_iqO(%esp),%xmm3
        mulps  nb211_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps  %xmm3,nb211_qqO(%esp)
        movaps  %xmm4,nb211_qqH(%esp)

        movl nb211_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb211_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb211_ntia(%esp),%edi
        addl %edi,%eax
        addl %edi,%ebx
        addl %edi,%ecx
        addl %edi,%edx

        movlps (%esi,%eax,4),%xmm6
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm6
        movhps (%esi,%edx,4),%xmm7

        movaps %xmm6,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm6 ## constant 11011101

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movd  %mm2,%ecx
        movd  %mm3,%edx

        movaps %xmm4,nb211_c6(%esp)
        movaps %xmm6,nb211_c12(%esp)

        movl nb211_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx     ## replace jnr with j3 
        leal  (%edx,%edx,2),%edx

        ## move four coordinates to xmm0-xmm2   
        movlps (%esi,%eax,4),%xmm4
        movlps (%esi,%ecx,4),%xmm5
        movss 8(%esi,%eax,4),%xmm2
        movss 8(%esi,%ecx,4),%xmm6

        movhps (%esi,%ebx,4),%xmm4
        movhps (%esi,%edx,4),%xmm5

        movss 8(%esi,%ebx,4),%xmm0
        movss 8(%esi,%edx,4),%xmm1

        shufps $0,%xmm0,%xmm2
        shufps $0,%xmm1,%xmm6

        movaps %xmm4,%xmm0
        movaps %xmm4,%xmm1

        shufps $136,%xmm6,%xmm2 ## constant 10001000

        shufps $136,%xmm5,%xmm0 ## constant 10001000
        shufps $221,%xmm5,%xmm1 ## constant 11011101            

        ## move ixO-izO to xmm4-xmm6 
        movaps nb211_ixO(%esp),%xmm4
        movaps nb211_iyO(%esp),%xmm5
        movaps nb211_izO(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb211_dxO(%esp)
        movaps %xmm5,nb211_dyO(%esp)
        movaps %xmm6,nb211_dzO(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        movaps %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movaps nb211_ixH1(%esp),%xmm4
        movaps nb211_iyH1(%esp),%xmm5
        movaps nb211_izH1(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb211_dxH1(%esp)
        movaps %xmm5,nb211_dyH1(%esp)
        movaps %xmm6,nb211_dzH1(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm6
        addps %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movaps nb211_ixH2(%esp),%xmm3
        movaps nb211_iyH2(%esp),%xmm4
        movaps nb211_izH2(%esp),%xmm5

        ## calc dr 
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        ## store dr 
        movaps %xmm3,nb211_dxH2(%esp)
        movaps %xmm4,nb211_dyH2(%esp)
        movaps %xmm5,nb211_dzH2(%esp)
        ## square it 
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm4,%xmm5
        addps %xmm3,%xmm5
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        movaps %xmm5,%xmm0
        movaps %xmm6,%xmm1
        movaps %xmm7,%xmm2

        mulps  nb211_krf(%esp),%xmm0
        mulps  nb211_krf(%esp),%xmm1
        mulps  nb211_krf(%esp),%xmm2

        movaps %xmm0,nb211_krsqH2(%esp)
        movaps %xmm1,nb211_krsqH1(%esp)
        movaps %xmm2,nb211_krsqO(%esp)

        ## start with rsqO - seed in xmm2       
        rsqrtps %xmm7,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb211_three(%esp),%xmm4
        mulps   %xmm7,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb211_half(%esp),%xmm4
        movaps  %xmm4,%xmm7     ## rinvO in xmm7 
        ## rsqH1 - seed in xmm2 
        rsqrtps %xmm6,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb211_three(%esp),%xmm4
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb211_half(%esp),%xmm4
        movaps  %xmm4,%xmm6     ## rinvH1 in xmm6 
        ## rsqH2 - seed in xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb211_three(%esp),%xmm4
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb211_half(%esp),%xmm4
        movaps  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do O interactions 
        movaps  %xmm7,%xmm4
        mulps   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq 
        movaps %xmm4,%xmm1
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  nb211_c6(%esp),%xmm1
        mulps  nb211_c12(%esp),%xmm2
        movaps %xmm2,%xmm3
        subps  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addps  nb211_Vvdwtot(%esp),%xmm3
        mulps  nb211_six(%esp),%xmm1
        mulps  nb211_twelve(%esp),%xmm2
        subps  %xmm1,%xmm2      ## nb part of fs  

        movaps %xmm7,%xmm0
        movaps nb211_krsqO(%esp),%xmm1
        addps  %xmm1,%xmm0
        mulps  nb211_two(%esp),%xmm1
        subps  nb211_crf(%esp),%xmm0   ## xmm0=rinv+ krsq-crf 
        subps  %xmm1,%xmm7
        mulps  nb211_qqO(%esp),%xmm0
        mulps  nb211_qqO(%esp),%xmm7
        addps  %xmm7,%xmm2

        mulps  %xmm2,%xmm4      ## total fsO in xmm4 

        addps  nb211_vctot(%esp),%xmm0
        movaps %xmm3,nb211_Vvdwtot(%esp)
        movaps %xmm0,nb211_vctot(%esp)

        movaps nb211_dxO(%esp),%xmm0
        movaps nb211_dyO(%esp),%xmm1
        movaps nb211_dzO(%esp),%xmm2
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2

        ## update O forces 
        movaps nb211_fixO(%esp),%xmm3
        movaps nb211_fiyO(%esp),%xmm4
        movaps nb211_fizO(%esp),%xmm7
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm7
        movaps %xmm3,nb211_fixO(%esp)
        movaps %xmm4,nb211_fiyO(%esp)
        movaps %xmm7,nb211_fizO(%esp)
        ## update j forces with water O 
        movaps %xmm0,nb211_fjx(%esp)
        movaps %xmm1,nb211_fjy(%esp)
        movaps %xmm2,nb211_fjz(%esp)

        ## H1 interactions 
        movaps  %xmm6,%xmm4
        mulps   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        movaps  %xmm6,%xmm7
        movaps  nb211_krsqH1(%esp),%xmm0
        addps   %xmm0,%xmm6     ## xmm6=rinv+ krsq 
        mulps   nb211_two(%esp),%xmm0
        subps   nb211_crf(%esp),%xmm6
        subps   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulps   nb211_qqH(%esp),%xmm6   ## vcoul 
        mulps   nb211_qqH(%esp),%xmm7
        mulps  %xmm7,%xmm4              ## total fsH1 in xmm4 

        addps  nb211_vctot(%esp),%xmm6

        movaps nb211_dxH1(%esp),%xmm0
        movaps nb211_dyH1(%esp),%xmm1
        movaps nb211_dzH1(%esp),%xmm2
        movaps %xmm6,nb211_vctot(%esp)
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2

        ## update H1 forces 
        movaps nb211_fixH1(%esp),%xmm3
        movaps nb211_fiyH1(%esp),%xmm4
        movaps nb211_fizH1(%esp),%xmm7
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm7
        movaps %xmm3,nb211_fixH1(%esp)
        movaps %xmm4,nb211_fiyH1(%esp)
        movaps %xmm7,nb211_fizH1(%esp)
        ## update j forces with water H1 
        addps  nb211_fjx(%esp),%xmm0
        addps  nb211_fjy(%esp),%xmm1
        addps  nb211_fjz(%esp),%xmm2
        movaps %xmm0,nb211_fjx(%esp)
        movaps %xmm1,nb211_fjy(%esp)
        movaps %xmm2,nb211_fjz(%esp)

        ## H2 interactions 
        movaps  %xmm5,%xmm4
        mulps   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        movaps  %xmm5,%xmm7
        movaps  nb211_krsqH2(%esp),%xmm0
        addps   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        mulps   nb211_two(%esp),%xmm0
        subps   nb211_crf(%esp),%xmm5
        subps   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulps   nb211_qqH(%esp),%xmm5   ## vcoul 
        mulps   nb211_qqH(%esp),%xmm7
        mulps  %xmm7,%xmm4              ## total fsH2 in xmm4 

        addps  nb211_vctot(%esp),%xmm5

        movaps nb211_dxH2(%esp),%xmm0
        movaps nb211_dyH2(%esp),%xmm1
        movaps nb211_dzH2(%esp),%xmm2
        movaps %xmm5,nb211_vctot(%esp)
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2

        ## update H2 forces 
        movaps nb211_fixH2(%esp),%xmm3
        movaps nb211_fiyH2(%esp),%xmm4
        movaps nb211_fizH2(%esp),%xmm7
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm7
        movaps %xmm3,nb211_fixH2(%esp)
        movaps %xmm4,nb211_fiyH2(%esp)
        movaps %xmm7,nb211_fizH2(%esp)

        movl nb211_faction(%ebp),%edi
        ## update j forces 
        addps nb211_fjx(%esp),%xmm0
        addps nb211_fjy(%esp),%xmm1
        addps nb211_fjz(%esp),%xmm2

        movlps (%edi,%eax,4),%xmm4
        movlps (%edi,%ecx,4),%xmm7
        movhps (%edi,%ebx,4),%xmm4
        movhps (%edi,%edx,4),%xmm7

        movaps %xmm4,%xmm3
        shufps $136,%xmm7,%xmm3 ## constant 10001000
        shufps $221,%xmm7,%xmm4 ## constant 11011101                          
        ## xmm3 has fjx, xmm4 has fjy 
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        ## unpack them back for storing 
        movaps %xmm3,%xmm7
        unpcklps %xmm4,%xmm7
        unpckhps %xmm4,%xmm3
        movlps %xmm7,(%edi,%eax,4)
        movlps %xmm3,(%edi,%ecx,4)
        movhps %xmm7,(%edi,%ebx,4)
        movhps %xmm3,(%edi,%edx,4)
        ## finally z forces 
        movss  8(%edi,%eax,4),%xmm0
        movss  8(%edi,%ebx,4),%xmm1
        movss  8(%edi,%ecx,4),%xmm3
        movss  8(%edi,%edx,4),%xmm4
        subss  %xmm2,%xmm0
        shufps $229,%xmm2,%xmm2 ## constant 11100101
        subss  %xmm2,%xmm1
        shufps $234,%xmm2,%xmm2 ## constant 11101010
        subss  %xmm2,%xmm3
        shufps $255,%xmm2,%xmm2 ## constant 11111111
        subss  %xmm2,%xmm4
        movss  %xmm0,8(%edi,%eax,4)
        movss  %xmm1,8(%edi,%ebx,4)
        movss  %xmm3,8(%edi,%ecx,4)
        movss  %xmm4,8(%edi,%edx,4)

        ## should we do one more iteration? 
        subl $4,nb211_innerk(%esp)
        jl    _nb_kernel211_ia32_sse.nb211_odd_inner
        jmp   _nb_kernel211_ia32_sse.nb211_unroll_loop
_nb_kernel211_ia32_sse.nb211_odd_inner: 
        addl $4,nb211_innerk(%esp)
        jnz   _nb_kernel211_ia32_sse.nb211_odd_loop
        jmp   _nb_kernel211_ia32_sse.nb211_updateouterdata
_nb_kernel211_ia32_sse.nb211_odd_loop: 
        movl  nb211_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb211_innerjjnr(%esp)

        xorps %xmm4,%xmm4
        movss nb211_iqO(%esp),%xmm4
        movl nb211_charge(%ebp),%esi
        movhps nb211_iqH(%esp),%xmm4
        movss (%esi,%eax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb211_qqO(%esp)    ## use oxygen qq for storage 

        xorps %xmm6,%xmm6
        movl nb211_type(%ebp),%esi
        movl (%esi,%eax,4),%ebx
        movl nb211_vdwparam(%ebp),%esi
        shll %ebx
        addl nb211_ntia(%esp),%ebx
        movlps (%esi,%ebx,4),%xmm6
        movaps %xmm6,%xmm7
        shufps $252,%xmm6,%xmm6 ## constant 11111100
        shufps $253,%xmm7,%xmm7 ## constant 11111101
        movaps %xmm6,nb211_c6(%esp)
        movaps %xmm7,nb211_c12(%esp)

        movl nb211_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## move j coords to xmm0-xmm2 
        movss (%esi,%eax,4),%xmm0
        movss 4(%esi,%eax,4),%xmm1
        movss 8(%esi,%eax,4),%xmm2
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movss nb211_ixO(%esp),%xmm3
        movss nb211_iyO(%esp),%xmm4
        movss nb211_izO(%esp),%xmm5

        movlps nb211_ixH1(%esp),%xmm6
        movlps nb211_ixH2(%esp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm3
        movlps nb211_iyH1(%esp),%xmm6
        movlps nb211_iyH2(%esp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm4
        movlps nb211_izH1(%esp),%xmm6
        movlps nb211_izH2(%esp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm5

        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        movaps %xmm3,nb211_dxO(%esp)
        movaps %xmm4,nb211_dyO(%esp)
        movaps %xmm5,nb211_dzO(%esp)

        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5

        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        ## rsq in xmm4 

        movaps %xmm4,%xmm0
        mulps nb211_krf(%esp),%xmm0
        movaps %xmm0,nb211_krsqO(%esp)

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb211_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb211_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        ## a little trick to avoid NaNs: 
        ## positions 0,2,and 3 are valid, but not 1. 
        ## If it contains NaN it doesnt help to mult by 0, 
        ## So we shuffle it and copy pos 0 to pos1! 
        shufps $224,%xmm0,%xmm0 ## constant 11100000

        movaps %xmm0,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvsq 
        movaps %xmm4,%xmm1
        mulss  %xmm4,%xmm1
        mulss  %xmm4,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  nb211_c6(%esp),%xmm1
        mulps  nb211_c12(%esp),%xmm2
        movaps %xmm2,%xmm5
        subss  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addps  nb211_Vvdwtot(%esp),%xmm5
        mulss  nb211_six(%esp),%xmm1
        mulss  nb211_twelve(%esp),%xmm2
        subss  %xmm1,%xmm2

        movaps %xmm0,%xmm1      ## xmm1=rinv 
        movaps nb211_krsqO(%esp),%xmm3
        addps  %xmm3,%xmm0      ## xmm0=rinv+ krsq 
        mulps  nb211_two(%esp),%xmm3
        subps  nb211_crf(%esp),%xmm0   ## xmm0=rinv+ krsq-crf 
        subps  %xmm3,%xmm1      ## xmm1=rinv-2*krsq 
        mulps  nb211_qqO(%esp),%xmm0    ## xmm0=vcoul 
        mulps  nb211_qqO(%esp),%xmm1    ## xmm1=coul part of fs 

        addps %xmm1,%xmm2       ## total fs 

        mulps  %xmm2,%xmm4      ## xmm4=total fscal 
        addps  nb211_vctot(%esp),%xmm0
        movaps %xmm0,nb211_vctot(%esp)

        movaps nb211_dxO(%esp),%xmm0
        movaps nb211_dyO(%esp),%xmm1
        movaps nb211_dzO(%esp),%xmm2

        movaps %xmm5,nb211_Vvdwtot(%esp)

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        movss  nb211_fixO(%esp),%xmm3
        movss  nb211_fiyO(%esp),%xmm4
        movss  nb211_fizO(%esp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,nb211_fixO(%esp)
        movss  %xmm4,nb211_fiyO(%esp)
        movss  %xmm5,nb211_fizO(%esp)   ## updated the O force now do the H's 
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        shufps $230,%xmm3,%xmm3 ## constant 11100110     ;# shift right 
        shufps $230,%xmm4,%xmm4 ## constant 11100110
        shufps $230,%xmm5,%xmm5 ## constant 11100110
        addss  nb211_fixH1(%esp),%xmm3
        addss  nb211_fiyH1(%esp),%xmm4
        addss  nb211_fizH1(%esp),%xmm5
        movss  %xmm3,nb211_fixH1(%esp)
        movss  %xmm4,nb211_fiyH1(%esp)
        movss  %xmm5,nb211_fizH1(%esp)          ## updated the H1 force 

        movl nb211_faction(%ebp),%edi
        shufps $231,%xmm3,%xmm3 ## constant 11100111     ;# shift right 
        shufps $231,%xmm4,%xmm4 ## constant 11100111
        shufps $231,%xmm5,%xmm5 ## constant 11100111
        addss  nb211_fixH2(%esp),%xmm3
        addss  nb211_fiyH2(%esp),%xmm4
        addss  nb211_fizH2(%esp),%xmm5
        movss  %xmm3,nb211_fixH2(%esp)
        movss  %xmm4,nb211_fiyH2(%esp)
        movss  %xmm5,nb211_fizH2(%esp)          ## updated the H2 force 

        ## the fj's - start by accumulating the tx/ty/tz force in xmm0, xmm1 
        xorps  %xmm5,%xmm5
        movaps %xmm0,%xmm3
        movlps (%edi,%eax,4),%xmm6
        movss  8(%edi,%eax,4),%xmm7
        unpcklps %xmm1,%xmm3
        movlhps  %xmm5,%xmm3
        unpckhps %xmm1,%xmm0
        addps    %xmm3,%xmm0
        movhlps  %xmm0,%xmm3
        addps    %xmm3,%xmm0    ## x,y sum in xmm0 

        movhlps  %xmm2,%xmm1
        addss    %xmm1,%xmm2
        shufps  $1,%xmm1,%xmm1
        addss    %xmm1,%xmm2   ## z sum in xmm2 
        subps    %xmm0,%xmm6
        subss    %xmm2,%xmm7

        movlps %xmm6,(%edi,%eax,4)
        movss  %xmm7,8(%edi,%eax,4)

        decl nb211_innerk(%esp)
        jz    _nb_kernel211_ia32_sse.nb211_updateouterdata
        jmp   _nb_kernel211_ia32_sse.nb211_odd_loop
_nb_kernel211_ia32_sse.nb211_updateouterdata: 
        movl  nb211_ii3(%esp),%ecx
        movl  nb211_faction(%ebp),%edi
        movl  nb211_fshift(%ebp),%esi
        movl  nb211_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movaps nb211_fixO(%esp),%xmm0
        movaps nb211_fiyO(%esp),%xmm1
        movaps nb211_fizO(%esp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addps  %xmm3,%xmm0
        addps  %xmm4,%xmm1
        addps  %xmm5,%xmm2 ## sum is in 1/2 in xmm0-xmm2 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        shufps $1,%xmm3,%xmm3
        shufps $1,%xmm4,%xmm4
        shufps $1,%xmm5,%xmm5
        addss  %xmm3,%xmm0
        addss  %xmm4,%xmm1
        addss  %xmm5,%xmm2      ## xmm0-xmm2 has single force in pos0 

        ## increment i force 
        movss  (%edi,%ecx,4),%xmm3
        movss  4(%edi,%ecx,4),%xmm4
        movss  8(%edi,%ecx,4),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,(%edi,%ecx,4)
        movss  %xmm4,4(%edi,%ecx,4)
        movss  %xmm5,8(%edi,%ecx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        movaps %xmm0,%xmm6
        movss %xmm2,%xmm7
        movlhps %xmm1,%xmm6
        shufps $8,%xmm6,%xmm6 ## constant 00001000      

        ## accumulate H1i forces in xmm0, xmm1, xmm2 
        movaps nb211_fixH1(%esp),%xmm0
        movaps nb211_fiyH1(%esp),%xmm1
        movaps nb211_fizH1(%esp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addps  %xmm3,%xmm0
        addps  %xmm4,%xmm1
        addps  %xmm5,%xmm2 ## sum is in 1/2 in xmm0-xmm2 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        shufps $1,%xmm3,%xmm3
        shufps $1,%xmm4,%xmm4
        shufps $1,%xmm5,%xmm5
        addss  %xmm3,%xmm0
        addss  %xmm4,%xmm1
        addss  %xmm5,%xmm2      ## xmm0-xmm2 has single force in pos0 

        ## increment i force 
        movss  12(%edi,%ecx,4),%xmm3
        movss  16(%edi,%ecx,4),%xmm4
        movss  20(%edi,%ecx,4),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,12(%edi,%ecx,4)
        movss  %xmm4,16(%edi,%ecx,4)
        movss  %xmm5,20(%edi,%ecx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        addss %xmm2,%xmm7
        movlhps %xmm1,%xmm0
        shufps $8,%xmm0,%xmm0 ## constant 00001000      
        addps   %xmm0,%xmm6

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movaps nb211_fixH2(%esp),%xmm0
        movaps nb211_fiyH2(%esp),%xmm1
        movaps nb211_fizH2(%esp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addps  %xmm3,%xmm0
        addps  %xmm4,%xmm1
        addps  %xmm5,%xmm2 ## sum is in 1/2 in xmm0-xmm2 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        shufps $1,%xmm3,%xmm3
        shufps $1,%xmm4,%xmm4
        shufps $1,%xmm5,%xmm5
        addss  %xmm3,%xmm0
        addss  %xmm4,%xmm1
        addss  %xmm5,%xmm2      ## xmm0-xmm2 has single force in pos0 

        ## increment i force 
        movss  24(%edi,%ecx,4),%xmm3
        movss  28(%edi,%ecx,4),%xmm4
        movss  32(%edi,%ecx,4),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,24(%edi,%ecx,4)
        movss  %xmm4,28(%edi,%ecx,4)
        movss  %xmm5,32(%edi,%ecx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        addss %xmm2,%xmm7
        movlhps %xmm1,%xmm0
        shufps $8,%xmm0,%xmm0 ## constant 00001000      
        addps   %xmm0,%xmm6

        ## increment fshift force  
        movlps  (%esi,%edx,4),%xmm3
        movss  8(%esi,%edx,4),%xmm4
        addps  %xmm6,%xmm3
        addss  %xmm7,%xmm4
        movlps  %xmm3,(%esi,%edx,4)
        movss  %xmm4,8(%esi,%edx,4)

        ## get n from stack
        movl nb211_n(%esp),%esi
        ## get group index for i particle 
        movl  nb211_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb211_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb211_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb211_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb211_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb211_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel211_ia32_sse.nb211_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb211_n(%esp)
        jmp _nb_kernel211_ia32_sse.nb211_outer
_nb_kernel211_ia32_sse.nb211_outerend: 
        ## check if more outer neighborlists remain
        movl  nb211_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel211_ia32_sse.nb211_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel211_ia32_sse.nb211_threadloop
_nb_kernel211_ia32_sse.nb211_end: 
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



.globl nb_kernel211nf_ia32_sse
.globl _nb_kernel211nf_ia32_sse
nb_kernel211nf_ia32_sse:        
_nb_kernel211nf_ia32_sse:       
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


        movl nb211nf_argkrf(%ebp),%esi
        movl nb211nf_argcrf(%ebp),%edi
        movss (%esi),%xmm5
        movss (%edi),%xmm6
        shufps $0,%xmm5,%xmm5
        shufps $0,%xmm6,%xmm6
        movaps %xmm5,nb211nf_krf(%esp)
        movaps %xmm6,nb211nf_crf(%esp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb211nf_half(%esp)
        movss nb211nf_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb211nf_half(%esp)
        movaps %xmm3,nb211nf_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb211nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb211nf_charge(%ebp),%edx
        movss (%edx,%ebx,4),%xmm3
        movss 4(%edx,%ebx,4),%xmm4
        movl nb211nf_p_facel(%ebp),%esi
        movss (%esi),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb211nf_iqO(%esp)
        movaps %xmm4,nb211nf_iqH(%esp)

        movl  nb211nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb211nf_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb211nf_ntia(%esp)

_nb_kernel211nf_ia32_sse.nb211nf_threadloop: 
        movl  nb211nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel211nf_ia32_sse.nb211nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel211nf_ia32_sse.nb211nf_spinlock

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
        jg  _nb_kernel211nf_ia32_sse.nb211nf_outerstart
        jmp _nb_kernel211nf_ia32_sse.nb211nf_end

_nb_kernel211nf_ia32_sse.nb211nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb211nf_nouter(%esp),%ebx
        movl %ebx,nb211nf_nouter(%esp)

_nb_kernel211nf_ia32_sse.nb211nf_outer: 
        movl  nb211nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb211nf_is3(%esp)            ## store is3 

        movl  nb211nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb211nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb211nf_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb211nf_ii3(%esp)

        addss (%eax,%ebx,4),%xmm3
        addss 4(%eax,%ebx,4),%xmm4
        addss 8(%eax,%ebx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb211nf_ixO(%esp)
        movaps %xmm4,nb211nf_iyO(%esp)
        movaps %xmm5,nb211nf_izO(%esp)

        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 12(%eax,%ebx,4),%xmm0
        addss 16(%eax,%ebx,4),%xmm1
        addss 20(%eax,%ebx,4),%xmm2
        addss 24(%eax,%ebx,4),%xmm3
        addss 28(%eax,%ebx,4),%xmm4
        addss 32(%eax,%ebx,4),%xmm5

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm0,nb211nf_ixH1(%esp)
        movaps %xmm1,nb211nf_iyH1(%esp)
        movaps %xmm2,nb211nf_izH1(%esp)
        movaps %xmm3,nb211nf_ixH2(%esp)
        movaps %xmm4,nb211nf_iyH2(%esp)
        movaps %xmm5,nb211nf_izH2(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb211nf_vctot(%esp)
        movaps %xmm4,nb211nf_Vvdwtot(%esp)

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
        subl  $4,%edx
        addl  nb211nf_ninner(%esp),%ecx
        movl  %ecx,nb211nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb211nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel211nf_ia32_sse.nb211nf_unroll_loop
        jmp   _nb_kernel211nf_ia32_sse.nb211nf_odd_inner
_nb_kernel211nf_ia32_sse.nb211nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb211nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 

        addl $16,nb211nf_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb211nf_charge(%ebp),%esi     ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4           ## and in xmm4 
        mulps  nb211nf_iqO(%esp),%xmm3
        mulps  nb211nf_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps  %xmm3,nb211nf_qqO(%esp)
        movaps  %xmm4,nb211nf_qqH(%esp)

        movl nb211nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb211nf_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb211nf_ntia(%esp),%edi
        addl %edi,%eax
        addl %edi,%ebx
        addl %edi,%ecx
        addl %edi,%edx

        movlps (%esi,%eax,4),%xmm6
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm6
        movhps (%esi,%edx,4),%xmm7

        movaps %xmm6,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm6 ## constant 11011101

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movd  %mm2,%ecx
        movd  %mm3,%edx

        movaps %xmm4,nb211nf_c6(%esp)
        movaps %xmm6,nb211nf_c12(%esp)

        movl nb211nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx     ## replace jnr with j3 
        leal  (%edx,%edx,2),%edx

        ## move four coordinates to xmm0-xmm2   
        movlps (%esi,%eax,4),%xmm4
        movlps (%esi,%ecx,4),%xmm5
        movss 8(%esi,%eax,4),%xmm2
        movss 8(%esi,%ecx,4),%xmm6

        movhps (%esi,%ebx,4),%xmm4
        movhps (%esi,%edx,4),%xmm5

        movss 8(%esi,%ebx,4),%xmm0
        movss 8(%esi,%edx,4),%xmm1

        shufps $0,%xmm0,%xmm2
        shufps $0,%xmm1,%xmm6

        movaps %xmm4,%xmm0
        movaps %xmm4,%xmm1

        shufps $136,%xmm6,%xmm2 ## constant 10001000

        shufps $136,%xmm5,%xmm0 ## constant 10001000
        shufps $221,%xmm5,%xmm1 ## constant 11011101            

        ## move ixO-izO to xmm4-xmm6 
        movaps nb211nf_ixO(%esp),%xmm4
        movaps nb211nf_iyO(%esp),%xmm5
        movaps nb211nf_izO(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        movaps %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movaps nb211nf_ixH1(%esp),%xmm4
        movaps nb211nf_iyH1(%esp),%xmm5
        movaps nb211nf_izH1(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm6
        addps %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movaps nb211nf_ixH2(%esp),%xmm3
        movaps nb211nf_iyH2(%esp),%xmm4
        movaps nb211nf_izH2(%esp),%xmm5

        ## calc dr 
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        ## square it 
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm4,%xmm5
        addps %xmm3,%xmm5
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        movaps %xmm5,%xmm0
        movaps %xmm6,%xmm1
        movaps %xmm7,%xmm2

        mulps  nb211nf_krf(%esp),%xmm0
        mulps  nb211nf_krf(%esp),%xmm1
        mulps  nb211nf_krf(%esp),%xmm2

        movaps %xmm0,nb211nf_krsqH2(%esp)
        movaps %xmm1,nb211nf_krsqH1(%esp)
        movaps %xmm2,nb211nf_krsqO(%esp)

        ## start with rsqO - seed in xmm2       
        rsqrtps %xmm7,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb211nf_three(%esp),%xmm4
        mulps   %xmm7,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb211nf_half(%esp),%xmm4
        movaps  %xmm4,%xmm7     ## rinvO in xmm7 
        ## rsqH1 - seed in xmm2 
        rsqrtps %xmm6,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb211nf_three(%esp),%xmm4
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb211nf_half(%esp),%xmm4
        movaps  %xmm4,%xmm6     ## rinvH1 in xmm6 
        ## rsqH2 - seed in xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb211nf_three(%esp),%xmm4
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb211nf_half(%esp),%xmm4
        movaps  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do O interactions 
        movaps  %xmm7,%xmm4
        mulps   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq 
        movaps %xmm4,%xmm1
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  nb211nf_c6(%esp),%xmm1
        mulps  nb211nf_c12(%esp),%xmm2
        movaps %xmm2,%xmm3
        subps  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addps  nb211nf_Vvdwtot(%esp),%xmm3

        movaps %xmm7,%xmm0
        movaps nb211nf_krsqO(%esp),%xmm1
        addps  %xmm1,%xmm0
        subps  nb211nf_crf(%esp),%xmm0   ## xmm0=rinv+ krsq-crf 
        subps  %xmm1,%xmm7
        mulps  nb211nf_qqO(%esp),%xmm0
        addps  nb211nf_vctot(%esp),%xmm0
        movaps %xmm3,nb211nf_Vvdwtot(%esp)
        movaps %xmm0,nb211nf_vctot(%esp)

        ## H1 interactions 
        movaps  nb211nf_krsqH1(%esp),%xmm0
        addps   %xmm0,%xmm6     ## xmm6=rinv+ krsq 
        subps   nb211nf_crf(%esp),%xmm6
        mulps   nb211nf_qqH(%esp),%xmm6   ## vcoul 
        addps  nb211nf_vctot(%esp),%xmm6
        movaps %xmm6,nb211nf_vctot(%esp)

        ## H2 interactions 
        movaps  %xmm5,%xmm7 ## rinv 
        movaps  nb211nf_krsqH2(%esp),%xmm0
        addps   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        subps   nb211nf_crf(%esp),%xmm5
        mulps   nb211nf_qqH(%esp),%xmm5   ## vcoul 
        addps   nb211nf_vctot(%esp),%xmm5
        movaps %xmm5,nb211nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $4,nb211nf_innerk(%esp)
        jl    _nb_kernel211nf_ia32_sse.nb211nf_odd_inner
        jmp   _nb_kernel211nf_ia32_sse.nb211nf_unroll_loop
_nb_kernel211nf_ia32_sse.nb211nf_odd_inner: 
        addl $4,nb211nf_innerk(%esp)
        jnz   _nb_kernel211nf_ia32_sse.nb211nf_odd_loop
        jmp   _nb_kernel211nf_ia32_sse.nb211nf_updateouterdata
_nb_kernel211nf_ia32_sse.nb211nf_odd_loop: 
        movl  nb211nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb211nf_innerjjnr(%esp)

        xorps %xmm4,%xmm4
        movss nb211nf_iqO(%esp),%xmm4
        movl nb211nf_charge(%ebp),%esi
        movhps nb211nf_iqH(%esp),%xmm4
        movss (%esi,%eax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb211nf_qqO(%esp)          ## use oxygen qq for storage 

        xorps %xmm6,%xmm6
        movl nb211nf_type(%ebp),%esi
        movl (%esi,%eax,4),%ebx
        movl nb211nf_vdwparam(%ebp),%esi
        shll %ebx
        addl nb211nf_ntia(%esp),%ebx
        movlps (%esi,%ebx,4),%xmm6
        movaps %xmm6,%xmm7
        shufps $252,%xmm6,%xmm6 ## constant 11111100
        shufps $253,%xmm7,%xmm7 ## constant 11111101
        movaps %xmm6,nb211nf_c6(%esp)
        movaps %xmm7,nb211nf_c12(%esp)

        movl nb211nf_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## move j coords to xmm0-xmm2 
        movss (%esi,%eax,4),%xmm0
        movss 4(%esi,%eax,4),%xmm1
        movss 8(%esi,%eax,4),%xmm2
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movss nb211nf_ixO(%esp),%xmm3
        movss nb211nf_iyO(%esp),%xmm4
        movss nb211nf_izO(%esp),%xmm5

        movlps nb211nf_ixH1(%esp),%xmm6
        movlps nb211nf_ixH2(%esp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm3
        movlps nb211nf_iyH1(%esp),%xmm6
        movlps nb211nf_iyH2(%esp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm4
        movlps nb211nf_izH1(%esp),%xmm6
        movlps nb211nf_izH2(%esp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm5

        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5

        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        ## rsq in xmm4 

        movaps %xmm4,%xmm0
        mulps nb211nf_krf(%esp),%xmm0
        movaps %xmm0,nb211nf_krsqO(%esp)

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb211nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb211nf_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 

        ## a little trick to avoid NaNs: 
        ## positions 0,2,and 3 are valid, but not 1. 
        ## If it contains NaN it doesnt help to mult by 0, 
        ## So we shuffle it and copy pos 0 to pos1! 
        shufps $224,%xmm0,%xmm0 ## constant 11100000     

        movaps %xmm0,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvsq 
        movaps %xmm4,%xmm1
        mulss  %xmm4,%xmm1
        mulss  %xmm4,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  nb211nf_c6(%esp),%xmm1
        mulps  nb211nf_c12(%esp),%xmm2
        movaps %xmm2,%xmm5
        subss  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addps  nb211nf_Vvdwtot(%esp),%xmm5
        movaps %xmm0,%xmm1      ## xmm1=rinv 
        movaps nb211nf_krsqO(%esp),%xmm3
        addps  %xmm3,%xmm0      ## xmm0=rinv+ krsq 
        subps  nb211nf_crf(%esp),%xmm0   ## xmm0=rinv+ krsq-crf 
        mulps  nb211nf_qqO(%esp),%xmm0          ## xmm0=vcoul 
        addps  nb211nf_vctot(%esp),%xmm0
        movaps %xmm0,nb211nf_vctot(%esp)
        movaps %xmm5,nb211nf_Vvdwtot(%esp)

        decl nb211nf_innerk(%esp)
        jz    _nb_kernel211nf_ia32_sse.nb211nf_updateouterdata
        jmp   _nb_kernel211nf_ia32_sse.nb211nf_odd_loop
_nb_kernel211nf_ia32_sse.nb211nf_updateouterdata: 
        ## get n from stack
        movl nb211nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb211nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb211nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb211nf_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb211nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb211nf_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb211nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel211nf_ia32_sse.nb211nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb211nf_n(%esp)
        jmp _nb_kernel211nf_ia32_sse.nb211nf_outer
_nb_kernel211nf_ia32_sse.nb211nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb211nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel211nf_ia32_sse.nb211nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel211nf_ia32_sse.nb211nf_threadloop
_nb_kernel211nf_ia32_sse.nb211nf_end: 
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


