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



.globl nb_kernel111_ia32_sse
.globl _nb_kernel111_ia32_sse
nb_kernel111_ia32_sse:  
_nb_kernel111_ia32_sse: 
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
.set nb111_p_krf, 56
.set nb111_p_crf, 60
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
        ## bottom of stack is cache-aligned for sse use 
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


        ## assume we have at least one i particle - start directly 
        movl  nb111_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb111_charge(%ebp),%edx
        movss (%edx,%ebx,4),%xmm3
        movss 4(%edx,%ebx,4),%xmm4
        movl nb111_p_facel(%ebp),%esi
        movss (%esi),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb111_iqO(%esp)
        movaps %xmm4,nb111_iqH(%esp)

        movl  nb111_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb111_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb111_ntia(%esp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb111_half(%esp)
        movss nb111_half(%esp),%xmm1
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
        movaps %xmm1,nb111_half(%esp)
        movaps %xmm3,nb111_three(%esp)
        movaps %xmm4,nb111_six(%esp)
        movaps %xmm5,nb111_twelve(%esp)

_nb_kernel111_ia32_sse.nb111_threadloop: 
        movl  nb111_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel111_ia32_sse.nb111_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel111_ia32_sse.nb111_spinlock

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
        jg  _nb_kernel111_ia32_sse.nb111_outerstart
        jmp _nb_kernel111_ia32_sse.nb111_end

_nb_kernel111_ia32_sse.nb111_outerstart: 
        ## ebx contains number of outer iterations
        addl nb111_nouter(%esp),%ebx
        movl %ebx,nb111_nouter(%esp)

_nb_kernel111_ia32_sse.nb111_outer: 
        movl  nb111_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb111_is3(%esp)      ## store is3 

        movl  nb111_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb111_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb111_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb111_ii3(%esp)

        addss (%eax,%ebx,4),%xmm3
        addss 4(%eax,%ebx,4),%xmm4
        addss 8(%eax,%ebx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb111_ixO(%esp)
        movaps %xmm4,nb111_iyO(%esp)
        movaps %xmm5,nb111_izO(%esp)

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
        movaps %xmm0,nb111_ixH1(%esp)
        movaps %xmm1,nb111_iyH1(%esp)
        movaps %xmm2,nb111_izH1(%esp)
        movaps %xmm3,nb111_ixH2(%esp)
        movaps %xmm4,nb111_iyH2(%esp)
        movaps %xmm5,nb111_izH2(%esp)

        ## clear potentials and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb111_vctot(%esp)
        movaps %xmm4,nb111_Vvdwtot(%esp)
        movaps %xmm4,nb111_fixO(%esp)
        movaps %xmm4,nb111_fiyO(%esp)
        movaps %xmm4,nb111_fizO(%esp)
        movaps %xmm4,nb111_fixH1(%esp)
        movaps %xmm4,nb111_fiyH1(%esp)
        movaps %xmm4,nb111_fizH1(%esp)
        movaps %xmm4,nb111_fixH2(%esp)
        movaps %xmm4,nb111_fiyH2(%esp)
        movaps %xmm4,nb111_fizH2(%esp)

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
        subl  $4,%edx
        addl  nb111_ninner(%esp),%ecx
        movl  %ecx,nb111_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb111_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel111_ia32_sse.nb111_unroll_loop
        jmp   _nb_kernel111_ia32_sse.nb111_odd_inner
_nb_kernel111_ia32_sse.nb111_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb111_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 

        addl $16,nb111_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb111_charge(%ebp),%esi     ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4           ## and in xmm4 
        mulps  nb111_iqO(%esp),%xmm3
        mulps  nb111_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps  %xmm3,nb111_qqO(%esp)
        movaps  %xmm4,nb111_qqH(%esp)

        movl nb111_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb111_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb111_ntia(%esp),%edi
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

        movaps %xmm4,nb111_c6(%esp)
        movaps %xmm6,nb111_c12(%esp)

        movl nb111_pos(%ebp),%esi        ## base of pos[] 

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
        movaps nb111_ixO(%esp),%xmm4
        movaps nb111_iyO(%esp),%xmm5
        movaps nb111_izO(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb111_dxO(%esp)
        movaps %xmm5,nb111_dyO(%esp)
        movaps %xmm6,nb111_dzO(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        movaps %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movaps nb111_ixH1(%esp),%xmm4
        movaps nb111_iyH1(%esp),%xmm5
        movaps nb111_izH1(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb111_dxH1(%esp)
        movaps %xmm5,nb111_dyH1(%esp)
        movaps %xmm6,nb111_dzH1(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm6
        addps %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movaps nb111_ixH2(%esp),%xmm3
        movaps nb111_iyH2(%esp),%xmm4
        movaps nb111_izH2(%esp),%xmm5

        ## calc dr 
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        ## store dr 
        movaps %xmm3,nb111_dxH2(%esp)
        movaps %xmm4,nb111_dyH2(%esp)
        movaps %xmm5,nb111_dzH2(%esp)
        ## square it 
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm4,%xmm5
        addps %xmm3,%xmm5
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## start with rsqO - seed in xmm2       
        rsqrtps %xmm7,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb111_three(%esp),%xmm4
        mulps   %xmm7,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb111_half(%esp),%xmm4
        movaps  %xmm4,%xmm7     ## rinvO in xmm7 
        ## rsqH1 - seed in xmm2 
        rsqrtps %xmm6,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb111_three(%esp),%xmm4
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb111_half(%esp),%xmm4
        movaps  %xmm4,%xmm6     ## rinvH1 in xmm6 
        ## rsqH2 - seed in xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb111_three(%esp),%xmm4
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb111_half(%esp),%xmm4
        movaps  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do O interactions 
        movaps  %xmm7,%xmm4
        mulps   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq 
        movaps %xmm4,%xmm1
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  nb111_qqO(%esp),%xmm7    ## xmm7=vcoul 

        mulps  nb111_c6(%esp),%xmm1
        mulps  nb111_c12(%esp),%xmm2
        movaps %xmm2,%xmm3
        subps  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addps  nb111_Vvdwtot(%esp),%xmm3
        mulps  nb111_six(%esp),%xmm1
        mulps  nb111_twelve(%esp),%xmm2
        subps  %xmm1,%xmm2
        addps  %xmm7,%xmm2
        mulps  %xmm2,%xmm4      ## total fsO in xmm4 

        addps  nb111_vctot(%esp),%xmm7

        movaps %xmm3,nb111_Vvdwtot(%esp)
        movaps %xmm7,nb111_vctot(%esp)

        movaps nb111_dxO(%esp),%xmm0
        movaps nb111_dyO(%esp),%xmm1
        movaps nb111_dzO(%esp),%xmm2
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2

        ## update O forces 
        movaps nb111_fixO(%esp),%xmm3
        movaps nb111_fiyO(%esp),%xmm4
        movaps nb111_fizO(%esp),%xmm7
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm7
        movaps %xmm3,nb111_fixO(%esp)
        movaps %xmm4,nb111_fiyO(%esp)
        movaps %xmm7,nb111_fizO(%esp)
        ## update j forces with water O 
        movaps %xmm0,nb111_fjx(%esp)
        movaps %xmm1,nb111_fjy(%esp)
        movaps %xmm2,nb111_fjz(%esp)

        ## H1 interactions 
        movaps  %xmm6,%xmm4
        mulps   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        mulps  nb111_qqH(%esp),%xmm6    ## xmm6=vcoul 
        mulps  %xmm6,%xmm4              ## total fsH1 in xmm4 

        addps  nb111_vctot(%esp),%xmm6

        movaps nb111_dxH1(%esp),%xmm0
        movaps nb111_dyH1(%esp),%xmm1
        movaps nb111_dzH1(%esp),%xmm2
        movaps %xmm6,nb111_vctot(%esp)
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2

        ## update H1 forces 
        movaps nb111_fixH1(%esp),%xmm3
        movaps nb111_fiyH1(%esp),%xmm4
        movaps nb111_fizH1(%esp),%xmm7
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm7
        movaps %xmm3,nb111_fixH1(%esp)
        movaps %xmm4,nb111_fiyH1(%esp)
        movaps %xmm7,nb111_fizH1(%esp)
        ## update j forces with water H1 
        addps  nb111_fjx(%esp),%xmm0
        addps  nb111_fjy(%esp),%xmm1
        addps  nb111_fjz(%esp),%xmm2
        movaps %xmm0,nb111_fjx(%esp)
        movaps %xmm1,nb111_fjy(%esp)
        movaps %xmm2,nb111_fjz(%esp)

        ## H2 interactions 
        movaps  %xmm5,%xmm4
        mulps   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        mulps  nb111_qqH(%esp),%xmm5    ## xmm5=vcoul 
        mulps  %xmm5,%xmm4              ## total fsH1 in xmm4 

        addps  nb111_vctot(%esp),%xmm5

        movaps nb111_dxH2(%esp),%xmm0
        movaps nb111_dyH2(%esp),%xmm1
        movaps nb111_dzH2(%esp),%xmm2
        movaps %xmm5,nb111_vctot(%esp)
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2

        ## update H2 forces 
        movaps nb111_fixH2(%esp),%xmm3
        movaps nb111_fiyH2(%esp),%xmm4
        movaps nb111_fizH2(%esp),%xmm7
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm7
        movaps %xmm3,nb111_fixH2(%esp)
        movaps %xmm4,nb111_fiyH2(%esp)
        movaps %xmm7,nb111_fizH2(%esp)

        movl nb111_faction(%ebp),%edi
        ## update j forces 
        addps nb111_fjx(%esp),%xmm0
        addps nb111_fjy(%esp),%xmm1
        addps nb111_fjz(%esp),%xmm2

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
        subl $4,nb111_innerk(%esp)
        jl    _nb_kernel111_ia32_sse.nb111_odd_inner
        jmp   _nb_kernel111_ia32_sse.nb111_unroll_loop
_nb_kernel111_ia32_sse.nb111_odd_inner: 
        addl $4,nb111_innerk(%esp)
        jnz   _nb_kernel111_ia32_sse.nb111_odd_loop
        jmp   _nb_kernel111_ia32_sse.nb111_updateouterdata
_nb_kernel111_ia32_sse.nb111_odd_loop: 
        movl  nb111_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb111_innerjjnr(%esp)

        xorps %xmm4,%xmm4
        movss nb111_iqO(%esp),%xmm4
        movl nb111_charge(%ebp),%esi
        movhps nb111_iqH(%esp),%xmm4
        movss (%esi,%eax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb111_qqO(%esp)    ## use oxygen qq for storage 

        xorps %xmm6,%xmm6
        movl nb111_type(%ebp),%esi
        movl (%esi,%eax,4),%ebx
        movl nb111_vdwparam(%ebp),%esi
        shll %ebx
        addl nb111_ntia(%esp),%ebx
        movlps (%esi,%ebx,4),%xmm6
        movaps %xmm6,%xmm7
        shufps $252,%xmm6,%xmm6 ## constant 11111100
        shufps $253,%xmm7,%xmm7 ## constant 11111101
        movaps %xmm6,nb111_c6(%esp)
        movaps %xmm7,nb111_c12(%esp)

        movl nb111_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## move j coords to xmm0-xmm2 
        movss (%esi,%eax,4),%xmm0
        movss 4(%esi,%eax,4),%xmm1
        movss 8(%esi,%eax,4),%xmm2
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movss nb111_ixO(%esp),%xmm3
        movss nb111_iyO(%esp),%xmm4
        movss nb111_izO(%esp),%xmm5

        movlps nb111_ixH1(%esp),%xmm6
        movlps nb111_ixH2(%esp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm3
        movlps nb111_iyH1(%esp),%xmm6
        movlps nb111_iyH2(%esp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm4
        movlps nb111_izH1(%esp),%xmm6
        movlps nb111_izH2(%esp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm5

        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        movaps %xmm3,nb111_dxO(%esp)
        movaps %xmm4,nb111_dyO(%esp)
        movaps %xmm5,nb111_dzO(%esp)

        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5

        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        ## rsq in xmm4 

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb111_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb111_half(%esp),%xmm0
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
        movaps nb111_qqO(%esp),%xmm3
        mulss  %xmm4,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  %xmm0,%xmm3      ## xmm3=vcoul 
        mulps  nb111_c6(%esp),%xmm1
        mulps  nb111_c12(%esp),%xmm2
        movaps %xmm2,%xmm5
        subss  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addps  nb111_Vvdwtot(%esp),%xmm5
        mulss  nb111_six(%esp),%xmm1
        mulss  nb111_twelve(%esp),%xmm2
        subss  %xmm1,%xmm2
        addps  %xmm3,%xmm2
        mulps  %xmm2,%xmm4      ## xmm4=total fscal 
        addps  nb111_vctot(%esp),%xmm3

        movaps nb111_dxO(%esp),%xmm0
        movaps nb111_dyO(%esp),%xmm1
        movaps nb111_dzO(%esp),%xmm2

        movaps %xmm3,nb111_vctot(%esp)
        movaps %xmm5,nb111_Vvdwtot(%esp)

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        movss  nb111_fixO(%esp),%xmm3
        movss  nb111_fiyO(%esp),%xmm4
        movss  nb111_fizO(%esp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,nb111_fixO(%esp)
        movss  %xmm4,nb111_fiyO(%esp)
        movss  %xmm5,nb111_fizO(%esp)   ## updated the O force now do the H's 
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        shufps $230,%xmm3,%xmm3 ## constant 11100110     ;# shift right 
        shufps $230,%xmm4,%xmm4 ## constant 11100110
        shufps $230,%xmm5,%xmm5 ## constant 11100110
        addss  nb111_fixH1(%esp),%xmm3
        addss  nb111_fiyH1(%esp),%xmm4
        addss  nb111_fizH1(%esp),%xmm5
        movss  %xmm3,nb111_fixH1(%esp)
        movss  %xmm4,nb111_fiyH1(%esp)
        movss  %xmm5,nb111_fizH1(%esp)          ## updated the H1 force 

        movl nb111_faction(%ebp),%edi
        shufps $231,%xmm3,%xmm3 ## constant 11100111     ;# shift right 
        shufps $231,%xmm4,%xmm4 ## constant 11100111
        shufps $231,%xmm5,%xmm5 ## constant 11100111
        addss  nb111_fixH2(%esp),%xmm3
        addss  nb111_fiyH2(%esp),%xmm4
        addss  nb111_fizH2(%esp),%xmm5
        movss  %xmm3,nb111_fixH2(%esp)
        movss  %xmm4,nb111_fiyH2(%esp)
        movss  %xmm5,nb111_fizH2(%esp)          ## updated the H2 force 

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

        decl nb111_innerk(%esp)
        jz    _nb_kernel111_ia32_sse.nb111_updateouterdata
        jmp   _nb_kernel111_ia32_sse.nb111_odd_loop
_nb_kernel111_ia32_sse.nb111_updateouterdata: 
        movl  nb111_ii3(%esp),%ecx
        movl  nb111_faction(%ebp),%edi
        movl  nb111_fshift(%ebp),%esi
        movl  nb111_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movaps nb111_fixO(%esp),%xmm0
        movaps nb111_fiyO(%esp),%xmm1
        movaps nb111_fizO(%esp),%xmm2

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
        movaps nb111_fixH1(%esp),%xmm0
        movaps nb111_fiyH1(%esp),%xmm1
        movaps nb111_fizH1(%esp),%xmm2

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
        movaps nb111_fixH2(%esp),%xmm0
        movaps nb111_fiyH2(%esp),%xmm1
        movaps nb111_fizH2(%esp),%xmm2

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
        movl nb111_n(%esp),%esi
        ## get group index for i particle 
        movl  nb111_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb111_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb111_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb111_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb111_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb111_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel111_ia32_sse.nb111_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb111_n(%esp)
        jmp _nb_kernel111_ia32_sse.nb111_outer
_nb_kernel111_ia32_sse.nb111_outerend: 
        ## check if more outer neighborlists remain
        movl  nb111_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel111_ia32_sse.nb111_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel111_ia32_sse.nb111_threadloop
_nb_kernel111_ia32_sse.nb111_end: 
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



.globl nb_kernel111nf_ia32_sse
.globl _nb_kernel111nf_ia32_sse
nb_kernel111nf_ia32_sse:        
_nb_kernel111nf_ia32_sse:       
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
.set nb111nf_p_krf, 56
.set nb111nf_p_crf, 60
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


        ## assume we have at least one i particle - start directly 
        movl  nb111nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb111nf_charge(%ebp),%edx
        movss (%edx,%ebx,4),%xmm3
        movss 4(%edx,%ebx,4),%xmm4
        movl nb111nf_p_facel(%ebp),%esi
        movss (%esi),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb111nf_iqO(%esp)
        movaps %xmm4,nb111nf_iqH(%esp)

        movl  nb111nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb111nf_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb111nf_ntia(%esp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb111nf_half(%esp)
        movss nb111nf_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb111nf_half(%esp)
        movaps %xmm3,nb111nf_three(%esp)

_nb_kernel111nf_ia32_sse.nb111nf_threadloop: 
        movl  nb111nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel111nf_ia32_sse.nb111nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel111nf_ia32_sse.nb111nf_spinlock

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
        jg  _nb_kernel111nf_ia32_sse.nb111nf_outerstart
        jmp _nb_kernel111nf_ia32_sse.nb111nf_end

_nb_kernel111nf_ia32_sse.nb111nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb111nf_nouter(%esp),%ebx
        movl %ebx,nb111nf_nouter(%esp)

_nb_kernel111nf_ia32_sse.nb111nf_outer: 
        movl  nb111nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb111nf_is3(%esp)            ## store is3 

        movl  nb111nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb111nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb111nf_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb111nf_ii3(%esp)

        addss (%eax,%ebx,4),%xmm3
        addss 4(%eax,%ebx,4),%xmm4
        addss 8(%eax,%ebx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb111nf_ixO(%esp)
        movaps %xmm4,nb111nf_iyO(%esp)
        movaps %xmm5,nb111nf_izO(%esp)

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
        movaps %xmm0,nb111nf_ixH1(%esp)
        movaps %xmm1,nb111nf_iyH1(%esp)
        movaps %xmm2,nb111nf_izH1(%esp)
        movaps %xmm3,nb111nf_ixH2(%esp)
        movaps %xmm4,nb111nf_iyH2(%esp)
        movaps %xmm5,nb111nf_izH2(%esp)

        ## clear vctot and Vvdwtot
        xorps %xmm4,%xmm4
        movaps %xmm4,nb111nf_vctot(%esp)
        movaps %xmm4,nb111nf_Vvdwtot(%esp)

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
        subl  $4,%edx
        addl  nb111nf_ninner(%esp),%ecx
        movl  %ecx,nb111nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb111nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel111nf_ia32_sse.nb111nf_unroll_loop
        jmp   _nb_kernel111nf_ia32_sse.nb111nf_odd_inner
_nb_kernel111nf_ia32_sse.nb111nf_unroll_loop: 

        ## quad-unroll innerloop here 
        movl  nb111nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 

        addl $16,nb111nf_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb111nf_charge(%ebp),%esi     ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4           ## and in xmm4 
        mulps  nb111nf_iqO(%esp),%xmm3
        mulps  nb111nf_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps  %xmm3,nb111nf_qqO(%esp)
        movaps  %xmm4,nb111nf_qqH(%esp)

        movl nb111nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb111nf_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb111nf_ntia(%esp),%edi
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

        movaps %xmm4,nb111nf_c6(%esp)
        movaps %xmm6,nb111nf_c12(%esp)

        movl nb111nf_pos(%ebp),%esi        ## base of pos[] 

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
        movaps nb111nf_ixO(%esp),%xmm4
        movaps nb111nf_iyO(%esp),%xmm5
        movaps nb111nf_izO(%esp),%xmm6

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
        movaps nb111nf_ixH1(%esp),%xmm4
        movaps nb111nf_iyH1(%esp),%xmm5
        movaps nb111nf_izH1(%esp),%xmm6

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
        movaps nb111nf_ixH2(%esp),%xmm3
        movaps nb111nf_iyH2(%esp),%xmm4
        movaps nb111nf_izH2(%esp),%xmm5

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

        ## start with rsqO - seed in xmm2       
        rsqrtps %xmm7,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb111nf_three(%esp),%xmm4
        mulps   %xmm7,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb111nf_half(%esp),%xmm4
        movaps  %xmm4,%xmm7     ## rinvO in xmm7 
        ## rsqH1 - seed in xmm2 
        rsqrtps %xmm6,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb111nf_three(%esp),%xmm4
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb111nf_half(%esp),%xmm4
        movaps  %xmm4,%xmm6     ## rinvH1 in xmm6 
        ## rsqH2 - seed in xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb111nf_three(%esp),%xmm4
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb111nf_half(%esp),%xmm4
        movaps  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do O interactions 
        movaps  %xmm7,%xmm4
        mulps   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq 
        movaps %xmm4,%xmm1
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  nb111nf_qqO(%esp),%xmm7          ## xmm7=vcoul 

        mulps  nb111nf_c6(%esp),%xmm1
        mulps  nb111nf_c12(%esp),%xmm2
        movaps %xmm2,%xmm3
        subps  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addps  nb111nf_Vvdwtot(%esp),%xmm3
        addps  nb111nf_vctot(%esp),%xmm7
        movaps %xmm3,nb111nf_Vvdwtot(%esp)
        movaps %xmm7,nb111nf_vctot(%esp)

        ## H1 & H2 interactions 
        addps  %xmm5,%xmm6          ## add H2 rinv 
        mulps  nb111nf_qqH(%esp),%xmm6          ## xmm6=vcoul 
        addps  nb111nf_vctot(%esp),%xmm6
        movaps %xmm6,nb111nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $4,nb111nf_innerk(%esp)
        jl    _nb_kernel111nf_ia32_sse.nb111nf_odd_inner
        jmp   _nb_kernel111nf_ia32_sse.nb111nf_unroll_loop
_nb_kernel111nf_ia32_sse.nb111nf_odd_inner: 
        addl $4,nb111nf_innerk(%esp)
        jnz   _nb_kernel111nf_ia32_sse.nb111nf_odd_loop
        jmp   _nb_kernel111nf_ia32_sse.nb111nf_updateouterdata
_nb_kernel111nf_ia32_sse.nb111nf_odd_loop: 
        movl  nb111nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb111nf_innerjjnr(%esp)

        xorps %xmm4,%xmm4
        movss nb111nf_iqO(%esp),%xmm4
        movl nb111nf_charge(%ebp),%esi
        movhps nb111nf_iqH(%esp),%xmm4
        movss (%esi,%eax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb111nf_qqO(%esp)          ## use oxygen qq for storage 

        xorps %xmm6,%xmm6
        movl nb111nf_type(%ebp),%esi
        movl (%esi,%eax,4),%ebx
        movl nb111nf_vdwparam(%ebp),%esi
        shll %ebx
        addl nb111nf_ntia(%esp),%ebx
        movlps (%esi,%ebx,4),%xmm6
        movaps %xmm6,%xmm7
        shufps $252,%xmm6,%xmm6 ## constant 11111100
        shufps $253,%xmm7,%xmm7 ## constant 11111101
        movaps %xmm6,nb111nf_c6(%esp)
        movaps %xmm7,nb111nf_c12(%esp)

        movl nb111nf_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## move j coords to xmm0-xmm2 
        movss (%esi,%eax,4),%xmm0
        movss 4(%esi,%eax,4),%xmm1
        movss 8(%esi,%eax,4),%xmm2
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movss nb111nf_ixO(%esp),%xmm3
        movss nb111nf_iyO(%esp),%xmm4
        movss nb111nf_izO(%esp),%xmm5

        movlps nb111nf_ixH1(%esp),%xmm6
        movlps nb111nf_ixH2(%esp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm3
        movlps nb111nf_iyH1(%esp),%xmm6
        movlps nb111nf_iyH2(%esp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm4
        movlps nb111nf_izH1(%esp),%xmm6
        movlps nb111nf_izH2(%esp),%xmm7
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

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb111nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb111nf_half(%esp),%xmm0
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
        movaps nb111nf_qqO(%esp),%xmm3
        mulss  %xmm4,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  %xmm0,%xmm3      ## xmm3=vcoul 
        mulps  nb111nf_c6(%esp),%xmm1
        mulps  nb111nf_c12(%esp),%xmm2
        movaps %xmm2,%xmm5
        subss  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addps  nb111nf_Vvdwtot(%esp),%xmm5
        addps  nb111nf_vctot(%esp),%xmm3
        movaps %xmm3,nb111nf_vctot(%esp)
        movaps %xmm5,nb111nf_Vvdwtot(%esp)

        decl nb111nf_innerk(%esp)
        jz    _nb_kernel111nf_ia32_sse.nb111nf_updateouterdata
        jmp   _nb_kernel111nf_ia32_sse.nb111nf_odd_loop
_nb_kernel111nf_ia32_sse.nb111nf_updateouterdata: 
        ## get n from stack
        movl nb111nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb111nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        movaps nb111nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb111nf_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb111nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb111nf_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb111nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel111nf_ia32_sse.nb111nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb111nf_n(%esp)
        jmp _nb_kernel111nf_ia32_sse.nb111nf_outer
_nb_kernel111nf_ia32_sse.nb111nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb111nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel111nf_ia32_sse.nb111nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel111nf_ia32_sse.nb111nf_threadloop
_nb_kernel111nf_ia32_sse.nb111nf_end: 
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




