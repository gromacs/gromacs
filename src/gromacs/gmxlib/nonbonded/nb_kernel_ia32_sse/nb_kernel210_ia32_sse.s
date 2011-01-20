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



.globl nb_kernel210_ia32_sse
.globl _nb_kernel210_ia32_sse
nb_kernel210_ia32_sse:  
_nb_kernel210_ia32_sse: 
.set nb210_p_nri, 8
.set nb210_iinr, 12
.set nb210_jindex, 16
.set nb210_jjnr, 20
.set nb210_shift, 24
.set nb210_shiftvec, 28
.set nb210_fshift, 32
.set nb210_gid, 36
.set nb210_pos, 40
.set nb210_faction, 44
.set nb210_charge, 48
.set nb210_p_facel, 52
.set nb210_argkrf, 56
.set nb210_argcrf, 60
.set nb210_Vc, 64
.set nb210_type, 68
.set nb210_p_ntype, 72
.set nb210_vdwparam, 76
.set nb210_Vvdw, 80
.set nb210_p_tabscale, 84
.set nb210_VFtab, 88
.set nb210_invsqrta, 92
.set nb210_dvda, 96
.set nb210_p_gbtabscale, 100
.set nb210_GBtab, 104
.set nb210_p_nthreads, 108
.set nb210_count, 112
.set nb210_mtx, 116
.set nb210_outeriter, 120
.set nb210_inneriter, 124
.set nb210_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb210_ix, 0
.set nb210_iy, 16
.set nb210_iz, 32
.set nb210_iq, 48
.set nb210_dx, 64
.set nb210_dy, 80
.set nb210_dz, 96
.set nb210_c6, 112
.set nb210_c12, 128
.set nb210_six, 144
.set nb210_twelve, 160
.set nb210_vctot, 176
.set nb210_Vvdwtot, 192
.set nb210_fix, 208
.set nb210_fiy, 224
.set nb210_fiz, 240
.set nb210_half, 256
.set nb210_three, 272
.set nb210_two, 288
.set nb210_krf, 304
.set nb210_crf, 320
.set nb210_is3, 336
.set nb210_ii3, 340
.set nb210_ntia, 344
.set nb210_innerjjnr, 348
.set nb210_innerk, 352
.set nb210_n, 356
.set nb210_nn1, 360
.set nb210_nri, 364
.set nb210_facel, 368
.set nb210_ntype, 372
.set nb210_nouter, 376
.set nb210_ninner, 380
.set nb210_salign, 384
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $400,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb210_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb210_p_nri(%ebp),%ecx
        movl nb210_p_facel(%ebp),%esi
        movl nb210_p_ntype(%ebp),%edi
        movl (%ecx),%ecx
        movl (%esi),%esi
        movl (%edi),%edi
        movl %ecx,nb210_nri(%esp)
        movl %esi,nb210_facel(%esp)
        movl %edi,nb210_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb210_nouter(%esp)
        movl %eax,nb210_ninner(%esp)


        movl nb210_argkrf(%ebp),%esi
        movl nb210_argcrf(%ebp),%edi
        movss (%esi),%xmm5
        movss (%edi),%xmm6
        shufps $0,%xmm5,%xmm5
        shufps $0,%xmm6,%xmm6
        movaps %xmm5,nb210_krf(%esp)
        movaps %xmm6,nb210_crf(%esp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb210_half(%esp)
        movss nb210_half(%esp),%xmm1
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
        movaps %xmm1,nb210_half(%esp)
        movaps %xmm2,nb210_two(%esp)
        movaps %xmm3,nb210_three(%esp)
        movaps %xmm4,nb210_six(%esp)
        movaps %xmm5,nb210_twelve(%esp)

_nb_kernel210_ia32_sse.nb210_threadloop: 
        movl  nb210_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel210_ia32_sse.nb210_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel210_ia32_sse.nb210_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb210_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb210_n(%esp)
        movl %ebx,nb210_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel210_ia32_sse.nb210_outerstart
        jmp _nb_kernel210_ia32_sse.nb210_end

_nb_kernel210_ia32_sse.nb210_outerstart: 
        ## ebx contains number of outer iterations
        addl nb210_nouter(%esp),%ebx
        movl %ebx,nb210_nouter(%esp)

_nb_kernel210_ia32_sse.nb210_outer: 
        movl  nb210_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb210_is3(%esp)      ## store is3 

        movl  nb210_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb210_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 

        movl  nb210_charge(%ebp),%edx
        movss (%edx,%ebx,4),%xmm3
        mulss nb210_facel(%esp),%xmm3
        shufps $0,%xmm3,%xmm3

        movl  nb210_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb210_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb210_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb210_pos(%ebp),%eax      ## eax = base of pos[]  

        addss (%eax,%ebx,4),%xmm0
        addss 4(%eax,%ebx,4),%xmm1
        addss 8(%eax,%ebx,4),%xmm2

        movaps %xmm3,nb210_iq(%esp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb210_ix(%esp)
        movaps %xmm1,nb210_iy(%esp)
        movaps %xmm2,nb210_iz(%esp)

        movl  %ebx,nb210_ii3(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb210_vctot(%esp)
        movaps %xmm4,nb210_Vvdwtot(%esp)
        movaps %xmm4,nb210_fix(%esp)
        movaps %xmm4,nb210_fiy(%esp)
        movaps %xmm4,nb210_fiz(%esp)

        movl  nb210_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb210_pos(%ebp),%esi
        movl  nb210_faction(%ebp),%edi
        movl  nb210_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb210_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb210_ninner(%esp),%ecx
        movl  %ecx,nb210_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb210_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel210_ia32_sse.nb210_unroll_loop
        jmp   _nb_kernel210_ia32_sse.nb210_finish_inner
_nb_kernel210_ia32_sse.nb210_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb210_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 
        addl $16,nb210_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb210_charge(%ebp),%esi     ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        movaps nb210_iq(%esp),%xmm2
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movl nb210_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb210_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb210_ntia(%esp),%edi
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

        movaps %xmm4,nb210_c6(%esp)
        movaps %xmm6,nb210_c12(%esp)

        movl nb210_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        mulps %xmm2,%xmm3
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

        ## move ix-iz to xmm4-xmm6 
        movaps nb210_ix(%esp),%xmm4
        movaps nb210_iy(%esp),%xmm5
        movaps nb210_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb210_dx(%esp)
        movaps %xmm5,nb210_dy(%esp)
        movaps %xmm6,nb210_dz(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        ## rsq in xmm4 

        movaps nb210_krf(%esp),%xmm7
        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb210_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb210_half(%esp),%xmm0
        mulps  %xmm4,%xmm7      ## xmm7=krsq 
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        movaps %xmm0,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvsq 
        movaps %xmm0,%xmm6
        addps  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        movaps %xmm4,%xmm1
        subps  nb210_crf(%esp),%xmm6
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  %xmm3,%xmm6      ## xmm6=vcoul=qq*(rinv+ krsq) 
        mulps  nb210_two(%esp),%xmm7
        mulps  nb210_c6(%esp),%xmm1
        mulps  nb210_c12(%esp),%xmm2
        movaps %xmm2,%xmm5
        subps  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addps  nb210_Vvdwtot(%esp),%xmm5
        mulps  nb210_six(%esp),%xmm1
        mulps  nb210_twelve(%esp),%xmm2
        subps  %xmm1,%xmm2
        subps  %xmm7,%xmm0
        mulps  %xmm0,%xmm3
        addps  %xmm3,%xmm2
        mulps  %xmm2,%xmm4      ## xmm4=total fscal 
        addps  nb210_vctot(%esp),%xmm6

        movaps nb210_dx(%esp),%xmm0
        movaps nb210_dy(%esp),%xmm1
        movaps nb210_dz(%esp),%xmm2

        movaps %xmm6,nb210_vctot(%esp)
        movaps %xmm5,nb210_Vvdwtot(%esp)

        movl   nb210_faction(%ebp),%edi
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movaps nb210_fix(%esp),%xmm3
        movaps nb210_fiy(%esp),%xmm4
        movaps nb210_fiz(%esp),%xmm5
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm5
        movaps %xmm3,nb210_fix(%esp)
        movaps %xmm4,nb210_fiy(%esp)
        movaps %xmm5,nb210_fiz(%esp)
        ## the fj's - start by accumulating x & y forces from memory 
        movlps (%edi,%eax,4),%xmm4
        movlps (%edi,%ecx,4),%xmm6
        movhps (%edi,%ebx,4),%xmm4
        movhps (%edi,%edx,4),%xmm6

        movaps %xmm4,%xmm3
        shufps $136,%xmm6,%xmm3 ## constant 10001000
        shufps $221,%xmm6,%xmm4 ## constant 11011101                          

        ## now xmm3-xmm5 contains fjx, fjy, fjz 
        subps  %xmm0,%xmm3
        subps  %xmm1,%xmm4

        ## unpack them back so we can store them - first x & y in xmm3/xmm4 

        movaps %xmm3,%xmm6
        unpcklps %xmm4,%xmm6
        unpckhps %xmm4,%xmm3
        ## xmm6(l)=x & y for j1, (h) for j2 
        ## xmm3(l)=x & y for j3, (h) for j4 
        movlps %xmm6,(%edi,%eax,4)
        movlps %xmm3,(%edi,%ecx,4)

        movhps %xmm6,(%edi,%ebx,4)
        movhps %xmm3,(%edi,%edx,4)

        ## and the z forces 
        movss  8(%edi,%eax,4),%xmm4
        movss  8(%edi,%ebx,4),%xmm5
        movss  8(%edi,%ecx,4),%xmm6
        movss  8(%edi,%edx,4),%xmm7
        subss  %xmm2,%xmm4
        shufps $229,%xmm2,%xmm2 ## constant 11100101
        subss  %xmm2,%xmm5
        shufps $234,%xmm2,%xmm2 ## constant 11101010
        subss  %xmm2,%xmm6
        shufps $255,%xmm2,%xmm2 ## constant 11111111
        subss  %xmm2,%xmm7
        movss  %xmm4,8(%edi,%eax,4)
        movss  %xmm5,8(%edi,%ebx,4)
        movss  %xmm6,8(%edi,%ecx,4)
        movss  %xmm7,8(%edi,%edx,4)

        ## should we do one more iteration? 
        subl $4,nb210_innerk(%esp)
        jl    _nb_kernel210_ia32_sse.nb210_finish_inner
        jmp   _nb_kernel210_ia32_sse.nb210_unroll_loop
_nb_kernel210_ia32_sse.nb210_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb210_innerk(%esp)
        movl  nb210_innerk(%esp),%edx
        andl  $2,%edx
        jnz   _nb_kernel210_ia32_sse.nb210_dopair
        jmp   _nb_kernel210_ia32_sse.nb210_checksingle
_nb_kernel210_ia32_sse.nb210_dopair: 
        movl nb210_charge(%ebp),%esi

    movl  nb210_innerjjnr(%esp),%ecx

        movl  (%ecx),%eax
        movl  4(%ecx),%ebx
        addl $8,nb210_innerjjnr(%esp)

        xorps %xmm3,%xmm3
        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ebx,4),%xmm6
        shufps $12,%xmm6,%xmm3 ## constant 00001100 
        shufps $88,%xmm3,%xmm3 ## constant 01011000 ;# xmm3(0,1) has the charges 

        movl nb210_type(%ebp),%esi
        movl  %eax,%ecx
        movl  %ebx,%edx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb210_vdwparam(%ebp),%esi
        shll %ecx
        shll %edx
        movl nb210_ntia(%esp),%edi
        addl %edi,%ecx
        addl %edi,%edx
        movlps (%esi,%ecx,4),%xmm6
        movhps (%esi,%edx,4),%xmm6
        movl nb210_pos(%ebp),%edi
        xorps  %xmm7,%xmm7
        movaps %xmm6,%xmm4
        shufps $8,%xmm4,%xmm4 ## constant 00001000       
        shufps $13,%xmm6,%xmm6 ## constant 00001101
        movlhps %xmm7,%xmm4
        movlhps %xmm7,%xmm6

        movaps %xmm4,nb210_c6(%esp)
        movaps %xmm6,nb210_c12(%esp)

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        ## move coordinates to xmm0-xmm2 
        movlps (%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2
        movhps (%edi,%ebx,4),%xmm1
        movss 8(%edi,%ebx,4),%xmm0

        mulps  nb210_iq(%esp),%xmm3

        movlhps %xmm7,%xmm3

        shufps $0,%xmm0,%xmm2

        movaps %xmm1,%xmm0

        shufps $136,%xmm2,%xmm2 ## constant 10001000

        shufps $136,%xmm0,%xmm0 ## constant 10001000
        shufps $221,%xmm1,%xmm1 ## constant 11011101

        movl   nb210_faction(%ebp),%edi
        ## move ix-iz to xmm4-xmm6 
        xorps   %xmm7,%xmm7

        movaps nb210_ix(%esp),%xmm4
        movaps nb210_iy(%esp),%xmm5
        movaps nb210_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb210_dx(%esp)
        movaps %xmm5,nb210_dy(%esp)
        movaps %xmm6,nb210_dz(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        ## rsq in xmm4 

        movaps nb210_krf(%esp),%xmm7
        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb210_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb210_half(%esp),%xmm0
        mulps  %xmm4,%xmm7      ## xmm7=krsq 
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        movaps %xmm0,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvsq 
        movaps %xmm0,%xmm6
        addps  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        movaps %xmm4,%xmm1
        subps  nb210_crf(%esp),%xmm6
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  %xmm3,%xmm6      ## xmm6=vcoul=qq*(rinv+ krsq-crf) 
        mulps  nb210_two(%esp),%xmm7
        mulps  nb210_c6(%esp),%xmm1
        mulps  nb210_c12(%esp),%xmm2
        movaps %xmm2,%xmm5
        subps  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addps  nb210_Vvdwtot(%esp),%xmm5
        mulps  nb210_six(%esp),%xmm1
        mulps  nb210_twelve(%esp),%xmm2
        subps  %xmm1,%xmm2
        subps  %xmm7,%xmm0
        mulps  %xmm0,%xmm3
        addps  %xmm3,%xmm2
        mulps  %xmm2,%xmm4      ## xmm4=total fscal 
        addps  nb210_vctot(%esp),%xmm6

        movaps nb210_dx(%esp),%xmm0
        movaps nb210_dy(%esp),%xmm1
        movaps nb210_dz(%esp),%xmm2

        movaps %xmm6,nb210_vctot(%esp)
        movaps %xmm5,nb210_Vvdwtot(%esp)

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movaps nb210_fix(%esp),%xmm3
        movaps nb210_fiy(%esp),%xmm4
        movaps nb210_fiz(%esp),%xmm5
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm5
        movaps %xmm3,nb210_fix(%esp)
        movaps %xmm4,nb210_fiy(%esp)
        movaps %xmm5,nb210_fiz(%esp)
        ## update the fj's 
        movss   (%edi,%eax,4),%xmm3
        movss   4(%edi,%eax,4),%xmm4
        movss   8(%edi,%eax,4),%xmm5
        subss   %xmm0,%xmm3
        subss   %xmm1,%xmm4
        subss   %xmm2,%xmm5
        movss   %xmm3,(%edi,%eax,4)
        movss   %xmm4,4(%edi,%eax,4)
        movss   %xmm5,8(%edi,%eax,4)

        shufps $225,%xmm0,%xmm0 ## constant 11100001
        shufps $225,%xmm1,%xmm1 ## constant 11100001
        shufps $225,%xmm2,%xmm2 ## constant 11100001

        movss   (%edi,%ebx,4),%xmm3
        movss   4(%edi,%ebx,4),%xmm4
        movss   8(%edi,%ebx,4),%xmm5
        subss   %xmm0,%xmm3
        subss   %xmm1,%xmm4
        subss   %xmm2,%xmm5
        movss   %xmm3,(%edi,%ebx,4)
        movss   %xmm4,4(%edi,%ebx,4)
        movss   %xmm5,8(%edi,%ebx,4)

_nb_kernel210_ia32_sse.nb210_checksingle:       
        movl  nb210_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel210_ia32_sse.nb210_dosingle
        jmp    _nb_kernel210_ia32_sse.nb210_updateouterdata
_nb_kernel210_ia32_sse.nb210_dosingle:  
        movl nb210_charge(%ebp),%esi
        movl nb210_pos(%ebp),%edi
        movl  nb210_innerjjnr(%esp),%ecx
        xorps %xmm3,%xmm3
        movl  (%ecx),%eax
        movss (%esi,%eax,4),%xmm3       ## xmm3(0) has the charge       

        movl nb210_type(%ebp),%esi
        movl %eax,%ecx
        movl (%esi,%ecx,4),%ecx
        movl nb210_vdwparam(%ebp),%esi
        shll %ecx
        addl nb210_ntia(%esp),%ecx
        xorps  %xmm6,%xmm6
        movlps (%esi,%ecx,4),%xmm6
        movaps %xmm6,%xmm4
        shufps $252,%xmm4,%xmm4 ## constant 11111100    
        shufps $253,%xmm6,%xmm6 ## constant 11111101    

        movaps %xmm4,nb210_c6(%esp)
        movaps %xmm6,nb210_c12(%esp)

        leal  (%eax,%eax,2),%eax

        ## move coordinates to xmm0-xmm2 
        movss (%edi,%eax,4),%xmm0
        movss 4(%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2

        mulps  nb210_iq(%esp),%xmm3

        xorps   %xmm7,%xmm7

        movaps nb210_ix(%esp),%xmm4
        movaps nb210_iy(%esp),%xmm5
        movaps nb210_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb210_dx(%esp)
        movaps %xmm5,nb210_dy(%esp)
        movaps %xmm6,nb210_dz(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        ## rsq in xmm4 

        movaps nb210_krf(%esp),%xmm7
        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb210_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb210_half(%esp),%xmm0
        mulps  %xmm4,%xmm7      ## xmm7=krsq 
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        movaps %xmm0,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvsq 
        movaps %xmm0,%xmm6
        addps  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        movaps %xmm4,%xmm1
        subps  nb210_crf(%esp),%xmm6
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  %xmm3,%xmm6      ## xmm6=vcoul 
        mulps  nb210_two(%esp),%xmm7
        mulps  nb210_c6(%esp),%xmm1
        mulps  nb210_c12(%esp),%xmm2
        movaps %xmm2,%xmm5
        subps  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addss  nb210_Vvdwtot(%esp),%xmm5
        mulps  nb210_six(%esp),%xmm1
        mulps  nb210_twelve(%esp),%xmm2
        subps  %xmm1,%xmm2
        subps  %xmm7,%xmm0
        mulps  %xmm0,%xmm3
        addps  %xmm3,%xmm2
        mulps  %xmm2,%xmm4      ## xmm4=total fscal 
        addss  nb210_vctot(%esp),%xmm6

        movl   nb210_faction(%ebp),%edi

        movaps nb210_dx(%esp),%xmm0
        movaps nb210_dy(%esp),%xmm1
        movaps nb210_dz(%esp),%xmm2

        movss %xmm6,nb210_vctot(%esp)
        movss %xmm5,nb210_Vvdwtot(%esp)

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movaps nb210_fix(%esp),%xmm3
        movaps nb210_fiy(%esp),%xmm4
        movaps nb210_fiz(%esp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movaps %xmm3,nb210_fix(%esp)
        movaps %xmm4,nb210_fiy(%esp)
        movaps %xmm5,nb210_fiz(%esp)
        ## update fj 

        movss   (%edi,%eax,4),%xmm3
        movss   4(%edi,%eax,4),%xmm4
        movss   8(%edi,%eax,4),%xmm5
        subss   %xmm0,%xmm3
        subss   %xmm1,%xmm4
        subss   %xmm2,%xmm5
        movss   %xmm3,(%edi,%eax,4)
        movss   %xmm4,4(%edi,%eax,4)
        movss   %xmm5,8(%edi,%eax,4)
_nb_kernel210_ia32_sse.nb210_updateouterdata: 
        movl  nb210_ii3(%esp),%ecx
        movl  nb210_faction(%ebp),%edi
        movl  nb210_fshift(%ebp),%esi
        movl  nb210_is3(%esp),%edx

        ## accumulate i forces in xmm0, xmm1, xmm2 
        movaps nb210_fix(%esp),%xmm0
        movaps nb210_fiy(%esp),%xmm1
        movaps nb210_fiz(%esp),%xmm2

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

        ## increment fshift force  
        movss  (%esi,%edx,4),%xmm3
        movss  4(%esi,%edx,4),%xmm4
        movss  8(%esi,%edx,4),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,(%esi,%edx,4)
        movss  %xmm4,4(%esi,%edx,4)
        movss  %xmm5,8(%esi,%edx,4)

        ## get n from stack
        movl nb210_n(%esp),%esi
        ## get group index for i particle 
        movl  nb210_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb210_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb210_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb210_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb210_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb210_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel210_ia32_sse.nb210_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb210_n(%esp)
        jmp _nb_kernel210_ia32_sse.nb210_outer
_nb_kernel210_ia32_sse.nb210_outerend: 
        ## check if more outer neighborlists remain
        movl  nb210_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel210_ia32_sse.nb210_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel210_ia32_sse.nb210_threadloop
_nb_kernel210_ia32_sse.nb210_end: 
        emms

        movl nb210_nouter(%esp),%eax
        movl nb210_ninner(%esp),%ebx
        movl nb210_outeriter(%ebp),%ecx
        movl nb210_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb210_salign(%esp),%eax
        addl %eax,%esp
        addl $400,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret







.globl nb_kernel210nf_ia32_sse
.globl _nb_kernel210nf_ia32_sse
nb_kernel210nf_ia32_sse:        
_nb_kernel210nf_ia32_sse:       
.set nb210nf_p_nri, 8
.set nb210nf_iinr, 12
.set nb210nf_jindex, 16
.set nb210nf_jjnr, 20
.set nb210nf_shift, 24
.set nb210nf_shiftvec, 28
.set nb210nf_fshift, 32
.set nb210nf_gid, 36
.set nb210nf_pos, 40
.set nb210nf_faction, 44
.set nb210nf_charge, 48
.set nb210nf_p_facel, 52
.set nb210nf_argkrf, 56
.set nb210nf_argcrf, 60
.set nb210nf_Vc, 64
.set nb210nf_type, 68
.set nb210nf_p_ntype, 72
.set nb210nf_vdwparam, 76
.set nb210nf_Vvdw, 80
.set nb210nf_p_tabscale, 84
.set nb210nf_VFtab, 88
.set nb210nf_invsqrta, 92
.set nb210nf_dvda, 96
.set nb210nf_p_gbtabscale, 100
.set nb210nf_GBtab, 104
.set nb210nf_p_nthreads, 108
.set nb210nf_count, 112
.set nb210nf_mtx, 116
.set nb210nf_outeriter, 120
.set nb210nf_inneriter, 124
.set nb210nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb210nf_ix, 0
.set nb210nf_iy, 16
.set nb210nf_iz, 32
.set nb210nf_iq, 48
.set nb210nf_c6, 64
.set nb210nf_c12, 80
.set nb210nf_vctot, 96
.set nb210nf_Vvdwtot, 112
.set nb210nf_half, 128
.set nb210nf_three, 144
.set nb210nf_krf, 160
.set nb210nf_crf, 176
.set nb210nf_is3, 192
.set nb210nf_ii3, 196
.set nb210nf_ntia, 200
.set nb210nf_innerjjnr, 204
.set nb210nf_innerk, 208
.set nb210nf_n, 212
.set nb210nf_nn1, 216
.set nb210nf_nri, 220
.set nb210nf_facel, 224
.set nb210nf_ntype, 228
.set nb210nf_nouter, 232
.set nb210nf_ninner, 236
.set nb210nf_salign, 240
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $256,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb210nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb210nf_p_nri(%ebp),%ecx
        movl nb210nf_p_facel(%ebp),%esi
        movl nb210nf_p_ntype(%ebp),%edi
        movl (%ecx),%ecx
        movl (%esi),%esi
        movl (%edi),%edi
        movl %ecx,nb210nf_nri(%esp)
        movl %esi,nb210nf_facel(%esp)
        movl %edi,nb210nf_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb210nf_nouter(%esp)
        movl %eax,nb210nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb210nf_half(%esp)
        movss nb210nf_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb210nf_half(%esp)
        movaps %xmm3,nb210nf_three(%esp)
        movl nb210nf_argkrf(%ebp),%esi
        movl nb210nf_argcrf(%ebp),%edi
        movss (%esi),%xmm5
        movss (%edi),%xmm6
        shufps $0,%xmm5,%xmm5
        shufps $0,%xmm6,%xmm6
        movaps %xmm5,nb210nf_krf(%esp)
        movaps %xmm6,nb210nf_crf(%esp)


_nb_kernel210nf_ia32_sse.nb210nf_threadloop: 
        movl  nb210nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel210nf_ia32_sse.nb210nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel210nf_ia32_sse.nb210nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb210nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb210nf_n(%esp)
        movl %ebx,nb210nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel210nf_ia32_sse.nb210nf_outerstart
        jmp _nb_kernel210nf_ia32_sse.nb210nf_end

_nb_kernel210nf_ia32_sse.nb210nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb210nf_nouter(%esp),%ebx
        movl %ebx,nb210nf_nouter(%esp)

_nb_kernel210nf_ia32_sse.nb210nf_outer: 
        movl  nb210nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb210nf_is3(%esp)            ## store is3 

        movl  nb210nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb210nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 

        movl  nb210nf_charge(%ebp),%edx
        movss (%edx,%ebx,4),%xmm3
        mulss nb210nf_facel(%esp),%xmm3
        shufps $0,%xmm3,%xmm3

        movl  nb210nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb210nf_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb210nf_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb210nf_pos(%ebp),%eax      ## eax = base of pos[]  

        addss (%eax,%ebx,4),%xmm0
        addss 4(%eax,%ebx,4),%xmm1
        addss 8(%eax,%ebx,4),%xmm2

        movaps %xmm3,nb210nf_iq(%esp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb210nf_ix(%esp)
        movaps %xmm1,nb210nf_iy(%esp)
        movaps %xmm2,nb210nf_iz(%esp)

        movl  %ebx,nb210nf_ii3(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb210nf_vctot(%esp)
        movaps %xmm4,nb210nf_Vvdwtot(%esp)

        movl  nb210nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb210nf_pos(%ebp),%esi
        movl  nb210nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb210nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb210nf_ninner(%esp),%ecx
        movl  %ecx,nb210nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb210nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel210nf_ia32_sse.nb210nf_unroll_loop
        jmp   _nb_kernel210nf_ia32_sse.nb210nf_finish_inner
_nb_kernel210nf_ia32_sse.nb210nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb210nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 
        addl $16,nb210nf_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb210nf_charge(%ebp),%esi     ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        movaps nb210nf_iq(%esp),%xmm2
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movl nb210nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb210nf_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb210nf_ntia(%esp),%edi
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

        movaps %xmm4,nb210nf_c6(%esp)
        movaps %xmm6,nb210nf_c12(%esp)

        movl nb210nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        mulps %xmm2,%xmm3
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

        ## move ix-iz to xmm4-xmm6 
        movaps nb210nf_ix(%esp),%xmm4
        movaps nb210nf_iy(%esp),%xmm5
        movaps nb210nf_iz(%esp),%xmm6

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
        ## rsq in xmm4 

        movaps nb210nf_krf(%esp),%xmm7
        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb210nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb210nf_half(%esp),%xmm0
        mulps  %xmm4,%xmm7      ## xmm7=krsq 
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        movaps %xmm0,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvsq 
        movaps %xmm0,%xmm6
        addps  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        movaps %xmm4,%xmm1
        subps  nb210nf_crf(%esp),%xmm6
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  %xmm3,%xmm6      ## xmm6=vcoul=qq*(rinv+ krsq) 
        mulps  nb210nf_c6(%esp),%xmm1
        mulps  nb210nf_c12(%esp),%xmm2
        movaps %xmm2,%xmm5
        subps  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addps  nb210nf_Vvdwtot(%esp),%xmm5
        addps  nb210nf_vctot(%esp),%xmm6
        movaps %xmm6,nb210nf_vctot(%esp)
        movaps %xmm5,nb210nf_Vvdwtot(%esp)

        ## should we do one more iteration? 
        subl $4,nb210nf_innerk(%esp)
        jl    _nb_kernel210nf_ia32_sse.nb210nf_finish_inner
        jmp   _nb_kernel210nf_ia32_sse.nb210nf_unroll_loop
_nb_kernel210nf_ia32_sse.nb210nf_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb210nf_innerk(%esp)
        movl  nb210nf_innerk(%esp),%edx
        andl  $2,%edx
        jnz   _nb_kernel210nf_ia32_sse.nb210nf_dopair
        jmp   _nb_kernel210nf_ia32_sse.nb210nf_checksingle
_nb_kernel210nf_ia32_sse.nb210nf_dopair: 
        movl nb210nf_charge(%ebp),%esi

        movl  nb210nf_innerjjnr(%esp),%ecx

        movl  (%ecx),%eax
        movl  4(%ecx),%ebx
        addl $8,nb210nf_innerjjnr(%esp)

        xorps %xmm3,%xmm3
        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ebx,4),%xmm6
        shufps $12,%xmm6,%xmm3 ## constant 00001100 
        shufps $88,%xmm3,%xmm3 ## constant 01011000 ;# xmm3(0,1) has the charges 

        movl nb210nf_type(%ebp),%esi
        movl  %eax,%ecx
        movl  %ebx,%edx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb210nf_vdwparam(%ebp),%esi
        shll %ecx
        shll %edx
        movl nb210nf_ntia(%esp),%edi
        addl %edi,%ecx
        addl %edi,%edx
        movlps (%esi,%ecx,4),%xmm6
        movhps (%esi,%edx,4),%xmm6
        movl nb210nf_pos(%ebp),%edi
        xorps  %xmm7,%xmm7
        movaps %xmm6,%xmm4
        shufps $8,%xmm4,%xmm4 ## constant 00001000       
        shufps $13,%xmm6,%xmm6 ## constant 00001101
        movlhps %xmm7,%xmm4
        movlhps %xmm7,%xmm6

        movaps %xmm4,nb210nf_c6(%esp)
        movaps %xmm6,nb210nf_c12(%esp)

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        ## move coordinates to xmm0-xmm2 
        movlps (%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2
        movhps (%edi,%ebx,4),%xmm1
        movss 8(%edi,%ebx,4),%xmm0

        mulps  nb210nf_iq(%esp),%xmm3

        movlhps %xmm7,%xmm3

        shufps $0,%xmm0,%xmm2

        movaps %xmm1,%xmm0

        shufps $136,%xmm2,%xmm2 ## constant 10001000

        shufps $136,%xmm0,%xmm0 ## constant 10001000
        shufps $221,%xmm1,%xmm1 ## constant 11011101

        ## move ix-iz to xmm4-xmm6 
        xorps   %xmm7,%xmm7

        movaps nb210nf_ix(%esp),%xmm4
        movaps nb210nf_iy(%esp),%xmm5
        movaps nb210nf_iz(%esp),%xmm6

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
        ## rsq in xmm4 

        movaps nb210nf_krf(%esp),%xmm7
        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb210nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb210nf_half(%esp),%xmm0
        mulps  %xmm4,%xmm7      ## xmm7=krsq 
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        movaps %xmm0,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvsq 
        movaps %xmm0,%xmm6
        addps  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        movaps %xmm4,%xmm1
        subps  nb210nf_crf(%esp),%xmm6
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  %xmm3,%xmm6      ## xmm6=vcoul=qq*(rinv+ krsq-crf) 
        mulps  nb210nf_c6(%esp),%xmm1
        mulps  nb210nf_c12(%esp),%xmm2
        movaps %xmm2,%xmm5
        subps  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addps  nb210nf_Vvdwtot(%esp),%xmm5
        addps  nb210nf_vctot(%esp),%xmm6
        movaps %xmm6,nb210nf_vctot(%esp)
        movaps %xmm5,nb210nf_Vvdwtot(%esp)

_nb_kernel210nf_ia32_sse.nb210nf_checksingle:   
        movl  nb210nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel210nf_ia32_sse.nb210nf_dosingle
        jmp    _nb_kernel210nf_ia32_sse.nb210nf_updateouterdata
_nb_kernel210nf_ia32_sse.nb210nf_dosingle: 
        movl nb210nf_charge(%ebp),%esi
        movl nb210nf_pos(%ebp),%edi
        movl  nb210nf_innerjjnr(%esp),%ecx
        xorps %xmm3,%xmm3
        movl  (%ecx),%eax
        movss (%esi,%eax,4),%xmm3       ## xmm3(0) has the charge       

        movl nb210nf_type(%ebp),%esi
        movl %eax,%ecx
        movl (%esi,%ecx,4),%ecx
        movl nb210nf_vdwparam(%ebp),%esi
        shll %ecx
        addl nb210nf_ntia(%esp),%ecx
        xorps  %xmm6,%xmm6
        movlps (%esi,%ecx,4),%xmm6
        movaps %xmm6,%xmm4
        shufps $252,%xmm4,%xmm4 ## constant 11111100    
        shufps $253,%xmm6,%xmm6 ## constant 11111101    

        movaps %xmm4,nb210nf_c6(%esp)
        movaps %xmm6,nb210nf_c12(%esp)

        leal  (%eax,%eax,2),%eax

        ## move coordinates to xmm0-xmm2 
        movss (%edi,%eax,4),%xmm0
        movss 4(%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2

        mulps  nb210nf_iq(%esp),%xmm3

        xorps   %xmm7,%xmm7

        movaps nb210nf_ix(%esp),%xmm4
        movaps nb210nf_iy(%esp),%xmm5
        movaps nb210nf_iz(%esp),%xmm6

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
        ## rsq in xmm4 

        movaps nb210nf_krf(%esp),%xmm7
        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb210nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb210nf_half(%esp),%xmm0
        mulps  %xmm4,%xmm7      ## xmm7=krsq 
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        movaps %xmm0,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvsq 
        movaps %xmm0,%xmm6
        addps  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        movaps %xmm4,%xmm1
        subps  nb210nf_crf(%esp),%xmm6
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  %xmm3,%xmm6      ## xmm6=vcoul 
        mulps  nb210nf_c6(%esp),%xmm1
        mulps  nb210nf_c12(%esp),%xmm2
        movaps %xmm2,%xmm5
        subps  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addss  nb210nf_Vvdwtot(%esp),%xmm5
        addss  nb210nf_vctot(%esp),%xmm6
        movss %xmm6,nb210nf_vctot(%esp)
        movss %xmm5,nb210nf_Vvdwtot(%esp)

_nb_kernel210nf_ia32_sse.nb210nf_updateouterdata: 
        ## get n from stack
        movl nb210nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb210nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb210nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb210nf_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb210nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb210nf_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb210nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel210nf_ia32_sse.nb210nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb210nf_n(%esp)
        jmp _nb_kernel210nf_ia32_sse.nb210nf_outer
_nb_kernel210nf_ia32_sse.nb210nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb210nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel210nf_ia32_sse.nb210nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel210nf_ia32_sse.nb210nf_threadloop
_nb_kernel210nf_ia32_sse.nb210nf_end: 
        emms

        movl nb210nf_nouter(%esp),%eax
        movl nb210nf_ninner(%esp),%ebx
        movl nb210nf_outeriter(%ebp),%ecx
        movl nb210nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb210nf_salign(%esp),%eax
        addl %eax,%esp
        addl $256,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret


