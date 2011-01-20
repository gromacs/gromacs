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




.globl nb_kernel110_ia32_sse
.globl _nb_kernel110_ia32_sse
nb_kernel110_ia32_sse:  
_nb_kernel110_ia32_sse: 
.set nb110_p_nri, 8
.set nb110_iinr, 12
.set nb110_jindex, 16
.set nb110_jjnr, 20
.set nb110_shift, 24
.set nb110_shiftvec, 28
.set nb110_fshift, 32
.set nb110_gid, 36
.set nb110_pos, 40
.set nb110_faction, 44
.set nb110_charge, 48
.set nb110_p_facel, 52
.set nb110_p_krf, 56
.set nb110_p_crf, 60
.set nb110_Vc, 64
.set nb110_type, 68
.set nb110_p_ntype, 72
.set nb110_vdwparam, 76
.set nb110_Vvdw, 80
.set nb110_p_tabscale, 84
.set nb110_VFtab, 88
.set nb110_invsqrta, 92
.set nb110_dvda, 96
.set nb110_p_gbtabscale, 100
.set nb110_GBtab, 104
.set nb110_p_nthreads, 108
.set nb110_count, 112
.set nb110_mtx, 116
.set nb110_outeriter, 120
.set nb110_inneriter, 124
.set nb110_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb110_ix, 0
.set nb110_iy, 16
.set nb110_iz, 32
.set nb110_iq, 48
.set nb110_dx, 64
.set nb110_dy, 80
.set nb110_dz, 96
.set nb110_c6, 112
.set nb110_c12, 128
.set nb110_six, 144
.set nb110_twelve, 160
.set nb110_vctot, 176
.set nb110_Vvdwtot, 192
.set nb110_fix, 208
.set nb110_fiy, 224
.set nb110_fiz, 240
.set nb110_half, 256
.set nb110_three, 272
.set nb110_is3, 288
.set nb110_ii3, 292
.set nb110_ntia, 296
.set nb110_innerjjnr, 300
.set nb110_innerk, 304
.set nb110_n, 308
.set nb110_nn1, 312
.set nb110_nri, 316
.set nb110_facel, 320
.set nb110_ntype, 324
.set nb110_nouter, 328
.set nb110_ninner, 332
.set nb110_salign, 336
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $340,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb110_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb110_p_nri(%ebp),%ecx
        movl nb110_p_facel(%ebp),%esi
        movl nb110_p_ntype(%ebp),%edi
        movl (%ecx),%ecx
        movl (%esi),%esi
        movl (%edi),%edi
        movl %ecx,nb110_nri(%esp)
        movl %esi,nb110_facel(%esp)
        movl %edi,nb110_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb110_nouter(%esp)
        movl %eax,nb110_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb110_half(%esp)
        movss nb110_half(%esp),%xmm1
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
        movaps %xmm1,nb110_half(%esp)
        movaps %xmm3,nb110_three(%esp)
        movaps %xmm4,nb110_six(%esp)
        movaps %xmm5,nb110_twelve(%esp)

_nb_kernel110_ia32_sse.nb110_threadloop: 
        movl  nb110_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel110_ia32_sse.nb110_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel110_ia32_sse.nb110_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb110_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb110_n(%esp)
        movl %ebx,nb110_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel110_ia32_sse.nb110_outerstart
        jmp _nb_kernel110_ia32_sse.nb110_end

_nb_kernel110_ia32_sse.nb110_outerstart: 
        ## ebx contains number of outer iterations
        addl nb110_nouter(%esp),%ebx
        movl %ebx,nb110_nouter(%esp)

_nb_kernel110_ia32_sse.nb110_outer: 
        movl  nb110_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb110_is3(%esp)      ## store is3 

        movl  nb110_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb110_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movl  nb110_charge(%ebp),%edx
        movss (%edx,%ebx,4),%xmm3
        mulss nb110_facel(%esp),%xmm3
        shufps $0,%xmm3,%xmm3

        movl  nb110_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb110_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb110_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb110_pos(%ebp),%eax      ## eax = base of pos[]  

        addss (%eax,%ebx,4),%xmm0
        addss 4(%eax,%ebx,4),%xmm1
        addss 8(%eax,%ebx,4),%xmm2

        movaps %xmm3,nb110_iq(%esp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb110_ix(%esp)
        movaps %xmm1,nb110_iy(%esp)
        movaps %xmm2,nb110_iz(%esp)

        movl  %ebx,nb110_ii3(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb110_vctot(%esp)
        movaps %xmm4,nb110_Vvdwtot(%esp)
        movaps %xmm4,nb110_fix(%esp)
        movaps %xmm4,nb110_fiy(%esp)
        movaps %xmm4,nb110_fiz(%esp)

        movl  nb110_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb110_pos(%ebp),%esi
        movl  nb110_faction(%ebp),%edi
        movl  nb110_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb110_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb110_ninner(%esp),%ecx
        movl  %ecx,nb110_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb110_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel110_ia32_sse.nb110_unroll_loop
        jmp   _nb_kernel110_ia32_sse.nb110_finish_inner
_nb_kernel110_ia32_sse.nb110_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb110_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 
        addl $16,nb110_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb110_charge(%ebp),%esi     ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        movaps nb110_iq(%esp),%xmm2
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movl nb110_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb110_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb110_ntia(%esp),%edi
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

        movaps %xmm4,nb110_c6(%esp)
        movaps %xmm6,nb110_c12(%esp)

        movl nb110_pos(%ebp),%esi        ## base of pos[] 

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
        movaps nb110_ix(%esp),%xmm4
        movaps nb110_iy(%esp),%xmm5
        movaps nb110_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb110_dx(%esp)
        movaps %xmm5,nb110_dy(%esp)
        movaps %xmm6,nb110_dz(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        ## rsq in xmm4 

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb110_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb110_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        movaps %xmm0,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvsq 
        movaps %xmm4,%xmm1
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  %xmm0,%xmm3      ## xmm3=vcoul 
        mulps  nb110_c6(%esp),%xmm1
        mulps  nb110_c12(%esp),%xmm2
        movaps %xmm2,%xmm5
        subps  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addps  nb110_Vvdwtot(%esp),%xmm5
        mulps  nb110_six(%esp),%xmm1
        mulps  nb110_twelve(%esp),%xmm2
        subps  %xmm1,%xmm2
        addps  %xmm3,%xmm2
        mulps  %xmm2,%xmm4      ## xmm4=total fscal 
        addps  nb110_vctot(%esp),%xmm3

        movaps nb110_dx(%esp),%xmm0
        movaps nb110_dy(%esp),%xmm1
        movaps nb110_dz(%esp),%xmm2

        movaps %xmm3,nb110_vctot(%esp)
        movaps %xmm5,nb110_Vvdwtot(%esp)

        movl   nb110_faction(%ebp),%edi
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movaps nb110_fix(%esp),%xmm3
        movaps nb110_fiy(%esp),%xmm4
        movaps nb110_fiz(%esp),%xmm5
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm5
        movaps %xmm3,nb110_fix(%esp)
        movaps %xmm4,nb110_fiy(%esp)
        movaps %xmm5,nb110_fiz(%esp)
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
        subl $4,nb110_innerk(%esp)
        jl    _nb_kernel110_ia32_sse.nb110_finish_inner
        jmp   _nb_kernel110_ia32_sse.nb110_unroll_loop
_nb_kernel110_ia32_sse.nb110_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb110_innerk(%esp)
        movl  nb110_innerk(%esp),%edx
        andl  $2,%edx
        jnz   _nb_kernel110_ia32_sse.nb110_dopair
        jmp   _nb_kernel110_ia32_sse.nb110_checksingle
_nb_kernel110_ia32_sse.nb110_dopair: 
        movl nb110_charge(%ebp),%esi

    movl  nb110_innerjjnr(%esp),%ecx

        movl  (%ecx),%eax
        movl  4(%ecx),%ebx
        addl $8,nb110_innerjjnr(%esp)

        xorps %xmm3,%xmm3
        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ebx,4),%xmm6
        shufps $12,%xmm6,%xmm3 ## constant 00001100 
        shufps $88,%xmm3,%xmm3 ## constant 01011000 ;# xmm3(0,1) has the charges 

        movl nb110_type(%ebp),%esi
        movl  %eax,%ecx
        movl  %ebx,%edx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb110_vdwparam(%ebp),%esi
        shll %ecx
        shll %edx
        movl nb110_ntia(%esp),%edi
        addl %edi,%ecx
        addl %edi,%edx
        movlps (%esi,%ecx,4),%xmm6
        movhps (%esi,%edx,4),%xmm6
        movl nb110_pos(%ebp),%edi
        xorps  %xmm7,%xmm7
        movaps %xmm6,%xmm4
        shufps $8,%xmm4,%xmm4 ## constant 00001000       
        shufps $13,%xmm6,%xmm6 ## constant 00001101
        movlhps %xmm7,%xmm4
        movlhps %xmm7,%xmm6

        movaps %xmm4,nb110_c6(%esp)
        movaps %xmm6,nb110_c12(%esp)

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        ## move coordinates to xmm0-xmm2 
        movlps (%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2
        movhps (%edi,%ebx,4),%xmm1
        movss 8(%edi,%ebx,4),%xmm0

        mulps  nb110_iq(%esp),%xmm3

        movlhps %xmm7,%xmm3

        shufps $0,%xmm0,%xmm2

        movaps %xmm1,%xmm0

        shufps $136,%xmm2,%xmm2 ## constant 10001000

        shufps $136,%xmm0,%xmm0 ## constant 10001000
        shufps $221,%xmm1,%xmm1 ## constant 11011101

        movl   nb110_faction(%ebp),%edi
        ## move ix-iz to xmm4-xmm6 
        xorps   %xmm7,%xmm7

        movaps nb110_ix(%esp),%xmm4
        movaps nb110_iy(%esp),%xmm5
        movaps nb110_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb110_dx(%esp)
        movaps %xmm5,nb110_dy(%esp)
        movaps %xmm6,nb110_dz(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        ## rsq in xmm4 

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb110_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb110_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        movaps %xmm0,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvsq 
        movaps %xmm4,%xmm1
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 

        mulps  %xmm0,%xmm3      ## xmm3=vcoul 
        mulps  nb110_c6(%esp),%xmm1
        mulps  nb110_c12(%esp),%xmm2
        movaps %xmm2,%xmm5
        subps  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addps  nb110_Vvdwtot(%esp),%xmm5
        mulps  nb110_six(%esp),%xmm1
        mulps  nb110_twelve(%esp),%xmm2
        subps  %xmm1,%xmm2
        addps  %xmm3,%xmm2
        mulps  %xmm2,%xmm4      ## xmm4=total fscal 
        addps  nb110_vctot(%esp),%xmm3

        movaps nb110_dx(%esp),%xmm0
        movaps nb110_dy(%esp),%xmm1
        movaps nb110_dz(%esp),%xmm2

        movaps %xmm3,nb110_vctot(%esp)
        movaps %xmm5,nb110_Vvdwtot(%esp)

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movaps nb110_fix(%esp),%xmm3
        movaps nb110_fiy(%esp),%xmm4
        movaps nb110_fiz(%esp),%xmm5
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm5
        movaps %xmm3,nb110_fix(%esp)
        movaps %xmm4,nb110_fiy(%esp)
        movaps %xmm5,nb110_fiz(%esp)
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

_nb_kernel110_ia32_sse.nb110_checksingle:       
        movl  nb110_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel110_ia32_sse.nb110_dosingle
        jmp    _nb_kernel110_ia32_sse.nb110_updateouterdata
_nb_kernel110_ia32_sse.nb110_dosingle:  
        movl nb110_charge(%ebp),%esi
        movl nb110_pos(%ebp),%edi
        movl  nb110_innerjjnr(%esp),%ecx
        xorps %xmm3,%xmm3
        movl  (%ecx),%eax
        movss (%esi,%eax,4),%xmm3       ## xmm3(0) has the charge       

        movl nb110_type(%ebp),%esi
        movl %eax,%ecx
        movl (%esi,%ecx,4),%ecx
        movl nb110_vdwparam(%ebp),%esi
        shll %ecx
        addl nb110_ntia(%esp),%ecx
        xorps  %xmm6,%xmm6
        movlps (%esi,%ecx,4),%xmm6
        movaps %xmm6,%xmm4
        shufps $252,%xmm4,%xmm4 ## constant 11111100    
        shufps $253,%xmm6,%xmm6 ## constant 11111101    

        movaps %xmm4,nb110_c6(%esp)
        movaps %xmm6,nb110_c12(%esp)

        leal  (%eax,%eax,2),%eax

        ## move coordinates to xmm0-xmm2 
        movss (%edi,%eax,4),%xmm0
        movss 4(%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2

        mulps  nb110_iq(%esp),%xmm3

        xorps   %xmm7,%xmm7

        movaps nb110_ix(%esp),%xmm4
        movaps nb110_iy(%esp),%xmm5
        movaps nb110_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb110_dx(%esp)
        movaps %xmm5,nb110_dy(%esp)
        movaps %xmm6,nb110_dz(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        ## rsq in xmm4 

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb110_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb110_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        movaps %xmm0,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvsq 
        movaps %xmm4,%xmm1
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  %xmm0,%xmm3      ## xmm3=vcoul 
        mulps  nb110_c6(%esp),%xmm1
        mulps  nb110_c12(%esp),%xmm2
        movaps %xmm2,%xmm5
        subps  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addss  nb110_Vvdwtot(%esp),%xmm5
        mulps  nb110_six(%esp),%xmm1
        mulps  nb110_twelve(%esp),%xmm2
        subps  %xmm1,%xmm2
        addps  %xmm3,%xmm2
        mulps  %xmm2,%xmm4      ## xmm4=total fscal 
        addss  nb110_vctot(%esp),%xmm3

        movl   nb110_faction(%ebp),%edi

        movaps nb110_dx(%esp),%xmm0
        movaps nb110_dy(%esp),%xmm1
        movaps nb110_dz(%esp),%xmm2

        movss %xmm3,nb110_vctot(%esp)
        movss %xmm5,nb110_Vvdwtot(%esp)

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movaps nb110_fix(%esp),%xmm3
        movaps nb110_fiy(%esp),%xmm4
        movaps nb110_fiz(%esp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movaps %xmm3,nb110_fix(%esp)
        movaps %xmm4,nb110_fiy(%esp)
        movaps %xmm5,nb110_fiz(%esp)
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
_nb_kernel110_ia32_sse.nb110_updateouterdata: 
        movl  nb110_ii3(%esp),%ecx
        movl  nb110_faction(%ebp),%edi
        movl  nb110_fshift(%ebp),%esi
        movl  nb110_is3(%esp),%edx

        ## accumulate i forces in xmm0, xmm1, xmm2 
        movaps nb110_fix(%esp),%xmm0
        movaps nb110_fiy(%esp),%xmm1
        movaps nb110_fiz(%esp),%xmm2

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
        movl nb110_n(%esp),%esi
        ## get group index for i particle 
        movl  nb110_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb110_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb110_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb110_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb110_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb110_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel110_ia32_sse.nb110_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb110_n(%esp)
        jmp _nb_kernel110_ia32_sse.nb110_outer
_nb_kernel110_ia32_sse.nb110_outerend: 
        ## check if more outer neighborlists remain
        movl  nb110_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel110_ia32_sse.nb110_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel110_ia32_sse.nb110_threadloop
_nb_kernel110_ia32_sse.nb110_end: 
        emms

        movl nb110_nouter(%esp),%eax
        movl nb110_ninner(%esp),%ebx
        movl nb110_outeriter(%ebp),%ecx
        movl nb110_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb110_salign(%esp),%eax
        addl %eax,%esp
        addl $340,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret



.globl nb_kernel110nf_ia32_sse
.globl _nb_kernel110nf_ia32_sse
nb_kernel110nf_ia32_sse:        
_nb_kernel110nf_ia32_sse:       
.set nb110nf_p_nri, 8
.set nb110nf_iinr, 12
.set nb110nf_jindex, 16
.set nb110nf_jjnr, 20
.set nb110nf_shift, 24
.set nb110nf_shiftvec, 28
.set nb110nf_fshift, 32
.set nb110nf_gid, 36
.set nb110nf_pos, 40
.set nb110nf_faction, 44
.set nb110nf_charge, 48
.set nb110nf_p_facel, 52
.set nb110nf_p_krf, 56
.set nb110nf_p_crf, 60
.set nb110nf_Vc, 64
.set nb110nf_type, 68
.set nb110nf_p_ntype, 72
.set nb110nf_vdwparam, 76
.set nb110nf_Vvdw, 80
.set nb110nf_p_tabscale, 84
.set nb110nf_VFtab, 88
.set nb110nf_invsqrta, 92
.set nb110nf_dvda, 96
.set nb110nf_p_gbtabscale, 100
.set nb110nf_GBtab, 104
.set nb110nf_p_nthreads, 108
.set nb110nf_count, 112
.set nb110nf_mtx, 116
.set nb110nf_outeriter, 120
.set nb110nf_inneriter, 124
.set nb110nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb110nf_ix, 0
.set nb110nf_iy, 16
.set nb110nf_iz, 32
.set nb110nf_iq, 48
.set nb110nf_c6, 64
.set nb110nf_c12, 80
.set nb110nf_vctot, 96
.set nb110nf_Vvdwtot, 112
.set nb110nf_half, 128
.set nb110nf_three, 144
.set nb110nf_is3, 160
.set nb110nf_ii3, 164
.set nb110nf_ntia, 168
.set nb110nf_innerjjnr, 172
.set nb110nf_innerk, 176
.set nb110nf_n, 180
.set nb110nf_nn1, 184
.set nb110nf_nri, 188
.set nb110nf_facel, 192
.set nb110nf_ntype, 196
.set nb110nf_nouter, 200
.set nb110nf_ninner, 204
.set nb110nf_salign, 208
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $216,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb110nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb110nf_p_nri(%ebp),%ecx
        movl nb110nf_p_facel(%ebp),%esi
        movl nb110nf_p_ntype(%ebp),%edi
        movl (%ecx),%ecx
        movl (%esi),%esi
        movl (%edi),%edi
        movl %ecx,nb110nf_nri(%esp)
        movl %esi,nb110nf_facel(%esp)
        movl %edi,nb110nf_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb110nf_nouter(%esp)
        movl %eax,nb110nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb110nf_half(%esp)
        movss nb110nf_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb110nf_half(%esp)
        movaps %xmm3,nb110nf_three(%esp)

_nb_kernel110nf_ia32_sse.nb110nf_threadloop: 
        movl  nb110nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel110nf_ia32_sse.nb110nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel110nf_ia32_sse.nb110nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb110nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb110nf_n(%esp)
        movl %ebx,nb110nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel110nf_ia32_sse.nb110nf_outerstart
        jmp _nb_kernel110nf_ia32_sse.nb110nf_end

_nb_kernel110nf_ia32_sse.nb110nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb110nf_nouter(%esp),%ebx
        movl %ebx,nb110nf_nouter(%esp)

_nb_kernel110nf_ia32_sse.nb110nf_outer: 
        movl  nb110nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb110nf_is3(%esp)            ## store is3 

        movl  nb110nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb110nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movl  nb110nf_charge(%ebp),%edx
        movss (%edx,%ebx,4),%xmm3
        mulss nb110nf_facel(%esp),%xmm3
        shufps $0,%xmm3,%xmm3

        movl  nb110nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb110nf_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb110nf_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb110nf_pos(%ebp),%eax      ## eax = base of pos[]  

        addss (%eax,%ebx,4),%xmm0
        addss 4(%eax,%ebx,4),%xmm1
        addss 8(%eax,%ebx,4),%xmm2

        movaps %xmm3,nb110nf_iq(%esp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb110nf_ix(%esp)
        movaps %xmm1,nb110nf_iy(%esp)
        movaps %xmm2,nb110nf_iz(%esp)

        movl  %ebx,nb110nf_ii3(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb110nf_vctot(%esp)
        movaps %xmm4,nb110nf_Vvdwtot(%esp)

        movl  nb110nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb110nf_pos(%ebp),%esi
        movl  nb110nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb110nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb110nf_ninner(%esp),%ecx
        movl  %ecx,nb110nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb110nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel110nf_ia32_sse.nb110nf_unroll_loop
        jmp   _nb_kernel110nf_ia32_sse.nb110nf_finish_inner
_nb_kernel110nf_ia32_sse.nb110nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb110nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 
        addl $16,nb110nf_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb110nf_charge(%ebp),%esi     ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        movaps nb110nf_iq(%esp),%xmm2
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movl nb110nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb110nf_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb110nf_ntia(%esp),%edi
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

        movaps %xmm4,nb110nf_c6(%esp)
        movaps %xmm6,nb110nf_c12(%esp)

        movl nb110nf_pos(%ebp),%esi        ## base of pos[] 

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
        movaps nb110nf_ix(%esp),%xmm4
        movaps nb110nf_iy(%esp),%xmm5
        movaps nb110nf_iz(%esp),%xmm6

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

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb110nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb110nf_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        movaps %xmm0,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvsq 
        movaps %xmm4,%xmm1
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  %xmm0,%xmm3      ## xmm3=vcoul 
        mulps  nb110nf_c6(%esp),%xmm1
        mulps  nb110nf_c12(%esp),%xmm2
        movaps %xmm2,%xmm5
        subps  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addps  nb110nf_Vvdwtot(%esp),%xmm5
        addps  nb110nf_vctot(%esp),%xmm3
        movaps %xmm3,nb110nf_vctot(%esp)
        movaps %xmm5,nb110nf_Vvdwtot(%esp)

        ## should we do one more iteration? 
        subl $4,nb110nf_innerk(%esp)
        jl    _nb_kernel110nf_ia32_sse.nb110nf_finish_inner
        jmp   _nb_kernel110nf_ia32_sse.nb110nf_unroll_loop
_nb_kernel110nf_ia32_sse.nb110nf_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb110nf_innerk(%esp)
        movl  nb110nf_innerk(%esp),%edx
        andl  $2,%edx
        jnz   _nb_kernel110nf_ia32_sse.nb110nf_dopair
        jmp   _nb_kernel110nf_ia32_sse.nb110nf_checksingle
_nb_kernel110nf_ia32_sse.nb110nf_dopair: 
        movl nb110nf_charge(%ebp),%esi

        movl  nb110nf_innerjjnr(%esp),%ecx

        movl  (%ecx),%eax
        movl  4(%ecx),%ebx
        addl $8,nb110nf_innerjjnr(%esp)

        xorps %xmm3,%xmm3
        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ebx,4),%xmm6
        shufps $12,%xmm6,%xmm3 ## constant 00001100 
        shufps $88,%xmm3,%xmm3 ## constant 01011000 ;# xmm3(0,1) has the charges 

        movl nb110nf_type(%ebp),%esi
        movl  %eax,%ecx
        movl  %ebx,%edx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb110nf_vdwparam(%ebp),%esi
        shll %ecx
        shll %edx
        movl nb110nf_ntia(%esp),%edi
        addl %edi,%ecx
        addl %edi,%edx
        movlps (%esi,%ecx,4),%xmm6
        movhps (%esi,%edx,4),%xmm6
        movl nb110nf_pos(%ebp),%edi
        xorps  %xmm7,%xmm7
        movaps %xmm6,%xmm4
        shufps $8,%xmm4,%xmm4 ## constant 00001000       
        shufps $13,%xmm6,%xmm6 ## constant 00001101
        movlhps %xmm7,%xmm4
        movlhps %xmm7,%xmm6

        movaps %xmm4,nb110nf_c6(%esp)
        movaps %xmm6,nb110nf_c12(%esp)

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        ## move coordinates to xmm0-xmm2 
        movlps (%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2
        movhps (%edi,%ebx,4),%xmm1
        movss 8(%edi,%ebx,4),%xmm0

        mulps  nb110nf_iq(%esp),%xmm3

        movlhps %xmm7,%xmm3

        shufps $0,%xmm0,%xmm2

        movaps %xmm1,%xmm0

        shufps $136,%xmm2,%xmm2 ## constant 10001000

        shufps $136,%xmm0,%xmm0 ## constant 10001000
        shufps $221,%xmm1,%xmm1 ## constant 11011101

        ## move ix-iz to xmm4-xmm6 
        xorps   %xmm7,%xmm7

        movaps nb110nf_ix(%esp),%xmm4
        movaps nb110nf_iy(%esp),%xmm5
        movaps nb110nf_iz(%esp),%xmm6

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

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb110nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb110nf_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        movaps %xmm0,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvsq 
        movaps %xmm4,%xmm1
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 

        mulps  %xmm0,%xmm3      ## xmm3=vcoul 
        mulps  nb110nf_c6(%esp),%xmm1
        mulps  nb110nf_c12(%esp),%xmm2
        movaps %xmm2,%xmm5
        subps  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addps  nb110nf_Vvdwtot(%esp),%xmm5
        addps  nb110nf_vctot(%esp),%xmm3
        movaps %xmm3,nb110nf_vctot(%esp)
        movaps %xmm5,nb110nf_Vvdwtot(%esp)

_nb_kernel110nf_ia32_sse.nb110nf_checksingle:   
        movl  nb110nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel110nf_ia32_sse.nb110nf_dosingle
        jmp    _nb_kernel110nf_ia32_sse.nb110nf_updateouterdata
_nb_kernel110nf_ia32_sse.nb110nf_dosingle: 
        movl nb110nf_charge(%ebp),%esi
        movl nb110nf_pos(%ebp),%edi
        movl  nb110nf_innerjjnr(%esp),%ecx
        xorps %xmm3,%xmm3
        movl  (%ecx),%eax
        movss (%esi,%eax,4),%xmm3       ## xmm3(0) has the charge       

        movl nb110nf_type(%ebp),%esi
        movl %eax,%ecx
        movl (%esi,%ecx,4),%ecx
        movl nb110nf_vdwparam(%ebp),%esi
        shll %ecx
        addl nb110nf_ntia(%esp),%ecx
        xorps  %xmm6,%xmm6
        movlps (%esi,%ecx,4),%xmm6
        movaps %xmm6,%xmm4
        shufps $252,%xmm4,%xmm4 ## constant 11111100    
        shufps $253,%xmm6,%xmm6 ## constant 11111101    

        movaps %xmm4,nb110nf_c6(%esp)
        movaps %xmm6,nb110nf_c12(%esp)

        leal  (%eax,%eax,2),%eax

        ## move coordinates to xmm0-xmm2 
        movss (%edi,%eax,4),%xmm0
        movss 4(%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2

        mulps  nb110nf_iq(%esp),%xmm3

        xorps   %xmm7,%xmm7

        movaps nb110nf_ix(%esp),%xmm4
        movaps nb110nf_iy(%esp),%xmm5
        movaps nb110nf_iz(%esp),%xmm6

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

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb110nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb110nf_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        movaps %xmm0,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvsq 
        movaps %xmm4,%xmm1
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  %xmm0,%xmm3      ## xmm3=vcoul 
        mulps  nb110nf_c6(%esp),%xmm1
        mulps  nb110nf_c12(%esp),%xmm2
        movaps %xmm2,%xmm5
        subps  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addss  nb110nf_Vvdwtot(%esp),%xmm5
        addss  nb110nf_vctot(%esp),%xmm3
        movss %xmm3,nb110nf_vctot(%esp)
        movss %xmm5,nb110nf_Vvdwtot(%esp)

_nb_kernel110nf_ia32_sse.nb110nf_updateouterdata: 
        ## get n from stack
        movl nb110nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb110nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb110nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb110nf_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb110nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb110nf_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb110nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel110nf_ia32_sse.nb110nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb110nf_n(%esp)
        jmp _nb_kernel110nf_ia32_sse.nb110nf_outer
_nb_kernel110nf_ia32_sse.nb110nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb110nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel110nf_ia32_sse.nb110nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel110nf_ia32_sse.nb110nf_threadloop
_nb_kernel110nf_ia32_sse.nb110nf_end: 
        emms

        movl nb110nf_nouter(%esp),%eax
        movl nb110nf_ninner(%esp),%ebx
        movl nb110nf_outeriter(%ebp),%ecx
        movl nb110nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb110nf_salign(%esp),%eax
        addl %eax,%esp
        addl $216,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret


