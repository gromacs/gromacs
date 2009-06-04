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





.globl nb_kernel100_ia32_sse
.globl _nb_kernel100_ia32_sse
nb_kernel100_ia32_sse:  
_nb_kernel100_ia32_sse: 
.set nb100_p_nri, 8
.set nb100_iinr, 12
.set nb100_jindex, 16
.set nb100_jjnr, 20
.set nb100_shift, 24
.set nb100_shiftvec, 28
.set nb100_fshift, 32
.set nb100_gid, 36
.set nb100_pos, 40
.set nb100_faction, 44
.set nb100_charge, 48
.set nb100_p_facel, 52
.set nb100_p_krf, 56
.set nb100_p_crf, 60
.set nb100_Vc, 64
.set nb100_type, 68
.set nb100_p_ntype, 72
.set nb100_vdwparam, 76
.set nb100_Vvdw, 80
.set nb100_p_tabscale, 84
.set nb100_VFtab, 88
.set nb100_invsqrta, 92
.set nb100_dvda, 96
.set nb100_p_gbtabscale, 100
.set nb100_GBtab, 104
.set nb100_p_nthreads, 108
.set nb100_count, 112
.set nb100_mtx, 116
.set nb100_outeriter, 120
.set nb100_inneriter, 124
.set nb100_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb100_ix, 0
.set nb100_iy, 16
.set nb100_iz, 32
.set nb100_iq, 48
.set nb100_dx, 64
.set nb100_dy, 80
.set nb100_dz, 96
.set nb100_vctot, 112
.set nb100_fix, 128
.set nb100_fiy, 144
.set nb100_fiz, 160
.set nb100_half, 176
.set nb100_three, 192
.set nb100_is3, 208
.set nb100_ii3, 212
.set nb100_innerjjnr, 216
.set nb100_innerk, 220
.set nb100_n, 224
.set nb100_nn1, 228
.set nb100_nri, 232
.set nb100_facel, 236
.set nb100_nouter, 240
.set nb100_ninner, 244
.set nb100_salign, 248
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $252,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb100_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb100_p_nri(%ebp),%ecx
        movl nb100_p_facel(%ebp),%esi
        movl (%ecx),%ecx
        movl (%esi),%esi
        movl %ecx,nb100_nri(%esp)
        movl %esi,nb100_facel(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb100_nouter(%esp)
        movl %eax,nb100_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb100_half(%esp)
        movss nb100_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb100_half(%esp)
        movaps %xmm3,nb100_three(%esp)

_nb_kernel100_ia32_sse.nb100_threadloop: 
        movl  nb100_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel100_ia32_sse.nb100_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel100_ia32_sse.nb100_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb100_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb100_n(%esp)
        movl %ebx,nb100_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel100_ia32_sse.nb100_outerstart
        jmp _nb_kernel100_ia32_sse.nb100_end

_nb_kernel100_ia32_sse.nb100_outerstart: 
        ## ebx contains number of outer iterations
        addl nb100_nouter(%esp),%ebx
        movl %ebx,nb100_nouter(%esp)

_nb_kernel100_ia32_sse.nb100_outer: 
        movl  nb100_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb100_is3(%esp)      ## store is3 

        movl  nb100_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb100_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 

        movl  nb100_charge(%ebp),%edx
        movss (%edx,%ebx,4),%xmm3
        mulss nb100_facel(%esp),%xmm3
        shufps $0,%xmm3,%xmm3

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb100_pos(%ebp),%eax      ## eax = base of pos[]  

        addss (%eax,%ebx,4),%xmm0
        addss 4(%eax,%ebx,4),%xmm1
        addss 8(%eax,%ebx,4),%xmm2

        movaps %xmm3,nb100_iq(%esp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb100_ix(%esp)
        movaps %xmm1,nb100_iy(%esp)
        movaps %xmm2,nb100_iz(%esp)

        movl  %ebx,nb100_ii3(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb100_vctot(%esp)
        movaps %xmm4,nb100_fix(%esp)
        movaps %xmm4,nb100_fiy(%esp)
        movaps %xmm4,nb100_fiz(%esp)

        movl  nb100_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb100_pos(%ebp),%esi
        movl  nb100_faction(%ebp),%edi
        movl  nb100_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb100_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb100_ninner(%esp),%ecx
        movl  %ecx,nb100_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb100_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel100_ia32_sse.nb100_unroll_loop
        jmp   _nb_kernel100_ia32_sse.nb100_finish_inner
_nb_kernel100_ia32_sse.nb100_unroll_loop: 
        ## quad-unrolled innerloop here 
        movl  nb100_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 
        addl $16,nb100_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb100_charge(%ebp),%esi     ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        movaps nb100_iq(%esp),%xmm5
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000          
        movl nb100_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        mulps %xmm5,%xmm3
        leal  (%ecx,%ecx,2),%ecx     ## replace jnr with j3 
        leal  (%edx,%edx,2),%edx

        ## move four coordinates to xmm0-xmm2   

        movlps (%esi,%eax,4),%xmm4      ## x1 y1 - - 
        movlps (%esi,%ecx,4),%xmm5      ## x3 y3 - - 
        movss 8(%esi,%eax,4),%xmm2      ## z1 -  - - 
        movss 8(%esi,%ecx,4),%xmm6      ## z3 -  - - 

        movhps (%esi,%ebx,4),%xmm4      ## x1 y1 x2 y2 
        movhps (%esi,%edx,4),%xmm5      ## x3 y3 x4 y4 

        movss 8(%esi,%ebx,4),%xmm0      ## z2 - - - 
        movss 8(%esi,%edx,4),%xmm1      ## z4 - - - 

        shufps $0,%xmm0,%xmm2          ## z1 z1 z2 z2 
        shufps $0,%xmm1,%xmm6          ## z3 z3 z4 z4 

        movaps %xmm4,%xmm0              ## x1 y1 x2 y2  
        movaps %xmm4,%xmm1              ## x1 y1 x2 y2 

        shufps $136,%xmm6,%xmm2 ## constant 10001000    ;# z1 z2 z3 z4 

        shufps $136,%xmm5,%xmm0 ## constant 10001000    ;# x1 x2 x3 x4 
        shufps $221,%xmm5,%xmm1 ## constant 11011101    ;# y1 y2 y3 y4          

        movl   nb100_faction(%ebp),%edi

        ## move nb100_ix-iz to xmm4-xmm6 
        movaps nb100_ix(%esp),%xmm4
        movaps nb100_iy(%esp),%xmm5
        movaps nb100_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb100_dx(%esp)
        movaps %xmm5,nb100_dy(%esp)
        movaps %xmm6,nb100_dz(%esp)
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
        movaps nb100_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb100_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        movaps %xmm0,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvsq 

        movaps nb100_vctot(%esp),%xmm5
        mulps  %xmm0,%xmm3      ## xmm3=vcoul 
        mulps  %xmm3,%xmm4      ## xmm4=fscal 
        addps  %xmm3,%xmm5

        movaps nb100_dx(%esp),%xmm0
        movaps nb100_dy(%esp),%xmm1
        movaps nb100_dz(%esp),%xmm2

        movaps %xmm5,nb100_vctot(%esp)

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movaps nb100_fix(%esp),%xmm3
        movaps nb100_fiy(%esp),%xmm4
        movaps nb100_fiz(%esp),%xmm5
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm5
        movaps %xmm3,nb100_fix(%esp)
        movaps %xmm4,nb100_fiy(%esp)
        movaps %xmm5,nb100_fiz(%esp)
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
        subl $4,nb100_innerk(%esp)
        jl    _nb_kernel100_ia32_sse.nb100_finish_inner
        jmp   _nb_kernel100_ia32_sse.nb100_unroll_loop
_nb_kernel100_ia32_sse.nb100_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb100_innerk(%esp)
        movl  nb100_innerk(%esp),%edx
        andl  $2,%edx
        jnz   _nb_kernel100_ia32_sse.nb100_dopair
        jmp   _nb_kernel100_ia32_sse.nb100_checksingle
_nb_kernel100_ia32_sse.nb100_dopair: 
        movl nb100_charge(%ebp),%esi
        movl nb100_pos(%ebp),%edi
        movl  nb100_innerjjnr(%esp),%ecx

        movl  (%ecx),%eax
        movl  4(%ecx),%ebx
        addl $8,nb100_innerjjnr(%esp)

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ebx,4),%xmm6
        shufps $0,%xmm6,%xmm3
        shufps $8,%xmm3,%xmm3 ## constant 00001000 ;# xmm3(0,1) has the charges 

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        ## move coordinates to xmm0-xmm2 
        movlps (%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2
        movhps (%edi,%ebx,4),%xmm1
        movss 8(%edi,%ebx,4),%xmm0

        mulps  nb100_iq(%esp),%xmm3
        xorps  %xmm7,%xmm7
        movlhps %xmm7,%xmm3

        shufps $0,%xmm0,%xmm2

        movaps %xmm1,%xmm0

        shufps $136,%xmm2,%xmm2 ## constant 10001000

        shufps $136,%xmm0,%xmm0 ## constant 10001000
        shufps $221,%xmm1,%xmm1 ## constant 11011101

        movl   nb100_faction(%ebp),%edi
        ## move nb100_ix-iz to xmm4-xmm6 
        xorps   %xmm7,%xmm7

        movaps nb100_ix(%esp),%xmm4
        movaps nb100_iy(%esp),%xmm5
        movaps nb100_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb100_dx(%esp)
        movaps %xmm5,nb100_dy(%esp)
        movaps %xmm6,nb100_dz(%esp)
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
        movaps nb100_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb100_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        movaps %xmm0,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvsq 

        movaps nb100_vctot(%esp),%xmm5
        mulps  %xmm0,%xmm3      ## xmm3=vcoul 
        mulps  %xmm3,%xmm4      ## xmm4=fscal 
        addps  %xmm3,%xmm5

        movaps nb100_dx(%esp),%xmm0
        movaps nb100_dy(%esp),%xmm1
        movaps nb100_dz(%esp),%xmm2

        movaps %xmm5,nb100_vctot(%esp)

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movaps nb100_fix(%esp),%xmm3
        movaps nb100_fiy(%esp),%xmm4
        movaps nb100_fiz(%esp),%xmm5
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm5
        movaps %xmm3,nb100_fix(%esp)
        movaps %xmm4,nb100_fiy(%esp)
        movaps %xmm5,nb100_fiz(%esp)
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
_nb_kernel100_ia32_sse.nb100_checksingle:       
        movl  nb100_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel100_ia32_sse.nb100_dosingle
        jmp    _nb_kernel100_ia32_sse.nb100_updateouterdata
_nb_kernel100_ia32_sse.nb100_dosingle:  
        movl nb100_charge(%ebp),%esi
        movl nb100_pos(%ebp),%edi
        movl  nb100_innerjjnr(%esp),%ecx
        movl  (%ecx),%eax
        movss (%esi,%eax,4),%xmm3       ## xmm3(0) has the charge       

        leal  (%eax,%eax,2),%eax

        ## move coordinates to xmm0-xmm2 
        movss (%edi,%eax,4),%xmm0
        movss 4(%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2

        mulps  nb100_iq(%esp),%xmm3

        xorps   %xmm7,%xmm7

        movaps nb100_ix(%esp),%xmm4
        movaps nb100_iy(%esp),%xmm5
        movaps nb100_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb100_dx(%esp)
        movaps %xmm5,nb100_dy(%esp)
        movaps %xmm6,nb100_dz(%esp)
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
        movaps nb100_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb100_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        movaps %xmm0,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvsq 
        movl   nb100_faction(%ebp),%edi
        movaps nb100_vctot(%esp),%xmm5
        mulps  %xmm0,%xmm3      ## xmm3=vcoul 
        mulps  %xmm3,%xmm4      ## xmm4=fscal 
        addss  %xmm3,%xmm5

        movaps nb100_dx(%esp),%xmm0
        movaps nb100_dy(%esp),%xmm1
        movaps nb100_dz(%esp),%xmm2

        movaps %xmm5,nb100_vctot(%esp)

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movaps nb100_fix(%esp),%xmm3
        movaps nb100_fiy(%esp),%xmm4
        movaps nb100_fiz(%esp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movaps %xmm3,nb100_fix(%esp)
        movaps %xmm4,nb100_fiy(%esp)
        movaps %xmm5,nb100_fiz(%esp)
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
_nb_kernel100_ia32_sse.nb100_updateouterdata: 
        movl  nb100_ii3(%esp),%ecx
        movl  nb100_faction(%ebp),%edi
        movl  nb100_fshift(%ebp),%esi
        movl  nb100_is3(%esp),%edx

        ## accumulate i forces in xmm0, xmm1, xmm2 
        movaps nb100_fix(%esp),%xmm0
        movaps nb100_fiy(%esp),%xmm1
        movaps nb100_fiz(%esp),%xmm2

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
        movl nb100_n(%esp),%esi
        ## get group index for i particle 
        movl  nb100_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb100_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb100_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb100_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel100_ia32_sse.nb100_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb100_n(%esp)
        jmp _nb_kernel100_ia32_sse.nb100_outer
_nb_kernel100_ia32_sse.nb100_outerend: 
        ## check if more outer neighborlists remain
        movl  nb100_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel100_ia32_sse.nb100_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel100_ia32_sse.nb100_threadloop
_nb_kernel100_ia32_sse.nb100_end: 
        emms

        movl nb100_nouter(%esp),%eax
        movl nb100_ninner(%esp),%ebx
        movl nb100_outeriter(%ebp),%ecx
        movl nb100_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb100_salign(%esp),%eax
        addl %eax,%esp
        addl $252,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret




.globl nb_kernel100nf_ia32_sse
.globl _nb_kernel100nf_ia32_sse
nb_kernel100nf_ia32_sse:        
_nb_kernel100nf_ia32_sse:       
.set nb100nf_p_nri, 8
.set nb100nf_iinr, 12
.set nb100nf_jindex, 16
.set nb100nf_jjnr, 20
.set nb100nf_shift, 24
.set nb100nf_shiftvec, 28
.set nb100nf_fshift, 32
.set nb100nf_gid, 36
.set nb100nf_pos, 40
.set nb100nf_faction, 44
.set nb100nf_charge, 48
.set nb100nf_p_facel, 52
.set nb100nf_p_krf, 56
.set nb100nf_p_crf, 60
.set nb100nf_Vc, 64
.set nb100nf_type, 68
.set nb100nf_p_ntype, 72
.set nb100nf_vdwparam, 76
.set nb100nf_Vvdw, 80
.set nb100nf_p_tabscale, 84
.set nb100nf_VFtab, 88
.set nb100nf_invsqrta, 92
.set nb100nf_dvda, 96
.set nb100nf_p_gbtabscale, 100
.set nb100nf_GBtab, 104
.set nb100nf_p_nthreads, 108
.set nb100nf_count, 112
.set nb100nf_mtx, 116
.set nb100nf_outeriter, 120
.set nb100nf_inneriter, 124
.set nb100nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb100nf_ix, 0
.set nb100nf_iy, 16
.set nb100nf_iz, 32
.set nb100nf_iq, 48
.set nb100nf_vctot, 64
.set nb100nf_half, 80
.set nb100nf_three, 96
.set nb100nf_is3, 112
.set nb100nf_ii3, 116
.set nb100nf_innerjjnr, 120
.set nb100nf_innerk, 124
.set nb100nf_n, 128
.set nb100nf_nn1, 132
.set nb100nf_nri, 136
.set nb100nf_facel, 140
.set nb100nf_nouter, 144
.set nb100nf_ninner, 148
.set nb100nf_salign, 152
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $156,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb100nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb100nf_p_nri(%ebp),%ecx
        movl nb100nf_p_facel(%ebp),%esi
        movl (%ecx),%ecx
        movl (%esi),%esi
        movl %ecx,nb100nf_nri(%esp)
        movl %esi,nb100nf_facel(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb100nf_nouter(%esp)
        movl %eax,nb100nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb100nf_half(%esp)
        movss nb100nf_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb100nf_half(%esp)
        movaps %xmm3,nb100nf_three(%esp)

_nb_kernel100nf_ia32_sse.nb100nf_threadloop: 
        movl  nb100nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel100nf_ia32_sse.nb100nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel100nf_ia32_sse.nb100nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb100nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb100nf_n(%esp)
        movl %ebx,nb100nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel100nf_ia32_sse.nb100nf_outerstart
        jmp _nb_kernel100nf_ia32_sse.nb100nf_end

_nb_kernel100nf_ia32_sse.nb100nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb100nf_nouter(%esp),%ebx
        movl %ebx,nb100nf_nouter(%esp)

_nb_kernel100nf_ia32_sse.nb100nf_outer: 
        movl  nb100nf_shift(%ebp),%eax          ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb100nf_is3(%esp)            ## store is3 

        movl  nb100nf_shiftvec(%ebp),%eax       ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb100nf_iinr(%ebp),%ecx           ## ecx = pointer into iinr[]    
        movl  (%ecx,%esi,4),%ebx                ## ebx =ii 

        movl  nb100nf_charge(%ebp),%edx
        movss (%edx,%ebx,4),%xmm3
        mulss nb100nf_facel(%esp),%xmm3
        shufps $0,%xmm3,%xmm3

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb100nf_pos(%ebp),%eax      ## eax = base of pos[]  

        addss (%eax,%ebx,4),%xmm0
        addss 4(%eax,%ebx,4),%xmm1
        addss 8(%eax,%ebx,4),%xmm2

        movaps %xmm3,nb100nf_iq(%esp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb100nf_ix(%esp)
        movaps %xmm1,nb100nf_iy(%esp)
        movaps %xmm2,nb100nf_iz(%esp)

        movl  %ebx,nb100nf_ii3(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb100nf_vctot(%esp)

        movl  nb100nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb100nf_pos(%ebp),%esi

        movl  nb100nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb100nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb100nf_ninner(%esp),%ecx
        movl  %ecx,nb100nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb100nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel100nf_ia32_sse.nb100nf_unroll_loop
        jmp   _nb_kernel100nf_ia32_sse.nb100nf_finish_inner
_nb_kernel100nf_ia32_sse.nb100nf_unroll_loop: 
        ## quad-unrolled innerloop here 
        movl  nb100nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 
        addl $16,nb100nf_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb100nf_charge(%ebp),%esi     ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        movaps nb100nf_iq(%esp),%xmm5
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000          
        movl nb100nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        mulps %xmm5,%xmm3
        leal  (%ecx,%ecx,2),%ecx     ## replace jnr with j3 
        leal  (%edx,%edx,2),%edx

        ## move four coordinates to xmm0-xmm2   

        movlps (%esi,%eax,4),%xmm4      ## x1 y1 - - 
        movlps (%esi,%ecx,4),%xmm5      ## x3 y3 - - 
        movss 8(%esi,%eax,4),%xmm2      ## z1 -  - - 
        movss 8(%esi,%ecx,4),%xmm6      ## z3 -  - - 

        movhps (%esi,%ebx,4),%xmm4      ## x1 y1 x2 y2 
        movhps (%esi,%edx,4),%xmm5      ## x3 y3 x4 y4 

        movss 8(%esi,%ebx,4),%xmm0      ## z2 - - - 
        movss 8(%esi,%edx,4),%xmm1      ## z4 - - - 

        shufps $0,%xmm0,%xmm2          ## z1 z1 z2 z2 
        shufps $0,%xmm1,%xmm6          ## z3 z3 z4 z4 

        movaps %xmm4,%xmm0              ## x1 y1 x2 y2  
        movaps %xmm4,%xmm1              ## x1 y1 x2 y2 

        shufps $136,%xmm6,%xmm2 ## constant 10001000    ;# z1 z2 z3 z4 

        shufps $136,%xmm5,%xmm0 ## constant 10001000    ;# x1 x2 x3 x4 
        shufps $221,%xmm5,%xmm1 ## constant 11011101    ;# y1 y2 y3 y4          

        ## move nb100nf_ix-iz to xmm4-xmm6 
        movaps nb100nf_ix(%esp),%xmm4
        movaps nb100nf_iy(%esp),%xmm5
        movaps nb100nf_iz(%esp),%xmm6

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
        movaps nb100nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb100nf_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 

        movaps nb100nf_vctot(%esp),%xmm5
        mulps  %xmm0,%xmm3      ## xmm3=vcoul 
        addps  %xmm3,%xmm5
        movaps %xmm5,nb100nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $4,nb100nf_innerk(%esp)
        jl    _nb_kernel100nf_ia32_sse.nb100nf_finish_inner
        jmp   _nb_kernel100nf_ia32_sse.nb100nf_unroll_loop
_nb_kernel100nf_ia32_sse.nb100nf_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb100nf_innerk(%esp)
        movl  nb100nf_innerk(%esp),%edx
        andl  $2,%edx
        jnz   _nb_kernel100nf_ia32_sse.nb100nf_dopair
        jmp   _nb_kernel100nf_ia32_sse.nb100nf_checksingle
_nb_kernel100nf_ia32_sse.nb100nf_dopair: 
        movl nb100nf_charge(%ebp),%esi
        movl nb100nf_pos(%ebp),%edi
        movl  nb100nf_innerjjnr(%esp),%ecx

        movl  (%ecx),%eax
        movl  4(%ecx),%ebx
        addl $8,nb100nf_innerjjnr(%esp)

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ebx,4),%xmm6
        shufps $0,%xmm6,%xmm3
        shufps $8,%xmm3,%xmm3 ## constant 00001000 ;# xmm3(0,1) has the charges 

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        ## move coordinates to xmm0-xmm2 
        movlps (%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2
        movhps (%edi,%ebx,4),%xmm1
        movss 8(%edi,%ebx,4),%xmm0

        mulps  nb100nf_iq(%esp),%xmm3
        xorps  %xmm7,%xmm7
        movlhps %xmm7,%xmm3

        shufps $0,%xmm0,%xmm2

        movaps %xmm1,%xmm0

        shufps $136,%xmm2,%xmm2 ## constant 10001000

        shufps $136,%xmm0,%xmm0 ## constant 10001000
        shufps $221,%xmm1,%xmm1 ## constant 11011101

        ## move nb100nf_ix-iz to xmm4-xmm6 
        xorps   %xmm7,%xmm7

        movaps nb100nf_ix(%esp),%xmm4
        movaps nb100nf_iy(%esp),%xmm5
        movaps nb100nf_iz(%esp),%xmm6

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
        movaps nb100nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb100nf_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        movaps %xmm0,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvsq 

        movaps nb100nf_vctot(%esp),%xmm5
        mulps  %xmm0,%xmm3      ## xmm3=vcoul 
        addps  %xmm3,%xmm5
        movaps %xmm5,nb100nf_vctot(%esp)
_nb_kernel100nf_ia32_sse.nb100nf_checksingle:   
        movl  nb100nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel100nf_ia32_sse.nb100nf_dosingle
        jmp    _nb_kernel100nf_ia32_sse.nb100nf_updateouterdata
_nb_kernel100nf_ia32_sse.nb100nf_dosingle: 
        movl nb100nf_charge(%ebp),%esi
        movl nb100nf_pos(%ebp),%edi
        movl  nb100nf_innerjjnr(%esp),%ecx
        movl  (%ecx),%eax
        movss (%esi,%eax,4),%xmm3       ## xmm3(0) has the charge       

        leal  (%eax,%eax,2),%eax

        ## move coordinates to xmm0-xmm2 
        movss (%edi,%eax,4),%xmm0
        movss 4(%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2

        mulps  nb100nf_iq(%esp),%xmm3

        xorps   %xmm7,%xmm7

        movaps nb100nf_ix(%esp),%xmm4
        movaps nb100nf_iy(%esp),%xmm5
        movaps nb100nf_iz(%esp),%xmm6

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
        movaps nb100nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb100nf_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        movaps %xmm0,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvsq 
        movaps nb100nf_vctot(%esp),%xmm5
        mulps  %xmm0,%xmm3      ## xmm3=vcoul 
        addss  %xmm3,%xmm5
        movaps %xmm5,nb100nf_vctot(%esp)

_nb_kernel100nf_ia32_sse.nb100nf_updateouterdata: 
        ## get n from stack
        movl nb100nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb100nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb100nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb100nf_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb100nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel100nf_ia32_sse.nb100nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb100nf_n(%esp)
        jmp _nb_kernel100nf_ia32_sse.nb100nf_outer
_nb_kernel100nf_ia32_sse.nb100nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb100nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel100nf_ia32_sse.nb100nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel100nf_ia32_sse.nb100nf_threadloop
_nb_kernel100nf_ia32_sse.nb100nf_end: 
        emms

        movl nb100nf_nouter(%esp),%eax
        movl nb100nf_ninner(%esp),%ebx
        movl nb100nf_outeriter(%ebp),%ecx
        movl nb100nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb100nf_salign(%esp),%eax
        addl %eax,%esp
        addl $156,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret

