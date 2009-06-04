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



.globl nb_kernel400_ia32_sse
.globl _nb_kernel400_ia32_sse
nb_kernel400_ia32_sse:  
_nb_kernel400_ia32_sse: 
.set nb400_p_nri, 8
.set nb400_iinr, 12
.set nb400_jindex, 16
.set nb400_jjnr, 20
.set nb400_shift, 24
.set nb400_shiftvec, 28
.set nb400_fshift, 32
.set nb400_gid, 36
.set nb400_pos, 40
.set nb400_faction, 44
.set nb400_charge, 48
.set nb400_p_facel, 52
.set nb400_argkrf, 56
.set nb400_argcrf, 60
.set nb400_Vc, 64
.set nb400_type, 68
.set nb400_p_ntype, 72
.set nb400_vdwparam, 76
.set nb400_Vvdw, 80
.set nb400_p_tabscale, 84
.set nb400_VFtab, 88
.set nb400_invsqrta, 92
.set nb400_dvda, 96
.set nb400_p_gbtabscale, 100
.set nb400_GBtab, 104
.set nb400_p_nthreads, 108
.set nb400_count, 112
.set nb400_mtx, 116
.set nb400_outeriter, 120
.set nb400_inneriter, 124
.set nb400_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb400_ix, 0
.set nb400_iy, 16
.set nb400_iz, 32
.set nb400_iq, 48
.set nb400_dx, 64
.set nb400_dy, 80
.set nb400_dz, 96
.set nb400_two, 112
.set nb400_gbtsc, 128
.set nb400_qq, 144
.set nb400_r, 160
.set nb400_vctot, 176
.set nb400_fix, 192
.set nb400_fiy, 208
.set nb400_fiz, 224
.set nb400_half, 240
.set nb400_three, 256
.set nb400_isai, 272
.set nb400_isaprod, 288
.set nb400_dvdasum, 304
.set nb400_gbscale, 320
.set nb400_is3, 336
.set nb400_ii3, 340
.set nb400_ii, 344
.set nb400_innerjjnr, 348
.set nb400_innerk, 352
.set nb400_n, 356
.set nb400_nn1, 360
.set nb400_jnra, 364
.set nb400_jnrb, 368
.set nb400_jnrc, 372
.set nb400_jnrd, 376
.set nb400_nri, 380
.set nb400_facel, 384
.set nb400_nouter, 388
.set nb400_ninner, 392
.set nb400_salign, 396
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
        movl %eax,nb400_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb400_p_nri(%ebp),%ecx
        movl nb400_p_facel(%ebp),%esi
        movl (%ecx),%ecx
        movl (%esi),%esi
        movl %ecx,nb400_nri(%esp)
        movl %esi,nb400_facel(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb400_nouter(%esp)
        movl %eax,nb400_ninner(%esp)


        movl nb400_p_gbtabscale(%ebp),%eax
        movss (%eax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb400_gbtsc(%esp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb400_half(%esp)
        movss nb400_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb400_half(%esp)
        movaps %xmm2,nb400_two(%esp)
        movaps %xmm3,nb400_three(%esp)

_nb_kernel400_ia32_sse.nb400_threadloop: 
        movl  nb400_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel400_ia32_sse.nb400_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel400_ia32_sse.nb400_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb400_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb400_n(%esp)
        movl %ebx,nb400_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel400_ia32_sse.nb400_outerstart
        jmp _nb_kernel400_ia32_sse.nb400_end

_nb_kernel400_ia32_sse.nb400_outerstart: 
        ## ebx contains number of outer iterations
        addl nb400_nouter(%esp),%ebx
        movl %ebx,nb400_nouter(%esp)

_nb_kernel400_ia32_sse.nb400_outer: 
        movl  nb400_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb400_is3(%esp)      ## store is3 

        movl  nb400_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb400_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 
        movl  %ebx,nb400_ii(%esp)

        movl  nb400_charge(%ebp),%edx
        movss (%edx,%ebx,4),%xmm3
        mulss nb400_facel(%esp),%xmm3
        shufps $0,%xmm3,%xmm3


        movl  nb400_invsqrta(%ebp),%edx         ## load invsqrta[ii]
        movss (%edx,%ebx,4),%xmm4
        shufps $0,%xmm4,%xmm4

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb400_pos(%ebp),%eax      ## eax = base of pos[]  

        addss (%eax,%ebx,4),%xmm0
        addss 4(%eax,%ebx,4),%xmm1
        addss 8(%eax,%ebx,4),%xmm2

        movaps %xmm3,nb400_iq(%esp)
        movaps %xmm4,nb400_isai(%esp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb400_ix(%esp)
        movaps %xmm1,nb400_iy(%esp)
        movaps %xmm2,nb400_iz(%esp)

        movl  %ebx,nb400_ii3(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb400_vctot(%esp)
        movaps %xmm4,nb400_dvdasum(%esp)
        movaps %xmm4,nb400_fix(%esp)
        movaps %xmm4,nb400_fiy(%esp)
        movaps %xmm4,nb400_fiz(%esp)

        movl  nb400_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb400_pos(%ebp),%esi
        movl  nb400_faction(%ebp),%edi
        movl  nb400_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb400_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb400_ninner(%esp),%ecx
        movl  %ecx,nb400_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb400_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel400_ia32_sse.nb400_unroll_loop
        jmp   _nb_kernel400_ia32_sse.nb400_finish_inner
_nb_kernel400_ia32_sse.nb400_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb400_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 
        addl $16,nb400_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        ## load isaj
        movl nb400_invsqrta(%ebp),%esi
        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7
        movaps nb400_isai(%esp),%xmm2
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all isaj in xmm3 
        mulps  %xmm3,%xmm2

        movaps %xmm2,nb400_isaprod(%esp)
        movaps %xmm2,%xmm1
        mulps nb400_gbtsc(%esp),%xmm1
        movaps %xmm1,nb400_gbscale(%esp)

        movl nb400_charge(%ebp),%esi     ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        mulps nb400_iq(%esp),%xmm2
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        mulps  %xmm2,%xmm3
        movaps %xmm3,nb400_qq(%esp)


        movl nb400_pos(%ebp),%esi        ## base of pos[] 

        movl %eax,nb400_jnra(%esp)
        movl %ebx,nb400_jnrb(%esp)
        movl %ecx,nb400_jnrc(%esp)
        movl %edx,nb400_jnrd(%esp)

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx
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
        movaps nb400_ix(%esp),%xmm4
        movaps nb400_iy(%esp),%xmm5
        movaps nb400_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb400_dx(%esp)
        movaps %xmm5,nb400_dy(%esp)
        movaps %xmm6,nb400_dz(%esp)
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
        movaps nb400_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb400_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r
        movaps %xmm4,nb400_r(%esp)
        mulps nb400_gbscale(%esp),%xmm4

        movhlps %xmm4,%xmm5
        cvttps2pi %xmm4,%mm6
        cvttps2pi %xmm5,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        cvtpi2ps %mm7,%xmm5
        movlhps %xmm5,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $2,%mm6
        pslld $2,%mm7

        movd %eax,%mm0
        movd %ebx,%mm1
        movd %ecx,%mm2
        movd %edx,%mm3

        movl nb400_GBtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        ## load coulomb table
        movaps (%esi,%eax,4),%xmm4
        movaps (%esi,%ebx,4),%xmm5
        movaps (%esi,%ecx,4),%xmm6
        movaps (%esi,%edx,4),%xmm7
        ## transpose, using xmm3 for scratch
        movaps %xmm6,%xmm3
        shufps $0xEE,%xmm7,%xmm3
        shufps $0x44,%xmm7,%xmm6
        movaps %xmm4,%xmm7
        shufps $0xEE,%xmm5,%xmm7
        shufps $0x44,%xmm5,%xmm4
        movaps %xmm4,%xmm5
        shufps $0xDD,%xmm6,%xmm5
        shufps $0x88,%xmm6,%xmm4
        movaps %xmm7,%xmm6
        shufps $0x88,%xmm3,%xmm6
        shufps $0xDD,%xmm3,%xmm7
        ## coulomb table ready, in xmm4-xmm7    

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  nb400_two(%esp),%xmm7    ## two*Heps2 
        movaps nb400_qq(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC

        ## get jnr from stack
        movl nb400_jnra(%esp),%eax
        movl nb400_jnrb(%esp),%ebx
        movl nb400_jnrc(%esp),%ecx
        movl nb400_jnrd(%esp),%edx

        movl nb400_dvda(%ebp),%esi

        ## Calculate dVda
        xorps  %xmm7,%xmm7
        mulps nb400_gbscale(%esp),%xmm3
        movaps %xmm3,%xmm6
        mulps  nb400_r(%esp),%xmm6
        addps  %xmm5,%xmm6
        addps  nb400_vctot(%esp),%xmm5
        movaps %xmm5,nb400_vctot(%esp)

        ## xmm6=(vcoul+fijC*r)
        subps  %xmm6,%xmm7
        movaps %xmm7,%xmm6

        ## update dvdasum
        addps  nb400_dvdasum(%esp),%xmm7
        movaps %xmm7,nb400_dvdasum(%esp)

        ## update j atoms dvdaj
        movhlps %xmm6,%xmm7
        movaps  %xmm6,%xmm5
        movaps  %xmm7,%xmm4
        shufps $0x1,%xmm5,%xmm5
        shufps $0x1,%xmm4,%xmm4
        ## xmm6=dvdaj1 xmm5=dvdaj2 xmm7=dvdaj3 xmm4=dvdaj4
        addss  (%esi,%eax,4),%xmm6
        addss  (%esi,%ebx,4),%xmm5
        addss  (%esi,%ecx,4),%xmm7
        addss  (%esi,%edx,4),%xmm4
        movss  %xmm6,(%esi,%eax,4)
        movss  %xmm5,(%esi,%ebx,4)
        movss  %xmm7,(%esi,%ecx,4)
        movss  %xmm4,(%esi,%edx,4)

        xorps  %xmm4,%xmm4
        mulps %xmm0,%xmm3
        subps  %xmm3,%xmm4

        movaps nb400_dx(%esp),%xmm0
        movaps nb400_dy(%esp),%xmm1
        movaps nb400_dz(%esp),%xmm2

        movd %mm0,%eax
        movd %mm1,%ebx
        movd %mm2,%ecx
        movd %mm3,%edx

        movl   nb400_faction(%ebp),%edi
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movaps nb400_fix(%esp),%xmm3
        movaps nb400_fiy(%esp),%xmm4
        movaps nb400_fiz(%esp),%xmm5
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm5
        movaps %xmm3,nb400_fix(%esp)
        movaps %xmm4,nb400_fiy(%esp)
        movaps %xmm5,nb400_fiz(%esp)
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
        subl $4,nb400_innerk(%esp)
        jl    _nb_kernel400_ia32_sse.nb400_finish_inner
        jmp   _nb_kernel400_ia32_sse.nb400_unroll_loop
_nb_kernel400_ia32_sse.nb400_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb400_innerk(%esp)
        movl  nb400_innerk(%esp),%edx
        andl  $2,%edx
        jnz   _nb_kernel400_ia32_sse.nb400_dopair
        jmp   _nb_kernel400_ia32_sse.nb400_checksingle
_nb_kernel400_ia32_sse.nb400_dopair: 
        movl  nb400_innerjjnr(%esp),%ecx

        movl  (%ecx),%eax
        movl  4(%ecx),%ebx
        addl $8,nb400_innerjjnr(%esp)

        xorps %xmm2,%xmm2
        movaps %xmm2,%xmm6

        ## load isaj
        movl nb400_invsqrta(%ebp),%esi
        movss (%esi,%eax,4),%xmm2
        movss (%esi,%ebx,4),%xmm3
        unpcklps %xmm3,%xmm2    ## isaj in xmm2(0,1)
        mulps  nb400_isai(%esp),%xmm2
        movaps %xmm2,nb400_isaprod(%esp)
        movaps %xmm2,%xmm1
        mulps nb400_gbtsc(%esp),%xmm1
        movaps %xmm1,nb400_gbscale(%esp)

        movl nb400_charge(%ebp),%esi     ## base of charge[]    
        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ebx,4),%xmm6
        unpcklps %xmm6,%xmm3 ## constant 00001000 ;# xmm3(0,1) has the charges 

        mulps  nb400_iq(%esp),%xmm2
        mulps  %xmm2,%xmm3
        movaps %xmm3,nb400_qq(%esp)

        movl nb400_pos(%ebp),%edi

        movd  %eax,%mm0         ## copy jnr to mm0/mm1
        movd  %ebx,%mm1

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        ## move coordinates to xmm0-xmm2 
        movlps (%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2
        movhps (%edi,%ebx,4),%xmm1
        movss 8(%edi,%ebx,4),%xmm0

        movlhps %xmm7,%xmm3

        shufps $0,%xmm0,%xmm2

        movaps %xmm1,%xmm0

        shufps $136,%xmm2,%xmm2 ## constant 10001000

        shufps $136,%xmm0,%xmm0 ## constant 10001000
        shufps $221,%xmm1,%xmm1 ## constant 11011101

        movl   nb400_faction(%ebp),%edi
        ## move ix-iz to xmm4-xmm6 
        xorps   %xmm7,%xmm7

        movaps nb400_ix(%esp),%xmm4
        movaps nb400_iy(%esp),%xmm5
        movaps nb400_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb400_dx(%esp)
        movaps %xmm5,nb400_dy(%esp)
        movaps %xmm6,nb400_dz(%esp)
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
        movaps nb400_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb400_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r 
        movaps %xmm4,nb400_r(%esp)
        mulps nb400_gbscale(%esp),%xmm4

        cvttps2pi %xmm4,%mm6    ## mm6 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6

        movl nb400_GBtab(%ebp),%esi
        movd %mm6,%ecx
        psrlq $32,%mm6
        movd %mm6,%edx

        ## load coulomb table
        movaps (%esi,%ecx,4),%xmm4
        movaps (%esi,%edx,4),%xmm7
        ## transpose, using xmm3 for scratch
        movaps %xmm4,%xmm6
        unpcklps %xmm7,%xmm4    ## Y1 Y2 F1 F2 
        unpckhps %xmm7,%xmm6    ## G1 G2 H1 H2
        movhlps  %xmm4,%xmm5    ## F1 F2 
        movhlps  %xmm6,%xmm7    ## H1 H2
        ## coulomb table ready, in xmm4-xmm7    

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  nb400_two(%esp),%xmm7    ## two*Heps2 
        movaps nb400_qq(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC

        ## get jnr from mm0/mm1
        movd %mm0,%ecx
        movd %mm1,%edx

        movl nb400_dvda(%ebp),%esi

        ## Calculate dVda
        xorps  %xmm7,%xmm7
        mulps nb400_gbscale(%esp),%xmm3
        movaps %xmm3,%xmm6
        mulps  nb400_r(%esp),%xmm6
        addps  %xmm5,%xmm6
        addps  nb400_vctot(%esp),%xmm5
        movaps %xmm5,nb400_vctot(%esp)

        ## xmm6=(vcoul+fijC*r)
        subps  %xmm6,%xmm7
        movaps %xmm7,%xmm6

        ## update dvdasum
        addps  nb400_dvdasum(%esp),%xmm7
        movaps %xmm7,nb400_dvdasum(%esp)

        ## update j atoms dvdaj
        movaps %xmm6,%xmm7
        shufps $0x1,%xmm7,%xmm7
        addss  (%esi,%ecx,4),%xmm6
        addss  (%esi,%edx,4),%xmm7
        movss  %xmm6,(%esi,%ecx,4)
        movss  %xmm7,(%esi,%edx,4)

        xorps  %xmm4,%xmm4
        mulps %xmm0,%xmm3
        subps  %xmm3,%xmm4

        movaps nb400_dx(%esp),%xmm0
        movaps nb400_dy(%esp),%xmm1
        movaps nb400_dz(%esp),%xmm2

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movaps nb400_fix(%esp),%xmm3
        movaps nb400_fiy(%esp),%xmm4
        movaps nb400_fiz(%esp),%xmm5
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm5
        movaps %xmm3,nb400_fix(%esp)
        movaps %xmm4,nb400_fiy(%esp)
        movaps %xmm5,nb400_fiz(%esp)
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

_nb_kernel400_ia32_sse.nb400_checksingle:       
        movl  nb400_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel400_ia32_sse.nb400_dosingle
        jmp    _nb_kernel400_ia32_sse.nb400_updateouterdata
_nb_kernel400_ia32_sse.nb400_dosingle: 
        movl nb400_charge(%ebp),%esi
        movl nb400_invsqrta(%ebp),%edx
        movl nb400_pos(%ebp),%edi
        movl  nb400_innerjjnr(%esp),%ecx
        movl  (%ecx),%eax
        xorps  %xmm2,%xmm2
        movaps %xmm2,%xmm6
        movss (%edx,%eax,4),%xmm2       ## isaj
        mulss nb400_isai(%esp),%xmm2
        movss %xmm2,nb400_isaprod(%esp)
        movss %xmm2,%xmm1
        mulss nb400_gbtsc(%esp),%xmm1
        movss %xmm1,nb400_gbscale(%esp)

        mulss  nb400_iq(%esp),%xmm2
        movss (%esi,%eax,4),%xmm6       ## xmm6(0) has the charge       
        mulss  %xmm2,%xmm6
        movss %xmm6,nb400_qq(%esp)

        movd  %eax,%mm0
        leal  (%eax,%eax,2),%eax

        ## move coordinates to xmm0-xmm2 
        movss (%edi,%eax,4),%xmm0
        movss 4(%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2

        movss nb400_ix(%esp),%xmm4
        movss nb400_iy(%esp),%xmm5
        movss nb400_iz(%esp),%xmm6

        ## calc dr 
        subss %xmm0,%xmm4
        subss %xmm1,%xmm5
        subss %xmm2,%xmm6

        ## store dr 
        movss %xmm4,nb400_dx(%esp)
        movss %xmm5,nb400_dy(%esp)
        movss %xmm6,nb400_dz(%esp)
        ## square it 
        mulss %xmm4,%xmm4
        mulss %xmm5,%xmm5
        mulss %xmm6,%xmm6
        addss %xmm5,%xmm4
        addss %xmm6,%xmm4
        ## rsq in xmm4 

        rsqrtss %xmm4,%xmm5
        ## lookup seed in xmm5 
        movss %xmm5,%xmm2
        mulss %xmm5,%xmm5
        movss nb400_three(%esp),%xmm1
        mulss %xmm4,%xmm5       ## rsq*lu*lu                    
        movss nb400_half(%esp),%xmm0
        subss %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulss %xmm2,%xmm1
        mulss %xmm1,%xmm0       ## xmm0=rinv 

        mulss %xmm0,%xmm4       ## xmm4=r 
        movss %xmm4,nb400_r(%esp)
        mulss nb400_gbscale(%esp),%xmm4

        cvttss2si %xmm4,%ebx    ## mm6 contain lu indices 
        cvtsi2ss %ebx,%xmm6
        subss %xmm6,%xmm4
        movss %xmm4,%xmm1       ## xmm1=eps 
        movss %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%ebx

        movl nb400_GBtab(%ebp),%esi

        movaps (%esi,%ebx,4),%xmm4
        movhlps %xmm4,%xmm6
        movaps %xmm4,%xmm5
        movaps %xmm6,%xmm7
        shufps $1,%xmm5,%xmm5
        shufps $1,%xmm7,%xmm7
        ## table ready in xmm4-xmm7 

        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp      
        mulss  nb400_two(%esp),%xmm7    ## two*Heps2 
        movss nb400_qq(%esp),%xmm3
        addss  %xmm6,%xmm7
        addss  %xmm5,%xmm7 ## xmm7=FF 
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 
        mulss  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulss  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 

        movd %mm0,%ebx
        movl nb400_dvda(%ebp),%esi

        ## Calculate dVda
        xorps %xmm7,%xmm7
        mulss nb400_gbscale(%esp),%xmm3
        movaps %xmm3,%xmm6
        mulss  nb400_r(%esp),%xmm6
        addss  %xmm5,%xmm6
        addss  nb400_vctot(%esp),%xmm5
        movss %xmm5,nb400_vctot(%esp)

        ## xmm6=(vcoul+fijC*r)
        subps  %xmm6,%xmm7
        movaps %xmm7,%xmm6

        ## update dvdasum
        addps  nb400_dvdasum(%esp),%xmm7
        movaps %xmm7,nb400_dvdasum(%esp)

        ## update j atoms dvdaj
        addss  (%esi,%ebx,4),%xmm6
        movss  %xmm6,(%esi,%ebx,4)

        xorps  %xmm4,%xmm4
        mulss %xmm0,%xmm3
        subss  %xmm3,%xmm4

        movl   nb400_faction(%ebp),%edi

        movss nb400_dx(%esp),%xmm0
        movss nb400_dy(%esp),%xmm1
        movss nb400_dz(%esp),%xmm2

        mulss  %xmm4,%xmm0
        mulss  %xmm4,%xmm1
        mulss  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movss nb400_fix(%esp),%xmm3
        movss nb400_fiy(%esp),%xmm4
        movss nb400_fiz(%esp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss %xmm3,nb400_fix(%esp)
        movss %xmm4,nb400_fiy(%esp)
        movss %xmm5,nb400_fiz(%esp)
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
_nb_kernel400_ia32_sse.nb400_updateouterdata: 
        movl  nb400_ii3(%esp),%ecx
        movl  nb400_faction(%ebp),%edi
        movl  nb400_fshift(%ebp),%esi
        movl  nb400_is3(%esp),%edx

        ## accumulate i forces in xmm0, xmm1, xmm2 
        movaps nb400_fix(%esp),%xmm0
        movaps nb400_fiy(%esp),%xmm1
        movaps nb400_fiz(%esp),%xmm2

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
        movl nb400_n(%esp),%esi
        ## get group index for i particle 
        movl  nb400_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb400_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb400_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate dVda and update it 
        movaps nb400_dvdasum(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        movl nb400_ii(%esp),%edx
        movl nb400_dvda(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb400_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel400_ia32_sse.nb400_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb400_n(%esp)
        jmp _nb_kernel400_ia32_sse.nb400_outer
_nb_kernel400_ia32_sse.nb400_outerend: 
        ## check if more outer neighborlists remain
        movl  nb400_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel400_ia32_sse.nb400_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel400_ia32_sse.nb400_threadloop
_nb_kernel400_ia32_sse.nb400_end: 
        emms

        movl nb400_nouter(%esp),%eax
        movl nb400_ninner(%esp),%ebx
        movl nb400_outeriter(%ebp),%ecx
        movl nb400_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb400_salign(%esp),%eax
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





.globl nb_kernel400nf_ia32_sse
.globl _nb_kernel400nf_ia32_sse
nb_kernel400nf_ia32_sse:        
_nb_kernel400nf_ia32_sse:       
.set nb400nf_p_nri, 8
.set nb400nf_iinr, 12
.set nb400nf_jindex, 16
.set nb400nf_jjnr, 20
.set nb400nf_shift, 24
.set nb400nf_shiftvec, 28
.set nb400nf_fshift, 32
.set nb400nf_gid, 36
.set nb400nf_pos, 40
.set nb400nf_faction, 44
.set nb400nf_charge, 48
.set nb400nf_p_facel, 52
.set nb400nf_argkrf, 56
.set nb400nf_argcrf, 60
.set nb400nf_Vc, 64
.set nb400nf_type, 68
.set nb400nf_p_ntype, 72
.set nb400nf_vdwparam, 76
.set nb400nf_Vvdw, 80
.set nb400nf_p_tabscale, 84
.set nb400nf_VFtab, 88
.set nb400nf_invsqrta, 92
.set nb400nf_dvda, 96
.set nb400nf_p_gbtabscale, 100
.set nb400nf_GBtab, 104
.set nb400nf_p_nthreads, 108
.set nb400nf_count, 112
.set nb400nf_mtx, 116
.set nb400nf_outeriter, 120
.set nb400nf_inneriter, 124
.set nb400nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb400nf_ix, 0
.set nb400nf_iy, 16
.set nb400nf_iz, 32
.set nb400nf_iq, 48
.set nb400nf_gbtsc, 64
.set nb400nf_qq, 80
.set nb400nf_vctot, 96
.set nb400nf_half, 112
.set nb400nf_three, 128
.set nb400nf_isai, 144
.set nb400nf_isaprod, 160
.set nb400nf_gbscale, 176
.set nb400nf_is3, 192
.set nb400nf_ii3, 196
.set nb400nf_innerjjnr, 200
.set nb400nf_innerk, 204
.set nb400nf_n, 208
.set nb400nf_nn1, 212
.set nb400nf_nri, 216
.set nb400nf_facel, 220
.set nb400nf_nouter, 224
.set nb400nf_ninner, 228
.set nb400nf_salign, 232
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $236,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb400nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb400nf_p_nri(%ebp),%ecx
        movl nb400nf_p_facel(%ebp),%esi
        movl (%ecx),%ecx
        movl (%esi),%esi
        movl %ecx,nb400nf_nri(%esp)
        movl %esi,nb400nf_facel(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb400nf_nouter(%esp)
        movl %eax,nb400nf_ninner(%esp)


        movl nb400nf_p_gbtabscale(%ebp),%eax
        movss (%eax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb400nf_gbtsc(%esp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb400nf_half(%esp)
        movss nb400nf_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb400nf_half(%esp)
        movaps %xmm3,nb400nf_three(%esp)

_nb_kernel400nf_ia32_sse.nb400nf_threadloop: 
        movl  nb400nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel400nf_ia32_sse.nb400nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel400nf_ia32_sse.nb400nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb400nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb400nf_n(%esp)
        movl %ebx,nb400nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel400nf_ia32_sse.nb400nf_outerstart
        jmp _nb_kernel400nf_ia32_sse.nb400nf_end

_nb_kernel400nf_ia32_sse.nb400nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb400nf_nouter(%esp),%ebx
        movl %ebx,nb400nf_nouter(%esp)

_nb_kernel400nf_ia32_sse.nb400nf_outer: 
        movl  nb400nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb400nf_is3(%esp)            ## store is3 

        movl  nb400nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb400nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 

        movl  nb400nf_charge(%ebp),%edx
        movss (%edx,%ebx,4),%xmm3
        mulss nb400nf_facel(%esp),%xmm3
        shufps $0,%xmm3,%xmm3

        movl  nb400nf_invsqrta(%ebp),%edx       ## load invsqrta[ii]
        movss (%edx,%ebx,4),%xmm4
        shufps $0,%xmm4,%xmm4

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb400nf_pos(%ebp),%eax      ## eax = base of pos[]  

        addss (%eax,%ebx,4),%xmm0
        addss 4(%eax,%ebx,4),%xmm1
        addss 8(%eax,%ebx,4),%xmm2

        movaps %xmm3,nb400nf_iq(%esp)
        movaps %xmm4,nb400nf_isai(%esp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb400nf_ix(%esp)
        movaps %xmm1,nb400nf_iy(%esp)
        movaps %xmm2,nb400nf_iz(%esp)

        movl  %ebx,nb400nf_ii3(%esp)

        ## clear vctot 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb400nf_vctot(%esp)

        movl  nb400nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb400nf_pos(%ebp),%esi
        movl  nb400nf_faction(%ebp),%edi
        movl  nb400nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb400nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb400nf_ninner(%esp),%ecx
        movl  %ecx,nb400nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb400nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel400nf_ia32_sse.nb400nf_unroll_loop
        jmp   _nb_kernel400nf_ia32_sse.nb400nf_finish_inner
_nb_kernel400nf_ia32_sse.nb400nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb400nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 
        addl $16,nb400nf_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        ## load isa2
        movl nb400nf_invsqrta(%ebp),%esi
        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7
        movaps nb400nf_isai(%esp),%xmm2
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        mulps  %xmm3,%xmm2

        movaps %xmm2,nb400nf_isaprod(%esp)
        movaps %xmm2,%xmm1
        mulps nb400nf_gbtsc(%esp),%xmm1
        movaps %xmm1,nb400nf_gbscale(%esp)

        movl nb400nf_charge(%ebp),%esi     ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        mulps nb400nf_iq(%esp),%xmm2
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        mulps  %xmm2,%xmm3
        movaps %xmm3,nb400nf_qq(%esp)


        movl nb400nf_pos(%ebp),%esi        ## base of pos[] 

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

        ## move ix-iz to xmm4-xmm6 
        movaps nb400nf_ix(%esp),%xmm4
        movaps nb400nf_iy(%esp),%xmm5
        movaps nb400nf_iz(%esp),%xmm6

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
        movaps nb400nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb400nf_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r
        mulps nb400nf_gbscale(%esp),%xmm4

        movhlps %xmm4,%xmm5
        cvttps2pi %xmm4,%mm6
        cvttps2pi %xmm5,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        cvtpi2ps %mm7,%xmm5
        movlhps %xmm5,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $2,%mm6
        pslld $2,%mm7

        movd %eax,%mm0
        movd %ebx,%mm1
        movd %ecx,%mm2
        movd %edx,%mm3

        movl nb400nf_GBtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        ## load coulomb table
        movaps (%esi,%eax,4),%xmm4
        movaps (%esi,%ebx,4),%xmm5
        movaps (%esi,%ecx,4),%xmm6
        movaps (%esi,%edx,4),%xmm7
        ## transpose, using xmm3 for scratch
        movaps %xmm6,%xmm3
        shufps $0xEE,%xmm7,%xmm3
        shufps $0x44,%xmm7,%xmm6
        movaps %xmm4,%xmm7
        shufps $0xEE,%xmm5,%xmm7
        shufps $0x44,%xmm5,%xmm4
        movaps %xmm4,%xmm5
        shufps $0xDD,%xmm6,%xmm5
        shufps $0x88,%xmm6,%xmm4
        movaps %xmm7,%xmm6
        shufps $0x88,%xmm3,%xmm6
        shufps $0xDD,%xmm3,%xmm7
        ## coulomb table ready, in xmm4-xmm7    

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        movaps nb400nf_qq(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        addps  nb400nf_vctot(%esp),%xmm5
        movaps %xmm5,nb400nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $4,nb400nf_innerk(%esp)
        jl    _nb_kernel400nf_ia32_sse.nb400nf_finish_inner
        jmp   _nb_kernel400nf_ia32_sse.nb400nf_unroll_loop
_nb_kernel400nf_ia32_sse.nb400nf_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb400nf_innerk(%esp)
        movl  nb400nf_innerk(%esp),%edx
        andl  $2,%edx
        jnz   _nb_kernel400nf_ia32_sse.nb400nf_dopair
        jmp   _nb_kernel400nf_ia32_sse.nb400nf_checksingle
_nb_kernel400nf_ia32_sse.nb400nf_dopair: 
        movl  nb400nf_innerjjnr(%esp),%ecx

        movl  (%ecx),%eax
        movl  4(%ecx),%ebx
        addl $8,nb400nf_innerjjnr(%esp)

        xorps %xmm2,%xmm2
        movaps %xmm2,%xmm6

        ## load isa2
        movl nb400nf_invsqrta(%ebp),%esi
        movss (%esi,%eax,4),%xmm2
        movss (%esi,%ebx,4),%xmm3
        unpcklps %xmm3,%xmm2    ## isa2 in xmm3(0,1)
        mulps  nb400nf_isai(%esp),%xmm2
        movaps %xmm2,nb400nf_isaprod(%esp)
        movaps %xmm2,%xmm1
        mulps nb400nf_gbtsc(%esp),%xmm1
        movaps %xmm1,nb400nf_gbscale(%esp)

        movl nb400nf_charge(%ebp),%esi     ## base of charge[]  
        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ebx,4),%xmm6
        unpcklps %xmm6,%xmm3 ## constant 00001000 ;# xmm3(0,1) has the charges 

        mulps  nb400nf_iq(%esp),%xmm2
        mulps  %xmm2,%xmm3
        movaps %xmm3,nb400nf_qq(%esp)

        movl nb400nf_pos(%ebp),%edi

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        ## move coordinates to xmm0-xmm2 
        movlps (%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2
        movhps (%edi,%ebx,4),%xmm1
        movss 8(%edi,%ebx,4),%xmm0

        movlhps %xmm7,%xmm3

        shufps $0,%xmm0,%xmm2

        movaps %xmm1,%xmm0

        shufps $136,%xmm2,%xmm2 ## constant 10001000

        shufps $136,%xmm0,%xmm0 ## constant 10001000
        shufps $221,%xmm1,%xmm1 ## constant 11011101

        movl   nb400nf_faction(%ebp),%edi
        ## move ix-iz to xmm4-xmm6 
        xorps   %xmm7,%xmm7

        movaps nb400nf_ix(%esp),%xmm4
        movaps nb400nf_iy(%esp),%xmm5
        movaps nb400nf_iz(%esp),%xmm6

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
        movaps nb400nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb400nf_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb400nf_gbscale(%esp),%xmm4

        cvttps2pi %xmm4,%mm6    ## mm6 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6

        movl nb400nf_GBtab(%ebp),%esi
        movd %mm6,%ecx
        psrlq $32,%mm6
        movd %mm6,%edx

        ## load coulomb table
        movaps (%esi,%ecx,4),%xmm4
        movaps (%esi,%edx,4),%xmm7
        ## transpose, using xmm3 for scratch
        movaps %xmm4,%xmm6
        unpcklps %xmm7,%xmm4    ## Y1 Y2 F1 F2 
        unpckhps %xmm7,%xmm6    ## G1 G2 H1 H2
        movhlps  %xmm4,%xmm5    ## F1 F2 
        movhlps  %xmm6,%xmm7    ## H1 H2
        ## coulomb table ready, in xmm4-xmm7    

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        movaps nb400nf_qq(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        addps  nb400nf_vctot(%esp),%xmm5
        movaps %xmm5,nb400nf_vctot(%esp)

_nb_kernel400nf_ia32_sse.nb400nf_checksingle:   
        movl  nb400nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel400nf_ia32_sse.nb400nf_dosingle
        jmp    _nb_kernel400nf_ia32_sse.nb400nf_updateouterdata
_nb_kernel400nf_ia32_sse.nb400nf_dosingle: 
        movl nb400nf_charge(%ebp),%esi
        movl nb400nf_invsqrta(%ebp),%edx
        movl nb400nf_pos(%ebp),%edi
        movl  nb400nf_innerjjnr(%esp),%ecx
        movl  (%ecx),%eax
        xorps  %xmm2,%xmm2
        movaps %xmm2,%xmm6
        movss (%edx,%eax,4),%xmm2       ## isa2
        mulss nb400nf_isai(%esp),%xmm2
        movss %xmm2,nb400nf_isaprod(%esp)
        movss %xmm2,%xmm1
        mulss nb400nf_gbtsc(%esp),%xmm1
        movss %xmm1,nb400nf_gbscale(%esp)

        mulss  nb400nf_iq(%esp),%xmm2
        movss (%esi,%eax,4),%xmm6       ## xmm6(0) has the charge       
        mulss  %xmm2,%xmm6
        movss %xmm6,nb400nf_qq(%esp)

        leal  (%eax,%eax,2),%eax

        ## move coordinates to xmm0-xmm2 
        movss (%edi,%eax,4),%xmm0
        movss 4(%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2

        movss nb400nf_ix(%esp),%xmm4
        movss nb400nf_iy(%esp),%xmm5
        movss nb400nf_iz(%esp),%xmm6

        ## calc dr 
        subss %xmm0,%xmm4
        subss %xmm1,%xmm5
        subss %xmm2,%xmm6

        ## square it 
        mulss %xmm4,%xmm4
        mulss %xmm5,%xmm5
        mulss %xmm6,%xmm6
        addss %xmm5,%xmm4
        addss %xmm6,%xmm4
        ## rsq in xmm4 

        rsqrtss %xmm4,%xmm5
        ## lookup seed in xmm5 
        movss %xmm5,%xmm2
        mulss %xmm5,%xmm5
        movss nb400nf_three(%esp),%xmm1
        mulss %xmm4,%xmm5       ## rsq*lu*lu                    
        movss nb400nf_half(%esp),%xmm0
        subss %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulss %xmm2,%xmm1
        mulss %xmm1,%xmm0       ## xmm0=rinv 

        mulss %xmm0,%xmm4       ## xmm4=r 
        mulss nb400nf_gbscale(%esp),%xmm4

        cvttss2si %xmm4,%ebx    ## mm6 contain lu indices 
        cvtsi2ss %ebx,%xmm6
        subss %xmm6,%xmm4
        movss %xmm4,%xmm1       ## xmm1=eps 
        movss %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%ebx

        movl nb400nf_GBtab(%ebp),%esi

        movaps (%esi,%ebx,4),%xmm4
        movhlps %xmm4,%xmm6
        movaps %xmm4,%xmm5
        movaps %xmm6,%xmm7
        shufps $1,%xmm5,%xmm5
        shufps $1,%xmm7,%xmm7
        ## table ready in xmm4-xmm7 

        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp      
        movss nb400nf_qq(%esp),%xmm3
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 
        mulss  %xmm3,%xmm5 ## vcoul=qq*VV  
        addss  nb400nf_vctot(%esp),%xmm5
        movss %xmm5,nb400nf_vctot(%esp)
_nb_kernel400nf_ia32_sse.nb400nf_updateouterdata: 
        ## get n from stack
        movl nb400nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb400nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb400nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb400nf_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb400nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel400nf_ia32_sse.nb400nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb400nf_n(%esp)
        jmp _nb_kernel400nf_ia32_sse.nb400nf_outer
_nb_kernel400nf_ia32_sse.nb400nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb400nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel400nf_ia32_sse.nb400nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel400nf_ia32_sse.nb400nf_threadloop
_nb_kernel400nf_ia32_sse.nb400nf_end: 
        emms

        movl nb400nf_nouter(%esp),%eax
        movl nb400nf_ninner(%esp),%ebx
        movl nb400nf_outeriter(%ebp),%ecx
        movl nb400nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb400nf_salign(%esp),%eax
        addl %eax,%esp
        addl $236,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret

