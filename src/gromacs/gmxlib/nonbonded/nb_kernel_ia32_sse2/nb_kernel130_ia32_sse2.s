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

.globl nb_kernel130_ia32_sse2
.globl _nb_kernel130_ia32_sse2
nb_kernel130_ia32_sse2: 
_nb_kernel130_ia32_sse2:        
.set nb130_p_nri, 8
.set nb130_iinr, 12
.set nb130_jindex, 16
.set nb130_jjnr, 20
.set nb130_shift, 24
.set nb130_shiftvec, 28
.set nb130_fshift, 32
.set nb130_gid, 36
.set nb130_pos, 40
.set nb130_faction, 44
.set nb130_charge, 48
.set nb130_p_facel, 52
.set nb130_argkrf, 56
.set nb130_argcrf, 60
.set nb130_Vc, 64
.set nb130_type, 68
.set nb130_p_ntype, 72
.set nb130_vdwparam, 76
.set nb130_Vvdw, 80
.set nb130_p_tabscale, 84
.set nb130_VFtab, 88
.set nb130_invsqrta, 92
.set nb130_dvda, 96
.set nb130_p_gbtabscale, 100
.set nb130_GBtab, 104
.set nb130_p_nthreads, 108
.set nb130_count, 112
.set nb130_mtx, 116
.set nb130_outeriter, 120
.set nb130_inneriter, 124
.set nb130_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb130_ix, 0
.set nb130_iy, 16
.set nb130_iz, 32
.set nb130_iq, 48
.set nb130_dx, 64
.set nb130_dy, 80
.set nb130_dz, 96
.set nb130_c6, 112
.set nb130_c12, 128
.set nb130_tsc, 144
.set nb130_fstmp, 160
.set nb130_vctot, 176
.set nb130_Vvdwtot, 192
.set nb130_fix, 208
.set nb130_fiy, 224
.set nb130_fiz, 240
.set nb130_half, 256
.set nb130_three, 272
.set nb130_two, 288
.set nb130_is3, 336
.set nb130_ii3, 340
.set nb130_ntia, 344
.set nb130_innerjjnr, 348
.set nb130_innerk, 352
.set nb130_n, 356
.set nb130_nn1, 360
.set nb130_nri, 364
.set nb130_facel, 368                         ## uses 8 bytes
.set nb130_ntype, 376
.set nb130_nouter, 380
.set nb130_ninner, 384
.set nb130_salign, 388
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $392,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb130_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb130_p_nri(%ebp),%ecx
        movl nb130_p_facel(%ebp),%esi
        movl nb130_p_ntype(%ebp),%edi
        movl (%ecx),%ecx
        movsd (%esi),%xmm7
        movl (%edi),%edi
        movl %ecx,nb130_nri(%esp)
        movsd %xmm7,nb130_facel(%esp)
        movl %edi,nb130_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb130_nouter(%esp)
        movl %eax,nb130_ninner(%esp)

        movl nb130_p_tabscale(%ebp),%eax
        movsd (%eax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb130_tsc(%esp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb130_half(%esp)
        movl %ebx,nb130_half+4(%esp)
        movsd nb130_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb130_half(%esp)
        movapd %xmm2,nb130_two(%esp)
        movapd %xmm3,nb130_three(%esp)

_nb_kernel130_ia32_sse2.nb130_threadloop: 
        movl  nb130_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel130_ia32_sse2.nb130_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel130_ia32_sse2.nb130_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb130_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb130_n(%esp)
        movl %ebx,nb130_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel130_ia32_sse2.nb130_outerstart
        jmp _nb_kernel130_ia32_sse2.nb130_end

_nb_kernel130_ia32_sse2.nb130_outerstart: 
        ## ebx contains number of outer iterations
        addl nb130_nouter(%esp),%ebx
        movl %ebx,nb130_nouter(%esp)

_nb_kernel130_ia32_sse2.nb130_outer: 
        movl  nb130_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb130_is3(%esp)      ## store is3 

        movl  nb130_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb130_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movl  nb130_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        mulsd nb130_facel(%esp),%xmm3
        shufpd $0,%xmm3,%xmm3

        movl  nb130_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb130_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb130_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb130_pos(%ebp),%eax      ## eax = base of pos[]  

        addsd (%eax,%ebx,8),%xmm0
        addsd 8(%eax,%ebx,8),%xmm1
        addsd 16(%eax,%ebx,8),%xmm2

        movapd %xmm3,nb130_iq(%esp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb130_ix(%esp)
        movapd %xmm1,nb130_iy(%esp)
        movapd %xmm2,nb130_iz(%esp)

        movl  %ebx,nb130_ii3(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb130_vctot(%esp)
        movapd %xmm4,nb130_Vvdwtot(%esp)
        movapd %xmm4,nb130_fix(%esp)
        movapd %xmm4,nb130_fiy(%esp)
        movapd %xmm4,nb130_fiz(%esp)

        movl  nb130_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb130_pos(%ebp),%esi
        movl  nb130_faction(%ebp),%edi
        movl  nb130_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb130_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb130_ninner(%esp),%ecx
        movl  %ecx,nb130_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb130_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel130_ia32_sse2.nb130_unroll_loop
        jmp   _nb_kernel130_ia32_sse2.nb130_checksingle
_nb_kernel130_ia32_sse2.nb130_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb130_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        addl $8,nb130_innerjjnr(%esp)                   ## advance pointer (unrolled 2) 

        movl nb130_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3

        movapd nb130_iq(%esp),%xmm5
        mulpd %xmm5,%xmm3               ## qq 

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movl nb130_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb130_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb130_ntia(%esp),%edi
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
        movapd %xmm4,nb130_c6(%esp)
        movapd %xmm6,nb130_c12(%esp)

        movl nb130_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move two coordinates to xmm0-xmm2    
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2
        movhpd (%esi,%ebx,8),%xmm0
        movhpd 8(%esi,%ebx,8),%xmm1
        movhpd 16(%esi,%ebx,8),%xmm2

        ## move ix-iz to xmm4-xmm6 
        movapd nb130_ix(%esp),%xmm4
        movapd nb130_iy(%esp),%xmm5
        movapd nb130_iz(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb130_dx(%esp)
        movapd %xmm5,nb130_dy(%esp)
        movapd %xmm6,nb130_dz(%esp)
        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4
        ## rsq in xmm4 

        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm5,%xmm2    ## lu in low xmm2 

        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulpd %xmm2,%xmm2       ## lu*lu 
        movapd nb130_three(%esp),%xmm1
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb130_half(%esp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb130_three(%esp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb130_half(%esp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=rinv 
        movapd %xmm0,%xmm6
        movapd %xmm4,%xmm1
        mulpd  %xmm3,%xmm6  ## vcoul = rinv*qq
        movapd %xmm6,%xmm3
        mulpd  %xmm0,%xmm3

        movapd %xmm3,nb130_fstmp(%esp)

        addpd  nb130_vctot(%esp),%xmm6
        movapd %xmm6,nb130_vctot(%esp)

        ## LJ table interaction. xmm0=rinv, xmm4=rsq

        mulpd %xmm0,%xmm4       ## xmm4=r 
        mulpd nb130_tsc(%esp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $3,%mm6           ## idx *= 8 

        movd %eax,%mm0
        movd %ebx,%mm1

        movl nb130_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx

        ## dispersion 
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
        ## dispersion table ready, in xmm4-xmm7         
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb130_two(%esp),%xmm7    ## two*Heps2 
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb130_c6(%esp),%xmm4
        mulpd  %xmm4,%xmm7       ## fijD 
        mulpd  %xmm4,%xmm5       ## Vvdw6 

        ## put scalar force on stack Update Vvdwtot directly 
        addpd  nb130_Vvdwtot(%esp),%xmm5
        movapd nb130_fstmp(%esp),%xmm3
        mulpd  nb130_tsc(%esp),%xmm7
        subpd  %xmm7,%xmm3
        movapd %xmm3,nb130_fstmp(%esp)
        movapd %xmm5,nb130_Vvdwtot(%esp)

        ## repulsion 
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

        ## table ready, in xmm4-xmm7    
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb130_two(%esp),%xmm7    ## two*Heps2 
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb130_c12(%esp),%xmm4
        mulpd  %xmm4,%xmm7
        mulpd  %xmm4,%xmm5

        addpd  nb130_Vvdwtot(%esp),%xmm5
        movapd nb130_fstmp(%esp),%xmm3
        mulpd  nb130_tsc(%esp),%xmm7
        subpd  %xmm7,%xmm3
        movapd %xmm5,nb130_Vvdwtot(%esp)

        mulpd  %xmm0,%xmm3

        movapd nb130_dx(%esp),%xmm0
        movapd nb130_dy(%esp),%xmm1
        movapd nb130_dz(%esp),%xmm2

        movd %mm0,%eax
        movd %mm1,%ebx

        movl   nb130_faction(%ebp),%edi
        mulpd  %xmm3,%xmm0
        mulpd  %xmm3,%xmm1
        mulpd  %xmm3,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movapd nb130_fix(%esp),%xmm3
        movapd nb130_fiy(%esp),%xmm4
        movapd nb130_fiz(%esp),%xmm5
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm5
        movapd %xmm3,nb130_fix(%esp)
        movapd %xmm4,nb130_fiy(%esp)
        movapd %xmm5,nb130_fiz(%esp)
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
        subl $2,nb130_innerk(%esp)
        jl    _nb_kernel130_ia32_sse2.nb130_checksingle
        jmp   _nb_kernel130_ia32_sse2.nb130_unroll_loop

_nb_kernel130_ia32_sse2.nb130_checksingle:      
        movl  nb130_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel130_ia32_sse2.nb130_dosingle
        jmp    _nb_kernel130_ia32_sse2.nb130_updateouterdata
_nb_kernel130_ia32_sse2.nb130_dosingle: 
        movl nb130_charge(%ebp),%esi
        movl nb130_pos(%ebp),%edi
        movl  nb130_innerjjnr(%esp),%ecx
        xorpd %xmm3,%xmm3
        movl  (%ecx),%eax

        movlpd (%esi,%eax,8),%xmm3
        movapd nb130_iq(%esp),%xmm5
        mulpd %xmm5,%xmm3               ## qq 

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movl nb130_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb130_vdwparam(%ebp),%esi
        shll %eax
        movl nb130_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 
        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movapd %xmm4,nb130_c6(%esp)
        movapd %xmm6,nb130_c12(%esp)

        movl nb130_pos(%ebp),%esi        ## base of pos[] 

        leal (%eax,%eax,2),%eax    ## replace jnr with j3 

        ## move two coordinates to xmm0-xmm2    
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ix-iz to xmm4-xmm6 
        movapd nb130_ix(%esp),%xmm4
        movapd nb130_iy(%esp),%xmm5
        movapd nb130_iz(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb130_dx(%esp)
        movapd %xmm5,nb130_dy(%esp)
        movapd %xmm6,nb130_dz(%esp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4
        ## rsq in xmm4 

        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm5,%xmm2    ## lu in low xmm2 

        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulsd %xmm2,%xmm2       ## lu*lu 
        movapd nb130_three(%esp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb130_half(%esp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb130_three(%esp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb130_half(%esp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=rinv 
        movapd %xmm0,%xmm6
        movapd %xmm4,%xmm1
        mulsd  %xmm3,%xmm6      ## xmm6=vcoul=qq*rinv
        movapd %xmm0,%xmm1
        movapd %xmm6,%xmm3
        mulsd  %xmm0,%xmm3

        movsd %xmm3,nb130_fstmp(%esp)

        addsd  nb130_vctot(%esp),%xmm6
        movsd %xmm6,nb130_vctot(%esp)

        ## LJ table interaction. xmm0=rinv, cmm4=rsq

        mulsd %xmm0,%xmm4       ## xmm4=r 
        mulsd nb130_tsc(%esp),%xmm4

        cvttsd2si %xmm4,%ebx    ## mm6 = lu idx 
        cvtsi2sd %ebx,%xmm5
        subsd %xmm5,%xmm4
        movsd %xmm4,%xmm1       ## xmm1=eps 
        movsd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll  $3,%ebx

        movl nb130_VFtab(%ebp),%esi

        ## dispersion 
        movlpd (%esi,%ebx,8),%xmm4      ## Y1   
        movhpd 8(%esi,%ebx,8),%xmm4     ## Y1 F1        
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%ebx,8),%xmm6    ## G1
        movhpd 24(%esi,%ebx,8),%xmm6    ## G1 H1        
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## dispersion table ready, in xmm4-xmm7         
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb130_two(%esp),%xmm7    ## two*Heps2 
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movsd nb130_c6(%esp),%xmm4
        mulsd  %xmm4,%xmm7       ## fijD 
        mulsd  %xmm4,%xmm5       ## Vvdw6 

        ## put scalar force on stack Update Vvdwtot directly 
        addsd  nb130_Vvdwtot(%esp),%xmm5
        movsd nb130_fstmp(%esp),%xmm3
        mulsd  nb130_tsc(%esp),%xmm7
        subsd  %xmm7,%xmm3
        movsd %xmm3,nb130_fstmp(%esp)
        movsd %xmm5,nb130_Vvdwtot(%esp)

        ## repulsion 
        movlpd 32(%esi,%ebx,8),%xmm4    ## Y1   
        movhpd 40(%esi,%ebx,8),%xmm4    ## Y1 F1        

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 48(%esi,%ebx,8),%xmm6    ## G1
        movhpd 56(%esi,%ebx,8),%xmm6    ## G1 H1        

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 

        ## table ready, in xmm4-xmm7    
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb130_two(%esp),%xmm7    ## two*Heps2 
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movsd nb130_c12(%esp),%xmm4
        mulsd  %xmm4,%xmm7
        mulsd  %xmm4,%xmm5

        addsd  nb130_Vvdwtot(%esp),%xmm5
        movsd nb130_fstmp(%esp),%xmm3
        mulsd  nb130_tsc(%esp),%xmm7
        subsd  %xmm7,%xmm3
        movsd %xmm5,nb130_Vvdwtot(%esp)

        mulsd  %xmm0,%xmm3

        movsd nb130_dx(%esp),%xmm0
        movsd nb130_dy(%esp),%xmm1
        movsd nb130_dz(%esp),%xmm2

        movl   nb130_faction(%ebp),%edi
        mulsd  %xmm3,%xmm0
        mulsd  %xmm3,%xmm1
        mulsd  %xmm3,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movlpd nb130_fix(%esp),%xmm3
        movlpd nb130_fiy(%esp),%xmm4
        movlpd nb130_fiz(%esp),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movlpd %xmm3,nb130_fix(%esp)
        movlpd %xmm4,nb130_fiy(%esp)
        movlpd %xmm5,nb130_fiz(%esp)
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

_nb_kernel130_ia32_sse2.nb130_updateouterdata: 
        movl  nb130_ii3(%esp),%ecx
        movl  nb130_faction(%ebp),%edi
        movl  nb130_fshift(%ebp),%esi
        movl  nb130_is3(%esp),%edx

        ## accumulate i forces in xmm0, xmm1, xmm2 
        movapd nb130_fix(%esp),%xmm0
        movapd nb130_fiy(%esp),%xmm1
        movapd nb130_fiz(%esp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addpd  %xmm3,%xmm0
        addpd  %xmm4,%xmm1
        addpd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

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

        ## increment fshift force  
        movsd  (%esi,%edx,8),%xmm3
        movsd  8(%esi,%edx,8),%xmm4
        movsd  16(%esi,%edx,8),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movsd  %xmm3,(%esi,%edx,8)
        movsd  %xmm4,8(%esi,%edx,8)
        movsd  %xmm5,16(%esi,%edx,8)

        ## get n from stack
        movl nb130_n(%esp),%esi
        ## get group index for i particle 
        movl  nb130_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb130_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb130_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb130_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb130_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb130_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel130_ia32_sse2.nb130_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb130_n(%esp)
        jmp _nb_kernel130_ia32_sse2.nb130_outer
_nb_kernel130_ia32_sse2.nb130_outerend: 
        ## check if more outer neighborlists remain
        movl  nb130_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel130_ia32_sse2.nb130_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel130_ia32_sse2.nb130_threadloop
_nb_kernel130_ia32_sse2.nb130_end: 
        emms

        movl nb130_nouter(%esp),%eax
        movl nb130_ninner(%esp),%ebx
        movl nb130_outeriter(%ebp),%ecx
        movl nb130_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb130_salign(%esp),%eax
        addl %eax,%esp
        addl $392,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret


.globl nb_kernel130nf_ia32_sse2
.globl _nb_kernel130nf_ia32_sse2
nb_kernel130nf_ia32_sse2:       
_nb_kernel130nf_ia32_sse2:      
.set nb130nf_p_nri, 8
.set nb130nf_iinr, 12
.set nb130nf_jindex, 16
.set nb130nf_jjnr, 20
.set nb130nf_shift, 24
.set nb130nf_shiftvec, 28
.set nb130nf_fshift, 32
.set nb130nf_gid, 36
.set nb130nf_pos, 40
.set nb130nf_faction, 44
.set nb130nf_charge, 48
.set nb130nf_p_facel, 52
.set nb130nf_argkrf, 56
.set nb130nf_argcrf, 60
.set nb130nf_Vc, 64
.set nb130nf_type, 68
.set nb130nf_p_ntype, 72
.set nb130nf_vdwparam, 76
.set nb130nf_Vvdw, 80
.set nb130nf_p_tabscale, 84
.set nb130nf_VFtab, 88
.set nb130nf_invsqrta, 92
.set nb130nf_dvda, 96
.set nb130nf_p_gbtabscale, 100
.set nb130nf_GBtab, 104
.set nb130nf_p_nthreads, 108
.set nb130nf_count, 112
.set nb130nf_mtx, 116
.set nb130nf_outeriter, 120
.set nb130nf_inneriter, 124
.set nb130nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb130nf_ix, 0
.set nb130nf_iy, 16
.set nb130nf_iz, 32
.set nb130nf_iq, 48
.set nb130nf_dx, 64
.set nb130nf_dy, 80
.set nb130nf_dz, 96
.set nb130nf_c6, 112
.set nb130nf_c12, 128
.set nb130nf_tsc, 144
.set nb130nf_vctot, 176
.set nb130nf_Vvdwtot, 192
.set nb130nf_half, 256
.set nb130nf_three, 272
.set nb130nf_is3, 336
.set nb130nf_ii3, 340
.set nb130nf_ntia, 344
.set nb130nf_innerjjnr, 348
.set nb130nf_innerk, 352
.set nb130nf_n, 356
.set nb130nf_nn1, 360
.set nb130nf_nri, 364
.set nb130nf_facel, 368                         ## uses 8 bytes
.set nb130nf_ntype, 376
.set nb130nf_nouter, 380
.set nb130nf_ninner, 384
.set nb130nf_salign, 388
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $392,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb130nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb130nf_p_nri(%ebp),%ecx
        movl nb130nf_p_facel(%ebp),%esi
        movl nb130nf_p_ntype(%ebp),%edi
        movl (%ecx),%ecx
        movsd (%esi),%xmm7
        movl (%edi),%edi
        movl %ecx,nb130nf_nri(%esp)
        movsd %xmm7,nb130nf_facel(%esp)
        movl %edi,nb130nf_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb130nf_nouter(%esp)
        movl %eax,nb130nf_ninner(%esp)

        movl nb130nf_p_tabscale(%ebp),%eax
        movsd (%eax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb130nf_tsc(%esp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb130nf_half(%esp)
        movl %ebx,nb130nf_half+4(%esp)
        movsd nb130nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb130nf_half(%esp)
        movapd %xmm3,nb130nf_three(%esp)

_nb_kernel130nf_ia32_sse2.nb130nf_threadloop: 
        movl  nb130nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel130nf_ia32_sse2.nb130nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel130nf_ia32_sse2.nb130nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb130nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb130nf_n(%esp)
        movl %ebx,nb130nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel130nf_ia32_sse2.nb130nf_outerstart
        jmp _nb_kernel130nf_ia32_sse2.nb130nf_end

_nb_kernel130nf_ia32_sse2.nb130nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb130nf_nouter(%esp),%ebx
        movl %ebx,nb130nf_nouter(%esp)

_nb_kernel130nf_ia32_sse2.nb130nf_outer: 
        movl  nb130nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb130nf_is3(%esp)            ## store is3 

        movl  nb130nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb130nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movl  nb130nf_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        mulsd nb130nf_facel(%esp),%xmm3
        shufpd $0,%xmm3,%xmm3

        movl  nb130nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb130nf_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb130nf_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb130nf_pos(%ebp),%eax      ## eax = base of pos[]  

        addsd (%eax,%ebx,8),%xmm0
        addsd 8(%eax,%ebx,8),%xmm1
        addsd 16(%eax,%ebx,8),%xmm2

        movapd %xmm3,nb130nf_iq(%esp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb130nf_ix(%esp)
        movapd %xmm1,nb130nf_iy(%esp)
        movapd %xmm2,nb130nf_iz(%esp)

        movl  %ebx,nb130nf_ii3(%esp)

        ## clear vctot
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb130nf_vctot(%esp)
        movapd %xmm4,nb130nf_Vvdwtot(%esp)

        movl  nb130nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb130nf_pos(%ebp),%esi
        movl  nb130nf_faction(%ebp),%edi
        movl  nb130nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb130nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb130nf_ninner(%esp),%ecx
        movl  %ecx,nb130nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb130nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel130nf_ia32_sse2.nb130nf_unroll_loop
        jmp   _nb_kernel130nf_ia32_sse2.nb130nf_checksingle
_nb_kernel130nf_ia32_sse2.nb130nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb130nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        addl $8,nb130nf_innerjjnr(%esp)                 ## advance pointer (unrolled 2) 

        movl nb130nf_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3

        movapd nb130nf_iq(%esp),%xmm5
        mulpd %xmm5,%xmm3               ## qq 

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movl nb130nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb130nf_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb130nf_ntia(%esp),%edi
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
        movapd %xmm4,nb130nf_c6(%esp)
        movapd %xmm6,nb130nf_c12(%esp)

        movl nb130nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move two coordinates to xmm0-xmm2    
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2
        movhpd (%esi,%ebx,8),%xmm0
        movhpd 8(%esi,%ebx,8),%xmm1
        movhpd 16(%esi,%ebx,8),%xmm2

        ## move ix-iz to xmm4-xmm6 
        movapd nb130nf_ix(%esp),%xmm4
        movapd nb130nf_iy(%esp),%xmm5
        movapd nb130nf_iz(%esp),%xmm6

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
        ## rsq in xmm4 

        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm5,%xmm2    ## lu in low xmm2 

        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulpd %xmm2,%xmm2       ## lu*lu 
        movapd nb130nf_three(%esp),%xmm1
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb130nf_half(%esp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb130nf_three(%esp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb130nf_half(%esp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=rinv 
        movapd %xmm0,%xmm6
        mulpd  %xmm3,%xmm6  ## vcoul

        addpd  nb130nf_vctot(%esp),%xmm6
        movapd %xmm6,nb130nf_vctot(%esp)

        ## LJ table interaction. xmm0=rinv, xmm4=rsq

        mulpd %xmm0,%xmm4       ## xmm4=r 
        mulpd nb130nf_tsc(%esp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $3,%mm6           ## idx *= 8 

        movl nb130nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx

        ## dispersion 
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
        ## dispersion table ready, in xmm4-xmm7         
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb130nf_c6(%esp),%xmm4
        mulpd  %xmm4,%xmm5       ## Vvdw6 

        ##Update Vvdwtot directly 
        addpd  nb130nf_Vvdwtot(%esp),%xmm5
        movapd %xmm5,nb130nf_Vvdwtot(%esp)

        ## repulsion 
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

        ## table ready, in xmm4-xmm7    
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb130nf_c12(%esp),%xmm4
        mulpd  %xmm4,%xmm5

        addpd  nb130nf_Vvdwtot(%esp),%xmm5
        movapd %xmm5,nb130nf_Vvdwtot(%esp)

        ## should we do one more iteration? 
        subl $2,nb130nf_innerk(%esp)
        jl    _nb_kernel130nf_ia32_sse2.nb130nf_checksingle
        jmp   _nb_kernel130nf_ia32_sse2.nb130nf_unroll_loop

_nb_kernel130nf_ia32_sse2.nb130nf_checksingle:  
        movl  nb130nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel130nf_ia32_sse2.nb130nf_dosingle
        jmp    _nb_kernel130nf_ia32_sse2.nb130nf_updateouterdata
_nb_kernel130nf_ia32_sse2.nb130nf_dosingle: 
        movl nb130nf_charge(%ebp),%esi
        movl nb130nf_pos(%ebp),%edi
        movl  nb130nf_innerjjnr(%esp),%ecx
        xorpd %xmm3,%xmm3
        movl  (%ecx),%eax

        movlpd (%esi,%eax,8),%xmm3
        movapd nb130nf_iq(%esp),%xmm5
        mulpd %xmm5,%xmm3               ## qq 

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movl nb130nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb130nf_vdwparam(%ebp),%esi
        shll %eax
        movl nb130nf_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 
        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movapd %xmm4,nb130nf_c6(%esp)
        movapd %xmm6,nb130nf_c12(%esp)

        movl nb130nf_pos(%ebp),%esi        ## base of pos[] 

        leal (%eax,%eax,2),%eax    ## replace jnr with j3 

        ## move two coordinates to xmm0-xmm2    
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ix-iz to xmm4-xmm6 
        movapd nb130nf_ix(%esp),%xmm4
        movapd nb130nf_iy(%esp),%xmm5
        movapd nb130nf_iz(%esp),%xmm6

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
        ## rsq in xmm4 

        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm5,%xmm2    ## lu in low xmm2 

        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulsd %xmm2,%xmm2       ## lu*lu 
        movapd nb130nf_three(%esp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb130nf_half(%esp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb130nf_three(%esp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb130nf_half(%esp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=rinv 
        movapd %xmm0,%xmm6
        movapd %xmm4,%xmm1
        mulsd  %xmm3,%xmm6      ## xmm6=vcoul=qq*rinv

        addsd  nb130nf_vctot(%esp),%xmm6
        movsd %xmm6,nb130nf_vctot(%esp)

        ## LJ table interaction. xmm0=rinv, cmm4=rsq

        mulsd %xmm0,%xmm4       ## xmm4=r 
        mulsd nb130nf_tsc(%esp),%xmm4

        cvttsd2si %xmm4,%ebx    ## mm6 = lu idx 
        cvtsi2sd %ebx,%xmm5
        subsd %xmm5,%xmm4
        movsd %xmm4,%xmm1       ## xmm1=eps 
        movsd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll  $3,%ebx

        movl nb130nf_VFtab(%ebp),%esi

        ## dispersion 
        movlpd (%esi,%ebx,8),%xmm4      ## Y1   
        movhpd 8(%esi,%ebx,8),%xmm4     ## Y1 F1        
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%ebx,8),%xmm6    ## G1
        movhpd 24(%esi,%ebx,8),%xmm6    ## G1 H1        
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## dispersion table ready, in xmm4-xmm7         
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movsd nb130nf_c6(%esp),%xmm4
        mulsd  %xmm4,%xmm5       ## Vvdw6 

        ## put scalar force on stack Update Vvdwtot directly 
        addsd  nb130nf_Vvdwtot(%esp),%xmm5
        movsd %xmm5,nb130nf_Vvdwtot(%esp)

        ## repulsion 
        movlpd 32(%esi,%ebx,8),%xmm4    ## Y1   
        movhpd 40(%esi,%ebx,8),%xmm4    ## Y1 F1        

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 48(%esi,%ebx,8),%xmm6    ## G1
        movhpd 56(%esi,%ebx,8),%xmm6    ## G1 H1        

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 

        ## table ready, in xmm4-xmm7    
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movsd nb130nf_c12(%esp),%xmm4
        mulsd  %xmm4,%xmm5

        addsd  nb130nf_Vvdwtot(%esp),%xmm5
        movsd %xmm5,nb130nf_Vvdwtot(%esp)

_nb_kernel130nf_ia32_sse2.nb130nf_updateouterdata: 
        ## get n from stack
        movl nb130nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb130nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb130nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb130nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb130nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb130nf_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb130nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel130nf_ia32_sse2.nb130nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb130nf_n(%esp)
        jmp _nb_kernel130nf_ia32_sse2.nb130nf_outer
_nb_kernel130nf_ia32_sse2.nb130nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb130nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel130nf_ia32_sse2.nb130nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel130nf_ia32_sse2.nb130nf_threadloop
_nb_kernel130nf_ia32_sse2.nb130nf_end: 
        emms

        movl nb130nf_nouter(%esp),%eax
        movl nb130nf_ninner(%esp),%ebx
        movl nb130nf_outeriter(%ebp),%ecx
        movl nb130nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb130nf_salign(%esp),%eax
        addl %eax,%esp
        addl $392,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret


