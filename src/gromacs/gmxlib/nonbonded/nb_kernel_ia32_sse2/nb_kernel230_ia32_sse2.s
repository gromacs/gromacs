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



.globl nb_kernel230_ia32_sse2
.globl _nb_kernel230_ia32_sse2
nb_kernel230_ia32_sse2: 
_nb_kernel230_ia32_sse2:        
.set nb230_p_nri, 8
.set nb230_iinr, 12
.set nb230_jindex, 16
.set nb230_jjnr, 20
.set nb230_shift, 24
.set nb230_shiftvec, 28
.set nb230_fshift, 32
.set nb230_gid, 36
.set nb230_pos, 40
.set nb230_faction, 44
.set nb230_charge, 48
.set nb230_p_facel, 52
.set nb230_argkrf, 56
.set nb230_argcrf, 60
.set nb230_Vc, 64
.set nb230_type, 68
.set nb230_p_ntype, 72
.set nb230_vdwparam, 76
.set nb230_Vvdw, 80
.set nb230_p_tabscale, 84
.set nb230_VFtab, 88
.set nb230_invsqrta, 92
.set nb230_dvda, 96
.set nb230_p_gbtabscale, 100
.set nb230_GBtab, 104
.set nb230_p_nthreads, 108
.set nb230_count, 112
.set nb230_mtx, 116
.set nb230_outeriter, 120
.set nb230_inneriter, 124
.set nb230_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb230_ix, 0
.set nb230_iy, 16
.set nb230_iz, 32
.set nb230_iq, 48
.set nb230_dx, 64
.set nb230_dy, 80
.set nb230_dz, 96
.set nb230_c6, 112
.set nb230_c12, 128
.set nb230_tsc, 144
.set nb230_fstmp, 160
.set nb230_vctot, 176
.set nb230_Vvdwtot, 192
.set nb230_fix, 208
.set nb230_fiy, 224
.set nb230_fiz, 240
.set nb230_half, 256
.set nb230_three, 272
.set nb230_two, 288
.set nb230_krf, 304
.set nb230_crf, 320
.set nb230_is3, 336
.set nb230_ii3, 340
.set nb230_ntia, 344
.set nb230_innerjjnr, 348
.set nb230_innerk, 352
.set nb230_n, 356
.set nb230_nn1, 360
.set nb230_nri, 364
.set nb230_facel, 368                         ## uses 8 bytes
.set nb230_ntype, 376
.set nb230_nouter, 380
.set nb230_ninner, 384
.set nb230_salign, 388
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
        movl %eax,nb230_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb230_p_nri(%ebp),%ecx
        movl nb230_p_facel(%ebp),%esi
        movl nb230_p_ntype(%ebp),%edi
        movl (%ecx),%ecx
        movsd (%esi),%xmm7
        movl (%edi),%edi
        movl %ecx,nb230_nri(%esp)
        movsd %xmm7,nb230_facel(%esp)
        movl %edi,nb230_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb230_nouter(%esp)
        movl %eax,nb230_ninner(%esp)

        movl nb230_p_tabscale(%ebp),%eax
        movsd (%eax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb230_tsc(%esp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb230_half(%esp)
        movl %ebx,nb230_half+4(%esp)
        movsd nb230_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb230_half(%esp)
        movapd %xmm2,nb230_two(%esp)
        movapd %xmm3,nb230_three(%esp)

        movl nb230_argkrf(%ebp),%esi
        movl nb230_argcrf(%ebp),%edi
        movsd (%esi),%xmm5
        movsd (%edi),%xmm6
        shufpd $0,%xmm5,%xmm5
        shufpd $0,%xmm6,%xmm6
        movapd %xmm5,nb230_krf(%esp)
        movapd %xmm6,nb230_crf(%esp)

_nb_kernel230_ia32_sse2.nb230_threadloop: 
        movl  nb230_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel230_ia32_sse2.nb230_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel230_ia32_sse2.nb230_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb230_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb230_n(%esp)
        movl %ebx,nb230_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel230_ia32_sse2.nb230_outerstart
        jmp _nb_kernel230_ia32_sse2.nb230_end

_nb_kernel230_ia32_sse2.nb230_outerstart: 
        ## ebx contains number of outer iterations
        addl nb230_nouter(%esp),%ebx
        movl %ebx,nb230_nouter(%esp)

_nb_kernel230_ia32_sse2.nb230_outer: 
        movl  nb230_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb230_is3(%esp)      ## store is3 

        movl  nb230_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb230_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movl  nb230_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        mulsd nb230_facel(%esp),%xmm3
        shufpd $0,%xmm3,%xmm3

        movl  nb230_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb230_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb230_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb230_pos(%ebp),%eax      ## eax = base of pos[]  

        addsd (%eax,%ebx,8),%xmm0
        addsd 8(%eax,%ebx,8),%xmm1
        addsd 16(%eax,%ebx,8),%xmm2

        movapd %xmm3,nb230_iq(%esp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb230_ix(%esp)
        movapd %xmm1,nb230_iy(%esp)
        movapd %xmm2,nb230_iz(%esp)

        movl  %ebx,nb230_ii3(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb230_vctot(%esp)
        movapd %xmm4,nb230_Vvdwtot(%esp)
        movapd %xmm4,nb230_fix(%esp)
        movapd %xmm4,nb230_fiy(%esp)
        movapd %xmm4,nb230_fiz(%esp)

        movl  nb230_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb230_pos(%ebp),%esi
        movl  nb230_faction(%ebp),%edi
        movl  nb230_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb230_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb230_ninner(%esp),%ecx
        movl  %ecx,nb230_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb230_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel230_ia32_sse2.nb230_unroll_loop
        jmp   _nb_kernel230_ia32_sse2.nb230_checksingle
_nb_kernel230_ia32_sse2.nb230_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb230_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        addl $8,nb230_innerjjnr(%esp)                   ## advance pointer (unrolled 2) 

        movl nb230_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3

        movapd nb230_iq(%esp),%xmm5
        mulpd %xmm5,%xmm3               ## qq 

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movl nb230_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb230_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb230_ntia(%esp),%edi
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
        movapd %xmm4,nb230_c6(%esp)
        movapd %xmm6,nb230_c12(%esp)

        movl nb230_pos(%ebp),%esi        ## base of pos[] 

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
        movapd nb230_ix(%esp),%xmm4
        movapd nb230_iy(%esp),%xmm5
        movapd nb230_iz(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb230_dx(%esp)
        movapd %xmm5,nb230_dy(%esp)
        movapd %xmm6,nb230_dz(%esp)
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

        movapd nb230_krf(%esp),%xmm7
        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulpd %xmm2,%xmm2       ## lu*lu 
        movapd nb230_three(%esp),%xmm1
        mulpd %xmm4,%xmm7       ## krsq 
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb230_half(%esp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb230_three(%esp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb230_half(%esp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=rinv 
        movapd %xmm0,%xmm6
        addpd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        movapd %xmm4,%xmm1
        subpd  nb230_crf(%esp),%xmm6
        mulpd  %xmm3,%xmm6
        mulpd  nb230_two(%esp),%xmm7
        movapd %xmm0,%xmm1
        subpd  %xmm7,%xmm1  ## rinv-2*krsq
        mulpd  %xmm0,%xmm1  ## (rinv-2*krsq)*rinv
        mulpd  %xmm3,%xmm1  ## qq*(rinv-2*krsq)*rinv

        movapd %xmm1,nb230_fstmp(%esp)

        addpd  nb230_vctot(%esp),%xmm6
        movapd %xmm6,nb230_vctot(%esp)

        ## LJ table interaction. xmm0=rinv, xmm4=rsq

        mulpd %xmm0,%xmm4       ## xmm4=r 
        mulpd nb230_tsc(%esp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $3,%mm6           ## idx *= 8 

        movd %eax,%mm0
        movd %ebx,%mm1

        movl nb230_VFtab(%ebp),%esi
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
        mulpd  nb230_two(%esp),%xmm7    ## two*Heps2 
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb230_c6(%esp),%xmm4
        mulpd  %xmm4,%xmm7       ## fijD 
        mulpd  %xmm4,%xmm5       ## Vvdw6 

        ## put scalar force on stack Update Vvdwtot directly 
        addpd  nb230_Vvdwtot(%esp),%xmm5
        movapd nb230_fstmp(%esp),%xmm3
        mulpd  nb230_tsc(%esp),%xmm7
        subpd  %xmm7,%xmm3
        movapd %xmm3,nb230_fstmp(%esp)
        movapd %xmm5,nb230_Vvdwtot(%esp)

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
        mulpd  nb230_two(%esp),%xmm7    ## two*Heps2 
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb230_c12(%esp),%xmm4
        mulpd  %xmm4,%xmm7
        mulpd  %xmm4,%xmm5

        addpd  nb230_Vvdwtot(%esp),%xmm5
        movapd nb230_fstmp(%esp),%xmm3
        mulpd  nb230_tsc(%esp),%xmm7
        subpd  %xmm7,%xmm3
        movapd %xmm5,nb230_Vvdwtot(%esp)

        mulpd  %xmm0,%xmm3

        movapd nb230_dx(%esp),%xmm0
        movapd nb230_dy(%esp),%xmm1
        movapd nb230_dz(%esp),%xmm2

        movd %mm0,%eax
        movd %mm1,%ebx

        movl   nb230_faction(%ebp),%edi
        mulpd  %xmm3,%xmm0
        mulpd  %xmm3,%xmm1
        mulpd  %xmm3,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movapd nb230_fix(%esp),%xmm3
        movapd nb230_fiy(%esp),%xmm4
        movapd nb230_fiz(%esp),%xmm5
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm5
        movapd %xmm3,nb230_fix(%esp)
        movapd %xmm4,nb230_fiy(%esp)
        movapd %xmm5,nb230_fiz(%esp)
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
        subl $2,nb230_innerk(%esp)
        jl    _nb_kernel230_ia32_sse2.nb230_checksingle
        jmp   _nb_kernel230_ia32_sse2.nb230_unroll_loop

_nb_kernel230_ia32_sse2.nb230_checksingle:      
        movl  nb230_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel230_ia32_sse2.nb230_dosingle
        jmp    _nb_kernel230_ia32_sse2.nb230_updateouterdata
_nb_kernel230_ia32_sse2.nb230_dosingle: 
        movl nb230_charge(%ebp),%esi
        movl nb230_pos(%ebp),%edi
        movl  nb230_innerjjnr(%esp),%ecx
        xorpd %xmm3,%xmm3
        movl  (%ecx),%eax

        movlpd (%esi,%eax,8),%xmm3
        movapd nb230_iq(%esp),%xmm5
        mulpd %xmm5,%xmm3               ## qq 

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movl nb230_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb230_vdwparam(%ebp),%esi
        shll %eax
        movl nb230_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 
        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movapd %xmm4,nb230_c6(%esp)
        movapd %xmm6,nb230_c12(%esp)

        movl nb230_pos(%ebp),%esi        ## base of pos[] 

        leal (%eax,%eax,2),%eax    ## replace jnr with j3 

        ## move two coordinates to xmm0-xmm2    
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ix-iz to xmm4-xmm6 
        movapd nb230_ix(%esp),%xmm4
        movapd nb230_iy(%esp),%xmm5
        movapd nb230_iz(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb230_dx(%esp)
        movapd %xmm5,nb230_dy(%esp)
        movapd %xmm6,nb230_dz(%esp)
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

        movapd nb230_krf(%esp),%xmm7
        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulsd %xmm2,%xmm2       ## lu*lu 
        movapd nb230_three(%esp),%xmm1
        mulsd %xmm4,%xmm7       ## krsq 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb230_half(%esp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb230_three(%esp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb230_half(%esp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=rinv 
        movapd %xmm0,%xmm6
        addsd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        movapd %xmm4,%xmm1
        subsd  nb230_crf(%esp),%xmm6
        mulsd  %xmm3,%xmm6      ## xmm6=vcoul=qq*(rinv+ krsq) 
        mulsd  nb230_two(%esp),%xmm7
        movapd %xmm0,%xmm1
        subsd  %xmm7,%xmm1  ## rinv-2*krsq
        mulsd  %xmm0,%xmm1  ## (rinv-2*krsq)*rinv
        mulsd  %xmm3,%xmm1
        movsd %xmm1,nb230_fstmp(%esp)

        addsd  nb230_vctot(%esp),%xmm6
        movsd %xmm6,nb230_vctot(%esp)

        ## LJ table interaction. xmm0=rinv, cmm4=rsq

        mulsd %xmm0,%xmm4       ## xmm4=r 
        mulsd nb230_tsc(%esp),%xmm4

        cvttsd2si %xmm4,%ebx    ## mm6 = lu idx 
        cvtsi2sd %ebx,%xmm5
        subsd %xmm5,%xmm4
        movsd %xmm4,%xmm1       ## xmm1=eps 
        movsd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll  $3,%ebx

        movl nb230_VFtab(%ebp),%esi

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
        mulsd  nb230_two(%esp),%xmm7    ## two*Heps2 
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movsd nb230_c6(%esp),%xmm4
        mulsd  %xmm4,%xmm7       ## fijD 
        mulsd  %xmm4,%xmm5       ## Vvdw6 

        ## put scalar force on stack Update Vvdwtot directly 
        addsd  nb230_Vvdwtot(%esp),%xmm5
        movsd nb230_fstmp(%esp),%xmm3
        mulsd  nb230_tsc(%esp),%xmm7
        subsd  %xmm7,%xmm3
        movsd %xmm3,nb230_fstmp(%esp)
        movsd %xmm5,nb230_Vvdwtot(%esp)

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
        mulsd  nb230_two(%esp),%xmm7    ## two*Heps2 
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movsd nb230_c12(%esp),%xmm4
        mulsd  %xmm4,%xmm7
        mulsd  %xmm4,%xmm5

        addsd  nb230_Vvdwtot(%esp),%xmm5
        movsd nb230_fstmp(%esp),%xmm3
        mulsd  nb230_tsc(%esp),%xmm7
        subsd  %xmm7,%xmm3
        movsd %xmm5,nb230_Vvdwtot(%esp)

        mulsd  %xmm0,%xmm3

        movsd nb230_dx(%esp),%xmm0
        movsd nb230_dy(%esp),%xmm1
        movsd nb230_dz(%esp),%xmm2

        movl   nb230_faction(%ebp),%edi
        mulsd  %xmm3,%xmm0
        mulsd  %xmm3,%xmm1
        mulsd  %xmm3,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movlpd nb230_fix(%esp),%xmm3
        movlpd nb230_fiy(%esp),%xmm4
        movlpd nb230_fiz(%esp),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movlpd %xmm3,nb230_fix(%esp)
        movlpd %xmm4,nb230_fiy(%esp)
        movlpd %xmm5,nb230_fiz(%esp)
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

_nb_kernel230_ia32_sse2.nb230_updateouterdata: 
        movl  nb230_ii3(%esp),%ecx
        movl  nb230_faction(%ebp),%edi
        movl  nb230_fshift(%ebp),%esi
        movl  nb230_is3(%esp),%edx

        ## accumulate i forces in xmm0, xmm1, xmm2 
        movapd nb230_fix(%esp),%xmm0
        movapd nb230_fiy(%esp),%xmm1
        movapd nb230_fiz(%esp),%xmm2

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
        movl nb230_n(%esp),%esi
        ## get group index for i particle 
        movl  nb230_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb230_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb230_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb230_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb230_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb230_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel230_ia32_sse2.nb230_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb230_n(%esp)
        jmp _nb_kernel230_ia32_sse2.nb230_outer
_nb_kernel230_ia32_sse2.nb230_outerend: 
        ## check if more outer neighborlists remain
        movl  nb230_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel230_ia32_sse2.nb230_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel230_ia32_sse2.nb230_threadloop
_nb_kernel230_ia32_sse2.nb230_end: 
        emms

        movl nb230_nouter(%esp),%eax
        movl nb230_ninner(%esp),%ebx
        movl nb230_outeriter(%ebp),%ecx
        movl nb230_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb230_salign(%esp),%eax
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


.globl nb_kernel230nf_ia32_sse2
.globl _nb_kernel230nf_ia32_sse2
nb_kernel230nf_ia32_sse2:       
_nb_kernel230nf_ia32_sse2:      
.set nb230nf_p_nri, 8
.set nb230nf_iinr, 12
.set nb230nf_jindex, 16
.set nb230nf_jjnr, 20
.set nb230nf_shift, 24
.set nb230nf_shiftvec, 28
.set nb230nf_fshift, 32
.set nb230nf_gid, 36
.set nb230nf_pos, 40
.set nb230nf_faction, 44
.set nb230nf_charge, 48
.set nb230nf_p_facel, 52
.set nb230nf_argkrf, 56
.set nb230nf_argcrf, 60
.set nb230nf_Vc, 64
.set nb230nf_type, 68
.set nb230nf_p_ntype, 72
.set nb230nf_vdwparam, 76
.set nb230nf_Vvdw, 80
.set nb230nf_p_tabscale, 84
.set nb230nf_VFtab, 88
.set nb230nf_invsqrta, 92
.set nb230nf_dvda, 96
.set nb230nf_p_gbtabscale, 100
.set nb230nf_GBtab, 104
.set nb230nf_p_nthreads, 108
.set nb230nf_count, 112
.set nb230nf_mtx, 116
.set nb230nf_outeriter, 120
.set nb230nf_inneriter, 124
.set nb230nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb230nf_ix, 0
.set nb230nf_iy, 16
.set nb230nf_iz, 32
.set nb230nf_iq, 48
.set nb230nf_dx, 64
.set nb230nf_dy, 80
.set nb230nf_dz, 96
.set nb230nf_c6, 112
.set nb230nf_c12, 128
.set nb230nf_tsc, 144
.set nb230nf_vctot, 176
.set nb230nf_Vvdwtot, 192
.set nb230nf_half, 256
.set nb230nf_three, 272
.set nb230nf_krf, 304
.set nb230nf_crf, 320
.set nb230nf_is3, 336
.set nb230nf_ii3, 340
.set nb230nf_ntia, 344
.set nb230nf_innerjjnr, 348
.set nb230nf_innerk, 352
.set nb230nf_n, 356
.set nb230nf_nn1, 360
.set nb230nf_nri, 364
.set nb230nf_facel, 368                         ## uses 8 bytes
.set nb230nf_ntype, 376
.set nb230nf_nouter, 380
.set nb230nf_ninner, 384
.set nb230nf_salign, 388
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
        movl %eax,nb230nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb230nf_p_nri(%ebp),%ecx
        movl nb230nf_p_facel(%ebp),%esi
        movl nb230nf_p_ntype(%ebp),%edi
        movl (%ecx),%ecx
        movsd (%esi),%xmm7
        movl (%edi),%edi
        movl %ecx,nb230nf_nri(%esp)
        movsd %xmm7,nb230nf_facel(%esp)
        movl %edi,nb230nf_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb230nf_nouter(%esp)
        movl %eax,nb230nf_ninner(%esp)

        movl nb230nf_p_tabscale(%ebp),%eax
        movsd (%eax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb230nf_tsc(%esp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb230nf_half(%esp)
        movl %ebx,nb230nf_half+4(%esp)
        movsd nb230nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb230nf_half(%esp)
        movapd %xmm3,nb230nf_three(%esp)

        movl nb230nf_argkrf(%ebp),%esi
        movl nb230nf_argcrf(%ebp),%edi
        movsd (%esi),%xmm5
        movsd (%edi),%xmm6
        shufpd $0,%xmm5,%xmm5
        shufpd $0,%xmm6,%xmm6
        movapd %xmm5,nb230nf_krf(%esp)
        movapd %xmm6,nb230nf_crf(%esp)

_nb_kernel230nf_ia32_sse2.nb230nf_threadloop: 
        movl  nb230nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel230nf_ia32_sse2.nb230nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel230nf_ia32_sse2.nb230nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb230nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb230nf_n(%esp)
        movl %ebx,nb230nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel230nf_ia32_sse2.nb230nf_outerstart
        jmp _nb_kernel230nf_ia32_sse2.nb230nf_end

_nb_kernel230nf_ia32_sse2.nb230nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb230nf_nouter(%esp),%ebx
        movl %ebx,nb230nf_nouter(%esp)

_nb_kernel230nf_ia32_sse2.nb230nf_outer: 
        movl  nb230nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb230nf_is3(%esp)            ## store is3 

        movl  nb230nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb230nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movl  nb230nf_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        mulsd nb230nf_facel(%esp),%xmm3
        shufpd $0,%xmm3,%xmm3

        movl  nb230nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb230nf_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb230nf_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb230nf_pos(%ebp),%eax      ## eax = base of pos[]  

        addsd (%eax,%ebx,8),%xmm0
        addsd 8(%eax,%ebx,8),%xmm1
        addsd 16(%eax,%ebx,8),%xmm2

        movapd %xmm3,nb230nf_iq(%esp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb230nf_ix(%esp)
        movapd %xmm1,nb230nf_iy(%esp)
        movapd %xmm2,nb230nf_iz(%esp)

        movl  %ebx,nb230nf_ii3(%esp)

        ## clear vctot
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb230nf_vctot(%esp)
        movapd %xmm4,nb230nf_Vvdwtot(%esp)

        movl  nb230nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb230nf_pos(%ebp),%esi
        movl  nb230nf_faction(%ebp),%edi
        movl  nb230nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb230nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb230nf_ninner(%esp),%ecx
        movl  %ecx,nb230nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb230nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel230nf_ia32_sse2.nb230nf_unroll_loop
        jmp   _nb_kernel230nf_ia32_sse2.nb230nf_checksingle
_nb_kernel230nf_ia32_sse2.nb230nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb230nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        addl $8,nb230nf_innerjjnr(%esp)                 ## advance pointer (unrolled 2) 

        movl nb230nf_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3

        movapd nb230nf_iq(%esp),%xmm5
        mulpd %xmm5,%xmm3               ## qq 

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movl nb230nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb230nf_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb230nf_ntia(%esp),%edi
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
        movapd %xmm4,nb230nf_c6(%esp)
        movapd %xmm6,nb230nf_c12(%esp)

        movl nb230nf_pos(%ebp),%esi        ## base of pos[] 

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
        movapd nb230nf_ix(%esp),%xmm4
        movapd nb230nf_iy(%esp),%xmm5
        movapd nb230nf_iz(%esp),%xmm6

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

        movapd nb230nf_krf(%esp),%xmm7
        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulpd %xmm2,%xmm2       ## lu*lu 
        movapd nb230nf_three(%esp),%xmm1
        mulpd %xmm4,%xmm7       ## krsq 
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb230nf_half(%esp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb230nf_three(%esp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb230nf_half(%esp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=rinv 
        movapd %xmm0,%xmm6
        addpd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        movapd %xmm4,%xmm1
        subpd  nb230nf_crf(%esp),%xmm6
        mulpd  %xmm3,%xmm6

        addpd  nb230nf_vctot(%esp),%xmm6
        movapd %xmm6,nb230nf_vctot(%esp)

        ## LJ table interaction. xmm0=rinv, xmm4=rsq

        mulpd %xmm0,%xmm4       ## xmm4=r 
        mulpd nb230nf_tsc(%esp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $3,%mm6           ## idx *= 8 

        movl nb230nf_VFtab(%ebp),%esi
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

        movapd nb230nf_c6(%esp),%xmm4
        mulpd  %xmm4,%xmm5       ## Vvdw6 

        ##Update Vvdwtot directly 
        addpd  nb230nf_Vvdwtot(%esp),%xmm5
        movapd %xmm5,nb230nf_Vvdwtot(%esp)

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

        movapd nb230nf_c12(%esp),%xmm4
        mulpd  %xmm4,%xmm5

        addpd  nb230nf_Vvdwtot(%esp),%xmm5
        movapd %xmm5,nb230nf_Vvdwtot(%esp)

        ## should we do one more iteration? 
        subl $2,nb230nf_innerk(%esp)
        jl    _nb_kernel230nf_ia32_sse2.nb230nf_checksingle
        jmp   _nb_kernel230nf_ia32_sse2.nb230nf_unroll_loop

_nb_kernel230nf_ia32_sse2.nb230nf_checksingle:  
        movl  nb230nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel230nf_ia32_sse2.nb230nf_dosingle
        jmp    _nb_kernel230nf_ia32_sse2.nb230nf_updateouterdata
_nb_kernel230nf_ia32_sse2.nb230nf_dosingle: 
        movl nb230nf_charge(%ebp),%esi
        movl nb230nf_pos(%ebp),%edi
        movl  nb230nf_innerjjnr(%esp),%ecx
        xorpd %xmm3,%xmm3
        movl  (%ecx),%eax

        movlpd (%esi,%eax,8),%xmm3
        movapd nb230nf_iq(%esp),%xmm5
        mulpd %xmm5,%xmm3               ## qq 

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movl nb230nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb230nf_vdwparam(%ebp),%esi
        shll %eax
        movl nb230nf_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 
        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movapd %xmm4,nb230nf_c6(%esp)
        movapd %xmm6,nb230nf_c12(%esp)

        movl nb230nf_pos(%ebp),%esi        ## base of pos[] 

        leal (%eax,%eax,2),%eax    ## replace jnr with j3 

        ## move two coordinates to xmm0-xmm2    
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ix-iz to xmm4-xmm6 
        movapd nb230nf_ix(%esp),%xmm4
        movapd nb230nf_iy(%esp),%xmm5
        movapd nb230nf_iz(%esp),%xmm6

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

        movapd nb230nf_krf(%esp),%xmm7
        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulsd %xmm2,%xmm2       ## lu*lu 
        movapd nb230nf_three(%esp),%xmm1
        mulsd %xmm4,%xmm7       ## krsq 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb230nf_half(%esp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb230nf_three(%esp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb230nf_half(%esp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=rinv 
        movapd %xmm0,%xmm6
        addsd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        movapd %xmm4,%xmm1
        subsd  nb230nf_crf(%esp),%xmm6
        mulsd  %xmm3,%xmm6      ## xmm6=vcoul=qq*(rinv+ krsq) 

        addsd  nb230nf_vctot(%esp),%xmm6
        movsd %xmm6,nb230nf_vctot(%esp)

        ## LJ table interaction. xmm0=rinv, cmm4=rsq

        mulsd %xmm0,%xmm4       ## xmm4=r 
        mulsd nb230nf_tsc(%esp),%xmm4

        cvttsd2si %xmm4,%ebx    ## mm6 = lu idx 
        cvtsi2sd %ebx,%xmm5
        subsd %xmm5,%xmm4
        movsd %xmm4,%xmm1       ## xmm1=eps 
        movsd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll  $3,%ebx

        movl nb230nf_VFtab(%ebp),%esi

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

        movsd nb230nf_c6(%esp),%xmm4
        mulsd  %xmm4,%xmm5       ## Vvdw6 

        ## put scalar force on stack Update Vvdwtot directly 
        addsd  nb230nf_Vvdwtot(%esp),%xmm5
        movsd %xmm5,nb230nf_Vvdwtot(%esp)

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

        movsd nb230nf_c12(%esp),%xmm4
        mulsd  %xmm4,%xmm5

        addsd  nb230nf_Vvdwtot(%esp),%xmm5
        movsd %xmm5,nb230nf_Vvdwtot(%esp)

_nb_kernel230nf_ia32_sse2.nb230nf_updateouterdata: 
        ## get n from stack
        movl nb230nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb230nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb230nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb230nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb230nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb230nf_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb230nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel230nf_ia32_sse2.nb230nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb230nf_n(%esp)
        jmp _nb_kernel230nf_ia32_sse2.nb230nf_outer
_nb_kernel230nf_ia32_sse2.nb230nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb230nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel230nf_ia32_sse2.nb230nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel230nf_ia32_sse2.nb230nf_threadloop
_nb_kernel230nf_ia32_sse2.nb230nf_end: 
        emms

        movl nb230nf_nouter(%esp),%eax
        movl nb230nf_ninner(%esp),%ebx
        movl nb230nf_outeriter(%ebp),%ecx
        movl nb230nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb230nf_salign(%esp),%eax
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


