# mark_description "Intel(R) C Intel(R) 64 Compiler for applications running on Intel(R) 64, Version 18.0.0.128 Build 20170811";
# mark_description "-S -o src/gromacs/hardware/identifyavx512fmaunits_asm.s -fasm-blocks";
	.file "identifyavx512fmaunits_asm.c"
	.text
..TXTST0:
# -- Begin  executeFmaAndShuffleLoop
	.text
# mark_begin;
       .align    16,0x90
	.globl executeFmaAndShuffleLoop
# --- executeFmaAndShuffleLoop(uint64_t)
executeFmaAndShuffleLoop:
# parameter 1: %rdi
..B1.1:                         # Preds ..B1.0
                                # Execution count [1.00e+00]
	.cfi_startproc
..___tag_value_executeFmaAndShuffleLoop.1:
..L2:
                                                          #48.1
        pushq     %rbp                                          #48.1
	.cfi_def_cfa_offset 16
        movq      %rsp, %rbp                                    #48.1
	.cfi_def_cfa 6, 16
	.cfi_offset 6, -16
        andq      $-64, %rsp                                    #48.1
        subq      $192, %rsp                                    #48.1
        movups    one_vec.4.0.0.1(%rip), %xmm0                  #50.45
        movups    16+one_vec.4.0.0.1(%rip), %xmm1               #50.45
        movups    32+one_vec.4.0.0.1(%rip), %xmm2               #50.45
        movups    48+one_vec.4.0.0.1(%rip), %xmm3               #50.45
        movq      %rdi, 128(%rsp)                               #49.20
        movups    %xmm0, 64(%rsp)                               #50.45
        movups    %xmm1, 80(%rsp)                               #50.45
        movups    %xmm2, 96(%rsp)                               #50.45
        movups    %xmm3, 112(%rsp)                              #50.45
                                # LOE rbx r12 r13 r14 r15
..B1.2:                         # Preds ..B1.1
                                # Execution count [5.00e-01]
        movups    shuf_vec.4.0.0.1(%rip), %xmm0                 #51.45
        movups    16+shuf_vec.4.0.0.1(%rip), %xmm1              #51.45
        movups    32+shuf_vec.4.0.0.1(%rip), %xmm2              #51.45
        movups    48+shuf_vec.4.0.0.1(%rip), %xmm3              #51.45
        movups    %xmm0, (%rsp)                                 #51.45
        movups    %xmm1, 16(%rsp)                               #51.45
        movups    %xmm2, 32(%rsp)                               #51.45
        movups    %xmm3, 48(%rsp)                               #51.45
                                # LOE rbx r12 r13 r14 r15
..B1.3:                         # Preds ..B1.2 ..B1.3
                                # Execution count [0.00e+00]
# Begin ASM
        vmovups   64(%rsp), %zmm0                               #54.9
        vmovups   64(%rsp), %zmm1                               #55.9
        vmovups   64(%rsp), %zmm2                               #56.9
        vmovups   64(%rsp), %zmm3                               #57.9
        vmovups   64(%rsp), %zmm4                               #58.9
        vmovups   64(%rsp), %zmm5                               #59.9
        vmovups   64(%rsp), %zmm6                               #60.9
        vmovups   64(%rsp), %zmm7                               #61.9
        vmovups   64(%rsp), %zmm8                               #62.9
        vmovups   64(%rsp), %zmm9                               #63.9
        vmovups   64(%rsp), %zmm10                              #64.9
        vmovups   64(%rsp), %zmm11                              #65.9
        vmovups   (%rsp), %zmm12                                #66.9
        vmovups   (%rsp), %zmm13                                #67.9
        vmovups   (%rsp), %zmm14                                #68.9
        vmovups   (%rsp), %zmm15                                #69.9
        vmovups   (%rsp), %zmm16                                #70.9
        vmovups   (%rsp), %zmm17                                #71.9
        vmovups   (%rsp), %zmm18                                #72.9
        vmovups   (%rsp), %zmm19                                #73.9
        vmovups   (%rsp), %zmm20                                #74.9
        vmovups   (%rsp), %zmm21                                #75.9
        vmovups   (%rsp), %zmm22                                #76.9
        vmovups   (%rsp), %zmm23                                #77.9
        vmovups   (%rsp), %zmm30                                #78.9
        movq      128(%rsp), %rdx                               #79.9
.L_2TAG_PACKET_0.0.1:                                           #80.1
        vfmadd231pd %zmm0, %zmm0, %zmm0                         #81.9
        vfmadd231pd %zmm1, %zmm1, %zmm1                         #82.9
        vfmadd231pd %zmm2, %zmm2, %zmm2                         #83.9
        vfmadd231pd %zmm3, %zmm3, %zmm3                         #84.9
        vfmadd231pd %zmm4, %zmm4, %zmm4                         #85.9
        vfmadd231pd %zmm5, %zmm5, %zmm5                         #86.9
        vfmadd231pd %zmm6, %zmm6, %zmm6                         #87.9
        vfmadd231pd %zmm7, %zmm7, %zmm7                         #88.9
        vfmadd231pd %zmm8, %zmm8, %zmm8                         #89.9
        vfmadd231pd %zmm9, %zmm9, %zmm9                         #90.9
        vfmadd231pd %zmm10, %zmm10, %zmm10                      #91.9
        vfmadd231pd %zmm11, %zmm11, %zmm11                      #92.9
        vpermd    %zmm30, %zmm30, %zmm12                        #93.9
        vpermd    %zmm30, %zmm30, %zmm13                        #94.9
        vpermd    %zmm30, %zmm30, %zmm14                        #95.9
        vpermd    %zmm30, %zmm30, %zmm15                        #96.9
        vpermd    %zmm30, %zmm30, %zmm16                        #97.9
        vpermd    %zmm30, %zmm30, %zmm17                        #98.9
        vpermd    %zmm30, %zmm30, %zmm18                        #99.9
        vpermd    %zmm30, %zmm30, %zmm19                        #100.9
        vpermd    %zmm30, %zmm30, %zmm20                        #101.9
        vpermd    %zmm30, %zmm30, %zmm21                        #102.9
        vpermd    %zmm30, %zmm30, %zmm22                        #103.9
        vpermd    %zmm30, %zmm30, %zmm23                        #104.9
        decq      %rdx                                          #105.9
        jg        .L_2TAG_PACKET_0.0.1                          #106.9
# End ASM
                                # LOE rbx r12 r13 r14 r15
..B1.4:                         # Preds ..B1.3
                                # Execution count [1.00e+00]
        movq      %rbp, %rsp                                    #108.1
        popq      %rbp                                          #108.1
	.cfi_def_cfa 7, 8
	.cfi_restore 6
        ret                                                     #108.1
        .align    16,0x90
                                # LOE
	.cfi_endproc
# mark_end;
	.type	executeFmaAndShuffleLoop,@function
	.size	executeFmaAndShuffleLoop,.-executeFmaAndShuffleLoop
	.section .rodata, "a"
	.align 64
	.align 64
one_vec.4.0.0.1:
	.long	0x00000000,0x3ff00000
	.long	0x00000000,0x3ff00000
	.long	0x00000000,0x3ff00000
	.long	0x00000000,0x3ff00000
	.long	0x00000000,0x3ff00000
	.long	0x00000000,0x3ff00000
	.long	0x00000000,0x3ff00000
	.long	0x00000000,0x3ff00000
	.align 64
shuf_vec.4.0.0.1:
	.long	0
	.long	1
	.long	2
	.long	3
	.long	4
	.long	5
	.long	6
	.long	7
	.long	8
	.long	9
	.long	10
	.long	11
	.long	12
	.long	13
	.long	14
	.long	15
	.data
# -- End  executeFmaAndShuffleLoop
	.text
# -- Begin  executeFmaOnlyLoop
	.text
# mark_begin;
       .align    16,0x90
	.globl executeFmaOnlyLoop
# --- executeFmaOnlyLoop(int)
executeFmaOnlyLoop:
# parameter 1: %edi
..B2.1:                         # Preds ..B2.0
                                # Execution count [5.00e-01]
	.cfi_startproc
..___tag_value_executeFmaOnlyLoop.9:
..L10:
                                                         #111.1
        pushq     %rbp                                          #111.1
	.cfi_def_cfa_offset 16
        movq      %rsp, %rbp                                    #111.1
	.cfi_def_cfa 6, 16
	.cfi_offset 6, -16
        andq      $-64, %rsp                                    #111.1
        subq      $128, %rsp                                    #111.1
        movups    one_vec.5.0.0.2(%rip), %xmm0                  #113.45
        movups    16+one_vec.5.0.0.2(%rip), %xmm1               #113.45
        movups    32+one_vec.5.0.0.2(%rip), %xmm2               #113.45
        movups    48+one_vec.5.0.0.2(%rip), %xmm3               #113.45
        movslq    %edi, %rdi                                    #111.1
        movq      %rdi, 64(%rsp)                                #112.20
        movups    %xmm0, (%rsp)                                 #113.45
        movups    %xmm1, 16(%rsp)                               #113.45
        movups    %xmm2, 32(%rsp)                               #113.45
        movups    %xmm3, 48(%rsp)                               #113.45
                                # LOE rbx r12 r13 r14 r15
..B2.2:                         # Preds ..B2.1 ..B2.2
                                # Execution count [0.00e+00]
# Begin ASM
        vmovups   (%rsp), %zmm0                                 #116.9
        vmovups   (%rsp), %zmm1                                 #117.9
        vmovups   (%rsp), %zmm2                                 #118.9
        vmovups   (%rsp), %zmm3                                 #119.9
        vmovups   (%rsp), %zmm4                                 #120.9
        vmovups   (%rsp), %zmm5                                 #121.9
        vmovups   (%rsp), %zmm6                                 #122.9
        vmovups   (%rsp), %zmm7                                 #123.9
        vmovups   (%rsp), %zmm8                                 #124.9
        vmovups   (%rsp), %zmm9                                 #125.9
        vmovups   (%rsp), %zmm10                                #126.9
        vmovups   (%rsp), %zmm11                                #127.9
        movq      64(%rsp), %rdx                                #128.9
.L_2TAG_PACKET_1.0.2:                                           #129.1
        vfmadd231pd %zmm0, %zmm0, %zmm0                         #130.9
        vfmadd231pd %zmm1, %zmm1, %zmm1                         #131.9
        vfmadd231pd %zmm2, %zmm2, %zmm2                         #132.9
        vfmadd231pd %zmm3, %zmm3, %zmm3                         #133.9
        vfmadd231pd %zmm4, %zmm4, %zmm4                         #134.9
        vfmadd231pd %zmm5, %zmm5, %zmm5                         #135.9
        vfmadd231pd %zmm6, %zmm6, %zmm6                         #136.9
        vfmadd231pd %zmm7, %zmm7, %zmm7                         #137.9
        vfmadd231pd %zmm8, %zmm8, %zmm8                         #138.9
        vfmadd231pd %zmm9, %zmm9, %zmm9                         #139.9
        vfmadd231pd %zmm10, %zmm10, %zmm10                      #140.9
        vfmadd231pd %zmm11, %zmm11, %zmm11                      #141.9
        decq      %rdx                                          #142.9
        jg        .L_2TAG_PACKET_1.0.2                          #143.9
# End ASM
                                # LOE rbx r12 r13 r14 r15
..B2.3:                         # Preds ..B2.2
                                # Execution count [1.00e+00]
        movq      %rbp, %rsp                                    #145.1
        popq      %rbp                                          #145.1
	.cfi_def_cfa 7, 8
	.cfi_restore 6
        ret                                                     #145.1
        .align    16,0x90
                                # LOE
	.cfi_endproc
# mark_end;
	.type	executeFmaOnlyLoop,@function
	.size	executeFmaOnlyLoop,.-executeFmaOnlyLoop
	.section .rodata, "a"
	.align 64
one_vec.5.0.0.2:
	.long	0x00000000,0x3ff00000
	.long	0x00000000,0x3ff00000
	.long	0x00000000,0x3ff00000
	.long	0x00000000,0x3ff00000
	.long	0x00000000,0x3ff00000
	.long	0x00000000,0x3ff00000
	.long	0x00000000,0x3ff00000
	.long	0x00000000,0x3ff00000
	.data
# -- End  executeFmaOnlyLoop
	.data
	.section .note.GNU-stack, ""
// -- Begin DWARF2 SEGMENT .eh_frame
	.section .eh_frame,"a",@progbits
.eh_frame_seg:
	.align 8
# End
