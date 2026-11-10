# Vendored Zstd, unmodified. Do not edit lib/; a re-vendor replaces it wholesale.
#
#   Upstream: https://github.com/facebook/zstd/releases/download/v1.5.6/zstd-1.5.6.tar.gz
#   SHA256:   8c29e06cf42aacc1eafc4077ae2ec6c6fcb96a626157e0593d5e82a34fd403c1
#   License:  LICENSE beside this file (BSD-3-Clause; upstream dual-licenses GPLv2)
#
# Not 1.5.7: it made level 3, the one Lossless_zstd uses, 26-30% slower for under 0.6% of size.
# Its correctness fix needs a ZSTD_CCtx reused across calls; SZ3 makes one per ZSTD_compress.

set(_zstd_lib "${CMAKE_CURRENT_SOURCE_DIR}/lib")

# Listed rather than globbed: a glob would silently start compiling whatever a re-vendor drops in.
add_library(sz3_zstd STATIC
        "${_zstd_lib}/common/debug.c"
        "${_zstd_lib}/common/entropy_common.c"
        "${_zstd_lib}/common/error_private.c"
        "${_zstd_lib}/common/fse_decompress.c"
        "${_zstd_lib}/common/pool.c"
        "${_zstd_lib}/common/threading.c"
        "${_zstd_lib}/common/xxhash.c"
        "${_zstd_lib}/common/zstd_common.c"

        "${_zstd_lib}/compress/fse_compress.c"
        "${_zstd_lib}/compress/hist.c"
        "${_zstd_lib}/compress/huf_compress.c"
        "${_zstd_lib}/compress/zstd_compress.c"
        "${_zstd_lib}/compress/zstd_compress_literals.c"
        "${_zstd_lib}/compress/zstd_compress_sequences.c"
        "${_zstd_lib}/compress/zstd_compress_superblock.c"
        "${_zstd_lib}/compress/zstd_double_fast.c"
        "${_zstd_lib}/compress/zstd_fast.c"
        "${_zstd_lib}/compress/zstd_lazy.c"
        "${_zstd_lib}/compress/zstd_ldm.c"
        "${_zstd_lib}/compress/zstd_opt.c"
        "${_zstd_lib}/compress/zstdmt_compress.c"

        "${_zstd_lib}/decompress/huf_decompress.c"
        "${_zstd_lib}/decompress/zstd_ddict.c"
        "${_zstd_lib}/decompress/zstd_decompress.c"
        "${_zstd_lib}/decompress/zstd_decompress_block.c"
)

# sz3_zstd, not zstd: a project that adds SZ3 as a subdirectory may have a zstd target of its own,
# and libzstd would overwrite the one already in the install prefix.
# Hidden visibility keeps ZSTD_* out of libhdf5sz3.so, where it interposes on the consumer's Zstd.
set_target_properties(sz3_zstd PROPERTIES
        EXPORT_NAME zstd
        POSITION_INDEPENDENT_CODE ON
        C_VISIBILITY_PRESET hidden)

# ZSTDLIB_VISIBILITY=: Zstd's headers mark its API visibility("default"), past C_VISIBILITY_PRESET.
# ZSTD_DISABLE_ASM=1: huf_decompress_amd64.S is not in the source list above.
# ZSTD_LEGACY_SUPPORT=0: lib/legacy is not vendored.
target_compile_definitions(sz3_zstd
        PRIVATE ZSTDLIB_VISIBILITY= ZSTDERRORLIB_VISIBILITY= ZSTD_DISABLE_ASM=1 ZSTD_LEGACY_SUPPORT=0)

# PRIVATE only: a public include directory would put this zstd.h ahead of a consumer's own.
target_include_directories(sz3_zstd
        PRIVATE
        "${_zstd_lib}"
        "${_zstd_lib}/common"
        "${_zstd_lib}/compress"
        "${_zstd_lib}/decompress"
)

if (SZ3_INSTALL)
    install(TARGETS sz3_zstd EXPORT SZ3Targets
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif ()

export(TARGETS sz3_zstd FILE zstd.cmake)
