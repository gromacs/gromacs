#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2014,2015,2016,2017,2018, by the GROMACS development team, led by
# Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
# and including many others, as listed in the AUTHORS file in the
# top-level source directory and at http://www.gromacs.org.
#
# GROMACS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# GROMACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with GROMACS; if not, see
# http://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at http://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out http://www.gromacs.org.

# This directory provides a unified place for building all kinds of
# GROMACS documentation. This includes some "static" content (Doxygen
# code documentation, reference manual, install guide, old online HTML
# images), and content generated from the gmx program for the various
# tools (man and HTML pages). It also provides the "webpage" target,
# that combines all of the above (except man pages in man format) into
# a form suitable for automated deployment to the GROMACS website. It
# also provides the INSTALL file for the tarball.
#
# The webpage is mostly built by Sphinx.  Variable values for Sphinx
# substitutions are configured by CMake (for things like version numbers),
# using gmx_configure_version_file().  This happens during build time instead
# of configure time, because 1) some of the version variables are only
# available during build time, and 2) we don't want to do all the Sphinx setup
# during configuration to save some time when not building the content.
# All the generated values get put into conf-vars.py (generated from
# conf-vars.py.cmakein), which in turn is included by the Sphinx configuration
# file conf.py.

set(SOURCE_MD5SUM "unknown" CACHE STRING
    "MD5 sum of the source tarball, normally used only for the pre-release webpage build")
# REGRESSIONTEST_MD5SUM is set in cmake/gmxVersionInfo.cmake because it is used also in tests/CMakeLists.txt
mark_as_advanced(SOURCE_MD5SUM)

set(EXPECTED_DOXYGEN_VERSION 1.8.5)

set(EXPECTED_SPHINX_VERSION 1.6.1)

if (DEFINED PYTHON_EXECUTABLE)
    # Keep quiet on subsequent runs of cmake
    set(PythonInterp_FIND_QUIETLY ON)
endif()
find_package(PythonInterp 2.7)
find_package(Sphinx ${EXPECTED_SPHINX_VERSION} QUIET COMPONENTS pygments)

# Even if we aren't going to make the full webpage, set up to put all
# the documentation output in the same place, for convenience
set(HTML_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/html")
file(MAKE_DIRECTORY ${HTML_OUTPUT_DIR})

# Prepare directories for pdf/tex output
set(TEX_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/manual")
set(LATEX_INPUT_FILES_DIR "${TEX_OUTPUT_DIR}")
file(MAKE_DIRECTORY ${TEX_OUTPUT_DIR})

# The directory from which man pages will be installed; if it remains
# empty, they will be silently skipped.
set(MAN_PAGE_DIR)
if (SOURCE_IS_SOURCE_DISTRIBUTION)
    # When building from the tarball, install the bundled man pages
    # (unless overridden).
    set(MAN_PAGE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
endif()

add_subdirectory(doxygen)

if (SPHINX_FOUND)
    # We need to have all the Sphinx input files in a single directory, and
    # since some of them are generated, we copy everything into the build tree,
    # to this directory.
    set(SPHINX_INPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/sphinx-input)
    # need to move manual part here because we anyhow need the sphinx information for
    # everything else
    add_subdirectory(manual)
    set(SPHINX_EXTENSION_PATH ${CMAKE_CURRENT_SOURCE_DIR})
    if (SOURCE_MD5SUM STREQUAL "unknown")
        # But for testing the webpage build (e.g. from the repo) we
        # need a default value.
        set(REGRESSIONTEST_MD5SUM_STRING "unknown")
    else()
        # The real build of the webpage happens from the tarball, and
        # this should be set to the matching MD5 sum.
        set(REGRESSIONTEST_MD5SUM_STRING "${REGRESSIONTEST_MD5SUM}")
    endif()
    set(WEBMANUAL_SPHINX_FILES
        web-manual/algorithms.rst
        web-manual/analysis.rst
        web-manual/averages.rst
        web-manual/definitions.rst
        web-manual/details.rst
        web-manual/file-formats.rst
        web-manual/functions.rst
        web-manual/index.rst
        web-manual/introduction.rst
        web-manual/preface.rst
        web-manual/references.rst
        web-manual/run-parameters.rst
        web-manual/special.rst
        web-manual/topologies.rst
        web-manual/plots/angle.eps
        web-manual/plots/awh-invN.eps
        web-manual/plots/awh-pmfs.eps
        web-manual/plots/awh-sampleweights.eps
        web-manual/plots/awh-traj.eps
        web-manual/plots/bstretch.eps
        web-manual/plots/chain.eps
        web-manual/plots/compelsetup.pdf
        web-manual/plots/dd-cells.eps
        web-manual/plots/dd-tric.eps
        web-manual/plots/decomp.eps
        web-manual/plots/dih-def.eps
        web-manual/plots/dih.eps
        web-manual/plots/distm.eps
        web-manual/plots/drift-all.eps
        web-manual/plots/dssp.eps
        web-manual/plots/dumaro.eps
        web-manual/plots/dummies.eps
        web-manual/plots/dumtypes.eps
        web-manual/plots/equipotential.pdf
        web-manual/plots/f-angle.eps
        web-manual/plots/f-bham.eps
        web-manual/plots/f-bond.eps
        web-manual/plots/fbposres.eps
        web-manual/plots/f-dih.eps
        web-manual/plots/f-dr.eps
        web-manual/plots/field.eps
        web-manual/plots/fig-02.eps
        web-manual/plots/fig-04.eps
        web-manual/plots/f-lj.eps
        web-manual/plots/flowchart.eps
        web-manual/plots/f-morse.eps
        web-manual/plots/fp-highres.eps
        web-manual/plots/f-imps.pdf
        web-manual/plots/f-pr.eps
        web-manual/plots/f-rbs.eps
        web-manual/plots/free1.eps
        web-manual/plots/free2.eps
        web-manual/plots/gaussians.pdf
        web-manual/plots/hbond.eps
        web-manual/plots/hbond-insert.eps
        web-manual/plots/hpr-wheel.eps
        web-manual/plots/int-mat.eps
        web-manual/plots/leapfrog.eps
        web-manual/plots/lincs.eps
        web-manual/plots/maxwell.eps
        web-manual/plots/mdpar.eps
        web-manual/plots/mpmd-pme.eps
        web-manual/plots/msdwater.eps
        web-manual/plots/ngmxdump.eps
        web-manual/plots/nstric.eps
        web-manual/plots/par-lincs2.eps
        web-manual/plots/parsort.eps
        web-manual/plots/pbctric.eps
        web-manual/plots/peregrine.png
        web-manual/plots/phipsi.eps
        web-manual/plots/pulldirrel.eps
        web-manual/plots/pull.eps
        web-manual/plots/pullref.eps
        web-manual/plots/rama.eps
        web-manual/plots/rdf.eps
        web-manual/plots/rdfO-O.eps
        web-manual/plots/rhododec.eps
        web-manual/plots/ring.eps
        web-manual/plots/ring-imp.eps
        web-manual/plots/rotation.pdf
        web-manual/plots/sgangle.pdf
        web-manual/plots/shiftf.eps
        web-manual/plots/softcore.eps
        web-manual/plots/subst-im.eps
        web-manual/plots/tetra-im.eps
        web-manual/plots/truncoct.eps
        web-manual/plots/vcrf.eps
        web-manual/plots/verlet-drift.eps
        web-manual/plots/vsite-4fdn.eps
        web-manual/plots/adress.png
        web-manual/plots/plotje.png
        web-manual/plots/xvgr.png
        )
    set(SPHINX_SOURCE_FILES
        index.rst
        download.rst
        conf.py
        links.dat
        dev-manual/build-system.rst
        dev-manual/commitstyle.rst
        dev-manual/documentation-generation.rst
        dev-manual/contribute.rst
        dev-manual/doxygen.rst
        dev-manual/error-handling.rst
        dev-manual/formatting.rst
        dev-manual/gmxtree.rst
        dev-manual/includestyle.rst
        dev-manual/index.rst
        dev-manual/jenkins.rst
        dev-manual/language-features.rst
        dev-manual/naming.rst
        dev-manual/overview.rst
        dev-manual/redmine-states.png
        dev-manual/relocatable-binaries.rst
        dev-manual/reportstyle.rst
        dev-manual/style.rst
        dev-manual/testutils.rst
        dev-manual/tools.rst
        dev-manual/uncrustify.rst
        fragments/doxygen-links.rst
        install-guide/index.rst
        release-notes/index.rst
        release-notes/2018/2018.1.rst
        release-notes/2018/major/highlights.rst
        release-notes/2018/major/features.rst
        release-notes/2018/major/performance.rst
        release-notes/2018/major/tools.rst
        release-notes/2018/major/bugs-fixed.rst
        release-notes/2018/major/removed-features.rst
        release-notes/2018/major/portability.rst
        release-notes/2018/major/miscellaneous.rst
        release-notes/2016/2016.5.rst
        release-notes/2016/2016.4.rst
        release-notes/2016/2016.3.rst
        release-notes/2016/2016.2.rst
        release-notes/2016/2016.1.rst
        release-notes/2016/major/highlights.rst
        release-notes/2016/major/new-features.rst
        release-notes/2016/major/performance.rst
        release-notes/2016/major/tools.rst
        release-notes/2016/major/bugs-fixed.rst
        release-notes/2016/major/removed-features.rst
        release-notes/2016/major/miscellaneous.rst
        release-notes/older/index.rst
        user-guide/index.rst
        user-guide/cmdline.rst
        user-guide/cutoff-schemes.rst
        user-guide/environment-variables.rst
        user-guide/faq.rst
        user-guide/floating-point.rst
        user-guide/system-preparation.rst
        user-guide/force-fields.rst
        user-guide/getting-started.rst
        user-guide/managing-simulations.rst
        user-guide/mdp-options.rst
        user-guide/mdrun-features.rst
        user-guide/mdrun-performance.rst
        user-guide/run-time-errors.rst
        user-guide/system-preparation.rst
        user-guide/terminology.rst
        )

    include(SphinxMacros.cmake)
    gmx_init_sphinx_setup(${SPHINX_INPUT_DIR})

    set(SPHINX_CONFIG_VARS_FILE ${SPHINX_INPUT_DIR}/conf-vars.py)
    gmx_configure_version_file(conf-vars.py.cmakein ${SPHINX_CONFIG_VARS_FILE}
        EXTRA_VARS
            SPHINX_EXTENSION_PATH RELENG_PATH
            EXPECTED_DOXYGEN_VERSION
            EXPECTED_SPHINX_VERSION
            CMAKE_MINIMUM_REQUIRED_VERSION REQUIRED_CUDA_VERSION
            REQUIRED_OPENCL_MIN_VERSION
            REQUIRED_CUDA_COMPUTE_CAPABILITY REGRESSIONTEST_VERSION
            SOURCE_MD5SUM REGRESSIONTEST_MD5SUM_STRING
            GMX_TNG_MINIMUM_REQUIRED_VERSION
            GMX_LMFIT_MINIMUM_REQUIRED_VERSION
        COMMENT "Configuring Sphinx configuration file")
    gmx_add_sphinx_input_file(${SPHINX_CONFIG_VARS_FILE})
    gmx_add_sphinx_source_files(FILES ${SPHINX_SOURCE_FILES})
    if (EXISTS ${RELENG_PATH}/docs/FileList.cmake)
        include(${RELENG_PATH}/docs/FileList.cmake)
        gmx_add_sphinx_source_files(
            FROM ${RELENG_PATH}/docs TO dev-manual/releng PREFIX releng/docs/
            FILES ${RELENG_SPHINX_FILES})
    else()
        gmx_add_sphinx_source_files(FILES
            dev-manual/releng/index.rst
            dev-manual/releng/jenkins-howto.rst
            dev-manual/releng/jenkins-ui.rst
            )
    endif()
    if (IMAGE_CONVERT_POSSIBLE)
        gmx_add_sphinx_source_files(
            FILES
            ${WEBMANUAL_SPHINX_FILES}
            user-guide/index.rst
            user-guide/flow.rst)
        gmx_add_sphinx_source_files(
            FROM ${CMAKE_CURRENT_SOURCE_DIR}/manual TO ../manual
            FILES
            picins.sty
            fancy.sty
            here.sty
            underscore.sty)
    else()
        gmx_add_sphinx_source_files(
            FROM ${CMAKE_CURRENT_SOURCE_DIR}/web-manual TO user-guide
            FILES
            file-formats.rst
            plots/plotje.png
            plots/xvgr.png)
        gmx_add_sphinx_source_files(
            FROM ${CMAKE_CURRENT_SOURCE_DIR}/web-manual-dummy TO web-manual
            FILES
            index.rst)
        gmx_add_sphinx_source_files(
            FROM ${CMAKE_CURRENT_SOURCE_DIR}/user-guide-old TO user-guide
            FILES
            index.rst
            flow.rst)
    endif()
    gmx_add_sphinx_input_target(sphinx-input)
    # Remove other rst files from the build tree, since they confuse Sphinx.
    # Skip generated files in onlinehelp/, and fragments.
    # The latter do not cause issues with obsolete files, as they
    # are not considered as Sphinx input files, but will only be
    # included using an explicit .. include::.
    gmx_remove_obsolete_sphinx_input_files("^(onlinehelp|fragments)/.*\\\\.rst$")

    # TODO: Make this remove obsolete .rst files.
    # TODO: This does not work in cross-compilation scenarios; disable up to
    # the necessary level.
    gmx_add_custom_output_target(sphinx-programs OUTPUT STAMP
        COMMAND ${CMAKE_COMMAND} -E make_directory onlinehelp
        COMMAND gmx -quiet help -export rst
        DEPENDS gmx
        WORKING_DIRECTORY ${SPHINX_INPUT_DIR}
        COMMENT "Generating reStructuredText help")
    # This dependency ensures that the directories exist before the
    # executable tries to write things there.
    add_dependencies(sphinx-programs sphinx-input)

    # Make the INSTALL file for CPack for the tarball. This gets put
    # into the tarball via the CPack rules below, which requires that
    # the INSTALL file is in a separate directory by itself.
    set(TEXT_INSTALL_GUIDE_OUTPUT_DIR "install-guide/text")
    add_custom_target(install-guide
        COMMAND
            ${SPHINX_EXECUTABLE}
            -q -E -b text
            -w sphinx-install.log
            -d ${CMAKE_CURRENT_BINARY_DIR}/install-guide/_doctrees
            -c ${SPHINX_INPUT_DIR}
            "${SPHINX_INPUT_DIR}/install-guide"
            "${TEXT_INSTALL_GUIDE_OUTPUT_DIR}"
        COMMAND
            ${CMAKE_COMMAND} -E rename
            ${TEXT_INSTALL_GUIDE_OUTPUT_DIR}/index.txt
            ${TEXT_INSTALL_GUIDE_OUTPUT_DIR}/INSTALL
        WORKING_DIRECTORY
            ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Building INSTALL with Sphinx"
        VERBATIM
        )
    add_dependencies(install-guide sphinx-input)
    gmx_cpack_add_generated_source_directory(install-guide/text DESTINATION /)

    # Sphinx cache with pickled ReST documents
    set(SPHINX_CACHE_DIR "${CMAKE_CURRENT_BINARY_DIR}/_doctrees")
    add_custom_target(webpage-sphinx
        DEPENDS sphinx-programs
        DEPENDS sphinx-input
        DEPENDS sphinx-image-conv-png
        COMMAND
            ${CMAKE_COMMAND} -E make_directory ${SPHINX_INPUT_DIR}/_static
        COMMAND
            ${SPHINX_EXECUTABLE}
            -q -E -b html
            -w sphinx-html.log
            -d "${SPHINX_CACHE_DIR}"
            "${SPHINX_INPUT_DIR}"
            "${HTML_OUTPUT_DIR}"
        WORKING_DIRECTORY
            ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Building HTML documentation with Sphinx"
        VERBATIM
        )
    add_dependencies(webpage-sphinx sphinx-input sphinx-programs)

    add_custom_target(man
        COMMAND
            ${SPHINX_EXECUTABLE}
            -q -E -b man
            -w sphinx-man.log
            -d ${SPHINX_CACHE_DIR}
            -t do_man
            ${SPHINX_INPUT_DIR}
            ${CMAKE_CURRENT_BINARY_DIR}/man
        COMMENT "Building man pages with Sphinx"
        VERBATIM)
    add_dependencies(man sphinx-input sphinx-programs)
    if (GMX_BUILD_HELP)
        # If requested, install the man pages built by the 'man' target
        # created above.  Nothing will be installed if the user did not
        # manually build the target.
        set(MAN_PAGE_DIR ${CMAKE_CURRENT_BINARY_DIR})
    endif()

    add_custom_target(sphinx-create-texman
        DEPENDS sphinx-programs
        DEPENDS sphinx-input
        DEPENDS sphinx-image-conv-pdf
        COMMAND
            ${SPHINX_EXECUTABLE}
            -q -E -b latex
            -w sphinx-latex.log
            -d "${SPHINX_CACHE_DIR}"
            "${SPHINX_INPUT_DIR}"
            "${TEX_OUTPUT_DIR}"
        WORKING_DIRECTORY
            ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Preparing PDF manual LaTeX input files"
        VERBATIM
        )
    add_dependencies(sphinx-create-texman sphinx-input sphinx-programs)

else()
    set(MANUAL_BUILD_IS_POSSIBLE OFF)
    set(MANUAL_BUILD_NOT_POSSIBLE_REASON "Sphinx version ${EXPECTED_SPHINX_VERSION} is not available")

    add_custom_target(webpage-sphinx
        COMMAND ${CMAKE_COMMAND} -E echo
            "HTML pages cannot be built because Sphinx version ${EXPECTED_SPHINX_VERSION} is not available"
        VERBATIM)
    add_custom_target(install-guide
        COMMAND ${CMAKE_COMMAND} -E echo
            "INSTALL cannot be built because Sphinx version ${EXPECTED_SPHINX_VERSION} is not available"
        VERBATIM)
    add_custom_target(man
        COMMAND ${CMAKE_COMMAND} -E echo
            "man pages cannot be built because Sphinx version ${EXPECTED_SPHINX_VERSION} is not available"
        VERBATIM)
    add_custom_target(sphinx-create-texman
        COMMAND ${CMAKE_COMMAND} -E echo
            "Cannot prepare LaTeX input files because Sphinx version ${EXPECTED_SPHINX_VERSION} is not available"
        VERBATIM)
endif()

if (MAN_PAGE_DIR)
    set(MAN_PAGE_DIR ${MAN_PAGE_DIR}/man)
    # Trailing slash on directory is significant for
    # install(DIRECTORY). See CMake docs.
    install(DIRECTORY ${MAN_PAGE_DIR}/
        DESTINATION ${MAN_INSTALL_DIR}/man1
        COMPONENT man OPTIONAL
        FILES_MATCHING PATTERN "*.1")
endif()
gmx_cpack_add_generated_source_directory(man)

# Determine whether we can build all the HTML pages and content linked from
# there.  If not, construct an informative message if the user tries to
# build the target; most people never need to know, unless they've asked for
# the webpage build.
set(HTML_BUILD_IS_POSSIBLE ON)
set(HTML_BUILD_NOT_POSSIBLE_REASON)
set(HTML_BUILD_WARNINGS)

# Next, turn it off if any of the preconditions are unsatisified
if (NOT PYTHON_EXECUTABLE)
    set(HTML_BUILD_IS_POSSIBLE OFF)
    set(HTML_BUILD_NOT_POSSIBLE_REASON "Python is required")
elseif (NOT SPHINX_FOUND)
    # Hardly anything gets built if Sphinx is not available, so don't bother.
    set(HTML_BUILD_IS_POSSIBLE OFF)
    set(HTML_BUILD_NOT_POSSIBLE_REASON "Sphinx version ${EXPECTED_SPHINX_VERSION} is required")
endif()
if (NOT MANUAL_BUILD_IS_POSSIBLE)
    list(APPEND HTML_BUILD_WARNINGS
         "Reference PDF manual was not built, so links to it do not work")
endif()
if (NOT DOXYGEN_EXECUTABLE)
    list(APPEND HTML_BUILD_WARNINGS
        "Doxygen was not available, so links to Doxygen do not work")
endif()
if (NOT DOXYGEN_DOT_EXECUTABLE)
    list(APPEND HTML_BUILD_WARNINGS
        "dot/graphviz was not found, so some graphs are missing")
endif()

if (HTML_BUILD_IS_POSSIBLE)
    set(_webpage_target_properties)
    if (HTML_BUILD_WARNINGS)
        list(APPEND _webpage_target_properties
             COMMAND ${CMAKE_COMMAND} -E echo
                 "webpage was built, but with the following limitations:")
        foreach(_warning ${HTML_BUILD_WARNINGS})
        list(APPEND _webpage_target_properties
             COMMAND ${CMAKE_COMMAND} -E echo " - ${_warning}")
        endforeach()
    endif()

    if (MANUAL_BUILD_IS_POSSIBLE)
        # Make the PDF reference guide
        # TODO Try to make the PDF arrive directly in ${HTML_OUTPUT_DIR}
        # TODO Make this depend on the output of the manual build, so that the
        # file actually gets copied multiple times.
        set(_manual_target_location ${HTML_OUTPUT_DIR}/manual-${GMX_VERSION_STRING}.pdf)
        add_custom_command(
            OUTPUT ${_manual_target_location}
            COMMAND ${CMAKE_COMMAND}
                -E remove -f ${_manual_target_location}
            COMMAND ${CMAKE_COMMAND}
                -E copy ${CMAKE_CURRENT_BINARY_DIR}/manual/Gromacs.pdf ${_manual_target_location}
            DEPENDS manual
            VERBATIM)
        list(APPEND _webpage_target_properties
             DEPENDS ${_manual_target_location})
    endif()

    # The Doxygen configuration in doxygen/Doxyfile-common.cmakein
    # makes all the Doxygen output directly in
    # ${HTML_OUTPUT_DIR}/doxygen (and makes the directory if it needs
    # to).

    # Add a top-level target that builds everything related to the webpage,
    # for Jenkins (and possibly others) to use
    add_custom_target(webpage ${_webpage_target_properties}
        COMMENT "Building webpage"
        VERBATIM)
    add_dependencies(webpage webpage-sphinx doxygen-all)
else()
    add_custom_target(webpage
        COMMAND ${CMAKE_COMMAND} -E echo
            "Cannot build webpage because ${HTML_BUILD_NOT_POSSIBLE_REASON}"
        COMMENT "Webpage build not possible"
        VERBATIM)
endif()
