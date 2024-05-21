#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2014- The GROMACS Authors
# and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
# Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
# https://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at https://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out https://www.gromacs.org.

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
# All of the generated values get put into conf.py (generated from
# conf.cmakein.py).

set(SOURCE_MD5SUM "unknown" CACHE STRING
    "MD5 sum of the source tarball, normally used only for the pre-release webpage build")
# REGRESSIONTEST_MD5SUM is set in cmake/gmxVersionInfo.cmake because it is used also in tests/CMakeLists.txt
mark_as_advanced(SOURCE_MD5SUM)

set(EXPECTED_DOXYGEN_VERSION 1.8.5)

set(EXPECTED_SPHINX_VERSION 4.0.0)

# By default, suppress output after first configuration.
if(SPHINX_ALREADY_SEARCHED)
    set(Sphinx_FIND_QUIETLY ON)
endif()
find_package(Sphinx ${EXPECTED_SPHINX_VERSION} COMPONENTS pygments)
set(SPHINX_ALREADY_SEARCHED TRUE CACHE BOOL "True if a search for Sphinx has already been done")
mark_as_advanced(SPHINX_ALREADY_SEARCHED)

# Even if we aren't going to make the full webpage, set up to put all
# the documentation output in the same place, for convenience
set(HTML_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/html")
file(MAKE_DIRECTORY ${HTML_OUTPUT_DIR})

# Prepare directories for pdf/tex output
set(TEX_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/manual")
set(SPHINX_LATEX_FILE "${TEX_OUTPUT_DIR}/gromacs.tex")
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

option(GMX_BUILD_MANUAL "Whether to try to configure to build the PDF manual" OFF)
mark_as_advanced(GMX_BUILD_MANUAL)
gmx_add_cache_dependency(GMX_BUILD_MANUAL BOOL GMX_DEVELOPER_BUILD ON)

set(MANUAL_BUILD_IS_POSSIBLE ON)
set(MANUAL_BUILD_NOT_POSSIBLE_REASON)
set(NO_IMAGE_CONVERT_REASON)
if (NOT GMX_BUILD_MANUAL)
    # Make sure we only do detection of manual-building dependencies
    # when the user opted in for that.
    set(MANUAL_BUILD_IS_POSSIBLE OFF)
    set(MANUAL_BUILD_NOT_POSSIBLE_REASON "GMX_BUILD_MANUAL is not ON")
elseif (BUILD_IS_INSOURCE)
    # UseLATEX.cmake does not work with an in-source build
    set(MANUAL_BUILD_IS_POSSIBLE OFF)
    set(MANUAL_BUILD_NOT_POSSIBLE_REASON "the build is in-source")
else()
    include(manual/UseLATEX.cmake)
    if(NOT SPHINX_FOUND)
        set(MANUAL_BUILD_IS_POSSIBLE OFF)
        set(MANUAL_BUILD_NOT_POSSIBLE_REASON "Sphinx has not been found and is needed to create the LaTex input files")
    elseif(NOT PDFLATEX_COMPILER)
        set(MANUAL_BUILD_IS_POSSIBLE OFF)
        set(MANUAL_BUILD_NOT_POSSIBLE_REASON "pdflatex is not available")
    elseif(NOT IMAGE_CONVERT_POSSIBLE)
        set(MANUAL_BUILD_IS_POSSIBLE OFF)
        set(MANUAL_BUILD_NOT_POSSIBLE_REASON "a working form of ImageMagick convert is not available")
        set(NO_IMAGE_CONVERT_REASON "Can not convert files for online or pdf manual")
        # TODO Later, identify other dependencies like bibtex,
        # make_index, date, some graphics conversion program,
        # etc. Perhaps patch UseLATEX.cmake and contribute upstream.
    endif()

endif()

if (SPHINX_FOUND)
    # We need to have all the Sphinx input files in a single directory, and
    # since some of them are generated, we copy everything into the build tree,
    # to this directory.
    set(SPHINX_INPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/sphinx-input)
    set(SPHINX_EXTENSION_PATH ${CMAKE_CURRENT_SOURCE_DIR})
    # As the manual build now depends also on Sphinx, the inclusion path needs
    # to be set after we know the basic information for Sphinx.
    add_subdirectory(manual)
    if (SOURCE_MD5SUM STREQUAL "unknown")
        # But for testing the webpage build (e.g. from the repo) we
        # need a default value.
        set(REGRESSIONTEST_MD5SUM_STRING "unknown")
    else()
        # The real build of the webpage happens from the tarball, and
        # this should be set to the matching MD5 sum.
        set(REGRESSIONTEST_MD5SUM_STRING "${REGRESSIONTEST_MD5SUM}")
    endif()
    # The reference manual still contains the individual sections from the
    # LaTeX document, with the files below grouped and annotated by chapter.
    set(REFERENCEMANUAL_SPHINX_FILES_GENERAL
        # Main index file, preface and introduction.
        reference-manual/index.rst
        reference-manual/preface.rst
        reference-manual/introduction.rst
        # Definitions and Units chapter
        reference-manual/definitions.rst
        # Topologies chapter
        reference-manual/topologies/topologies.rst
        reference-manual/topologies/particle-type.rst
        reference-manual/topologies/parameter-files.rst
        reference-manual/topologies/molecule-definition.rst
        reference-manual/topologies/constraint-algorithm-section.rst
        reference-manual/topologies/pdb2gmx-input-files.rst
        reference-manual/topologies/topology-file-formats.rst
        reference-manual/topologies/force-field-organization.rst
        # File formats chapter
        reference-manual/file-formats.rst
        # Run parameters and programs chapter
        reference-manual/run-parameters.rst
        # Implementation details chapter
        reference-manual/details.rst
        # Averages and fluctations chapter
        reference-manual/averages.rst
        # References
        reference-manual/references.rst
        # PNG formated plot files that don't need to be converted into PNG
        # for the web page.
        reference-manual/plots/GMX_logos/gmx_falcon_blue.png
        reference-manual/plots/plotje.png
        reference-manual/plots/xvgr.png
        )
    set(REFERENCEMANUAL_SPHINX_FILES_WITH_IMAGES
        # Algorithms chapter
        reference-manual/algorithms/algorithms.rst
        reference-manual/algorithms/periodic-boundary-conditions.rst
        reference-manual/algorithms/group-concept.rst
        reference-manual/algorithms/molecular-dynamics.rst
        reference-manual/algorithms/shell-molecular-dynamics.rst
        reference-manual/algorithms/constraint-algorithms.rst
        reference-manual/algorithms/simulated-annealing.rst
        reference-manual/algorithms/stochastic-dynamics.rst
        reference-manual/algorithms/brownian-dynamics.rst
        reference-manual/algorithms/energy-minimization.rst
        reference-manual/algorithms/normal-mode-analysis.rst
        reference-manual/algorithms/free-energy-calculations.rst
        reference-manual/algorithms/replica-exchange.rst
        reference-manual/algorithms/essential-dynamics.rst
        reference-manual/algorithms/expanded-ensemble.rst
        reference-manual/algorithms/parallelization-domain-decomp.rst
        # Interaction functions and force fields chapter
        reference-manual/functions/functions.rst
        reference-manual/functions/bonded-interactions.rst
        reference-manual/functions/force-field.rst
        reference-manual/functions/free-energy-interactions.rst
        reference-manual/functions/interaction-methods.rst
        reference-manual/functions/long-range-electrostatics.rst
        reference-manual/functions/long-range-vdw.rst
        reference-manual/functions/nonbonded-interactions.rst
        reference-manual/functions/polarization.rst
        reference-manual/functions/restraints.rst
        # Special topics chapter
        reference-manual/special/special.rst
        reference-manual/special/free-energy-implementation.rst
        reference-manual/special/pulling.rst
        reference-manual/special/awh.rst
        reference-manual/special/enforced-rotation.rst
        reference-manual/special/electric-fields.rst
        reference-manual/special/comp-electrophys.rst
        reference-manual/special/free-energy-pmf.rst
        reference-manual/special/remove-fast-dgf.rst
        reference-manual/special/viscosity-calculation.rst
        reference-manual/special/shearing.rst
        reference-manual/special/tabulated-interaction-functions.rst
        reference-manual/special/qmmm.rst
        reference-manual/special/vmd-imd.rst
        reference-manual/special/membrane-embedding.rst
        reference-manual/special/mimic-qmmm.rst
        reference-manual/special/density-guided-simulation.rst
        reference-manual/special/colvars.rst
        # Analysis chapter
        reference-manual/analysis/analysis.rst
        reference-manual/analysis/using-groups.rst
        reference-manual/analysis/looking-at-trajectory.rst
        reference-manual/analysis/general-properties.rst
        reference-manual/analysis/radial-distribution-function.rst
        reference-manual/analysis/correlation-function.rst
        reference-manual/analysis/curve-fitting.rst
        reference-manual/analysis/mean-square-displacement.rst
        reference-manual/analysis/bond-angle-dihedral.rst
        reference-manual/analysis/radius-of-gyration.rst
        reference-manual/analysis/rmsd.rst
        reference-manual/analysis/covariance-analysis.rst
        reference-manual/analysis/dihedral-pca.rst
        reference-manual/analysis/hydrogen-bonds.rst
        reference-manual/analysis/protein-related.rst
        reference-manual/analysis/interface-related.rst)
    # The image files have also been ordered by the respective
    # chapter they are included in in the reference manual
    set(REFERENCEMANUAL_IMAGE_FILES
        # General folder
        reference-manual/plots/decomp.pdf
        reference-manual/plots/dih.pdf
        reference-manual/plots/drift-all.pdf
        reference-manual/plots/f-angle.pdf
        reference-manual/plots/f-bond.pdf
        reference-manual/plots/fp-highres.pdf
        reference-manual/plots/int-mat.pdf
        reference-manual/plots/mdpar.pdf
        reference-manual/plots/parsort.pdf
        reference-manual/plots/ring.pdf
        reference-manual/plots/shiftf.pdf
        # Algorithms chapter
        reference-manual/algorithms/plots/dd-cells.pdf
        reference-manual/algorithms/plots/dd-tric.pdf
        reference-manual/algorithms/plots/flowchart.pdf
        reference-manual/algorithms/plots/free1.pdf
        reference-manual/algorithms/plots/free2.pdf
        reference-manual/algorithms/plots/leapfrog.pdf
        reference-manual/algorithms/plots/lincs.pdf
        reference-manual/algorithms/plots/maxwell.pdf
        reference-manual/algorithms/plots/mpmd-pme.pdf
        reference-manual/algorithms/plots/nstric.pdf
        reference-manual/algorithms/plots/par-lincs2.pdf
        reference-manual/algorithms/plots/pbctric.pdf
        reference-manual/algorithms/plots/rhododec.pdf
        reference-manual/algorithms/plots/truncoct.pdf
        reference-manual/algorithms/plots/verlet-drift.pdf
        # Interaction functions chapter
        reference-manual/functions/plots/angle.pdf
        reference-manual/functions/plots/bstretch.pdf
        reference-manual/functions/plots/chain.pdf
        reference-manual/functions/plots/dummies.pdf
        reference-manual/functions/plots/f-bham.pdf
        reference-manual/functions/plots/fbposres.pdf
        reference-manual/functions/plots/f-dih.pdf
        reference-manual/functions/plots/f-dr.pdf
        reference-manual/functions/plots/fig-02.pdf
        reference-manual/functions/plots/fig-04.pdf
        reference-manual/functions/plots/f-imps.pdf
        reference-manual/functions/plots/f-lj.pdf
        reference-manual/functions/plots/f-morse.pdf
        reference-manual/functions/plots/f-pr.pdf
        reference-manual/functions/plots/f-rbs.pdf
        reference-manual/functions/plots/ring-imp.pdf
        reference-manual/functions/plots/softcore.pdf
        reference-manual/functions/plots/gapsys-sc.pdf
        reference-manual/functions/plots/subst-im.pdf
        reference-manual/functions/plots/tetra-im.pdf
        reference-manual/functions/plots/vcrf.pdf
        reference-manual/functions/plots/vsite-4fdn.pdf
        # Special topics chapter
        reference-manual/special/plots/awh-invN.pdf
        reference-manual/special/plots/awh-pmfs.pdf
        reference-manual/special/plots/awh-sampleweights.pdf
        reference-manual/special/plots/awh-traj.pdf
        reference-manual/special/plots/compelsetup.pdf
        reference-manual/special/plots/dumaro.pdf
        reference-manual/special/plots/dumtypes.pdf
        reference-manual/special/plots/equipotential.pdf
        reference-manual/special/plots/field.pdf
        reference-manual/special/plots/gaussians.pdf
        reference-manual/special/plots/lambda-values.pdf
        reference-manual/special/plots/pulldirrel.pdf
        reference-manual/special/plots/pull.pdf
        reference-manual/special/plots/pullref.pdf
        reference-manual/special/plots/rotation.pdf
        # Analysis chapter
        reference-manual/analysis/plots/dih-def.pdf
        reference-manual/analysis/plots/distm.pdf
        reference-manual/analysis/plots/hbond-insert.pdf
        reference-manual/analysis/plots/hbond.pdf
        reference-manual/analysis/plots/hpr-wheel.pdf
        reference-manual/analysis/plots/msdwater.pdf
        reference-manual/analysis/plots/phipsi.pdf
        reference-manual/analysis/plots/rama.pdf
        reference-manual/analysis/plots/rdfO-O.pdf
        reference-manual/analysis/plots/rdf.pdf
        reference-manual/analysis/plots/sgangle.pdf 
        )
    set(SPHINX_SOURCE_FILES
        index.rst
        texindex.rst
        download.rst
        links.dat
        api/gmxlibs.rst
        api/index.rst
        dev-manual/build-system.rst
        dev-manual/change-management.rst
        dev-manual/commitstyle.rst
        dev-manual/documentation-generation.rst
        dev-manual/contribute.rst
        dev-manual/doxygen.rst
        dev-manual/error-handling.rst
        dev-manual/formatting.rst
        dev-manual/gitlab-ci.rst
        dev-manual/gmxtree.rst
        dev-manual/includestyle.rst
        dev-manual/index.rst
        dev-manual/infrastructure.rst
        dev-manual/known-issues.rst
        dev-manual/language-features.rst
        dev-manual/naming.rst
        dev-manual/overview.rst
        dev-manual/physical_validation.rst
        dev-manual/redmine-states.png
        dev-manual/relocatable-binaries.rst
        dev-manual/reportstyle.rst
        dev-manual/style.rst
        dev-manual/testutils.rst
        dev-manual/tools.rst
        dev-manual/code-formatting.rst
        fragments/doxygen-links.rst
        how-to/index.rst
        how-to/beginners.rst
        how-to/topology.rst
        how-to/special.rst
        how-to/visualize.rst
        install-guide/index.rst
        install-guide/exotic.rst
        reference-manual/plots/GMX_logos/gmx_logo_blue.svg
        release-notes/index.rst
        release-notes/2025/major/highlights.rst
        release-notes/2025/major/features.rst
        release-notes/2025/major/performance.rst
        release-notes/2025/major/tools.rst
        release-notes/2025/major/bugs-fixed.rst
        release-notes/2025/major/removed-functionality.rst
        release-notes/2025/major/deprecated-functionality.rst
        release-notes/2025/major/portability.rst
        release-notes/2025/major/miscellaneous.rst
        release-notes/2025/major/api.rst
        release-notes/2024/2024.1.rst
        release-notes/2024/2024.2.rst
        release-notes/2024/2024.3.rst
        release-notes/2024/major/highlights.rst
        release-notes/2024/major/features.rst
        release-notes/2024/major/performance.rst
        release-notes/2024/major/tools.rst
        release-notes/2024/major/bugs-fixed.rst
        release-notes/2024/major/removed-functionality.rst
        release-notes/2024/major/deprecated-functionality.rst
        release-notes/2024/major/portability.rst
        release-notes/2024/major/miscellaneous.rst
        release-notes/2024/major/api.rst
        release-notes/2023/2023.1.rst
        release-notes/2023/2023.2.rst
        release-notes/2023/2023.3.rst
        release-notes/2023/2023.4.rst
        release-notes/2023/2023.5.rst
        release-notes/2023/2023.6.rst
        release-notes/2023/major/highlights.rst
        release-notes/2023/major/features.rst
        release-notes/2023/major/performance.rst
        release-notes/2023/major/tools.rst
        release-notes/2023/major/bugs-fixed.rst
        release-notes/2023/major/removed-functionality.rst
        release-notes/2023/major/deprecated-functionality.rst
        release-notes/2023/major/portability.rst
        release-notes/2023/major/miscellaneous.rst
        release-notes/2023/major/api.rst
        release-notes/2022/2022.1.rst
        release-notes/2022/2022.2.rst
        release-notes/2022/2022.3.rst
        release-notes/2022/2022.4.rst
        release-notes/2022/2022.5.rst
        release-notes/2022/2022.6.rst
        release-notes/2022/major/highlights.rst
        release-notes/2022/major/features.rst
        release-notes/2022/major/performance.rst
        release-notes/2022/major/tools.rst
        release-notes/2022/major/bugs-fixed.rst
        release-notes/2022/major/removed-functionality.rst
        release-notes/2022/major/deprecated-functionality.rst
        release-notes/2022/major/portability.rst
        release-notes/2022/major/miscellaneous.rst
        release-notes/2022/major/api.rst
        release-notes/2021/2021.7.rst
        release-notes/2021/2021.6.rst
        release-notes/2021/2021.5.rst
        release-notes/2021/2021.4.rst
        release-notes/2021/2021.3.rst
        release-notes/2021/2021.2.rst
        release-notes/2021/2021.1.rst
        release-notes/2021/major/highlights.rst
        release-notes/2021/major/features.rst
        release-notes/2021/major/performance.rst
        release-notes/2021/major/tools.rst
        release-notes/2021/major/bugs-fixed.rst
        release-notes/2021/major/removed-functionality.rst
        release-notes/2021/major/deprecated-functionality.rst
        release-notes/2021/major/portability.rst
        release-notes/2021/major/miscellaneous.rst
        release-notes/2020/2020.1.rst
        release-notes/2020/2020.2.rst
        release-notes/2020/2020.3.rst
        release-notes/2020/2020.4.rst
        release-notes/2020/2020.5.rst
        release-notes/2020/2020.6.rst
        release-notes/2020/2020.7.rst
        release-notes/2020/major/highlights.rst
        release-notes/2020/major/features.rst
        release-notes/2020/major/performance.rst
        release-notes/2020/major/tools.rst
        release-notes/2020/major/bugs-fixed.rst
        release-notes/2020/major/removed-functionality.rst
        release-notes/2020/major/deprecated-functionality.rst
        release-notes/2020/major/portability.rst
        release-notes/2020/major/miscellaneous.rst
        release-notes/2019/2019.6.rst
        release-notes/2019/2019.5.rst
        release-notes/2019/2019.4.rst
        release-notes/2019/2019.3.rst
        release-notes/2019/2019.2.rst
        release-notes/2019/2019.1.rst
        release-notes/2019/major/highlights.rst
        release-notes/2019/major/features.rst
        release-notes/2019/major/performance.rst
        release-notes/2019/major/tools.rst
        release-notes/2019/major/bugs-fixed.rst
        release-notes/2019/major/removed-functionality.rst
        release-notes/2019/major/deprecated-functionality.rst
        release-notes/2019/major/portability.rst
        release-notes/2019/major/miscellaneous.rst
        release-notes/2018/2018.7.rst
        release-notes/2018/2018.6.rst
        release-notes/2018/2018.5.rst
        release-notes/2018/2018.4.rst
        release-notes/2018/2018.3.rst
        release-notes/2018/2018.2.rst
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
        # the entry for user-guide/index.rst should not appear here,
        # as it will be included conditionally further down depending on
        # if the documentation will be build with the full reference
        # manual or without.
        user-guide/cmdline.rst
        user-guide/deprecation-policy.rst
        user-guide/environment-variables.rst
        user-guide/faq.rst
        user-guide/floating-point.rst
        user-guide/flow.rst
        user-guide/force-fields.rst
        user-guide/getting-started.rst
        user-guide/index.rst
        user-guide/known-issues.rst
        user-guide/managing-simulations.rst
        user-guide/mdp-options.rst
        user-guide/mdrun-features.rst
        user-guide/mdrun-performance.rst
        user-guide/run-time-errors.rst
        user-guide/security.rst
        user-guide/system-preparation.rst
        user-guide/terminology.rst
        )

    include(SphinxMacros.cmake)
    gmx_init_sphinx_setup(${SPHINX_INPUT_DIR})

    # set temporary variables for doi inclusion
    # into the manual, plain string + some wrapping
    # for release builds, and dummy string for non-release
    # builds
    if("${GMX_MANUAL_DOI}" STREQUAL "")
      # empty string means no doi, set dummy text
      set(GMX_MANUAL_DOI_STRING "This is not a release build of GROMACS, so please reference")
      set(GMX_MANUAL_DOI_STRING "${GMX_MANUAL_DOI_STRING} one of the GROMACS papers and the base release of the manual.")
    else()
      # release version, set plain old boring string
      set(GMX_MANUAL_DOI_STRING "Please reference this documentation as https://doi.org/${GMX_MANUAL_DOI}.")
    endif()
    # same for source doi, but modify the text
    if("${GMX_SOURCE_DOI}" STREQUAL "")
      # empty string means no release build
      set(GMX_SOURCE_DOI_STRING "This is not a release build of GROMACS. Please reference one of the")
      set(GMX_SOURCE_DOI_STRING "${GMX_SOURCE_DOI_STRING} GROMACS papers, as well as the base release that this version is built from.")
      set(GMX_SOURCE_DOI_STRING "${GMX_SOURCE_DOI_STRING} Also, please state what modifcations have been performed or where the version")
      set(GMX_SOURCE_DOI_STRING "${GMX_SOURCE_DOI_STRING} was sourced from.")
    else()
      # release version, give them a doi url string
      set(GMX_SOURCE_DOI_STRING "To cite the source code for this release, please cite")
      set(GMX_SOURCE_DOI_STRING "${GMX_SOURCE_DOI_STRING} https://doi.org/${GMX_SOURCE_DOI}.")
    endif()

    if(IMAGE_CONVERT_POSSIBLE)
        set(IMAGE_CONVERT_STRING "gmx_image_convert_possible")
    else()
        set(IMAGE_CONVERT_STRING "gmx_image_convert_impossible")
    endif()

    set(SPHINX_CONFIG_FILE ${SPHINX_INPUT_DIR}/conf.py)
    if (GMX_PYTHON_PACKAGE)
        set(GMXAPI_PYTHON_STAGING_DIR ${CMAKE_BINARY_DIR}/python_packaging/gmxapi/gmxapi_staging)
        # TODO: Resolve circular reference. We would like to get the CMake build-time directory for
        # the gmxapi Python package from the _gmxapi target, as we do when building sample_restraint.
        # Instead of the above hard-coded path, how can we do
        # get_target_property(GMXAPI_PYTHON_STAGING_DIR _gmxapi staging_dir)
        # in this context?
    endif ()

    gmx_configure_version_file(
        conf.cmakein.py ${SPHINX_CONFIG_FILE}
        EXTRA_VARS
            CMAKE_MINIMUM_REQUIRED_VERSION
            EXPECTED_DOXYGEN_VERSION
            EXPECTED_SPHINX_VERSION
            GMX_ADMIN_DIR
            GMX_LMFIT_REQUIRED_VERSION
            GMX_MANUAL_DOI_STRING
            GMX_TNG_MINIMUM_REQUIRED_VERSION
            GMX_SOURCE_DOI_STRING
            GMXAPI_PYTHON_STAGING_DIR
            IMAGE_CONVERT_STRING
            REGRESSIONTEST_VERSION
            REQUIRED_CUDA_COMPUTE_CAPABILITY
            REQUIRED_CUDA_VERSION
            REQUIRED_OPENCL_MIN_VERSION
            REGRESSIONTEST_MD5SUM_STRING
            RELENG_PATH
            SOURCE_MD5SUM
            SPHINX_EXTENSION_PATH
            GMX_CURRENT_CONTRIBUTORS_STRING
            GMX_PREVIOUS_CONTRIBUTORS_STRING
            GMX_CURRENT_PROJECT_LEADERS_STRING
        COMMENT "Configuring Sphinx configuration file")
    gmx_add_sphinx_input_file(${SPHINX_CONFIG_FILE})
    gmx_add_sphinx_source_files(FILES ${SPHINX_SOURCE_FILES})
    if (EXISTS ${RELENG_PATH}/docs/FileList.cmake)
        include(${RELENG_PATH}/docs/FileList.cmake)
        gmx_add_sphinx_source_files(
            FROM ${RELENG_PATH}/docs TO dev-manual/releng PREFIX releng/docs/
            FILES ${RELENG_SPHINX_FILES})
    else()
        gmx_add_sphinx_source_files(FILES
            dev-manual/releng/index.rst
            )
    endif()

    gmx_add_sphinx_source_files(
            FILES
            gmxapi/index.rst
            gmxapi/userguide/install.rst
            gmxapi/userguide/usage.rst
    )

    if (GMX_PYTHON_PACKAGE)
        gmx_add_sphinx_source_files(
            FILES
            gmxapi/userguide/pythonreference.rst
            )
    else()
        gmx_add_sphinx_source_files(
            FROM ${CMAKE_CURRENT_SOURCE_DIR}/gmxapi/userguide-stub
            TO gmxapi/userguide/ 
            FILES pythonreference.rst 
            )
    endif ()

    gmx_add_sphinx_source_files(
            FILES
            nblib/index.rst
            nblib/guide-to-writing-MD-programs.rst
    )

    gmx_add_sphinx_source_files(
        FILES
        ${REFERENCEMANUAL_SPHINX_FILES_GENERAL})
    if (IMAGE_CONVERT_POSSIBLE)
        gmx_add_sphinx_source_files(
            FILES
            ${REFERENCEMANUAL_SPHINX_FILES_WITH_IMAGES}
            ${REFERENCEMANUAL_IMAGE_FILES})
        gmx_add_sphinx_image_conversion_files(
            FILES
            ${REFERENCEMANUAL_IMAGE_FILES})
    endif()
    gmx_add_sphinx_input_target(sphinx-input-rst)
    gmx_add_sphinx_image_conversion_target(sphinx-image-conversion)
    add_custom_target(sphinx-input)
    add_dependencies(sphinx-input sphinx-input-rst sphinx-image-conversion)
    if (GMX_PYTHON_PACKAGE)
        add_dependencies(sphinx-input _gmxapi)
    endif()
    # Remove other rst files from the build tree, since they confuse Sphinx.
    # Skip generated files in onlinehelp/, and fragments.
    # The latter do not cause issues with obsolete files, as they
    # are not considered as Sphinx input files, but will only be
    # included using an explicit .. include::.
    gmx_remove_obsolete_sphinx_input_files("^(onlinehelp|fragments)/.*\\\\.rst$")

    # TODO: Make this remove obsolete .rst files.
    # TODO: This does not work in cross-compilation scenarios; disable up to
    # the necessary level.
    gmx_add_custom_output_target(sphinx-programs
        OUTPUT ${SPHINX_INPUT_DIR}/conf-man.py
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
    # Note: the automatically generated builder tags are not available until after
    # conf.py is processed, so we use additional custom tags in our sphinx targets.
    add_custom_target(install-guide
        COMMAND
            ${SPHINX_EXECUTABLE}
            -q -b text
            -w sphinx-install.log
            -d ${CMAKE_CURRENT_BINARY_DIR}/install-guide/_doctrees
            -c ${SPHINX_INPUT_DIR}
            -t do_text
            "${SPHINX_INPUT_DIR}/install-guide"
            "${TEXT_INSTALL_GUIDE_OUTPUT_DIR}"
            "${SPHINX_INPUT_DIR}/install-guide/index.rst"
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

    set(SPHINX_CONFIG_OVERRIDES "")
    if (GMX_DEVELOPER_BUILD)
        set(SPHINX_CONFIG_OVERRIDES "-Dtodo_include_todos=1")
    endif()
    add_custom_target(webpage-sphinx
        DEPENDS sphinx-programs
        DEPENDS sphinx-input
        DEPENDS sphinx-image-conversion
        DEPENDS manual
        # Note: the automatically generated builder tags are not available until after
        # conf.py is processed, so we use additional custom tags in our sphinx targets.
        COMMAND
            ${SPHINX_EXECUTABLE}
            -q -b html
            -w sphinx-html.log
            -d "${CMAKE_CURRENT_BINARY_DIR}/_html_doctrees"
            -t do_html
            ${SPHINX_CONFIG_OVERRIDES}
            "${SPHINX_INPUT_DIR}"
            "${HTML_OUTPUT_DIR}"
        WORKING_DIRECTORY
            ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Building HTML documentation with Sphinx"
        VERBATIM
        )

    # Note: the automatically generated builder tags are not available until after
    # conf.py is processed, so we use additional custom tags in our sphinx targets.
    add_custom_target(man
        DEPENDS ${SPHINX_INPUT_DIR}/conf-man.py
        COMMAND
            ${SPHINX_EXECUTABLE}
            -q -b man
            -w sphinx-man.log
            -d "${CMAKE_CURRENT_BINARY_DIR}/_man_doctrees"
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

else()
    set(MANUAL_BUILD_IS_POSSIBLE OFF)
    set(MANUAL_BUILD_NOT_POSSIBLE_REASON "Sphinx expected minimum version ${EXPECTED_SPHINX_VERSION} is not available")

    add_custom_target(webpage-sphinx
        COMMAND ${CMAKE_COMMAND} -E echo
            "HTML pages cannot be built because Sphinx expected minimum version ${EXPECTED_SPHINX_VERSION} is not available"
        VERBATIM)
    add_custom_target(install-guide
        COMMAND ${CMAKE_COMMAND} -E echo
            "INSTALL cannot be built because Sphinx expected minimum version ${EXPECTED_SPHINX_VERSION} is not available"
        VERBATIM)
    add_custom_target(man
        COMMAND ${CMAKE_COMMAND} -E echo
            "man pages cannot be built because Sphinx expected minimum version ${EXPECTED_SPHINX_VERSION} is not available"
        VERBATIM)
    add_custom_target(sphinx-create-texman
        COMMAND ${CMAKE_COMMAND} -E echo
            "Cannot prepare LaTeX input files because Sphinx expected minimum version ${EXPECTED_SPHINX_VERSION} is not available"
        VERBATIM)
    add_custom_target(manual
        COMMAND ${CMAKE_COMMAND} -E echo
            "manual cannot be built because Sphinx expected minimum version ${EXPECTED_SPHINX_VERSION} is not available")
endif()

if (MAN_PAGE_DIR)
    set(MAN_PAGE_DIR ${MAN_PAGE_DIR}/man)
    # Trailing slash on directory is significant for
    # install(DIRECTORY). See CMake docs.
    install(DIRECTORY ${MAN_PAGE_DIR}/
        DESTINATION ${CMAKE_INSTALL_MANDIR}/man1
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

# Next, turn it off if any of the preconditions are unsatisfied
if (NOT Python3_Interpreter_FOUND)
    set(HTML_BUILD_IS_POSSIBLE OFF)
    set(HTML_BUILD_NOT_POSSIBLE_REASON "Python is required")
elseif (NOT SPHINX_FOUND)
    # Hardly anything gets built if Sphinx is not available, so don't bother.
    set(HTML_BUILD_IS_POSSIBLE OFF)
    set(HTML_BUILD_NOT_POSSIBLE_REASON "Sphinx expected minimum version ${EXPECTED_SPHINX_VERSION} is required")
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
                -E copy ${CMAKE_CURRENT_BINARY_DIR}/manual/gromacs.pdf ${_manual_target_location}
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
    # for CI (and possibly others) to use
    add_custom_target(webpage ${_webpage_target_properties}
        COMMAND ${CMAKE_COMMAND} -E echo
            "Webpage output: file://${HTML_OUTPUT_DIR}/index.html"
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
