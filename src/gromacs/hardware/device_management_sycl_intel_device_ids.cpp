/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2022- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 *  \brief Defines the mapping from Intel PCIE device ID to hardware version.
 *
 *  \author Andrey Alekseenko <al42and@gmail.com>
 *
 * \ingroup module_hardware
 */
#include "gmxpre.h"

#include "device_management_sycl_intel_device_ids.h"

#include <optional>
#include <unordered_map>
#include <unordered_set>

// Based on PRODUCT_CONFIG enum from
// https://github.com/intel/compute-runtime/blob/d75eccc0269be71ebb359d93ca7fff36cb03b9d1/third_party/aot_config_headers/platforms.h
enum class IntelProductConfig : unsigned int
{
    BDW        = 0x02000000, //! Broadwell (Gen8), 8.0.0
    SKL        = 0x02400009, //! Skylake (Gen9), 9.0.9
    KBL        = 0x02404009, //! Kaby Lake (Gen9p5), 9.1.9
    CFL        = 0x02408009, //! Coffee Lake, 9.2.9
    APL        = 0x0240c000, //! Apollo Lake (Gen9LP) 9.3.0
    GLK        = 0x02410000, //! Gemini Lake (Gen9LP) 9.4.0
    WHL        = 0x02414000,
    AML        = 0x02418000,
    CML        = 0x0241c000,
    ICL        = 0x02c00000, //! Ice Lake (Gen11), 11.0.0
    LKF        = 0x02c04000,
    EHL        = 0x02c08000,
    TGL        = 0x03000000, //! Tiger Lake (Gen12), 12.0.0
    RKL        = 0x03004000,
    RPL_S      = 0x03008000,
    ADL_S      = 0x03008000, //! Alder Lake S (Gen12), 12.2.0
    ADL_P      = 0x0300c000,
    ADL_N      = 0x03010000,
    DG1        = 0x03028000, //! DG1 (XeLP), 12.10.0
    XEHP_SDV   = 0x030c8004, //! XeHP SDV, 12.50.4
    DG2_G10_A0 = 0x030dc000, //! DG2 G10 (Arc, Xe-HPG), 12.55.0
    DG2_G10_A1 = 0x030dc001,
    DG2_G10_B0 = 0x030dc004,
    DG2_G10_C0 = 0x030dc008,
    DG2_G11_A0 = 0x030e0000, //! DG2 G11 (Arc, Xe-HPG), 12.56.0
    DG2_G11_B0 = 0x030e0004,
    DG2_G11_B1 = 0x030e0005,
    DG2_G12_A0 = 0x030e4000, //! DG2 G12(Arc, Xe-HPG), 12.57.0
    PVC_XL_A0  = 0x030f0000, //! Ponte Veccio (Xe-HPC), 12.60.0
    PVC_XL_A0P = 0x030f0001,
    PVC_XT_A0  = 0x030f0003,
    PVC_XT_B0  = 0x030f0005,
    PVC_XT_B1  = 0x030f0006,
    PVC_XT_C0  = 0x030f0007,
};

// Map from IntelProductConfig to list of PCIE IDs
/*
 * Check out https://github.com/intel/compute-runtime, then use the following script:
 *
#!/bin/bash
set -euo pipefail
IFS=$'\n\t'
BIN="$(mktemp)"
SRC="${BIN}.cpp"
cat >"${SRC}" <<EOF
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
std::ostream& operator<<(std::ostream& os, const std::vector<short unsigned int>& vec) {
  for (const auto e: vec) {
    os << "0x" << std::setw(4) << std::setfill('0') << std::hex << e << ", ";
  }
  return os;
}
#define DEVICE_CONFIG(type, a, list, ...) \
  std::cout << "{ IntelProductConfig::" << #type << ", { " << NEO::list << "} }," << std::endl;
EOF
for i in $(find ./shared/source/ -name 'device_ids_configs_*.h'); do
   grep -v 'include\|pragma' $i >> "${SRC}"
   echo >> "${SRC}"
done
echo 'int main() {' >> "${SRC}"
grep ^DEVICE_CONFIG shared/source/dll/devices/product_config_base.inl >> "${SRC}"
echo 'return 0;}' >> "${SRC}"
g++ "${SRC}" -o ${BIN}
${BIN}
rm ${SRC} ${BIN}
 *
 * to generate the list below
 */
static const std::unordered_map<IntelProductConfig, std::unordered_set<unsigned int>> c_pciExpressIdsForProduct = {
    { IntelProductConfig::PVC_XL_A0,
      {
              0x0bd0,
      } },
    { IntelProductConfig::PVC_XL_A0P,
      {
              0x0bd0,
      } },
    { IntelProductConfig::PVC_XT_A0,
      {
              0x0bd5,
              0x0bd6,
              0x0bd7,
              0x0bd8,
              0x0bd9,
              0x0bda,
              0x0bdb,
      } },
    { IntelProductConfig::PVC_XT_B0,
      {
              0x0bd5,
              0x0bd6,
              0x0bd7,
              0x0bd8,
              0x0bd9,
              0x0bda,
              0x0bdb,
      } },
    { IntelProductConfig::PVC_XT_B1,
      {
              0x0bd5,
              0x0bd6,
              0x0bd7,
              0x0bd8,
              0x0bd9,
              0x0bda,
              0x0bdb,
      } },
    { IntelProductConfig::PVC_XT_C0,
      {
              0x0bd5,
              0x0bd6,
              0x0bd7,
              0x0bd8,
              0x0bd9,
              0x0bda,
              0x0bdb,
      } },
    { IntelProductConfig::DG2_G10_A0,
      {
              0x4f80,
              0x4f81,
              0x4f82,
              0x4f83,
              0x4f84,
              0x5690,
              0x5691,
              0x5692,
              0x56a0,
              0x56a1,
              0x56a2,
              0x56c0,
      } },
    { IntelProductConfig::DG2_G10_A1,
      {
              0x4f80,
              0x4f81,
              0x4f82,
              0x4f83,
              0x4f84,
              0x5690,
              0x5691,
              0x5692,
              0x56a0,
              0x56a1,
              0x56a2,
              0x56c0,
      } },
    { IntelProductConfig::DG2_G10_B0,
      {
              0x4f80,
              0x4f81,
              0x4f82,
              0x4f83,
              0x4f84,
              0x5690,
              0x5691,
              0x5692,
              0x56a0,
              0x56a1,
              0x56a2,
              0x56c0,
      } },
    { IntelProductConfig::DG2_G10_C0,
      {
              0x4f80,
              0x4f81,
              0x4f82,
              0x4f83,
              0x4f84,
              0x5690,
              0x5691,
              0x5692,
              0x56a0,
              0x56a1,
              0x56a2,
              0x56c0,
      } },
    { IntelProductConfig::DG2_G11_A0,
      {
              0x4f87,
              0x4f88,
              0x5693,
              0x5694,
              0x5695,
              0x56a5,
              0x56a6,
              0x56b0,
              0x56b1,
              0x56c1,
      } },
    { IntelProductConfig::DG2_G11_B0,
      {
              0x4f87,
              0x4f88,
              0x5693,
              0x5694,
              0x5695,
              0x56a5,
              0x56a6,
              0x56b0,
              0x56b1,
              0x56c1,
      } },
    { IntelProductConfig::DG2_G11_B1,
      {
              0x4f87,
              0x4f88,
              0x5693,
              0x5694,
              0x5695,
              0x56a5,
              0x56a6,
              0x56b0,
              0x56b1,
              0x56c1,
      } },
    { IntelProductConfig::DG2_G12_A0,
      {
              0x5696,
              0x5697,
              0x56a3,
              0x56a4,
              0x56b2,
              0x56b3,
              0x4f85,
              0x4f86,
      } },
    { IntelProductConfig::XEHP_SDV,
      {
              0x0201,
              0x0202,
              0x0203,
              0x0204,
              0x0205,
              0x0206,
              0x0207,
              0x0208,
              0x0209,
              0x020a,
              0x020b,
              0x020c,
              0x020d,
              0x020e,
              0x020f,
              0x0210,
      } },
    { IntelProductConfig::TGL,
      {
              0xff20,
              0x9a49,
              0x9a40,
              0x9a59,
              0x9a60,
              0x9a68,
              0x9a70,
              0x9a78,
      } },
    { IntelProductConfig::DG1,
      {
              0x4905,
              0x4906,
              0x4907,
              0x4908,
      } },
    { IntelProductConfig::RKL,
      {
              0x4c80,
              0x4c8a,
              0x4c8b,
              0x4c8c,
              0x4c90,
              0x4c9a,
      } },
    { IntelProductConfig::ADL_S,
      {
              0x4680,
              0x4682,
              0x4688,
              0x468a,
              0x4690,
              0x4692,
              0x4693,
              0xa780,
              0xa781,
              0xa782,
              0xa783,
              0xa788,
              0xa789,
              0xa78b,
      } },
    { IntelProductConfig::ADL_P,
      {
              0x46a0, 0x46b0, 0x46a1, 0x46a2, 0x46a3, 0x46a6, 0x46a8, 0x46aa,
              0x462a, 0x4626, 0x4628, 0x46b1, 0x46b2, 0x46b3, 0x46c0, 0x46c1,
              0x46c2, 0x46c3, 0xa7a0, 0xa720, 0xa7a8, 0xa7a1, 0xa721, 0xa7a9,
      } },
    { IntelProductConfig::ADL_N,
      {
              0x46d0,
              0x46d1,
              0x46d2,
      } },
    { IntelProductConfig::ICL,
      {
              0xff05,
              0x8a56,
              0x8a58,
              0x8a5c,
              0x8a5a,
              0x8a50,
              0x8a52,
              0x8a51,
      } },
    { IntelProductConfig::EHL,
      {
              0x4500,
              0x4541,
              0x4551,
              0x4571,
              0x4555,
              0x4e51,
              0x4e61,
              0x4e71,
              0x4e55,
      } },
    { IntelProductConfig::LKF,
      {
              0x9840,
      } },
    { IntelProductConfig::SKL,
      {
              0x1902, 0x190b, 0x190a, 0x1906, 0x190e, 0x1917, 0x1913, 0x1915, 0x1912,
              0x191b, 0x191a, 0x1916, 0x191e, 0x191d, 0x1921, 0x9905, 0x192b, 0x192d,
              0x192a, 0x1923, 0x1926, 0x1927, 0x1932, 0x193b, 0x193a, 0x193d,
      } },
    { IntelProductConfig::KBL,
      {
              0x5902, 0x590b, 0x590a, 0x5906, 0x590e, 0x5908, 0x5913, 0x5915, 0x5912,
              0x591b, 0x5917, 0x591a, 0x5916, 0x591e, 0x591d, 0x5921, 0x5926, 0x5927,
              0x592b, 0x592a, 0x5923, 0x5932, 0x593b, 0x593a, 0x593d,
      } },
    { IntelProductConfig::AML,
      {
              0x591c,
      } },
    { IntelProductConfig::CFL,
      {
              0x3e90, 0x3e93, 0x3e99, 0x3e92, 0x3e9b, 0x3e94, 0x3e91, 0x3e96, 0x3e9a, 0x3ea9,
              0x3e98, 0x3e95, 0x3ea6, 0x3ea7, 0x3ea8, 0x3ea5, 0x9bab, 0x9ba0, 0x9bc0,
      } },
    { IntelProductConfig::CML,
      {
              0x9b21,
              0x9b41,
              0x9ba2,
              0x9ba4,
              0x9ba5,
              0x9ba8,
              0x9baa,
              0x9bac,
              0x9bc2,
              0x9bc4,
              0x9bc5,
              0x9bc6,
              0x9bc8,
              0x9bca,
              0x9bcb,
              0x9bcc,
              0x9be6,
              0x9bf6,
      } },
    { IntelProductConfig::WHL,
      {
              0x3ea1,
              0x3ea3,
              0x3ea4,
              0x3ea0,
              0x3ea2,
      } },
    { IntelProductConfig::GLK,
      {
              0x3184,
              0x3185,
      } },
    { IntelProductConfig::APL,
      {
              0x9906,
              0x9907,
              0x0a84,
              0x5a84,
              0x5a85,
              0x1a85,
              0x1a84,
              0x9908,
      } },
    { IntelProductConfig::BDW,
      {
              0x1602,
              0x160a,
              0x1606,
              0x160e,
              0x160d,
              0x1612,
              0x161a,
              0x1616,
              0x161e,
              0x161d,
              0x1622,
              0x162a,
              0x1626,
              0x162b,
              0x162e,
              0x162d,
      } },
};


static constexpr std::tuple<int, int, int> getHardwareVersionFromIntelProductConfig(const IntelProductConfig productConfig)
{
    // Convert IntelProductConfig value into (major, minor) tuple.
    // PRODUCT_CONFIG layout is described in https://github.com/intel/compute-runtime/blob/d75eccc0269be71ebb359d93ca7fff36cb03b9d1/shared/source/helpers/product_config_helper.h#L31
    const int major = ((static_cast<unsigned>(productConfig) & 0xffc00000) >> 22); // first 10 bits
    const int minor = ((static_cast<unsigned>(productConfig) & 0x3fc000) >> 14);   // next 8 bits
    // next 8 bits are reserved
    const int revision = (static_cast<unsigned>(productConfig) & 0x3f); // last 6 bits
    return std::make_tuple(major, minor, revision);
}

static std::optional<IntelProductConfig> getProductConfigFromPciExpressID(unsigned int pciExpressID)
{
    for (const auto& [intelProduct, pciExpressIds] : c_pciExpressIdsForProduct)
    {
        if (pciExpressIds.find(pciExpressID) != pciExpressIds.cend())
        {
            return intelProduct;
        }
    }
    return std::nullopt;
}

std::optional<std::tuple<int, int, int>> getIntelHardwareVersionFromPciExpressID(unsigned int pciExpressID)
{
    const auto productConfig = getProductConfigFromPciExpressID(pciExpressID);
    if (productConfig.has_value())
    {
        return getHardwareVersionFromIntelProductConfig(productConfig.value());
    }
    else
    {
        return std::nullopt;
    }
}
