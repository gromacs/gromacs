import numpy as np
from pathlib import Path
from pysz import SZ
import platform

# prepare your data in numpy array format
HOME = str(Path.home())
data = np.fromfile(HOME + '/data/hurricane-100x500x500/Uf48.bin.dat', dtype=np.float32)
data = np.reshape(data, (100, 500, 500))


# init SZ (both SZ2 and SZ3 are supported)
# Please change the path to the SZ dynamic library file in your system
lib_extention = "so" if platform.system() == 'Linux' else "dylib"
sz = SZ("../../build/tools/sz3c/libSZ3c.{}".format(lib_extention))
# sz = SZ("../../../sz2/build/sz/libSZ.{}".format(lib_extention))

# compress, both input and output data are numpy array
data_cmpr, cmpr_ratio = sz.compress(data, 0, 1e-3, 0, 0)
print("compression ratio = {:5G}".format(cmpr_ratio))

# decompress, both input and output data are numpy array
data_dec = sz.decompress(data_cmpr, data.shape, data.dtype)

# verify
sz.verify(data, data_dec)