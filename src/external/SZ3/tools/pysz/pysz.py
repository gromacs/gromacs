import ctypes
from ctypes.util import find_library
import numpy as np

"""
Python API for SZ2/SZ3
"""


class SZ:
    def __init__(self, szpath):
        """
        init SZ
        :param szpath: the path to SZ dynamic library
        """
        self.sz = ctypes.cdll.LoadLibrary(szpath)

        self.sz.SZ_compress_args.argtypes = (ctypes.c_int, ctypes.c_void_p, ctypes.POINTER(ctypes.c_size_t),
                                             ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double,
                                             ctypes.c_size_t, ctypes.c_size_t, ctypes.c_size_t, ctypes.c_size_t,
                                             ctypes.c_size_t)
        self.sz.SZ_compress_args.restype = ctypes.POINTER(ctypes.c_ubyte)

        self.sz.SZ_decompress.argtypes = (ctypes.c_int, ctypes.POINTER(ctypes.c_ubyte), ctypes.c_size_t,
                                          ctypes.c_size_t, ctypes.c_size_t, ctypes.c_size_t, ctypes.c_size_t,
                                          ctypes.c_size_t)

        self.libc = ctypes.CDLL(ctypes.util.find_library('c'))
        self.libc.free.argtypes = (ctypes.c_void_p,)

    def __sz_datatype(self, dtype, data=None):
        if dtype == np.float32:
            return 0, data.ctypes.data_as(ctypes.POINTER(ctypes.c_float)) if data is not None else None
        elif dtype == np.float64:
            return 1, data.ctypes.data_as(ctypes.POINTER(ctypes.c_double)) if data is not None else None
        else:
            print('SZ currently supports float32 and float64\n')
            exit(0)

    def verify(self, src_data, dec_data):
        """
        Compare the decompressed data with original data
        :param src_data: original data, numpy array
        :param dec_data: decompressed data, numpy array
        :return: max_diff, psnr, nrmse
        """
        data_range = np.max(src_data) - np.min(src_data)
        diff = src_data - dec_data
        max_diff = np.max(abs(diff))
        print("abs err={:.8G}".format(max_diff))
        mse = np.mean(diff ** 2)
        nrmse = np.sqrt(mse) / data_range
        psnr = 20 * np.log10(data_range) - 10 * np.log10(mse)
        return max_diff, psnr, nrmse

    def decompress(self, data_cmpr, original_shape, original_dtype):
        """
        Decompress data with SZ
        :param data_cmpr: compressed data, numpy array format, dtype should be np.uint8
        :param original_shape: the shape of original data
        :param original_dtype: the dtype of original data
        :return: decompressed data,numpy array format
        """

        r5, r4, r3, r2, r1 = [0] * (5 - len(original_shape)) + list(original_shape)
        ori_type, ori_null = self.__sz_datatype(original_dtype)
        self.sz.SZ_decompress.restype = ctypes.POINTER(
            ctypes.c_float if original_dtype == np.float32 else ctypes.c_double)
        data_dec_c = self.sz.SZ_decompress(ori_type,
                                           data_cmpr.ctypes.data_as(ctypes.POINTER(ctypes.c_ubyte)),
                                           data_cmpr.size,
                                           r5, r4, r3, r2, r1)

        data_dec = np.array(data_dec_c[:np.prod(original_shape)]).reshape(original_shape)
        self.libc.free(data_dec_c)
        return data_dec

    def compress(self, data, eb_mode, eb_abs, eb_rel, eb_pwr):
        """
        Compress data with SZ
        :param data: original data, numpy array format, dtype is FP32 or FP64
        :param eb_mode:# error bound mode, integer (0: ABS, 1:REL, 2:ABS_AND_REL, 3:ABS_OR_REL, 4:PSNR, 5:NORM, 10:PW_REL)
        :param eb_abs: optional, abs error bound, double
        :param eb_rel: optional, rel error bound, double
        :param eb_pwr: optional, pwr error bound, double
        :return: compressed data, numpy array format, dtype is np.uint8
                 compression ratio
        """
        assert len(data.shape) <= 5, "SZ only supports 1D to 5D input data"
        cmpr_size = ctypes.c_size_t()
        r5, r4, r3, r2, r1 = [0] * (5 - len(data.shape)) + list(data.shape)
        datatype, datap = self.__sz_datatype(data.dtype, data)
        data_cmpr_c = self.sz.SZ_compress_args(datatype, datap,
                                               ctypes.byref(cmpr_size),
                                               eb_mode, eb_abs, eb_rel, eb_pwr,
                                               r5, r4, r3, r2, r1)

        cmpr_ratio = data.size * data.itemsize / cmpr_size.value

        data_cmpr = np.array(data_cmpr_c[:cmpr_size.value], dtype=np.uint8)
        self.libc.free(data_cmpr_c)
        return data_cmpr, cmpr_ratio
