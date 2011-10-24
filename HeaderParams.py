import ctypes as C

header_keys = {
    "HEADER_START":None,
    "HEADER_END":None,
    "filename": 'str',
    "telescope_id": 'I',
    "machine_id": 'I',
    "data_type": 'I',
    "rawdatafile": 'str',
    "source_name": 'str',
    "barycentric": 'I',
    "pulsarcentric": 'I',
    "az_start": 'd',
    "za_start": 'd',
    "src_raj": 'd',
    "src_dej": 'd',
    "tstart": 'd',
    "tsamp": 'd',
    "nbits": 'I',
    "nsamples": 'I',
    "fch1": 'd',
    "foff": 'd',
    "fchannel": 'd',
    "nchans": 'I',
    "nifs": 'I',
    "refdm": 'd',
    "period": 'd',
    "nbeams": 'I',
    "ibeam": 'I',
    "hdrlen": 'I',
    "orig_hdrlen": 'I',
    "new_hdrlen": 'I',
    "sampsize": 'I',
    "bandwidth": 'd',
    "fbottom": 'd',
    "date": 'str',
    "signed": 'b'}

header_units = {
    "filename": "",
    "telescope_id": "",
    "machine_id": "",
    "data_type": "",
    "rawdatafile": "",
    "source_name": "",
    "barycentric": "",
    "pulsarcentric": "",
    "az_start": "(Degrees)",
    "za_start": "(Degrees)",
    "src_raj": "(HH:MM:SS.ss)",
    "src_dej": "(DD:MM:SS.ss)",
    "tstart": '(MJD)',
    "tsamp": '(s)',
    "nbits": "",
    "nsamples": "",
    "fch1": "(MHz)",
    "foff": "(MHz)",
    "fchannel": "(MHz)",
    "nchans": "",
    "nifs": "",
    "refdm": "(pccm^-3)",
    "period": "(s)",
    "nbeams": "",
    "ibeam": "",
    "hdrlen": "(Bytes)",
    "orig_hdrlen": "(Bytes)",
    "new_hdrlen": "(Bytes)",
    "sampsize": "(Bytes)",
    "bandwidth": "(MHz)",
    "fbottom": "(MHz)",
    "date":"",
    "signed":""}

data_types = {
    1:"Filterbank file",
    2:"Timeseries file"}

struct_to_numpy = {
    "I":"uint",
    "d":"float",
    "str":"S256"}

telescope_ids = {
    "Fake": 0,
    "Arecibo": 1,
    "Ooty": 2,
    "Nancay": 3,
    "Parkes": 4, 
    "Jodrell": 5,
    "GBT": 6,
    "GMRT": 7,
    "Effelsberg": 8,
    "Effelsberg LOFAR":9,
    "Other": 10}

ids_to_telescope = dict(zip(telescope_ids.values(), telescope_ids.keys()))

machine_ids = {
    "FAKE": 0,
    "PSPM": 1,
    "Wapp": 2,
    "AOFTM": 3,
    "BCPM1": 4,
    "OOTY": 5,
    "SCAMP": 6,
    "GBT Pulsar Spigot": 7,
    "PFFTS": 8}

ids_to_machine = dict(zip(machine_ids.values(), machine_ids.keys()))

telescope_lats_longs = {
    "Effelsberg":(50.52485,6.883593)
    }

nptypes_to_ctypes = {"|b1":C.c_bool,
                     "|S1":C.c_char,
                     "<i1":C.c_byte,
                     "<u1":C.c_ubyte,
                     "<i2":C.c_short,
                     "<u2":C.c_ushort,
                     "<i4":C.c_int,
                     "<u4":C.c_uint,
                     "<i8":C.c_long,
                     "<u8":C.c_ulong,
                     "<f4":C.c_float,
                     "<f8":C.c_double}
