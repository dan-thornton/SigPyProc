import os,sys,struct,numpy,string,math
import HeaderParams as conf

class Header:
    def __init__(self,indata):
        if isinstance(indata,str):
            if os.path.isfile(indata):
                self.file = indata
                self.basename = self.file.strip(".fil")
                self.header = self._header_from_file(indata)
                self.basename = self.file.strip(".fil")
            else:
                raise ValueError,"File not found..."

        elif isinstance(indata,dict):
            self.header = self._header_from_dict(indata)
            self.file   = None
            self.hdrlen = len(self.write_header()) 
            if "filename" in indata.keys():
                self.basename = indata["filename"].strip("fil")
        
        else:
            raise ValueError,"Bad input data type to Header.__init__"
        
        self.original_header = self.header.copy()
        self.info = {}
        self._get_info()

    def modify_header(self,key,value):
        if key in conf.header_keys.keys():
            self.header[key] = value 
            self._get_info()
        else:
            print "key data type undefined"

    def remove_from_header(self,key):
        del self.header[key]

    def reset_header(self):
        self.header = self.original_header.copy()
        self._get_info()

    def add_to_mjd(self,time_in_secs):
        if not self.header['tstart'] == 0.0:
            val = self.header['tstart']+(time_in_secs/86400.)
            self.modify_header("tstart",val)

    def make_inf(self,outfile=None):
         inf = (" Data file name without suffix          =  %s\n"%(self.basename),
                " Telescope used                         =  Effelsberg\n",
                " Instrument used                        =  PFFTS\n",
                " Object being observed                  =  %s\n"%(self.info["source_name"]),
                " J2000 Right Ascension (hh:mm:ss.ssss)  =  %s\n"%(self._reformat_radec(self.info["src_raj"])),
                " J2000 Declination     (dd:mm:ss.ssss)  =  %s\n"%(self._reformat_radec(self.info["src_dej"])),
                " Data observed by                       =  Ewan Barr\n",
                " Epoch of observation (MJD)             =  %.09f\n"%(self.info["tstart"]),
                " Barycentered?           (1=yes, 0=no)  =  %d\n"%(self.info.get("barycentric",0)),
                " Number of bins in the time series      =  %d\n"%(self.info["nsamples"]),     
                " Width of each time series bin (sec)    =  %.17gf\n"%(self.info["tsamp"]),
                " Any breaks in the data? (1=yes, 0=no)  =  0\n",
                " Type of observation (EM band)          =  Radio\n",
                " Beam diameter (arcsec)                 =  9.22\n",
                " Dispersion measure (cm-3 pc)           =  %.03f\n"%(self.info.get("refdm",0)),
                " Central freq of low channel (Mhz)      =  %.05f\n"%(self.info.get("fbottom",self.info["fch1"])),
                " Total bandwidth (Mhz)                  =  %.05f\n"%(self.info.get("bandwidth",0)),
                " Number of channels                     =  %d\n"%(self.info.get("nchans",1)),
                " Channel bandwidth (Mhz)                =  %.09f\n"%(abs(self.info.get("foff",0))),
                " Data analyzed by                       =  Ewan Barr\n",
                " Any additional notes:",
                "    Input filterbank samples have %d bits."%(self.info["nbits"]))
         
         if outfile==None:
             return "".join(inf)
         else:
             f = open(outfile,"w+")
             f.write("".join(inf))
             f.close()
                 
                
    def display(self):
        self._get_info()
        self._get_formatted_keys()

        for key in self.info.keys():
            if key in ["src_dej","src_raj"]:
                print "%s %s"%(self.formatted_keys[key],self._reformat_radec(self.info[key]))
            elif key == "telescope_id":
                print "%s %s"%(self.formatted_keys[key],conf.ids_to_telescope.get(self.info[key],"unknown"))
            elif key == "machine_id":
                print "%s %s"%(self.formatted_keys[key],conf.ids_to_machine.get(self.info[key],"unknown"))
            elif key == "data_type":
                if self.info[key] in conf.data_types.keys():
                    print "%s %s"%(self.formatted_keys[key],conf.data_types[self.info[key]])
                else:
                    print "%s Unknown"%(self.formatted_keys[key])
            else:
                print "%s %s"%(self.formatted_keys[key],self.info[key])

    def write_header(self):
        hstart  = "HEADER_START"
        hend    = "HEADER_END"
        header  = "".join([struct.pack("I",len(hstart)),hstart])
        for key in self.header.keys():
            if conf.header_keys[key] == "str":
                header = "".join([header,self._write_string(key,self.header[key])])
            elif conf.header_keys[key] == "I":
                header = "".join([header,self._write_int(key,self.header[key])])
            elif conf.header_keys[key] == "d":
                header = "".join([header,self._write_double(key,self.header[key])])
            elif conf.header_keys[key] == "b":
                header = "".join([header,self._write_char(key,self.header[key])])
        return "".join([header,struct.pack("I",len(hend)),hend])


    def _header_from_dict(self,infodict):
        header = {}
        for key in infodict.keys():
            if key in conf.header_keys:
                header[key] = infodict[key]
        return header

    def _header_from_file(self,file):
        self._f = open(file,"r")
        header = {}
        while True:
            keylen = struct.unpack("I",self._f.read(4))[0]
            key = self._f.read(keylen)
            if not key in conf.header_keys: 
                print "'%s' not recognised header key"%(key)
                return None
            
            if conf.header_keys[key] == "str":
                value = self._read_string()
                header[key] = value
            elif conf.header_keys[key] == "I":
                value = self._read_int()
                header[key] = value
            elif conf.header_keys[key] == "b":
                value = self._read_char()
                header[key] = value
            elif conf.header_keys[key] == "d":
                value = self._read_double()
                header[key] = value
            
            if key == "HEADER_END":
                break

        self.hdrlen = self._f.tell()
        self._f.seek(0,2)
        self.filelen = self._f.tell()
        self._f.close()
        return header

    def _read_char(self):
        return struct.unpack("b",self._f.read(1))[0]

    def _read_string(self):
        strlen = struct.unpack("I",self._f.read(4))[0]
        return self._f.read(strlen)
    
    def _read_int(self):
        return struct.unpack("I",self._f.read(4))[0]

    def _read_double(self):
        return struct.unpack("d",self._f.read(8))[0]
    
    def _write_string(self,key,value):
        return "".join([struct.pack("I", len(key)),
                        key,struct.pack('I',len(value)),
                        value])

    def _write_int(self,key,value):
        return "".join([struct.pack('I',len(key)),
                        key,struct.pack('I',value)])

    def _write_double(self,key,value):
        return "".join([struct.pack('I',len(key)),
                        key,struct.pack('d',value)])
    def _write_char(self,key,value):
        return "".join([struct.pack('I',len(key)),
                        key,struct.pack('b',value)])
    def _get_info(self):
        
        for key,val in self.header.items():
            self.info[key] = val
        
        self.info['orig_hdrlen']    = self.hdrlen
        self.info['new_hdrlen']     = len(self.write_header())
        if self.header['data_type'] == 1:
            self.info['bandwidth']  = abs(self.header['foff'])*self.header['nchans']
            self.info['fbottom']    = self.header['fch1']-abs(self.info['bandwidth'])

        if self.file: 
            self.info['filename']   = self.file
            self.info['nsamples']   = ((os.path.getsize(self.file)-self.hdrlen)
                                       /self.header['nchans'])*8/self.header['nbits']

        self.ctype = conf.nbits_to_ctypes[self.info["nbits"]]
        self.dtype = conf.ctypes_to_nptypes[self.ctype]
        self.info['ctype'] = self.ctype
        self.info['dtype'] = self.dtype


    def _reformat_radec(self,val):
        if val < 0:
            sign = -1
        else:
            sign = 1

        fractional,integral = math.modf(abs(val))
        xx = (integral-(integral%10000))/10000
        yy = ((integral-(integral%100))/100)-xx*100
        zz = integral - 100*yy - 10000*xx + fractional
        zz = string.zfill("%.4f"%(zz),7)
        return "%02d:%02d:%s"%(sign*xx,yy,zz)

    def _get_formatted_keys(self):
        combi_keys = self.header.keys()+self.info.keys()
        maxstrlen = len(max(combi_keys,key=len))+\
            len(max(conf.header_units.keys(), key=len))+3
        
        formatted_keys = {}
        
        for key in combi_keys:
            prstr = "%s %s"%(key,conf.header_units[key]) 
            prstr ="".join([prstr," "*(maxstrlen-len(prstr))]) 
            formatted_keys[key] = "%s:"%(prstr)
        self.formatted_keys = formatted_keys
                
    def npheader(self):
        dtypes = [(key,conf.struct_to_numpy[dtype]) for key,dtype\
                      in zip(conf.header_keys.keys(),conf.header_keys.values())\
                      if not dtype in [None,"b"]]

        array = numpy.empty(shape=[1],dtype=dtypes)
        for key in self.header.keys():
            array[key] = self.header[key]
        for key in self.info.keys():
            array[key] = self.info[key]
        return(array)


class MultiHeader:
    def __init__(self,files):
        if not isinstance(files,list):
            raise ValueError,"Bad value passed to MultiHeader.__init__()"

        self.headers = [Header(filename) for filename in files]
        self.info = {}
        ref = self.headers[0]
        for key in ref.info.keys():
            if self.testAttribute(key):
                self.info[key] = ref.info[key]
                
    def testAttribute(self,key,mode="same"):
        if mode == "same":
            return (len(set([h.info[key] for h in self.headers]))==1)
        elif mode == "different":
            return (len(set([h.info[key] for h in self.headers]))==len(self.headers))
            
            
        
            
            
        
                
