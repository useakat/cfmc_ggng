OS 10.5: 09/01/06

For OS 10.5, some machines give the following error message on executing td_mac_intel (td_g4 will work using the rosetta emulator which is already install on the intel macs.)

"dyld: Library not loaded: /usr/local/lib/libg2c.0.dylib"

Whether you get this error or not depends on exactly how you have installed gcc and g77.  There is a simple work around, install the file 
libg2c.0.dylib
which is included in the tarball in /usr/local/lib/  as follows

sudo cp libg2c.0.dylib /usr/local/lib/

(the adminstrator password will be required)
If the directory /usr/local/lib/ doesn't exist it will need to be created first using sudo mkdir.

With this work around, td will work even if gcc and g77 are not installed on your machine. 

notes:
The origin of this error is associated with dynamic verses static link when td is compiled.  Someday, I'll figure out a more elegant solution.

The file libg2c.0.dylib comes from 
http://prdownloads.sourceforge.net/hpc/g77-intel-bin.tar.gz?download
see /usr/local/lib/ in this tarball.

