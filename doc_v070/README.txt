
General
-------
DFDC and its plot library should compile on any Unix system 
with normal Fortran-77, C, and X-Windows support.

DFDC can also be compiled on Windows systems (Win32) using a Win32 
version of the Unix plot library.  This version works on Win32 but is 
not a native "Windows" program (i.e. limited to console operation with
a single graphic window for display). 


Unix build sequence
--------------

1) Build Xplot11 library in  ./plotlib  ...

 % cd plotlib
 % edit config.make (or choose/edit config.make.xx and copy to config.make)
 % make           (creates libPlt.a)


2) Build dfdc in ./bin  ...

 % cd bin
 % edit Makefile.xxx  (set compiler flags for your system)
 % make 


The executable will appear in the ./bin directory.

NOTE: two versions of userio.f exist in the src-dll directory.  The issue is
that g77 objects to certain format statements in the IO in userio2.f.  
See the file src-dll/userio.readme for details.


Win32 build sequence
--------------

1) Build Xplot11 library in  ./plotlib  ...

 % cd plotlib\win32
 % edit Makefile.NT       (set compiler flags for your system)
 % nmake -f Makefile.NT   (creates libPlt.lib)


2) Build dfdc in .\bin  ...

 % cd bin
 % edit Makefile.xxx  (use version xxx closest to your system and edit)
 % make -f Makefile.xxx

3) Build dfdcdll.dll in .\src-dll  ...

 % cd src-dll
 % edit Makefile.NT (to match your system compilers, etc)
 % make -f Makefile.NT


The test executable and dll library will appear in the .\src-dll directory.


Documentation
-------------
User Guide is in the  dfdc_doc.txt  file.  If impatient, you can just
run DFDC in the runs/ directory, which contains a few input files:

 % cd runs
 % ../bin/dfdc sample1.case

The files  session1.txt, session2.txt  contain keyboard input examples.



