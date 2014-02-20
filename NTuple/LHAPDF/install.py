import os
import sys
import commands
import glob
import multiprocessing

YELLOW = '\x1b[33m'
RED = '\x1b[31m'
PINK = '\x1b[35m'
NORMAL = '\x1b[0m'
CYAN = '\x1b[36m'

ncores = 4

def CompileLHAPDF():

    # YAML-CPP
    print YELLOW+"Building yaml-cpp for LHAPDF ..."+NORMAL
    print " - Configuring ..."
    os.system("cd install/yaml-cpp/build/ && cmake -DBUILD_SHARED_LIBS=ON ..")
    print " - Compiling ..."
    os.system("cd install/yaml-cpp/build/ && make")
    print " - Moving the headers folder..."
    os.system("cd install/yaml-cpp/ && cp -rv include/* ../../include/")
    print " - Moving the dynamic library ..."
    os.system("cd install/yaml-cpp/build/ && cp *.so ../../../lib/")
    print " - Checking the presence of the library 'libyaml-cpp.so' ..."
    if os.path.isfile('lib/libyaml-cpp.so'):
        print CYAN+" -> OK"+NORMAL
    else:
        print RED +" -> BAD"+NORMAL
        print "Stop now the recipe!"
        sys.exit()
    print ""

    # LHAPDF
    print YELLOW+"Building LHAPDF6 ..."+NORMAL
    print " - Configuring ..."
    mycmd = './configure --prefix='+os.getcwd()+'/'
    mycmd += ' --disable-lhaglue --disable-python'
    mycmd += ' --with-boost=/opt/exp_soft/cms/slc5_amd64_gcc462/external/boost/1.50.0-cms'
    mycmd += ' --with-yaml-cpp-lib='+os.getcwd()+'/lib/'
    mycmd += ' --with-yaml-cpp-inc='+os.getcwd()+'/include/'
    print mycmd
    os.system("cd install/LHAPDF-6.0.2 && "+mycmd)
    print " - Compiling ..."
    os.system("cd install/LHAPDF-6.0.2 && make -j"+str(ncores))
    print " - Moving the dynamic library and the headers ..."
    os.system("cd install/LHAPDF-6.0.2 && make install")
    print " - Checking the presence of the library 'libyaml-cpp.so' ..."
    if os.path.isfile('lib/libLHAPDF.so'):
        print CYAN+" -> OK"+NORMAL
    else:
        print RED +" -> BAD"+NORMAL
        print "Stop now the recipe!"
        sys.exit()
    print ""

def InstallPDF(pdf):

    print YELLOW+"Downloading the PDF sets ..."+NORMAL
    begin = "cd share/LHAPDF && " +\
            "wget --no-check-certificate "
    for item in pdf:
        if not item.endswith('.tar.gz'):
            print RED+"skip the file with bad extension: "+item+NORMAL
            continue
        os.system(begin+'"'+item+'"')
        words = item.split('/')
        if len(words)==0:
            continue
        if not words[-1].endswith('.tar.gz'):
            print RED+"skip the file with bad extension: "+words[-1]+NORMAL
            continue
        filename = words[-1]
        os.system("cd share/LHAPDF && tar xvzf "+filename)
        os.system("cd share/LHAPDF && rm -f "+filename)
        print ""



# Compiling LHAPDF
CompileLHAPDF()

# Downloading PDF sets
print YELLOW+"Download PDF sets ..."+NORMAL
pdf = [ 'http://www.hepforge.org/archive/lhapdf/pdfsets/6.0/cteq6l1.tar.gz',\
        'http://www.hepforge.org/archive/lhapdf/pdfsets/6.0/unvalidated/cteq66.tar.gz',\
        'http://www.hepforge.org/archive/lhapdf/pdfsets/6.0/unvalidated/MSTW2008lo68cl.tar.gz',\
        'http://www.hepforge.org/archive/lhapdf/pdfsets/6.0/unvalidated/MSTW2008nlo68cl_asmz+68cl.tar.gz',\
        'http://www.hepforge.org/archive/lhapdf/pdfsets/6.0/unvalidated/MSTW2008nlo68cl_asmz-68cl.tar.gz',\
        'http://www.hepforge.org/archive/lhapdf/pdfsets/6.0/NNPDF23_nlo_as_0119.tar.gz',\
        'http://www.hepforge.org/archive/lhapdf/pdfsets/6.0/NNPDF23_nlo_as_0118.tar.gz',\
        'http://www.hepforge.org/archive/lhapdf/pdfsets/6.0/NNPDF23_nlo_as_0120.tar.gz']
InstallPDF(pdf)















