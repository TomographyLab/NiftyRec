#!/usr/bin/python

import sys
from os.path import walk as walk
from os import listdir as listdir
from os.path import isdir as isdir
from os.path import isfile as isfile




class Counter:
    def __init__(self, exclude_dir=None):
        self.directory = ""
        self.exclude_dir = exclude_dir
        self.clear()

    def process_dir(self,arg,dirname,fnames):
        if self.exclude_dir:
            for excl_dir in self.exclude_dir:
                if dirname.startswith(self.directory+excl_dir):
                    return 
        self.dirs.append(dirname)
        self.n_dirs += 1
        for f in fnames:
            the_file = dirname+"/"+f
            #print the_file
            if isfile(the_file):
                if f.endswith(arg):
                    #print dirname+"/"+f
                    self.files.append(the_file)
                    self.n_files += 1
                    self.n_lines += self.count_lines_in_file(the_file)

    def count_lines_in_directory_tree(self,directory,extension):
        self.directory = directory
        self.clear()
        walk(directory,self.process_dir,extension)
        n_dirs  = self.n_dirs
        n_files = self.n_files
        n_lines = self.n_lines
        self.clear()
        return n_dirs, n_files, n_lines

    def count_lines_in_directory(self,directory,extension):
        n_lines = 0
        n_files = 0
        for a in listdir(directory):
            the_file = directory+"/"+a
            if isfile(the_file):
                if the_file.endswith(extension):
                    n_files += 1
                    n_lines += self.count_lines_in_file(the_file)
        return n_files, n_lines
                

    def count_lines_in_file(self,the_file):
        f = open(the_file,"r")
        n_lines = 0
        while 1:
            d = f.readline()
            if not d:
                break
            n_lines += 1
        f.close()
        return n_lines

    def clear(self):
        self.n_dirs = 0
        self.n_files = 0
        self.n_lines = 0
        self.files = []
        self.dirs = []



if __name__ == "__main__":
    directory = sys.argv[1]
    exclude_dir = []
    if len(sys.argv) >= 3:
        for i in range(len(sys.argv)-2):
            exclude_dir.append(sys.argv[2+i])
    C = Counter(exclude_dir)
    py_d,  py_f,  py_l  = C.count_lines_in_directory_tree(directory,".py")
    c_d,   c_f,   c_l   = C.count_lines_in_directory_tree(directory,".c")
    cpp_d, cpp_f, cpp_l = C.count_lines_in_directory_tree(directory,".cpp")
    cu_d,  cu_f,  cu_l  = C.count_lines_in_directory_tree(directory,".cu")
    h_d, h_f, h_l       = C.count_lines_in_directory_tree(directory,".h")
    hpp_d, hpp_f, hpp_l = C.count_lines_in_directory_tree(directory,".hpp")
    gl_d, gl_f, gl_l    = C.count_lines_in_directory_tree(directory,".glade")
    xml_d, xml_f, xml_l = C.count_lines_in_directory_tree(directory,".xml")
    m_d, m_f, m_l       = C.count_lines_in_directory_tree(directory,".m")
#    mex_d, mex_f, mex_l = C.count_lines_in_directory_tree(directory,"mex.cpp")
#    meh_d, meh_f, meh_l = C.count_lines_in_directory_tree(directory,"mex.h"); mex_f+=meh_f; mex_l+=meh_l;
    sh_d, sh_f, sh_l    = C.count_lines_in_directory_tree(directory,".sh")
    iss_d, iss_f, iss_l = C.count_lines_in_directory_tree(directory,".iss")
    mk_d, mk_f, mk_l    = C.count_lines_in_directory_tree(directory,"makefile")
    cm_d, cm_f, cm_l    = C.count_lines_in_directory_tree(directory,"CMakeLists.txt")

    print "----------------------------------------------------------------------"   
    print "Lines of code in directory",directory,"(including subdirectories)"
    print "Excluding ",exclude_dir
    print ""
    print "   LANGUAGE   SUBDIRS    FILES    LINES"
    print "" 
    print "   C:            %-8d %-8d %-8d" %(c_d,c_f,c_l)
    print "   h:            %-8d %-8d %-8d" %(h_d,h_f,h_l)
    print "   C++:          %-8d %-8d %-8d" %(cpp_d,cpp_f,cpp_l)
    print "   h++:          %-8d %-8d %-8d" %(hpp_d,hpp_f,hpp_l)
    print "   Cuda:         %-8d %-8d %-8d" %(cu_d,cu_f,cu_l)
#    print "   mex:          %-8d %-8d %-8d" %(mex_d,mex_f,mex_l)
    print "   Matlab:       %-8d %-8d %-8d" %(m_d,m_f,m_l)
    print "   Python:       %-8d %-8d %-8d" %(py_d,py_f,py_l)
    print "   Glade:        %-8d %-8d %-8d" %(gl_d,gl_f,gl_l)
    print "   XML:          %-8d %-8d %-8d" %(xml_d,xml_f,xml_l)
    print "   CMake:        %-8d %-8d %-8d" %(cm_d,cm_f,cm_l)
    print "   Shell script: %-8d %-8d %-8d" %(sh_d,sh_f,sh_l)
    print "   Makefiles:    %-8d %-8d %-8d" %(mk_d,mk_f,mk_l)
    print "   Inno Setup:   %-8d %-8d %-8d" %(iss_d,iss_f,iss_l)
    print "   TOTAL:                          %-8d" %(py_l+c_l+h_l+cpp_l+hpp_l+cu_l+gl_l+xml_l+m_l+sh_l+iss_l+mk_l+cm_l)
    print ""
    print "----------------------------------------------------------------------"   
