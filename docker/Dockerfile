FROM centos:centos7
MAINTAINER Remy Zyngfogel <remy@metocean.co.nz>

ADD hdf5-1.10.1.tar.gz /

RUN yum -y install epel-release &&\
    yum -y install git vim sudo tar unzip subversion wget &&\
    yum -y update &&\
    yum clean all

# install base packages through yum
RUN yum install -y expat-devel libtool tcl-devel libjpeg-devel libpng-devel gcc-gfortran cmake && \
    yum install -y gcc-c++ ruby flex bison libxslt tk tk-devel make file deltarpm

RUN yum -y install glibc-static

RUN yum -y groupinstall 'Development Tools'

RUN yum -y install libcurl-devel

# Get HDF5, NETCDF, NETDF-FORTRAN AND MPICH2 via ftp/http
RUN echo "-----------------GET Netcdf MPICH and HDF5"-----------------"" &&\ 
#    wget ftp://ftp.hdfgroup.org/HDF5/releases/hdf5-1.10/hdf5-1.10.1/src/hdf5-1.10.1.tar.gz && \
    wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.5.0.tar.gz && \
    wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-fortran-4.4.4.tar.gz && \
    wget http://www.mpich.org/static/downloads/3.2/mpich-3.2.tar.gz 



# Install HDF5
#RUN tar zxvf hdf5-1.10.1.tar.gz && \
RUN cd hdf5-1.10.1 && \
    ./configure --prefix=/usr/local \
                --enable-shared \
                --enable-hl \
                FC=gfortran && \
    make -j4 && \
    make install && \
    cd .. 

# Install NETCDF
RUN tar zxvf netcdf-4.5.0.tar.gz && \
    cd netcdf-4.5.0 && \
    ./configure --enable-netcdf-4 \
                --enable-dap \
                --enable-shared \
                --prefix=/usr/local \
                FC=gfortran \
                CPPFLAGS=-I/usr/local/include \
                LDFLAGS=-L/usr/local/lib && \
    make -j4 && \
    make install && \
    cd .. 
#    rm -rf /netcdf-4.5.0 /netcdf-4.5.0.tar.gz



# Install NETCDF-FORTRAN
RUN tar zxvf netcdf-fortran-4.4.4.tar.gz && \
    cd netcdf-fortran-4.4.4 && \
    export LD_LIBRARY_PATH=/usr/local/lib:${LD_LIBRARY_PATH} && \
    ./configure --prefix=/usr/local \
                FC=gfortran \
                CPPFLAGS=-I/usr/local/include \
                 LDFLAGS=-L/usr/local/lib && \
    make -j4 && \
    make install && \
    cd .. 
#    rm -rf /netcdf-fortran-4.4.3 /netcdf-fortran-4.4.3.tar.gz
#
## Install MPICH2
RUN tar zxvf mpich-3.2.tar.gz && \
    cd mpich-3.2 && \
    ./configure --prefix=/usr/local \
                --enable-shared \
                --enable-lib-depend \
                --enable-sharedlibs=gcc \
                F77=gfortran \
                FC=gfortran \
                CC=gcc \
                CXX=g++ \
                FCFLAGS=-fPIC \
                CFLAGS=-fPIC \
                CCFLAGS=-fPIC \
                CXXFLAGS=-fPIC && \
    make -j4 && \
    make install && \
    cd .. && \
    rm -rf /mpich-3.2 /mpich-3.2.tar.gz

RUN cp /usr/lib/gcc/x86_64-redhat-linux/4.8.2/libquadmath.so /usr/local/lib/libquadmath.a

RUN useradd -ms /bin/bash user



## Downloading schism
RUN cd /home/user/ && \
    svn co http://columbia.vims.edu/schism/tags/v5.5.0

## edit the makefile
RUN cd /home/user/v5.5.0/mk &&\
    cp Make.defs.centos.gnu.ompi Make.defs.local &&\
    sed -i "s|ENV = CentOS.gnu.ompi|ENV = HYDRO|g" Make.defs.local && \
    sed -i "s|FCS = gfortran -ffree-line-length-none|FCS = /usr/local/bin/mpif90 -f90=gfortran -ffree-line-length-none|g" Make.defs.local &&\
    sed -i "s|FCP = mpif90 -ffree-line-length-none|FCP = /usr/local/bin/mpif90 -ffree-line-length-none|g" Make.defs.local &&\
    sed -i "s|#EXEC   := othername.ex|EXEC   := pschism_5.5.0_$(ENV)|g" Make.defs.local

## Compile SCHISM
RUN cd /home/user/v5.5.0/src/ &&\
    make clean &&\
    make &&\
    sudo cp /home/user/v5.5.0/src/pschism_5.5.0__VL /usr/local/bin/schism


## Compile combine code
RUN cd /home/user/v5.5.0/src/Utility/Combining_Scripts/ &&\
    gfortran -ffree-line-length-none -cpp -O2 -o combine_output10 ../UtilLib/argparse.F90 combine_output10.f90 ../UtilLib/schism_geometry.f90 -I/usr/local/include -L/usr/local/lib -lnetcdf -lnetcdff &&\
    sudo cp combine_output10 /usr/local/bin/


## Download git lfs 
RUN cd /home/user/ &&\
    wget 'https://github.com/git-lfs/git-lfs/releases/download/v2.3.4/git-lfs-linux-amd64-2.3.4.tar.gz' &&\
    tar zxvf git-lfs-linux-amd64-2.3.4.tar.gz

USER root
    
# install git-lfs
RUN cd /home/user/git-lfs-2.3.4 &&\
    sudo ./install.sh &&\
    git lfs install

# download all the code from github
RUN cd /home/user &&\
#    git lfs clone git://github.com/metocean/SSC.git
git lfs clone https://github.com/metocean/SSC.git

## unzip initial files
RUN cd /home/user/SSC/initial_files &&\
    unzip input_files.zip &&\
    mv input_files/* /home/user/SSC/initial_files &&\
    rm /home/user/SSC/initial_files/input_files.zip &&\
    rm -rf /home/user/SSC/initial_files/input_files

   
## ADD libspatial for Rtree package
RUN cd /home/user/ &&\
    wget 'http://download.osgeo.org/libspatialindex/spatialindex-src-1.8.5.tar.gz' &&\
    tar zxvf spatialindex-src-1.8.5.tar.gz &&\
    cd /home/user/spatialindex-src-1.8.5/ &&\
    ./configure &&\
    make &&\
    sudo make install

ENV LD_LIBRARY_PATH=/usr/local/lib/:$LD_LIBRARY_PATH


## ADD more python tools
RUN sudo yum -y install python-pip python-devel

## ADD python lib
RUN sudo pip install netCDF4 pyyaml mako numpy==1.8.2 rtree pyproj scipy

USER user
