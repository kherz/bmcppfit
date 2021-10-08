# latest gcc image
FROM gcc:latest

# install dependencies
RUN apt-get update && apt-get install -y git cmake 

# make working folder
RUN mkdir -p /opt/code

#install vcpkg
RUN cd /opt/code && \ 
git clone https://github.com/microsoft/vcpkg vcpkg && \
vcpkg/bootstrap-vcpkg.sh -disableMetrics

# clone bmcppfit and prepare folders
RUN cd /opt/code && \ 
git clone https://github.com/kherz/bmcppfit && \
cd bmcppfit && \
mkdir build && \
mkdir install

# configure
RUN cd /opt/code/bmcppfit/build && \
cmake ../ -DCMAKE_TOOLCHAIN_FILE=/opt/code/vcpkg/scripts/buildsystems/vcpkg.cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../install

# build and install
RUN cd /opt/code/bmcppfit/build && \
make && \
make install 

#remove unnecessary depenencies
RUN rm -rf /opt/code/vcpkg/packages && \
rm -rf /opt/code/vcpkg/buildtrees && \
rm -rf /opt/code/vcpkg/downloads && \
rm -rf /opt/code/bmcppfit/build
