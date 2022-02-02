# latest gcc image
FROM gcc:latest

# install dependencies
RUN apt-get update && apt-get install -y git cmake zip  

# make working folder
RUN mkdir -p /opt/code

# buil everything and clean up afterwards
RUN cd /opt/code && \ 
git clone https://github.com/microsoft/vcpkg vcpkg && \ 
vcpkg/bootstrap-vcpkg.sh -disableMetrics && \
git clone https://github.com/kherz/bmcppfit && \
cd bmcppfit && \
mkdir install && \
mkdir build && cd build && \
cmake ../ -DCMAKE_TOOLCHAIN_FILE=/opt/code/vcpkg/scripts/buildsystems/vcpkg.cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../install && \
make && \
make install && \
rm -rf /opt/code/vcpkg && \
rm -rf /opt/code/bmcppfit/build
