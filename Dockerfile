FROM ubuntu:hirsute
RUN apt update
RUN apt upgrade -y
RUN DEBIAN_FRONTEND=noninteractive apt install -y cmake git g++ doxygen graphviz bash rsync curl libblas-dev liblapack-dev libeigen3-dev liblapacke-dev ninja-build