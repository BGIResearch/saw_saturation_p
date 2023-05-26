
# change the libPath to your specific path
libPath="/root/lib"
python3="python3"


cmakePath="$libPath/cmake-3.17.2"
gccPath="$libPath/gcc-9.1.0"
CLI11Path="$libPath/CLI11-1.9.0"
hdf5Path="$libPath/hdf5-1.12.1"
libxPath="$libPath/libx-1.1"


binPath="/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/bin"
export PATH="$gccPath/bin:$cmakePath/bin:$binPath"

#export LD_LIBRARY_PATH="$gccPath/lib64:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH=$gccPath/lib64:$hdf5Path/lib
export LIBRARY_PATH=$LD_LIBRARY_PATH

export C_INCLUDE_PATH="$libxPath/include:$hdf5Path/include:$CLI11Path/include:$gccPath/include:$C_INCLUDE_PATH"
export CPLUS_INCLUDE_PATH=$C_INCLUDE_PATH

export CC="$gccPath/bin/gcc"
export CXX="$gccPath/bin/g++"

absPath(){
    relativePath=$1
    mkdir -p $relativePath
    if [ ${relativePath:0:1} == "/" ]
    then
        echo $relativePath
        return 0
    fi
    ap="$(cd $relativePath; pwd)"
    echo "$ap"
}

install(){
    libFile="$1"
    tlib="$installPath/lib"
    if [ -d "$tlib" ]
    then
        if [ -e "$libFile" ]
        then
            echo "Installing: $libFile into $tlib"
            cp $libFile $tlib
        else
            echo "[WARN] File not exists: $libFile"
        fi
    fi
}

srcPath="$(cd $(dirname $(dirname $0)); pwd)"
buildPath="build"
buildPath="$(absPath $buildPath)"
mkdir -p $buildPath
cd $buildPath
installPath="$buildPath/install"
mkdir -p $installPath/bin $installPath/lib

timeStart=$(date +%s)

cmake $srcPath -DINSTALL_PATH=$installPath -DHDF5PATH=$hdf5Path/lib

thread=$(grep -c ^processor /proc/cpuinfo)
make -j $thread #VERBOSE=1

if [[ $? == 0 ]];
then
    install $gccPath/lib64/libstdc++.so.6
    install $gccPath/lib64/libgomp.so.1
    install $gccPath/lib64/libgcc_s.so.1
    install $hdf5Path/lib/libhdf5.so.200
    # install $hdf5Path/lib/libhdf5_cpp.so.200

    find $srcPath/src -name "*.py"| xargs -I {} $python3 -c "import py_compile;py_compile.compile('{}',cfile='{}c')"
    cp $srcPath/src/plot.pyc $installPath/bin/

    cd $srcPath
    cp -R $installPath $srcPath
    echo -n "Success! "
else
    echo -n "Fail! "
fi

timeEnd=$(date +%s)
secs=$(($timeEnd - $timeStart))
mins=$(($secs/60))
secs=$(($secs%60))
echo "Cost: $mins Mins $secs Secs"
