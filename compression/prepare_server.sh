#From build dir
base="/home/lamsade/lgnecco/stage/gedlib/compression"
echo -ne '\007'
eval "cd" $base
eval "rm -rf build"
eval "mkdir build"
eval "cd build"
echo "cmake"
eval "cmake .."
echo "Adapt to server"
eval "cd .."
eval "python3 adapt_to_server.py"
echo "make"
eval "cd build"
eval "make"
echo -ne '\007'


