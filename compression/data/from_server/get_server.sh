#From build dir
from="lgnecco@ssh.lamsade.dauphine.fr:/home/lamsade/lgnecco/stage/gedlib/compression/data/output/"
to="/home/lucas/Documents/stage/gedlib/compression/data/from_server/"
echo -ne '\007'
eval "rsync -az -e ssh" $from $to
echo -ne '\007'


