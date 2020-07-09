#From build dir
from="/home/lucas/Documents/stage/gedlib/compression/"
to="lgnecco@ssh.lamsade.dauphine.fr:/home/lamsade/lgnecco/stage/gedlib/compression"
echo -ne '\007'
eval "rsync -az -e ssh" $from $to
echo -ne '\007'


