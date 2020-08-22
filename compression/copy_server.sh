#From build dir
#from="lgnecco@ssh.lamsade.dauphine.fr:/home/lamsade/lgnecco/stage/gedlib/compression/data/output/compressed/"
#to="/home/lucas/Documents/stage/gedlib/compression/data/from_server"
from="/home/lucas/Documents/stage/gedlib/compression/data/output/compressed/"
to="lgnecco@ssh.lamsade.dauphine.fr:/home/lamsade/lgnecco/stage/gedlib/compression/data/output/compressed"
echo -ne '\007'
eval "rsync -az -e ssh" $from $to
echo -ne '\007'


