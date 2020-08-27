#From build dir
from="lgnecco@ssh.lamsade.dauphine.fr:/home/lamsade/lgnecco/stage/gedlib/compression/data/output/compressed/results_compression_new.csv"
to="/home/lucas/Documents/stage/gedlib/compression/data/from_server/results_compression_new_1.csv"
#from="/home/lucas/Documents/stage/gedlib/compression/tests/"
#to="lgnecco@ssh.lamsade.dauphine.fr:/home/lamsade/lgnecco/stage/gedlib/compression/tests"
echo -ne '\007'
eval "rsync -az -e ssh" $from $to
echo -ne '\007'


