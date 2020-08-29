#From build dir
#from="lgnecco@ssh.lamsade.dauphine.fr:/home/lamsade/lgnecco/stage/gedlib/compression/data/output/compressed/results_refinement_2.csv"
#to="/home/lucas/Documents/stage/gedlib/compression/data/from_server/results_refinement_server_1.csv"

from="/home/lucas/Documents/stage/gedlib/compression/src/"
to="lgnecco@ssh.lamsade.dauphine.fr:/home/lamsade/lgnecco/stage/gedlib/compression/src"

#from="/home/lucas/Documents/stage/gedlib/run.sh"
#to="lgnecco@ssh.lamsade.dauphine.fr:/home/lamsade/lgnecco/stage/gedlib/run.sh"

echo -ne '\007'
eval "rsync -az -e ssh" $from $to
echo -ne '\007'

from="/home/lucas/Documents/stage/gedlib/compression/tests/"
to="lgnecco@ssh.lamsade.dauphine.fr:/home/lamsade/lgnecco/stage/gedlib/compression/tests"

echo -ne '\007'
eval "rsync -az -e ssh" $from $to
echo -ne '\007'
