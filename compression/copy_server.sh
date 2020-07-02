#From build dir
from="/home/lucas/Documents/stage/"
to="lgnecco@ssh.lamsade.dauphine.fr:/home/lamsade/lgnecco/stage"
echo -ne '\007'
eval "rsync -az -e ssh" $from $to
echo -ne '\007'

