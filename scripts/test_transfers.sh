

CONF0=$1
CONF1=$2

L=4
T=4

DIRHALOSIZE=$((L*L*L*2/2*9))
DIRSIZE=$((L*L*L*(T+4)/2*9 ))


diffcheck(){

    DIR=$1

    ENDLINE0=$((DIRSIZE*DIR+1+2*DIRHALOSIZE))
    ENDLINE1=$((DIRSIZE*(DIR+1)+1))


    echo $((ENDLINE0-DIRHALOSIZE))-$ENDLINE0 vs$((ENDLINE1-DIRHALOSIZE))-$ENDLINE1\
        $(paste <(head -$ENDLINE0 $CONF1 | tail -$DIRHALOSIZE) <(head -$ENDLINE1 $CONF0 | tail -$DIRHALOSIZE) | awk '{a+=($1-$3)^2 + ($2-$4)^2 }END{print a}')\
        $(paste <(head -$ENDLINE0 $CONF0 | tail -$DIRHALOSIZE) <(head -$ENDLINE1 $CONF1 | tail -$DIRHALOSIZE) | awk '{a+=($1-$3)^2 + ($2-$4)^2 }END{print a}')

    ENDLINE0=$((DIRSIZE*DIR+1+DIRHALOSIZE))
    ENDLINE1=$((DIRSIZE*(DIR+1)+1-DIRHALOSIZE))

    echo $((ENDLINE0-DIRHALOSIZE))-$ENDLINE0 vs$((ENDLINE1-DIRHALOSIZE))-$ENDLINE1\
        $(paste <(head -$ENDLINE0 $CONF1 | tail -$DIRHALOSIZE) <(head -$ENDLINE1 $CONF0 | tail -$DIRHALOSIZE) | awk '{a+=($1-$3)^2 + ($2-$4)^2 }END{print a}')\
        $(paste <(head -$ENDLINE0 $CONF0 | tail -$DIRHALOSIZE) <(head -$ENDLINE1 $CONF1 | tail -$DIRHALOSIZE) | awk '{a+=($1-$3)^2 + ($2-$4)^2 }END{print a}')





}

diffcheck 0
diffcheck 1
diffcheck 2
diffcheck 3
diffcheck 4
diffcheck 5
diffcheck 6
diffcheck 7



