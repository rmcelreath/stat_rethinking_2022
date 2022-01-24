library(rethinking)
library(animation)
library(dagitty)
library(plotrix)


# ani.saveqz(dpi=150)

# convert -alpha remove -background white -delay 10 -loop 0 frame*.png exdag.gif

# convert -delay 10 exdag.gif exdag.gif

####################################
# TESTS
if ( FALSE ) {

####################################
# complex example

exdag <- dagitty( "dag {
    X -> Y
    X <- C -> Y
    X <- Z -> Y
    X <- A -> Z
    Y <- B -> Z
}")
coordinates( exdag ) <- list( 
    x=c(X=0,Y=2,C=1,Z=1,A=0,B=2) , 
    y=c(X=0,Y=0,C=1,Z=-1,A=-1,B=-1) )

SCM <- list(
    list( var="C" , f=function(sv) rbern(1,0.5) ),
    list( var="A" , f=function(sv) rbern(1,0.5) ),
    list( var="B" , f=function(sv) rbern(1,0.5) ),
    list( var="Z" , f=function(sv) rbern(1, inv_logit( 2*sv['A'] - 2*sv['B'] ) ) ),
    list( var="X" , f=function(sv) rbern(1, inv_logit( 2*sv['A'] - 2*sv['Z'] + 2*sv['C'] ) ) ),
    list( var="Y" , f=function(sv) rbern(1, inv_logit( sv['Z'] + sv['X'] - sv['C'] - sv['B'] ) ) )
)

dagfx_anim_forward( exdag , Y="Y" , fade_rate=1 , n_loops=7 , path_frames=30 , color_list=list(X=2,C=4,Z=5,A=3,B="white") , buffer=c(0,0.4,0,0) , arc_lwd=15 , angle_gap=pi/4 , angle1_inc=pi/12 , accel=0.5 , vel0=0 , fix="M" , scm=SCM )

oopts = ani.options(interval = 0.1)
ani.replay()

####################################
# table 2 fallacy smoking DAG

# without U
exdag <- dagitty( "dag {
    X -> Y
    S -> Y
    S -> X
    A -> S
    A -> X
    A -> Y
}")
coordinates( exdag ) <- list( x=c(A=0,X=1,S=0,Y=2) , y=c(A=0,X=-1,S=-2,Y=-1) )

SCM <- list(
    list( var="A" , f=function(sv) rbern(1,0.5) ),
    list( var="S" , f=function(sv) rbern(1, inv_logit( 2*sv['A'] ) ) ),
    list( var="X" , f=function(sv) rbern(1, inv_logit( 2*sv['A'] - 2*sv['S'] ) ) ),
    list( var="Y" , f=function(sv) rbern(1, inv_logit( sv['X'] + sv['A'] - sv['S'] ) ) )
)

# with U
exdag <- dagitty( "dag {
    X -> Y
    S -> Y
    S -> X
    S <- U -> Y
    A -> S
    A -> X
    A -> Y
}")
coordinates( exdag ) <- list( x=c(A=0,X=1,S=0,U=1,Y=2) , y=c(A=0,X=-1,S=-2,Y=-1,U=-2.5) )

SCM <- list(
    list( var="A" , f=function(sv) rbern(1,0.5) ),
    list( var="U" , f=function(sv) rbern(1,0.5) ),
    list( var="S" , f=function(sv) rbern(1, inv_logit( 2*sv['A'] - 2*sv['U'] ) ) ),
    list( var="X" , f=function(sv) rbern(1, inv_logit( 2*sv['A'] - 2*sv['S'] ) ) ),
    list( var="Y" , f=function(sv) rbern(1, inv_logit( sv['X'] + sv['A'] - sv['S'] - sv['U'] ) ) )
)

dagfx_anim_forward( exdag , Y="Y" , fade_rate=1 , n_loops=5 , path_frames=50 , color_list=list(X=2,U=5,Z=3,Y=1,A=4,B=5,V=4,M=3,H=2,S=3) , buffer=c(0,0.4,0,0) , arc_lwd=15 , angle_gap=pi/4 , angle1_inc=pi/12 , accel=0.5 , vel0=0 , fix="" , force_on="X" , scm=SCM )

oopts = ani.options(interval = 0.07)
ani.replay()

}

