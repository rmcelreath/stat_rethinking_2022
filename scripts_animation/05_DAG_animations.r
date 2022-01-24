library(rethinking)
library(animation)
library(dagitty)
library(plotrix)


# ani.saveqz(dpi=120)

# convert -alpha remove -background white -delay 10 -loop 0 frame*.png exdag.gif

# convert -delay 10 exdag.gif exdag.gif

####################################
# TESTS
if ( FALSE ) {

##########
# FORK EXAMPLE

exdag <- dagitty( "dag {
    X <- Z -> Y
}")
coordinates( exdag ) <- list( x=c(X=0,Z=1,Y=2) , y=c(X=0,Z=1,Y=0) )

SCM <- list(
    list( var="Z" , f=function(sv) rbern(1,0.5) ),
    list( var="X" , f=function(sv) rbern(1, inv_logit(2*sv['Z']-1) ) ),
    list( var="Y" , f=function(sv) rbern(1, inv_logit(2*sv['Z']-1) ) )
)

dagfx_anim_forward( exdag , Y="" , fade_rate=1 , n_loops=5 , path_frames=50 , color_list=list(X=2,Z=4) , buffer=c(0,0.4,0,0) , arc_lwd=15 , angle_gap=pi/4 , angle1_inc=pi/12 , accel=0.5 , vel0=0 , fix="" , scm=SCM )

oopts = ani.options(interval = 0.05)
ani.replay()

##########
# FORK EXAMPLE with NOISE sources

exdag <- dagitty( "dag {
    X <- Z -> Y
    e_X -> X
    e_Y -> Y
}")
coordinates( exdag ) <- list( x=c(X=0,Z=1,Y=2,e_X=-1,e_Y=3) , y=c(X=0,Z=1,Y=0,e_X=0,e_Y=0) )

SCM <- list(
    list( var="Z" , f=function(sv) rbern(1,0.5) ),
    list( var="e_X" , f=function(sv) rbern(1,0.5) ),
    list( var="e_Y" , f=function(sv) rbern(1,0.5) ),
    list( var="X" , f=function(sv) rbern(1, inv_logit(2*(sv['Z']-sv['e_X'])-1) ) ),
    list( var="Y" , f=function(sv) rbern(1, inv_logit(2*(sv['Z']-sv['e_Y'])-1) ) )
)

dagfx_anim_forward( exdag , Y="" , fade_rate=1 , n_loops=10 , path_frames=20 , color_list=list(X=2,Z=4,e_Y=5,e_X=2) , buffer=c(0,0.4,0,0) , arc_lwd=15 , angle_gap=pi/4 , angle1_inc=pi/12 , accel=0.5 , vel0=0 , fix="" , scm=SCM )

oopts = ani.options(interval = 0.05)
ani.replay()


##########
# PIPE EXAMPLE 

exdag <- dagitty( "dag {
    X -> Z -> Y
}")
coordinates( exdag ) <- list( x=c(X=0,Z=1,Y=2,e_X=-1,e_Y=3) , y=c(X=0,Z=1,Y=0,e_X=0,e_Y=0) )

SCM <- list(
    list( var="X" , f=function(sv) rbern(1, 0.5 ) ),
    list( var="Z" , f=function(sv) rbern(1, inv_logit(2*sv['X']-1) ) ),
    list( var="Y" , f=function(sv) rbern(1, inv_logit(2*sv['Z']-1) ) )
)

dagfx_anim_forward( exdag , Y="" , fade_rate=1 , n_loops=10 , path_frames=20 , color_list=list(X=2,Z=4,e_Y=5,e_X="pink") , buffer=c(0,0.4,0,0) , arc_lwd=15 , angle_gap=pi/4 , angle1_inc=pi/12 , accel=0.5 , vel0=0 , fix="" , scm=SCM )

oopts = ani.options(interval = 0.05)
ani.replay()

##########
# COLLIDER EXAMPLE 

exdag <- dagitty( "dag {
    X -> Z <- Y
}")
coordinates( exdag ) <- list( x=c(X=0,Z=1,Y=2) , y=c(X=0,Z=1,Y=0) )

SCM <- list(
    list( var="X" , f=function(sv) rbern(1, 0.5 ) ),
    list( var="Y" , f=function(sv) rbern(1, 0.5 ) ),
    list( var="Z" , f=function(sv) rbern(1, inv_logit(2*(sv['X']+sv['Y'])-2) ) )
)

dagfx_anim_forward( exdag , Y="" , fade_rate=1 , n_loops=10 , path_frames=20 , color_list=list(X=2,Z=4,Y=4) , buffer=c(0,0.4,0,0) , arc_lwd=15 , angle_gap=pi/4 , angle1_inc=pi/12 , accel=0.5 , vel0=0 , fix="" , scm=SCM )

oopts = ani.options(interval = 0.05)
ani.replay()

####################################
# instrument
exdag <- dagitty( "dag {
    U [unobserved]
    Z -> X -> Y
    X <- U -> Y
}")
coordinates( exdag ) <- list( x=c(Z=0,X=1,Y=2,U=1.5) , y=c(Z=0,X=0,Y=0,U=-1) )

SCM <- list(
    list( var="U" , f=function(sv) rbern(1,0.5) ),
    list( var="Z" , f=function(sv) rbern(1,0.5) ),
    list( var="X" , f=function(sv) rbern(1, inv_logit(sv['Z']-sv['U']) ) ),
    list( var="Y" , f=function(sv) rbern(1, inv_logit(sv['X']-sv['U']) ) )
)

test <- dagfx_anim_forward( exdag , Y="Y" , fade_rate=1 , n_loops=10 , path_frames=20 , color_list=list(X=2,U=5,Z=3,Y=1,A=4,B=5,V=4,M=3,H=2,S=3) , buffer=c(0,0.4,0,0) , arc_lwd=15 , angle_gap=pi/4 , angle1_inc=pi/12 , accel=0.5 , vel0=0 , fix="" , scm=SCM )

oopts = ani.options(interval = 0.05)
ani.replay()

####################################
# m-bias
exdag <- dagitty( "dag {
    X -> Y
    X <- A -> Z
    Y <- B -> Z
}")
coordinates( exdag ) <- list( x=c(X=0,A=0,Z=1,Y=2,B=2) , y=c(X=0,A=-1,Z=-0.7,B=-1,Y=0) )

##########
# height, weight ,sex

exdag <- dagitty( "dag {
    S -> H
    H -> W
    S -> W
}")
coordinates( exdag ) <- list( x=c(H=0,S=0,W=1) , y=c(H=0,S=1,W=0) )

dagfx_anim_forward( exdag , Y="W" , fade_rate=1 , n_loops=10 , path_frames=50 , color_list=list(X=2,U=5,Z=3,Y=1,A=4,B=5,V=4,M=3,H=2,S=4) , buffer=c(0,0.4,0,0) , arc_lwd=15 , angle_gap=pi/4 , angle1_inc=pi/12 , accel=0.5 , vel0=0 , fix="" , scm=NULL )

oopts = ani.options(interval = 0.05)
ani.replay()

####################################
# mediator

exdag <- dagitty( "dag {
    X -> Y
    X -> M -> Y
    M <- U -> Y
}")
coordinates( exdag ) <- list( x=c(X=0,Y=2,M=1,U=2) , y=c(X=0,Y=0,M=-1,U=-1) )

dagfx_anim_forward( exdag , Y="Y" , fade_rate=1 , n_loops=4 , path_frames=50 , color_list=list(X=2,U=5,Z=3,Y=1,A=4,B=5,V=4,M=3,H=2,S=3) , buffer=c(0,0.4,0,0) , arc_lwd=15 , angle_gap=pi/4 , angle1_inc=pi/12 , accel=0.5 , vel0=0 , fix="" , scm=NULL )

oopts = ani.options(interval = 0.05)
ani.replay()


####################################
# mediator 2 confounds
exdag <- dagitty( "dag {
    X -> Y
    X -> M -> Y
    M <- U -> Y
    X <- V -> M
}")
coordinates( exdag ) <- list( x=c(X=0,Y=2,M=1,U=2,V=0) , y=c(X=0,Y=0,M=-1,U=-1,V=-1) )

do_anim( exdag , Y="Y" , fade_rate=1 , n_loops=4 , path_frames=50 , color_list=list(X=2,U=5,Z=3,Y=1,A=4,B=5,V=4,M=3,H=2,S=3) , buffer=c(0,0.4,0,0) , arc_lwd=15 , angle_gap=pi/4 , angle1_inc=pi/12 , accel=0.5 , vel0=0 , fix="M" )

oopts = ani.options(interval = 0.1)
ani.replay()

####################################
# table 2 fallacy smoking DAG
exdag <- dagitty( "dag {
    X -> Y
    S -> Y
    S -> X
    S <- U -> Y
    A -> S
    A -> X
    A -> Y
}")
coordinates( exdag ) <- list( x=c(A=0,X=1,S=0,U=1,Y=2) , y=c(A=0,X=1,S=2,Y=1,U=2.5) )

do_anim( exdag , Y="Y" , fade_rate=1 , n_loops=5 , path_frames=50 , color_list=list(X=2,U=5,Z=3,Y=1,A=4,B=5,V=4,M=3,H=2,S=3) , buffer=c(0,0.4,0,0) , arc_lwd=15 , angle_gap=pi/4 , angle1_inc=pi/12 , accel=0.5 , vel0=0 , fix=c("A","S") , force_on="X" )

oopts = ani.options(interval = 0.07)
ani.replay()

}
