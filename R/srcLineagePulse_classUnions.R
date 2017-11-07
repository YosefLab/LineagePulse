###############################################################################
#####################     Class unions for S4 slots     #######################
###############################################################################

# Define class unions for slots
setClassUnion('numericORNULL', members=c('numeric', 'NULL'))
setClassUnion('matrixORNULL', members=c('matrix', 'NULL'))
setClassUnion('characterORNULL', members=c('character', 'NULL'))
setClassUnion('listORNULL', members=c('list', 'NULL'))
setClassUnion('data.frameORNULL', members=c('data.frame', 'NULL'))
