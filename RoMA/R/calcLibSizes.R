calcLibSizes <-
function(e){
size<-e$samples$lib.size*e$samples$norm.factors

e$lib.size<-size

e
}
