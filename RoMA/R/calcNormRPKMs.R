calcNormRPKMs <-
function(e){

ScaleFactors <- e$lib.size
Count.normalized <- t(t(e$counts+0.5)/ScaleFactors) * mean(ScaleFactors)

GeneLen<-e$length

Count.normalized.lenadj<-Count.normalized/GeneLen
Count.normalized.total<-colSums(Count.normalized)
Count.normalized.lenadj<-Count.normalized.lenadj/median(Count.normalized.total)
Count.normalized.lenadj<-Count.normalized.lenadj*10^9

e$rpkms<-Count.normalized.lenadj

e
}
