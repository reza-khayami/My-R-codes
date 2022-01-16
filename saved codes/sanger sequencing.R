f2 <- readsangerseq("f2/166_Senseprimer.ab1")
m <- readsangerseq("m/187_Senseprimer.ab1")
p <- readsangerseq("p/186_Senseprimer.ab1")
p2 <- readsangerseq("p/186_Antisenseprimer.ab1")
ref <- readDNAStringSet("b.fasta")
ref <- toString(ref)
ref <- subseq(primarySeq(m, string=TRUE), start=30, width=600)

chromatogram(p2, width=200, height=2, trim5=100, trim3=100, 
             showcalls='both', filename="p2.pdf")
hetcallsf2 <- makeBaseCalls(f2, ratio=0.33)

chromatogram(hetcallsm, width=100, height=2, trim5=100, trim3=100, 
             showcalls='both', filename="chromatogrampm.pdf")

hetseqallelesp <- setAllelePhase(hetcallsp1, ref, trim5=50, trim3=50)
pa <- pairwiseAlignment(primarySeq(hetseqallelesp), 
                        secondarySeq(hetseqallelesp), 
                        type="global-local")
writePairwiseAlignments(pa)


a <- readDNAStringSet("a.fasta")
b <- readDNAStringSet("b.fasta")

hetseqallelesp <- setAllelePhase(toString(b), toString(a), trim5=50, trim3=50)
a
b
pa <- pairwiseAlignment(toString(b), 
                        toString(a), 
                        type="global-local")
