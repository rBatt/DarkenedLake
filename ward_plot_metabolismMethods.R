

# ================================================
# = Plot Ward 2010 Epi Metabolism across Methods =
# ================================================
# dev.new(width=3.5, height=7)
pdf("~/Desktop/ward2010_epilimnion.pdf", width=3.5, height=7)
par(mfcol=c(5,2), mar=c(2, 2, 0.5, 0.5), oma=c(0.1, 0.1, 2, 0.1), ps=9, mgp=c(1, 0.25, 0), tcl=-0.35, family="Times", cex=1)

plot(w10e.r4[w10e.r4[,"method"]=="bookkeep","R"], w10e.r4[w10e.r4[,"method"]=="ols","R"], xlab="BK", ylab="OLS"); abline(a=0, b=1)
plot(w10e.r4[w10e.r4[,"method"]=="bookkeep","R"], w10e.r4[w10e.r4[,"method"]=="mle","R"], xlab="BK", ylab="MLE"); abline(a=0, b=1)
plot(w10e.r4[w10e.r4[,"method"]=="bookkeep","R"], w10e.r4[w10e.r4[,"method"]=="kalman","R"], xlab="BK", ylab="Kalman"); abline(a=0, b=1)
plot(w10e.r4[w10e.r4[,"method"]=="bookkeep","R"], w10e.r4[w10e.r4[,"method"]=="bayes","R"], xlab="BK", ylab="Bayes"); abline(a=0, b=1)

plot(w10e.r4[w10e.r4[,"method"]=="ols","R"], w10e.r4[w10e.r4[,"method"]=="mle","R"], xlab="OLS", ylab="MLE"); abline(a=0, b=1)
plot(w10e.r4[w10e.r4[,"method"]=="ols","R"], w10e.r4[w10e.r4[,"method"]=="kalman","R"], xlab="OLS", ylab="Kalman"); abline(a=0, b=1)
plot(w10e.r4[w10e.r4[,"method"]=="ols","R"], w10e.r4[w10e.r4[,"method"]=="bayes","R"], xlab="OLS", ylab="Bayes"); abline(a=0, b=1)

plot(w10e.r4[w10e.r4[,"method"]=="mle","R"], w10e.r4[w10e.r4[,"method"]=="kalman","R"], xlab="MLE", ylab="Kalman"); abline(a=0, b=1)
plot(w10e.r4[w10e.r4[,"method"]=="mle","R"], w10e.r4[w10e.r4[,"method"]=="bayes","R"], xlab="MLE", ylab="Bayes"); abline(a=0, b=1)

plot(w10e.r4[w10e.r4[,"method"]=="kalman","R"], w10e.r4[w10e.r4[,"method"]=="bayes","R"], xlab="Kalman", ylab="Bayes"); abline(a=0, b=1)

mtext("Ward 2010 Epi\nRespiration", outer=TRUE, line=0, side=3, font=2, cex=1.2)


# dev.new(width=3.5, height=7)
par(mfcol=c(5,2), mar=c(2, 2, 0.5, 0.5), oma=c(0.1, 0.1, 2, 0.1), ps=9, mgp=c(1, 0.25, 0), tcl=-0.35, family="Times", cex=1)

plot(w10e.r4[w10e.r4[,"method"]=="bookkeep","GPP"], w10e.r4[w10e.r4[,"method"]=="ols","GPP"], xlab="BK", ylab="OLS"); abline(a=0, b=1)
plot(w10e.r4[w10e.r4[,"method"]=="bookkeep","GPP"], w10e.r4[w10e.r4[,"method"]=="mle","GPP"], xlab="BK", ylab="MLE"); abline(a=0, b=1)
plot(w10e.r4[w10e.r4[,"method"]=="bookkeep","GPP"], w10e.r4[w10e.r4[,"method"]=="kalman","GPP"], xlab="BK", ylab="Kalman"); abline(a=0, b=1)
plot(w10e.r4[w10e.r4[,"method"]=="bookkeep","GPP"], w10e.r4[w10e.r4[,"method"]=="bayes","GPP"], xlab="BK", ylab="Bayes"); abline(a=0, b=1)

plot(w10e.r4[w10e.r4[,"method"]=="ols","GPP"], w10e.r4[w10e.r4[,"method"]=="mle","GPP"], xlab="OLS", ylab="MLE"); abline(a=0, b=1)
plot(w10e.r4[w10e.r4[,"method"]=="ols","GPP"], w10e.r4[w10e.r4[,"method"]=="kalman","GPP"], xlab="OLS", ylab="Kalman"); abline(a=0, b=1)
plot(w10e.r4[w10e.r4[,"method"]=="ols","GPP"], w10e.r4[w10e.r4[,"method"]=="bayes","GPP"], xlab="OLS", ylab="Bayes"); abline(a=0, b=1)

plot(w10e.r4[w10e.r4[,"method"]=="mle","GPP"], w10e.r4[w10e.r4[,"method"]=="kalman","GPP"], xlab="MLE", ylab="Kalman"); abline(a=0, b=1)
plot(w10e.r4[w10e.r4[,"method"]=="mle","GPP"], w10e.r4[w10e.r4[,"method"]=="bayes","GPP"], xlab="MLE", ylab="Bayes"); abline(a=0, b=1)

plot(w10e.r4[w10e.r4[,"method"]=="kalman","GPP"], w10e.r4[w10e.r4[,"method"]=="bayes","GPP"], xlab="Kalman", ylab="Bayes"); abline(a=0, b=1)

mtext("Ward 2010 Epi\nGross Primary Production", outer=TRUE, line=0, side=3, font=2, cex=1.2)


# dev.new(width=3.5, height=7)
par(mfcol=c(5,2), mar=c(2, 2, 0.5, 0.5), oma=c(0.1, 0.1, 2, 0.1), ps=9, mgp=c(1, 0.25, 0), tcl=-0.35, family="Times", cex=1)

plot(w10e.r4[w10e.r4[,"method"]=="bookkeep","NEP"], w10e.r4[w10e.r4[,"method"]=="ols","NEP"], xlab="BK", ylab="OLS"); abline(a=0, b=1)
plot(w10e.r4[w10e.r4[,"method"]=="bookkeep","NEP"], w10e.r4[w10e.r4[,"method"]=="mle","NEP"], xlab="BK", ylab="MLE"); abline(a=0, b=1)
plot(w10e.r4[w10e.r4[,"method"]=="bookkeep","NEP"], w10e.r4[w10e.r4[,"method"]=="kalman","NEP"], xlab="BK", ylab="Kalman"); abline(a=0, b=1)
plot(w10e.r4[w10e.r4[,"method"]=="bookkeep","NEP"], w10e.r4[w10e.r4[,"method"]=="bayes","NEP"], xlab="BK", ylab="Bayes"); abline(a=0, b=1)

plot(w10e.r4[w10e.r4[,"method"]=="ols","NEP"], w10e.r4[w10e.r4[,"method"]=="mle","NEP"], xlab="OLS", ylab="MLE"); abline(a=0, b=1)
plot(w10e.r4[w10e.r4[,"method"]=="ols","NEP"], w10e.r4[w10e.r4[,"method"]=="kalman","NEP"], xlab="OLS", ylab="Kalman"); abline(a=0, b=1)
plot(w10e.r4[w10e.r4[,"method"]=="ols","NEP"], w10e.r4[w10e.r4[,"method"]=="bayes","NEP"], xlab="OLS", ylab="Bayes"); abline(a=0, b=1)

plot(w10e.r4[w10e.r4[,"method"]=="mle","NEP"], w10e.r4[w10e.r4[,"method"]=="kalman","NEP"], xlab="MLE", ylab="Kalman"); abline(a=0, b=1)
plot(w10e.r4[w10e.r4[,"method"]=="mle","NEP"], w10e.r4[w10e.r4[,"method"]=="bayes","NEP"], xlab="MLE", ylab="Bayes"); abline(a=0, b=1)

plot(w10e.r4[w10e.r4[,"method"]=="kalman","NEP"], w10e.r4[w10e.r4[,"method"]=="bayes","NEP"], xlab="Kalman", ylab="Bayes"); abline(a=0, b=1)

mtext("Ward 2010 Epi\nNet Ecosystem Production", outer=TRUE, line=0, side=3, font=2, cex=1.2)

dev.off()





# ================================================
# = Plot Ward 2010 Epi Metabolism across Methods =
# ================================================
pdf("~/Desktop/ward2010_metalimnion.pdf", width=3.5, height=7)
par(mfcol=c(5,2), mar=c(2, 2, 0.5, 0.5), oma=c(0.1, 0.1, 2, 0.1), ps=9, mgp=c(1, 0.25, 0), tcl=-0.35, family="Times", cex=1)

plot(w10m.r4[w10m.r4[,"method"]=="bookkeep","R"], w10m.r4[w10m.r4[,"method"]=="ols","R"], xlab="BK", ylab="OLS"); abline(a=0, b=1)
plot(w10m.r4[w10m.r4[,"method"]=="bookkeep","R"], w10m.r4[w10m.r4[,"method"]=="mle","R"], xlab="BK", ylab="MLE"); abline(a=0, b=1)
plot(w10m.r4[w10m.r4[,"method"]=="bookkeep","R"], w10m.r4[w10m.r4[,"method"]=="kalman","R"], xlab="BK", ylab="Kalman"); abline(a=0, b=1)
plot(w10m.r4[w10m.r4[,"method"]=="bookkeep","R"], w10m.r4[w10m.r4[,"method"]=="bayes","R"], xlab="BK", ylab="Bayes"); abline(a=0, b=1)

plot(w10m.r4[w10m.r4[,"method"]=="ols","R"], w10m.r4[w10m.r4[,"method"]=="mle","R"], xlab="OLS", ylab="MLE"); abline(a=0, b=1)
plot(w10m.r4[w10m.r4[,"method"]=="ols","R"], w10m.r4[w10m.r4[,"method"]=="kalman","R"], xlab="OLS", ylab="Kalman"); abline(a=0, b=1)
plot(w10m.r4[w10m.r4[,"method"]=="ols","R"], w10m.r4[w10m.r4[,"method"]=="bayes","R"], xlab="OLS", ylab="Bayes"); abline(a=0, b=1)

plot(w10m.r4[w10m.r4[,"method"]=="mle","R"], w10m.r4[w10m.r4[,"method"]=="kalman","R"], xlab="MLE", ylab="Kalman"); abline(a=0, b=1)
plot(w10m.r4[w10m.r4[,"method"]=="mle","R"], w10m.r4[w10m.r4[,"method"]=="bayes","R"], xlab="MLE", ylab="Bayes"); abline(a=0, b=1)

plot(w10m.r4[w10m.r4[,"method"]=="kalman","R"], w10m.r4[w10m.r4[,"method"]=="bayes","R"], xlab="Kalman", ylab="Bayes"); abline(a=0, b=1)

mtext("Ward 2010 Meta\n Respiration", outer=TRUE, line=0, side=3, font=2, cex=1.2)



par(mfcol=c(5,2), mar=c(2, 2, 0.5, 0.5), oma=c(0.1, 0.1, 2, 0.1), ps=9, mgp=c(1, 0.25, 0), tcl=-0.35, family="Times", cex=1)

plot(w10m.r4[w10m.r4[,"method"]=="bookkeep","GPP"], w10m.r4[w10m.r4[,"method"]=="ols","GPP"], xlab="BK", ylab="OLS"); abline(a=0, b=1)
plot(w10m.r4[w10m.r4[,"method"]=="bookkeep","GPP"], w10m.r4[w10m.r4[,"method"]=="mle","GPP"], xlab="BK", ylab="MLE"); abline(a=0, b=1)
plot(w10m.r4[w10m.r4[,"method"]=="bookkeep","GPP"], w10m.r4[w10m.r4[,"method"]=="kalman","GPP"], xlab="BK", ylab="Kalman"); abline(a=0, b=1)
plot(w10m.r4[w10m.r4[,"method"]=="bookkeep","GPP"], w10m.r4[w10m.r4[,"method"]=="bayes","GPP"], xlab="BK", ylab="Bayes"); abline(a=0, b=1)

plot(w10m.r4[w10m.r4[,"method"]=="ols","GPP"], w10m.r4[w10m.r4[,"method"]=="mle","GPP"], xlab="OLS", ylab="MLE"); abline(a=0, b=1)
plot(w10m.r4[w10m.r4[,"method"]=="ols","GPP"], w10m.r4[w10m.r4[,"method"]=="kalman","GPP"], xlab="OLS", ylab="Kalman"); abline(a=0, b=1)
plot(w10m.r4[w10m.r4[,"method"]=="ols","GPP"], w10m.r4[w10m.r4[,"method"]=="bayes","GPP"], xlab="OLS", ylab="Bayes"); abline(a=0, b=1)

plot(w10m.r4[w10m.r4[,"method"]=="mle","GPP"], w10m.r4[w10m.r4[,"method"]=="kalman","GPP"], xlab="MLE", ylab="Kalman"); abline(a=0, b=1)
plot(w10m.r4[w10m.r4[,"method"]=="mle","GPP"], w10m.r4[w10m.r4[,"method"]=="bayes","GPP"], xlab="MLE", ylab="Bayes"); abline(a=0, b=1)

plot(w10m.r4[w10m.r4[,"method"]=="kalman","GPP"], w10m.r4[w10m.r4[,"method"]=="bayes","GPP"], xlab="Kalman", ylab="Bayes"); abline(a=0, b=1)

mtext("Ward 2010 Meta\nGross Primary Production", outer=TRUE, line=0, side=3, font=2, cex=1.2)


# dev.new(width=3.5, height=7)
par(mfcol=c(5,2), mar=c(2, 2, 0.5, 0.5), oma=c(0.1, 0.1, 2, 0.1), ps=9, mgp=c(1, 0.25, 0), tcl=-0.35, family="Times", cex=1)

plot(w10m.r4[w10m.r4[,"method"]=="bookkeep","NEP"], w10m.r4[w10m.r4[,"method"]=="ols","NEP"], xlab="BK", ylab="OLS"); abline(a=0, b=1)
plot(w10m.r4[w10m.r4[,"method"]=="bookkeep","NEP"], w10m.r4[w10m.r4[,"method"]=="mle","NEP"], xlab="BK", ylab="MLE"); abline(a=0, b=1)
plot(w10m.r4[w10m.r4[,"method"]=="bookkeep","NEP"], w10m.r4[w10m.r4[,"method"]=="kalman","NEP"], xlab="BK", ylab="Kalman"); abline(a=0, b=1)
plot(w10m.r4[w10m.r4[,"method"]=="bookkeep","NEP"], w10m.r4[w10m.r4[,"method"]=="bayes","NEP"], xlab="BK", ylab="Bayes"); abline(a=0, b=1)

plot(w10m.r4[w10m.r4[,"method"]=="ols","NEP"], w10m.r4[w10m.r4[,"method"]=="mle","NEP"], xlab="OLS", ylab="MLE"); abline(a=0, b=1)
plot(w10m.r4[w10m.r4[,"method"]=="ols","NEP"], w10m.r4[w10m.r4[,"method"]=="kalman","NEP"], xlab="OLS", ylab="Kalman"); abline(a=0, b=1)
plot(w10m.r4[w10m.r4[,"method"]=="ols","NEP"], w10m.r4[w10m.r4[,"method"]=="bayes","NEP"], xlab="OLS", ylab="Bayes"); abline(a=0, b=1)

plot(w10m.r4[w10m.r4[,"method"]=="mle","NEP"], w10m.r4[w10m.r4[,"method"]=="kalman","NEP"], xlab="MLE", ylab="Kalman"); abline(a=0, b=1)
plot(w10m.r4[w10m.r4[,"method"]=="mle","NEP"], w10m.r4[w10m.r4[,"method"]=="bayes","NEP"], xlab="MLE", ylab="Bayes"); abline(a=0, b=1)

plot(w10m.r4[w10m.r4[,"method"]=="kalman","NEP"], w10m.r4[w10m.r4[,"method"]=="bayes","NEP"], xlab="Kalman", ylab="Bayes"); abline(a=0, b=1)

mtext("Ward 2010 Meta\nNet Ecosystem Production", outer=TRUE, line=0, side=3, font=2, cex=1.2)

dev.off()






# ================================================
# = Plot Ward 2012 Epi Metabolism across Methods =
# ================================================
# dev.new(width=3.5, height=7)
pdf("~/Desktop/ward2012_epilimnion.pdf", width=3.5, height=7)
par(mfcol=c(5,2), mar=c(2, 2, 0.5, 0.5), oma=c(0.1, 0.1, 2, 0.1), ps=9, mgp=c(1, 0.25, 0), tcl=-0.35, family="Times", cex=1)

plot(w12e.r4[w12e.r4[,"method"]=="bookkeep","R"], w12e.r4[w12e.r4[,"method"]=="ols","R"], xlab="BK", ylab="OLS"); abline(a=0, b=1)
plot(w12e.r4[w12e.r4[,"method"]=="bookkeep","R"], w12e.r4[w12e.r4[,"method"]=="mle","R"], xlab="BK", ylab="MLE"); abline(a=0, b=1)
plot(w12e.r4[w12e.r4[,"method"]=="bookkeep","R"], w12e.r4[w12e.r4[,"method"]=="kalman","R"], xlab="BK", ylab="Kalman"); abline(a=0, b=1)
plot(w12e.r4[w12e.r4[,"method"]=="bookkeep","R"], w12e.r4[w12e.r4[,"method"]=="bayes","R"], xlab="BK", ylab="Bayes"); abline(a=0, b=1)

plot(w12e.r4[w12e.r4[,"method"]=="ols","R"], w12e.r4[w12e.r4[,"method"]=="mle","R"], xlab="OLS", ylab="MLE"); abline(a=0, b=1)
plot(w12e.r4[w12e.r4[,"method"]=="ols","R"], w12e.r4[w12e.r4[,"method"]=="kalman","R"], xlab="OLS", ylab="Kalman"); abline(a=0, b=1)
plot(w12e.r4[w12e.r4[,"method"]=="ols","R"], w12e.r4[w12e.r4[,"method"]=="bayes","R"], xlab="OLS", ylab="Bayes"); abline(a=0, b=1)

plot(w12e.r4[w12e.r4[,"method"]=="mle","R"], w12e.r4[w12e.r4[,"method"]=="kalman","R"], xlab="MLE", ylab="Kalman"); abline(a=0, b=1)
plot(w12e.r4[w12e.r4[,"method"]=="mle","R"], w12e.r4[w12e.r4[,"method"]=="bayes","R"], xlab="MLE", ylab="Bayes"); abline(a=0, b=1)

plot(w12e.r4[w12e.r4[,"method"]=="kalman","R"], w12e.r4[w12e.r4[,"method"]=="bayes","R"], xlab="Kalman", ylab="Bayes"); abline(a=0, b=1)

mtext("Ward 2012 Epi\nRespiration", outer=TRUE, line=0, side=3, font=2, cex=1.2)


# dev.new(width=3.5, height=7)
par(mfcol=c(5,2), mar=c(2, 2, 0.5, 0.5), oma=c(0.1, 0.1, 2, 0.1), ps=9, mgp=c(1, 0.25, 0), tcl=-0.35, family="Times", cex=1)

plot(w12e.r4[w12e.r4[,"method"]=="bookkeep","GPP"], w12e.r4[w12e.r4[,"method"]=="ols","GPP"], xlab="BK", ylab="OLS"); abline(a=0, b=1)
plot(w12e.r4[w12e.r4[,"method"]=="bookkeep","GPP"], w12e.r4[w12e.r4[,"method"]=="mle","GPP"], xlab="BK", ylab="MLE"); abline(a=0, b=1)
plot(w12e.r4[w12e.r4[,"method"]=="bookkeep","GPP"], w12e.r4[w12e.r4[,"method"]=="kalman","GPP"], xlab="BK", ylab="Kalman"); abline(a=0, b=1)
plot(w12e.r4[w12e.r4[,"method"]=="bookkeep","GPP"], w12e.r4[w12e.r4[,"method"]=="bayes","GPP"], xlab="BK", ylab="Bayes"); abline(a=0, b=1)

plot(w12e.r4[w12e.r4[,"method"]=="ols","GPP"], w12e.r4[w12e.r4[,"method"]=="mle","GPP"], xlab="OLS", ylab="MLE"); abline(a=0, b=1)
plot(w12e.r4[w12e.r4[,"method"]=="ols","GPP"], w12e.r4[w12e.r4[,"method"]=="kalman","GPP"], xlab="OLS", ylab="Kalman"); abline(a=0, b=1)
plot(w12e.r4[w12e.r4[,"method"]=="ols","GPP"], w12e.r4[w12e.r4[,"method"]=="bayes","GPP"], xlab="OLS", ylab="Bayes"); abline(a=0, b=1)

plot(w12e.r4[w12e.r4[,"method"]=="mle","GPP"], w12e.r4[w12e.r4[,"method"]=="kalman","GPP"], xlab="MLE", ylab="Kalman"); abline(a=0, b=1)
plot(w12e.r4[w12e.r4[,"method"]=="mle","GPP"], w12e.r4[w12e.r4[,"method"]=="bayes","GPP"], xlab="MLE", ylab="Bayes"); abline(a=0, b=1)

plot(w12e.r4[w12e.r4[,"method"]=="kalman","GPP"], w12e.r4[w12e.r4[,"method"]=="bayes","GPP"], xlab="Kalman", ylab="Bayes"); abline(a=0, b=1)

mtext("Ward 2012 Epi\nGross Primary Production", outer=TRUE, line=0, side=3, font=2, cex=1.2)


# dev.new(width=3.5, height=7)
par(mfcol=c(5,2), mar=c(2, 2, 0.5, 0.5), oma=c(0.1, 0.1, 2, 0.1), ps=9, mgp=c(1, 0.25, 0), tcl=-0.35, family="Times", cex=1)

plot(w12e.r4[w12e.r4[,"method"]=="bookkeep","NEP"], w12e.r4[w12e.r4[,"method"]=="ols","NEP"], xlab="BK", ylab="OLS"); abline(a=0, b=1)
plot(w12e.r4[w12e.r4[,"method"]=="bookkeep","NEP"], w12e.r4[w12e.r4[,"method"]=="mle","NEP"], xlab="BK", ylab="MLE"); abline(a=0, b=1)
plot(w12e.r4[w12e.r4[,"method"]=="bookkeep","NEP"], w12e.r4[w12e.r4[,"method"]=="kalman","NEP"], xlab="BK", ylab="Kalman"); abline(a=0, b=1)
plot(w12e.r4[w12e.r4[,"method"]=="bookkeep","NEP"], w12e.r4[w12e.r4[,"method"]=="bayes","NEP"], xlab="BK", ylab="Bayes"); abline(a=0, b=1)

plot(w12e.r4[w12e.r4[,"method"]=="ols","NEP"], w12e.r4[w12e.r4[,"method"]=="mle","NEP"], xlab="OLS", ylab="MLE"); abline(a=0, b=1)
plot(w12e.r4[w12e.r4[,"method"]=="ols","NEP"], w12e.r4[w12e.r4[,"method"]=="kalman","NEP"], xlab="OLS", ylab="Kalman"); abline(a=0, b=1)
plot(w12e.r4[w12e.r4[,"method"]=="ols","NEP"], w12e.r4[w12e.r4[,"method"]=="bayes","NEP"], xlab="OLS", ylab="Bayes"); abline(a=0, b=1)

plot(w12e.r4[w12e.r4[,"method"]=="mle","NEP"], w12e.r4[w12e.r4[,"method"]=="kalman","NEP"], xlab="MLE", ylab="Kalman"); abline(a=0, b=1)
plot(w12e.r4[w12e.r4[,"method"]=="mle","NEP"], w12e.r4[w12e.r4[,"method"]=="bayes","NEP"], xlab="MLE", ylab="Bayes"); abline(a=0, b=1)

plot(w12e.r4[w12e.r4[,"method"]=="kalman","NEP"], w12e.r4[w12e.r4[,"method"]=="bayes","NEP"], xlab="Kalman", ylab="Bayes"); abline(a=0, b=1)

mtext("Ward 2012 Epi\nNet Ecosystem Production", outer=TRUE, line=0, side=3, font=2, cex=1.2)

dev.off()

