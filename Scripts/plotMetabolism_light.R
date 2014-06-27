
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/metab.light.RData")

# =============
# = Plot Ward =
# =============
dev.new(width=8, height=8)
par(mfcol=c(3,2), mar=c(3.5,3.5,0.1,0.1), mgp=c(1.5, 0.5, 0), tcl=-0.25, ps=12, family="Times", cex=1)


w10.ind <- metab.light[,"lake"]=="Ward" & metab.light[,"year"]==2010
w10.gpp.sonde <- metab.light[w10.ind,"gpp.sonde"]
w10.gpp.bottom <- metab.light[w10.ind,"gpp.bot"]
w10.gpp.total <- w10.gpp.sonde + w10.gpp.bottom

w12.ind <- metab.light[,"lake"]=="Ward" & metab.light[,"year"]==2012
w12.gpp.sonde <- metab.light[w12.ind,"gpp.sonde"]
w12.gpp.bottom <- metab.light[w12.ind, "gpp.bot"]
w12.gpp.total <- w12.gpp.sonde + w12.gpp.bottom

#plot ward 2010
plot(metab.light[w10.ind, "doy"], w10.gpp.sonde, type="o", xlab="", ylab=bquote(Top~~GPP~~(g~O[2]~m^-2~d^-1)), ylim=range(c(w10.gpp.sonde, w12.gpp.sonde)))
plot(metab.light[w10.ind, "doy"], w10.gpp.bottom, type="o", xlab="", ylab=bquote(Bottom~~GPP~~(g~O[2]~m^-2~d^-1)), ylim=range(c(w10.gpp.bottom, w12.gpp.bottom)))
plot(metab.light[w10.ind, "doy"], w10.gpp.total, type="o", xlab="", ylab=bquote(Total~~GPP~~(g~O[2]~m^-2~d^-1)), ylim=range(c(w10.gpp.total, w12.gpp.total)))


#plot ward 2012
plot(metab.light[w12.ind, "doy"], w12.gpp.sonde, type="o", xlab="", ylab=bquote(Top~~GPP~~(g~O[2]~m^-2~d^-1)), ylim=range(c(w10.gpp.sonde, w12.gpp.sonde)))
plot(metab.light[w12.ind, "doy"], w12.gpp.bottom, type="o", xlab="", ylab=bquote(Bottom~~GPP~~(g~O[2]~m^-2~d^-1)), ylim=range(c(w10.gpp.bottom, w12.gpp.bottom)))
plot(metab.light[w12.ind, "doy"], w12.gpp.total, type="o", xlab="", ylab=bquote(Total~~GPP~~(g~O[2]~m^-2~d^-1)), ylim=range(c(w10.gpp.total, w12.gpp.total)))


# =============
# = Plot Paul =
# =============
dev.new(width=8, height=8)
par(mfcol=c(3,2), mar=c(3.5,3.5,0.1,0.1), mgp=c(1.5, 0.5, 0), tcl=-0.25, ps=12, family="Times", cex=1)
l10.ind <- metab.light[,"lake"]=="Paul" & metab.light[,"year"]==2010
l10.gpp.sonde <- metab.light[l10.ind,"gpp.sonde"]
l10.gpp.bottom <- metab.light[l10.ind,"gpp.bot"]
l10.gpp.total <- l10.gpp.sonde + l10.gpp.bottom

l12.ind <- metab.light[,"lake"]=="Paul" & metab.light[,"year"]==2012
l12.gpp.sonde <- metab.light[l12.ind,"gpp.sonde"]
l12.gpp.bottom <- metab.light[l12.ind, "gpp.bot"]
l12.gpp.total <- l12.gpp.sonde + l12.gpp.bottom

#plot lard 2010
plot(metab.light[l10.ind, "doy"], l10.gpp.sonde, type="o", xlab="", ylab=bquote(Top~~GPP~~(g~O[2]~m^-2~d^-1)), ylim=range(c(l10.gpp.sonde, l12.gpp.sonde)))
plot(metab.light[l10.ind, "doy"], l10.gpp.bottom, type="o", xlab="", ylab=bquote(Bottom~~GPP~~(g~O[2]~m^-2~d^-1)), ylim=range(c(l10.gpp.bottom, l12.gpp.bottom)))
plot(metab.light[l10.ind, "doy"], l10.gpp.total, type="o", xlab="", ylab=bquote(Total~~GPP~~(g~O[2]~m^-2~d^-1)), ylim=range(c(l10.gpp.total, l12.gpp.total)))


#plot lard 2012
plot(metab.light[l12.ind, "doy"], l12.gpp.sonde, type="o", xlab="", ylab=bquote(Top~~GPP~~(g~O[2]~m^-2~d^-1)), ylim=range(c(l10.gpp.sonde, l12.gpp.sonde)))
plot(metab.light[l12.ind, "doy"], l12.gpp.bottom, type="o", xlab="", ylab=bquote(Bottom~~GPP~~(g~O[2]~m^-2~d^-1)), ylim=range(c(l10.gpp.bottom, l12.gpp.bottom)))
plot(metab.light[l12.ind, "doy"], l12.gpp.total, type="o", xlab="", ylab=bquote(Total~~GPP~~(g~O[2]~m^-2~d^-1)), ylim=range(c(l10.gpp.total, l12.gpp.total)))

