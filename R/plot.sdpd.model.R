plot.sdpd.model <- function (res.fit, n.units="all", n.vars="all", which=c(1,2)[1], t.axis=list(t.labels=NULL, t.points=NULL), xlimit=NULL, ylimit=NULL, max.col=5, col.punti=1)
{
  if(which==1)
    plot1.sdpd.model(res.fit=res.fit, n.units=n.units, n.vars=n.vars, t.axis=t.axis, xlimit=xlimit, ylimit=ylimit, max.col=max.col, col.punti=col.punti)
}
