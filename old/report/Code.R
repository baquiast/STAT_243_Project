par(mfrow=c(2,2))
y = ars(n=5,dnorm,mean=0,sd=1,x.start=c(-2,2), plot.type="bounds",
        pdf.name="Normal Distribution(n=5): ")  
y = ars(n=10,dnorm,mean=0,sd=1,x.start=c(-2,2), plot.type="bounds",
        pdf.name="Normal Distribution(n=10): ")  
y = ars(n=30,dnorm,mean=0,sd=1,x.start=c(-2,2), plot.type="bounds",
        pdf.name="Normal Distribution(n=30): ")  
y = ars(n=100,dnorm,mean=0,sd=1,x.start=c(-2,2), plot.type="bounds",
        pdf.name="Normal Distribution(n=100): ")  
