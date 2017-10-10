
logit <- data.frame(response=c("What's that?", "I've heard of it", "I can use it or apply it", "I understand it well"), Frequency=c(0,5,13,0))
probit <- data.frame(response=c("What's that?", "I've heard of it", "I can use it or apply it", "I understand it well"), Frequency=c(6,6,6,0))
ml <- data.frame(response=c("What's that?", "I've heard of it", "I can use it or apply it", "I understand it well"), Frequency=c(2,11,4,1))
R <- data.frame(response=c("What's that?", "I've heard of it", "I can use it or apply it", "I understand it well"), Frequency=c(2,4,12,0))
latex <- data.frame(response=c("What's that?", "I've heard of it", "I can use it or apply it", "I understand it well"), Frequency=c(8,8,3,0))

logit$response <- factor(logit$response, levels=c("What's that?", "I've heard of it", "I can use it or apply it", "I understand it well"))
probit$response <- factor(probit$response, levels=c("What's that?", "I've heard of it", "I can use it or apply it", "I understand it well"))
ml$response <- factor(ml$response, levels=c("What's that?", "I've heard of it", "I can use it or apply it", "I understand it well"))
R$response <- factor(R$response, levels=c("What's that?", "I've heard of it", "I can use it or apply it", "I understand it well"))
latex$response <- factor(latex$response, levels=c("What's that?", "I've heard of it", "I can use it or apply it", "I understand it well"))


logitplot <- ggplot(data=logit, aes(x=response, y=Frequency))+geom_bar(stat="identity", fill="steelblue", width=0.5)+
  geom_text(aes(label=Frequency, vjust=-0.3, size=0.5))+theme(plot.title = element_text(hjust = 0.5), legend.position="none")+scale_y_continuous(limits=c(0, 15))+labs(title="Logit", x="Response", y="Frequency")

logitplot

probitplot <- ggplot(data=probit, aes(x=response, y=Frequency))+geom_bar(stat="identity", fill="steelblue", width=0.5)+
  geom_text(aes(label=Frequency, vjust=-0.3, size=0.5))+theme(plot.title = element_text(hjust = 0.5), legend.position="none")+scale_y_continuous(limits=c(0, 15))+labs(title="Probit", x="Response", y="Frequency")

probitplot

mlplot <- ggplot(data=ml, aes(x=response, y=Frequency))+geom_bar(stat="identity", fill="steelblue", width=0.5)+
  geom_text(aes(label=Frequency, vjust=-0.3, size=0.5))+theme(plot.title = element_text(hjust = 0.5), legend.position="none")+scale_y_continuous(limits=c(0, 15))+labs(title="Maximum Likelihood", x="Response", y="Frequency")

mlplot

Rplot <- ggplot(data=ml, aes(x=response, y=Frequency))+geom_bar(stat="identity", fill="steelblue", width=0.5)+
  geom_text(aes(label=Frequency, vjust=-0.3, size=0.5))+theme(plot.title = element_text(hjust = 0.5), legend.position="none")+scale_y_continuous(limits=c(0, 15))+labs(title="R", x="Response", y="Frequency")

Rplot


latexplot <- ggplot(data=latex, aes(x=response, y=Frequency))+geom_bar(stat="identity", fill="steelblue", width=0.5)+
  geom_text(aes(label=Frequency, vjust=-0.3, size=0.5))+theme(plot.title = element_text(hjust = 0.5), legend.position="none")+scale_y_continuous(limits=c(0, 15))+labs(title="Latex", x="Response", y="Frequency")

latexplot

