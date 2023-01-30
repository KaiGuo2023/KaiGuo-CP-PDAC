


library(GEOquery)

eSet <- getGEO("GSE143754",
               destdir = '.',
               getGPL = F)  # ��ȡƽ̨��Ϣ
exp <- exprs(eSet[[1]])  # �������
GPL <- eSet[[1]]@annotation  # ƽ̨��Ϣ������ȡоƬƽ̨���
GPL

GPL=getGEO(filename = 'GSE143754_family.soft.gz')
gpl=GPL@gpls[[1]]@dataTable@table
colnames(gpl)

write.csv(gpl[,c(1,8)],file="gpl.csv")
###########



gpl=read.csv("gpl1.csv",header = T)

ids=gpl


#install.packages("dplyr")
library(dplyr)
colnames(ids) = c("probe_id" ,"symbol")
exp=as.data.frame(exp)
exp$probe_id=rownames(exp) # ��������Ϊ����Ϊprobe_id��һ��
# exp��ԭ���ı������
exp2= merge(exp,ids,by.x="probe_id", by.y="probe_id") # �ϲ�����
#����ظ���������Ĺ߳�����
library(limma)
exp2[,1]=exp2[,ncol(exp2)]
exp2=exp2[,-ncol(exp2)]

rt=as.matrix(exp2)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=t(avereps(data))
exp2=t(avereps(data))
#end

write.table(exp2,file = "ids_exprs.txt",sep = "\t",row.names=T,col.names = T)
write.csv(exp2,file = "ids_exprs.csv")
