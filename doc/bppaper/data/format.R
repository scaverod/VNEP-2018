library("xtable")

for (file in c('fdensevir1', 'fdensevir4', 'fsparsevir1', 'fsparsevir4', 'fhiervir1', 'fhiervir4', 'ftsvir1', 'ftsvir4')) {
  print(file)
  a <- read.table(gettextf('%s.dat', file), sep='&')
  ag <- aggregate(a, by=list(a$V2), FUN=mean, na.rm=T)
  print(xtable(ag))
}
