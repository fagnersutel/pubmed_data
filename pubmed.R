#install.packages("pubmed.mineR", dependencies = T)
#install.packages('wordcloud', dependencies = T)
#install.packages('tm', dependencies = T)
#install.packages('fpc', dependencies = T)
#install.packages('cluster', dependencies = T)
#install.packages("rentrez", dependencies = T)
library(rentrez)
library (wordcloud)
library (tm)
library (RISmed)
library (cluster)
library(pubmed.mineR)
library(ggplot2)

pubtator_output <- pubtator_function(29656905)
head(pubtator_output)


#query <- "(exome OR whole OR deep OR high-throughput OR (next AND generation) OR (massively AND parallel)) AND sequencing"
query <- "bioethics AND medical AND Brazil"

search <- EUtilsSummary(query, type="esearch",db = "pubmed",mindate=2000, maxdate=2018, retmax=30000)
QueryCount(search)
records <- EUtilsGet(search)
years <- RISmed::YearPubmed(records)
pubs_count <- as.data.frame(table(years))
pubs_count

total <- NULL
for (i in 2000:2018){
  peryear <- EUtilsSummary("", type="esearch", db="pubmed", mindate=i, maxdate=i)
  total[i] <- QueryCount(peryear)
}

year <- 2000:2018
total_pubs_count<- as.data.frame(cbind(year,total[year]))
names(total_pubs_count) <- c("year","Total_publications")
names(pubs_count) <-  c("year","publications")
pubs_year <-  merge(pubs_count,total_pubs_count,by="year")
pubs_year$publications_normalized <-  pubs_year$publications *100000 / pubs_year$Total_publications
write.table(pubs_year,"publications_per_year.txt",quote=F,sep="\t",row.names=F)

journal <- MedlineTA(records)
journal_count <- as.data.frame(table(journal))
journal_count_top25 <- journal_count[order(-journal_count[,2]),][1:25,]

journal_names <- paste(journal_count_top25$journal,"[jo]",sep="")

total_journal <- NULL
for (i in journal_names){
  perjournal <- EUtilsSummary(i, type='esearch', db='pubmed',mindate=1980, maxdate=2013)
  total_journal[i] <- QueryCount(perjournal)
}

journal_total <- cbind(journal_count_top25,total_journal)
names(journal_total) <- c("journal","publications","Total_publications")
journal_total$publications_normalized <- journal_total$publications / journal_total$Total_publications

write.table(journal_total,"publications_per_journal.txt",quote=F,sep="\t",row.names=F)


pubs_per_year <- read.table("publications_per_year.txt",header = T,sep="\t")
pubs_per_journal <- read.table("publications_per_journal.txt",header = T,sep="\t")

ggplot(pubs_per_year,aes(year, publications_normalized)) + geom_line (colour="blue",size=2) +
  xlab("Ano") +
  ylab("Artigos/100000 artigos")+
  ggtitle("Artigos {Bioetica, Medical, Brazil}PubMed")

ggplot(pubs_per_journal,aes(journal, publications,fill=journal)) + geom_bar(stat="identity")+
  coord_flip()+
  theme(legend.position="none")

ggplot(pubs_per_journal ,aes(journal, publications_normalized,fill=journal)) + geom_bar(stat="identity")+
  coord_flip()+
  theme(legend.position="none")




#Procurar dados de autor
res <- EUtilsSummary('Wolnei Caumo', type='esearch', db='pubmed')

summary(res)
QueryId(res)


res2 <- EUtilsSummary('Wolnei Caumo', type='esearch', db='pubmed', mindate='2012', maxdate='2018')

QueryId(res2)




res3 <- EUtilsSummary('chronobiology', type='esearch', db='pubmed')

summary(res3)


QueryCount(res3)

#tally each year beginning at 1970
#In order not to overload the E-utility servers, NCBI recommends that users post no more than three
#URL requests per second and limit large jobs to either weekends or between 9:00 PM and 5:00 AM
#Eastern time during weekdays. Failure to comply with this policy may result in an IP address being
#blocked from accessing NCBI.

tally <- array()
x <- 1
for (i in 2013:2017){
  Sys.sleep(1)
  r <- EUtilsSummary('Pain Neuromodulation', type='esearch', db='pubmed', mindate=i, maxdate=i)
  tally[x] <- QueryCount(r)
  x <- x + 1
}

names(tally) <- 2013:2017
max(tally)

barplot(tally, las=2, ylim=c(0,50000), main="Num artigos PUBMED contendo Pain Neuromodulatio")



transposon <- array()
x <- 1
for (i in 2013:2017){
  Sys.sleep(1)
  r <- EUtilsSummary('pain and therapeutic', type='esearch', db='pubmed', mindate=i, maxdate=i)
  transposon[x] <- QueryCount(r)
  x <- x + 1
}

names(transposon) <- 2013:2017
max(transposon)


barplot(transposon, las=2, ylim=c(0,20000), main="Num artigos PUBMED contendo pain and therapeutic")


trna <- array()
x <- 1
for (i in 2013:2017){
  Sys.sleep(1)
  r <- EUtilsSummary('Pharmacological modulation', type='esearch', db='pubmed', mindate=i, maxdate=i)
  trna[x] <- QueryCount(r)
  x <- x + 1
}

names(trna) <- 2013:2017
max(trna)

barplot(trna, las=2, ylim=c(0,5000), main="Num artigos PUBMED contendo Pharmacological modulation")



test <- EUtilsSummary('', type='esearch', db='pubmed', mindate=2013, maxdate=2017)
summary(test)

#artigos por ano
total <- array()
x <- 1
for (i in 2013:2017){
  Sys.sleep(1)
  r <- EUtilsSummary('', type='esearch', db='pubmed', mindate=i, maxdate=i)
  total[x] <- QueryCount(r)
  x <- x + 1
}

names(total) <- 2013:2017
max(total)

options(scipen=999)
barplot(total, las=2, ylim=c(0,1200000), main="Numero de artigos do PubMed por ano")


#var um pelo total
tally_norm <- tally / total
#var dois pelo total
transposon_norm <- transposon / total
#var tres pelo total
trna_norm <- trna / total

par(mfrow=c(2,2))
barplot(tally_norm, las=2)
title(main = "Pain Neuromodulation", font.main = 4)
barplot(transposon_norm, las=2)
title(main = "Pain and Therapeutic", font.main = 4)
barplot(trna_norm, las=2)
title(main = "Pharmacological modulation", font.main = 4)
barplot(total, las=2)
title(main = "Universal", font.main = 4)
#reset
par(mfrow=c(1,1))






query <- 'Capp E'
query_level2 <- EUtilsSummary(query, retmax=100, mindate=2016, maxdate=2016)
query_level3<- EUtilsGet(query_level2)
class(query_level3)
pubmed_data <- data.frame('Abstract'= AbstractText(query_level3))
for (Abs in 1:100){
    doc1 <- data.frame(pubmed_data[Abs, ])
    doc2 <- file.path("c:/Users/fsmoura/Documents/R-files/abs/", paste0(Abs, ".txt"))
    write.table(doc1, file = doc2, sep = "", row.names = FALSE, col.names = FALSE, quote = FALSE,
                 append = FALSE)
}

source <- DirSource("c:/Users/fsmoura/Documents/R-files/abs/")
testdoc <- Corpus(source)
testdoc1 <- tm_map(testdoc, removeWords, c("may","are","use","can","the", "then", "this", "is", "a", "well", stopwords ("english")))
testdoc1
testdoc2 <- TermDocumentMatrix (testdoc1, control = list(tokenize=scan_tokenizer, stopwords =
                                                            TRUE, removePunctuation = TRUE, stripWhitespace = TRUE, stemming = TRUE, removeNumbers=
                                                            TRUE
))
testdoc2
testdoc3 <- as.matrix(testdoc2)
View(testdoc3)
testdoc4 <- sort(rowSums(testdoc3),decreasing=TRUE)
head(testdoc4, 15)
testdoc5 <- data.frame(word = names(testdoc4),freq=testdoc4)
head(testdoc5, 15)

findAssocs(x=testdoc2, term="ovarian", corlimit=0.6)

set.seed(1234)
wordcloud(words = testdoc5$word, freq = testdoc5$freq, min.freq = 1,
           max.words=200, random.order=FALSE, rot.per=0.35,
           colors=brewer.pal(8, "Dark2"))


testdoc5 <- removeSparseTerms(testdoc2, 0.70)
c1 <- as.matrix(testdoc5)
c1
c2 <- dist(c1)
c2
c3 <- hclust(c2, method="ward.D")
plot(c1, hang=-1) 


plot.new()
plot(c3, hang=-1) 


km1 <- kmeans(c2, 2)
clusplot(as.matrix(c2), km1$cluster, color=T, shade=T, labels=2, lines=0)


doctest <- DocumentTermMatrix (testdoc1, control = list (tokenize=scan_tokenizer, stopwords =
                                                           TRUE, removePunctuation = TRUE,
                                                         stripWhitespace = TRUE,
                                                         stemming = TRUE,
                                                         removeNumbers= TRUE
                                                         
))
c1 <- as.matrix(testdoc5)
c2 <- dist(c1)
c3 <- hclust(c2, method="ward.D") 

plot(c3, hang=-1)
plot.new()
plot(c3, hang=-1) 




#[AFFL]    Affiliation
#[ALL]     All Fields
#[AUTH]    Author
#[FAUT]    Author - First
#[LAUT]    Author - Last
#[PDAT]    Date - Publication
#[FILT]    Filter
#[JOUR]    Journal
#[LANG]    Language
#[MAJR]    MeSH Major Topic
#[SUBH]    MeSH Subheading
#[MESH]    MeSH Terms
#[PTYP]    Publication Type
#[WORD]    Text Word
#[TITL]    Title
#[TIAB]    Title/Abstract
#[UID]     UID


pubmed_search <- entrez_search(db="pubmed", term="(Transcranial direct current[TITL]) AND 2017[PDATT]",use_history=TRUE)
pkg_paper_summs <- entrez_summary(db="pubmed", web_history=pubmed_search$web_history)
pkg_paper_summs

journals <- extract_from_esummary(pkg_paper_summs, "fulljournalname")
journals_by_R_pkgs <- sort(table(journals), decreasing = TRUE)
head(journals_by_R_pkgs,10)
View(head(journals_by_R_pkgs,10))


pubmed_search <- entrez_search(db="pubmed", term="(Transcranial direct current[TITL]) AND 2017[PDATT] AND CAUMO W[AUTH]",use_history=TRUE)
pkg_paper_summs <- entrez_summary(db="pubmed", web_history=pubmed_search$web_history)
pkg_paper_summs

journals <- extract_from_esummary(pkg_paper_summs, "fulljournalname")
journals_by_R_pkgs <- sort(table(journals), decreasing = TRUE)
head(journals_by_R_pkgs,10)
View(head(journals_by_R_pkgs,10))
pkg_paper_summs$`29274881`$title
pkg_paper_summs$`29274881`$fulljournalname
pkg_paper_summs$`29274881`$authors
pkg_paper_summs$`29274881`$references
pkg_paper_summs$`29274881`$sortfirstauthor



app_gene <- entrez_search(db="gene", term="(Homo sapiens[ORGN]) AND APP[GENE]")
app_gene
nuc_links <- entrez_link(dbfrom="gene", id=app_gene$ids, db="nuccore")
nuc_links$links
raw_recs <- entrez_fetch(db="nuccore",
                         id=nuc_links$links$gene_nuccore_refseqrna,
                         rettype="fasta")
cat(substr(raw_recs, 1,303), "...")
cat(raw_recs, file="APP_transcripts.fasta")
tf <- tempfile()
cat(raw_recs, file=tf)
ape::read.dna(tf, format="fasta")



#install.packages("easyPubMed", dependencies = TRUE )
library(easyPubMed)
first_PM_query <- "T cell receptor"                  # text of the query
first_PM_records <- get_pubmed_ids(first_PM_query)   # submit the query to PubMed
fetch_pubmed_data(first_PM_records, retmax = 1)      # retrieve the first output record


library(easyPubMed)
my_query <- "Wolnei Caumo[AU]"
my_entrez_id <- get_pubmed_ids(my_query)
my_abstracts_txt <- fetch_pubmed_data(my_entrez_id, format = "abstract")
my_abstracts_txt[1:10]


my_abstracts_xml <- fetch_pubmed_data(my_entrez_id)
class(my_abstracts_xml)
my_titles <- unlist(xpathApply(my_abstracts_xml, "//ArticleTitle", saveXML))
#Para não pesar pedi os primeiros 10 apenas
my_titles <- gsub("(^.{5,10}Title>)|(<\\/.*$)", "", my_titles[1:10])
my_titles[nchar(my_titles)>75] <- paste(substr(my_titles[nchar(my_titles)>75], 1, 70),
                                        "...", sep = "")
print(my_titles)



new_query <- "(Transcranial direc[TI] OR Neurobiological[TI]) AND (2012[PDAT]:2016[PDAT])"
out.A <- batch_pubmed_download(pubmed_query_string = new_query,
                               format = "xml",
                               batch_size = 150,
                               dest_file_prefix = "easyPM_example")
out.A
my_PM_list <- articles_to_list(my_abstracts_xml)
class(my_PM_list[[4]])
cat(substr(my_PM_list[[4]], 1, 984))


curr_PM_record <- my_PM_list[[4]]
custom_grep(curr_PM_record, tag = "ArticleTitle", format = "char")
custom_grep(curr_PM_record, tag = "EIdType")
custom_grep(curr_PM_record, tag = "Title")
custom_grep(curr_PM_record, tag = "DateRevised")
custom_grep(curr_PM_record, tag = "Year", format = "char")
custom_grep(curr_PM_record, tag = "LastName", format = "char")


my.df <- article_to_df(curr_PM_record, max_chars = 18)
colnames(my.df)
my.df$title <- substr(my.df$title, 1, 15)
my.df$address <- substr(my.df$address, 1, 19)
my.df$jabbrv <- substr(my.df$jabbrv, 1, 10)
my.df[,c("pmid", "title", "jabbrv", "firstname", "address")]



my.df2 <- article_to_df(curr_PM_record, autofill = TRUE)
my.df2$title <- substr(my.df2$title, 1, 15)
my.df2$jabbrv <- substr(my.df2$jabbrv, 1, 10)

my.df2$address <- substr(my.df2$address, 1, 19)
my.df2[,c("pmid", "title", "jabbrv", "firstname", "address")]



new_PM_query <- "(Transcranial direc[TI] OR Neurobiological[TI]) AND (2012[PDAT]:2016[PDAT])"
out.B <- batch_pubmed_download(pubmed_query_string = new_PM_query, dest_file_prefix = "apex1_sample")
new_PM_file <- out.B[1]
new_PM_df <- table_articles_byAuth(pubmed_data = new_PM_file, included_authors = "first", max_chars = 0)
new_PM_df$address <- substr(new_PM_df$address, 1, 28)
new_PM_df$jabbrv <- substr(new_PM_df$jabbrv, 1, 9)
print(new_PM_df[1:10, c("pmid", "year", "jabbrv", "lastname", "address")]) 




# Example 01: retrieve data in XML format
dami_query_string <- "Damiano Fantini[AU]"
dami_on_pubmed <- get_pubmed_ids(dami_query_string)
dami_papers <- fetch_pubmed_data(dami_on_pubmed)
titles <- unlist(xpathApply(dami_papers, "//ArticleTitle", saveXML))
title_pos <- regexpr('<ArticleTitle>.*<\\/ArticleTitle>', titles)
titles <- substr(titles, title_pos + 14, title_pos + attributes(title_pos)$match.length - 16)
print(titles)
#
# Example 02: retrieve data in TXT format
dami_query_string <- "Damiano Fantini[AU]"
dami_on_pubmed <- get_pubmed_ids(dami_query_string)
dami_papers <- fetch_pubmed_data(dami_on_pubmed, format = "abstract")
dami_papers[dami_papers == ""] <- "\n"
cat(paste(dami_papers[1:10], collapse = ""))
#
## Not run: 
# Example 03: retrieve data from PubMed and save as XML file
ml_query <- "Machine Learning[TI] AND 2016[PD]"
out1 <- batch_pubmed_download(pubmed_query_string = ml_query, batch_size = 10)
XML::xmlParse(out1[1])
#
# Example 04: retrieve data from PubMed and save as TXT file
ml_query <- "Machine Learning[TI] AND 2016[PD]"
out2 <- batch_pubmed_download(pubmed_query_string = ml_query, batch_size = 180, format = "medline")
readLines(out2[1])[1:30]
#
# Example 05: extract information from a single PubMed record 
ml_query <- "Machine Learning[TI] AND 2016[PD]"
out3 <- batch_pubmed_download(pubmed_query_string = ml_query, batch_size = 180)
PM_data <- articles_to_list(out3[1])
PM_record_df <- article_to_df(PM_data[[100]])
print(PM_record_df[1,])
print(PM_record_df[,"address"])
#
# Example 06: query PubMed and extract information from multiple records in one step 
ml_query <- "Machine Learning[TI] AND 2016[PD]"
out4 <- batch_pubmed_download(pubmed_query_string = ml_query, batch_size = 10)
PM_tab <- table_articles_byAuth(out4[1], autofill = TRUE, included_authors = "last")
PM_tab$address <- substr(PM_tab$address, 1, 15)
PM_tab[50:70,c("pmid", "jabbrv", "year", "lastname", "address")]
#



library(purrr)
Query <- c('rituximab OR bevacizumab','meningitis OR headache')
Heading <- c('A','B')
x <- as.data.frame(cbind(Heading,Query),stringsAsFactors = F)
x$PMID<- ""
ids <- map(x[,"Query"],get_pubmed_ids)
for (i in 1:length(ids)) {
  x[i,"PMID"]<- paste(ids[[i]][["IdList"]],collapse = ",")
}
x







search_topic <- 'copd'
search_query <- EUtilsSummary(search_topic, retmax=100, mindate=2012,maxdate=2012)
summary(search_query)

# see the ids of our returned query
QueryId(search_query)

# get actual data from PubMed
records<- EUtilsGet(search_query)
class(records)

# store it
pubmed_data <- data.frame('Title'=ArticleTitle(records),'Abstract'=AbstractText(records))
head(pubmed_data,1)

pubmed_data$Abstract <- as.character(pubmed_data$Abstract)
pubmed_data$Abstract <- gsub(",", " ", pubmed_data$Abstract, fixed = TRUE)

# see what we have
str(pubmed_data)

table(pubmed_data)

