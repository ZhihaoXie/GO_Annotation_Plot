#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author:    Zhihao Xie  \(#^o^)/
# Date:      2017.07.28
# Version:   v1.0.0
# CopyRight: Copyright ©Zhihao Xie, All rights reserved.

import sys
import re
import os
import argparse
from numpy import unique
Script_path = os.path.split(os.path.realpath(sys.argv[0]))[0]
import obo_parser as obop

def get_params():
    usage = "python3 %(prog)s -i <gene_go_file> -o <out_prefix> --obo <go-basic.obo>"
    parser = argparse.ArgumentParser(usage=usage, description="%(prog)s v1.0.0")
    parser.add_argument('-i', '--input', dest='go_file', help='a gene2go file as input')
    parser.add_argument('-o', '--out', dest='prefix', help='output prefix')
    parser.add_argument('-l', '--level', dest='level', default=2, type=int, help='go level, default is 2')
    parser.add_argument('--obo', dest='obo', default=Script_path+os.sep+"go-basic.obo", help='go-basic.obo, default is %s. you can download this file from http://geneontology.org/page/download-ontology' % (os.path.join(Script_path, "go-basic.obo")))

    args = parser.parse_args()
    if not args.go_file or not args.prefix:
        parser.print_help()
        sys.exit()
    if not os.path.exists(args.obo):
        sys.stderr.write("The %s doesn't exists. you can download this file from http://geneontology.org/page/download-ontology\n")
        sys.exit()
    return args

def main():
    args = get_params()
    go_file = os.path.abspath(args.go_file)
    prefix = args.prefix
    obo_file = os.path.abspath(args.obo)
    level = args.level

    # parse input file
    tmp_go_file = prefix+".temp_go_list.txt"
    with open(go_file, 'r') as fh:
        with open(tmp_go_file, 'w') as outfh:
            for line in fh:
                if not re.search(r"^\s*$", line):
                    gene_id, go_term = line.strip().split("\t", 1)
                    go_term = go_term.strip()
                    if re.search(r',', go_term):
                        go_lists = re.split(r",\s*", go_term)
                    elif re.search(r"\s", go_term):
                        go_lists = re.split(r'\s+', go_term)
                    else:
                        go_lists = [go_term]
                    for iterm in go_lists:
                        outfh.write("%s\t%s\n" % (gene_id, iterm))

    # get go parents
    g = obop.GODag(obo_file)
    
    go_level_out_file = prefix+".parseAllGO.AddLevel.txt"
    all_go_level_lists = []
    with open(tmp_go_file, 'r') as fh:
        for line in fh:
            if re.search(r'^\s*$', line):
                continue
            gene_id, go_term = line.strip().split("\t")
            rec = g.query_term(go_term)
            all_parents = rec.get_all_parents()
            all_go_lists = [x for x in all_parents]
            all_go_lists.append(go_term)

            for go_id in all_go_lists:
                temp_rec = g.query_term(go_id)
                go_name = temp_rec.name
                go_namespace = temp_rec.namespace
                go_level = int(temp_rec.level) + 1
                tmp_str = "%s\t%s\t%s\t%s\t%s\n" % (go_id, go_name, go_namespace, go_level, gene_id)
                all_go_level_lists.append(tmp_str)

    all_go_level_lists = unique(all_go_level_lists)
    with open(go_level_out_file, 'w') as fh:
        fh.write("#go_id\tname\tnamespace\tlevel\tgene\n")
        for i in all_go_level_lists:
            fh.write(i)
    
    # add level to gene
    tmp_hash_go = {}
    with open(go_level_out_file, 'r') as fh:
        for line in fh:
            if not re.search(r"^#|^\s*$", line):
                fields = line.strip().split("\t")
                tmp_hash_go.setdefault(fields[0], "\t".join(fields[1:4]))

    add_level2gene_fp = prefix+".gene2GO.AddLevel.txt"
    with open(add_level2gene_fp, 'w') as outfh:
        outfh.write("\t".join(["gene", "go_id", "name", "namespace", "level"])+"\n")
        with open(tmp_go_file, 'r') as fh:
            for line in fh:
                if not re.search("^\s*$", line):
                    fields = line.strip().split("\t")
                    if fields[1] in tmp_hash_go:
                        outfh.write(line.strip()+"\t"+tmp_hash_go[fields[1]]+"\n")
                    else:
                        outfh.write(line)

    # go statis
    go_statis_file = prefix+".parseAllGO.stat.xls"
    count = {}
    with open(go_level_out_file, 'r') as fh:
        for line in fh:
            if len(line)==0:
                break
            elif re.search(r'^\s*$|^#', line):
                continue
            arrays = line.strip().split("\t")
            consent = "\t".join(arrays[1:3])
            if int(arrays[3]) == level:
                count.setdefault(consent, 0)
                count[consent] += 1
    
    with open(go_statis_file, 'w') as fh:
        for key in count:
            if count[key] > 10:
                fh.write("%s\t%s\n" % (key, count[key]))

    # go statis plot Rscript
    plot_out_file = prefix+".GO_distribution_under_level"+str(level)+".pdf"
    with open("GOlevel2Grapher.R", 'w') as fh:
        cmd_str = """# Rscript GOlevel2Grapher.R project_GO.parseAllGO.stat.xls
argv <- commandArgs(TRUE)
inputfile <- argv[1]
outputfile <- "%s"
library(ggplot2)
z<-read.table(inputfile,sep="\\t",header=F)
names(z)<-c("Class","Ontology","Counts")
z<-z[order(z$Counts,decreasing=F),]
p <- ggplot(z,aes(Class,Counts)) + scale_fill_brewer(palette="Set1")
p <- p + facet_grid(.~Ontology,scales="free_x",space="free")+geom_bar(aes(fill=Ontology,weight=Counts,position="dodge"),stat="identity")+theme(legend.position="none",axis.text.x=element_text(angle=75,hjust=1.0,size=8,colour="black"))+scale_y_log10()+ylab("Number of Gene Hits")+xlab("GO term")
p <- p + theme(panel.background = element_blank(),panel.border = element_rect(colour="black",fill=NA,size=0.1),panel.grid.minor = element_line(colour = "grey", linetype = "dotted", size=0.1),panel.grid.major = element_line(colour = "grey", linetype = "dotted", size=0.1))
pdf(file=outputfile)
plot(p)
dev.off()
q()
""" % (plot_out_file)
        fh.write(cmd_str)

    cmd_str = """Rscript GOlevel2Grapher.R %s &&\\
convert -density 300 -resize 50%% %s %s
""" % (go_statis_file, plot_out_file, plot_out_file.replace(".pdf", ".png"))
    status = os.system(cmd_str)
    if status == 0:
        sys.stderr.write("good. go plot finished.\n")
        os.remove(tmp_go_file)
    else:
        sys.stderr.write("bad. go plot failed.\n")


# <-- main --> #
if __name__ == "__main__":
    main()
