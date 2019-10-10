using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using Mono.Options;
using Bio;
using Bio.IO;
using System.Drawing.Drawing2D;
using System.Drawing;

namespace SELDLA{
    class Chain{
        public struct chr_match{
            public string chr;
            public List<double> match;
        }
        public struct chain{
            public string chr1;
            public string chr2;
            public string pos1;
            public string pos2;
            public double match;
        }
        public struct stranded_chr{
            public string chr;
            public int order;
            public double match;
        }
        public struct edge{
            public string chr;
            public string pos;
            public string prechr;
            public double match;
        }
        public struct scafpos{
            public string chr;
            public string order;
            public long pos1;
            public long pos2;
            public double cm1;
            public double cm2;
        }
        public struct linkagescaf{
            public long length;
            public List<scafpos> scafs;
        }
        int nbp=10*1000;
        public void run(Dictionary<string, string> refseqs2, Dictionary<string, SortedDictionary<int, Dictionary<string, int[]>>> datas, double opt_s, double opt_nonzerophase, int num_member, string opt_o){
            List<string> croschr=new List<string>();
            Dictionary<string, double> matchcroschr = new Dictionary<string, double>();
            foreach(KeyValuePair<string, SortedDictionary<int, Dictionary<string, int[]>>> pair in datas){
                double tempmatch = match_rate(pair.Value, pair.Value, "start", "end", opt_nonzerophase, num_member);
                if(tempmatch!=1 && tempmatch>0){
                    croschr.Add(pair.Key);
                    matchcroschr.Add(pair.Key, tempmatch);
                }
            }
            List<string> extrachr=new List<string>();
            foreach(string tempchr in datas.Keys){
                if(!croschr.Contains(tempchr)){
                    extrachr.Add(tempchr);
                }
            }

            //scaffoldの両端のフェーズが異なるscaffold間の一致率計算
            List<SortableMatch> sortedmatch = new List<SortableMatch>();
            List<chr_match> mainmatchss = croschr
                                                .AsParallel()
                                                .Select(f => calc_match_rate(f, "start", "start", croschr, datas, opt_nonzerophase, num_member)).ToList();
            foreach(chr_match tempmatch in mainmatchss){
                int i=0;
                foreach(string tempchr in croschr){
                    if (tempmatch.chr != tempchr)
                    {
                        SortableMatch temp = new SortableMatch(refseqs2[tempmatch.chr].Length + refseqs2[tempchr].Length
                         , tempmatch.match[i], tempmatch.chr, tempchr, "start", "start");
                        sortedmatch.Add(temp);
                        // if((temp.chr1=="scaffold_1" && temp.chr2=="scaffold_3")||(temp.chr1=="scaffold_3" && temp.chr2=="scaffold_1")){
                        //     Console.WriteLine(temp.matchrate+", "+temp.chr1+", "+temp.chr2);
                        // }
                    }
                    i++;
                }
            }
            List<chr_match> mainmatchse = croschr
                                                .AsParallel()
                                                .Select(f => calc_match_rate(f, "start", "end", croschr, datas, opt_nonzerophase, num_member)).ToList();
            foreach(chr_match tempmatch in mainmatchse){
                int i=0;
                foreach(string tempchr in croschr){
                    if (tempmatch.chr != tempchr)
                    {
                        SortableMatch temp = new SortableMatch(refseqs2[tempmatch.chr].Length + refseqs2[tempchr].Length
                         , tempmatch.match[i], tempmatch.chr, tempchr, "start", "end");
                        sortedmatch.Add(temp);
                    }
                    i++;
                }
            }
            List<chr_match> mainmatchee = croschr
                                                .AsParallel()
                                                .Select(f => calc_match_rate(f, "end", "end", croschr, datas, opt_nonzerophase, num_member)).ToList();
            foreach(chr_match tempmatch in mainmatchee){
                int i=0;
                foreach(string tempchr in croschr){
                    if (tempmatch.chr != tempchr)
                    {
                        SortableMatch temp = new SortableMatch(refseqs2[tempmatch.chr].Length + refseqs2[tempchr].Length
                         , tempmatch.match[i], tempmatch.chr, tempchr, "end", "end");
                        sortedmatch.Add(temp);
                    }
                    i++;
                }
            }

            sortedmatch.Sort();
            sortedmatch.Reverse();
            List<chain> mainchains = new List<chain>();
            List<string> flagmainstart = new List<string>();
            List<string> flagmainend = new List<string>();
            Dictionary<string, List<string>> flagmain = new Dictionary<string, List<string>>();
            flagmain.Add("start", flagmainstart);
            flagmain.Add("end", flagmainend);
            for(int i=0; i<sortedmatch.Count; i++){
                // if((sortedmatch[i].chr1=="scaffold_5" && sortedmatch[i].chr2=="scaffold_70") || (sortedmatch[i].chr2=="scaffold_5" && sortedmatch[i].chr1=="scaffold_70")){
                //     int temp=10;
                // }
                if(sortedmatch[i].matchrate<opt_s){
                    break;
                }
                bool isdo = true;
                if(flagmain[sortedmatch[i].pos1].Contains(sortedmatch[i].chr1)){
                    isdo=false;
                }else if(flagmain[sortedmatch[i].pos2].Contains(sortedmatch[i].chr2)){
                    isdo=false;
                }
                if(isdo){
                    chain[] tempall = mainchains.Where(f=> (f.chr1 == sortedmatch[i].chr1 && f.chr2 == sortedmatch[i].chr2)
                                          || (f.chr1 == sortedmatch[i].chr2 && f.chr2 == sortedmatch[i].chr1)).ToArray();
                    if(tempall.Length>0){
                        isdo=false;
                    }
                }
                if(isdo){
                    flagmain[sortedmatch[i].pos1].Add(sortedmatch[i].chr1);
                    flagmain[sortedmatch[i].pos2].Add(sortedmatch[i].chr2);
                    chain tempchain1 = new chain();
                    tempchain1.chr1=sortedmatch[i].chr1;
                    tempchain1.chr2=sortedmatch[i].chr2;
                    tempchain1.pos1=sortedmatch[i].pos1;
                    tempchain1.pos2=sortedmatch[i].pos2;
                    tempchain1.match=sortedmatch[i].matchrate;
                    mainchains.Add(tempchain1);
                    chain tempchain2 = new chain();
                    tempchain2.chr2=sortedmatch[i].chr1;
                    tempchain2.chr1=sortedmatch[i].chr2;
                    tempchain2.pos2=sortedmatch[i].pos1;
                    tempchain2.pos1=sortedmatch[i].pos2;
                    tempchain2.match=sortedmatch[i].matchrate;
                    mainchains.Add(tempchain2);

                    // if((sortedmatch[i].chr1=="scaffold_5" && sortedmatch[i].chr2=="scaffold_70") || (sortedmatch[i].chr2=="scaffold_5" && sortedmatch[i].chr1=="scaffold_70")){
                    //     Console.WriteLine(sortedmatch[i].matchrate+" "+sortedmatch[i].chr1+" "+sortedmatch[i].chr2+" "+sortedmatch[i].length+" "+sortedmatch[i].pos1+" "+sortedmatch[i].pos2);
                    // }
                }
            }

            //scaffoldの両端のフェーズが一致するextra scaffoldとメインのscaffold間の一致率計算
            // List<SortableMatch> sortedextramatch = new List<SortableMatch>();
            // ParallelQuery<chr_match> extramatchss = extrachr
            //                                     .AsParallel()
            //                                     .Select(f => calc_match_rate(f, "start", "start", croschr, datas));
            // foreach(chr_match tempmatch in extramatchss){
            //     int i=0;
            //     foreach(string tempchr in croschr){
            //         if (tempmatch.chr != tempchr)
            //         {
            //             SortableMatch temp = new SortableMatch(refseqs2[tempmatch.chr].Length + refseqs2[tempchr].Length
            //              , tempmatch.match[i], tempmatch.chr, tempchr, "start", "start");
            //             sortedextramatch.Add(temp);
            //         }
            //         i++;
            //     }
            // }
            List<edge> extramatchse = extrachr
                                            .AsParallel()
                                            .Select(f => get_max_edge(f, croschr, datas, opt_nonzerophase, num_member)).ToList();

            //scaffoldの並び順を探索する
            List<string> flagpassed = new List<string>();
            int num_ls=0;
            List<linkagescaf> finallinks = new List<linkagescaf>();
            foreach(string chr in croschr){
                if(!flagpassed.Contains(chr)){
                    flagpassed.Add(chr);
                    List<stranded_chr> stchrs = new List<stranded_chr>();
                    List<stranded_chr> tempstchrs = searchRev(chr, "start", mainchains, flagpassed);
                    double oldmatch=1;
                    foreach(stranded_chr temp in tempstchrs){
                        stranded_chr newtemp = new stranded_chr();
                        newtemp.chr=temp.chr;
                        newtemp.order=temp.order;
                        newtemp.match=oldmatch;
                        stchrs.Add(newtemp);
                        oldmatch=temp.match;
                    }
                    stranded_chr tempstchr = new stranded_chr();
                    tempstchr.chr=chr;
                    tempstchr.order=1;
                    tempstchr.match=oldmatch;
                    stchrs.Add(tempstchr);
                    stchrs.AddRange(searchFor(chr, "end", mainchains, flagpassed));

                    //最終的な連鎖地図の作成
                    //Console.WriteLine("###"+chr);
                    num_ls++;
                    //Console.WriteLine("linkage_scaffold_"+num_ls);
                    int num_ordered_sca = 0;
                    double tempdistcm = 0;
                    long tempdist = 0;
                    string temp_prev_chr="";
                    int temp_prev_order=0;
                    double temp_prev_distcm=0;
                    linkagescaf finallink = new linkagescaf();
                    finallink.scafs = new List<scafpos>();
                    foreach(stranded_chr stchr in stchrs){
                        num_ordered_sca++;
                        if(num_ordered_sca==1){
                            double tempmatch=-1;
                            List<edge> tempedges = new List<edge>();
                            if (stchr.order == 1)
                            {
                                tempedges = extramatchse.Where(x => x.chr == stchr.chr && x.match >= opt_s && x.pos == "start").OrderBy(x => x.match).ToList();
                            }
                            else
                            {
                                tempedges = extramatchse.Where(x => x.chr == stchr.chr && x.match >= opt_s && x.pos == "end").OrderBy(x => x.match).ToList();
                            }
                            foreach (edge tempedge in tempedges)
                            {
                                if (tempmatch == -1)
                                {
                                    tempdistcm = 0;
                                    tempdist=0;
                                }
                                else
                                {
                                    tempdistcm+= (tempedge.match - tempmatch) * 100;
                                    tempdist+=nbp;
                                }
                                finallink.length=(tempdist+refseqs2[tempedge.prechr].Length);
                                finallink.scafs.Add(createScafPos(tempedge.prechr, "na", tempdist, (tempdist+refseqs2[tempedge.prechr].Length), tempdistcm, tempdistcm));
                                //Console.WriteLine(num_ls + "\t" + tempedge.prechr + "\tna\t" + tempdist+"\t"+ (tempdist+refseqs2[tempedge.prechr].Length) + "\t" + tempdistcm + "\t" + tempdistcm
                                // + "\t" + tempedge.chr + " " + tempedge.pos + " " + tempedge.match);
                                tempmatch = tempedge.match;
                                tempdist+=refseqs2[tempedge.prechr].Length;
                            }
                            if(tempmatch==-1){
                                tempdistcm=0;
                                tempdist=0;
                            }else{
                                tempdistcm+=(1-tempmatch)*100;
                                tempdist+=nbp;
                            }
                            string temporder;
                            if(stchr.order==1){
                                temporder="+";
                            }else{
                                temporder="-";
                            }
                            finallink.length=(tempdist+refseqs2[stchr.chr].Length);
                            finallink.scafs.Add(createScafPos(stchr.chr, temporder, tempdist, (tempdist+refseqs2[stchr.chr].Length), tempdistcm, (tempdistcm+(1-matchcroschr[stchr.chr])*100)));
                            //Console.WriteLine(num_ls + "\t"+ stchr.chr+"\t"+temporder+"\t"+tempdist + "\t"+ (tempdist+refseqs2[stchr.chr].Length)+"\t"+tempdistcm+"\t"+(tempdistcm+(1-matchcroschr[stchr.chr])*100)
                            //+"\tmatch_rate_of_previous_ordered_scaffold:"+stchr.match+" match_rate_between_both_edge_of_this_scaffold"+matchcroschr[stchr.chr]);
                            tempdist+=refseqs2[stchr.chr].Length;
                            tempdistcm+=(1-matchcroschr[stchr.chr])*100;
                            temp_prev_distcm=tempdistcm;
                        }else if (num_ordered_sca > 1)
                        {
                            List<edge> tempedgespre = new List<edge>();
                            List<edge> tempedges = new List<edge>();
                            if(temp_prev_order==1){
                                tempedgespre = extramatchse.Where(x => x.chr == temp_prev_chr && x.match >= opt_s && x.pos == "end").OrderByDescending(x => x.match).ToList();
                            }else{
                                tempedgespre = extramatchse.Where(x => x.chr == temp_prev_chr && x.match >= opt_s && x.pos == "start").OrderByDescending(x => x.match).ToList();
                            }
                            if(stchr.order==1){
                                tempedges=extramatchse.Where(x => x.chr == stchr.chr && x.match >= opt_s && x.pos == "start").OrderBy(x => x.match).ToList();
                            }else{
                                tempedges=extramatchse.Where(x => x.chr == stchr.chr && x.match >= opt_s && x.pos == "end").OrderBy(x => x.match).ToList();
                            }
                            double totalmatch=0;
                            double matchnorm=1;
                            if(tempedgespre.Count>0){
                                totalmatch+=1-tempedgespre[tempedgespre.Count-1].match;
                            }
                            if(tempedges.Count>0){
                                totalmatch+=1-tempedges[0].match;
                            }
                            if(totalmatch>0 && totalmatch>(1-stchr.match)){
                                matchnorm=(1-stchr.match)/totalmatch;
                            }
                            foreach (edge tempedge in tempedgespre)
                            {
                                tempdistcm=temp_prev_distcm+(1-tempedge.match)*100*matchnorm;
                                tempdist+=nbp;
                                finallink.length=(tempdist+refseqs2[tempedge.prechr].Length);
                                finallink.scafs.Add(createScafPos(tempedge.prechr, "na", tempdist, (tempdist+refseqs2[tempedge.prechr].Length), tempdistcm, tempdistcm));
                                //Console.WriteLine(num_ls + "\t" + tempedge.prechr + "\tna\t" + tempdist+"\t"+ (tempdist+refseqs2[tempedge.prechr].Length) + "\t" + tempdistcm + "\t" + tempdistcm
                                // + "\t" + tempedge.chr + " " + tempedge.pos + " " + tempedge.match);
                                tempdist+=refseqs2[tempedge.prechr].Length;
                            }
                            double newtempdistcm = temp_prev_distcm+(1-stchr.match)*100;
                            foreach (edge tempedge in tempedges)
                            {
                                tempdistcm=newtempdistcm-(1-tempedge.match)*100*matchnorm;
                                tempdist+=nbp;
                                finallink.length=(tempdist+refseqs2[tempedge.prechr].Length);
                                finallink.scafs.Add(createScafPos(tempedge.prechr, "na", tempdist, (tempdist+refseqs2[tempedge.prechr].Length), tempdistcm, tempdistcm));
                                //Console.WriteLine(num_ls + "\t" + tempedge.prechr + "\tna\t" + tempdist+"\t"+ (tempdist+refseqs2[tempedge.prechr].Length) + "\t" + tempdistcm + "\t" + tempdistcm
                                // + "\t" + tempedge.chr + " " + tempedge.pos + " " + tempedge.match);
                                tempdist+=refseqs2[tempedge.prechr].Length;
                            }
                            tempdistcm=temp_prev_distcm+(1-stchr.match)*100;
                            tempdist+=nbp;
                            string temporder;
                            if(stchr.order==1){
                                temporder="+";
                            }else{
                                temporder="-";
                            }
                            finallink.length=(tempdist+refseqs2[stchr.chr].Length);
                            finallink.scafs.Add(createScafPos(stchr.chr, temporder, tempdist, (tempdist+refseqs2[stchr.chr].Length), tempdistcm, (tempdistcm+(1-matchcroschr[stchr.chr])*100)));
                            //Console.WriteLine(num_ls + "\t"+ stchr.chr+"\t"+temporder+"\t"+tempdist + "\t"+ (tempdist+refseqs2[stchr.chr].Length)+"\t"+tempdistcm+"\t"+(tempdistcm+(1-matchcroschr[stchr.chr])*100)
                            //+"\tmatch_rate_of_previous_ordered_scaffold:"+stchr.match+" match_rate_between_both_edge_of_this_scaffold"+matchcroschr[stchr.chr]);
                            tempdist+=refseqs2[stchr.chr].Length;
                            tempdistcm+=(1-matchcroschr[stchr.chr])*100;
                            temp_prev_distcm=tempdistcm;
                        }
                        temp_prev_chr=stchr.chr;
                        temp_prev_order=stchr.order;

                        if(num_ordered_sca == stchrs.Count){
                            double tempmatch=-1;
                            List<edge> tempedges = new List<edge>();
                            if (stchr.order == 1)
                            {
                                tempedges = extramatchse.Where(x => x.chr == stchr.chr && x.match >= opt_s && x.pos == "end").OrderByDescending(x => x.match).ToList();
                            }
                            else
                            {
                                tempedges = extramatchse.Where(x => x.chr == stchr.chr && x.match >= opt_s && x.pos == "start").OrderByDescending(x => x.match).ToList();
                            }
                            foreach (edge tempedge in tempedges)
                            {
                                if (tempmatch == -1)
                                {
                                    tempdistcm+=(1-tempedge.match)*100;
                                    tempdist+=nbp;
                                }
                                else
                                {
                                    tempdistcm+= (tempmatch - tempedge.match) * 100;
                                    tempdist+=nbp;
                                }
                                finallink.length=(tempdist+refseqs2[tempedge.prechr].Length);
                                finallink.scafs.Add(createScafPos(tempedge.prechr, "na", tempdist, (tempdist+refseqs2[tempedge.prechr].Length), tempdistcm, tempdistcm));
                                //Console.WriteLine(num_ls + "\t" + tempedge.prechr + "\tna\t" + tempdist+"\t"+ (tempdist+refseqs2[tempedge.prechr].Length) + "\t" + tempdistcm + "\t" + tempdistcm
                                // + "\t" + tempedge.chr + " " + tempedge.pos + " " + tempedge.match);
                                tempmatch = tempedge.match;
                                tempdist+=refseqs2[tempedge.prechr].Length;
                            }
                        }
                    }
                    finallinks.Add(finallink);
                }
            }
            if(finallinks.Count==0){
                throw new Exception("There were no crossed chromosomes.");
            }
            List<linkagescaf> listlinkages = finallinks.OrderByDescending(x=> x.length).ToList();
            double maxcm = listlinkages.Max(x=> x.scafs[x.scafs.Count-1].cm2);
            long maxbp = listlinkages.Max(x=> x.length);
            int num_big_ls = listlinkages.Where(x => x.length>1000*1000).Count();
            num_ls=0;
            StreamWriter writer = new StreamWriter(opt_o+"_chain.txt");
            StreamWriter writerfa = new StreamWriter(opt_o+"_extended.fasta");
            StreamWriter writerunloc = new StreamWriter(opt_o+"_unordered.fasta");
            StreamWriter writerfanainchr = new StreamWriter(opt_o+"_include_unordered_in_chr.fasta");
            int res_num_scaf = 0;
            int res_num_scaf_loc =0;
            int res_num_scaf_order =0;
            int res_num_scaf_order_big = 0;
            long res_bp_scaf=0;
            long res_bp_scaf_loc=0;
            long res_bp_scaf_order=0;
            long res_bp_scaf_order_big=0;
            List<string> inChain = new List<string>();
            List<string> inChainWithOriented = new List<string>(); 
            foreach(linkagescaf scafs in listlinkages){
                num_ls++;
                Console.WriteLine("#"+num_ls+" "+scafs.length);
                writer.WriteLine("#"+num_ls+" "+scafs.length);
                writerfa.WriteLine(">linkage_scaffold_"+num_ls);
                writerfanainchr.WriteLine(">linkage_scaffold_"+num_ls);
                int num_scaf=0;
                // var image = new Bitmap(500, 2000);
                // var g = Graphics.FromImage(image);
                // if (scafs.length > 1000 * 1000)
                // {
                //     g.SmoothingMode = SmoothingMode.AntiAlias;
                //     g.InterpolationMode = InterpolationMode.HighQualityBicubic;
                //     g.PixelOffsetMode = PixelOffsetMode.HighQuality;
                //     drawChrBase(g, maxbp, maxcm, 2, num_ls, scafs.length, scafs.scafs[scafs.scafs.Count-1].cm2, 0, 0);
                // }
                foreach(scafpos scaf in scafs.scafs){
                    num_scaf++;
                    res_num_scaf++;
                    res_bp_scaf+=refseqs2[scaf.chr].Length;
                    res_num_scaf_loc++;
                    res_bp_scaf_loc+=refseqs2[scaf.chr].Length;
                    //Console.WriteLine(num_ls+"\t"+ scaf.chr+"\t"+scaf.order+"\t"+scaf.pos1+"\t"+scaf.pos2+"\t"+scaf.cm1+"\t"+scaf.cm2);
                    writer.WriteLine(num_ls+"\t"+ scaf.chr+"\t"+scaf.order+"\t"+scaf.pos1+"\t"+scaf.pos2+"\t"+scaf.cm1+"\t"+scaf.cm2);
                    // if (scafs.length > 1000 * 1000)
                    // {
                    //     drawBpCm(g, 2, scaf, maxbp, maxcm, 0, 0);
                    // }
                    if(num_scaf>1){
                        for(int i=0;i<nbp;i++){
                            writerfa.Write("N");
                            writerfanainchr.Write("N");
                        }
                    }
                    if(scaf.order=="+"){
                        writerfa.Write(refseqs2[scaf.chr]);
                        writerfanainchr.Write(refseqs2[scaf.chr]);
                    }else if(scaf.order=="-"){
                        writerfa.Write(getRevComp(refseqs2[scaf.chr]));
                        writerfanainchr.Write(getRevComp(refseqs2[scaf.chr]));
                    }else{
                        for(int i=0;i<refseqs2[scaf.chr].Length;i++){
                            writerfa.Write("N");
                        }
                        writerfanainchr.Write(refseqs2[scaf.chr]);
                    }
                    if(scaf.order=="na"){
                        writerunloc.WriteLine(">linkage_scaffold_"+num_ls+"_unordered_pos"+scaf.pos1+"_old_"+scaf.chr);
                        writerunloc.WriteLine(refseqs2[scaf.chr]);
                        inChain.Add(scaf.chr);
                    }else{
                        res_num_scaf_order++;
                        res_bp_scaf_order+=refseqs2[scaf.chr].Length;
                        if(scafs.length>=1000*1000){
                            res_num_scaf_order_big++;
                            res_bp_scaf_order_big+=refseqs2[scaf.chr].Length;
                        }
                        inChain.Add(scaf.chr);
                        inChainWithOriented.Add(scaf.chr);
                    }
                }
                writerfa.WriteLine("");
                writerfanainchr.WriteLine("");
                // if (scafs.length > 1000 * 1000)
                // {
                //     image.Save(opt_o + "_map_"+num_ls+".png");
                // }
            }
            foreach(KeyValuePair<string, string> Pair in refseqs2){
                if(!inChain.Contains(Pair.Key)){
                    res_num_scaf++;
                    res_bp_scaf+=refseqs2[Pair.Key].Length;
                    writerfa.WriteLine(">old_"+Pair.Key);
                    writerfa.WriteLine(Pair.Value);
                    writerfanainchr.WriteLine(">old_"+Pair.Key);
                    writerfanainchr.WriteLine(Pair.Value);
                }
            }
            writer.Close();
            writerfa.Close();
            writerunloc.Close();
            writerfanainchr.Close();
            StreamWriter txtwrite = new StreamWriter(opt_o+"_log.txt");
            string str="Total: "+res_num_scaf.ToString("N0")+" scaffolds, Located: "+res_num_scaf_loc.ToString("N0")
            +" ("+((double)res_num_scaf_loc/res_num_scaf*100).ToString("F3")+"), Ordered: "
            +res_num_scaf_order.ToString("N0")+" ("+((double)res_num_scaf_order/res_num_scaf*100).ToString("F3")+"), Ordered in >=1Mbp chr ("+num_big_ls+"): "
            +res_num_scaf_order_big.ToString("N0")+" ("+((double)res_num_scaf_order_big/res_num_scaf*100).ToString("F3")+")";
            Console.WriteLine(str);
            txtwrite.WriteLine(str);

            str="Total: "+res_bp_scaf.ToString("N0")+" bp, Located: "+res_bp_scaf_loc.ToString("N0")
            +" ("+((double)res_bp_scaf_loc/res_bp_scaf*100).ToString("F3")+"), Ordered: "
            +res_bp_scaf_order+" ("+((double)res_bp_scaf_order/res_bp_scaf*100).ToString("F3")+"), Ordered in >=1Mbp chr ("+num_big_ls+"): "
            +res_bp_scaf_order_big+" ("+((double)res_bp_scaf_order_big/res_bp_scaf*100).ToString("F3")+")";
            Console.WriteLine(str);
            txtwrite.WriteLine(str);
            txtwrite.Close();

            int per_width=250;
            int per_height=1000;
            int all_per_num=8;
            try{
                var imageall = new Bitmap(per_width * all_per_num, per_height * ((num_big_ls - 1) / all_per_num + 1));
                var gall = Graphics.FromImage(imageall);
                gall.SmoothingMode = SmoothingMode.AntiAlias;
                gall.InterpolationMode = InterpolationMode.HighQualityBicubic;
                gall.PixelOffsetMode = PixelOffsetMode.HighQuality;
                for (int i = 1; i <= num_big_ls; i++)
                {
                    int j = (i - 1) % all_per_num + 1;
                    int k = (i - 1) / all_per_num + 1;
                    drawChrBase(gall, maxbp, maxcm, 1, i, listlinkages[i - 1].length, listlinkages[i - 1].scafs[listlinkages[i - 1].scafs.Count - 1].cm2, (j - 1) * per_width, (k - 1) * per_height);
                    foreach (scafpos scaf in listlinkages[i - 1].scafs)
                    {
                        drawBpCm(gall, 1, scaf, maxbp, maxcm, (j - 1) * per_width, (k - 1) * per_height);
                    }
                    drawChrFin(gall, maxbp, maxcm, 1, i, listlinkages[i - 1].length, listlinkages[i - 1].scafs[listlinkages[i - 1].scafs.Count - 1].cm2, (j - 1) * per_width, (k - 1) * per_height);
                }
                imageall.Save(opt_o + "_map_all.png");
            }
            catch (Exception e)
            {
                Console.WriteLine(e);
            }
        }

        public void drawChrFin(Graphics g, long maxbp, double maxcm, int fold, int num_ls, long each_maxbp, double each_maxcm, int startx, int starty){
            float y0=50;
            float y1=50+900*each_maxbp/maxbp;
            g.FillPolygon(Brushes.White, new[] {
                        new PointF(startx+25*fold,starty+(float)(y0-5)*fold),
                        new PointF(startx+25*fold,starty+(float)(y0+7.5-0)*fold),
                        new PointF(startx+30*fold,starty+(float)(y0+7.5-3)*fold),
                        new PointF(startx+35*fold,starty+(float)(y0+7.5-5.5)*fold),
                        new PointF(startx+40*fold,starty+(float)(y0+7.5-6.5)*fold),
                        new PointF(startx+45*fold,starty+(float)(y0+7.5-7)*fold),
                        new PointF(startx+50*fold,starty+(float)(y0+7.5-7.5)*fold),
                        new PointF(startx+55*fold,starty+(float)(y0+7.5-7)*fold),
                        new PointF(startx+60*fold,starty+(float)(y0+7.5-6.5)*fold),
                        new PointF(startx+65*fold,starty+(float)(y0+7.5-5.5)*fold),
                        new PointF(startx+70*fold,starty+(float)(y0+7.5-3)*fold),
                        new PointF(startx+75*fold,starty+(float)(y0+7.5-0)*fold),
                        new PointF(startx+75*fold,starty+(float)(y0-5)*fold),
            });
            g.FillPolygon(Brushes.White, new[] {
                        new PointF(startx+75*fold,starty+(float)(y1+5)*fold),
                        new PointF(startx+75*fold,starty+(float)(y1-7.5+0)*fold),
                        new PointF(startx+70*fold,starty+(float)(y1-7.5+3)*fold),
                        new PointF(startx+65*fold,starty+(float)(y1-7.5+5.5)*fold),
                        new PointF(startx+60*fold,starty+(float)(y1-7.5+6.5)*fold),
                        new PointF(startx+55*fold,starty+(float)(y1-7.5+7)*fold),
                        new PointF(startx+50*fold,starty+(float)(y1-7.5+7.5)*fold),
                        new PointF(startx+45*fold,starty+(float)(y1-7.5+7)*fold),
                        new PointF(startx+40*fold,starty+(float)(y1-7.5+6.5)*fold),
                        new PointF(startx+35*fold,starty+(float)(y1-7.5+5.5)*fold),
                        new PointF(startx+30*fold,starty+(float)(y1-7.5+3)*fold),
                        new PointF(startx+25*fold,starty+(float)(y1-7.5+0)*fold),
                        new PointF(startx+25*fold,starty+(float)(y1+5)*fold),
                    });
            var pen = new Pen(Color.Black, (float)1.5*fold);
            g.DrawLines(pen, new[] {
                        new PointF(startx+25*fold,starty+(float)(y0+7.5-0)*fold),
                        new PointF(startx+30*fold,starty+(float)(y0+7.5-3)*fold),
                        new PointF(startx+35*fold,starty+(float)(y0+7.5-5.5)*fold),
                        new PointF(startx+40*fold,starty+(float)(y0+7.5-6.5)*fold),
                        new PointF(startx+45*fold,starty+(float)(y0+7.5-7)*fold),
                        new PointF(startx+50*fold,starty+(float)(y0+7.5-7.5)*fold),
                        new PointF(startx+55*fold,starty+(float)(y0+7.5-7)*fold),
                        new PointF(startx+60*fold,starty+(float)(y0+7.5-6.5)*fold),
                        new PointF(startx+65*fold,starty+(float)(y0+7.5-5.5)*fold),
                        new PointF(startx+70*fold,starty+(float)(y0+7.5-3)*fold),
                        new PointF(startx+75*fold,starty+(float)(y0+7.5-0)*fold),

                        new PointF(startx+75*fold,starty+(float)(y1-7.5+0)*fold),
                        new PointF(startx+70*fold,starty+(float)(y1-7.5+3)*fold),
                        new PointF(startx+65*fold,starty+(float)(y1-7.5+5.5)*fold),
                        new PointF(startx+60*fold,starty+(float)(y1-7.5+6.5)*fold),
                        new PointF(startx+55*fold,starty+(float)(y1-7.5+7)*fold),
                        new PointF(startx+50*fold,starty+(float)(y1-7.5+7.5)*fold),
                        new PointF(startx+45*fold,starty+(float)(y1-7.5+7)*fold),
                        new PointF(startx+40*fold,starty+(float)(y1-7.5+6.5)*fold),
                        new PointF(startx+35*fold,starty+(float)(y1-7.5+5.5)*fold),
                        new PointF(startx+30*fold,starty+(float)(y1-7.5+3)*fold),
                        new PointF(startx+25*fold,starty+(float)(y1-7.5+0)*fold),

                        new PointF(startx+25*fold,starty+(float)(y0+7.5-0)*fold)
                    });
            Font myfont = new Font("Tahoma", 15*fold);
            g.DrawString("Linkage " + num_ls, myfont, Brushes.Black, new Point(startx+50*fold, starty+5*fold));
            myfont = new Font("Tahoma", 15*fold);
            g.DrawString("0 bp", myfont, Brushes.Black, new Point(startx+5*fold, starty+25*fold));
            g.DrawString("0 cM", myfont, Brushes.Black, new Point(startx+155*fold, starty+25*fold));
            g.DrawString(each_maxbp.ToString("N0") + " bp", myfont, Brushes.Black, new PointF(startx+5*fold, starty+(60 + 900 * (float)each_maxbp / maxbp)*fold));
            g.DrawString(each_maxcm.ToString("F1") + " cM", myfont, Brushes.Black, new PointF(startx+155*fold, starty+(60 + 900 * (float)(each_maxcm / maxcm))*fold));
        }
        public void drawBpCm(Graphics g, int fold, scafpos scaf, long maxbp, double maxcm, int startx, int starty){
            var pen = new Pen(Color.Gray, fold);
            var brush = Brushes.LightGray;
            g.DrawLine(pen, new PointF(startx+75*fold, starty+(50 + 900 * ((float)scaf.pos1 / maxbp))*fold), new PointF(startx+175*fold, starty+(50 + 900 * (float)(scaf.cm1 / maxcm))*fold));
            g.DrawLine(pen, new PointF(startx+75*fold, starty+(50 + 900 * ((float)scaf.pos2 / maxbp))*fold), new PointF(startx+175*fold, starty+(50 + 900 * (float)(scaf.cm2 / maxcm))*fold));
            if (scaf.order == "na")
            {
                pen = new Pen(Color.Blue, (float)1.5*fold);
                brush = Brushes.White;
            }
            else
            {
                pen = new Pen(Color.Red, (float)1.5*fold);
                brush = Brushes.LightGray;
            }
            g.FillRectangle(brush, new RectangleF(startx+26*fold, starty+(50 + 900 * ((float)scaf.pos1 / maxbp))*fold, 48*fold,900*(scaf.pos2-scaf.pos1)/(float)maxbp*fold));
            //g.DrawLine(pen, new PointF(startx+100*fold, starty+(50 + 900 * ((float)scaf.pos1 / maxbp))*fold), new PointF(startx+100*fold, starty+(50 + 900 *((float)scaf.pos2 / maxbp))*fold));
            g.DrawLine(pen, new PointF(startx+175*fold, starty+(50 + 900 * (float)(scaf.cm1 / maxcm))*fold), new PointF(startx+175*fold, starty+(50 + 900 * (float)(scaf.cm2 / maxcm))*fold));
            pen = new Pen(Color.Black, fold*2);
            g.DrawLine(pen, new PointF(startx+165*fold, starty+(50 + 900 * (float)(scaf.cm1 / maxcm))*fold), new PointF(startx+185*fold, starty+(50 + 900 * (float)(scaf.cm1 / maxcm))*fold));
            g.DrawLine(pen, new PointF(startx+165*fold, starty+(50 + 900 * (float)(scaf.cm2 / maxcm))*fold), new PointF(startx+185*fold, starty+(50 + 900 * (float)(scaf.cm2 / maxcm))*fold));
        }
        public void drawChrBase(Graphics g, long maxbp, double maxcm, int fold, int num_ls, long each_maxbp, double each_maxcm, int startx, int starty){
            g.FillRectangle(Brushes.White, startx, starty, 250*fold, 1000*fold);
            var pen = new Pen(Color.Black, (float)1.5*fold);
            // g.DrawLines(pen, new[] {
            //             new PointF(startx+25*fold,starty+50*fold),
            //             new PointF(startx+30*fold,starty+47*fold),
            //             new PointF(startx+35*fold,starty+(float)44.5*fold),
            //             new PointF(startx+40*fold,starty+(float)43.5*fold),
            //             new PointF(startx+45*fold,starty+(float)43*fold),
            //             new PointF(startx+50*fold,starty+(float)42.5*fold),
            //             new PointF(startx+55*fold,starty+(float)43*fold),
            //             new PointF(startx+60*fold,starty+(float)43.5*fold),
            //             new PointF(startx+65*fold,starty+(float)44.5*fold),
            //             new PointF(startx+70*fold,starty+47*fold),
            //             new PointF(startx+75*fold,starty+50*fold),
            //             new PointF(startx+75*fold,starty+(50+900*each_maxbp/maxbp)*fold),
            //             new PointF(startx+70*fold,starty+(53+900*each_maxbp/maxbp)*fold),
            //             new PointF(startx+65*fold,starty+((float)55.5+900*each_maxbp/maxbp)*fold),
            //             new PointF(startx+60*fold,starty+((float)56.5+900*each_maxbp/maxbp)*fold),
            //             new PointF(startx+55*fold,starty+((float)57+900*each_maxbp/maxbp)*fold),
            //             new PointF(startx+50*fold,starty+((float)57.5+900*each_maxbp/maxbp)*fold),
            //             new PointF(startx+45*fold,starty+((float)57+900*each_maxbp/maxbp)*fold),
            //             new PointF(startx+40*fold,starty+((float)56.5+900*each_maxbp/maxbp)*fold),
            //             new PointF(startx+35*fold,starty+((float)55.5+900*each_maxbp/maxbp)*fold),
            //             new PointF(startx+30*fold,starty+(53+900*each_maxbp/maxbp)*fold),
            //             new PointF(startx+25*fold,starty+(50+900*each_maxbp/maxbp)*fold),
            //             new PointF(startx+25*fold,starty+50*fold)
            //         });
            pen = new Pen(Color.Gray, fold);
            //g.DrawLine(pen, new PointF(startx+100*fold, starty+50*fold), new PointF(startx+100*fold, starty+(50 + 900 * (float)each_maxbp / maxbp)*fold));
            g.DrawLine(pen, new PointF(startx+175*fold, starty+50*fold), new PointF(startx+175*fold, starty+(50 + 900 * (float)(each_maxcm / maxcm))*fold));

        }
        public static string getRevCompBp(string bp){
            switch (bp.ToUpper())
            {
                case "A":
                    return "T";
                case "C":
                    return "G";
                case "G":
                    return "C";
                case "T":
                    return "A";
                default:
                    return "N";
            }
        }
        public static string getRevComp(string seq){
            StringBuilder sb1 = new StringBuilder();
            for(int i=seq.Length-1;i>=0;i--){
                sb1.Append(getRevCompBp(seq.Substring(i,1)));
            }
            return sb1.ToString();
        }
        public scafpos createScafPos(string chr, string order, long pos1, long pos2, double cm1, double cm2)
        {
            scafpos tempscafpos = new scafpos();
            tempscafpos.chr = chr;
            tempscafpos.order = order;
            tempscafpos.pos1 = pos1;
            tempscafpos.pos2 = pos2;
            tempscafpos.cm1 = cm1;
            tempscafpos.cm2 = cm2;
            return tempscafpos;
        }
        public List<stranded_chr> searchFor(string chr, string pos, List<chain> mainchains, List<string> flagpassed){
            List<stranded_chr> result = new List<stranded_chr>();
            chain[] tempchain = mainchains.Where(f => f.chr1==chr && f.pos1==pos).ToArray();
            // if(chr=="scaffold_5"){
            //     bool test = flagpassed.Contains("scaffold_70");
            //     int temp=10;
            // }
            if(tempchain.Length==1){
                if(!flagpassed.Contains(tempchain[0].chr2)){
                    flagpassed.Add(tempchain[0].chr2);
                    string temppos;
                    int temporder;
                    if(tempchain[0].pos2=="start"){
                        temppos="end";
                        temporder=1;
                    }else{
                        temppos="start";
                        temporder=-1;
                    }
                    stranded_chr tempstchr= new stranded_chr();
                    tempstchr.chr=tempchain[0].chr2;
                    tempstchr.order=temporder;
                    tempstchr.match=tempchain[0].match;
                    result.Add(tempstchr);
                    List<stranded_chr> tempres = searchFor(tempchain[0].chr2, temppos, mainchains, flagpassed);
                    result.AddRange(tempres);
                }
            }
            return result;
        }
        public List<stranded_chr> searchRev(string chr, string pos, List<chain> mainchains, List<string> flagpassed){
            List<stranded_chr> result = new List<stranded_chr>();
            chain[] tempchain = mainchains.Where(f => f.chr1==chr && f.pos1==pos).ToArray();
            if(tempchain.Length==1){
                if(!flagpassed.Contains(tempchain[0].chr2)){
                    flagpassed.Add(tempchain[0].chr2);
                    string temppos;
                    int temporder;
                    if(tempchain[0].pos2=="start"){
                        temppos="end";
                        temporder=-1;
                    }else{
                        temppos="start";
                        temporder=1;
                    }
                    stranded_chr tempstchr= new stranded_chr();
                    tempstchr.chr=tempchain[0].chr2;
                    tempstchr.order=temporder;
                    tempstchr.match=tempchain[0].match;
                    List<stranded_chr> tempres = searchRev(tempchain[0].chr2, temppos, mainchains, flagpassed);
                    result.AddRange(tempres);
                    result.Add(tempstchr);
                }
            }
            return result;
        }

        public static edge get_max_edge(string chr, List<string> croschr, Dictionary<string, SortedDictionary<int, Dictionary<string, int[]>>> datas, double opt_nonzerophase, int num_member){
            edge result = new edge();
            chr_match temp1 = calc_match_rate(chr, "start", "start", croschr, datas, opt_nonzerophase, num_member);
            chr_match temp2 = calc_match_rate(chr, "start", "end", croschr, datas, opt_nonzerophase, num_member);
            int i=0;
            double max=0;
            string maxchr="";
            string maxpos="";
            foreach(string tempchr in croschr){
                if(temp1.match[i]>max){
                    max=temp1.match[i];
                    maxchr=tempchr;
                    maxpos="start";
                }
                if(temp2.match[i]>max){
                    max=temp2.match[i];
                    maxchr=tempchr;
                    maxpos="end";
                }
                i++;
            }
            result.chr=maxchr;
            result.pos=maxpos;
            result.prechr=chr;
            result.match=max;
            return result;
        }
        public static chr_match calc_match_rate(string chr, string xpos, string ypos, List<string> croschr, Dictionary<string, SortedDictionary<int, Dictionary<string, int[]>>> datas, double opt_nonzerophase, int num_member){
            chr_match resmatch = new chr_match();
            List<double> res = new List<double>();
            SortedDictionary<int, Dictionary<string, int[]>> tempx = datas[chr];
            foreach(string tempchr in croschr){
                // if(chr=="scaffold_70" && tempchr=="scaffold_5"){
                //     double temptemp=match_rate(tempx, datas[tempchr], xpos, ypos, opt_nonzerophase, num_member);
                //     Console.Write(match_rate(tempx, datas[tempchr], xpos, ypos, opt_nonzerophase, num_member));
                // }
                res.Add(match_rate(tempx, datas[tempchr], xpos, ypos, opt_nonzerophase, num_member));
            }
            resmatch.chr=chr;
            resmatch.match=res;
            return resmatch;
        }
        public static double match_rate(SortedDictionary<int, Dictionary<string, int[]>> x, SortedDictionary<int, Dictionary<string, int[]>> y, string xpos, string ypos, double opt_nonzerophase, int num_member){
            int num_ind=0;
            int num_match=0;
            foreach(KeyValuePair<int, Dictionary<string, int[]>> pair in x){
                if(y.ContainsKey(pair.Key)){
                    int[] temp = get_dist(pair.Value[xpos], y[pair.Key][ypos]);
                    num_ind+= temp[1];
                    num_match+=temp[0];
                }
            }
            if(num_ind < num_member*opt_nonzerophase){
                return 0;
            }else{
                return num_match/(double)num_ind;
            }
        }
        public static int[] get_dist(int[] x, int[] y)
        {
            int n = 0;
            int nm = 0;
            int nn = 0;
            for (int i = 0; i < x.Length; i++)
            {
                if (x[i] != -1 && y[i] != -1)
                {
                    n++;
                    if (x[i] == y[i]) { nm++; } else { nn++; };
                }
            }
            if (n == 0)
            {
                return new int[]{0, 0};
            }
            else
            {
                if (nm >= nn)
                {
                    return new int[]{nm, n};
                }
                else
                {
                    return new int[]{nn, n};
                }
            }
        }

    }
}