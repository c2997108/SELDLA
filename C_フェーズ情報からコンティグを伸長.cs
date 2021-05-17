using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Drawing.Drawing2D;
using System.Drawing;

namespace SELDLA{
    class C_フェーズ情報からコンティグを伸長{
        public struct contig_match{
            public string contig;
            public List<double> I_LIST_一致率;
        }
        public struct S_コンティグ1StartOrEndと2StartOrEndの一致率{
            public string contig1;
            public string contig2;
            public string start_end_1;
            public string start_end_2;
            public double v_一致率;
        }
        public struct stranded_contig{
            public string contig;
            public int orient;
            public double match;
        }
        public struct edge{
            public string chr;
            public string pos;
            public string prechr;
            public double match;
        }
        public struct S_連鎖群中の向きと位置{
            public string V_連鎖群名;
            public string V_向き;
            public long pos1;
            public long pos2;
            public double cm1;
            public double cm2;
        }
        public struct linkagescaf{
            public long length;
            public List<S_連鎖群中の向きと位置> I_LIST_連鎖群中の向きと位置;
        }
        public struct S_伸長したcontigたちの平均フェーズ
        {
            public double V_不一致度;
            public double V_有効個体数割合;
            public SortedDictionary<int, long[]> I_伸長した合計塩基数;
            public SortedDictionary<int, double[]> I_家系別のコンティグの平均フェーズ配列;
        }

        int nbp=10*1000;
        public void run(Dictionary<string, string> I_コンティグ名から塩基配列への連想配列, Dictionary<string, SortedDictionary<int, Dictionary<string, int[]>>> I_コンティグ名から_家系IDから_StartEndから各個人のジェノタイプ情報,
            double opt_s, double opt_nonzerophase, int num_member, string opt_o, double opt_shiftFromCenter, double opt_nonzerophasetotal)
        {
            List<string> I_LIST_コンティグ両端のフェーズが異なるコンティグ名リスト=new List<string>();
            Dictionary<string, double> I_DICT_コンティグのStartEnd間の一致度 = new Dictionary<string, double>();
            //コンティグ両端のフェーズが異なるコンティグを抽出
            List<string> I_LIST_コンティグ両端のフェーズが一致するコンティグリスト = new List<string>();
            foreach (KeyValuePair<string, SortedDictionary<int, Dictionary<string, int[]>>> pair in I_コンティグ名から_家系IDから_StartEndから各個人のジェノタイプ情報){
                double tempmatch = match_rate(pair.Value, pair.Value, "start", "end", opt_nonzerophase, num_member);
                if(tempmatch!=1 && tempmatch>0){
                    I_LIST_コンティグ両端のフェーズが異なるコンティグ名リスト.Add(pair.Key);
                    I_DICT_コンティグのStartEnd間の一致度.Add(pair.Key, tempmatch);
                }
                else
                {
                    I_LIST_コンティグ両端のフェーズが一致するコンティグリスト.Add(pair.Key);
                }
            }

            //scaffoldの両端のフェーズが異なるscaffold間の一致率計算
            List<SortableMatch> sortedmatch = new List<SortableMatch>();
            List<contig_match> mainmatchss = I_LIST_コンティグ両端のフェーズが異なるコンティグ名リスト
                                                .AsParallel()
                                                .Select(f => calc_match_rate(f, "start", "start", I_LIST_コンティグ両端のフェーズが異なるコンティグ名リスト, I_コンティグ名から_家系IDから_StartEndから各個人のジェノタイプ情報, opt_nonzerophase, num_member)).ToList();
            foreach(contig_match tempmatch in mainmatchss){
                int i=0;
                foreach(string tempchr in I_LIST_コンティグ両端のフェーズが異なるコンティグ名リスト){
                    if (tempmatch.contig != tempchr)
                    {
                        SortableMatch temp = new SortableMatch(I_コンティグ名から塩基配列への連想配列[tempmatch.contig].Length + I_コンティグ名から塩基配列への連想配列[tempchr].Length
                         , tempmatch.I_LIST_一致率[i], tempmatch.contig, tempchr, "start", "start");
                        sortedmatch.Add(temp);
                    }
                    i++;
                }
            }
            List<contig_match> mainmatchse = I_LIST_コンティグ両端のフェーズが異なるコンティグ名リスト
                                                .AsParallel()
                                                .Select(f => calc_match_rate(f, "start", "end", I_LIST_コンティグ両端のフェーズが異なるコンティグ名リスト, I_コンティグ名から_家系IDから_StartEndから各個人のジェノタイプ情報, opt_nonzerophase, num_member)).ToList();
            foreach(contig_match tempmatch in mainmatchse){
                int i=0;
                foreach(string tempchr in I_LIST_コンティグ両端のフェーズが異なるコンティグ名リスト){
                    if (tempmatch.contig != tempchr)
                    {
                        SortableMatch temp = new SortableMatch(I_コンティグ名から塩基配列への連想配列[tempmatch.contig].Length + I_コンティグ名から塩基配列への連想配列[tempchr].Length
                         , tempmatch.I_LIST_一致率[i], tempmatch.contig, tempchr, "start", "end");
                        sortedmatch.Add(temp);
                    }
                    i++;
                }
            }
            List<contig_match> mainmatchee = I_LIST_コンティグ両端のフェーズが異なるコンティグ名リスト
                                                .AsParallel()
                                                .Select(f => calc_match_rate(f, "end", "end", I_LIST_コンティグ両端のフェーズが異なるコンティグ名リスト, I_コンティグ名から_家系IDから_StartEndから各個人のジェノタイプ情報, opt_nonzerophase, num_member)).ToList();
            foreach(contig_match tempmatch in mainmatchee){
                int i=0;
                foreach(string tempchr in I_LIST_コンティグ両端のフェーズが異なるコンティグ名リスト){
                    if (tempmatch.contig != tempchr)
                    {
                        SortableMatch temp = new SortableMatch(I_コンティグ名から塩基配列への連想配列[tempmatch.contig].Length + I_コンティグ名から塩基配列への連想配列[tempchr].Length
                         , tempmatch.I_LIST_一致率[i], tempmatch.contig, tempchr, "end", "end");
                        sortedmatch.Add(temp);
                    }
                    i++;
                }
            }

            sortedmatch.Sort();
            sortedmatch.Reverse();

            Dictionary<(string, string, string, string), double> I_Dict_コンティグ1名前とStartEndとコンティグ2名前とStartEndから一致度 = new Dictionary<(string, string, string, string), double>();
            sortedmatch.ForEach(x=>I_Dict_コンティグ1名前とStartEndとコンティグ2名前とStartEndから一致度.Add((x.contig1,x.start_end_1,x.contig2,x.start_end_2),x.V_一致率)); //contig1:end -> contig2:start方向は欠損したデータになっている
            
            List<S_コンティグ1StartOrEndと2StartOrEndの一致率> mainchains = new List<S_コンティグ1StartOrEndと2StartOrEndの一致率>();
            List<string> flagmainstart = new List<string>();
            List<string> flagmainend = new List<string>();
            Dictionary<string, List<string>> flagmain = new Dictionary<string, List<string>>();
            flagmain.Add("start", flagmainstart);
            flagmain.Add("end", flagmainend);
            for(int i=0; i<sortedmatch.Count; i++){
                if(sortedmatch[i].V_一致率<opt_s){
                    break;
                }
                bool isdo = true;
                if(flagmain[sortedmatch[i].start_end_1].Contains(sortedmatch[i].contig1)){
                    isdo=false;
                }else if(flagmain[sortedmatch[i].start_end_2].Contains(sortedmatch[i].contig2)){
                    isdo=false;
                }
                if(isdo){
                    S_コンティグ1StartOrEndと2StartOrEndの一致率[] tempall = mainchains.Where(f=> (f.contig1 == sortedmatch[i].contig1 && f.contig2 == sortedmatch[i].contig2)
                                          || (f.contig1 == sortedmatch[i].contig2 && f.contig2 == sortedmatch[i].contig1)).ToArray();
                    if(tempall.Length>0){
                        isdo=false;
                    }
                }
                if(isdo){
                    flagmain[sortedmatch[i].start_end_1].Add(sortedmatch[i].contig1);
                    flagmain[sortedmatch[i].start_end_2].Add(sortedmatch[i].contig2);
                    S_コンティグ1StartOrEndと2StartOrEndの一致率 tempchain1 = new S_コンティグ1StartOrEndと2StartOrEndの一致率();
                    tempchain1.contig1=sortedmatch[i].contig1;
                    tempchain1.contig2=sortedmatch[i].contig2;
                    tempchain1.start_end_1=sortedmatch[i].start_end_1;
                    tempchain1.start_end_2=sortedmatch[i].start_end_2;
                    tempchain1.v_一致率=sortedmatch[i].V_一致率;
                    mainchains.Add(tempchain1);
                    S_コンティグ1StartOrEndと2StartOrEndの一致率 tempchain2 = new S_コンティグ1StartOrEndと2StartOrEndの一致率();
                    tempchain2.contig2=sortedmatch[i].contig1;
                    tempchain2.contig1=sortedmatch[i].contig2;
                    tempchain2.start_end_2=sortedmatch[i].start_end_1;
                    tempchain2.start_end_1=sortedmatch[i].start_end_2;
                    tempchain2.v_一致率=sortedmatch[i].V_一致率;
                    mainchains.Add(tempchain2);
                }
            }

            List<edge> extramatchse = I_LIST_コンティグ両端のフェーズが一致するコンティグリスト
                                            .AsParallel()
                                            .Select(f => get_max_edge(f, I_LIST_コンティグ両端のフェーズが異なるコンティグ名リスト, I_コンティグ名から_家系IDから_StartEndから各個人のジェノタイプ情報, opt_nonzerophase, num_member)).ToList();

            //scaffoldの並び順を探索する
            List<string> flagpassed = new List<string>();
            int num_ls=0;
            List<linkagescaf> finallinks = new List<linkagescaf>();
            foreach(string temp_contig_name in I_LIST_コンティグ両端のフェーズが異なるコンティグ名リスト){
                if(!flagpassed.Contains(temp_contig_name)){
                    flagpassed.Add(temp_contig_name);
                    List<stranded_contig> stchrs = new List<stranded_contig>();

                    S_伸長したcontigたちの平均フェーズ I_平均フェーズ = new S_伸長したcontigたちの平均フェーズ();
                    I_平均フェーズ.I_伸長した合計塩基数 = new SortedDictionary<int, long[]>();
                    I_平均フェーズ.I_家系別のコンティグの平均フェーズ配列 = new SortedDictionary<int, double[]>();
                    foreach (KeyValuePair<int, Dictionary<string, int[]>> I_計算前の家系ごとのコンティグの平均フェーズ配列 in I_コンティグ名から_家系IDから_StartEndから各個人のジェノタイプ情報[temp_contig_name])
                    {
                        int num_家系内の個体数 = I_計算前の家系ごとのコンティグの平均フェーズ配列.Value["start"].Length;
                        I_平均フェーズ.I_伸長した合計塩基数.Add(I_計算前の家系ごとのコンティグの平均フェーズ配列.Key, new long[num_家系内の個体数]);
                        I_平均フェーズ.I_家系別のコンティグの平均フェーズ配列.Add(I_計算前の家系ごとのコンティグの平均フェーズ配列.Key, new double[num_家系内の個体数]);
                    }
                    //Console.WriteLine(temp_contig_name);
                    List<stranded_contig> tempstchrs = searchRev(temp_contig_name, "start", mainchains, flagpassed, I_コンティグ名から塩基配列への連想配列, I_コンティグ名から_家系IDから_StartEndから各個人のジェノタイプ情報, I_平均フェーズ, opt_shiftFromCenter, opt_nonzerophasetotal);
                    double oldmatch = 1; 
                    if (tempstchrs != null)
                    {
                        foreach (stranded_contig temp in tempstchrs) {
                            stranded_contig newtemp = new stranded_contig();
                            newtemp.contig = temp.contig;
                            newtemp.orient = temp.orient;
                            newtemp.match = oldmatch;
                            stchrs.Add(newtemp);
                            oldmatch = temp.match;
                        }
                    }
                    stranded_contig tempstchr = new stranded_contig();
                    tempstchr.contig = temp_contig_name;
                    tempstchr.orient = 1;
                    tempstchr.match = oldmatch;
                    stchrs.Add(tempstchr);

                    I_平均フェーズ = new S_伸長したcontigたちの平均フェーズ();
                    I_平均フェーズ.I_伸長した合計塩基数 = new SortedDictionary<int, long[]>();
                    I_平均フェーズ.I_家系別のコンティグの平均フェーズ配列 = new SortedDictionary<int, double[]>();
                    foreach (KeyValuePair<int, Dictionary<string, int[]>> I_計算前の家系ごとのコンティグの平均フェーズ配列 in I_コンティグ名から_家系IDから_StartEndから各個人のジェノタイプ情報[temp_contig_name])
                    {
                        int num_家系内の個体数 = I_計算前の家系ごとのコンティグの平均フェーズ配列.Value["start"].Length;
                        I_平均フェーズ.I_伸長した合計塩基数.Add(I_計算前の家系ごとのコンティグの平均フェーズ配列.Key, new long[num_家系内の個体数]);
                        I_平均フェーズ.I_家系別のコンティグの平均フェーズ配列.Add(I_計算前の家系ごとのコンティグの平均フェーズ配列.Key, new double[num_家系内の個体数]);
                    }
                    //Console.WriteLine("for:");
                    List<stranded_contig> tempstchrs2 = searchFor(temp_contig_name, "end", mainchains, flagpassed, I_コンティグ名から塩基配列への連想配列, I_コンティグ名から_家系IDから_StartEndから各個人のジェノタイプ情報, I_平均フェーズ, opt_shiftFromCenter, opt_nonzerophasetotal);
                    if (tempstchrs2 != null)
                    {
                        stchrs.AddRange(tempstchrs2);
                    }

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
                    finallink.I_LIST_連鎖群中の向きと位置 = new List<S_連鎖群中の向きと位置>();
                    foreach(stranded_contig stchr in stchrs){
                        num_ordered_sca++;
                        if(num_ordered_sca==1){
                            double tempmatch=-1;
                            List<edge> tempedges = new List<edge>();
                            if (stchr.orient == 1)
                            {
                                tempedges = extramatchse.Where(x => x.chr == stchr.contig && x.match >= opt_s && x.pos == "start").OrderBy(x => x.match).ToList();
                            }
                            else
                            {
                                tempedges = extramatchse.Where(x => x.chr == stchr.contig && x.match >= opt_s && x.pos == "end").OrderBy(x => x.match).ToList();
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
                                finallink.length=(tempdist+I_コンティグ名から塩基配列への連想配列[tempedge.prechr].Length);
                                finallink.I_LIST_連鎖群中の向きと位置.Add(createScafPos(tempedge.prechr, "na", tempdist, (tempdist+I_コンティグ名から塩基配列への連想配列[tempedge.prechr].Length), tempdistcm, tempdistcm));
                                //Console.WriteLine(num_ls + "\t" + tempedge.prechr + "\tna\t" + tempdist+"\t"+ (tempdist+refseqs2[tempedge.prechr].Length) + "\t" + tempdistcm + "\t" + tempdistcm
                                // + "\t" + tempedge.chr + " " + tempedge.pos + " " + tempedge.match);
                                tempmatch = tempedge.match;
                                tempdist+=I_コンティグ名から塩基配列への連想配列[tempedge.prechr].Length;
                            }
                            if(tempmatch==-1){
                                tempdistcm=0;
                                tempdist=0;
                            }else{
                                tempdistcm+=(1-tempmatch)*100;
                                tempdist+=nbp;
                            }
                            string temporder;
                            if(stchr.orient==1){
                                temporder="+";
                            }else{
                                temporder="-";
                            }
                            finallink.length=(tempdist+I_コンティグ名から塩基配列への連想配列[stchr.contig].Length);
                            finallink.I_LIST_連鎖群中の向きと位置.Add(createScafPos(stchr.contig, temporder, tempdist, (tempdist+I_コンティグ名から塩基配列への連想配列[stchr.contig].Length), tempdistcm, (tempdistcm+(1-I_DICT_コンティグのStartEnd間の一致度[stchr.contig])*100)));
                            //Console.WriteLine(num_ls + "\t"+ stchr.chr+"\t"+temporder+"\t"+tempdist + "\t"+ (tempdist+refseqs2[stchr.chr].Length)+"\t"+tempdistcm+"\t"+(tempdistcm+(1-matchcroschr[stchr.chr])*100)
                            //+"\tmatch_rate_of_previous_ordered_scaffold:"+stchr.match+" match_rate_between_both_edge_of_this_scaffold"+matchcroschr[stchr.chr]);
                            tempdist+=I_コンティグ名から塩基配列への連想配列[stchr.contig].Length;
                            tempdistcm+=(1-I_DICT_コンティグのStartEnd間の一致度[stchr.contig])*100;
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
                            if(stchr.orient==1){
                                tempedges=extramatchse.Where(x => x.chr == stchr.contig && x.match >= opt_s && x.pos == "start").OrderBy(x => x.match).ToList();
                            }else{
                                tempedges=extramatchse.Where(x => x.chr == stchr.contig && x.match >= opt_s && x.pos == "end").OrderBy(x => x.match).ToList();
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
                                finallink.length=(tempdist+I_コンティグ名から塩基配列への連想配列[tempedge.prechr].Length);
                                finallink.I_LIST_連鎖群中の向きと位置.Add(createScafPos(tempedge.prechr, "na", tempdist, (tempdist+I_コンティグ名から塩基配列への連想配列[tempedge.prechr].Length), tempdistcm, tempdistcm));
                                //Console.WriteLine(num_ls + "\t" + tempedge.prechr + "\tna\t" + tempdist+"\t"+ (tempdist+refseqs2[tempedge.prechr].Length) + "\t" + tempdistcm + "\t" + tempdistcm
                                // + "\t" + tempedge.chr + " " + tempedge.pos + " " + tempedge.match);
                                tempdist+=I_コンティグ名から塩基配列への連想配列[tempedge.prechr].Length;
                            }
                            double newtempdistcm = temp_prev_distcm+(1-stchr.match)*100;
                            foreach (edge tempedge in tempedges)
                            {
                                tempdistcm=newtempdistcm-(1-tempedge.match)*100*matchnorm;
                                tempdist+=nbp;
                                finallink.length=(tempdist+I_コンティグ名から塩基配列への連想配列[tempedge.prechr].Length);
                                finallink.I_LIST_連鎖群中の向きと位置.Add(createScafPos(tempedge.prechr, "na", tempdist, (tempdist+I_コンティグ名から塩基配列への連想配列[tempedge.prechr].Length), tempdistcm, tempdistcm));
                                //Console.WriteLine(num_ls + "\t" + tempedge.prechr + "\tna\t" + tempdist+"\t"+ (tempdist+refseqs2[tempedge.prechr].Length) + "\t" + tempdistcm + "\t" + tempdistcm
                                // + "\t" + tempedge.chr + " " + tempedge.pos + " " + tempedge.match);
                                tempdist+=I_コンティグ名から塩基配列への連想配列[tempedge.prechr].Length;
                            }
                            tempdistcm=temp_prev_distcm+(1-stchr.match)*100;
                            tempdist+=nbp;
                            string temporder;
                            if(stchr.orient==1){
                                temporder="+";
                            }else{
                                temporder="-";
                            }
                            finallink.length=(tempdist+I_コンティグ名から塩基配列への連想配列[stchr.contig].Length);
                            finallink.I_LIST_連鎖群中の向きと位置.Add(createScafPos(stchr.contig, temporder, tempdist, (tempdist+I_コンティグ名から塩基配列への連想配列[stchr.contig].Length), tempdistcm, (tempdistcm+(1-I_DICT_コンティグのStartEnd間の一致度[stchr.contig])*100)));
                            //Console.WriteLine(num_ls + "\t"+ stchr.chr+"\t"+temporder+"\t"+tempdist + "\t"+ (tempdist+refseqs2[stchr.chr].Length)+"\t"+tempdistcm+"\t"+(tempdistcm+(1-matchcroschr[stchr.chr])*100)
                            //+"\tmatch_rate_of_previous_ordered_scaffold:"+stchr.match+" match_rate_between_both_edge_of_this_scaffold"+matchcroschr[stchr.chr]);
                            tempdist+=I_コンティグ名から塩基配列への連想配列[stchr.contig].Length;
                            tempdistcm+=(1-I_DICT_コンティグのStartEnd間の一致度[stchr.contig])*100;
                            temp_prev_distcm=tempdistcm;
                        }
                        temp_prev_chr=stchr.contig;
                        temp_prev_order=stchr.orient;

                        if(num_ordered_sca == stchrs.Count){
                            double tempmatch=-1;
                            List<edge> tempedges = new List<edge>();
                            if (stchr.orient == 1)
                            {
                                tempedges = extramatchse.Where(x => x.chr == stchr.contig && x.match >= opt_s && x.pos == "end").OrderByDescending(x => x.match).ToList();
                            }
                            else
                            {
                                tempedges = extramatchse.Where(x => x.chr == stchr.contig && x.match >= opt_s && x.pos == "start").OrderByDescending(x => x.match).ToList();
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
                                finallink.length=(tempdist+I_コンティグ名から塩基配列への連想配列[tempedge.prechr].Length);
                                finallink.I_LIST_連鎖群中の向きと位置.Add(createScafPos(tempedge.prechr, "na", tempdist, (tempdist+I_コンティグ名から塩基配列への連想配列[tempedge.prechr].Length), tempdistcm, tempdistcm));
                                //Console.WriteLine(num_ls + "\t" + tempedge.prechr + "\tna\t" + tempdist+"\t"+ (tempdist+refseqs2[tempedge.prechr].Length) + "\t" + tempdistcm + "\t" + tempdistcm
                                // + "\t" + tempedge.chr + " " + tempedge.pos + " " + tempedge.match);
                                tempmatch = tempedge.match;
                                tempdist+=I_コンティグ名から塩基配列への連想配列[tempedge.prechr].Length;
                            }
                        }
                    }
                    finallinks.Add(finallink);
                }
            }
            if(finallinks.Count==0){
                throw new Exception("There were no crossed chromosomes.");
            }

            //結果ファイルの書き出し
            List<linkagescaf> listlinkages = finallinks.OrderByDescending(x=> x.length).ToList();
            double maxcm = listlinkages.Max(x=> x.I_LIST_連鎖群中の向きと位置[x.I_LIST_連鎖群中の向きと位置.Count-1].cm2);
            long maxbp = listlinkages.Max(x=> x.length);
            int num_big_ls = listlinkages.Where(x => x.length>1000*1000).Count();
            num_ls=0;
            StreamWriter writer = new StreamWriter(opt_o+"_chain.txt");
            StreamWriter writer_chph = new StreamWriter(opt_o + "_chain_ph.txt");
            writer_chph.Write("#chr\tname\tdirection\tstart_bp\tend_bp\tgenetical_start_cM\tgenetical_end_cM");
            I_コンティグ名から_家系IDから_StartEndから各個人のジェノタイプ情報.Take(1).ToList()[0].Value.ToList().ForEach(x => { for (int i = 0; i < x.Value["start"].Length; i++) { writer_chph.Write("\tfamily_" + x.Key); } });
            writer_chph.WriteLine();
            StreamWriter writerfa = new StreamWriter(opt_o+"_extended.fasta");
            StreamWriter writerunloc = new StreamWriter(opt_o+"_unoriented.fasta");
            StreamWriter writerfanainchr = new StreamWriter(opt_o+"_include_unoriented_in_chr.fasta");
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
            List<(int, string, string)> matchRateOrder = new List<(int, string, string)>();
            List<(int, string)> scafLenInChrList = new List<(int, string)>();
            int totalScafLen = 0;
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
                foreach (S_連鎖群中の向きと位置 scaf in scafs.I_LIST_連鎖群中の向きと位置) {
                    num_scaf++;
                    res_num_scaf++;
                    res_bp_scaf += I_コンティグ名から塩基配列への連想配列[scaf.V_連鎖群名].Length;
                    res_num_scaf_loc++;
                    res_bp_scaf_loc += I_コンティグ名から塩基配列への連想配列[scaf.V_連鎖群名].Length;
                    //Console.WriteLine(num_ls+"\t"+ scaf.chr+"\t"+scaf.order+"\t"+scaf.pos1+"\t"+scaf.pos2+"\t"+scaf.cm1+"\t"+scaf.cm2);
                    writer.WriteLine(num_ls + "\t" + scaf.V_連鎖群名 + "\t" + scaf.V_向き + "\t" + scaf.pos1 + "\t" + scaf.pos2 + "\t" + scaf.cm1 + "\t" + scaf.cm2);
                    //フェーズ付きチェインファイル作成
                    if (scaf.V_向き == "+" || scaf.V_向き == "na")
                    {
                        writer_chph.Write(num_ls + "\t" + scaf.V_連鎖群名 + "\t" + scaf.V_向き + "\t" + scaf.pos1 + "\t" + scaf.pos2 + "\t" + scaf.cm1 + "\t" + scaf.cm2);
                        I_コンティグ名から_家系IDから_StartEndから各個人のジェノタイプ情報[scaf.V_連鎖群名].ToList().ForEach(x => { x.Value["start"].ToList().ForEach(y => writer_chph.Write("\t" + y)); });
                        writer_chph.WriteLine();
                        writer_chph.Write(num_ls + "\t" + scaf.V_連鎖群名 + "\t" + scaf.V_向き + "\t" + scaf.pos1 + "\t" + scaf.pos2 + "\t" + scaf.cm1 + "\t" + scaf.cm2);
                        I_コンティグ名から_家系IDから_StartEndから各個人のジェノタイプ情報[scaf.V_連鎖群名].ToList().ForEach(x => { x.Value["end"].ToList().ForEach(y => writer_chph.Write("\t" + y)); });
                        writer_chph.WriteLine();
                    }
                    else
                    {
                        writer_chph.Write(num_ls + "\t" + scaf.V_連鎖群名 + "\t" + scaf.V_向き + "\t" + scaf.pos1 + "\t" + scaf.pos2 + "\t" + scaf.cm1 + "\t" + scaf.cm2);
                        I_コンティグ名から_家系IDから_StartEndから各個人のジェノタイプ情報[scaf.V_連鎖群名].ToList().ForEach(x => { x.Value["end"].ToList().ForEach(y => writer_chph.Write("\t" + y)); });
                        writer_chph.WriteLine();
                        writer_chph.Write(num_ls + "\t" + scaf.V_連鎖群名 + "\t" + scaf.V_向き + "\t" + scaf.pos1 + "\t" + scaf.pos2 + "\t" + scaf.cm1 + "\t" + scaf.cm2);
                        I_コンティグ名から_家系IDから_StartEndから各個人のジェノタイプ情報[scaf.V_連鎖群名].ToList().ForEach(x => { x.Value["start"].ToList().ForEach(y => writer_chph.Write("\t" + y)); });
                        writer_chph.WriteLine();
                    }
                    //一致率テーブルのためのインデックス作成
                    //ヒートマップのためのscaffoldの長さ保存
                    totalScafLen += I_コンティグ名から塩基配列への連想配列[scaf.V_連鎖群名].Length;
                    if (scaf.V_向き == "+" || scaf.V_向き=="-")
                    {
                        matchRateOrder.Add((num_ls, scaf.V_連鎖群名, scaf.V_向き));
                        scafLenInChrList.Add((I_コンティグ名から塩基配列への連想配列[scaf.V_連鎖群名].Length, "+-"));
                    }
                    else
                    {
                        scafLenInChrList.Add((I_コンティグ名から塩基配列への連想配列[scaf.V_連鎖群名].Length, "na"));
                    }
                    

                    // if (scafs.length > 1000 * 1000)
                    // {
                    //     drawBpCm(g, 2, scaf, maxbp, maxcm, 0, 0);
                    // }
                    if (num_scaf>1){
                        for(int i=0;i<nbp;i++){
                            writerfa.Write("N");
                            writerfanainchr.Write("N");
                        }
                    }
                    if(scaf.V_向き=="+"){
                        writerfa.Write(I_コンティグ名から塩基配列への連想配列[scaf.V_連鎖群名]);
                        writerfanainchr.Write(I_コンティグ名から塩基配列への連想配列[scaf.V_連鎖群名]);
                    }else if(scaf.V_向き=="-"){
                        writerfa.Write(getRevComp(I_コンティグ名から塩基配列への連想配列[scaf.V_連鎖群名]));
                        writerfanainchr.Write(getRevComp(I_コンティグ名から塩基配列への連想配列[scaf.V_連鎖群名]));
                    }else{
                        for(int i=0;i<I_コンティグ名から塩基配列への連想配列[scaf.V_連鎖群名].Length;i++){
                            writerfa.Write("N");
                        }
                        writerfanainchr.Write(I_コンティグ名から塩基配列への連想配列[scaf.V_連鎖群名]);
                    }
                    if(scaf.V_向き=="na"){
                        writerunloc.WriteLine(">linkage_scaffold_"+num_ls+"_unoriented_pos"+scaf.pos1+"_old_"+scaf.V_連鎖群名);
                        writerunloc.WriteLine(I_コンティグ名から塩基配列への連想配列[scaf.V_連鎖群名]);
                        inChain.Add(scaf.V_連鎖群名);
                    }else{
                        res_num_scaf_order++;
                        res_bp_scaf_order+=I_コンティグ名から塩基配列への連想配列[scaf.V_連鎖群名].Length;
                        if(scafs.length>=1000*1000){
                            res_num_scaf_order_big++;
                            res_bp_scaf_order_big+=I_コンティグ名から塩基配列への連想配列[scaf.V_連鎖群名].Length;
                        }
                        inChain.Add(scaf.V_連鎖群名);
                        inChainWithOriented.Add(scaf.V_連鎖群名);
                    }
                }
                writerfa.WriteLine("");
                writerfanainchr.WriteLine("");
                // if (scafs.length > 1000 * 1000)
                // {
                //     image.Save(opt_o + "_map_"+num_ls+".png");
                // }
            }
            foreach(KeyValuePair<string, string> Pair in I_コンティグ名から塩基配列への連想配列){
                if(!inChain.Contains(Pair.Key)){
                    res_num_scaf++;
                    res_bp_scaf+=I_コンティグ名から塩基配列への連想配列[Pair.Key].Length;
                    writerfa.WriteLine(">old_"+Pair.Key);
                    writerfa.WriteLine(Pair.Value);
                    writerfanainchr.WriteLine(">old_"+Pair.Key);
                    writerfanainchr.WriteLine(Pair.Value);
                }
            }
            writer.Close();
            writer_chph.Close();
            writerfa.Close();
            writerunloc.Close();
            writerfanainchr.Close();
            StreamWriter txtwrite = new StreamWriter(opt_o+"_log.txt");
            string str="Total: "+res_num_scaf.ToString("N0")+" scaffolds, Located: "+res_num_scaf_loc.ToString("N0")
            +" ("+((double)res_num_scaf_loc/res_num_scaf*100).ToString("F3")+"), Oriented: "
            +res_num_scaf_order.ToString("N0")+" ("+((double)res_num_scaf_order/res_num_scaf*100).ToString("F3")+"), Oriented in >=1Mbp chr ("+num_big_ls+"): "
            +res_num_scaf_order_big.ToString("N0")+" ("+((double)res_num_scaf_order_big/res_num_scaf*100).ToString("F3")+")";
            Console.WriteLine(str);
            txtwrite.WriteLine(str);

            str="Total: "+res_bp_scaf.ToString("N0")+" bp, Located: "+res_bp_scaf_loc.ToString("N0")
            +" ("+((double)res_bp_scaf_loc/res_bp_scaf*100).ToString("F3")+"), Oriented: "
            +res_bp_scaf_order+" ("+((double)res_bp_scaf_order/res_bp_scaf*100).ToString("F3")+"), Oriented in >=1Mbp chr ("+num_big_ls+"): "
            +res_bp_scaf_order_big+" ("+((double)res_bp_scaf_order_big/res_bp_scaf*100).ToString("F3")+")";
            Console.WriteLine(str);
            txtwrite.WriteLine(str);
            txtwrite.Close();

            //一致率テーブル作成
            StreamWriter writer_mtable = new StreamWriter(opt_o + "_match_rate");
            writer_mtable.Write("#chr_contig_pos");
            List<int> heatmap_chr_border = new List<int>();
            int heatmap_chr_border_temp = -1;
            int heatmap_chr_border_cnt = 0;
            double[,] heatmap_data = new double[matchRateOrder.Count*2, matchRateOrder.Count*2];
            matchRateOrder.ForEach(x =>
            {
                if (x.Item3 == "+")
                {
                    writer_mtable.Write("\t" + x.Item1 + "_" + x.Item2 + "_s\t" + x.Item1 + "_" + x.Item2 + "_e");
                }
                else //x.Item3 == "-"のみ
                {
                    writer_mtable.Write("\t" + x.Item1 + "_" + x.Item2 + "_e\t" + x.Item1 + "_" + x.Item2 + "_s");
                }
                if(x.Item1 != heatmap_chr_border_temp)
                {
                    heatmap_chr_border_temp = x.Item1;
                    heatmap_chr_border.Add(heatmap_chr_border_cnt);
                }
                heatmap_chr_border_cnt++;
            });
            writer_mtable.WriteLine();
            matchRateOrder.ForEach(x =>
            {
                if (x.Item3 == "+")
                {
                    writer_mtable.Write(x.Item1 + "_" + x.Item2 + "_s");
                    matchRateOrder.ForEach(y =>
                    {
                        if(x.Item2==y.Item2)
                        {
                            writer_mtable.Write("\t1\t"+ match_rate_for_ends(I_コンティグ名から_家系IDから_StartEndから各個人のジェノタイプ情報[x.Item2], opt_nonzerophase, num_member));
                        }
                        else if (y.Item3 == "+") {
                            double temp;
                            try { temp = I_Dict_コンティグ1名前とStartEndとコンティグ2名前とStartEndから一致度[(x.Item2, "start", y.Item2, "start")]; } catch (Exception) { temp = 0; }
                            writer_mtable.Write("\t" + temp);
                            try { temp = I_Dict_コンティグ1名前とStartEndとコンティグ2名前とStartEndから一致度[(x.Item2, "start", y.Item2, "end")]; } catch (Exception) { temp = 0; }
                            writer_mtable.Write("\t" + temp);
                        }
                        else
                        {
                            double temp;
                            try { temp = I_Dict_コンティグ1名前とStartEndとコンティグ2名前とStartEndから一致度[(x.Item2, "start", y.Item2, "end")]; } catch (Exception) { temp = 0; }
                            writer_mtable.Write("\t" + temp);
                            try { temp = I_Dict_コンティグ1名前とStartEndとコンティグ2名前とStartEndから一致度[(x.Item2, "start", y.Item2, "start")]; } catch (Exception) { temp = 0; }
                            writer_mtable.Write("\t" + temp);
                        }
                    });
                    writer_mtable.WriteLine();
                    writer_mtable.Write(x.Item1 + "_" + x.Item2 + "_e");
                    matchRateOrder.ForEach(y =>
                    {
                        if (x.Item2 == y.Item2)
                        {
                            writer_mtable.Write("\t"+ match_rate_for_ends(I_コンティグ名から_家系IDから_StartEndから各個人のジェノタイプ情報[x.Item2], opt_nonzerophase, num_member)+"\t1");
                        }
                        else if (y.Item3 == "+")
                        {
                            double temp;
                            try { temp = I_Dict_コンティグ1名前とStartEndとコンティグ2名前とStartEndから一致度[(y.Item2, "start", x.Item2, "end")]; } catch (Exception) { temp = 0; } //end->startはデータがないため
                            writer_mtable.Write("\t" + temp);
                            try { temp = I_Dict_コンティグ1名前とStartEndとコンティグ2名前とStartEndから一致度[(x.Item2, "end", y.Item2, "end")]; } catch (Exception) { temp = 0; }
                            writer_mtable.Write("\t" + temp);
                        }
                        else
                        {
                            double temp;
                            try { temp = I_Dict_コンティグ1名前とStartEndとコンティグ2名前とStartEndから一致度[(x.Item2, "end", y.Item2, "end")]; } catch (Exception) { temp = 0; }
                            writer_mtable.Write("\t" + temp);
                            try { temp = I_Dict_コンティグ1名前とStartEndとコンティグ2名前とStartEndから一致度[(y.Item2, "start", x.Item2, "end")]; } catch (Exception) { temp = 0; }
                            writer_mtable.Write("\t" + temp);
                        }
                    });
                    writer_mtable.WriteLine();
                }
                else //x.Item3 == "-"のみ
                {
                    writer_mtable.Write(x.Item1 + "_" + x.Item2 + "_e");
                    matchRateOrder.ForEach(y =>
                    {
                        if (x.Item2 == y.Item2)
                        {
                            writer_mtable.Write("\t1\t"+ match_rate_for_ends(I_コンティグ名から_家系IDから_StartEndから各個人のジェノタイプ情報[x.Item2], opt_nonzerophase, num_member));
                        }
                        else if (y.Item3 == "+")
                        {
                            double temp;
                            try { temp = I_Dict_コンティグ1名前とStartEndとコンティグ2名前とStartEndから一致度[(y.Item2, "start", x.Item2, "end")]; } catch (Exception) { temp = 0; }
                            writer_mtable.Write("\t" + temp);
                            try { temp = I_Dict_コンティグ1名前とStartEndとコンティグ2名前とStartEndから一致度[(x.Item2, "end", y.Item2, "end")]; } catch (Exception) { temp = 0; }
                            writer_mtable.Write("\t" + temp);
                        }
                        else
                        {
                            double temp;
                            try { temp = I_Dict_コンティグ1名前とStartEndとコンティグ2名前とStartEndから一致度[(x.Item2, "end", y.Item2, "end")]; } catch (Exception) { temp = 0; }
                            writer_mtable.Write("\t" + temp);
                            try { temp = I_Dict_コンティグ1名前とStartEndとコンティグ2名前とStartEndから一致度[(y.Item2, "start", x.Item2, "end")]; } catch (Exception) { temp = 0; }
                            writer_mtable.Write("\t" + temp);
                        }
                    });
                    writer_mtable.WriteLine();
                    writer_mtable.Write(x.Item1 + "_" + x.Item2 + "_s");
                    matchRateOrder.ForEach(y =>
                    {
                        if (x.Item2 == y.Item2)
                        {
                            writer_mtable.Write("\t"+ match_rate_for_ends(I_コンティグ名から_家系IDから_StartEndから各個人のジェノタイプ情報[x.Item2], opt_nonzerophase, num_member)+"\t1");
                        }
                        else if (y.Item3 == "+")
                        {
                            double temp;
                            try { temp = I_Dict_コンティグ1名前とStartEndとコンティグ2名前とStartEndから一致度[(x.Item2, "start", y.Item2, "start")]; } catch (Exception) { temp = 0; }
                            writer_mtable.Write("\t" + temp);
                            try { temp = I_Dict_コンティグ1名前とStartEndとコンティグ2名前とStartEndから一致度[(x.Item2, "start", y.Item2, "end")]; } catch (Exception) { temp = 0; }
                            writer_mtable.Write("\t" + temp);
                        }
                        else
                        {
                            double temp;
                            try { temp = I_Dict_コンティグ1名前とStartEndとコンティグ2名前とStartEndから一致度[(x.Item2, "start", y.Item2, "end")]; } catch (Exception) { temp = 0; }
                            writer_mtable.Write("\t" + temp);
                            try { temp = I_Dict_コンティグ1名前とStartEndとコンティグ2名前とStartEndから一致度[(x.Item2, "start", y.Item2, "start")]; } catch (Exception) { temp = 0; }
                            writer_mtable.Write("\t" + temp);
                        }
                    });
                    writer_mtable.WriteLine();
                }
            });
            writer_mtable.Close();

            //グラフ作成　連鎖地図
            Console.WriteLine("Drawing genetical and phisycal map...");
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
                    drawChrBase(gall, maxbp, maxcm, 1, i, listlinkages[i - 1].length, listlinkages[i - 1].I_LIST_連鎖群中の向きと位置[listlinkages[i - 1].I_LIST_連鎖群中の向きと位置.Count - 1].cm2, (j - 1) * per_width, (k - 1) * per_height);
                    foreach (S_連鎖群中の向きと位置 scaf in listlinkages[i - 1].I_LIST_連鎖群中の向きと位置)
                    {
                        drawBpCm(gall, 1, scaf, maxbp, maxcm, (j - 1) * per_width, (k - 1) * per_height);
                    }
                    drawChrFin(gall, maxbp, maxcm, 1, i, listlinkages[i - 1].length, listlinkages[i - 1].I_LIST_連鎖群中の向きと位置[listlinkages[i - 1].I_LIST_連鎖群中の向きと位置.Count - 1].cm2, (j - 1) * per_width, (k - 1) * per_height);
                }
                imageall.Save(opt_o + "_map_all.png");
            }
            catch (Exception e)
            {
                Console.WriteLine(e);
            }

            //グラフ作成　ヒートマップ
            Console.WriteLine("Drawing heatmap...");
            try
            {
                var imageall = new Bitmap(1200,1200);
                var gall = Graphics.FromImage(imageall);
                gall.SmoothingMode = SmoothingMode.AntiAlias;
                gall.InterpolationMode = InterpolationMode.HighQualityBicubic;
                gall.PixelOffsetMode = PixelOffsetMode.HighQuality;
                System.Drawing.SolidBrush myBrush = new System.Drawing.SolidBrush(Color.White);
                gall.FillRectangle(myBrush, 0, 0, 1200, 1200);
                int hmall = matchRateOrder.Count*2;
                float hmunit = 1000 / (float)hmall;

                File.ReadLines(opt_o + "_match_rate").Select((line, y) => (line,y)).ToList().ForEach(z=>
                {
                    if (z.y > 0)
                    {
                        z.line.Split("\t").Skip(1).Select((item,x)=>(item,x)).ToList().ForEach(w=>
                        {
                            double itemtmp = double.Parse(w.item);
                            if (itemtmp < 0.5) { itemtmp = 0.5; }
                            itemtmp = (itemtmp - 0.5) * 2;
                            myBrush = new System.Drawing.SolidBrush(Color.FromArgb((int)(itemtmp*255) , Color.Green));
                            //Console.WriteLine((int)(itemtmp * 255) + " " + (100 + hmunit * w.x) + " " + (100 + hmunit * z.y));
                            gall.FillRectangle(myBrush, 100 + hmunit * w.x, 1100 - hmunit * z.y, hmunit, hmunit);
                        });
                    }
                });
                var pen = new Pen(Color.Black);
                Font myfont = new Font("Tahoma", 20);
                heatmap_chr_border.Select((x,i)=>(x,i)).ToList().ForEach(z =>
                {
                    gall.DrawLine(pen, 100 + z.x * 2 * hmunit, 100, 100 + z.x * 2 * hmunit, 1100);
                    gall.DrawLine(pen, 100, 1100 - z.x * 2 * hmunit, 1100, 1100 - z.x * 2 * hmunit);
                    gall.DrawString((z.i+1).ToString(), myfont, Brushes.Black, new Point(30, (int)(1100 - z.x * 2 * hmunit) - 30));
                    gall.DrawString((z.i+1).ToString(), myfont, Brushes.Black, new Point((int)(100 + z.x * 2 * hmunit) + 10, 1130));
                });
                imageall.Save(opt_o + "_heatmap_phase.png");
            }
            catch (Exception e)
            {
                Console.WriteLine(e);
            }
            //グラフ作成　ヒートマップ物理距離対応版
            try
            {
                var imageall = new Bitmap(2200, 2200);
                var gall = Graphics.FromImage(imageall);
                gall.SmoothingMode = SmoothingMode.AntiAlias;
                gall.InterpolationMode = InterpolationMode.HighQualityBicubic;
                gall.PixelOffsetMode = PixelOffsetMode.HighQuality;
                System.Drawing.SolidBrush myBrush = new System.Drawing.SolidBrush(Color.White);
                System.Drawing.SolidBrush myBrush2 = new System.Drawing.SolidBrush(Color.LightGray);
                gall.FillRectangle(myBrush, 0, 0, 2200, 2200);
                float hmunit = 2000 / (float)totalScafLen;
                double[,] matchAllData = new double[matchRateOrder.Count * 2, matchRateOrder.Count * 2];
                
                File.ReadLines(opt_o + "_match_rate").Select((line, y) => (line, y)).ToList().ForEach(z =>
                {
                    if (z.y > 0)
                    {
                        z.line.Split("\t").Skip(1).Select((item, x) => (item, x)).ToList().ForEach(w =>
                        {
                            //Console.WriteLine(w.x + " " + z.y + " " + matchRateOrder.Count);
                            matchAllData[w.x, z.y - 1] = double.Parse(w.item);
                        });
                    }
                });
                var pen = new Pen(Color.Black);
                Font myfont = new Font("Tahoma", 20);
                int orientedScafNumY = 0;
                int cntScafLenY = 0;
                scafLenInChrList.ForEach(z =>
                {
                    int orientedScafNumX = 0;
                    int cntScafLenX = 0;
                    cntScafLenY += z.Item1;
                    if (z.Item2 == "+-")
                    {
                        orientedScafNumY++;
                        scafLenInChrList.ForEach(w =>
                        {
                            cntScafLenX += w.Item1;
                            if (w.Item2 == "+-")
                            {
                                orientedScafNumX++;
                                double itemtmp=0;
                                //Console.WriteLine((2 * orientedScafNumX - 2) + " " + (2 * orientedScafNumY - 2) + " " + scafLenInChrList.Count);
                                itemtmp = matchAllData[2 * orientedScafNumX - 2, 2 * orientedScafNumY - 2];
                                if (itemtmp < 0.5) { itemtmp = 0.5; }
                                itemtmp = (itemtmp - 0.5) * 2;
                                //myBrush = new System.Drawing.SolidBrush(Color.FromArgb((int)(255 * (1 - itemtmp)), 255, (int)(255 * (1 - itemtmp))));
                                myBrush = new System.Drawing.SolidBrush(Color.FromArgb((int)(itemtmp * 255), Color.Green));
                                gall.FillRectangle(myBrush, 100 + hmunit * (cntScafLenX - w.Item1), 2100 - hmunit * (cntScafLenY - z.Item1 / (float)2), hmunit * w.Item1/(float)2, hmunit * z.Item1/(float)2);

                                itemtmp = matchAllData[2 * orientedScafNumX - 1, 2 * orientedScafNumY - 2];
                                if (itemtmp < 0.5) { itemtmp = 0.5; }
                                itemtmp = (itemtmp - 0.5) * 2;
                                //myBrush = new System.Drawing.SolidBrush(Color.FromArgb((int)(255 * (1 - itemtmp)), 255, (int)(255 * (1 - itemtmp))));
                                myBrush = new System.Drawing.SolidBrush(Color.FromArgb((int)(itemtmp * 255), Color.Green));
                                gall.FillRectangle(myBrush, 100 + hmunit * (cntScafLenX - w.Item1 / (float)2), 2100 - hmunit * (cntScafLenY - z.Item1 / (float)2), hmunit * w.Item1 / (float)2, hmunit * z.Item1 / (float)2);

                                itemtmp = matchAllData[2 * orientedScafNumX - 2, 2 * orientedScafNumY - 1];
                                if (itemtmp < 0.5) { itemtmp = 0.5; }
                                itemtmp = (itemtmp - 0.5) * 2;
                                //myBrush = new System.Drawing.SolidBrush(Color.FromArgb((int)(255 * (1 - itemtmp)), 255, (int)(255 * (1 - itemtmp))));
                                myBrush = new System.Drawing.SolidBrush(Color.FromArgb((int)(itemtmp * 255), Color.Green));
                                gall.FillRectangle(myBrush, 100 + hmunit * (cntScafLenX - w.Item1), 2100 - hmunit * (cntScafLenY), hmunit * w.Item1 / (float)2, hmunit * z.Item1 / (float)2);

                                itemtmp = matchAllData[2 * orientedScafNumX - 1, 2 * orientedScafNumY - 1];
                                if (itemtmp < 0.5) { itemtmp = 0.5; }
                                itemtmp = (itemtmp - 0.5) * 2;
                                //myBrush = new System.Drawing.SolidBrush(Color.FromArgb((int)(255 * (1 - itemtmp)), 255, (int)(255 * (1 - itemtmp))));
                                myBrush = new System.Drawing.SolidBrush(Color.FromArgb((int)(itemtmp * 255), Color.Green));
                                gall.FillRectangle(myBrush, 100 + hmunit * (cntScafLenX - w.Item1 / (float)2), 2100 - hmunit * (cntScafLenY), hmunit * w.Item1 / (float)2, hmunit * z.Item1 / (float)2);
                            }
                            else
                            {
                                gall.FillRectangle(myBrush2, 100 + hmunit * (cntScafLenX - w.Item1), 2100 - hmunit * (cntScafLenY), hmunit * w.Item1, hmunit * z.Item1);
                            }
                        });
                    }
                    else
                    {
                        gall.FillRectangle(myBrush2, 100, 2100 - hmunit * (cntScafLenY), 2000, hmunit * z.Item1);
                    }
                });
                //線を引く
                orientedScafNumY = 0;
                cntScafLenY = 0;
                scafLenInChrList.ForEach(z =>
                {
                    cntScafLenY += z.Item1;
                    if (z.Item2 == "+-")
                    {
                        orientedScafNumY++;
                        heatmap_chr_border.Select((x, i) => (x, i)).ToList().ForEach(w =>
                        {
                            if (w.x == orientedScafNumY - 1)
                            {
                                gall.DrawLine(pen, 100 + hmunit * (cntScafLenY - z.Item1), 100, 100 + hmunit * (cntScafLenY - z.Item1), 2100);
                                gall.DrawLine(pen, 100, 2100 - hmunit * (cntScafLenY - z.Item1), 2100, 2100 - hmunit * (cntScafLenY - z.Item1));
                                gall.DrawString((w.i + 1).ToString(), myfont, Brushes.Black, new Point(30, (int)(2100 - hmunit * (cntScafLenY - z.Item1)) - 30));
                                gall.DrawString((w.i + 1).ToString(), myfont, Brushes.Black, new Point((int)(100 + hmunit * (cntScafLenY - z.Item1)) + 10, 2130));
                            }
                        });
                    }
                });

                imageall.Save(opt_o + "_heatmap_phase_physical.png");
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
        public void drawBpCm(Graphics g, int fold, S_連鎖群中の向きと位置 scaf, long maxbp, double maxcm, int startx, int starty){
            var pen = new Pen(Color.Gray, fold);
            var brush = Brushes.LightGray;
            g.DrawLine(pen, new PointF(startx+75*fold, starty+(50 + 900 * ((float)scaf.pos1 / maxbp))*fold), new PointF(startx+175*fold, starty+(50 + 900 * (float)(scaf.cm1 / maxcm))*fold));
            g.DrawLine(pen, new PointF(startx+75*fold, starty+(50 + 900 * ((float)scaf.pos2 / maxbp))*fold), new PointF(startx+175*fold, starty+(50 + 900 * (float)(scaf.cm2 / maxcm))*fold));
            if (scaf.V_向き == "na")
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
        public S_連鎖群中の向きと位置 createScafPos(string chr, string order, long pos1, long pos2, double cm1, double cm2)
        {
            S_連鎖群中の向きと位置 tempscafpos = new S_連鎖群中の向きと位置();
            tempscafpos.V_連鎖群名 = chr;
            tempscafpos.V_向き = order;
            tempscafpos.pos1 = pos1;
            tempscafpos.pos2 = pos2;
            tempscafpos.cm1 = cm1;
            tempscafpos.cm2 = cm2;
            return tempscafpos;
        }
        public List<stranded_contig> searchFor(string query_contig_Name, string query_contig_StartEnd, List<S_コンティグ1StartOrEndと2StartOrEndの一致率> mainchains, List<string> flagpassed,
            Dictionary<string, string> I_コンティグ名から塩基配列への連想配列, Dictionary<string, SortedDictionary<int, Dictionary<string, int[]>>> I_コンティグ名から_家系IDから_StartEndから各個人のジェノタイプ情報,
            S_伸長したcontigたちの平均フェーズ I_平均フェーズ, double opt_shiftFromCenter, double opt_nonzerophasetotal)
        {
            List<stranded_contig> result = new List<stranded_contig>();
            S_コンティグ1StartOrEndと2StartOrEndの一致率[] tempchain = mainchains.Where(f => f.contig1==query_contig_Name && f.start_end_1==query_contig_StartEnd).ToArray();
            if(tempchain.Length==1){
                if(!flagpassed.Contains(tempchain[0].contig2)){
                    S_伸長したcontigたちの平均フェーズ I_temp_平均フェーズ = F_平均フェーズを計算(I_平均フェーズ, query_contig_Name, query_contig_StartEnd, I_コンティグ名から塩基配列への連想配列, I_コンティグ名から_家系IDから_StartEndから各個人のジェノタイプ情報);
                    if (I_temp_平均フェーズ.V_不一致度 > opt_shiftFromCenter || I_temp_平均フェーズ.V_有効個体数割合 < opt_nonzerophasetotal) { return null; }
                    flagpassed.Add(tempchain[0].contig2);
                    string temppos;
                    int temporder;
                    if(tempchain[0].start_end_2=="start"){
                        temppos="end";
                        temporder=1;
                    }else{
                        temppos="start";
                        temporder=-1;
                    }
                    stranded_contig tempstchr= new stranded_contig();
                    tempstchr.contig=tempchain[0].contig2;
                    tempstchr.orient=temporder;
                    tempstchr.match=tempchain[0].v_一致率;
                    result.Add(tempstchr);
                    List<stranded_contig> tempres = searchFor(tempchain[0].contig2, temppos, mainchains, flagpassed, I_コンティグ名から塩基配列への連想配列, I_コンティグ名から_家系IDから_StartEndから各個人のジェノタイプ情報, I_temp_平均フェーズ, opt_shiftFromCenter, opt_nonzerophasetotal);
                    if (tempres != null) { result.AddRange(tempres); }
                }
            }
            return result;
        }
        public List<stranded_contig> searchRev(string query_contig_Name, string query_contig_StartEnd, List<S_コンティグ1StartOrEndと2StartOrEndの一致率> mainchains, List<string> flagpassed,
            Dictionary<string, string> I_コンティグ名から塩基配列への連想配列, Dictionary<string, SortedDictionary<int, Dictionary<string, int[]>>> I_コンティグ名から_家系IDから_StartEndから各個人のジェノタイプ情報,
            S_伸長したcontigたちの平均フェーズ I_平均フェーズ, double opt_shiftFromCenter, double opt_nonzerophasetotal)
        {
            List<stranded_contig> result = new List<stranded_contig>();
            S_コンティグ1StartOrEndと2StartOrEndの一致率[] tempchain = mainchains.Where(f => f.contig1==query_contig_Name && f.start_end_1==query_contig_StartEnd).ToArray();
            if(tempchain.Length==1){
                if(!flagpassed.Contains(tempchain[0].contig2)){
                    S_伸長したcontigたちの平均フェーズ I_temp_平均フェーズ = F_平均フェーズを計算(I_平均フェーズ, query_contig_Name, query_contig_StartEnd, I_コンティグ名から塩基配列への連想配列, I_コンティグ名から_家系IDから_StartEndから各個人のジェノタイプ情報);
                    if (I_temp_平均フェーズ.V_不一致度 > opt_shiftFromCenter || I_temp_平均フェーズ.V_有効個体数割合 < opt_nonzerophasetotal) { return null; }
                    flagpassed.Add(tempchain[0].contig2);
                    string temppos;
                    int temporder;
                    if(tempchain[0].start_end_2=="start"){
                        temppos="end";
                        temporder=-1;
                    }else{
                        temppos="start";
                        temporder=1;
                    }
                    stranded_contig tempstchr= new stranded_contig();
                    tempstchr.contig=tempchain[0].contig2;
                    tempstchr.orient=temporder;
                    tempstchr.match=tempchain[0].v_一致率;
                    List<stranded_contig> tempres = searchRev(tempchain[0].contig2, temppos, mainchains, flagpassed, I_コンティグ名から塩基配列への連想配列, I_コンティグ名から_家系IDから_StartEndから各個人のジェノタイプ情報, I_temp_平均フェーズ, opt_shiftFromCenter, opt_nonzerophasetotal);
                    if (tempres != null) { result.AddRange(tempres); }
                    result.Add(tempstchr);
                }
            }
            return result;
        }

        public S_伸長したcontigたちの平均フェーズ F_平均フェーズを計算(S_伸長したcontigたちの平均フェーズ oldフェーズ, string query_contig_Name, string StartEnd, 
            Dictionary<string, string> I_コンティグ名から塩基配列への連想配列, Dictionary<string, SortedDictionary<int, Dictionary<string, int[]>>> I_コンティグ名から_家系IDから_StartEndから各個人のジェノタイプ情報)
        {
            int V_コンティグ長 = I_コンティグ名から塩基配列への連想配列[query_contig_Name].Length;
            S_伸長したcontigたちの平均フェーズ I_平均フェーズ計算結果 = new S_伸長したcontigたちの平均フェーズ();
            I_平均フェーズ計算結果.I_伸長した合計塩基数 = new SortedDictionary<int, long[]>();
            I_平均フェーズ計算結果.I_家系別のコンティグの平均フェーズ配列 = new SortedDictionary<int, double[]>();
            int total_有効個体数 = 0;
            int total_個体数 = 0;
            double total_不一致度合 = 0;
            foreach (KeyValuePair<int, double[]> I_計算前の家系ごとのコンティグの平均フェーズ配列 in oldフェーズ.I_家系別のコンティグの平均フェーズ配列)
            {
                double V_temp_フェーズ一致時の計算結果 = 0;
                double V_temp_フェーズ反転時の計算結果 = 0;
                long[] A_temp_コンティグ長合計 = new long[I_計算前の家系ごとのコンティグの平均フェーズ配列.Value.Length]; //0で初期化
                double[] A_temp_フェーズ一致時の平均フェーズ = new double[I_計算前の家系ごとのコンティグの平均フェーズ配列.Value.Length]; //0で初期化
                double[] A_temp_フェーズ反転時の平均フェーズ = new double[I_計算前の家系ごとのコンティグの平均フェーズ配列.Value.Length]; //0で初期化
                int tempN = 0;
                total_個体数 += I_計算前の家系ごとのコンティグの平均フェーズ配列.Value.Length;
                for (int i=0; i<I_計算前の家系ごとのコンティグの平均フェーズ配列.Value.Length; i++)
                {
                    if(I_コンティグ名から_家系IDから_StartEndから各個人のジェノタイプ情報[query_contig_Name][I_計算前の家系ごとのコンティグの平均フェーズ配列.Key][StartEnd][i] != -1)
                    {
                        tempN++;
                        if (oldフェーズ.I_伸長した合計塩基数[I_計算前の家系ごとのコンティグの平均フェーズ配列.Key][i] == 0) 
                        { 
                            I_計算前の家系ごとのコンティグの平均フェーズ配列.Value[i] = I_コンティグ名から_家系IDから_StartEndから各個人のジェノタイプ情報[query_contig_Name][I_計算前の家系ごとのコンティグの平均フェーズ配列.Key][StartEnd][i]; 
                        }
                        V_temp_フェーズ一致時の計算結果 += Math.Abs(I_コンティグ名から_家系IDから_StartEndから各個人のジェノタイプ情報[query_contig_Name][I_計算前の家系ごとのコンティグの平均フェーズ配列.Key][StartEnd][i]
                                                            - I_計算前の家系ごとのコンティグの平均フェーズ配列.Value[i]);
                        V_temp_フェーズ反転時の計算結果 += Math.Abs(1 - I_コンティグ名から_家系IDから_StartEndから各個人のジェノタイプ情報[query_contig_Name][I_計算前の家系ごとのコンティグの平均フェーズ配列.Key][StartEnd][i]
                                                            - I_計算前の家系ごとのコンティグの平均フェーズ配列.Value[i]);
                        A_temp_フェーズ一致時の平均フェーズ[i] = (I_計算前の家系ごとのコンティグの平均フェーズ配列.Value[i] * oldフェーズ.I_伸長した合計塩基数[I_計算前の家系ごとのコンティグの平均フェーズ配列.Key][i] +
                                                                    I_コンティグ名から_家系IDから_StartEndから各個人のジェノタイプ情報[query_contig_Name][I_計算前の家系ごとのコンティグの平均フェーズ配列.Key][StartEnd][i]* V_コンティグ長)/
                                                                    (oldフェーズ.I_伸長した合計塩基数[I_計算前の家系ごとのコンティグの平均フェーズ配列.Key][i] + V_コンティグ長);
                        A_temp_フェーズ反転時の平均フェーズ[i] = (I_計算前の家系ごとのコンティグの平均フェーズ配列.Value[i] * oldフェーズ.I_伸長した合計塩基数[I_計算前の家系ごとのコンティグの平均フェーズ配列.Key][i] +
                                            (1 - I_コンティグ名から_家系IDから_StartEndから各個人のジェノタイプ情報[query_contig_Name][I_計算前の家系ごとのコンティグの平均フェーズ配列.Key][StartEnd][i]) * V_コンティグ長) /
                                            (oldフェーズ.I_伸長した合計塩基数[I_計算前の家系ごとのコンティグの平均フェーズ配列.Key][i] + V_コンティグ長);
                        A_temp_コンティグ長合計[i] += V_コンティグ長;
                    }
                }
                I_平均フェーズ計算結果.I_伸長した合計塩基数[I_計算前の家系ごとのコンティグの平均フェーズ配列.Key] = A_temp_コンティグ長合計;
                total_有効個体数 += tempN;
                if (Math.Abs(V_temp_フェーズ一致時の計算結果) <= Math.Abs(V_temp_フェーズ反転時の計算結果))
                {
                    I_平均フェーズ計算結果.I_家系別のコンティグの平均フェーズ配列[I_計算前の家系ごとのコンティグの平均フェーズ配列.Key] = A_temp_フェーズ一致時の平均フェーズ;
                    //System.Console.WriteLine(Math.Abs(V_temp_フェーズ一致時の計算結果 / tempN));
                    total_不一致度合 += V_temp_フェーズ一致時の計算結果;
                }
                else
                {
                    I_平均フェーズ計算結果.I_家系別のコンティグの平均フェーズ配列[I_計算前の家系ごとのコンティグの平均フェーズ配列.Key] = A_temp_フェーズ反転時の平均フェーズ;
                    //System.Console.WriteLine(Math.Abs(V_temp_フェーズ反転時の計算結果 / tempN));
                    total_不一致度合 += V_temp_フェーズ反転時の計算結果;
                }
            }
            I_平均フェーズ計算結果.V_有効個体数割合 = total_有効個体数 / (double) total_個体数;
            if (total_有効個体数 > 0)
            {
                I_平均フェーズ計算結果.V_不一致度 = total_不一致度合 / total_有効個体数;
                //System.Console.WriteLine("不一致度:" + total_不一致度合 / total_有効個体数);
            }
            else
            {
                I_平均フェーズ計算結果.V_不一致度 = 2;
                //System.Console.WriteLine("total_有効個体数=0");
            }
            return I_平均フェーズ計算結果;
        }

        public static edge get_max_edge(string chr, List<string> croschr, Dictionary<string, SortedDictionary<int, Dictionary<string, int[]>>> datas, double opt_nonzerophase, int num_member){
            edge result = new edge();
            contig_match temp1 = calc_match_rate(chr, "start", "start", croschr, datas, opt_nonzerophase, num_member);
            contig_match temp2 = calc_match_rate(chr, "start", "end", croschr, datas, opt_nonzerophase, num_member);
            int i=0;
            double max=0;
            string maxchr="";
            string maxpos="";
            foreach(string tempchr in croschr){
                if(temp1.I_LIST_一致率[i]>max){
                    max=temp1.I_LIST_一致率[i];
                    maxchr=tempchr;
                    maxpos="start";
                }
                if(temp2.I_LIST_一致率[i]>max){
                    max=temp2.I_LIST_一致率[i];
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
        public static contig_match calc_match_rate(string chr, string xpos, string ypos, List<string> croschr, Dictionary<string, SortedDictionary<int, Dictionary<string, int[]>>> datas, double opt_nonzerophase, int num_member){
            contig_match resmatch = new contig_match();
            List<double> res = new List<double>();
            SortedDictionary<int, Dictionary<string, int[]>> tempx = datas[chr];
            foreach(string tempchr in croschr){
                res.Add(match_rate(tempx, datas[tempchr], xpos, ypos, opt_nonzerophase, num_member));
            }
            resmatch.contig=chr;
            resmatch.I_LIST_一致率=res;
            return resmatch;
        }
        public static double match_rate_for_ends(SortedDictionary<int, Dictionary<string, int[]>> x, double opt_nonzerophase, int num_member)
        {
            int num_ind = 0;
            int num_match = 0;
            foreach (KeyValuePair<int, Dictionary<string, int[]>> pair in x)
            {
                    int[] temp = get_dist(pair.Value["start"], pair.Value["end"]);
                    num_ind += temp[1];
                    num_match += temp[0];
            }
            if (num_ind < num_member * opt_nonzerophase)
            {
                return 0;
            }
            else
            {
                return num_match / (double)num_ind;
            }
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