using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;

namespace SELDLA
{

    class Ld2Ph : Snp2Ld
    {
        public double opt_clust_match_rate = 0.8;
        public int opt_clust_size = 2;
        public double opt_split_match_rate = 0.7;
        public StreamWriter write_ph;
        public StreamWriter write_ldimp;
        public struct datapos
        {
            public List<int[]> mydata;
            public List<int> mypos;
            public List<int> mystart;
            public List<int> myend;
        }
        public struct adatapos
        {
            public int[] mydata;
            public int mypos;
            public int mystart;
            public int myend;
        }
        public struct enddatastruct
        {
            public adatapos mystartdata;
            public adatapos myenddata;
        }
        public void run2(string input, double opt_balance, double opt_cm, int opt_cs, double opt_sm, bool flagsplit, int opt_ldnum, bool maxLdClusterOnly, double rateOfNotNA, int opt_ldseqnum)
        {

            opt_clust_match_rate = opt_cm;
            opt_split_match_rate = opt_sm;
            opt_clust_size = opt_cs;
            int counter = 0;
            string line;
            List<int> start_data = new List<int>();
            List<int> end_data = new List<int>();

            // Read the file and display it line by line.  
            System.IO.StreamReader file = new System.IO.StreamReader(input);
            while ((line = file.ReadLine()) != null)
            {
                //System.Console.WriteLine(line);
                if (counter == 0)
                {
                    string[] header = line.Split("\t");
                    //Console.WriteLine(string.Join(", ", header));
                    writer = new StreamWriter(input + ".break");
                    writer.WriteLine("breakpoint\tstart\tend");
                    write_ph = new StreamWriter(input + ".ph");
                    write_ph.Write("#chr\tnum_of_split/quality\tLD_start_pos");
                    line.Split("\t").Skip(3).ToList().ForEach(x=>write_ph.Write("\t"+x));
                    write_ph.WriteLine();
                    write_ldimp = new StreamWriter(input + "imp");
                    write_ldimp.WriteLine(line);
                }
                else
                {
                    string[] values = line.Split("\t");
                    int n0 = 0;
                    int n1 = 0;
                    int num_people = values.Length - 3;
                    //string[] chrpos=values[0].Split(":");
                    if (oldchr != values[1] && oldchr != "" && num > 0)
                    {
                        Console.WriteLine("calculating the phases of "+oldchr);
                        clustsearch2(data, pos, oldchr, flagsplit, opt_ldnum, start_data, end_data, maxLdClusterOnly, rateOfNotNA, opt_ldseqnum);
                        initialize();
                        start_data = new List<int>();
                        end_data = new List<int>();
                    }
                    oldchr = values[1];
                    for (int i = 3; i < values.Length; i++)
                    {
                        if (values[i] == "0") { n0++; } else if (values[i] == "1") { n1++; }
                    }
                    if (n0 > num_people * opt_balance && n1 > num_people * opt_balance)
                    {
                        num++;
                        int[] tempdata = new int[num_people];
                        for (int i = 3; i < values.Length; i++)
                        {
                            tempdata[i - 3] = Int32.Parse(values[i]);
                        }
                        data.Add(tempdata);
                        pos.Add(Int32.Parse(values[2]));

                        string[] val2 = values[0].Split("#");
                        start_data.Add(Int32.Parse(val2[1]));
                        end_data.Add(Int32.Parse(val2[2]));
                    }
                }
                counter++;
            }
            //ファイルの最後のscaffoldのクラスターを検索
            if (num > 0)
            {
                Console.WriteLine("calculating the phases of " + oldchr);
                clustsearch2(data, pos, oldchr, flagsplit, opt_ldnum, start_data, end_data, maxLdClusterOnly, rateOfNotNA, opt_ldseqnum);
            }

            file.Close();
            writer.Close();
            write_ph.Close();
            write_ldimp.Close();
            //System.Console.WriteLine("There were {0} lines.", counter);
            // Suspend the screen.  
            // System.Console.ReadLine();
        }

        public void clustsearch2(List<int[]> cdata, List<int> cpos, string chr, bool flagsplit, int opt_ldnum,
         List<int> start_data, List<int> end_data, bool maxLdClusterOnly, double rateOfNotNA, int opt_ldseqnum)
        {
            //LD間の総当たり距離を計算し、クラスター化の準備を行う
            int numdata = cpos.Count;
            double[,] distmatrix = new double[numdata, numdata];
            int[,] corder = new int[numdata, numdata];
            bool[,] flagclust = new bool[numdata, numdata];
            bool[] skipflag = new bool[numdata];
            bool[,] flagclustat1 = new bool[numdata, numdata];
            for (int i = 0; i < numdata; i++)
            {
                for (int j = 0; j < numdata; j++)
                {
                    // i==jの場合は0で初期化されている
                    if (i > j)
                    {
                        distance tempdist = get_dist(cdata[i], cdata[j], rateOfNotNA);
                        distmatrix[i, j] = tempdist.match_rate;
                        distmatrix[j, i] = tempdist.match_rate;
                        corder[i, j] = tempdist.order;
                        corder[j, i] = tempdist.order;
                        if (tempdist.match_rate >= opt_clust_match_rate)
                        {
                            flagclust[i, j] = true;
                            flagclust[j, i] = true;
                            if (tempdist.match_rate == 1)
                            {
                                flagclustat1[i, j] = true;
                                flagclustat1[j, i] = true;
                            }
                        }
                    }
                }
            }
            //閾値以下のクラスター数のLDは無視するフラグを付ける
            for (int i = 0; i < numdata; i++)
            {
                if (getSameLdNum(i, flagclustat1, numdata) < opt_ldnum)
                {
                    skipflag[i] = true;
                }
            }

            //最大クラスターのみ使用する場合、最大クラスター以外のLDにskipフラグを付ける
            if (maxLdClusterOnly && !flagsplit)
            {
                bool[] forMaxClusterSearchFlag = skipflag.Clone() as bool[];
                int tempMaxClusterSize = 0;
                List<int> tempMaxClusterId = new List<int>();
                //最大クラスターの探索
                for (int i = 0; i < numdata; i++)
                {
                    if (!forMaxClusterSearchFlag[i])
                    {
                        forMaxClusterSearchFlag[i] = true;
                        List<int> tempIdArrForSearchMaxCluster = new List<int>();
                        tempIdArrForSearchMaxCluster.Add(i);
                        tempIdArrForSearchMaxCluster.AddRange(searchclust(flagclust, forMaxClusterSearchFlag, i));
                        if (tempIdArrForSearchMaxCluster.Count > tempMaxClusterSize)
                        {
                            tempMaxClusterSize = tempIdArrForSearchMaxCluster.Count;
                            tempMaxClusterId = tempIdArrForSearchMaxCluster;
                        }
                    }
                }
                //スキップフラグを付ける
                for (int i = 0; i < numdata; i++)
                {
                    if (!tempMaxClusterId.Contains(i))
                    {
                        skipflag[i] = true;
                    }
                }
            }

            //フェーズ決定のためのインピュテーション、ブレークポイント探しなど
            List<enddatastruct> listends = new List<enddatastruct>();
            for (int i = 0; i < numdata; i++)
            {
                if (!skipflag[i])
                {
                    //scaffold内のクラスターを探索
                    skipflag[i] = true;
                    List<int> temp = new List<int>();
                    temp.Add(i);
                    temp.AddRange(searchclust(flagclust, skipflag, i));
                    temp.Sort();

                    if (temp.Count >= opt_clust_size)
                    {
                        List<int[]> templist = new List<int[]>();
                        List<int> temppos = new List<int>();
                        List<int> tempstart = new List<int>();
                        List<int> tempend = new List<int>();
                        foreach (int j in temp)
                        {
                            templist.Add(cdata[j]);
                            temppos.Add(cpos[j]);
                            tempstart.Add(start_data[j]);
                            tempend.Add(end_data[j]);
                        }
                        datapos tempdatapos = get_imputed(templist, temppos, tempstart, tempend, rateOfNotNA, opt_ldseqnum);
                        foreach(var item in tempdatapos.mypos.Select((Value, Index) => new { Value, Index }))
                        {
                            write_ldimp.Write("info\t" + chr + "\t" + item.Value);
                            foreach(var item2 in tempdatapos.mydata[item.Index])
                            {
                                write_ldimp.Write("\t" + item2);
                            }
                            write_ldimp.WriteLine();
                        }

                        //フェーズの両端を探す
                        adatapos startdat = new adatapos();
                        adatapos enddat = new adatapos();
                        for (int j = 0; j < tempdatapos.mypos.Count; j++)
                        {
                            //フェーズの最初
                            if (tempdatapos.mypos[j] == tempdatapos.mypos.Min())
                            {
                                //Console.WriteLine("min: "+tempdatapos.mypos[j]);
                                startdat.mydata = tempdatapos.mydata[j];
                                startdat.mypos = tempdatapos.mypos[j];
                                startdat.mystart = tempdatapos.mystart[j];
                                startdat.myend = tempdatapos.myend[j];
                                break;
                            }
                        }
                        for (int j = 0; j < tempdatapos.mypos.Count; j++)
                        {
                            //フェーズの最後
                            if (tempdatapos.mypos[j] == tempdatapos.mypos.Max())
                            {
                                //Console.WriteLine("max: "+tempdatapos.mypos[j]);
                                enddat.mydata = tempdatapos.mydata[j];
                                enddat.mypos = tempdatapos.mypos[j];
                                enddat.mystart = tempdatapos.mystart[j];
                                enddat.myend = tempdatapos.myend[j];
                                break;
                            }
                        }
                        enddatastruct ends = new enddatastruct();
                        ends.mystartdata = startdat;
                        ends.myenddata = enddat;
                        //クラスターごとに両端のフェーズを既知のクラスター両端と比較し、新しい両端のほうが長ければ統合する
                        foreach (enddatastruct tempends in listends)
                        {
                            if (tempends.mystartdata.mypos < ends.mystartdata.mypos
                            && tempends.myenddata.mypos > ends.mystartdata.mypos)
                            {
                                ends.mystartdata = tempends.mystartdata;
                            }
                            if (tempends.mystartdata.mypos < ends.myenddata.mypos
                            && tempends.myenddata.mypos > ends.myenddata.mypos)
                            {
                                ends.myenddata = tempends.myenddata;
                            }
                        }
                        for (int j = 0; j < listends.Count; j++)
                        {
                            if (listends[j].mystartdata.mypos >= ends.mystartdata.mypos
                            && listends[j].myenddata.mypos <= ends.myenddata.mypos)
                            {
                                listends.RemoveAt(j);
                            }
                        }
                        listends.Add(ends);

                        //if(temp.Count>=maxmember.Count){
                        //    maxmember = temp;
                        //}
                    }
                }
            }
            if (listends.Count >= 1)
            {
                List<EndData> listends2 = new List<EndData>();
                foreach (enddatastruct enddata in listends)
                {
                    //Console.WriteLine("start: "+enddata.mystartdata.mypos);
                    //Console.WriteLine("end: "+enddata.myenddata.mypos);
                    EndData tempend = new EndData();
                    tempend.setStartPos(enddata.mystartdata.mystart);
                    tempend.setStartData(enddata.mystartdata.mydata);
                    tempend.setEndPos(enddata.myenddata.myend);
                    tempend.setEndData(enddata.myenddata.mydata);
                    listends2.Add(tempend);
                }
                listends2.Sort();
                if (listends2.Count > 1)
                {
                    //入力のscaffoldを分割する
                    if (flagsplit)
                    {
                        int tempstartpos = listends2[0].spos;
                        int[] tempstartdata = listends2[0].sdata;
                        int numsplit = 0;
                        for (int i = 1; i < listends2.Count; i++)
                        {
                            //double tempdist=get_dist(listends2[i - 1].edata, listends2[i].sdata).match_rate;
                            if (get_dist(listends2[i - 1].edata, listends2[i].sdata, rateOfNotNA).match_rate < opt_split_match_rate)
                            {
                                //Console.WriteLine(chr + "\t" + listends2[i - 1].epos + "\t" + listends2[i].spos);
                                writer.WriteLine(chr + "\t" + listends2[i - 1].epos + "\t" + listends2[i].spos);
                                numsplit++;
                                write_ph.WriteLine(chr + "\t" + numsplit + "\t" + tempstartpos + "\t" + arr2str(tempstartdata));
                                write_ph.WriteLine(chr + "\t" + numsplit + "\t" + listends2[i - 1].epos + "\t" + arr2str(listends2[i - 1].edata));
                                tempstartpos = listends2[i].spos;
                                tempstartdata = listends2[i].sdata;
                            }
                        }
                        numsplit++;
                        write_ph.WriteLine(chr + "\t" + numsplit + "\t" + tempstartpos + "\t" + arr2str(tempstartdata));
                        write_ph.WriteLine(chr + "\t" + numsplit + "\t" + listends2[listends2.Count - 1].epos + "\t" + arr2str(listends2[listends2.Count - 1].edata));
                    }
                    else
                    {
                        List<EndData> templistends2 = listends2.OrderByDescending(x => x.epos - x.spos).ToList();
                        write_ph.WriteLine(chr + "\t" + listends2.Count + "\t" + templistends2[0].spos + "\t" + arr2str(templistends2[0].sdata));
                        write_ph.WriteLine(chr + "\t" + listends2.Count + "\t" + templistends2[0].epos + "\t" + arr2str(templistends2[0].edata));
                    }
                }
                else
                {
                    write_ph.WriteLine(chr + "\t0\t" + listends2[0].spos + "\t" + arr2str(listends2[0].sdata));
                    write_ph.WriteLine(chr + "\t0\t" + listends2[0].epos + "\t" + arr2str(listends2[0].edata));
                }
            }
            else
            {
                int[] tempdata = get_consensus(cdata, rateOfNotNA);
                write_ph.WriteLine(chr + "\tlowqual\t" + cpos[0] + "\t" + arr2str(tempdata));
                write_ph.WriteLine(chr + "\tlowqual\t" + cpos[cpos.Count - 1] + "\t" + arr2str(tempdata));
            }
        }
        /// <summary>
        /// 同じクラスターに属するLDの数を数える
        /// </summary>
        /// <param name="key"></param>
        /// <param name="flagclustat1"></param>
        /// <param name="numdata"></param>
        /// <returns></returns>
        public int getSameLdNum(int key, bool[,] flagclustat1, int numdata)
        {
            int res = 1;
            for (int i = 0; i < numdata; i++)
            {
                if (flagclustat1[key, i])
                {
                    res++;
                }
            }
            return res;
        }
        public datapos get_imputed(List<int[]> cdata, List<int> cpos, List<int> start_data, List<int> end_data, double rateOfNotNA, int opt_ldseqnum)
        {
            int[] tempcons = new int[cdata[0].Length];
            for (int j = 0; j < tempcons.Length; j++) { tempcons[j] = -1; }

            datapos phasedatapos = get_phased(cdata, cpos, start_data, end_data, rateOfNotNA);
            List<int[]> phaseddata = phasedatapos.mydata;

            for (int j = 0; j < tempcons.Length; j++)
            {
                int[] tempPhValues = new int[phaseddata.Count];
                List<int> tempSamePh = new List<int>();
                int tempCurPh = -1;
                int[] tempNumPh = new int[2]; //フェーズは0 or 1の2値
                //閾値以上の連続する同じLDが存在するか調べる
                for (int i = 0; i < phaseddata.Count; i++) { tempPhValues[i] = -1; }
                for (int i = 0; i < phaseddata.Count; i++) 
                {
                    if (phaseddata[i][j] != -1)
                    {
                        tempNumPh[phaseddata[i][j]]++;
                        if (phaseddata[i][j] == tempCurPh)
                        {
                            tempSamePh.Add(i);
                        }
                        else
                        {
                            if (tempSamePh.Count >= opt_ldseqnum)
                            {
                                tempSamePh.ForEach(x => { tempPhValues[x] = tempCurPh; });
                            }
                            tempSamePh.Clear();
                            tempSamePh.Add(i);
                            tempCurPh = phaseddata[i][j];
                        }
                    }
                }
                if (tempSamePh.Count >= opt_ldseqnum)
                {
                    tempSamePh.ForEach(x => { tempPhValues[x] = tempCurPh; });
                }

                //インピュテーションを行う
                int startImp = -1;
                int endImp = -1;
                int cntOfPh = 0;
                for (int i = 0; i < phaseddata.Count; i++)
                {
                    if (startImp == -1 && tempPhValues[i] != -1)
                    {
                        startImp = i;
                        cntOfPh++;
                    }
                    else if (startImp > -1 && tempPhValues[i] != -1)
                    {
                        endImp = i;
                        if (endImp - startImp > 1)
                        {
                            for(int k = startImp + 1; k < endImp; k++)
                            {
                                if (phasedatapos.mypos[k] - phasedatapos.mypos[startImp] <= phasedatapos.mypos[endImp] - phasedatapos.mypos[k])
                                {
                                    tempPhValues[k] = tempPhValues[startImp];
                                }
                                else
                                {
                                    tempPhValues[k] = tempPhValues[endImp];
                                }
                            }
                        }

                        startImp = i;
                        endImp = -1;
                    }
                }
                //閾値以上のフェーズがない場合は、最も多いフェーズで埋めておく
                if (cntOfPh == 0) { 
                    for (int i = 0; i < phaseddata.Count; i++) {
                        if (tempNumPh[0] > tempNumPh[1]) 
                        {
                            tempPhValues[i] = 0;
                        }else if (tempNumPh[1] > tempNumPh[0])
                        {
                            tempPhValues[i] = 1;
                        }
                        else
                        {
                            tempPhValues[i] = -1;
                        }
                    }
                }
                else
                {
                    //末端のフェーズを補完しておく
                    int tempPhase = -1;
                    for (int i = 0; i < phaseddata.Count; i++)
                    {
                        if ( tempPhValues[i] != -1)
                        {
                            tempPhase = tempPhValues[i];
                        }
                        else
                        {
                            tempPhValues[i] = tempPhase;
                        }
                    }
                    tempPhase = -1;
                    for (int i = phaseddata.Count-1; i >=0; i--)
                    {
                        if (tempPhValues[i] != -1)
                        {
                            tempPhase = tempPhValues[i];
                        }
                        else
                        {
                            tempPhValues[i] = tempPhase;
                        }
                    }
                }

                for (int i = 0; i < phaseddata.Count; i++)
                {
                    phaseddata[i][j] = tempPhValues[i];
                }
            }


            //ver2.0.9までのインピュテーション
            //for (int i = 0; i < phaseddata.Count; i++)
            //{
            //    for (int j = 0; j < tempcons.Length; j++)
            //    {
            //        if (phaseddata[i][j] != -1)
            //        {
            //            tempcons[j] = phaseddata[i][j];
            //        }
            //        else if (phaseddata[i][j] == -1 && tempcons[j] != -1)
            //        {
            //            phaseddata[i][j] = tempcons[j];
            //        }
            //    }
            //}
            //for (int i = phaseddata.Count - 1; i >= 0; i--)
            //{
            //    for (int j = 0; j < tempcons.Length; j++)
            //    {
            //        if (phaseddata[i][j] != -1)
            //        {
            //            tempcons[j] = phaseddata[i][j];
            //        }
            //        else if (phaseddata[i][j] == -1 && tempcons[j] != -1)
            //        {
            //            phaseddata[i][j] = tempcons[j];
            //        }
            //    }
            //}
            phasedatapos.mydata = phaseddata;
            return phasedatapos;
        }
        /// <summary>
        /// 与えられたLDのフェーズをそろえる
        /// </summary>
        /// <param name="cdata"></param>
        /// <param name="cpos"></param>
        /// <param name="start_data"></param>
        /// <param name="end_data"></param>
        /// <param name="rateOfNotNA"></param>
        /// <returns></returns>
        public datapos get_phased(List<int[]> cdata, List<int> cpos, List<int> start_data, List<int> end_data, double rateOfNotNA)
        {
            datapos resultdatapos = new datapos();
            List<int[]> result = new List<int[]>();
            List<int> resultpos = new List<int>();
            List<int> resultstart = new List<int>();
            List<int> resultend = new List<int>();

            int[] keyphase = new int[cdata[0].Length];
            Array.Copy(cdata[0], keyphase, keyphase.Length);
            result.Add(keyphase);
            resultpos.Add(cpos[0]);
            resultstart.Add(start_data[0]);
            resultend.Add(end_data[0]);
            for (int i = 1; i < cdata.Count; i++)
            {
                int[] temp = new int[keyphase.Length];
                distance tempdist = get_dist(keyphase, cdata[i], rateOfNotNA);
                if (tempdist.order == 1)
                {
                    Array.Copy(cdata[i], temp, keyphase.Length);
                    result.Add(temp);
                    resultpos.Add(cpos[i]);
                    resultstart.Add(start_data[i]);
                    resultend.Add(end_data[i]);
                }
                else if (tempdist.order == -1)
                {
                    for (int j = 0; j < keyphase.Length; j++)
                    {
                        if (cdata[i][j] == 0)
                        {
                            temp[j] = 1;
                        }
                        else if (cdata[i][j] == 1)
                        {
                            temp[j] = 0;
                        }
                        else
                        {
                            temp[j] = -1;
                        }
                    }
                    result.Add(temp);
                    resultpos.Add(cpos[i]);
                    resultstart.Add(start_data[i]);
                    resultend.Add(end_data[i]);
                }
            }
            resultdatapos.mydata = result;
            resultdatapos.mypos = resultpos;
            resultdatapos.mystart = resultstart;
            resultdatapos.myend = resultend;
            return resultdatapos;
        }

    }


}