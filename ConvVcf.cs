using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace SELDLA
{
    class ConvVcf
    {
        public struct newpos
        {
            public string newchr;
            public string order;
            public int pos;
        }
        public void run(string inputvcf, string inputbreak, string inputmap, string opt_o, Dictionary<string, string> refseqs2)
        {

            System.IO.StreamReader file = new System.IO.StreamReader(inputbreak);
            string line;
            string old = "";
            Dictionary<string, List<int>> bpold2new = new Dictionary<string, List<int>>();
            while ((line = file.ReadLine()) != null)
            {
                string[] values = line.Split("\t");
                string temp_breaked_chr=values[0];
                int temp_breaked_pos=Int32.Parse(values[3]);
                if (temp_breaked_chr != old)
                {
                    List<int> breaked_position = new List<int>();
                    breaked_position.Add(0);
                    bpold2new.Add(temp_breaked_chr, breaked_position);
                }
                bpold2new[temp_breaked_chr].Add(temp_breaked_pos);
                old = temp_breaked_chr;
            }
            file.Close();

            System.IO.StreamReader file2 = new System.IO.StreamReader(inputmap);
            Dictionary<string, newpos> posold2new = new Dictionary<string, newpos>();
            Dictionary<string, newpos> posold2newincna = new Dictionary<string, newpos>();
            while ((line = file2.ReadLine()) != null)
            {
                if (!line.StartsWith("#"))
                {
                    string[] values = line.Split("\t");
                    newpos temp = new newpos();
                    if (values[2] == "na")
                    {
                        temp.newchr = "linkage_scaffold_" + values[0] + "_unoriented_pos" + values[3] +"_old_" + values[1];
                        temp.pos = 0;
                    }
                    else
                    {
                        temp.newchr = "linkage_scaffold_" + values[0];
                        temp.pos = Int32.Parse(values[3]);
                    }
                    temp.order = values[2];
                    posold2new.Add(values[1], temp);

                    //NAもchrに含める場合
                    newpos tempincna = new newpos();
                    tempincna.newchr = "linkage_scaffold_" + values[0];
                    tempincna.pos = Int32.Parse(values[3]);
                    tempincna.order = values[2];
                    posold2newincna.Add(values[1], tempincna);
                }
            }
            file2.Close();

            System.IO.StreamReader file3 = new System.IO.StreamReader(inputvcf);
            System.IO.StreamWriter writer = new System.IO.StreamWriter(opt_o + "_newpos.vcf");
            System.IO.StreamWriter writer2 = new System.IO.StreamWriter(opt_o + "_newpos_include_unoriented_in_chr.vcf");
            while ((line = file3.ReadLine()) != null)
            {
                if (line.StartsWith("#"))
                {
                    writer.WriteLine(line);
                    writer2.WriteLine(line);
                }
                else
                {
                    string[] values = line.Split("\t");
                    string oldchr = values[0];
                    int oldpos = Int32.Parse(values[1]);
                    if (bpold2new.ContainsKey(oldchr))
                    {
                        int newind = findbreaked(oldpos, bpold2new[oldchr]);
                        oldpos = oldpos - bpold2new[oldchr][newind - 1];
                        oldchr = oldchr + "_" + newind;
                    }
                    if (posold2new.ContainsKey(oldchr))
                    {
                        newpos temp = posold2new[oldchr];
                        values[0] = temp.newchr;
                        if (temp.order == "-")
                        {
                            //-側は変異を相補鎖に変換 ここは下のNAもchrに含める場合にも影響
                            string tempref = values[3];
                            if (tempref.Length > 1)
                            {
                                oldpos = oldpos + tempref.Length - 1;
                            }
                            values[3] = Chain.getRevComp(values[3]);
                            string[] tempalt = values[4].Split(",");
                            values[4] = Chain.getRevComp(tempalt[0]);
                            for (int i = 1; i < tempalt.Length; i++)
                            {
                                values[4] = values[4] + "," + Chain.getRevComp(tempalt[i]);
                            }
                            //以下はNAを含めない場合のみに影響
                            values[1] = (temp.pos + refseqs2[oldchr].Length - oldpos + 1).ToString();
                        }
                        else
                        {
                            values[1] = (temp.pos + oldpos).ToString();
                        }
                    }
                    else
                    {
                        values[0] = "old_"+oldchr;
                        values[1] = oldpos.ToString();
                    }
                    writer.WriteLine(string.Join("\t", values));

                    //NAもchrに含める場合
                    if (posold2newincna.ContainsKey(oldchr))
                    {
                        newpos temp = posold2newincna[oldchr];
                        values[0] = temp.newchr;
                        if (temp.order == "-")
                        {
                            values[1] = (temp.pos + refseqs2[oldchr].Length - oldpos + 1).ToString();
                        }
                        else
                        {
                            values[1] = (temp.pos + oldpos).ToString();
                        }
                    }
                    else
                    {
                        values[0] = "old_"+oldchr;
                        values[1] = oldpos.ToString();
                    }
                    writer2.WriteLine(string.Join("\t", values));
                }
            }
            file3.Close();
            writer.Close();
            writer2.Close();
        }

        public int findbreaked(int pos, List<int> list)
        {
            int res = 0;
            for (int i = 0; i < list.Count; i++)
            {
                if (pos > list[i])
                {
                    res = i;
                }
            }
            return res + 1;
        }
    }
}