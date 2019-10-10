using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace SELDLA{
    class ConvSNP{
        public struct newpos{
            public string newchr;
            public string order;
            public int pos;
        }
        public void run(int num_fam, string inputbreak, string inputmap, string opt_o, Dictionary<string, string> refseqs2){
            
            System.IO.StreamReader file = new System.IO.StreamReader(inputbreak);
            string line;
            int dup=0;
            string old="";
            Dictionary<string, List<int>> bpold2new = new Dictionary<string, List<int>>();
            while ((line = file.ReadLine()) != null){
                string[] values = line.Split("\t");
                if(values[0]!=old){
                    dup=2;
                    List<int> bps = new List<int>();
                    bps.Add(1);
                    bpold2new.Add(values[0],bps);
                }else{
                    dup++;
                }
                bpold2new[values[0]].Add(Int32.Parse(values[3]));
                old=values[0];
            }
            file.Close();

            System.IO.StreamReader file2 = new System.IO.StreamReader(inputmap);
            Dictionary<string, newpos> posold2new = new Dictionary<string, newpos>();
            while ((line = file2.ReadLine()) != null){
                if(!line.StartsWith("#")){
                    string[] values = line.Split("\t");
                    newpos temp = new newpos();
                    if(values[2]=="na"){
                        temp.newchr="linkage_scaffold_"+values[0];
                        temp.pos=Int32.Parse(values[3]);
                    }else{
                        temp.newchr="linkage_scaffold_"+values[0];
                        temp.pos=Int32.Parse(values[3]);
                    }
                    temp.order=values[2];
                    posold2new.Add(values[1], temp);
                }
            }
            file2.Close();

            for(int i=1;i<=num_fam;i++){
                System.IO.StreamReader file3 = new System.IO.StreamReader(opt_o + "_split_" + i + ".txt");
                System.IO.StreamWriter writer = new System.IO.StreamWriter(opt_o + "_split_" + i + "_newpos.txt");

                while ((line = file3.ReadLine()) != null)
                {
                    if (line.StartsWith("#"))
                    {
                        writer.WriteLine("#oldchr\t"+line);
                    }
                    else
                    {
                        string[] values = line.Split("\t");
                        string oldchr = values[0];
                        int oldpos = Int32.Parse(values[1]);
                        if (bpold2new.ContainsKey(oldchr))
                        {
                            int newind = findbreaked(oldpos, bpold2new[oldchr]);
                            oldpos = oldpos - bpold2new[oldchr][newind - 1] - 1;
                            oldchr = oldchr + "_" + newind;
                        }
                        if (posold2new.ContainsKey(oldchr))
                        {
                            newpos temp = posold2new[oldchr];
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
                            values[0] = oldchr;
                        }
                        if (posold2new.ContainsKey(oldchr) && posold2new[oldchr].order == "na")
                        {
                            writer.WriteLine("unordered_" + oldchr + "\t" + string.Join("\t", values));
                        }
                        else
                        {
                            writer.WriteLine(oldchr + "\t" + string.Join("\t", values));
                        }
                    }
                }
                file3.Close();
                writer.Close();
            }
        }

        public int findbreaked(int pos, List<int> list){
            int res=0;
            for(int i=0;i<list.Count;i++){
                if(pos>=list[i]){
                    res=i;
                }
            }
            return res+1;
        }
    }
}