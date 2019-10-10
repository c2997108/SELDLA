using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace SELDLA{
        public class EndData : System.IComparable
    {
        public int spos;
        public int epos;
        public int[] sdata;
        public int[] edata;
        public void setStartPos(int x)
        {
            spos = x;
        }
        public void setEndPos(int x)
        {
            epos = x;
        }
        public void setStartData(int[] x)
        {
            sdata = x;
        }
        public void setEndData(int[] x)
        {
            edata = x;
        }

        //自分自身がobjより小さいときはマイナスの数、大きいときはプラスの数、
        //同じときは0を返す
        public int CompareTo(object obj)
        {
            //nullより大きい
            if (obj == null)
            {
                return 1;
            }

            //違う型とは比較できない
            if (this.GetType() != obj.GetType())
            {
                throw new ArgumentException("別の型とは比較できません。", "obj");
            }
            //このクラスが継承されることが無い（構造体など）ならば、次のようにできる
            //if (!(other is TestClass)) { }

            //Priceを比較する
            return this.spos.CompareTo(((EndData)obj).spos);
            //または、次のようにもできる
            //return this.Price - ((Product)other).Price;
        }
    }
}