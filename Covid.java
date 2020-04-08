import java.io.*;
import java.util.*;
import java.net.*;
import java.util.zip.*;


public class Covid {

   public static String humangenome[] = new String[23];
   public static String sarscov2 = "";
   public static String blnk = "";

   public static int nnt = 10;
   public static int nword = 1 << (nnt*2);
   public static Vector ntindex[] = new Vector[nword];

   public static void main(String[] args) {

      for(int i=1; i<=20; i++) {
         humangenome[i-1] = "Human/Homo_sapiens.GRCh38.dna.chromosome."+i+".fa";
      }
      humangenome[20] = "Human/Homo_sapiens.GRCh38.dna.chromosome.X.fa";
      humangenome[21] = "Human/Homo_sapiens.GRCh38.dna.chromosome.Y.fa";
      humangenome[22] = "Human/Homo_sapiens.GRCh38.dna.chromosome.MT.fa";

      String file1 = "virus/Sars-Cov-2.txt";
      String file2 = "virus/Sars-Cov-1.txt";
      String file3 = "virus/RaTG13.txt";
      String file4 = "virus/pangolin.txt";
      String file5 = "BLNK/homosapiensvariant1.txt";
      String file6 = "BLNK/macacamulattatvariantX3.txt";
      String file7 = "virus/adintovirus.txt";
      String file8 = "virus/HIV.txt";
      String file9 = "TB/strain1.txt";
      String file10 = "FURIN/FASTA.txt";
      String file11 = "genes/SERPINE1.txt";
      String file12 = "virus/python_nidovirus_KJ541759.1.txt";
      String file13 = "genes/chimp-OTUD6A.txt";
      String file14 = "genes/human-OTUD6A.txt";
      String file15 = "virus/avian_bronchitis_NC_001451.1.txt";
      String file16 = "virus/cellulophaga_phage_NC_021800.1.txt";
      String file17 = "virus/dengue_NC_001477.1.txt";
      String file18 = "virus/Hepatitis_B_NC_003977.2.txt";
      String file19 = "virus/HIV1_NC_001802.1.txt";
      String file20 = "virus/human_coronavirus_HKU1_NC_006577.2.txt";
      String file21 = "virus/human_herpesvirus_4_NC_009334.1.txt";
      String file22 = "virus/Influenza_A_NC_026437.1.txt";
      String file23 = "virus/mumps_virus_NC_002200.1.txt";
      String file24 = "virus/human_rhinovirus_3_NC_038312.1.txt";

      String trans1 = "Human/gencode.v33.pc_transcripts.fa";
      String trans2 = "Human/gencode.v33.transcripts.fa";
      String trans3 = "Human/GCF_000001405.39_GRCh38.p13_rna.fna";

      String bieti1 = "Monkey/GCF_001698545.1_ASM169854v1_genomic.gbff";
      String bieti2 = "Monkey/GCF_001698545.1_ASM169854v1_genomic.fna";
      String bieti3 = "Monkey/GCF_001698545.1_ASM169854v1_rna.fna";
      String macaca = "Monkey/GCF_003339765.1_Mmul_10_rna.fna";

      String battr = "Bat/GCF_004115265.1_mRhiFer1_v1.p_rna.fna";

      //readgenome(bieti2);


      String covi1 = "Sars-Cov-2/20191226-LR757995.1.txt";
      String covi2 = "Sars-Cov-2/20191230-MN996527.1.txt";
      String covi3 = "Sars-Cov-2/20200101-LR757996.1.txt";
      String covi4 = "Sars-Cov-2/20200108-MT093631.2.txt";
      String covi5 = "Sars-Cov-2/20200110-MN938384.1.txt";
      String covi6 = "Sars-Cov-2/20200120-MT226610.1.txt";
      String covi7 = "Sars-Cov-2/20200122-MT192773.1.txt";
      String covi8 = "Sars-Cov-2/20200125-MT192759.1.txt";
      String covi9 = "Sars-Cov-2/20200217-MT159705.1.txt";
      String covi10 = "Sars-Cov-2/20200224-MT184913.1.txt";
      String covi11 = "Sars-Cov-2/20200229-MT163718.1.txt";
      String covi12 = "Sars-Cov-2/20200305-MT188341.1.txt";
      String covi13 = "Sars-Cov-2/20200309-MT188339.1.txt";
      String covi14 = "Sars-Cov-2/20200311-MT192765.1.txt";
      String covi15 = "Sars-Cov-2/20200210-LC528232.1.txt";

      //comparei(file1,bieti2);

      /*
      for(int i=0; i<23; i++) {
         System.err.println("");
         System.out.println("");
         comparei(file3, humangenome[i]);
      }
      */

      //System.out.println("reverse,complement,junk,length,fileA,lengthA,startA,endA,textA,fileB,lengthB,startB,endB,textB");
      //System.err.println("reverse,complement,junk,length,fileA,lengthA,startA,endA,textA,fileB,lengthB,startB,endB,textB");

 
      //comparef("coronavirus","coronavirus/pangolin_MT084071.1.txt");

      // *****************************************************

      //System.out.println("reverse,complement,junk,length,idA,titelA,lengthA,startA,endA,textA,fileB,lengthB,startB,endB,textB");
      //System.err.println("reverse,complement,junk,length,idA,titleA,lengthA,startA,endA,textA,fileB,lengthB,startB,endB,textB");

      sequences("sequences/sequences.fasta", file1);

   }

   public static void sequences(String fileS, String fileB) {

      String B = strip(readgenome(fileB));

      makeindex(B);

      StringBuffer buffer = new StringBuffer("");

      try { 
         BufferedReader br = new BufferedReader(new FileReader(fileS));

         String line = "";
         String id = null;
         String title = "";
         boolean more = true;
         while(more) {
            line = br.readLine();
            if(line == null || line.startsWith(">")) {
               //System.err.println(line); 
               if(id != null) {
                  String A = buffer.toString();
                  System.err.println(id+","+title+","+A.length());
                  String fileA = id+","+title;
                  matchi(false, false, A, B, fileA, fileB);
                  matchi(false, true,  A, B, fileA, fileB);
                  matchi(true,  false, A, B, fileA, fileB);
                  matchi(true,  true,  A, B, fileA, fileB);
               }
               if(line != null) {
                  int i = line.indexOf("|");
                  id = line.substring(1,i).trim();
                  title = line.substring(i+1).trim();
                  i = title.indexOf("|");
                  if(i > -1) title = title.substring(0,i);
                  buffer = new StringBuffer("");
               } else {
                  more = false;
               }
            } else {
               buffer.append(line);
            }
         }

         br.close();
        

     } catch (IOException ioe) { 
            System.out.println(ioe);
     } 

   }

   public static void comparei(String fileA, String fileB) {
      String files[] = {fileA};
      comparei(files,fileB);
   }

   public static void comparef(String foldername, String fileB) {

      String B = strip(readgenome(fileB));

      makeindex(B);

      File folder = new File(foldername);
      String[] files = folder.list();

      for(int f=0; f<files.length; f++) {
         if(!(folder+"/"+files[f]).equals(fileB)) {
            String A = readgenome(folder+"/"+files[f]).toUpperCase();

            matchi(false, false, A, B, files[f], fileB);
            matchi(false, true,  A, B, files[f], fileB);
            matchi(true,  false, A, B, files[f], fileB);
            matchi(true,  true,  A, B, files[f], fileB);
         }
      }
   }

   public static void comparex(String foldername, String fileB) {

      File folder = new File(foldername);
      String[] files = folder.list();

      for(int istart=0; istart<39000000; istart+=5000000) {
         int iend = istart+5000000;

         String B = readgenome(fileB,istart,iend);

         makeindex(B);

         for(int f=0; f<files.length; f++) {
            String A = readgenome(folder+"/"+files[f]).toUpperCase();

            matchi(false, false, A, B, files[f], fileB);
            matchi(false, true,  A, B, files[f], fileB);
            matchi(true,  false, A, B, files[f], fileB);
            matchi(true,  true,  A, B, files[f], fileB);
         }
      }
   }

   public static void comparei(String files[], String fileB) {

      String B = strip(readgenome(fileB));

      makeindex(B);

      for(int f=0; f<files.length; f++) {
         String A = readgenome(files[f]).toUpperCase();

         matchi(false, false, A, B, files[f], fileB);
         matchi(false, true,  A, B, files[f], fileB);
         matchi(true,  false, A, B, files[f], fileB);
         matchi(true,  true,  A, B, files[f], fileB);
      }
   }

   public static void matchi(boolean reverse, boolean complement, String A, String B, String fileA, String fileB) {

     String T = strip(A);
     if(reverse) T = reverse(T);
     if(complement) T = complement(T);

     // must index B before using this version

     int n = 15;   // must be at least nnt

     int ilast = 0;
     int h = 0;
     for(int k=0; k<T.length()-n-h; k++) {
      String text = T.substring(k,k+n+h);
      if(text.indexOf("NNNNNNNNNN") == -1) {
        int i = indexOf(text,B);
        if(i == -1) {
           h--;
           if(h < 0) h = 0;
        } else {
           int j = i;
           String text1 = text;
           while (j > -1) {
              i = j;
              text = text1;
              if(k+n+h+1 > T.length()) {
                 j = -1;
              } else {
                 text1 = T.substring(k,k+n+h+1);
                 j = indexOf(text1,B);
              }
              h++;
           }
           if(k+n+h > ilast) {
              String textB = B.substring(i,i+n+h-1);
              String textT = T.substring(k,k+n+h-1);
              if(reverse) textT = reverse(textT);
              if(complement) textT = complement(textT);
              boolean junk = false;
              if(textB.indexOf("NNNNNNNNNNNNNN") == -1) {
                 if(textB.indexOf("AAAAAAA") > -1 || 
                    textB.indexOf("TTTTTTT") > -1 || 
                    textB.indexOf("NNNN") > -1 ||
                    textB.indexOf("ATATATATATAT") > -1) junk = true;
                 int istart = k+1;
                 int iend = k+n+h-1;
                 if(reverse) {
                    int ist = T.length()+1-iend;
                    iend = T.length()+1-istart;
                    istart = ist;
                 }
                 String line = reverse+","+complement+","+junk+","+(n+h-1)+","+fileA+","+A.length()+","+istart+","+iend+","+textT+","+fileB+","+B.length()+","+(i+1)+","+(i+n+h-1)+","+textB;
                 System.out.println(line);
                 System.err.println(line);
              }
           }
           ilast = k+n+h;
           //k += (h-1);
        }
      }
     }

   }

   public static String strip(String A) {
      A = A.replaceAll(" ","");
      A = A.replaceAll("0","");
      A = A.replaceAll("1","");
      A = A.replaceAll("2","");
      A = A.replaceAll("3","");
      A = A.replaceAll("4","");
      A = A.replaceAll("5","");
      A = A.replaceAll("6","");
      A = A.replaceAll("7","");
      A = A.replaceAll("8","");
      A = A.replaceAll("9","");
      A = A.toUpperCase();

      return A;
   }

   public static int indexOf(String text, String genome) { 
      // index of genome must be made before using this function

      int n = text.length();
      String word = text.substring(0,nnt).toUpperCase();
      int w = wordtoint(word);
      int m = ntindex[w].size();
      for(int k=0; k<m; k++) {
         Integer I = (Integer) ntindex[w].elementAt(k);
         int i = I.intValue();
         if(i+n <= genome.length() && genome.substring(i,i+n).equalsIgnoreCase(text)) return i;
      }

      return -1;
   }

   public static void makeindex(String text) {

      for(int w=0; w<nword; w++) ntindex[w] = new Vector();

      for(int i=0; i<text.length()-nnt+1; i++) {
         if(i%1000000 == 0) System.err.println(i+"/"+(text.length()-nnt+1));
         String word = text.substring(i,i+nnt).toUpperCase();
         int w = wordtoint(word);
         ntindex[w].add(new Integer(i));
      }

   }

   public static int wordtoint(String word) {
      int w = 0;
      for(int i=0; i<word.length(); i++) {
         String C = word.substring(i,i+1);
         int n = 0;
         if(C.equals("A")) n = 0;
         if(C.equals("C")) n = 1;
         if(C.equals("G")) n = 2;
         if(C.equals("T")) n = 3;
         w = w*4+n;
      }

      return w;
   }
         

   public static void compare(String file1, String file2) {

      System.err.println("");
      System.out.println("");

      System.err.println("comparing "+file1+" to "+file2);
      System.out.println("comparing "+file1+" to "+file2);

      String A = readgenome(file1).toUpperCase();
      //String B = readgenome(file2).toUpperCase();
      String B = readgenome(file2);
      //blnk = readgenome("BLNK/homosapiensvariant1.txt").toUpperCase();

      System.out.println("direct");
      System.err.println("direct");
      match(A, B);
      System.out.println("complement");
      System.err.println("complement");
      match(complement(A), B);
      System.out.println("reverse");
      System.err.println("reverse");
      match(reverse(A), B);
      System.out.println("reverse complement");
      System.err.println("reverse complement");
      match(reverse(complement(A)), B);
   }

   public static void match(String A, String B) {

     int n = 14;

     for(int k=0; k<A.length()-n; k++) {
        String text = A.substring(k,k+n);
        int i = B.indexOf(text);
        if(i > -1) {
           int h = 0;
           int j = i;
           String text1 = text;
           while (j > -1) {
              i = j;
              text = text1;
              if(k+n+h+1 > A.length()) {
                 j = -1;
              } else {
                 text1 = A.substring(k,k+n+h+1);
                 j = B.indexOf(text1);
              }
              h++;
           }
           System.out.println(k+"-"+(k+n+h-2)+" "+i+"-"+(i+n+h-2)+" "+(n+h-1)+" "+text);
           System.err.println(k+"-"+(k+n+h-2)+" "+i+"-"+(i+n+h-2)+" "+(n+h-1)+" "+text);
           k += (h-1);
        }
     }

   }

   public static String complement(String A) {

      StringBuffer B = new StringBuffer();

      for(int i=0; i<A.length(); i++) {
         String C = A.substring(i,i+1);
         if(C.equals("C")) B.append("G");
         if(C.equals("G")) B.append("C");
         if(C.equals("A")) B.append("T");
         if(C.equals("T")) B.append("A");
         if(C.equals("N")) B.append("N");
      }
      if(A.length() != B.length()) System.err.println("letters missed");

      return B.toString();

   }

   public static String reverse(String A) {

      StringBuffer B = new StringBuffer();

      for(int i=0; i<A.length(); i++) {
         int j = A.length()-1-i;
         String C = A.substring(j,j+1);
         B.append(C);
      }
      if(A.length() != B.length()) System.err.println("letters missed");

      return B.toString();

   }


   public static String readgenome(String filename) {

      System.err.println("reading "+filename);

      StringBuffer buffer = new StringBuffer("");

      try { 
         BufferedReader br = new BufferedReader(new FileReader(filename));

         int k=0;
         int g=0;
         String line = null;
         while((line = br.readLine() ) != null) {
            buffer.append(line);
            //if(k< 100) System.err.println(line);
            //if(line.startsWith("LOCUS")) System.err.println(line);
           k++;
           //for(int i=0; i<line.length(); i++) {
              //String C = line.substring(i,i+1);
              //if(C.equals("C") || C.equals("G") || C.equals("A") || C.equals("T")) g++;
           //}
         }

         System.err.println("k = "+k+" g = "+ g);

         br.close();
        

     } catch (IOException ioe) { 
            System.out.println(ioe);
     } 

     return buffer.toString();     


  }
 
    public static String readgenome(String filename, int istart, int iend) {

      System.err.println("reading "+filename);

      StringBuffer buffer = new StringBuffer("");

      try { 
         BufferedReader br = new BufferedReader(new FileReader(filename));

         int k=0;
         String line = null;
         while((line = br.readLine() ) != null) {
            if(k%1000000 == 0) System.err.println(istart+" "+k);
            if(k >= istart && k < iend) {
               buffer.append(line);
            }
            k++;
         }

         br.close();
        

     } catch (IOException ioe) { 
            System.out.println(ioe);
     } 

     return buffer.toString();     


  }

  //public static void String filter( String line ) {


  //}
 

}