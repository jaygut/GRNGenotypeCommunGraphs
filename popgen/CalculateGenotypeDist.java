
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * 
 * @author jaysongutierrez
 *
 *A class holding methods to calculate the average Nei's genetic distance (Nei, Tajima & Tateno 1983) between two diploid ARO genotypes
 */


public class CalculateGenotypeDist {
    
    public static List<String> readStringsFromFile(String genomeFileName) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(genomeFileName));
        String line = reader.readLine();
        List<String> lines = new ArrayList<String>();
        while (line != null) {
            lines.add(line);
            line = reader.readLine();
        }
        reader.close();
        return lines;
    }
    
    /*
     * Method to fetch diploid genotypes (e.g. functional sequences) as an array of strings, e.g. {["AT",....,"GC], ["AT",....,"GC], ["AA",....,"CC]}
     */
    public static List<String[]> getGenotypedSitesList(String inputFileName) throws IOException{
        
        List<String[]> genotypedSitesList = new ArrayList<String[]>();
        
        List<String> genotypes = readStringsFromFile(inputFileName);
            
        for (String g : genotypes){
            
            String[] genotypedSites = g.split(",");
            
            genotypedSitesList.add(genotypedSites);
            
        }
        
        return genotypedSitesList;
    }
    
    /*
     * Method to calculate the Nei's genetic distance (Nei, Tajima & Tateno 1983) between two genotypes at a particular site.
     * See Nei's distance definition in: Flexible methods for estimating genetic distances from single nucleotide polymorphisms (2015).
     */
    public static double perSiteNeiDistance(String allelesSite0, String allelesSite1){
        
        double total = 0.0;
        
        char[] nucls = {'A','C','T','G'};
        
        //This is the number of alleles per site, e.g. 2 alleles for a diploid at any site: AA
        int numAlleles = allelesSite0.length();
        
        for (char n: nucls){
            
            //Assess frequency of a given nucleotide at the site for genotype 0
            long nuclCount0 = allelesSite0.chars().filter(ch -> ch == n).count();
            double freq0 = (double) nuclCount0/numAlleles;
            
            //Assess frequency of a given nucleotide at the site for genotype 1
            long nuclCount1 = allelesSite1.chars().filter(ch -> ch == n).count();
            double freq1 = (double) nuclCount1/numAlleles;
            
            //Sum up the sqrt of f0*f1 for each nucleotide
            double sqrtFreqProduct = Math.sqrt(freq0*freq1);
            
            total += sqrtFreqProduct;
        }
        
        //calculated of Nei's genetic distance as described in: Flexible methods for estimating genetic distances from single nucleotide polymorphisms (2015).
        return 1-total; 
    }
    
    /*
     * This method calculates the average Nei's genetic distance over all DNA positions/sites between two diploid genotypes
     */
    public static double calcNeiGeneticDistance(String[] genotypedSitesList0, String[] genotypedSitesList1){
        
        double numSites = genotypedSitesList0.length;
        double distNei = 0.0;
        
        for(int i=0; i<numSites; i++){
            
            distNei += perSiteNeiDistance(genotypedSitesList0[i], genotypedSitesList1[i]);
        }
        
        double avgNeiDistance = distNei/numSites;
        
        //return Nei distance averaged across all DNA positions/sites
        return avgNeiDistance;
    }
    
    /*
     * Method to dump to a file a matrix, e.g. the distance matrix
     */
    static void writeMatrix(String filename, double[][] matrix) {
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(filename));

            for (int i = 0; i < matrix.length; i++) {
                for (int j = 0; j < matrix[i].length; j++) {
                    if(j == matrix[i].length - 1){    
                        bw.write(String.format("%.4f",matrix[i][j]));
                    } else{
                        bw.write(String.format("%.4f",matrix[i][j]) + ",");
                    }        
                }
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (IOException e) {
            //why does the catch need its own curly?
        }
    }
    
    /*
     * Method to calculate the genetic distance among the individuals in a population
     */
    public static void calcGeneticDistMatrixPopWide(List<String[]> genotypedSitesList, String outputFileName){
        
        int rows = genotypedSitesList.size();
        int cols = genotypedSitesList.size();
        double matrix[][] = new double[rows][cols];
        
        
        for(int i=0; i<rows; i++){
            
            for(int j=0; j<cols; j++){
                
                matrix[i][j] = calcNeiGeneticDistance(genotypedSitesList.get(i), genotypedSitesList.get(j));
                
                }
           }
        
        writeMatrix(outputFileName,matrix);         
    }
    
    /*
     * main method to run calcGeneticDistMatrixPopWide
     */
    public static void main(String[] args) throws IOException {
        
        String composedInputFN = "./DerivedPop1_FromAnc1_FuncGenotypes.txt";
        String composedOutputFN = "./GeneticDistanceMatrix_DerivedPop1_FromAnc1.txt";
        
        List<String[]> genotypedSitesList = getGenotypedSitesList(composedInputFN);
        
        calcGeneticDistMatrixPopWide(genotypedSitesList, composedOutputFN);
                
    }
    
    
}
