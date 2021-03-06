import java.util.Arrays;
import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.LinkedList;
import java.util.ArrayList;

public class FastBuckets2 {
	int l_param;
	int tuples;
	int dims;
	int bucketSize; // = (tuples/l_param);
	Map<Short, ArrayList<Integer>> sTuples; //SA to tuples'index in map[].
	Map<Integer, ArrayList<Integer>> changeSA; //points of new SA per bucket.
	short[] sArray;			// [|SA|]
	short[][] map;			// [n][dims]
	short[][][] buckets;	// [l][(n/l)][dims]
	double cost;
	static byte[] cardinalities = {79, 2, 17, 6, 9, 10, 83, 51};
	double threshold = 0.0;
	

	public FastBuckets2(int l_param, int tuples, int dims, short[][] map,
				   short[][][] buckets, double th) {
		this.l_param = l_param;
		this.tuples = tuples;
		this.dims = dims;
		this.map = map;
		this.buckets = buckets;
		bucketSize = (tuples/l_param);
		this.threshold = th;
	}

	public int getL() {
		return l_param;
	}

	public int getSize() {
		return tuples;
	}
	
	public void quickSortSA(int left, int right) {
		int i = left, j = right;
		int ref = left + (right-left)/2; //i.e., (right+left)/2
		int pivot = sTuples.get(sArray[ref]).size();
		short temp2;
		while (i <= j) {
			//DESCENDING ORDER:
			while (pivot < sTuples.get(sArray[i]).size())
				i++;
			while (sTuples.get(sArray[j]).size() < pivot) 
				j--;
			if (i <= j) {
				temp2=sArray[i];
				sArray[i]=sArray[j];
				sArray[j]=temp2;
				i++;
				j--;
			}
		}
		// recursion
		if (left < j)
			quickSortSA(left, j);
		if (i < right) {
			quickSortSA(i, right);
		}
	}
	
	//returns the most empty bucket
	public int getEmptyBucket(int[] empty){
		int best=0;
		for(int i=1; i<l_param; i++){
			if(empty[i] > empty[best])
				best = i;
		}
		return best;
	}
	
	//returns the next most empty bucket, after currentB
	//or -1 if no other bucket is empty.
	public int getNextEmptyBucket(int[] empty, int currentB, int remains){
		int best=0;
				
		for(int i=1; i<l_param; i++) {
			if (empty[i] > empty[best]){
				if ( (empty[i] < empty[currentB]) ||
					( (empty[i] == empty[currentB]) && (i>currentB) ) ) {
					best = i;
				}
			}
		}
		
		if ((empty[best] == 0) || (empty[best] > empty[currentB])) {
			//System.out.println(" - ");
			return -1;
		}
		
		if (best == currentB){
			//System.out.println(" -- ");
			return -1;
		}
		
		if ( (empty[currentB] >= remains) && (empty[best] < remains) ){ //*0.9
			//System.out.println(" ---- ");
			return -1;
		}
	
		return best;
	}
	
	public double getAvgNCP(int[] empty, ArrayList<Integer> tupSet, int offset, int b){
		double cost = 0.0;
		int canFit = tupSet.size() - offset;
		if (empty[b] < canFit) {
			canFit = empty[b];
		}
		for(int i=offset; i<(offset+canFit); i++){
			for(int j=0; j<(bucketSize - empty[b]); j++){
				for(int a=0; a<(dims-1); a++) {
					//cost += (double)Math.abs(buckets[b][j][a] - map[tupSet.get(i)][a])/(double)(cardinalities[a]-1);
					if (a==0 || a==2){
						cost += (double)Math.abs(buckets[b][j][a] - map[tupSet.get(i)][a]) / (double)(cardinalities[a]-1);
					}else{
						if (buckets[b][j][a] != map[tupSet.get(i)][a])
							cost+=(double)(1.0)/(double)(cardinalities[a]-1);
					}
				}
			}
		}
		double denom = (double)((bucketSize - empty[b]) * canFit * (dims-1));
		return cost/denom;
	}

	public double getMaxNCP(int[] empty, ArrayList<Integer> tupSet, int offset, int b){
		double cost = 0.0;
		double maxcost = 0.0;
		int canFit = tupSet.size() - offset;
		if (empty[b] < canFit) {
			canFit = empty[b];
		}
		for(int i=offset; i<offset+canFit; i++){
			for(int j=0; j<(bucketSize - empty[b]); j++){
				cost = 0.0;
				for(int a=0; a<(dims-1); a++) {
					//cost = (double)Math.abs(buckets[b][j][a] - map[tupSet.get(i)][a])/(double)(cardinalities[a]-1);
					if (a==0 || a==2){
						cost += (double)Math.abs(buckets[b][j][a] - map[tupSet.get(i)][a])/(double)(cardinalities[a]-1);
					}else{
						if (buckets[b][j][a] != map[tupSet.get(i)][a])
							cost += (double)(1.0)/(double)(cardinalities[a]-1);
					}
				}
				if (cost > maxcost)
					maxcost = cost;
			}
		}
		return maxcost/((double)(dims-1));
	}

	public void printBuckets(){
		// Print tuples in buckets:
		for(int b=0; b<l_param; b++){
			System.out.println("bucket="+b);
			for(int t=0; t<bucketSize; t++){
				System.out.print("[");
				for(int a=0; a<dims; a++){
					System.out.print(" "+buckets[b][t][a]);
				}
				System.out.println(" ]");
			}
		}
	}
	
	public short[][][] bucketization(int SA){
		System.out.println("1. check if SA is in dims range");
		if ((SA<0)||(SA>=dims)) {
			//SA is out of range
			System.err.println("Error: "+SA+" is out of range [0-"+dims+"]");
			//SA = dims;
			System.exit(-1);
		}

		//count frequencies:
		System.out.println("2. count frequencies");
		sTuples = new HashMap<Short, ArrayList<Integer>>();
		int maxFreq = 1;
		for (int i = 0; i < tuples;i++) {
			short value = map[i][SA];
			if (sTuples.containsKey(value)) {
				int tempFreq = sTuples.get(value).size() + 1;
				sTuples.get(value).add(i);
				if (maxFreq < tempFreq){
					maxFreq = tempFreq;
				}
			}else{
				ArrayList<Integer> tl = new ArrayList<Integer>();
				tl.add(i);
				sTuples.put(value, tl);
			}
		}
		
		int frSize = sTuples.size();
		
		//l-eligibility:
		System.out.println("4. l-eligibility");
		if (frSize < l_param){
			System.err.println("Error: There are less than "+l_param+" SA values");
			System.exit(-1);
		}
		if (maxFreq > bucketSize){ //(tuples / l_param)
			System.err.println("Error: max SA frequency("+maxFreq
							   +") is > than n/l ("+bucketSize+")");
			System.exit(-1);
		}
		
		//Sort SA values wrt frequency:
		sArray = new short[frSize];
		int ii = 0;
		for (Map.Entry<Short, ArrayList<Integer>> entry : sTuples.entrySet()) {
			sArray[ii++] = entry.getKey();
		}
		quickSortSA(0, frSize-1);

		//buckets:
		int offset = 0;
		int[] empty = new int[l_param];
		int freq = 0;
		int bucketOffset = 0;
		changeSA = new HashMap<Integer, ArrayList<Integer>>();
		
		//first l most frequent SA values:
		System.out.println("5. first l most frequent SA values");
		for (int i=0; i<l_param; i++){
			ArrayList<Integer> tmp = sTuples.get(sArray[i]);
			freq = tmp.size();
			for (int j=0; j<freq; j++){
				buckets[i][j] = map[tmp.get(j)];
			}
			empty[i] = bucketSize - freq;
			offset++;
			/*
			System.out.println("SA="+(map[offset][SA])+"\tb="+i
								+" fits="+freq+"/"+freq);
			 */
			ArrayList<Integer> chg = new ArrayList<Integer>();
			chg.add(freq-1);
			changeSA.put(i, chg); //(bucket, last_indx_of_this_SA)
		}
		
		//remaining SA values:
		System.out.println("6. remaining SA values");
		int fits = 0;
		int remains = 0;
		int v = l_param; //the (l+1) SA value.
		while(v < frSize){
			ArrayList<Integer> tmp = sTuples.get(sArray[v++]);
			freq = tmp.size();
			remains = freq;
			offset = 0;
			while(remains > 0){
				int i = getEmptyBucket(empty);
				int j = i;
				int best = i;
				double cost = getMaxNCP(empty, tmp, offset, j); //or getAvgNCP
				double maxcost = cost;
				//System.out.println(j+" "+cost);
				while (cost < threshold){
					j = getNextEmptyBucket(empty, j, remains);
					if (j == -1){ //not found.
						j = i; //most empty bucket
						break;
					}
					cost = getMaxNCP(empty, tmp, offset, j); //or getAvgNCP
					//System.out.println(j+" "+cost);
					if (cost > maxcost){
						maxcost = cost;
						best = j;
					}
				}
				//i = j;
				i  = best;
				cost = getMaxNCP(empty, tmp, offset, i); //or getAvgNCP
				//System.out.println("finally: "+i+" "+cost);
				if(empty[i] < remains){
					fits = empty[i];
				}else{
					fits = remains;
				}
				bucketOffset = bucketSize - empty[i];
				for (int jj=0; jj<fits; jj++){
					buckets[i][bucketOffset + jj] = map[tmp.get(offset + jj)];
				}
				empty[i] -= fits;
				remains -= fits;
				offset += fits;
				/*
				System.out.println("SA="+(map[offset][SA])+"\tb="+i
								   +" fits="+fits+"/"+freq);
				*/
				changeSA.get(i).add(bucketOffset+fits-1); //(bucket, last_indx_of_this_SA)
			}
		}
		
		//printBuckets(); //For debugging.
		
		//sTuples.clear();
		
		
		System.out.println("7. return buckets");
		return buckets;
	}
}
